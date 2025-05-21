"""
Ryan Pranantyo
EOS, May 2025

mesh generation using OCSMesh
based on my old work at Reask

main input files = 
    - slab dep in .tif
    - bounding polygon
"""
import os, sys
import numpy as np
from ocsmesh import Raster, Geom, Hfun, Mesh, JigsawDriver, utils
from shapely.geometry import shape
import fiona
USE_PYGEOS=0
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon

#import xarray as xr
import rasterio
from osgeo import gdal
from copy import deepcopy
import jigsawpy
from shapely import geometry

from pathlib import Path

def create_raster_list(dem_paths):
    rast_list = list()
    for f in dem_paths:
        rast_list.append(Raster(f, crs='EPSG:4326'))
    return rast_list

### slab depth
slab_depth = ["/home/ignatius.pranantyo/DATA/USGS__Slab2/sum_slab2_dep_02_23_18.tif"]

### domain of mesh
bbox_f = "/home/ignatius.pranantyo/Tsunamis/Learnings/OCSMesh_examples/bbox__slab2__jawa.shp"
base_gdf = gpd.read_file(bbox_f)

###
rast_list_1 = create_raster_list(slab_depth)

###
base_geom = Geom(rast_list_1,
                 base_shape = base_gdf.union_all(),
                 base_shape_crs = base_gdf.crs,)

rast_list_2 = create_raster_list(slab_depth)

base_hfun = Hfun(rast_list_2,
                 base_shape = base_gdf.union_all(),
                 base_shape_crs = base_gdf.crs,
                 hmax = 25000,
                 hmin = 25000,
                 method = "fast")

print("... generating mesh ...")
driver = JigsawDriver(geom=base_geom, hfun=base_hfun)
base_mesh = driver.run()
utils.reproject(base_mesh.msh_t, "EPSG:4326")

### save basic_mesh
where_to_save = "/home/ignatius.pranantyo/Tsunamis/Learnings/OCSMesh_examples/"
base_mesh.write(os.path.join(where_to_save, "basic_mesh.gr3"), format="grd", overwrite = True)

### interpolate
base_mesh = Mesh.open(os.path.join(where_to_save, "basic_mesh.gr3"), crs="EPSG:4326")
utils.reproject(base_mesh.msh_t, "EPSG:4326")
#rast_list = create_raster_list(slab_depth)
rast_list = [Raster(slab_depth[0])]
base_mesh.interpolate(rast_list, method="nearest", nprocs=4)

### saving
base_mesh.write(os.path.join(where_to_save, "mesh.gr3"), format="grd", overwrite = True)
base_mesh.write(os.path.join(where_to_save, "mesh.2dm"), format="2dm", overwrite = True)




"""
writing fault mesh for fakequakes
NO CENTROID NODE1   NODE2   NODE3   MEAN_LENGTH AREA   STRIKE DIP
   (x,y,z)  (x,y,z) (x,y,z) (x,y,z) (km)        (km^2) (deg)  (deg)
"""
#Earth radius, km
R_EARTH = 6371.0

def geo_to_ecef(lat, lon, depth_km):
    """
    convert geographic coordinates to ECEF
    """
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)
    r = R_EARTH - depth_km

    x = r * np.cos(lat_rad) * np.cos(lon_rad)
    y = r * np.cos(lat_rad) * np.sin(lon_rad)
    z = r * np.sin(lat_rad)
    return np.array([x,y,z])

def ecef_to_geo(x, y, z):
    """
    Convert ECEF (x, y, z) back to (lat, lon, depth) in degrees and km.
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    lat = np.degrees(np.arcsin(z / r))
    lon = np.degrees(np.arctan2(y, x))
    depth = R_EARTH - r
    return lon, lat, depth

def compute_strike_dip_area_centroid(p1, p2, p3):
    """
    Strike and dip computations and visualization scripts were refined with the help of OpenAI’s ChatGPT (2025)
    OpenAI. (2025). ChatGPT (May 2025 version). Retrieved from https://chat.openai.com
    """
    P1 = geo_to_ecef(p1[1], p1[0], p1[2])
    P2 = geo_to_ecef(p2[1], p2[0], p2[2])
    P3 = geo_to_ecef(p3[1], p3[0], p3[2])

    # Vectors in ECEF
    v1 = P2 - P1
    v2 = P3 - P1
    n = np.cross(v1, v2)
    n_unit = n / np.linalg.norm(n)

    # Centroid for local ENU basis
    centroid = (P1 + P2 + P3) / 3
    lon_c, lat_c, depth_c = ecef_to_geo(*centroid)

    # ENU basis at centroid
    lat0 = np.radians(lat_c)
    lon0 = np.radians(lon_c)

    e = np.array([-np.sin(lon0), np.cos(lon0), 0])
    n_ = np.array([-np.sin(lat0)*np.cos(lon0),
                   -np.sin(lat0)*np.sin(lon0),
                    np.cos(lat0)])
    u = np.array([np.cos(lat0)*np.cos(lon0),
                  np.cos(lat0)*np.sin(lon0),
                  np.sin(lat0)])

    # Normal in ENU
    n_E = np.dot(n_unit, e)
    n_N = np.dot(n_unit, n_)
    n_U = np.dot(n_unit, u)

    # Flip normal to point downward
    if n_U > 0:
        n_E, n_N, n_U = -n_E, -n_N, -n_U

    # Dip: angle between normal and vertical
    dip = np.degrees(np.arccos(-n_U))
    dip = max(0.0, min(90.0, dip))  # Clamp to [0, 90]

    # Dip direction azimuth
    dip_az = (np.degrees(np.arctan2(n_E, n_N))) % 360

    # Strike is dip direction - 90
    strike = (dip_az - 90) % 360

    # Ensure strike follows your model convention (≥ 200°)
    if strike < 180:
        strike += 180

    # Area & edge lengths
    area = 0.5 * np.linalg.norm(np.cross(v1, v2))
    mean_len = (np.linalg.norm(v1) + np.linalg.norm(v2) + np.linalg.norm(P3 - P2)) / 3

    return strike, dip, area, mean_len, (lon_c, lat_c, depth_c)


###############################################################################

x,y,z = base_mesh.coord[:,0], base_mesh.coord[:,1], base_mesh.value.flatten()
triangles = base_mesh.triangles

data = list()
data2 = list()

polys = []
for ii, tri in enumerate(triangles):
    print(tri)
    x1,y1,z1 = x[tri[0]], y[tri[0]], z[tri[0]]
    x2,y2,z2 = x[tri[1]], y[tri[1]], z[tri[1]]
    x3,y3,z3 = x[tri[2]], y[tri[2]], z[tri[2]]
    
    
    ### going to estimate strike, dip, area, mean length, and centroid
    P1 = (x1,y1,z1)
    P2 = (x2,y2,z2)
    P3 = (x3,y3,z3)

    strike, dip, area, mean_length, centroid = compute_strike_dip_area_centroid(P1, P2, P3)
    xc = centroid[0]
    yc = centroid[1]
    zc = centroid[2]

    row = [ii+1, xc, yc, zc, x1, y1, z1, x2, y2, z2, x3, y3, z3, mean_length, area, strike, dip]
    data.append(row)

    poly = [(x1,y1), (x2,y2), (x3,y3)]
    polys.append(poly)

df = pd.DataFrame(
        columns = ["NO", "Lon_c", "Lat_c", "Depth_c", "Lon1", "Lat1", "Depth1", "Lon2", "Lat2", "Depth2", "Lon3", "Lat3", "Depth3", "Mean_Length", "Area", "Strike", "Dip" ],
        data = data)

### save the mesh for fakequakes
cols2convert = df.columns[1:]
df[cols2convert] = df[cols2convert].map(lambda x: f"{x:.6f}")   #f[cols2convert].astype("float32")

mshout_fout = os.path.join(where_to_save, f"mesh.mshout")
df.to_csv(mshout_fout, header = None, sep="\t", index=False)

###plotting
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize

# Create subplots
fig, axes = plt.subplots(1, 2, figsize=(20, 6), constrained_layout=True)

# Subplot 1: Dip colormap
norm_dip = Normalize(vmin=0, vmax=60)
coll1 = PolyCollection(polys, array=df['Dip'].values, cmap='viridis', edgecolor='k', norm=norm_dip)
coll1.vmin = 5
coll1.vmax = 60
axes[0].add_collection(coll1)
#axes[0].tripcolor(polys, facecolors=df['Dip'], cmap='viridis', edgecolors='k', vmin = 0, vmax=60)
axes[0].autoscale()
axes[0].set_title("Dip (°)")
axes[0].set_xlabel("Longitude")
axes[0].set_ylabel("Latitude")
fig.colorbar(coll1, ax=axes[0], label='Dip (°)')

# Subplot 2: Strike colormap
norm_strike = Normalize(vmin=0, vmax=360)
coll2 = PolyCollection(polys, array=df['Strike'].values, cmap='twilight', edgecolor='k', norm=norm_strike)
coll2.vmin = 0
coll2.vmax = 360
axes[1].add_collection(coll2)
axes[1].autoscale()
axes[1].set_title("Strike (°)")
axes[1].set_xlabel("Longitude")
axes[1].set_ylabel("Latitude")
fig.colorbar(coll2, ax=axes[1], label='strike')

fout = os.path.join(where_to_save, "mesh.png")
fig.savefig(fout, dpi=300)
plt.close()
