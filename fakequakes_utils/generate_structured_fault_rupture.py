"""
Ryan Pranantyo
EOS, May 2025

to create a rectangular mesh fault plane
with help from ChatGPT!
hehe.
thanks
"""


import geopandas as gpd
import numpy as np
import pandas as pd
import pygmsh
import meshio
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
from sklearn.linear_model import LinearRegression

# === PARAMETERS ===
contour_file = "your_contours.shp"  # <-- your input shapefile
output_csv = "fault_plane_mesh_pygmsh.csv"
dx = 10000  # 10 km
dy = 20000  # 20 km

# === STEP 1: LOAD CONTOURS ===
contours = gpd.read_file(contour_file)

# Reproject to UTM for metric units
utm_crs = contours.estimate_utm_crs()
contours = contours.to_crs(utm_crs)

# Get bounding box
minx, miny, maxx, maxy = contours.total_bounds

# === STEP 2: CREATE MESH WITH PYGMSH ===
with pygmsh.geo.Geometry() as geom:
    rectangle = geom.add_rectangle(xmin=minx, xmax=maxx, ymin=miny, ymax=maxy, z=0, mesh_size=dx)
    mesh = geom.generate_mesh()

# Extract mesh cell centers
points = mesh.points
cells = mesh.cells_dict.get("quad", mesh.cells_dict.get("triangle", None))  # quads preferred, fallback triangles

if cells is None:
    raise ValueError("Mesh does not have quadrilaterals or triangles.")

cell_centers = []
cell_polygons = []

for cell in cells:
    pts = points[cell]
    center = np.mean(pts[:, :2], axis=0)
    poly = Polygon(pts[:, :2])
    cell_centers.append(center)
    cell_polygons.append(poly)

cell_centers = np.array(cell_centers)

# === STEP 3: PROCESS EACH CELL ===
records = []
for center_xy, polygon in zip(cell_centers, cell_polygons):
    center_x, center_y = center_xy

    # Search points nearby
    search_radius = np.sqrt(dx**2 + dy**2) / 2

    nearby = []
    for geom, elev in zip(contours.geometry, contours['elevation']):
        if geom.distance(Point(center_x, center_y)) < search_radius:
            if geom.geom_type == 'LineString':
                for x, y in np.array(geom.coords):
                    nearby.append([x, y, elev])
            elif geom.geom_type == 'MultiLineString':
                for linestring in geom.geoms:
                    for x, y in np.array(linestring.coords):
                        nearby.append([x, y, elev])

    nearby = np.array(nearby)
    
    if len(nearby) < 3:
        continue  # Not enough points to fit a plane

    # Fit plane
    X = nearby[:, :2]
    y = nearby[:, 2]
    model = LinearRegression().fit(X, y)
    a, b = model.coef_

    dip_rad = np.arctan(np.sqrt(a**2 + b**2))
    dip_deg = np.degrees(dip_rad)

    strike_rad = np.arctan2(a, b)
    strike_deg = np.degrees(strike_rad)
    strike_deg = (90 - strike_deg) % 360  # compass convention

    mean_depth = np.mean(nearby[:, 2])

    # Reproject center back to lat/lon
    center_lonlat = gpd.GeoSeries([Point(center_x, center_y)], crs=utm_crs).to_crs(epsg=4326).geometry[0]

    records.append({
        'lon': center_lonlat.x,
        'lat': center_lonlat.y,
        'depth': mean_depth,
        'strike': strike_deg,
        'dip': dip_deg,
        'length': 10,  # km
        'width': 20,   # km
        'geometry': polygon
    })

# === STEP 4: SAVE OUTPUT ===
df = pd.DataFrame(records)
df_simple = df[['lon', 'lat', 'depth', 'strike', 'dip', 'length', 'width']]
df_simple.to_csv(output_csv, index=False)

print(f"Mesh saved as '{output_csv}'.")

# === STEP 5: PLOTTING FUNCTION ===
def plot_fault_mesh(df, utm_crs, color_by='dip', arrows=True):
    gdf = gpd.GeoDataFrame(df, geometry='geometry', crs=utm_crs)

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    gdf.plot(column=color_by, cmap='viridis', edgecolor='k', linewidth=0.5, legend=True, ax=ax)

    if arrows:
        for idx, row in gdf.iterrows():
            center = row.geometry.centroid
            strike_rad = np.radians(row['strike'])
            dx = np.cos(strike_rad)
            dy = np.sin(strike_rad)
            ax.arrow(center.x, center.y, dx*5000, dy*5000, head_width=2000, head_length=2000, fc='red', ec='red')

    ax.set_title(f"Fault Plane Mesh (colored by {color_by})")
    ax.set_xlabel('X (meters)')
    ax.set_ylabel('Y (meters)')
    plt.grid(True)
    plt.axis('equal')
    plt.show()

# === STEP 6: PLOT THE RESULT ===
plot_fault_mesh(df, utm_crs)

