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

import xarray as xr
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
slab_depth = ["/home/ryan/OneDrive_NTU_Projects/DATA/USGS/Slab2__SumatraJawa/sum_slab2_dep_02_23_18.tif"]

### domain of mesh
bbox_f = "/home/ryan/OneDrive_NTU_Projects/FaultMesh/input_files/bbox__slab2__jawa.shp"
base_gdf = gpd.read_file(bbox_f)

###
rast_list_1 = create_raster_list(slab_depth)

###
base_geom = Geom(rast_list_1,
                 base_shape = base_gdf.unary_union,
                 base_shape_crs = base_gdf.crs,)

rast_list_2 = create_raster_list(slab_depth)

base_hfun = Hfun(rast_list_2,
                 base_shape = base_gdf.unary_union,
                 base_shape_crs = base_gdf.crs,
                 hmax = 1000,
                 hmin = 1000,
                 method = "fast")

print("... generating mesh ...")
driver = JigsawDriver(geom=base_geom, hfun=base_hfun)
base_mesh = driver.run()
utils.reproject(base_mesh.msh_t, "EPSG:4326")

### interpolate to slab depth
rast_list_3 = create_raster_list(slab_depth)
base_mesh.interpolate(rast_lis)

where_to_save = "/home/ryan/OneDrive_NTU_Projects/FaultMesh/"
base_mesh.write(os.path.join(where_to_save, "mesh.gr3"), format="grd", overwrite = True)
base_mesh.write(os.path.join(where_to_save, "mesh.2dm"), format="2dm", overwrite = True)

### 
