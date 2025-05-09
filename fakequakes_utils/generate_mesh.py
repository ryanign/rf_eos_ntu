"""
Ryan Pranantyo
EOS, May 2025

script to generate 2d rectangular mesh, 
intended to build a finite-fault model for earthquake model
"""
import os, sys
import meshio
import gmsh
import pygmsh
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

"""
this is example
"""
slab_gdf = gpd.read_file('/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/input_files__SouthernJava/contours/PUSGEN2017__WestCentralJavaMegathrust.shp')

top_surface = slab_gdf[slab_gdf['level'] == 0.]
bottom_surface = slab_gdf[slab_gdf['level'] == 60.]

top_coord = top_surface.get_coordinates()
top_coord['z'] = 0.0
top_coord = top_coord.reset_index()
bottom_coord = bottom_surface.get_coordinates()
bottom_coord['z'] = 0.0 # 60.0
bottom_coord = bottom_coord.reset_index()

tmp_df = pd.concat([top_coord, bottom_coord], ignore_index=True)

### initiate pygmsh
geometry = pygmsh.geo.Geometry()
model = geometry.__enter__()

### add points
mesh_size = 0.1
### top points
top_points = []
for ii in top_coord.index:
    top_points.append(model.add_point((top_coord['x'][ii], top_coord['y'][ii], top_coord['z'][ii]) , mesh_size = mesh_size))

### bottom points
bottom_points = []
for ii in bottom_coord.index[::-1]:
    bottom_points.append(model.add_point((bottom_coord['x'][ii], bottom_coord['y'][ii], bottom_coord['z'][ii]) , mesh_size = mesh_size))

### collect all points
### later, we can add points from level 10 km, 20 km, etc.
points = np.hstack(( top_points , bottom_points ))

### creating area to mesh
channel_lines = []
# connect the top first
for ii in range(len(top_points) - 1):
    channel_lines.append(model.add_line(top_points[ii], top_points[ii + 1]))
# connect side
channel_lines.append(model.add_line(top_points[-1] , bottom_points[0]))
# connect the bottom
for ii in range(len(bottom_points) - 1):
    channel_lines.append(model.add_line(bottom_points[ii], bottom_points[ii + 1]))
# connect the second side
channel_lines.append(model.add_line(bottom_points[-1] , top_points[0]))

## create a line loop
channel_loop = model.add_curve_loop(channel_lines)

## create plane surface for meshing
plane_surface = model.add_plane_surface(channel_loop)

# Call gmsh kernel before add physical entities
model.synchronize()

geometry.generate_mesh(dim=2)
gmsh.write("mesh.msh")

#mesh = meshio.read("mesh.msh")

#gmsh.clear()
#geometry.__exit__()
    
