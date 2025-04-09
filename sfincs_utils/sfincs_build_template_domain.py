"""
Ryan Pranantyo
EOS, March 2025

a script to automate the process of building template domain model for SFINCS
this is based on my older script prepared at Reask in March 2024
-- thank you Reask :) --
"""
import os, sys
import numpy as np
import pandas as pd
import geopandas as gpd

from hydromt_sfincs import SfincsModel
from hydromt_sfincs import utils
from hydromt import DataCatalog

from pathlib import Path


### loading / reading input files using argparser
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--output_dir", type = str,
        default = "/home/ignatius.pranantyo/Tsunamis/Learnings/SFINCS_examples/domain_exercise",
        help = "where to save template domain")
parser.add_argument("--domain_name", type = str,
        default = "Pangandaran_example",
        help = "name of template domain")
parser.add_argument("--domain_area", type = str,
        default = "/home/ignatius.pranantyo/Tsunamis/Learnings/SFINCS_examples/input_exercise/domain_pangandaran.gpkg",
        help = "domain polygon in .gpkg format")
parser.add_argument("--model_resolution", type = int,
        default = 30,
        help = "model resolution to be generated, in m")
parser.add_argument("--dem_file", type = str,
        default = "/home/ignatius.pranantyo/Tsunamis/Learnings/SFINCS_examples/input_exercise/dem_pangandaran_small.tif",
        help = "DEM in .tif file")
parser.add_argument("--masking_area", type = str,
        default = "/home/ignatius.pranantyo/Tsunamis/Learnings/SFINCS_examples/input_exercise/masking_area__pangandaran_v2.gpkg",
        help = "polygon to define active cells for simulations, in .gpkg format")
parser.add_argument("--manning_land", type = float,
        default = 0.04,
        help = "Manning's roughness coefficient for land")
parser.add_argument("--manning_sea", type = float,
        default = 0.02,
        help = "Manning's roughness coefficient for sea")
parser.add_argument("--open_bc_line", type = str,
        default =  "/home/ignatius.pranantyo/Tsunamis/Learnings/SFINCS_examples/input_exercise/open_bc_line__pangandaran.gpkg",
        help = "line where open BC water level is")
parser.add_argument("--walltime", type = int, default = 24,
        help = "wall time requested, in hour and int")
parser.add_argument("--ncpus", type = int, default = 4,
        help = "number of cpus requested")
parser.add_argument("--sfincs_exe", type = str,
        default = "/home/ignatius.pranantyo/apps/sfincs/executable/bin/sfincs",
        help = "path to sfincs exe file, will be used to create a symbolic link to template domain")
parser.add_argument("--project_code", type=str, 
        default = "eos_luca.dalzilio",
        help = "Project code name on WildFly")
parser.add_argument("--qtype", type=str, 
        default = "qintel_wfly",
        help = "queue type in Wildfly: qintel_wfly, qamd_wfly, dev")

args = parser.parse_args()
### 

### below is main script ###
outdir = args.output_dir 
domain = args.domain_name
domain_poly = args.domain_area
model_res = args.model_resolution
dem_input = args.dem_file
open_bc_line = args.open_bc_line
mask_poly = args.masking_area
manning_land = args.manning_land
manning_sea = args.manning_sea


print(f"=========================================")
print(f"Building domain model = {domain}")
print(f"=========================================")
# where to save the template #
root_dir = os.path.join(outdir, domain, f"resolution__{model_res}_m")

sf = SfincsModel(root = root_dir, mode = "w+")
sf.setup_grid_from_region(
        region = {"geom" : domain_poly},
        res = model_res,
        rotated = False,
        #crs = "4326"       ### trying with epsg:4326 so that I do not need to make any projection conversion
        )

# load elevation data
print(f"  > Reading DEM = {dem_input}")
da = sf.data_catalog.get_rasterdataset(
        dem_input,
        variable = "elevtn",
        geom = sf.region,
        meta = {"version" : "1"},
        buffer = 1000,
        )

datasets_dep = [{"elevtn" : da}]
dep = sf.setup_dep(datasets_dep = datasets_dep)

### if masking area is available
#print(f"  > Polygon masking from {mask_poly}")
sf.setup_mask_active(include_mask = mask_poly, reset_mask = True)

### instead of using masking area, we can define this based on the zmin and zmax from the DEM

# setup manning's roughness
sf.setup_manning_roughness(
        manning_land = manning_land,
        manning_sea = manning_sea,
        rgh_lev_land = 0
        )

# read open BC polylin
# this is to show where the water level BC is
bc_waterlevel = sf.data_catalog.get_geodataframe(open_bc_line)


### if masking area is unavailable
#sf.setup_mask_bounds(btype="waterlevel", zmax = mask_zmin + 1, reset_bounds = True)
### if masking area is available
sf.setup_mask_bounds(btype="waterlevel", zmax = -5, reset_bounds = True)

# read open BC points
# seems add this when add water level forcing

# write template domain 
fig_name = f"Domain__{domain}__Resolution_{model_res}_m"
fig,ax = sf.plot_basemap(fn_out = fig_name, bmap = "sat", zoomlevel = 6)
sf.write()

# write sfincs.slurm script
f = open(os.path.join(root_dir, "sfincs.slurm"), "w")
f.write(f"#!/bin/bash\n")
f.write(f"#PBS -N SFINCS\n")
f.write(f"#PBS -P {args.project_code}\n")
f.write(f"#PBS -q {args.qtype}\n")
f.write(f"#PBS -l walltime={args.walltime:02d}:00:00\n")
f.write(f"#PBS -l select=1:ncpus={args.ncpus}\n")
f.write(f"module load gnu/gcc-12.3\n")
f.write(f"module load hdf5/1.14.3-intel2023-parallel\n")
f.write(f"cd $PBS_O_WORKDIR\n")
f.write(f"time ./sfincs\n")
f.close()

# create a symlink from sfincs.exe
cmd = f"ln -s {args.sfincs_exe} {root_dir}/sfincs"
os.system(cmd)
# ADD THIS LATER!

print(f"")
print(f"=========================================")
print(f" Chek = {root_dir}")
print(f"=========================================")

