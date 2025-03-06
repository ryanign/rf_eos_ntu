"""
Ryan Pranantyo
EOS, February 2025

this is a basic script to launch stochastic simulation for Method 1
1. Prepare slip models using plotting_stochastic_sources.py / convert_unit_sources__to__grd.py,
2. Once (1) is done, then execute this script.
"""
import os, sys
import numpy as np
import pandas as pd
import argparse
from pathlib import Path

###
jagurs_template_folder = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/jagurs_runs/trial_template"
displacement_grd_file  = Path("/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/OUTPUTS/displacement_ready2use/stochastic_sources__Mw_7.900000__Lon_104.104800__Lat_-7.122200__table_simplified_10.grd")
where_to_run = Path("/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/jagurs_runs/trial_runs")

###
curdir = os.getcwd()
where_to_run.mkdir(exist_ok = True)

### scenario name
model_name = displacement_grd_file.name[:-4]
print(f"Going to launch {model_name} ...")

outdir = Path(os.path.join(where_to_run, f"jagurs__{model_name}"))
outdir.mkdir(exist_ok = True)
os.chdir(outdir)
cmd = f"rm -rf *"
os.system(cmd)

cmd = f"cp -P {jagurs_template_folder}/*.nc ."
os.system(cmd)

cmd = f"cp -P {jagurs_template_folder}/*.grd ."
os.system(cmd)

cmd = f"ln -s {jagurs_template_folder}/tsun.par ."
os.system(cmd)

cmd = f"ln -s {jagurs_template_folder}/gridfile.dat ."
os.system(cmd)

cmd = f"cp -f {jagurs_template_folder}/basic_jagurs.pbs jagurs.pbs"
os.system(cmd)

cmd = f"cp -P {jagurs_template_folder}/jagurs_* ."
os.system(cmd)

cmd = f"ln -s {displacement_grd_file} disp__g00.grd"
os.system(cmd)

### update basic_jagurs.pbs

### submit job
cmd = f"qsub jagurs.pbs"
os.system(cmd)

os.chdir(curdir)
