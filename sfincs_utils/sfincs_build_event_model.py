"""
Ryan Pranantyo
EOS, March 2025

a script to build an event model,
I am following my old script developed when I was in Reask.

> Read SFINCS template domain then copy to a new path
> Read basic configuration then add water level boundary condition forcing from JAGURS
"""
import os, sys
import numpy as np
import pandas as pd
import geopandas as gpd

from hydromt_sfincs import SfincsModel
from hydromt_sfincs import utils
from hydromt import DataCatalog

from pathlib import Path

### READING INPUT FILES ... ###
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--sfincs_template", type = str,
        default = "/home/ignatius.pranantyo/Tsunamis/Learnings/SFINCS_examples/domain_exercise/Pangandaran_example/resolution__30_m",
        help = "Prepared SFINCS template domain model")
parser.add_argument("--target_dir", type = str,
        default = "/home/ignatius.pranantyo/Tsunamis/Learnings/SFINCS_examples/run_exercise",
        help = "Where to run this event")
parser.add_argument("--event_name", type = str,
        default = "Pangandaran2006_testing",
        help = "Name of this run model")
parser.add_argument("--domain_name", type = str,
        default = "Pangandaran-small_testing")

parser.add_argument("--bc_waterlevel", type = str,
        default = "/home/ignatius.pranantyo/Tsunamis/Learnings/SFINCS_examples/input_exercise/extracted__SD01/timeseries__at_points__open_bc_line__pangandaran.csv",
        help = "Elevation timeseries extracted from JAGURS")
parser.add_argument("--bc_points_csv", type = str,
        default = "/home/ignatius.pranantyo/Tsunamis/Learnings/SFINCS_examples/input_exercise/points__open_bc_line__pangandaran.csv",
        help = "Location of water level boundary condition points")
parser.add_argument("--cfl", type = float,
        default = 0.5,
        help = "CFL criteria for SFINCS")
parser.add_argument("--advection", type = int,
        default = 2,
        help = "advection flag in SFINCS, 0 = NO, 1 = 1d, 2 = 2d")
parser.add_argument("--scale_ratio", type = int,
        default = 10000,
        help = "scale ratio to 'normalise' the timeseries value due to the conversion from JAGURS simulation to save some space. Please Use the same scale ratio")
# all done #
args = parser.parse_args()
sfincs_template = Path(args.sfincs_template)
sfincs_target = args.target_dir
event = args.event_name
domain = args.domain_name

print(f"====================================")
print(f"PREPARING FOLDER TO RUN SFINCS ")

# where to run #
root_dir = Path(os.path.join(sfincs_target, f"{domain}__{event}__{sfincs_template.name}"))
#os.system(f"rm -rf {root_dir}")
root_dir.mkdir(exist_ok = True)
cmd = f"cp -rPf {sfincs_template}/* {root_dir}/"
os.system(cmd)
print(f"  working directory = {root_dir}")

# read pre-configure SFINCS template
sf = SfincsModel(root = root_dir, mode = "r+")
sf.config['epsg']
print(sf.crs.utm_zone)
#sf_utm_zone = sf.crs.utm_zone

#sys.exit()
# setup open BC locations
bc_points_f = args.bc_points_csv
df_bc = pd.read_csv(bc_points_f)
### need to reproject lon lat from epsg:4326 to utm following sfincs template
pnts = gpd.points_from_xy(df_bc.LON.values, df_bc.LAT.values)
open_bnd = gpd.GeoDataFrame(index = df_bc.index + 1, geometry = pnts, crs = "epsg:4326")   ### original projection
open_bnd_utm = open_bnd.to_crs(sf.crs)                                                     ### converted to sf projection

# read water level bc
df_ts = pd.read_csv(args.bc_waterlevel)
df_ts["time"] = pd.to_datetime(df_ts["time"])
df_ts_time = df_ts["time"].dt.strftime("%Y%m%d %H%M%S")
time_start = df_ts["time"][0]
time_end = df_ts["time"].iloc[-1]
dtmaxout = int((time_end - time_start).total_seconds())

# setup sf config time
sf.setup_config(**{"tref" : time_start.strftime("%Y%m%d %H%M%S"),
                   "tstart" : time_start.strftime("%Y%m%d %H%M%S"),
                   "tstop" : time_end.strftime("%Y%m%d %H%M%S"),
                   "dtmaxout" : f"{dtmaxout}"})

# setup other parameters inside sfincs.inp
sf.setup_config(**
        {"dtout" : "60",
         "advection" : f"{args.advection}"} )
                

# add water level bc to sfincs config
timeS = pd.date_range(
        start = utils.parse_datetime(sf.config["tstart"]),
        end = utils.parse_datetime(sf.config["tstop"]),
        periods = len(df_ts_time)
        )

### NEED TO DOUBLE CHECK THE SCALE RATIO!!!
bzs = df_ts.drop(columns=["time"])
bzspd = pd.DataFrame(index = timeS, columns = df_bc.index.values + 1, data = bzs.values / args.scale_ratio)

# update sf config
sf.setup_waterlevel_forcing(timeseries = bzspd, locations = open_bnd_utm)
sf.forcing.keys()

sf.write()

### execute
os.chdir(root_dir)
os.system("qsub sfincs.slurm")
