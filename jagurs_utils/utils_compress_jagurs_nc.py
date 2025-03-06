"""
Ryan Pranantyo
EOS, March 2025

a small script to compress jagurs output file:
    > convert data to float16
    > remove not needed variables
"""
import os, sys
import xarray as xr
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--jagurs_nc", type=str, default="SD00.nc",
        help="JAGURS nc output file")
args = parser.parse_args()

nc = xr.open_dataset(args.jagurs_nc)

### drop vars
nc = nc.drop_vars(["step", "max_velocity"])

var2compress = ["initial_displacement", "max_height", "wave_height"]
for var in var2compress:
    nc[var].data = nc[var].data.astype("float16")
nc.to_netcdf(f"{args.jagurs_nc[:-3]}_tmp.nc")
nc.close()
os.system(f"mv {args.jagurs_nc[:-3]}_tmp.nc {args.jagurs_nc}")

