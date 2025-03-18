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
parser.add_argument("--save_unit_as", type=str, default="cm",
        help="format of the values will be saved as 'm' (float16) or 'cm' and 'mm' (integer)")
args = parser.parse_args()

if args.save_unit_as.lower() == "cm":
    print(f"convert to cm and int32")
    scale = 100
    dtype = "int32"
elif args.save_unit_as.lower() == "mm":
    print(f"convert to mm and int32")
    scale = 1000
    dtype = "int32"
else:
    print(f"convert to float16")
    scale = 1.
    dtype = "float16"

nc = xr.open_dataset(args.jagurs_nc)

### drop vars
nc = nc.drop_vars(["step", "max_velocity"])

var2compress = ["initial_displacement", "max_height", "wave_height"]
for var in var2compress:
    nc[var].data = (nc[var].data * scale).astype(dtype)
    nc[var].attrs.update(description = 
        f"Ouput values from JAGURS converted into {args.save_unit_as.lower()} and {dtype} format by Ryan Pranantyo")

nc.to_netcdf(f"{args.jagurs_nc[:-3]}_tmp.nc")
nc.close()
os.system(f"mv {args.jagurs_nc[:-3]}_tmp.nc {args.jagurs_nc}")
