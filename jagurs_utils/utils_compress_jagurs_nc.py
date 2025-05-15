"""
Ryan Pranantyo
EOS, March 2025

a small script to compress jagurs output file:
    > convert data to int32
    > remove not needed variables
"""
import os, sys
import xarray as xr
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--jagurs_nc", type=str, default="SD00.nc",
        help = "JAGURS nc output file")
parser.add_argument("--scale_ratio", type=int, default=10000,
        help = "output values will be multipled by this scale and save it as int")
args = parser.parse_args()

scale = args.scale_ratio
dtype = 'int32'

nc = xr.open_dataset(args.jagurs_nc)

### variables to be compressed
var2compress = ["initial_displacement", "max_height", "wave_height", "max_velocity"]

for var in var2compress:
    nc[var] = nc[var] * scale

encoding = {'wave_height'          : {'dtype' : dtype, '_FillValue' : 0, 'zlib' : True, 'complevel' : 5},
            'max_height'           : {'dtype' : dtype, '_FillValue' : 0, 'zlib' : True, 'complevel' : 5},
            'initial_displacement' : {'dtype' : dtype, '_FillValue' : 0, 'zlib' : True, 'complevel' : 5},
            'max_velocity'         : {'dtype' : dtype, '_FillValue' : 0, 'zlib' : True, 'complevel' : 5}}

nc.to_netcdf(f"{args.jagurs_nc[:-3]}_tmp.nc", encoding = encoding, mode = 'w')
nc.close()
os.system(f"mv {args.jagurs_nc[:-3]}_tmp.nc {args.jagurs_nc}")




sys.exit()
