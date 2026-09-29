"""
Ryan Pranantyo
EOS, April 2025

going to take average sea surface height above sea level for Indonesia 
based on the ECMWF dataset
"""
import os, sys
import numpy as np
import xarray as xr
from pathlib import Path
from joblib import Parallel

bbox = [94,150,-12,8]

ecmwf_path = '/scratch/ignatius.pranantyo/DATA/ECMWF_GlobalSeaLevel'

where_to_save = Path('/home/ignatius.pranantyo/DATA/GlobalSeaLevel/Indonesia')
where_to_save.mkdir(exist_ok = True)

infile = os.path.join(ecmwf_path, 'dt_global_twosat_phy_l4_202312_vDT2024-M01.nc')

# load sea surface height above sea level only
nc = xr.open_dataset(infile)['sla']

# assign projection
nc = nc.rio.write_crs('epsg:4326')

# clip to bbox
nc = nc.rio.clip_box(minx=bbox[0],
        maxx=bbox[1],
        miny=bbox[2],
        maxy=bbox[3])

fout = os.path.join(where_to_save, 'Indonesia__202312.tif')
nc.rio.to_raster(fout, driver='GTiff', compress='LZW')
nc.close()

sys.exit()
msl = nc['sla']
nc.close()

