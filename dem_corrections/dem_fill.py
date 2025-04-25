"""
Ryan Pranantyo
EOS, April 2025

a script to filled in NaN from interpolated DEM,
also can extent the domain and make it a low-resolution DEM.
all will be done through GMT
"""
import os, sys
import xarray as xr
import rioxarray
from pathlib import Path

### background DEM used to fill in NaN
bg_dem = "/scratch/ignatius.pranantyo/DATA/TanahAirIndonesia/BATNAS_v1-5-0.grd"

### (foreground) DEM to be filled in
fg_dem = Path("/home/ignatius.pranantyo/DATA/working_deltadtm/blended__batnas_deltadtm__correctedby-0.0m__res1s/blend_dem.grd")
fg_dir = fg_dem.parent

### config for outputs
# if want to extent the domain and make it low-resolution
# maximum domain extension is following the bg_dem bbox
extent_domain = True
new_extent_domain = [100.000, 120.000, -14.000, -2.000]
new_resolution = '9s' ## '1s=~30m, 3s=~90m, 9s=~270m'

###main program
print("JUST FILL IN NaN without extent the domain and resample the resolution ...")
tmp_new_fg_dem = os.path.join(fg_dir, f'tmp__{fg_dem.name[:-4]}_filled.grd')
new_fg_dem = os.path.join(fg_dir, f'{fg_dem.name[:-4]}_filled.grd')
new_fg_dem_tif = os.path.join(fg_dir, f'{fg_dem.name[:-4]}_filled.tif')
cmd = f"gmt grdfill {fg_dem} -Ag{bg_dem} -G{tmp_new_fg_dem}=nf"
os.system(cmd)

### clip
cmd = f"gmt grdclip {tmp_new_fg_dem} -Sa30/30 -G{new_fg_dem}"
os.system(cmd)

dem = xr.open_dataset(new_fg_dem)
#dem = dem.rio.write_crs('epsg:4326')
#dem.rio.to_raster(new_fg_dem_tif, driver='GTiff', compress='LZW')

oxmin, oxmax, oymin, oymax = dem.lon.min().values, dem.lon.max().values, dem.lat.min().values, dem.lat.max().values

### clean-up
#os.system(f"rm -f {tmp_new_fg_dem}")

#sys.exit()
if extent_domain == True:
    print(f"FILL IN NaN, extent the domain, and resample the resolution ...")
    ### resample first
    tmp_res = os.path.join(fg_dir, 'tmp_res.grd')
    cmd = f"gmt grdsample {tmp_new_fg_dem} -R{oxmin}/{oxmax}/{oymin}/{oymax} -I{new_resolution}/{new_resolution} -G{tmp_res}"
    os.system(cmd)

    ### extent the domain
    xmin = new_extent_domain[0]
    xmax = new_extent_domain[1]
    ymin = new_extent_domain[2]
    ymax = new_extent_domain[3]

    tmp_ext = os.path.join(fg_dir, 'tmp_ext.grd')
    cmd = f"gmt grdcut {tmp_res} -G{tmp_ext} -R{xmin}/{xmax}/{ymin}/{ymax} -N"
    os.system(cmd)

    tmp_fg = os.path.join(fg_dir, f'tmp_fill.grd')
    new_fg_dem = os.path.join(fg_dir, f'{fg_dem.name[:-4]}_newextent_resolution{new_resolution}.grd')
    new_fg_dem_tif = os.path.join(fg_dir, f'{fg_dem.name[:-4]}_newextent_resolution{new_resolution}.tif')
    cmd = f"gmt grdfill {tmp_ext} -Ag{bg_dem} -G{tmp_fg}=nf"
    os.system(cmd)

    ###clip
    cmd = f"gmt grdclip {tmp_fg} -Sa30/30 -G{new_fg_dem}"
    os.system(cmd)

    dem = xr.open_dataset(new_fg_dem)
    dem = dem.rio.write_crs('epsg:4326')
    dem.rio.to_raster(new_fg_dem_tif, driver='GTiff', compress='LZW')

    ### clean-up
    os.system(f"rm -f {tmp_res} {tmp_ext} {tmp_fg}")

os.system(f"rm -f {tmp_new_fg_dem}")

