"""
Ryan Pranantyo
EOS, 12 May 2026

A script to convert JAGURS output.nc file to a GTiff file.
A this moment, I only convert:
    > Maximum elevation tsunami footprint

Useage:
    python convert__jagurs2gtiff.py --jagurs_nc <path/to/output.nc> --what_to_convert <1> --scale_ratio 10000

    what_to_convert: 1. tsunami maximum footprint
"""
import os
import sys
import numpy as np
import pandas as pd
import xarray as xr
import rioxarray as rio
import argparse
from pathlib import Path

def convert_maximum_footprint(infile, dem_file, gtiff_dir, scale_ratio, crs, bottom_deformation):
    print(f'='*60)
    print(f' Converting maximum tsunami footprint to a GTiff file')
    print(f'  JAGURS output file = {infile}')
    print(f'='*60)

    nc = xr.open_dataset(infile, decode_times=False)

    # load maximum footprint
    max_footprint = nc['max_height'] / scale_ratio

    # load initial displacement
    init_disp = nc['initial_displacement'].sum(axis=0) / scale_ratio
    init_disp_ori = init_disp

    # load dem
    dem = xr.open_dataset(dem_file)
    nx = int(dem['dimension'].values[0])
    ny = int(dem['dimension'].values[1])
    x0,x1 = dem['x_range'].values
    y0,y1 = dem['y_range'].values

    lon_arr = np.linspace(x0, x1, nx)
    lat_arr = np.linspace(y0, y1, ny)
    zz = dem['z'].values.reshape(ny, nx).astype(np.float32)

    dem = xr.DataArray(
            data   = zz[::-1],            ### convert water to negative and land to positive
            dims   = ['lat', 'lon'],
            coords = {
                'lon': lon_arr,
                'lat': lat_arr,
                },
            name = 'DEM',
            )
    dem = dem.assign_coords(
            lat = init_disp.lat,
            lon = init_disp.lon,
            )

    # mask
    land_mask = dem > 0
    max_footprint = max_footprint.where(~land_mask)
    dem_mask = dem.where(~land_mask)
    init_disp = init_disp.where(~land_mask)
    
    #if bottom_deformation:
    #    flow_depth = max_footprint + dem_mask - init_disp
    #else:
    #    flow_depth = max_footprint + dem_mask

    flow_depth = max_footprint + dem_mask
   
    #flow_depth = only the inundation above land

    nc = flow_depth.astype(np.float32) #max_footprint #, low_depth

    # assign a projection
    nc = nc.rio.set_spatial_dims(x_dim='lon', y_dim='lat')
    nc.rio.write_crs(crs, inplace=True)

    # convert to a GTiff
    print(f'  Saving max footprint ...')
    fname = infile.name[:-3]
    gtiff_name = os.path.join(gtiff_dir, f'{fname}__max_footprint.tiff')
    nc.rio.to_raster(gtiff_name, driver="GTiff", compress="ZSTD", zstd_level=9,
            tiled=True, blockxsize=256, blockysize=256)
    print(f'  > check {gtiff_name}')

    # initial displacement
    print(f'  Saving initial displacement ...')
    print(np.nanmin(init_disp_ori), np.nanmax(init_disp_ori))
    nc = init_disp_ori.astype(np.float32)
    nc = nc.rio.set_spatial_dims(x_dim='lon', y_dim='lat')                                              
    nc.rio.write_crs(crs, inplace=True)
    gtiff_name = os.path.join(gtiff_dir, f'{fname}__init_disp.tiff')
    nc.rio.to_raster(gtiff_name, driver="GTiff", compress="ZSTD", zstd_level=9,
            tiled=True, blockxsize=256, blockysize=256)
    print(f'  > check {gtiff_name}')


    return nc
    #"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '--jagurs_nc', type=str,
            #default = '/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/Benchmarking__Historical/JAGURS__vs__SFINCS/2006Java__FujiiSatake2006__mod-10GPa__0-5m/SD00.nc',
            default = '/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/Benchmarking__Historical/JAGURS__vs__SFINCS/Tile-5__LargeEvents/DeltaDTM_0-5m__SFFM__Mw_8.900000__Lon_108.909670__Lat_-8.697910__unit_source_slip__5/SD04.nc',
            help = 'JAGURS nc output file',
            )
    parser.add_argument(
            '--dem_file', type=str,
            default = '/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/Benchmarking__Historical/JAGURS__vs__SFINCS/Tile-5__LargeEvents/DeltaDTM_0-5m__SFFM__Mw_8.900000__Lon_108.909670__Lat_-8.697910__unit_source_slip__5/DEM__SD04__NEG.grd',
            #'/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/Benchmarking__Historical/JAGURS__vs__SFINCS/2006Java__FujiiSatake2006__mod-10GPa__0-5m/DEM__SD00__NEG.grd',
            help = 'DEM used',
            )
    parser.add_argument(
            '--what_to_convert', type=int,
            default = 1,
            help = '1: tsunami maximum footprint',
            )
    parser.add_argument(
            '--consider_bottom_deformation',
            type=bool,
            default = False,
            help = 'True: JAGURS consider bottom deformation, False: otherwise and will not included in the conversion'
            )
    parser.add_argument(
            '--scale_ratio', type=float,
            default = 10000.,
            help = 'In order to save some space, JAGURS outputs is saved in an integer format, multiplied it by 10,000'
            )
    parser.add_argument(
            '--crs', type=str,
            default = 'epsg:4326',
            help = 'CRS projection to use in an epsg: format',
            )
    args = parser.parse_args()

    if os.path.exists(args.jagurs_nc):
        infile = Path(args.jagurs_nc)
        parent_dir = infile.parent
        gtiff_dir  = Path(os.path.join(parent_dir, 'GTiff'))
        gtiff_dir.mkdir(exist_ok = True)

        if args.what_to_convert == 1:
            nc = convert_maximum_footprint(
                    infile, 
                    args.dem_file, gtiff_dir, 
                    args.scale_ratio, args.crs,
                    args.consider_bottom_deformation,
                    )

        nc.close()
        
    else:
        print(f'{args.jagurs_nc} is not available!')
        sys.exit()


