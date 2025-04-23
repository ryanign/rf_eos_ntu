"""
Ryan Pranantyo
EOS, April 2025

a script to shift topo by X-m then mask it by a polygon
- a manually prepared polygon
"""
import os, sys
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import rioxarray
from pathlib import Path
from joblib import Parallel, delayed

def main(topo_file, mask_land_shp, land_correction,
         bathy_file, mask_sea_shp,
         where_to_save):
    """
    to clip and correct topography inland
    to clip bathymetry on the sea followoing topography extent
    """
    topo_file = Path(topo_file)
    bathy_file = Path(bathy_file)

    ### inland stage
    print('working on topography ...')
    dem_topo = rioxarray.open_rasterio(topo_file)
    mask_land = gpd.read_file(mask_land_shp)
    
    #clip/mask dem
    dem_topo_c = dem_topo.rio.clip(mask_land.geometry, mask_land.crs)

    #shift dem by correction values
    dem_topo_c = dem_topo_c.where(dem_topo_c >= 30 , dem_topo_c + land_correction)

    #replace below 0 with NaN
    dem_topo_c = dem_topo_c.where(dem_topo_c >= 0, np.nan)

    #temp save dem_topo_c
    temp_topo_fout = os.path.join(where_to_save, f'temp_{topo_file.name}')
    dem_topo_c.rio.to_raster(temp_topo_fout, driver='GTiff', compress='LZW')
    #convert to xyz
    temp_topo_xyz = os.path.join(where_to_save, f'temp_{topo_file.name[:-4]}.pts')
    cmd = f'gdal_translate {temp_topo_fout} -of XYZ {temp_topo_xyz}'
    os.system(cmd)
    
    #reshape xyz
    df = pd.read_table(temp_topo_xyz, sep=' ', header=None)
    df = df.dropna()
    pts_topo_xyz = os.path.join(where_to_save, f'clip-adjust__{topo_file.name[:-4]}.pts')
    df.to_csv(pts_topo_xyz, header=None, index=False, sep=' ', float_format="%.5f")

    ### water stage
    print('working on bathymetry ...')
    dem_bathy = rioxarray.open_rasterio(bathy_file)
    mask_sea = gpd.read_file(mask_sea_shp)

    #clip/mask bathy based on dem_topo bbox
    dem_bathy_c = dem_bathy.rio.clip_box(minx = dem_topo.x.min(),
                                         maxx = dem_topo.x.max(),
                                         miny = dem_topo.y.min(),
                                         maxy = dem_topo.y.max())
    #clip/mask bathy based on mask_sea
    dem_bathy_c = dem_bathy_c.rio.clip(mask_sea.geometry, mask_sea.crs)

    #remove dem_bathy_c above -2m to avoid 'magic island'
    dem_bathy_c = dem_bathy_c.where(dem_bathy_c < -2, np.nan)

    #temp save dem_bathy_c
    temp_bathy_fout = os.path.join(where_to_save, f'temp_{topo_file.name[:-4]}_{bathy_file.name}')
    dem_bathy_c.rio.to_raster(temp_bathy_fout, driver='GTiff', compress='LZW')
    #convert to xyz
    temp_bathy_xyz = os.path.join(where_to_save, f'temp_{topo_file.name[:-4]}_{bathy_file.name[:-4]}.pts')
    cmd = f'gdal_translate {temp_bathy_fout} -of XYZ {temp_bathy_xyz}'
    os.system(cmd)

    #reshape xyz
    df = pd.read_table(temp_bathy_xyz, sep=' ', header=None)
    df = df.dropna()
    pts_bathy_xyz = os.path.join(where_to_save, f'clip__{bathy_file.name[:-4]}_at_{topo_file.name[:-4]}.pts')
    df.to_csv(pts_bathy_xyz, header=None, index=False, sep=' ', float_format = "%.5f")

    ### cleaning up
    cmd = f'rm -f {temp_bathy_fout} {temp_bathy_xyz} {temp_topo_fout} {temp_topo_xyz}'
    os.system(cmd)

    return dem_topo, mask_land, dem_topo_c, dem_bathy, dem_bathy_c

if __name__ == "__main__":
    topo_files = ['/home/ryan/OneDrive_NTU_Projects/DATA/DeltaDTM/DeltaDTM_v1_1_S08E108.tif',
                  '/home/ryan/OneDrive_NTU_Projects/DATA/DeltaDTM/DeltaDTM_v1_1_S08E109.tif']
    topo_file = '/home/ryan/OneDrive_NTU_Projects/DATA/DeltaDTM/DeltaDTM_v1_1_S08E108.tif'
    mask_land_shp = '/home/ryan/OneDrive_NTU_Projects/DATA/DEM_corrections/sedang_dikerjakan/vectors/polygon_garispantai_rbi__buffered-100m.shp'

    shift_by = -0.5   ### correction value

    bathy_file = '/home/ryan/OneDrive_NTU_Projects/DATA/BATNAS/BATNAS_v1-5-0.tif'
    mask_sea_shp = '/home/ryan/OneDrive_NTU_Projects/DATA/DEM_corrections/sedang_dikerjakan/vectors/mask_sea_jawabalilombok__buffered1000m.shp'

    where_to_save = Path('/home/ryan/OneDrive_NTU_Projects/DATA/DEM_corrections/sedang_dikerjakan/raster_deltadem_correction/')

    Parallel(n_jobs = 2)(delayed(main)(topo_files[ii], mask_land_shp, shift_by, bathy_file, mask_sea_shp, where_to_save) for ii in range(len(topo_files)))
    #dem_topo, mask_land, dem_topo_c, dem_bathy, dem_bathy_c = main(topo_file, mask_land_shp, shift_by,
    #                            bathy_file, mask_sea_shp, where_to_save)

    ### checking temp outputs
    #dem_topo_c.rio.to_raster('/home/ryan/OneDrive_NTU_Projects/DATA/DEM_corrections/sedang_dikerjakan/raster_deltadem_correction/clip-and-shift__DeltaDTM_v1_1_S08E108.tif', driver='GTiff', compress='LZW')

    #dem_bathy_c.rio.to_raster('/home/ryan/OneDrive_NTU_Projects/DATA/DEM_corrections/sedang_dikerjakan/raster_deltadem_correction/clip-and-shift__BATNAS_v1-5_S08E108.tif', driver='GTiff', compress='LZW')
