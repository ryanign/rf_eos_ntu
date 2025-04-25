"""
Ryan Pranantyo
EOS, April 2025

reinterpolate data
"""
import os, sys
import pandas as pd
import geopandas as gpd
import xarray as xr
import rioxarray
from pathlib import Path
from joblib import Parallel, delayed

def main_large_domain(points_path, coast, res, where_to_save, mask_shp,
        domain=[105.000,117.000,-9.000,-5.000]):
    """ 
    I am trying to reinterpolate the whole domain in one go,
    this to avoid discontinuity between tiles border
    """
    where_to_save = Path(where_to_save)
    where_to_save.mkdir(exist_ok = True)

    #create temporary concate points file
    pts_fin = os.path.join(points_path, 'clip*.pts')
    tmp_pts_fout = os.path.join(where_to_save, f'tmp_dem.pts')
    cmd = f'cat {pts_fin} > {tmp_pts_fout}'
    os.system(cmd)

    #start to reinterpolate
    xmin = domain[0]
    xmax = domain[1]
    ymin = domain[2]
    ymax = domain[3]
    fout_bm = os.path.join(where_to_save, f'tmp_bm_dem.pts')
    cmd = f'gmt blockmean {tmp_pts_fout} -I{res}/{res} -V -R{xmin}/{xmax}/{ymin}/{ymax} > {fout_bm}'
    os.system(cmd)

    tmp_fout_grd = os.path.join(where_to_save, f'tmp_blend_dem.grd')
    cmd = f'gmt surface {fout_bm} -I{res}/{res} -T0.73 -V -D{coast} -Lu30 -G{tmp_fout_grd} -R{xmin}/{xmax}/{ymin}/{ymax}'
    os.system(cmd)

    #clip by another polygon mask and convert to GTiff
    mask = gpd.read_file(mask_shp)
    fout_tif = os.path.join(where_to_save, f'blend_dem.tif')
    dem = xr.open_dataset(tmp_fout_grd)
    dem = dem.rio.write_crs('epsg:4326')
    dem = dem.rio.clip(mask.geometry, mask.crs)
    dem.rio.to_raster(fout_tif, driver='GTiff', compress='LZW')

    fout_grd = os.path.join(where_to_save, f'blend_dem.grd')
    dem.to_netcdf(fout_grd, mode='w', format='NETCDF4_CLASSIC', engine='netcdf4')

    """
    - load interpolated grd file
    - load mask shp to clip, to make sure 'nice'
    - clip grd file with the mask
    - write to raster GTiff
    """
    #cleaning up
    cmd = f'rm {tmp_pts_fout} {fout_bm} {tmp_fout_grd}'
    os.system(cmd)

def main(topo, bathy, coast, res, where_to_save):
    topo = Path(topo)
    bathy = Path(bathy)
    fname = topo.name.split('_')[-1][:-4]
    print(f'blending tile {fname}')

    where_to_save = Path(where_to_save)
    where_to_save.mkdir(exist_ok = True)

    #concate pts from topo and bathy
    fout = os.path.join(where_to_save, f'tmp_pts__{fname}.pts')
    cmd = f'cat {topo} {bathy} > {fout}'
    os.system(cmd)

    #collect extent
    df = pd.read_table(fout, header=None, sep=' ')
    xmin = df[0].min()
    xmax = df[0].max()
    ymin = df[1].min()
    ymax = df[1].max()

    #start to reinterpolate
    fout_bm = os.path.join(where_to_save, f'tmp_bm__{fname}.pts')
    cmd = f'gmt blockmean {fout} -I{res}/{res} -V -R{xmin+0.0001}/{xmax-0.0001}/{ymin+0.0001}/{ymax-0.0001} > {fout_bm}'
    os.system(cmd)

    fout_grd = os.path.join(where_to_save, f'blend_dem__{fname}.grd')
    cmd = f'gmt surface {fout_bm} -I{res}/{res} -T0.73 -V -D{coast} -Lu30 -G{fout_grd} -R{xmin+0.0001}/{xmax-0.0001}/{ymin+0.0001}/{ymax-0.0001}'
    os.system(cmd)

    #cleaning up
    cmd = f'rm {fout} {fout_bm}'
    os.system(cmd)

    return

if __name__ == "__main__":
    ### tiles from DeltaDTM
    tiles = pd.read_csv('DeltaDTM_tiles2correct.csv')
    
    ### configuration for blended DEM
    # resolution, '1s' = ~30m, '3s' = ~90m, '9s' = ~270m
    resolution = '1s'
    # dem extent
    domain = [105.000, 117.000, -9.000, -5.000]
    # masking shp to clean up the DEM file
    mask_shp = '/home/ignatius.pranantyo/DATA/working_deltadtm/vectors/Mask_DeltaDTM_tiles__Sumatra-Jawa-Bali-Lombok-Sumbawa.shp'

    ### points to be used
    # coastline breaklines
    coast_file = '/home/ignatius.pranantyo/DATA/working_deltadtm/vectors/gmt__softbreaklines_coastline__20250424.pts'

    ## correction configuration
    correction = -0.0 # -0.5, -1.0, -1.5
    pts_path = f'/home/ignatius.pranantyo/DATA/working_deltadtm/points__batnas_deltadem_correction_by{correction}m'
    where_to_save = Path(f'/home/ignatius.pranantyo/DATA/working_deltadtm/blended__batnas_deltadtm__correctedby{correction}m__res{resolution}')
    where_to_save.mkdir(exist_ok = True)
    Path(f'/home/ignatius.pranantyo/DATA/working_deltadtm/blended__batnas_deltadtm__correctedby-1.0m__res{resolution}')
    where_to_save.mkdir(exist_ok = True)

    main_large_domain(pts_path, coast_file, resolution, where_to_save, mask_shp, domain)    


    sys.exit()
    ####
    ####
    ####
    ####
    #pts_path = '/home/ignatius.pranantyo/DATA/working_deltadtm/points__batnas_deltadem_correction_by-1.0m'
    #topo_file = os.path.join(pts_path, 'clip-adjust__DeltaDTM_v1_1_S06E105.pts')
    #bathy_file = os.path.join(pts_path, 'clip__BATNAS_v1-5-0_at_DeltaDTM_v1_1_S06E105.pts')
    #coast_file = os.path.join(pts_path, 'softbreaklines_garpan_20250408.pts')
    

    ### coastline definition used as softbreaklines in surface GMT
    coast_file = '/home/ignatius.pranantyo/DATA/working_deltadtm/vectors/gmt__softbreaklines_coastline__20250424.pts'
    
    ### resolution reinterpolate dem, 1s~30m, 3s~90m
    #resolution = '3s'

    ### path to save
    where_to_save = Path(f'/home/ignatius.pranantyo/DATA/working_deltadtm/blended__batnas_deltadtm__correctedby-1.0m__res{resolution}')
    where_to_save.mkdir(exist_ok = True)
    
    #domain = [105.000, 117.000, -9.000, -5.000]
    #mask_shp = '/home/ignatius.pranantyo/DATA/working_deltadtm/vectors/Mask_DeltaDTM_tiles__Sumatra-Jawa-Bali-Lombok-Sumbawa.shp'
    main_large_domain(pts_path, coast_file, resolution, where_to_save, mask_shp, domain)

    ### correction 0.0 m
    ### resolution 3s
    #pts_path = '/home/ignatius.pranantyo/DATA/working_deltadtm/points__batnas_deltadem_correction_by-0.0m'
    #where_to_save = Path(f'/home/ignatius.pranantyo/DATA/working_deltadtm/blended__batnas_deltadtm__correctedby-0.0m__res{resolution}')


    sys.exit()
    topos = []
    bathys = []
    for ii in tiles.index:
        tile = tiles['tile'][ii]
        tname = tile[:-4]
        topos.append(os.path.join(pts_path, f'clip-adjust__{tname}.pts'))
        bathys.append(os.path.join(pts_path, f'clip__BATNAS_v1-5-0_at_{tname}.pts'))

    Parallel(n_jobs = 8)(delayed(main)(topos[ii], bathys[ii], coast_file, resolution, where_to_save) for ii in tiles.index) 
    
    
    sys.exit()
    
    #range(2))

    #main(topo_file, bathy_file, coast_file, resolution, pts_path)
