"""
Ryan Pranantyo
EOS, April 2025

script to launch SFINCS simulations for one event:
    1. Extract max elevation at open bc points, filter out if the values less than a defined threshold (e.g. 0.5m)
    2. Remap the points belong to which tile
    3. Extract offshore timeseries per-tile
    4. Build event for SFINCS and execute
"""
import os, sys
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pathlib import Path
from shapely.geometry import Point
from joblib import Parallel, delayed

def mapping_max_values(script, jagurs_nc, bc_open_f, scenario, where_to_save, tiles, threshold=1.0, scale_ratio=10000, plot=True):
    """
    main target of this function is to filter which tiles to be simulated for sfincs simulation
    1. extract jagurs_nc max elevation at bc_open_f
    2. remove points if max elevation below threshold (df_max)
    3. find the point (df_max) inside sfincs domain (tiles/df_tile)

    scale_ratio = JAGURS output file is saved in int32 format by multiplying the values by this scale ratio,
        by default, I use 10000. 
    Hence, the threshold to define which tiles to be activated will be multiplied by this scale_ratio.

    df_tile['tile_name'] = list of tiles to be simulated for coastal inundation, in this case using sfincs.
    """
    where_to_save = Path(os.path.join(where_to_save, scenario))
    where_to_save.mkdir(exist_ok = True)
    where_output_is = os.path.join(where_to_save, f'extracted__{jagurs_nc.name[:-3]}')

    ### extract maximum values at open bc points
    cmd = f'python -W ignore {script} --jagurs_nc {jagurs_nc} --list_of_coordinates {bc_open_f} --where_to_save {where_to_save} --var2extract 2'
    os.system(cmd)
    max_value = os.path.join(where_output_is, f'maxvalues_at_{bc_open_f.name}')
    df = pd.read_csv(max_value)
    df = df.T
    df = df.rename(columns = {0 : 'max_value'})
    df['NAME'] = df.index
    df = df.reset_index()
    bc_df = pd.read_csv(bc_open_f)
    df = df.merge(bc_df)

    ### filterout max values under threshold
    df_max = df[df['max_value'] >= threshold * scale_ratio]

    df_tile = tiles
    ### check which tile to run
    df_tile['RUN'] = 0
    for ii in df_max.index:
        point = Point(df_max['LON'][ii], df_max['LAT'][ii])
        con = df_tile.contains(point)
        for jj in df_tile.index:
            if con[jj] == True:
                df_tile.loc[jj, 'RUN'] += 1
    ### if more than df_tile['RUN'] is more than 10 nodes, then we run
    df_tile = df_tile[df_tile['RUN'] >  (1/3 * df_tile['NUM_PTS'])]

    ### save tiles to run
    fout = os.path.join(where_to_save, f'tiles2run__{scenario}.csv')
    tiles = df_tile['tile_name']
    tiles.to_csv(fout, index=False, header=None)

    print(df_tile)
    
    if plot == True:
        bc_df = pd.read_csv(bc_open_f)
        fig = plt.figure(constrained_layout = True)
        ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
        ax.coastlines()
        ax.scatter(df.LON, df.LAT, marker='o', fc='none', ec='gray')
        cb = ax.scatter(df_max.LON, df_max.LAT, marker='o', c=df_max['max_value'] / scale_ratio, 
                ec='none', vmin=threshold)
        df_tile.boundary.plot(ax=ax, linewidth=0.2)
        fig.colorbar(cb, pad=0, orientation='vertical', shrink=0.5, extend='max', 
                label='max values at SFINCS openbc, m')
        fout = os.path.join(where_output_is, 'max_values.png')
        fig.savefig(fout)
    return df, df_max, df_tile

def check_points_inside_tile(bc_pts, tiles):
    """
    df_tile = to count number of open bc points inside one domain
    bc_df = to show point X belong to tile Y
    """
    bc_df = pd.read_csv(bc_pts)
    df_tile = gpd.read_file(tiles)
    df_tile['NUM_PTS'] = 0
    for ii in bc_df.index:
        point = Point(bc_df['LON'][ii], bc_df['LAT'][ii])
        con = df_tile.contains(point)
        for jj in df_tile.index:
            if con[jj] == True:
                df_tile.loc[jj, 'NUM_PTS'] += 1
                bc_df.loc[ii, df_tile['tile_name'][jj]] = True

    #re-save bc_pts
    fname = f"{bc_pts.name[:-4]}__belong2whichdomain.csv"
    fout = os.path.join(bc_pts.parent, fname)
    bc_df.to_csv(fout, index=False)

    #re-split bc_pts to per-tile
    for ii in df_tile.index:
        tile_name = df_tile['tile_name'][ii]
        bc_tmp = bc_df[bc_df[tile_name] == True]
        fout = os.path.join(bc_pts.parent, f"{bc_pts.name[:-4]}__{tile_name}.csv")
        bc_tmp.to_csv(fout, index=False)
    
    return df_tile, bc_df


def prepare_openbc_pertile(script, jagurs_nc, bc_open_f, scenario, where_to_save, tile):
    """
    more and less similar with mapping_max_values() but this is for elevation time-series
    and will be done per-sfincs-tile

    bc_file = is input file for SFINCS
    """
    where_to_save = Path(os.path.join(where_to_save, scenario, f'sfincs-open-bc__at__{tile}'))
    where_to_save.mkdir(exist_ok = True)
    
    ### main part to extract elev time series
    cmd = f'python -W ignore {script} --jagurs_nc {jagurs_nc} --list_of_coordinates {bc_open_f} --where_to_save {where_to_save} --var2extract 1'
    print(cmd)
    os.system(cmd)

    bc_file = os.path.join(where_to_save, f'extracted__{jagurs_nc.name[:-3]}', f'timeseries__at_{bc_open_f.name}')

    ### can add a short function to clip the timeseries just around the peak.

    return cmd, bc_file

def check_individual_event(
        ii, 
        events_df, 
        jagurs_simulations_path,
        script2extract,
        points_openbc_line,
        main_where_to_save,
        sfincs_tiles,
        jagurs_nc_out='SD01.nc',
        threshold=1.0,
        scale_ratio=10000):
    """ 
    main part 
    """
    
    event_name = events_df['EVENT_NAME'][ii]
    event_id = events_df['EVENT_ID'][ii]
    scenario = f'{event_id:08d}__{event_name}'
    print(event_name, event_id)

    jagurs_nc = Path(os.path.join(jagurs_simulations_path, scenario, jagurs_nc_out))

    ### extract sfincs tiles to simulate
    df, df_max, df_tile = mapping_max_values(script2extract, jagurs_nc, points_openbc_line, scenario, main_where_to_save, sfincs_tiles, threshold, scale_ratio, plot=True)
    
    ### to update events_df by adding path where the openbc timeseries would be saved
    openbc_timeseries_path = os.path.join(main_where_to_save, scenario)
    tiles2run = os.path.join(main_where_to_save, scenario, f'tiles2run__{scenario}.csv')

    ### prepare open bc timeseries per-tile
    for jj in df_tile.index:
        tile_name = df_tile['tile_name'][jj]
        bc_open_tile_f = Path(os.path.join(points_openbc_line.parent, f"{points_openbc_line.name[:-4]}__{tile_name}.csv"))

        cmd, sfincs_open_bc_f = prepare_openbc_pertile(
                script2extract, 
                jagurs_nc, 
                bc_open_tile_f, 
                scenario, 
                main_where_to_save, 
                tile_name
                )

    return openbc_timeseries_path, tiles2run


if __name__ == "__main__":

    """
    this section should be fixed, update as needed
    """
    # script to extract elevation timeseries
    script2extract = "/home/ignatius.pranantyo/apps/rf_eos_ntu/tsunami_stochastic/main_extract_elev_timeseries.py"

    # base file where the open bc points are located in SFINCS
    points_openbc_line = Path("/home/ignatius.pranantyo/Tsunamis/AOGS2025__StressCaseScenarios/SFINCS_config/vectors/domains_tile/points_open_bc_line_jawa.csv")

    # where to save open bc timeseries for SFINCS
    main_where_to_save = Path("/home/ignatius.pranantyo/Tsunamis/AOGS2025__StressCaseScenarios/SFINCS_onshore_simulations/ELEV_AT_OPENBC_JAWA")
    main_where_to_save.mkdir(exist_ok = True)

    # SFINCS domain tiles in a polygon, as a base to determine which tiles to be simulated by SFINCS
    ##domain_tiles = "/home/ignatius.pranantyo/Tsunamis/AOGS2025__StressCaseScenarios/SFINCS_config/vectors/domain_tiles.shp"
    #JUST FOR AOGS
    domain_tiles = "/home/ignatius.pranantyo/Tsunamis/AOGS2025__StressCaseScenarios/SFINCS_config/vectors/domain_tiles__AOGS2025.shp"
    # JAGURS output file to use
    jagurs_nc_out = 'SD01.nc'
    # multiplication factor coef used to convert JAGURS output from float to int32
    scale_ratio = 10000

    # offshore minimum elevation to be considered for SFINCS simulation
    threshold = 0.05 # in metre

    sfincs_tiles, bc_df = check_points_inside_tile(points_openbc_line, domain_tiles)

    #sys.exit()


    """
    section below is the way to launch timeseries extrcation and determine which tiles to be activated for SFINCS simulations
    """
    ### section below follow scenario to check! ###
    #MODEL_CONFIGURATION_NAME = "TESTING__20250502"
    #events_catalogue = Path('../EventsCatalogue__Jawa.csv')

    MODEL_CONFIGURATION_NAME = "20250513"
    events_catalogue = Path('/home/ignatius.pranantyo/Tsunamis/AOGS2025__StressCaseScenarios/EventsCatalogue__Sample__AOGS2025.csv')


    events_df = pd.read_csv(events_catalogue)
    events_df['openbc_timeseries_path'] = None
    events_df['tiles2run'] = None
    jagurs_simulations_path = '/home/ignatius.pranantyo/Tsunamis/AOGS2025__StressCaseScenarios/JAGURS_offshore_simulations/'
    
    #ii = 0
    results = Parallel(n_jobs = 6)(
            delayed(check_individual_event)(
                ii,
                events_df,
                jagurs_simulations_path,
                script2extract,
                points_openbc_line,
                main_where_to_save,
                sfincs_tiles,
                jagurs_nc_out,
                threshold,
                scale_ratio) 
            for ii in events_df.index)
    
    open_bc_timeseries_path, tiles2run = zip(*results)
    
    events_df['openbc_timeseries_path'] = open_bc_timeseries_path
    events_df['tiles2run'] = tiles2run

    ### re-save events_catalogue dataframe
    fout = os.path.join(events_catalogue.parent, events_catalogue.name[:-4] + '__' + MODEL_CONFIGURATION_NAME + '.csv')
    events_df.to_csv(fout, index = False)

    print(events_df)
    
    sys.exit()
    """
    DONE
    DONE
    DONE
    DONE
    """




