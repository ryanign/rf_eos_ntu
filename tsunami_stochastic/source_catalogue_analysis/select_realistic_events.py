"""
Ryan Pranantyo
EOS, May 2025

-26 May 2025-
Because we generated too many SFFM, the previous step was not that efficient.
Therefore, here, I am trying to revise the process.

Once the SFFM have been generated, one file is for one target Mw and target epicentre with N realisation model
We need to check whether they are realistic or not first, then proceed to the next steps.
"""
import os, sys
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import argparse
from pathlib import Path
#from joblib import Parallel, delayed
#from joblib.externals.loky import get_reusable_executor
from shapely.geometry import Point

def fig_extent(nc):
    xmin = nc.lon.min().values
    xmax = nc.lon.max().values
    ymin = nc.lat.min().values
    ymax = nc.lat.max().values
    dlon = abs(xmax - xmin)
    dlat = abs(ymax - ymin)
    return xmin, xmax, ymin, ymax, dlon, dlat

def fig_size(dlon, dlat, fwidth = 10.):
    """ if dlat < dlon, figure height should be smaller than fwidth """
    ratio = dlon / dlat
    if ratio <= 1.1 and ratio >= 0.9:
        fheight = fwidth - 0.1
    elif ratio > 1.1:
        fheight = (fwidth / ratio) - 0.035
    else:
        fheight = (fwidth / ratio) + 0.035
    return fwidth, fheight

def check_sffm_realistic(grid_gdf): # , targetMw):
    """ checking if slip distribution is realistic for a given Mw """
    rupt_gdf = grid_gdf.dropna()
    slips = rupt_gdf['slip']
    #print(slips.max(), slips.mean(), slips.std())

    ### check if there is a slip above 40 m, then reject
    if slips.max() >= 40:
        print(f"  unrealistic model ... there is a slip up to {slips.max():.2f} m")
        return False

    ### histogram of slips distribution within 1m bin
    slip_bins = np.arange(-0.5, np.ceil(slips.max()), 1)
    hist , _ = np.histogram(slips, bins = slip_bins)

    ### check if the largest bin has at least 2 patches for large events with sli above 15m
    if slips.max() >= 15:
        if hist[-1] == 1:
            print(slips.max(), hist)
            print(f"  might be still a realistic model, but there is only one value within the highest slip bin")
            return False

    #print("SLIP MODEL IS REALISTIC")
    return True

def check_epi_inside_rupture_area(grid_gdf, target_lon, target_lat):
    P = Point(target_lon, target_lat)

    ### to check wheter there is a slip value at the epicentre
    ### if not, will not be used
    epi_slip = grid_gdf.intersects(P)
    epi_slip_idx = epi_slip[epi_slip == True].index[0]
    epi_nan = np.isnan(grid_gdf['slip'][epi_slip_idx])
    if epi_nan == True:
        epi_is_zero = True
    else:
        epi_is_zero = False

    ### to check wheter epicentre is inside convex hull of the slip model
    ### if it outside, will not be used
    rupt_gdf = grid_gdf.dropna()
    rupt_gdf = rupt_gdf.dissolve()
    convex_hull = rupt_gdf.convex_hull
    rupt_gdf = gpd.GeoDataFrame(geometry = convex_hull, crs = grid_gdf.crs)

    inside = rupt_gdf.intersects(P)

    if inside[0] == False:
        epi_is_outside = True
    else:
        epi_is_outside = False

    if epi_is_zero == True and epi_is_outside == True:
        #print(f" epicentre is outside slip model >> will not be used ...")
        return False
    elif epi_is_zero == False and epi_is_outside == False:
        #print(f" epicentre is inside slip model >> might be used ...")
        return True
    else:
        return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--SFFM_file", type = str,
                        default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/OUTPUTS__Slab2__Jawa/SourceCombinations__20250523/stochastic_slips__SLAB2__Jawa/stochastic_sources__Mw_8.600000__Lon_109.300140__Lat_-9.782850__table.csv",
                        help = "file of SFFM generated")
    parser.add_argument("--sourcename", type = str,
                        default = "SLAB2__Jawa",
                        help = "Name of the sourcename")
    parser.add_argument("--grid_source", type = str,
                        default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/OUTPUTS__Slab2__Jawa/unit_source_grid/SLAB2__Jawa.shp",
                        help = "a polygon shapefile generated when generating unit_source using RPTHA")
    #parser.add_argument("--ncpus", type = int,
    #                    default = 8,
    #                    help = "number of cpus to use for plotting")
    parser.add_argument("--SFFM_plot", type = bool,
                        default = False,
                        help = "plotting individual SFFM")
    args = parser.parse_args()


    ### main program ###
    sffm_f = Path(args.SFFM_file)

    #################################
    ###                           ###
    ### FROM BELOW, THE MAIN CORE ###
    ###                           ###
    #################################

    ### create a new folder to save realistic SFFM tables
    sffm_path_ready2use = Path(os.path.join(sffm_f.parent, "..", f"SFFM_realistic__{args.sourcename}"))
    sffm_path_ready2use.mkdir(exist_ok = True)
    sffm_path_figs = Path(os.path.join(sffm_path_ready2use, "figures"))
    sffm_path_figs.mkdir(exist_ok = True)

    ### unit source filename
    unit_f = os.path.join(sffm_f.parent, f"unit_source_grid_raster_filename_index__{args.sourcename}.csv")
    unit_df = pd.read_csv(unit_f)
    unit_df = unit_df.rename(
            columns = {
                'Unnamed: 0' : 'unit_source_index',
                'x'          : 'unit_source_filename'}
            )

    ### read SFFM grid polygon used
    grid_source = Path(args.grid_source)
    grid_gdf = gpd.read_file(grid_source)
    #re-order index of the SFFM to match SFFM tables
    grid_gdf = grid_gdf.sort_values(by=['alngst_', 'dwndp_n']).reset_index()
    #to make same column name
    grid_gdf['unit_source_index'] = (np.arange(len(grid_gdf))+1).astype(int)
    grid_gdf = grid_gdf.drop(columns=['index'])
    grid_gdf.set_index = grid_gdf['unit_source_index']
    grid_fout = os.path.join(sffm_path_ready2use, f"final__{grid_source.name}")
    if os.path.exists(grid_fout) == False:
        grid_gdf.to_file(grid_fout) 
    
    ### configuration for figsize, incase we plot the SFFM
    grid_xmin = grid_gdf.bounds.minx.min() - 1
    grid_xmax = grid_gdf.bounds.maxx.max() + 1
    grid_ymin = grid_gdf.bounds.miny.min() - 1
    grid_ymax = grid_gdf.bounds.maxy.max() + 1
    fwidth, fheight = fig_size( grid_xmax - grid_xmin , grid_ymax - grid_ymin, 10)

    ### building a DataFrame of SFFM
    sffm_df = pd.read_csv(sffm_f)
    unit_source_index = (np.arange(len(grid_gdf))+1).astype(int)
    slip_df = pd.DataFrame(data = unit_source_index,
                           columns = ['unit_source_index'])
    info_df = pd.DataFrame(data = ['target_lon', 
                                   'target_lat', 
                                   'Mw', 
                                   'physical_corner_wavenumber_x', 
                                   'physical_corner_wavenumber_y', 
                                   'sourcename'],
                           columns = ['unit_source_index'])

    #df['unit_source_index'] = np.hstack([unit_source_index, additional])
    for src in sffm_df.index:
        flt_index_str = sffm_df.event_index_string[src]
        flt_slip_str  = sffm_df.event_slip_string[src]
        flt_index = list(map(int, flt_index_str.split('-')[:-1]))
        flt_slip  = list(map(float, flt_slip_str.split('_')[:-1]))
        check_tmp_df = pd.DataFrame(data = {f'unit_source_slip__{src}' : flt_slip,
                                      'unit_source_index' : flt_index},
                              index = flt_index)

        grid_gdf['slip'] = check_tmp_df[f'unit_source_slip__{src}']

        target_lon = sffm_df['target_lon'][src]
        target_lat = sffm_df['target_lat'][src]
        target_Mw  = sffm_df['Mw'][src]
        
        inside = check_epi_inside_rupture_area(grid_gdf, target_lon, target_lat)
        realistic = check_sffm_realistic(grid_gdf)

        if inside == True and realistic == True:
            tmp_df = pd.DataFrame(data = {f'unit_source_slip__{src}' : flt_slip,
                                          f'unit_source_index' : flt_index},
                                  index = flt_index)

            info1 = [sffm_df['target_lon'][src],
                     sffm_df['target_lat'][src],
                     sffm_df['Mw'][src],
                     sffm_df['physical_corner_wavenumber_x'][src],
                     sffm_df['physical_corner_wavenumber_y'][src],
                     sffm_df['sourcename'][src]]
            info2 = ['target_lon',
                     'target_lat',
                     'Mw',
                     'physical_corner_wavenumber_x',
                     'physical_corner_wavenumber_y',
                     'sourcename']
            inftmp_df = pd.DataFrame(data = {f'unit_source_slip__{src}' : info1,
                                             f'unit_source_index' : info2})
            
            slip_df = slip_df.merge(tmp_df, on='unit_source_index', how = 'left')
            info_df = info_df.merge(inftmp_df, on='unit_source_index', how = 'left')

            if args.SFFM_plot == True: 
                fig = plt.figure(figsize = (fwidth, fheight), constrained_layout = True)
                ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
                ax.coastlines()
                ax.set_title(f"Mw = {target_Mw:.1f} -- Lon = {target_lon:.5f} -- Lat = {target_lat:.5f} -- Sample Num = {src}")
                ax.scatter(target_lon, target_lat, marker='*', c='green', zorder=99,
                               s=50)
                grid_gdf.plot(ax = ax,
                            column = 'slip',
                            ec = 'gray',
                            legend = True,
                            legend_kwds = {
                                "label" : "slip, m",
                                "orientation" : "vertical",
                                "shrink" : 0.8,
                                "pad" : 0},
                            cmap = 'inferno_r')
                grid_gdf.plot(
                            ax = ax,
                            fc = 'none', ec='gray',
                            linewidth = 0.1)
                ax.set_xlim(grid_xmin, grid_xmax)
                ax.set_ylim(grid_ymin, grid_ymax)
    
                fout = os.path.join(sffm_path_figs,
                        f'SFFM__{src:08d}__Mw_{target_Mw:.6f}__Lon_{target_lon:.6f}__Lat_{target_lat:.6f}.png')
                fig.savefig(fout)
                plt.close()

            ### constructing a DF
            #slip_df = pd.merge([slip_df, tmp_df], on='unit_source_index', how = 'left')
            #info_df = pd.merge([info_df, inftmp_df], on='unit_source_index', how = 'left')

            #trgt_idx = df.index[df['unit_source_index'] == 'target_lon'].to_list()[0]
            #tmp_df.loc[trgt_idx] = [sffm_df['target_lon'][src], 'target_lon']

            #trgt_idx = df.index[df['unit_source_index'] == 'target_lat'].to_list()[0]
            #tmp_df.loc[trgt_idx] = [sffm_df['target_lat'][src], 'target_lat']

            #trgt_idx = df.index[df['unit_source_index'] == 'Mw'].to_list()[0]
            #tmp_df.loc[trgt_idx] = [sffm_df['Mw'][src], 'Mw']

            #trgt_idx = df.index[df['unit_source_index'] == 'physical_corner_wavenumber_x'].to_list()[0]
            #tmp_df.loc[trgt_idx] = [sffm_df['physical_corner_wavenumber_x'][src], 'physical_corner_wavenumber_x']

            #trgt_idx = df.index[df['unit_source_index'] == 'physical_corner_wavenumber_y'].to_list()[0]
            #tmp_df.loc[trgt_idx] = [sffm_df['physical_corner_wavenumber_y'][src], 'physical_corner_wavenumber_y']
            #
            #trgt_idx = df.index[df['unit_source_index'] == 'sourcename'].to_list()[0]
            #tmp_df.loc[trgt_idx] = [sffm_df['sourcename'][src], 'sourcename']
            #df = pd.merge(df, tmp_df, on='unit_source_index', how='left')
    
    slip_df = slip_df.fillna(0)                                     ### NaN slip convert to 0.0
    df = pd.concat([slip_df, unit_df], axis=1)                      ### concate with unitsourcefilename
    df = df.loc[:, ~df.columns.duplicated()]                        ### remove duplicates columns
    df = pd.concat([df, info_df], axis=0)                           ### concate info_df
    print(df)

    ### saving realistic SFFM into a file 
    fout = os.path.join(sffm_path_ready2use, f'realistic__{sffm_f.name}')
    df.to_csv(fout, index = False)




    sys.exit()
