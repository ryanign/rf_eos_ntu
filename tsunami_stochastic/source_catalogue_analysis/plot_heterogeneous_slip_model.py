"""
Ryan Pranantyo
EOS, May 2025

grid source fault plane and random slip are based on the RPTHA code
"""
import os, sys
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import argparse
from pathlib import Path

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

def clean_up_sffm(sffm_df, grid_gdf):
    print("cleaning up SFFM model")

    df = pd.DataFrame()
    df['unit_source_index'] = (np.arange(len(grid_gdf))+1).astype(int)
    for src in sffm_df.index:
        flt_index_str = sffm_df.event_index_string[src]
        flt_slip_str  = sffm_df.event_slip_string[src]
        flt_index = list(map(int, flt_index_str.split('-')[:-1]))
        flt_slip  = list(map(float, flt_slip_str.split('_')[:-1]))
        tmp_df = pd.DataFrame(data = {'unit_source_index' : flt_index,
                                      f'unit_source_slip__{src}'  : flt_slip})
        df = pd.merge(df, tmp_df, on='unit_source_index', how='left')
    return df, grid_gdf

def main(args):
    print(args)

    ### load grid shp
    grid_gdf = gpd.read_file(args.grid_source)
    grid_xmin = grid_gdf.bounds.minx.min() - 1
    grid_xmax = grid_gdf.bounds.maxx.max() + 1
    grid_ymin = grid_gdf.bounds.miny.min() - 1
    grid_ymax = grid_gdf.bounds.maxy.max() + 1

    fwidth, fheight = fig_size( grid_xmax - grid_xmin , grid_ymax - grid_ymin, 10)
    
    grid_gdf = grid_gdf.sort_values(by=['alngst_', 'dwndp_n']).reset_index()

    ### cleaning up
    df = pd.read_csv(args.SFFM_model)
    if args.SFFM_model_uncombined == False:
        sffm_df, grid_gdf = clean_up_sffm(df, grid_gdf)
    
    print(df)

    for ii in range(args.NSample):
        grid_gdf['slip'] = sffm_df[f'unit_source_slip__{ii}']

        ###
        fig = plt.figure(figsize = (fwidth, fheight), constrained_layout = True)
        ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
        ax.coastlines()
        ax.scatter(df['target_lon'][ii], df['target_lat'][ii], marker='*', c='green', zorder=99,
                   s=50)
        grid_gdf.plot(ax = ax,
                      column = 'slip',
                      ec = 'gray',
                      legend = True,
                      legend_kwds = {"label" : "slip, m", 
                                     "orientation" : "vertical",
                                     "shrink" : 0.8,
                                     "pad" : 0},
                      cmap = 'inferno_r')
        grid_gdf.plot(ax = ax,
                      fc = 'none', ec='gray',
                      linewidth = 0.1)
        ax.set_xlim(grid_xmin, grid_xmax)
        ax.set_ylim(grid_ymin, grid_ymax)

        fout = os.path.join(args.where_to_save, f'scenario_{ii}.png')
        fig.savefig(fout)
        plt.close()

    return grid_gdf, df





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--grid_source", type = str,
                        default = "/home/ryan/OneDrive_NTU_Projects/Tsunamis_PTHA_inundation/EarthquakeCatalogue/Segmentations/SourcesGrid/SLAB2__Jawa.shp",
                        help = "a polygon shapefile generated when generating unit_source using RPTHA")
    parser.add_argument("--SFFM_model", type = str,
                        default = "/home/ryan/IO_for_github/earthquake_catalogue/stochastic_sources__Mw_7.700000__Lon_113.114200__Lat_-10.113100__table.csv",
                        #default = "/home/ryan/IO_for_github/earthquake_catalogue/stochastic_sources__Mw_7.500000__Lon_105.253600__Lat_-8.241500__table.csv",
                        help = "a list of stochastic finite fault mode generated from automate_stochastic_generation.py")
    parser.add_argument("--SFFM_model_uncombined", type = bool,
                        default = False,
                        help = "False: SFFM model need to be 'cleaned up' first; True = is ready to use")
    parser.add_argument("--NSample", type = int,
                        default = 10,
                        help = "number of samples to ploe")
    parser.add_argument("--where_to_save", type = str,
                        default = "/home/ryan/IO_for_github/earthquake_catalogue",
                        help = "where to save the figure")

    args = parser.parse_args()

    grid_gdf, df = main(args)
