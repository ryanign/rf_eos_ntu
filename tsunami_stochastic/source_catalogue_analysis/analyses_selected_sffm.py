"""
Ryan Pranantyo
EOS, 28 May 2025

After we obtained SFFM_filename and SFFM_src_id to use for the catalogue (from select_events_catalogue2use.py),
I would like to spatially check whether we have enough scenarios for a specific Mw bin.

What to check:
    - Whether SFFMs have filled up all the grid
    - Distribution of mean and max accross the grid

If visually we see some gap, potentiall we need to enlarge the catalogue
"""
import os, sys
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cmcrameri.cm as cm
import argparse
from pathlib import Path

### FUNCTIONS ###
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


### MAIN CORE ###
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalogue_f", type = str,
                        default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/input_files__SouthernJava/earthquake_catalogue__region/20250602__cat_6.5-8.7_100k_5samples.dat__EVENT_LIST__Mw6.95+.csv")
    parser.add_argument("--Mw", type = float,
                        default = 7.0)
    parser.add_argument("--grid_source", type = str,
                        default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/OUTPUTS__Slab2__Jawa/unit_source_grid/SLAB2__Jawa.shp")
    parser.add_argument("--plot_SFFMs", type = bool,
                        default = False,
                        help = "plotting individual SFFM up to the first 70 models")
    args = parser.parse_args()


    ### let's go!

    catalogue_f = Path(args.catalogue_f)
    cat_fname = catalogue_f.name[:-4]
    figures_path = Path(os.path.join(catalogue_f.parent, "figures"))
    figures_path.mkdir(exist_ok = True)
    
    ### load the catalogue
    df = pd.read_csv(catalogue_f)
    df = df[df["target_Mw"] == args.Mw]

    if len(df) == 0:
        print(f'you do not have event of Mw {args.Mw} in this catalogue!')
        sys.exit()

    ### load the grid source
    grid_gdf = gpd.read_file(args.grid_source)
    # make sure the grid in correct order following strike
    grid_gdf = grid_gdf.sort_values(by = ['alngst_', 'dwndp_n']).reset_index()
    grid_gdf['unit_source_index'] = np.arange(len(grid_gdf)).astype(int) + 1
    grid_xmin = grid_gdf.bounds.minx.min() - 1
    grid_xmax = grid_gdf.bounds.maxx.max() + 1
    grid_ymin = grid_gdf.bounds.miny.min() - 1
    grid_ymax = grid_gdf.bounds.maxy.max() + 1
    fwidth, fheight = fig_size( grid_xmax - grid_xmin , grid_ymax - grid_ymin, 10)

    ### start to collect slip from all model available
    data = list()
    info = list()
    colname = list()
    for ii in df.index:
        sffm_df = pd.read_csv(df["SFFM_filename"][ii])
        sffm_slip = sffm_df[df["SFFM_src_id"][ii]]
        slip = sffm_slip.iloc[:-6].astype(float)
        data.append(slip)
        info.append(sffm_slip.iloc[-6:])
        colname.append(f'realisastion_{ii}')
    slip_df = pd.DataFrame(data = np.vstack(data).T,
                           columns = colname)
    info_df = pd.DataFrame(data = np.vstack(info).T,
                           columns = colname)

    
    stat2checks = ['mean', 'max']

    for ii, stat in enumerate(stat2checks):
        print(ii, stat)
        title = f"Mw{args.Mw} -- Num SFFM = {len(df.index)} -- {stat}"

        if stat == 'mean':
            vals = slip_df.mean(axis=1)
        elif stat == 'max':
            vals = slip_df.max(axis=1)
        else:
            print(f'wrong variable')
            sys.exit()

        grid_gdf['values'] = np.where(vals == 0, np.nan, vals)

        fig = plt.figure(figsize = (fwidth, fheight), constrained_layout = True)
        ax = fig.add_subplot(1,1,1, projection = ccrs.PlateCarree())
        ax.set_title(title)
        ax.coastlines()

        grid_gdf.plot(
                ax = ax,
                column = 'values',
                legend = True,
                legend_kwds = {
                    'label' : f'{stat} - slip, m',
                    'orientation' : 'vertical',
                    'shrink' : 0.8,
                    'pad' : 0.
                    },
                cmap = 'inferno_r',
                )

        grid_gdf.plot(
                ax = ax,
                fc = 'none',
                ec = 'gray',
                linewidth = 0.1
                )

        for jj in range(len(colname)):
            ax.scatter(
                    float(info_df[colname[jj]][0]), float(info_df[colname[jj]][1]), 
                    marker='*',
                    c = 'green',
                    s = 50
                    )

        ax.set_xlim(grid_xmin, grid_xmax)
        ax.set_ylim(grid_ymin, grid_ymax)


        fout = os.path.join(figures_path, f'{cat_fname}__stat_{stat}__Mw{args.Mw}__N-SFFM_{len(df.index)}.png')
        fig.savefig(fout)
        plt.close()


    ### plot individual SFFMs
    if args.plot_SFFMs:
        fig = plt.figure(figsize = (50, 30), constrained_layout = True)
        fig.suptitle(f'Mw{args.Mw}')
        for ii, var in enumerate(colname[:70]):
            print(var)
            ax = fig.add_subplot(7, 10, ii+1, projection = ccrs.PlateCarree())
            ax.coastlines()
            vals = slip_df[var]
            grid_gdf['values'] = np.where(vals == 0, np.nan, vals)
            grid_gdf.plot(
                    ax = ax,
                    column = 'values',
                    legend = True,
                    legend_kwds = {
                        'label' : 'slip, m',
                        'orientation' : 'vertical',
                        'shrink' : 0.2,
                        'pad' : 0.,
                        },
                    cmap = 'inferno_r',
                    vmin = 0,
                    vmax = slip_df.max().max(),
                    )
            ax.scatter(
                    float(info_df[var][0]),
                    float(info_df[var][1]),
                    marker = '*',
                    c = 'green',
                    s = 50,
                    )

            ax.set_xlim(grid_xmin, grid_xmax)
            ax.set_ylim(grid_ymin, grid_ymax)

        fout = os.path.join(figures_path, f'individual__Mw{args.Mw}__N-SFFM_{len(df.index)}.png')
        fig.savefig(fout)
        plt.close(fig)






    




