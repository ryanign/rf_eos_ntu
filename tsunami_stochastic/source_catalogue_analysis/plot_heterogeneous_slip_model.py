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
from joblib import Parallel, delayed
from joblib.externals.loky import get_reusable_executor
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

    ### reading ...
    SFFM_f = Path(args.SFFM_model)
    where_to_save = Path(os.path.join(SFFM_f.parent, "figures"))
    where_to_save.mkdir(exist_ok = True)

    df = pd.read_csv(SFFM_f, low_memory = False)
    if args.SPECIFIC_ID is not None:
        """ just to check specific event """
        cols2take = list()
        cols2take.append(df.columns[0 ])
        cols2take.append(df.columns[1 ])
        cols2take.append(df.columns[-1])
        cols2take.append(f"unit_source_slip__{args.SPECIFIC_ID}")
        df = df[cols2take]

    sffm_df = df.iloc[:-6]
    info_df = df.iloc[-6:]
    
    print(df)
    if args.SPECIFIC_ID is not None:
        eventid2use = plot_individual_sffm(args.SPECIFIC_ID, grid_gdf, sffm_df, info_df, where_to_save)
    else:
        if args.NSample == 0:
            NSamples = len(df.columns) - 4
        else:
            NSamples = args.NSample
        eventid2use = Parallel(n_jobs = args.ncpus)(
                                delayed(
                                plot_individual_sffm)(
                                    ii, grid_gdf, sffm_df, info_df, where_to_save) 
                                for ii in range(NSamples)
                                )

        ### resaving SFFM model that is realistc
        cols2take = list()
        cols2take.append(df.columns[0 ])
        cols2take.append(df.columns[-1])
        for idx in eventid2use:
            if idx != -99999999:
                cols2take.append(f"unit_source_slip__{idx}")
        filt_df = df[cols2take]
        final_sffm_path = Path(os.path.join(SFFM_f.parent, "sffm_read2use"))
        final_sffm_path.mkdir(exist_ok = True)
        filt_df.to_csv(os.path.join(final_sffm_path, f"filtered__N{NSamples}__{SFFM_f.name}"), index = False)

    print(f"DONE ...")
    print(f"DONE ...")

    return grid_gdf, df, eventid2use

def plot_individual_sffm(ii, grid_gdf, sffm_df, info_df, where_to_save):
    ### figure configuration
    grid_xmin = grid_gdf.bounds.minx.min() - 1
    grid_xmax = grid_gdf.bounds.maxx.max() + 1
    grid_ymin = grid_gdf.bounds.miny.min() - 1
    grid_ymax = grid_gdf.bounds.maxy.max() + 1
    fwidth, fheight = fig_size( grid_xmax - grid_xmin , grid_ymax - grid_ymin, 10)

    ### assign slip to polygon
    grid_gdf['slip'] = (sffm_df[f'unit_source_slip__{ii}']).astype(float)
    target_lon = float(info_df.iloc[-6][f'unit_source_slip__{ii}'])
    target_lat = float(info_df.iloc[-5][f'unit_source_slip__{ii}'])
    target_Mw  = float(info_df.iloc[-4][f'unit_source_slip__{ii}'])

    ### check if target_lon and target_lat inside random slip model
    use_and_plot = check_epi_inside_rupture_area(grid_gdf, target_lon, target_lat)

    ### check if slip distribution is realistic, particularly the maximum slip
    realistic = check_sffm_realistic(grid_gdf) #, target_Mw)

    if use_and_plot == True and realistic == True:    
        ### plot
        fig = plt.figure(figsize = (fwidth, fheight), constrained_layout = True)
        ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
        ax.set_title(f'Mw = {target_Mw}', pad = 0)
        ax.coastlines()
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

        fout = os.path.join(where_to_save, 
                f'SFFM__{ii:08d}__{target_Mw:.6f}__Lon_{target_lon:.6f}__Lat_{target_lat:.6f}.png')
        fig.savefig(fout)
        plt.close()

        return ii
    else:
        return -99999999

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
        print(f"SLIP MODEL WILL NOT BE USED ...")
        use_and_plot = False
    elif epi_is_zero == False and epi_is_outside == False:
        print(f"SLIP MODEL WILL BE USED ...")
        use_and_plot = True
    else:
        use_and_plot = False

    return use_and_plot

def check_sffm_realistic(grid_gdf): # , targetMw):
    """ checking if slip distribution is realistic for a given Mw """
    rupt_gdf = grid_gdf.dropna()
    slips = rupt_gdf['slip']
    print(slips.max(), slips.mean(), slips.std())

    ### check if there is a slip above 40 m, then reject
    if slips.max() >= 40:
        print(f"  unrealistic model ... there is a slip up to {slips.max():.2f} m")
        return False

    ### histogram of slips distribution within 1m bin
    slip_bins = np.arange(-0.5, np.ceil(slips.max()), 1)
    hist , _ = np.histogram(slips, bins = slip_bins)
    
    ### check if the largest bin has at least 2 patches
    if hist[-1] == 1:
        print(f"  might be still a realistic model, but there is only one value within the highest slip bin")
        return False

    print("SLIP MODEL IS REALISTIC")
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--grid_source", type = str,
                        default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/OUTPUTS__Slab2__Jawa/SourceCombinations__20250516/stochastic_slips__SLAB2__Jawa/SFFM_tables__collection/final__SLAB2__Jawa.shp",
                        #default = "/home/ryan/OneDrive_NTU_Projects/Tsunamis_PTHA_inundation/EarthquakeCatalogue/Segmentations/SourcesGrid/SLAB2__Jawa.shp",
                        help = "a polygon shapefile generated when generating unit_source using RPTHA")
    parser.add_argument("--SFFM_model", type = str,
                        default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/OUTPUTS__Slab2__Jawa/SourceCombinations__20250516/stochastic_slips__SLAB2__Jawa/SFFM_tables__collection/stochastic_sources__Mw_8.700000__table.csv",
                        #default = "/home/ryan/IO_for_github/earthquake_catalogue/stochastic_sources__Mw_7.700000__Lon_113.114200__Lat_-10.113100__table.csv",
                        #default = "/home/ryan/IO_for_github/earthquake_catalogue/stochastic_sources__Mw_7.500000__Lon_105.253600__Lat_-8.241500__table.csv",
                        help = "a list of stochastic finite fault mode generated from automate_stochastic_generation.py")
    parser.add_argument("--SFFM_model_uncombined", type = bool,
                        default = False,
                        help = "False: SFFM model need to be 'cleaned up' first; True = is ready to use")
    parser.add_argument("--NSample", type = int,
                        default = 10,
                        help = "number of samples to plot, if 0 = plot all samples")
    parser.add_argument("--SPECIFIC_ID", type = int,
                        default = None, #300,
                        help = "analyses specific event number from the SFFM list")
    parser.add_argument("--ncpus", type = int,
                        default = 8,
                        help = "number of cpus to use for plotting")
    #parser.add_argument("--where_to_save", type = str,
    #                    default = "/home/ryan/IO_for_github/earthquake_catalogue",
    #                    help = "where to save the figure")

    args = parser.parse_args()

    grid_gdf, df, eventid2use = main(args)

    get_reusable_executor().shutdown(wait=True)
    sys.exit()
