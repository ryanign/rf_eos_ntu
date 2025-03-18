"""
Ryan Pranantyo
EOS, March 2025

main script to extract wave height timeseries or max value from JAGURS output file
inputs: - SD00.nc, SD01.nc, or any netcdf resulted from JAGURS
        - a list of coordinates where to extract
outputs: - wave height elevation timeseries
         - max wave height eleavation
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
import argparse
from pathlib import Path

def haversine(lon1, lat1, lon2, lat2):
    """ 
    calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    https://stackoverflow.com/q/29545704
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lon2 = map(np.radians, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
    c = 2 * np.asin(np.sqrt(a))
    dist = 6367 * c #in km
    return dist

def find_index(lons, lats, lon_target, lat_target):
    """
    find the closest index from data
    """
    dlon = abs(lons - lon_target)
    dlat = abs(lats - lat_target)
    ilon = np.where(dlon == dlon.min())[0][0]
    ilat = np.where(dlat == dlat.min())[0][0]
    return ilon, ilat

def extract(args):
    print(f"GOING TO EXTRACT ...")
    print(args)
    ### configuration ###
    infile = Path(args.jagurs_nc)
    fname = infile.name[:-3]
    main_where_to_save = Path(args.where_to_save)
    main_where_to_save.mkdir(exist_ok = True)

    stations = Path(args.list_of_coordinates)
    stations_name = stations.name[:-4]

    where_to_save = Path(os.path.join(main_where_to_save, f"extracted__{fname}"))
    where_to_save.mkdir(exist_ok = True)

    if args.var2extract == 1:
        fout = os.path.join(where_to_save, f"timeseries__at_{stations_name}.csv")
    elif args.var2extract == 2:
        fout = os.path.join(where_to_save, f"maxvalues_at_{stations_name}.csv")
    else:
        print(f"Wrong choice of variable to extract: STOP")
        sys.exit()
    ###

    ### main ###
    nc = xr.open_dataset(infile)
    station_df = pd.read_csv(stations)

    ### find index
    idx_lon = np.zeros(len(station_df)).astype(int)
    idx_lat = np.zeros(len(station_df)).astype(int)
    """
    STILL NEED TO MAKE SURE TO EXTRACT THE CORRECT INDEX
    """
    bc_df = pd.DataFrame()
    if args.var2extract == 1:
        data = nc.wave_height
        bc_df["time"] = data.time.astype("datetime64[s]")
    else:
        data = nc.max_height

    for idx in station_df.index:
        idx_lon[idx], idx_lat[idx] = find_index(nc.lon.data, nc.lat.data, 
                station_df["LON"][idx], station_df["LAT"][idx])
        #bc_df[f"bc_{station_df['ID'][idx]}"] = data[:, idx_lon[idx], idx_lat[idx]]
        
        ### if want to plot the timeseries extracted
        if args.var2extract == 1:
            bc_df[f"bc_{station_df['ID'][idx]}"] = data[:, idx_lon[idx], idx_lat[idx]]
            if args.plot_extracted:
                plot_timeseries(data, idx_lon[idx], idx_lat[idx], station_df['NAME'][idx], where_to_save)
        else:
            bc_df[f"bc_{station_df['ID'][idx]}"] = data[idx_lon[idx], idx_lat[idx]]

    ### save output
    bc_df.to_csv(fout, index=False)

    plot_map(station_df, where_to_save)

    return nc, station_df

def plot_timeseries(data, ilon, ilat, name, where_to_save):
    """
    plot elevation timeseries
    """
    import matplotlib.dates as mdates
    time = data.time.values
    fig = plt.figure(constrained_layout = True)
    ax = fig.add_subplot(1,1,1)
    ax.set_title(name)
    ax.plot(time, data[:,ilat, ilon]/10000, c='black')
    ax.set_xlabel('time')
    ax.set_ylabel('elevation [m]')
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%m"))
    fig_fout = os.path.join(where_to_save, f"timeseries__{name}.png")
    fig.savefig(fig_fout)
    plt.close()

def plot_map(station_df, where_to_save):
    ### plot map ###
    fig = plt.figure(constrained_layout = True)
    ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.scatter(station_df['LON'], station_df['LAT'], marker='o', c='red')
    ax.set_xlim(station_df['LON'].min() - 5, station_df['LON'].max() + 5)
    ax.set_ylim(station_df['LAT'].min() - 5, station_df['LAT'].max() + 5)
    fig_fout = os.path.join(where_to_save, f"map.png")
    fig.savefig(fig_fout, dpi=150)
    plt.close(fig)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--jagurs_nc", type=str, default = "/scratch/ignatius.pranantyo/tsunami_footprints__2700m_in_integer10000/footprints_stochastic_sources__Mw_7.900000__Lon_105.233600__Lat_-7.923100__table_simplified/footprints__sample_num__24.nc",
            help = "JAGURS nc output file")
    parser.add_argument("--list_of_coordinates", type=str, default = "./input_examples/coordinates_jakarta_sukabumi.csv",
            help = "a list of coordinates where to extract: ID, LON, LAT, (DEPTH)")
    parser.add_argument("--where_to_save", type=str, default = "./output_examples/timeseries",
            help = "where to save the extracted values")
    parser.add_argument("--var2extract", type=int, default = 1,
            help = "1: extract elevation timeseries, 2: extract maximum elevation")
    parser.add_argument("--plot_extracted", type=bool, default = False,
            help = "True or False, plotting the extracted data, for now only for the timeseries")
    args = parser.parse_args()

    nc, station_df = extract(args)
