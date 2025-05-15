"""
Ryan Pranantyo
EOS, March 2025

a basic script just to check JAGURS nc output file
    > plot initial displacement
    > plot max wave height footprint
    > plot footprint timeseries, if wanted
"""
import os, sys
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib as mpl
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

def initial_displacement(infile, vmin, vmax, scale_ratio, dem, XX, YY):
    print(infile, vmin, vmax)
    nc = xr.open_dataset(infile)
    nc = nc["initial_displacement"]
    fname = infile.name[:-3]
    parent_dir = infile.parent
    fig_dir = os.path.join(parent_dir, "figures")

    xmin, xmax, ymin, ymax, dlon, dlat = fig_extent(nc)
    fwidth, fheight = fig_size(dlon, dlat)

    fig = plt.figure(figsize = (fwidth, fheight), constrained_layout = True)
    ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
    ax.set_title(f"Initial displacement {fname} [m]", pad=0.01)
    ax.coastlines()
    disp = ax.imshow(nc.data[0] / scale_ratio, vmin=vmin, vmax=vmax,
            extent=[xmin, xmax, ymin, ymax], cmap="RdBu_r",
            origin="lower")
    fig.colorbar(disp, orientation="vertical", 
            extend="both", pad=0, shrink=0.8)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    if dem is not None and XX is not None and YY is not None:
        ax.contour(XX, YY, dem[::-1], [0], colors='white', zorder=100, linewidths=2)

    fout = os.path.join(fig_dir, f"init_disp__{fname}.png")
    fig.savefig(fout, dpi=300)
    plt.close()

    return nc

def tsunami_max_footprint(infile, scale_ratio, dem, XX, YY):
    print(infile)
    nc = xr.open_dataset(infile)
    nc = nc["max_height"]
    nc = nc.where(nc.data != 0 , np.nan)

    fname = infile.name[:-3]
    parent_dir = infile.parent
    fig_dir = os.path.join(parent_dir, "figures")

    xmin, xmax, ymin, ymax, dlon, dlat = fig_extent(nc)
    fwidth, fheight = fig_size(dlon, dlat)

    fig = plt.figure(figsize = (fwidth, fheight), constrained_layout = True)
    ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
    ax.set_title(f"Max tsunami elevation {fname} [m]", pad=0.01)
    #ax.coastlines(zorder=100)
    disp = ax.imshow(nc.data / scale_ratio, vmin=0, vmax=np.nanmax(nc.data / scale_ratio),
            extent=[xmin, xmax, ymin, ymax], cmap="hot_r",
            origin="lower")
    fig.colorbar(disp, orientation="vertical",
            extend="max", pad=0, shrink=0.8)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_facecolor("gray")

    if dem is not None and XX is not None and YY is not None:
        ax.contour(XX, YY, dem[::-1], [0], colors='white', zorder=100, linewidths=2)

    #print(scale_ratio)
    #print(nc.data/scale_ratio)
    
    fout = os.path.join(fig_dir, f"max_footprint__{fname}.png")
    fig.savefig(fout, dpi=300)
    plt.close()

    return nc

def tsunami_footprint_timeseries(infile, vmin, vmax, t0, t1, scale_ratio, dem, XX, YY):
    print(infile, vmin, vmax, t0, t1)
    nc = xr.open_dataset(infile)
    nc = nc["wave_height"]

    fname = infile.name[:-3]
    parent_dir = infile.parent
    fig_dir = os.path.join(parent_dir, "figures")

    xmin, xmax, ymin, ymax, dlon, dlat = fig_extent(nc)
    fwidth, fheight = fig_size(dlon, dlat)

    for ii, tt in enumerate(nc.time[t0:t1]):
        time = pd.to_datetime(tt.values)
        string = f"{time.hour:02d}:{time.minute:02d} __ {fname}"
        fig = plt.figure(figsize = (fwidth, fheight), constrained_layout = True)
        ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
        ax.set_title(f"{string} {fname} [m]", pad=0.01)
        ax.coastlines(zorder=100)
        disp = ax.imshow(nc.data[ii] / scale_ratio, 
                vmin=vmin, vmax=vmax,
                extent=[xmin, xmax, ymin, ymax], cmap="RdBu_r",
                origin="lower")
        fig.colorbar(disp, orientation="vertical",
            extend="both", pad=0, shrink=0.8)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        if dem is not None and XX is not None and YY is not None:   
            ax.contour(XX, YY, dem[::-1], [0], colors='white', zorder=100, linewidths=2)

        fout = os.path.join(fig_dir, f"footprint_{ii:04d}__{fname}.png")
        fig.savefig(fout)
        plt.close(fig)

    ### convert to gif
    inf = os.path.join(fig_dir, f"footprint_*__{fname}.png")
    cmd = f"magick -delay 20 -loop 0 {inf} {fig_dir}/animation__{fname}.gif"
    os.system(cmd)
    cmd = f"rm -f {inf}"
    os.system(cmd)

    return nc

def load_dem(dem_f):
    ds = xr.open_dataset(dem_f)
    res_x, res_y = ds.spacing.values
    nx, ny = ds.dimension.values
    xmin, xmax = ds.x_range.values.ravel()
    ymin, ymax = ds.y_range.values.ravel()
    x_range = np.linspace(xmin, xmax, nx)
    y_range = np.linspace(ymin, ymax, ny)
    z_values = ds.z.values
    XX , YY = np.meshgrid(x_range, y_range)
    dem = z_values.reshape(XX.shape)
    return dem, XX, YY


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--jagurs_nc", type=str,
            default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/jagurs_runs/uji_coba/tsunami_footprints/sample_jawa__2700m_to_0900m/footprints_stochastic_sources__Mw_8.400000__Lon_107.288800__Lat_-7.132900__table_simplified/footprints__sample_num__21.nc",
            help = "JAGURS nc output file")
    parser.add_argument("--what_to_check", type=int,
            default = 1,
            help = "1: initial displacement; 2: tsunami max footprint; 3: tsunami footprint timeseries")
    parser.add_argument("--vmin_vmax", type=float, nargs="+",
            default = [-1, 1],
            help = "vmin and vmax for colourmap")
    parser.add_argument("--t0_t1", type=int, nargs="+",
            default = [20, 30],
            help = "range of time step to plot footprint timeseries, required when what_to_check is 3")
    parser.add_argument("--scale_ratio", type=float,
            default = 10000.,
            help = "In order to save some space, JAGURS outputs is saved in an integer format, multiplied it by 10,000. In this case, --scale_ratio will be 10000")
    parser.add_argument("--dem_file", type=str,
            default = None,
            help = "to generate coastline for plotting")
    args = parser.parse_args()

    if args.dem_file is not None:
        dem, XX, YY = load_dem(args.dem_file)
    else:
        dem, XX, YY = None, None, None

    if os.path.exists(args.jagurs_nc):
        infile = Path(args.jagurs_nc)
        parent_dir = infile.parent
        fig_dir = Path(os.path.join(parent_dir, "figures"))
        fig_dir.mkdir(exist_ok = True)
        if args.what_to_check == 1:
            nc = initial_displacement(infile, args.vmin_vmax[0], args.vmin_vmax[1], args.scale_ratio, dem, XX, YY)
        elif args.what_to_check == 2:
            nc = tsunami_max_footprint(infile, args.scale_ratio, dem, XX, YY)
        elif args.what_to_check == 3:
            nc = tsunami_footprint_timeseries(infile, args.vmin_vmax[0], args.vmin_vmax[1], 
                    args.t0_t1[0], args.t0_t1[1], args.scale_ratio, dem, XX, YY)
        else:
            print("WRONG what_to_check CODE!")
            sys.exit()

        nc.close()
    else:
        print(f"{args.jagurs_nc} IS NOT AVAILABLE! STOP!")
        sys.exit()


