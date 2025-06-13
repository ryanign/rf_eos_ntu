"""
Ryan Pranantyo
EOS, March 2025

a basic script to check sfincs results
- based on my old script built at Reask
"""
import os, sys
import numpy as np
import pandas as pd
import geopandas as gpd

from hydromt_sfincs import SfincsModel
from hydromt_sfincs import utils
from hydromt import DataCatalog

import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import rasterio
import rioxarray

from pathlib import Path

### what to check
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--sfincs_model", type = str,
        default = "/home/ignatius.pranantyo/Tsunamis/Learnings/SFINCS_examples/run_exercise/Pangandaran-small_testing__Pangandaran2006_testing__resolution__30_m",
        help = "full path to sfincs model run")
parser.add_argument("--find_unstable", type = bool, 
        default = True,)
parser.add_argument("--steps", type = int,
        default = 6,
        help = "create snapshots every N steps")
args = parser.parse_args()

sfincs_root = args.sfincs_model
# where to save
fig_dir = Path(os.path.join(sfincs_root, "figs"))
fig_dir.mkdir(exist_ok = True)

### main core ###
sf = SfincsModel(sfincs_root, mode = "r")
sf.read_results()

utm_zone = sf.crs.utm_zone
utm = utm_zone[:-1]
hemisphere = utm_zone[-1]

# mask water depth
hmin = 0.05
da_h = sf.results["h"].copy()
da_h = da_h.where(da_h > hmin).drop("spatial_ref")
da_h.attrs.update(long_name = "flood_depth", unit = "m")

dem = sf.results["zb"].copy()
dem.attrs.update(long_name = "DEM", unit = "m")
X,Y = np.meshgrid(dem.x, dem.y)

width = 8
height = (da_h.y.shape[0] / da_h.x.shape[0]) * width

# plot maximum water depth
da_h = sf.results["hmax"].copy()
da_h = da_h.where(da_h > hmin).drop("spatial_ref")
da_h.attrs.update(long_name = "max_flood_depth", unit = "m")

plt.close()
fig = plt.figure(figsize = (width, height), constrained_layout = True)
if utm_zone.upper() == "N":
    ax = fig.add_subplot(1,1,1,
            projection = ccrs.UTM(zone = utm,
                southern_hemishpere = False))
else:
    ax = fig.add_subplot(1,1,1,
            projection = ccrs.UTM(zone = utm,
                southern_hemisphere = True))
ax.set_title("Maximum flood depth (m)", loc = "left")
ax.contour(X, Y, dem, levels = np.arange(-10, np.nanmax(dem), 10), colors = "gray",
        zorder = 10, linewidths = 0.3)
ax.contour(X, Y, dem, levels = [0], colors = "black",
        zorder = 10, linewidths = 0.5)

flood_depth = np.where(dem <=0, np.nan, da_h[0].values)
var = ax.imshow(flood_depth, origin = "lower",
        extent = [da_h.x.min(), da_h.x.max(), da_h.y.min(), da_h.y.max()],
        vmin = 0, vmax = 10.)
fig.colorbar(var, orientation = "vertical", shrink = 0.5, pad=0, extend='max')
fout = os.path.join(fig_dir, "flood_depth_maximum.png")
fig.savefig(fout, dpi=200)
plt.close()

### export to GTiff
sf.write_raster("results.hmax", compress="LZW")
land_flood = da_h.where(dem >=0, np.nan)
land_flood_epsg = land_flood.rio.reproject("EPSG:4326")
land_flood_epsg.rio.to_raster(f"{sfincs_root}/gis/hmax_onland__epsg4326.tif", driver="GTiff", compress="LZW")

def plot_snapshots(t, da_h, dem):
    print(t)
    fig = plt.figure(figsize = (width, height), constrained_layout = True)
    if utm_zone.upper() == "N":
        ax = fig.add_subplot(1,1,1,
                projection = ccrs.UTM(zone = utm,
                    southern_hemishpere = False))
    else:   
        ax = fig.add_subplot(1,1,1,
                projection = ccrs.UTM(zone = utm,
                    southern_hemisphere = True))

    ax.set_title(f"time step {t}", loc = "left")
    ax.contour(X, Y, dem, levels = np.arange(-10, np.nanmax(dem), 10), colors = "gray",
            zorder = 10, linewidths = 0.3)
    ax.contour(X, Y, dem, levels = [0], colors = "black",
            zorder = 10, linewidths = 0.5)

    water = da_h.values[t]
    water = np.where(dem <= 0, water + dem, water)
    var = ax.imshow(water, origin = "lower",
            extent = [da_h.x.min(), da_h.x.max(), da_h.y.min(), da_h.y.max()],
            vmin = -0.5, vmax = 0.5,
            cmap = "RdBu_r")
    fig.colorbar(var, shrink = 0.5, pad = 0, extend = "both")
    fout = os.path.join(fig_dir, f"flood_depth__{t:05d}.png")
    fig.savefig(fout, dpi=300)
    plt.close(fig)


if args.steps == 0:
    sys.exit()

### plot snapshots
da_h = sf.results["h"].copy()
da_h = da_h.where(da_h > hmin).drop("spatial_ref")
da_h.attrs.update(long_name = "flood_depth", unit = "m")
#da_h = np.where(dem <=0, da_h + dem, da_h)
from joblib import Parallel, delayed
Parallel(n_jobs = 8)(delayed(plot_snapshots)(t, da_h, dem) for t in range(0, da_h.time.size, args.steps))

inf = os.path.join(fig_dir, f"flood_depth__*.png")
cmd = f"magick -delay 20 -loop 0 {inf} {fig_dir}/animation.gif"
os.system(cmd)
