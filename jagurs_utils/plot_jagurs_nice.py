"""
Ryan Pranantyo
EOS, April 2025

a nice quick plot for presentation to show result from JAGURS

outputs: initial displacement and max sea surface height
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

to_plot = ['initial_displacement', 'max_height']

jagurs_nc = Path('/home/ignatius.pranantyo/Tsunamis/SensitivityTests__ModelResolutions/2006_SouthernJava/runs__DEM__DeltaDTM_GEBCO_BATNAS/Java_medRes__Inundation4grids__adjustedDeltaDTM__Supendi2022/SD00.nc')
where_to_save = Path(os.path.join(jagurs_nc.parent, 'figures'))
where_to_save.mkdir(exist_ok = True)

### loading basic config
ds = xr.open_dataset(jagurs_nc)
xmin = ds.lon.min().values
xmax = ds.lon.max().values
ymin = ds.lat.min().values
ymax = ds.lat.max().values
dlon = abs(xmax - xmin)
dlat = abs(ymax - ymin)

### figsize
fig_rat = dlon / dlat
fig_width = 10.
if fig_rat <= 1.1 and fig_rat >= 0.9:
    fig_height = fig_width - 0.1
elif fig_rat > 1.1:
    fig_height = (fig_width / fig_rat) - 0.035
else:
    fig_height = (fig_width / fig_rat) + 0.035

for var in to_plot:
    print(var)
    da = ds[var]
    fig = plt.figure(figsize = (fig_width, fig_height), constrained_layout = True)
    ax = fig.add_subplot(1,1,1, projection = ccrs.PlateCarree())

    if var == 'initial_displacement':
        title = 'initial displacement (m)'
        cmap = 'RdBu_r'
        vmin = -5; vmax = 5
        cb_extend = 'both'
        fout = os.path.join(where_to_save, f'initial_displacement_{jagurs_nc.name[:-3]}.png')
        var2plot = da.data[0]

        ax.coastlines()
    elif var == 'max_height':
        title = 'maximum sea surface elevation (m)'
        cmap = 'hot_r'
        vmin = 0; vmax = 10
        cb_extend = 'max'
        fout = os.path.join(where_to_save, f'max_footprint_{jagurs_nc.name[:-3]}.png')
        var2plot = da.data
    else:
        print('something wrong > exit.')
        sys.exit()

    ax.set_title(title, pad = 0)
    plot = ax.imshow(var2plot, vmin = vmin, vmax = vmax,
            cmap = cmap, origin = 'lower',
            extent = [xmin, xmax, ymin, ymax])
    fig.colorbar(plot, orientation = 'vertical', 
            extend = cb_extend, pad = 0, shrink = 0.85)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    fig.savefig(fout)
    plt.close()


