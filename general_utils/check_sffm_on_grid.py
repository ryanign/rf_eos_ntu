"""
it is still confusing to understand the 'fault index' of the SFFM obtained from the RPTHA code.

here, I would like to check whether I am mapping the correct slip on the grid source


OUTCOME:
    my way was correct!
    the fault index should 'sort' by along strike first followed by downdip
"""
import os, sys
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


### grid
grid_source_f = '/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/OUTPUTS__Slab2__Jawa/unit_source_grid/SLAB2__Jawa.shp'
grid_gdf = gpd.read_file(grid_source_f)

### sffm example
sffm_f = '/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/OUTPUTS__Slab2__Jawa/SourceCombinations__20250523/stochastic_slips__SLAB2__Jawa/stochastic_sources__Mw_8.000000__Lon_116.894850__Lat_-9.251680__table.csv'
sffm_df = pd.read_csv(sffm_f)


### grid sort by downdip
grid_downdip = grid_gdf.sort_values(by = ['dwndp_n', 'alngst_']).reset_index()
grid_downdip['unit_source_index'] = np.arange(len(grid_gdf)).astype(int) + 1

### grid sort by alongstrike
grid_strike = grid_gdf.sort_values(by = ['alngst_', 'dwndp_n']).reset_index()
grid_strike['unit_source_index'] = np.arange(len(grid_gdf)).astype(int) + 1

grid_xmin = grid_gdf.bounds.minx.min() - 1
grid_xmax = grid_gdf.bounds.maxx.max() + 1
grid_ymin = grid_gdf.bounds.miny.min() - 1
grid_ymax = grid_gdf.bounds.maxy.max() + 1

### quick check
where_to_save = '../../../sandpits'

df = pd.DataFrame()
df['unit_source_index'] = np.arange(len(grid_gdf)).astype(int) + 1

for src in sffm_df.index:
    print(src)

    flt_index_str = sffm_df.event_index_string[src]
    flt_slip_str = sffm_df.event_slip_string[src]
    flt_index = list(map(int, flt_index_str.split('-')[:-1]))
    flt_slip = list(map(float, flt_slip_str.split('_')[:-1]))
    
    tmp_df = pd.DataFrame(
            data = {f'unit_source_slip__{src}' : flt_slip,
                     'unit_source_index' : flt_index},
                index = flt_index)

    grid_downdip['slip'] = tmp_df[f'unit_source_slip__{src}']
    grid_strike['slip'] = tmp_df[f'unit_source_slip__{src}']

    target_lon = sffm_df['target_lon'][src]
    target_lat = sffm_df['target_lat'][src]

    fig = plt.figure(figsize=(15, 5), constrained_layout = True)

    for ii, var2plot in enumerate([grid_downdip, grid_strike]):
        ax = fig.add_subplot(1,2,ii+1, projection=ccrs.PlateCarree())
        ax.coastlines()
        ax.scatter(target_lon, target_lat, marker = '*', s= 50, c='green', zorder=99)
        var2plot.plot(
                ax = ax,
                column = 'slip',
                legend = True,
                legend_kwds = {
                    'label' : 'slip, m',
                    'orientation' : 'horizontal',
                    'shrink' : 0.5,
                    'pad' : 0},
                cmap = 'inferno_r'
                )

        var2plot.plot(
                ax = ax,
                fc = 'none',
                ec = 'gray',
                linewidth = 0.1)
        ax.set_xlim(grid_xmin, grid_xmax)
        ax.set_ylim(grid_ymin, grid_ymax)

    fout = os.path.join(where_to_save, f'sffm__{src:08d}.png')
    fig.savefig(fout)
    plt.close()
                    
