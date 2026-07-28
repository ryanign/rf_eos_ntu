"""
Ryan Pranantyo
EOS, 27 July 2026

basic post processing of ground motion field output from scenario based calcuation

usage:
    python post-processing_gmf-data.py --input_file ./outputs/gmf-data_22.csv

"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import argparse
import xarray as xr
import cartopy.crs as ccrs
from pathlib import Path
from datetime import datetime

#-----------------------------------------
# BASIC VARIABLES
#-----------------------------------------
PERCENTILES = [2.5, 5, 16, 25, 50, 75, 84, 95, 97.5]

MMI_COLORS = {
    1 : '#FFFFFF',   # I      — not felt
    2 : '#BFCCFF',   # II     — weak
    3 : '#BFCCFF',   # III    — weak
    4 : '#A0F0A0',   # IV     — light
    5 : '#FFFF00',   # V      — moderate
    6 : '#FFC800',   # VI     — strong
    7 : '#FF9100',   # VII    — very strong
    8 : '#FF0000',   # VIII   — severe
    9 : '#C80000',   # IX     — violent
    10: '#800000',   # X+     — extreme
}

MMI_LABELS = {
    1 : 'I',
    2 : 'II',
    3 : 'III',
    4 : 'IV',
    5 : 'V',
    6 : 'VI',
    7 : 'VII',
    8 : 'VIII',
    9 : 'IX',
    10: 'X+',
}


#-----------------------------------------
# FUNCTIONS
#-----------------------------------------

def parse_realizations(rlz_file):
    """
    Parse realizations CSV from OpenQuake.
    """
    df = pd.read_csv(rlz_file, comment='#')
    df.columns = df.columns.str.strip()

    # Extract GMPE name from branch_path
    df['gmpe_name'] = (df['branch_path']
            .str.replace('[','', regex=False)
            .str.replace(']','', regex=False)
            .str.strip())
    return df[['rlz_id', 'gmpe_name', 'weight']]

def build_event_map(rlz_df, n_events_total):
    """
    Map rlz_id -> (gmpe_name, weight, event_id_start, event_id_end).
    Assumes sequential event_id assignment per rlz_id.
    """
    n_gmpe = len(rlz_df)
    n_sim  = n_events_total // n_gmpe

    event_map = {}
    for i, row in rlz_df.iterrows():
        start = row['rlz_id'] * n_sim
        end   = start + n_sim - 1
        event_map[row['rlz_id']] = {
                'gmpe'   : row['gmpe_name'],
                'weight' : row['weight'],
                'start'  : start,
                'end'    : end,
                }
    return event_map, n_sim

def weighted_quantile(values, quantiles, weights):
    """
    Compute weighted quantiles
    """
    values    = np.asarray(values,  dtype=float)
    weights   = np.asarray(weights, dtype=float)
    quantiles = np.asarray(quantiles)

    # sort by values
    sorter  = np.argsort(values)
    values  = values[sorter]
    weights = weights[sorter]

    # cumulative weights normalised
    cumw = np.cumsum(weights)
    cumw = (cumw - 0.5 * weights) / cumw[-1]

    return np.interp(quantiles, cumw, values)

def compute_stats(x, percentiles, weights=None):
    """
    Statistic log-normal for 1 column of GMV.
    If weights provided, uses weighted mean and weighted quantiles.
    weights : per-realisation weights (same length as x)
    """
    x = np.asarray(x, dtype=float)
    ln_x = np.log(x)

    if weights is not None:
        weights = np.asarray(weights, dtype=float)
        weights = weights / weights.sum()
        mu_ln   = np.average(ln_x, weights=weights)
        var_ln  = np.average((ln_x - mu_ln)**2, weights=weights)
        sig_ln  = np.sqrt(var_ln)
        q_vals  = weighted_quantile(x, [p/100 for p in percentiles], weights)
    else:
        mu_ln   = ln_x.mean()
        sig_ln  = ln_x.std()
        q_vals  = np.percentile(x, percentiles)

    stats  = {
            'mean' : np.exp(mu_ln),
            'std'  : np.exp(sig_ln),
            'min'  : x.min(),
            'max'  : x.max(),
            }

    for p, q in zip(percentiles, q_vals):
        stats[f'p{str(p).replace(".","pt")}'] = q

    return pd.Series(stats)

def build_xarray(stats_dfs, site_coords, imts, imt_cols,
        percentiles, n_realizations, group_name, attrs_extra=None):
    """
    Build xarray
    """
    sites      = site_coords['custom_site_id'].values
    stat_names = (['mean', 'std', 'min', 'max'] + 
                  [f'p{str(p).replace(".","pt")}' for p in percentiles])
    data_vars  = {}

    # Coordinates
    data_vars['lon']   = xr.DataArray(site_coords['lon'].values,   dims='site')
    data_vars['lat']   = xr.DataArray(site_coords['lat'].values,   dims='site')
    data_vars['depth'] = xr.DataArray(site_coords['depth'].values, dims='site')
    data_vars['vs30']  = xr.DataArray(site_coords['vs30'].values,  dims='site')

    # IMT statistics
    for imt, col in zip(imts, imt_cols):
        df_stat   = stats_dfs[imt].sort_values('custom_site_id')
        imt_clean = (imt.replace('(','_').replace(')','').replace('.','p'))
        for stat in stat_names:
            varname = f'{imt_clean}_{stat}'
            data_vars[varname] = xr.DataArray(
                    df_stat[stat].values,
                    dims='site',
                    attrs={'imt': imt, 'statistic': stat, 'units': 'g'},
                    )

    attrs = {
            'group'          : group_name,
            'n_realizations' : n_realizations,
            'imts'           : ', '.join(imts),
            'percentiles'    : str(percentiles),
            'created'        : datetime.now().isoformat(),
            'source'         : 'OpenQuake Engine - Scenario calculation',
            'Conventions'    : 'CF-1.8',
            }

    if attrs_extra:
        attrs.update(attrs_extra)

    ds = xr.Dataset(data_vars, coords={'site': sites}, attrs=attrs)
    return ds

def get_pga_mmi_bounds():
    """
    PGA (%g) to MMI bounds from Wald et al. (1999) Table 1.
    Returns list of (threshold_%g, mmi) tuples, sorted ascending.
    """
    return [
        (0.17,  1),
        (1.4,   2),
        (3.9,   4),
        (9.2,   5),
        (18.0,  6),
        (34.0,  7),
        (65.0,  8),
        (124.0, 9),
    ]

def get_pgv_mmi_bounds():
    """
    PGV (cm/s) to MMI bounds from Wald et al. (1999) Table 1.
    Returns list of (threshold_cm/s, mmi) tuples, sorted ascending.
    """
    return [
        (0.1,   1),
        (1.1,   2),
        (3.4,   4),
        (8.1,   5),
        (16.0,  6),
        (31.0,  7),
        (60.0,  8),
        (116.0, 9),
    ]

def pga_to_mmi(pga_g):
    pga_pct = np.asarray(pga_g) * 100.0
    mmi     = np.ones_like(pga_pct, dtype=int)
    for threshold, intensity in get_pga_mmi_bounds():
        mmi[pga_pct >= threshold] = intensity
    mmi[pga_pct >= 124.0] = 10
    return mmi

def pgv_to_mmi(pgv_cms):
    pgv = np.asarray(pgv_cms)
    mmi = np.ones_like(pgv, dtype=int)
    for threshold, intensity in get_pgv_mmi_bounds():
        mmi[pgv >= threshold] = intensity
    mmi[pgv >= 116.0] = 10
    return mmi

def quick_plot(ds, imts, imt_cols, output_file='gmf_quickplot.png', title_suffix=''):
    """
    plotting all IMTs available
    """

    n_rows = len(imts)
    fig, axes = plt.subplots(n_rows, 3, figsize=(15, 4 * n_rows))

    if n_rows == 1:
        axes = axes[np.newaxis, :]

    col_stats  = ['p2pt5', 'mean', 'p97pt5']
    col_titles = ['Lower CI (2.5%)', 'Mean', 'Upper CI (97.5%)']

    lon = ds['lon'].values
    lat = ds['lat'].values

    for row, (imt, col) in enumerate(zip(imts, imt_cols)):
        imt_clean = (imt.replace('(', '_')
                        .replace(')', '')
                        .replace('.', 'p'))

        vals = [ds[f'{imt_clean}_{s}'].values for s in col_stats]
        vmin = min(v.min() for v in vals)
        vmax = max(v.max() for v in vals)

        for col_idx, (stat, title) in enumerate(zip(col_stats, col_titles)):
            ax  = axes[row, col_idx]
            z   = ds[f'{imt_clean}_{stat}'].values

            sc = ax.scatter(lon, lat, c=z, cmap='hot_r',
                            vmin=vmin, vmax=vmax, s=10)

            plt.colorbar(sc, ax=ax, label='g', shrink=0.8)
            ax.set_xlabel('Longitude', fontsize=8)
            ax.tick_params(labelsize=7)
            ax.grid(True, alpha=0.3, linewidth=0.5)

            if row == 0:
                ax.set_title(title, fontsize=10, fontweight='bold', pad=8)
            if col_idx == 0:
                ax.set_ylabel(f'{imt}\nLatitude', fontsize=8)

    plt.suptitle(f'Scenario GMF Statistics — 95% Confidence Interval {title_suffix}',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.show()
    print(f"Saved: {output_file}")

def quick_plot_pga_mmi(ds, output_file='gmf_PGA_MMI.png', title_suffix=''):
    """
    plotting PGA and MMI
    """

    if 'PGA_mean' not in ds:
        print('Skipping: PGA not found in dataset')
        return

    lon = ds['lon'].values
    lat = ds['lat'].values
    col_stats  = ['p2pt5', 'mean', 'p97pt5']
    col_titles = ['Lower CI (2.5%)', 'Mean', 'Upper CI (97.5%)']

    mmi_levels = list(range(1,11))
    cmap_mmi   = mcolors.ListedColormap([MMI_COLORS[m] for m in mmi_levels])
    norm_mmi   = mcolors.BoundaryNorm(
                    boundaries=[0.5 + i for i in range(11)],
                    ncolors=len(mmi_levels))

    vals = [ds[f'PGA_{s}'].values for s in col_stats]
    vmin = min(v.min() for v in vals)
    vmax = max(v.max() for v in vals)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))

    for col_idx, (stat, title) in enumerate(zip(col_stats, col_titles)):
        ax = axes[0, col_idx]
        sc = ax.scatter(lon, lat, c=ds[f'PGA_{stat}'].values,
                        cmap='hot_r', vmin=vmin, vmax=vmax, s=10)
        plt.colorbar(sc, ax=ax, label='g', shrink=0.8)
        ax.set_title(title, fontsize=10, fontweight='bold', pad=8)
        ax.set_xlabel('Longitude', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3, linewidth=0.5)
        if col_idx == 0:
            ax.set_ylabel('PGA\nLatitude', fontsize=8)

    for col_idx, stat in enumerate(col_stats):
        ax      = axes[1, col_idx]
        mmi_val = pga_to_mmi(ds[f'PGA_{stat}'].values)
        sc      = ax.scatter(lon, lat, c=mmi_val,
                             cmap=cmap_mmi, norm=norm_mmi, s=10)
        cbar    = plt.colorbar(sc, ax=ax, shrink=0.8, ticks=mmi_levels)
        cbar.set_ticklabels([MMI_LABELS[m] for m in mmi_levels])
        cbar.set_label('MMI', fontsize=9)
        ax.set_xlabel('Longitude', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3, linewidth=0.5)
        if col_idx == 0:
            ax.set_ylabel('MMI (PGA)\nLatitude', fontsize=8)

    plt.suptitle(f'Scenario GMF — PGA & MMI (Wald et al. 1999) {title_suffix}',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")

def quick_plot_pgv_mmi(ds, output_file='gmf_PGV_MMI.png', title_suffix=''):
    """
    quick plot PGV and MMI
    """
    if 'PGV_mean' not in ds:
        print("Skipping: PGV not found in dataset")
        return

    lon = ds['lon'].values
    lat = ds['lat'].values
    col_stats  = ['p2pt5', 'mean', 'p97pt5']
    col_titles = ['Lower CI (2.5%)', 'Mean', 'Upper CI (97.5%)']

    mmi_levels = list(range(1, 11))
    cmap_mmi   = mcolors.ListedColormap([MMI_COLORS[m] for m in mmi_levels])
    norm_mmi   = mcolors.BoundaryNorm(
                    boundaries=[0.5 + i for i in range(11)],
                    ncolors=len(mmi_levels))

    vals = [ds[f'PGV_{s}'].values for s in col_stats]
    vmin = min(v.min() for v in vals)
    vmax = max(v.max() for v in vals)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))

    for col_idx, (stat, title) in enumerate(zip(col_stats, col_titles)):
        ax = axes[0, col_idx]
        sc = ax.scatter(lon, lat, c=ds[f'PGV_{stat}'].values,
                        cmap='hot_r', vmin=vmin, vmax=vmax, s=10)
        plt.colorbar(sc, ax=ax, label='cm/s', shrink=0.8)
        ax.set_title(title, fontsize=10, fontweight='bold', pad=8)
        ax.set_xlabel('Longitude', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3, linewidth=0.5)
        if col_idx == 0:
            ax.set_ylabel('PGV\nLatitude', fontsize=8)

    for col_idx, stat in enumerate(col_stats):
        ax      = axes[1, col_idx]
        mmi_val = pgv_to_mmi(ds[f'PGV_{stat}'].values)
        sc      = ax.scatter(lon, lat, c=mmi_val,
                             cmap=cmap_mmi, norm=norm_mmi, s=10)
        cbar    = plt.colorbar(sc, ax=ax, shrink=0.8, ticks=mmi_levels)
        cbar.set_ticklabels([MMI_LABELS[m] for m in mmi_levels])
        cbar.set_label('MMI', fontsize=9)
        ax.set_xlabel('Longitude', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3, linewidth=0.5)
        if col_idx == 0:
            ax.set_ylabel('MMI (PGV)\nLatitude', fontsize=8)

    plt.suptitle(f'Scenario GMF — PGV & MMI (Wald et al. 1999) {title_suffix}',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Post-processing OpenQuake scenario GMF',
            )
    parser.add_argument(
            '--input_file',
            default = './outputs/gmf-data_34.csv',
            help = 'GMF output file from scenario based calculation.'
            )
    args = parser.parse_args()


    #-- Paths --
    input_f = Path(args.input_file)
    calc_id = int(input_f.stem.split('_')[-1])
    outdir  = input_f.parent
    rlz_f   = outdir / f'realizations_{calc_id}.csv'
    site_f  = outdir / f'site_model_{calc_id}.csv'

    print(f'cacl_id     : {calc_id}')
    print(f'GMF         : {input_f}')
    print(f'Realization : {rlz_f}'  )
    print(f'Sites       : {site_f}' )

    #-- Load data --
    df      = pd.read_csv(input_f, skiprows=1)
    site_df = pd.read_csv(site_f,  skiprows=1)
    rlz_df  = parse_realizations(rlz_f)

    print(f'\nGMPEs found ({len(rlz_df)}):')
    for _, row in rlz_df.iterrows():
        print(f'  rlz_id={row["rlz_id"]} w={row["weight"]:.4f} {row["gmpe_name"]}')

    #-- Merge site coordinates
    df = df.merge(
            site_df[['custom_site_id', 'lon', 'lat', 'depth', 'vs30']],
            on='custom_site_id', how='left'
            )

    #sys.exit()

    #-- Detect intensity meassured (IMTs) --
    imt_cols = [c for c in df.columns if c.startswith('gmv_')]
    imts     = [c.replace('gmv_','') for c in imt_cols]
    print(f'\nIMTs: {imts}')

    #-- Event mapping --
    n_events_total = df['event_id'].nunique()
    event_map, n_sim = build_event_map(rlz_df, n_events_total)

    print(f'\nEvent mapping (n_sim={n_sim} per-GMPE):')
    for rlz_id, info in event_map.items():
        print(f'  rlz_id={rlz_id} event_id {info["start"]}-{info["end"]} '
              f'w={info["weight"]:.4f} {info["gmpe"]}')

    #-- Site coordinates --
    site_coords = (df.groupby('custom_site_id')
            .agg(lon=('lon', 'first'),
                 lat=('lat', 'first'),
                 depth=('depth', 'first'),
                 vs30=('vs30', 'first'))
            .reset_index()
            .sort_values('custom_site_id'))
    sites = site_coords['custom_site_id'].values

    #-----------------------
    # compute stats per GMPE
    #-----------------------
    datasets = {} # {group_name: xr.Dataset}

    for rlz_id, info in event_map.items():
        gmpe   = info['gmpe']
        weight = info['weight']
        print(f'\n[{gmpe}] computing stats ...')

        mask    = (df['event_id'] >= info['start']) & (df['event_id'] <= info['end'])
        df_gmpe = df[mask].copy()

        # Equal weights within GMPE
        per_event_w = np.ones(n_sim) / n_sim

        stats_df = {}
        for col in imt_cols:
            imt = col.replace('gmv_', '')
            print(f'   {imt}')
            raw = (df_gmpe.groupby('custom_site_id')[col]
                    .apply(lambda g: compute_stats( 
                        g.values, 
                        PERCENTILES, 
                        weights=None))
                    .reset_index())
            stats_df[imt] = raw.pivot(
                index='custom_site_id',
                columns='level_1',
                values=col).reset_index()

        datasets[gmpe] = build_xarray(
            stats_df, site_coords, imts, imt_cols,
            PERCENTILES, n_sim, gmpe,
            attrs_extra={'gmpe': gmpe, 'weight': float(weight)})

    #-----------------------
    # compute weighted
    #-----------------------
    print('\n[all] computing weighted pooled stats ...')

    # assign per-realisation weight to each row
    weight_map = {}
    for rlz_id, info in event_map.items():
        per_event_w = info['weight'] / n_sim
        for eid in range(info['start'], info['end'] + 1):
            weight_map[eid] = per_event_w

    df['_w'] = df['event_id'].map(weight_map)

    stats_dfs_all = {}
    for col in imt_cols:
        imt = col.replace('gmv_', '')
        print(f'   {imt}')
        raw = (df.groupby('custom_site_id')
                .apply(lambda g: compute_stats(
                    g[col].values,
                    PERCENTILES,
                    weights=g['_w'].values),
                    include_groups=False)
                .reset_index())
        stats_dfs_all[imt] = raw

    datasets['all'] = build_xarray(
            stats_dfs_all, site_coords, imts, imt_cols,
            PERCENTILES, n_events_total, 'all',
            attrs_extra={
                'description' : 'Weighted pooled statistics across all GMPEs',
                'gmpes'       : ', '.join(event_map[r]['gmpe'] for r in event_map),
                'weights'     : ', '.join(str(event_map[r]['weight']) for r in event_map),
                })
    

    #--------------
    # saving
    #--------------
    output_nc = f'gmf_statistics__{calc_id}.nc'
    print(f'\n{"="*60}')
    print(f'Saving NetCDF: {output_nc}')

    # write each group
    for group_name, ds in datasets.items():
        mode = 'w' if group_name == list(datasets.keys())[0] else 'a'
        ds.to_netcdf(
                output_nc,
                group=group_name,
                mode=mode,
                encoding={v: {'zlib': True, 'complevel': 4} for v in ds.data_vars})
        print(f'  Written group: /{group_name}/')

    # CSV export
    output_csv = f'gmf_statistics__{calc_id}.csv'
    df_out = datasets['all'].to_dataframe().reset_index()
    df_out.to_csv(output_csv, index=False)
    print(f'  CSV (all): {output_csv}')

    #--------------
    # quick plots
    #--------------
    ds_all = datasets['all']
    print('\nGenerating plots ...')
    quick_plot(
            ds_all, imts, imt_cols, 
            output_file=f'gmf_quickplot__{calc_id}.png',
            title_suffix='(all GMPEs weighted)'
            )

    quick_plot_pga_mmi(
            ds_all,
            output_file=f'gmf_quickplot_pga-mmi__{calc_id}.png',
            title_suffix='(all GMPEs weighted)',
            )

    quick_plot_pgv_mmi(
            ds_all,
            output_file=f'gmf_quickplot_pgv-mmi__{calc_id}.png',
            title_suffix='(all GMPEs weighted)',
            )

    # per-GMPE
    for gmpe, ds_gmpe in datasets.items():
        if gmpe == 'all':
            continue
        gmpe_slug = gmpe.replace(' ','_')
        quick_plot(
                ds_gmpe, imts, imt_cols,
                output_file=f'gmf_quickplot__{calc_id}_{gmpe_slug}.png',
                title_suffix=f'({gmpe})'
                )

    print(f'\nDone!')
    print(f'  NetCDF : {output_nc}')
    print(f'  CSV    : {output_csv}')
    print(f'  Groups : {list(datasets.keys())}')


                        


    sys.exit()
    # saving

    # quick plot of imts
    quick_plot(ds, imts, imt_cols, output_file=f'gmf_quickplot__{calc_id}.png')

    # quick plot of mmi from pga
    quick_plot_pga_mmi(ds, output_file=f'gmf_quickplot_pga-mmi__{calc_id}.png')

    # quick plot of mmi from pgv
    quick_plot_pgv_mmi(ds, output_file=f'gmf_quickplot_pgv-mmi__{calc_id}.png')
    
    print(f' csv format    = {output_f}')












