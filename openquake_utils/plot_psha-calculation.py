"""
plot_psha.py
============
Quick plots for OpenQuake Classical PSHA results.

Reads from NetCDF built by build_psha_netcdf.py.

Plots:
  Plot 1 — Hazard curves  : N sites x M IMTs, mean + Q10/Q90 band
  Plot 2 — Hazard maps    : Q10 | mean | Q90 per IMT per RP (cartopy)
  Plot 3 — UHS            : 1 figure per site, all PoEs

Usage:
    python plot_psha.py --netcdf psha_summary_6.nc \
                        --sites  selected_sites.csv \
                        --calc_id 6
"""
import os
import re
import argparse
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pathlib import Path
from scipy.spatial import cKDTree

# ══════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════

def read_sites_csv(path):
    """Read selected sites CSV: site_label, lon, lat"""
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    return df


def find_nearest_sites(selected, ds):
    """
    For each selected site (lon, lat), find nearest site in dataset.
    Returns DataFrame with matched site indices and actual coordinates.
    """
    ds_lon = ds['lon'].values
    ds_lat = ds['lat'].values
    ds_ids = ds['site'].values

    tree = cKDTree(np.column_stack([ds_lon, ds_lat]))

    results = []
    for _, row in selected.iterrows():
        dist, idx = tree.query([row['lon'], row['lat']])
        results.append({
            'label'      : row['site_label'],
            'req_lon'    : row['lon'],
            'req_lat'    : row['lat'],
            'site_id'    : ds_ids[idx],
            'actual_lon' : ds_lon[idx],
            'actual_lat' : ds_lat[idx],
            'dist_km'    : dist * 111.0,
        })

    df = pd.DataFrame(results)
    print('\nSite matching:')
    for _, r in df.iterrows():
        print(f"  {r['label']:20s} "
              f"requested ({r['req_lon']:.4f}, {r['req_lat']:.4f}) "
              f"-> nearest ({r['actual_lon']:.4f}, {r['actual_lat']:.4f}) "
              f"dist={r['dist_km']:.2f} km")
    return df


def detect_imts_from_ds(ds):
    """Detect IMTs from variable names in hazard_curves group."""
    imts = []
    seen = set()
    for var in ds.data_vars:
        if var.endswith('_mean') and '_iml' not in var:
            imt_var = var.replace('_mean', '')
            if imt_var not in seen:
                imts.append(imt_var)
                seen.add(imt_var)
    return imts


def imt_var_to_label(imt_var):
    """SA1p0 -> SA(1.0), PGA -> PGA"""
    if imt_var == 'PGA':
        return 'PGA'
    m = re.match(r'SA(\d+)p(\d+)', imt_var)
    if m:
        return f'SA({m.group(1)}.{m.group(2)})'
    return imt_var


def build_rp_lines(return_periods, investigation_time=50.0):
    """
    Generate RP_LINES dict automatically from list of return periods.
    Returns: {rp: (linestyle, color, label)}
    """
    colors = ['#E67E22', '#E74C3C', '#8E44AD',
              '#27AE60', '#2980B9', '#F39C12']
    styles = ['--', ':', '-.', '--', ':', '-.']

    rp_lines = {}
    for i, rp in enumerate(sorted(return_periods)):
        poe_pct = round((1 - np.exp(-investigation_time / rp)) * 100, 1)
        label   = f'{rp}yr ({poe_pct}%/{int(investigation_time)}yr)'
        color   = colors[i % len(colors)]
        ls      = styles[i % len(styles)]
        rp_lines[rp] = (ls, color, label)

    return rp_lines


# ══════════════════════════════════════════════════════════════════════
# PLOT 1 — HAZARD CURVES
# ══════════════════════════════════════════════════════════════════════

def plot_hazard_curves(ds_curves, matched_sites, imts,
                       return_periods, calc_id, output_prefix,
                       where2save):
    """
    N rows (sites) x M cols (IMTs).
    Each panel: mean + Q10/Q90 shaded band + RP horizontal lines.
    Y-axis: PoE (linear), X-axis: IML (log).
    """
    n_sites      = len(matched_sites)
    n_imts       = len(imts)
    inv_time     = 50.0
    rp_lines     = build_rp_lines(return_periods, inv_time)

    fig, axes = plt.subplots(n_sites, n_imts,
                             figsize=(5 * n_imts, 4 * n_sites),
                             squeeze=False)

    for row, (_, site) in enumerate(matched_sites.iterrows()):
        sid = site['site_id']

        for col, imt_var in enumerate(imts):
            ax    = axes[row, col]
            label = imt_var_to_label(imt_var)

            # IML axis
            iml_dim = f'{imt_var}_iml'
            if iml_dim not in ds_curves:
                ax.set_visible(False)
                continue
            iml = ds_curves[iml_dim].values

            # Extract PoE per statistic
            poe = {}
            for stat in ['mean', 'q10', 'q90']:
                var = f'{imt_var}_{stat}'
                if var in ds_curves:
                    poe[stat] = ds_curves[var].sel(site=sid).values

            # Q10/Q90 band
            if 'q10' in poe and 'q90' in poe:
                ax.fill_between(iml, poe['q10'], poe['q90'],
                                alpha=0.25, color='#3498DB',
                                label='Q10–Q90')

            # Mean
            if 'mean' in poe:
                ax.plot(iml, poe['mean'],
                        color='#2C3E50', lw=2.0, label='Mean')

            # Return period horizontal lines
            for rp, (ls, color, rp_label) in rp_lines.items():
                poe_rp = 1 - np.exp(-inv_time / rp)
                ax.axhline(poe_rp, ls=ls, color=color,
                           lw=1.2, label=rp_label)

            ax.set_xscale('log')
            ax.set_xlim(iml.min(), iml.max())
            ax.set_ylim(0, 1)
            ax.grid(True, which='both', alpha=0.3, lw=0.5)
            ax.tick_params(labelsize=8)

            if row == 0:
                ax.set_title(label, fontsize=11, fontweight='bold')
            if col == 0:
                ax.set_ylabel(f"{site['label']}\nPoE in {int(inv_time)}yr",
                              fontsize=9)
            if row == n_sites - 1:
                ax.set_xlabel('IML (g)', fontsize=9)
            if row == 0 and col == n_imts - 1:
                ax.legend(fontsize=7, loc='upper right')

    plt.suptitle(f'Hazard Curves — calc #{calc_id}',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    out = os.path.join(where2save, f'{output_prefix}_hazard_curves.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {out}')


# ══════════════════════════════════════════════════════════════════════
# PLOT 2 — HAZARD MAPS (cartopy)
# ══════════════════════════════════════════════════════════════════════

def plot_hazard_maps(ds_maps, imts, return_periods,
                     matched_sites, calc_id, output_prefix,
                     where2save):
    """
    Per IMT per return period: 3 columns (Q10 | mean | Q90).
    Uses cartopy. Overlays selected sites as star markers.
    """
    stats       = ['q10', 'mean', 'q90']
    stat_titles = ['Q10 (lower)', 'Mean', 'Q90 (upper)']

    lon  = ds_maps['lon'].values
    lat  = ds_maps['lat'].values
    proj = ccrs.PlateCarree()

    # Map extent with buffer
    buf     = 0.5
    lon_min = lon.min() - buf
    lon_max = lon.max() + buf
    lat_min = lat.min() - buf
    lat_max = lat.max() + buf

    for imt_var in imts:
        label  = imt_var_to_label(imt_var)
        n_rows = len(return_periods)
        n_cols = 3

        fig = plt.figure(figsize=(6 * n_cols, 5 * n_rows))

        for row, rp in enumerate(return_periods):

            # Consistent vmin/vmax across stats for this RP
            all_vals = []
            for stat in stats:
                var = f'{imt_var}_{stat}_{rp}yr'
                if var in ds_maps:
                    all_vals.append(ds_maps[var].values)
            if not all_vals:
                continue
            vmin = min(v.min() for v in all_vals)
            vmax = max(v.max() for v in all_vals)

            for col, (stat, stitle) in enumerate(zip(stats, stat_titles)):
                ax_idx = row * n_cols + col + 1
                ax = fig.add_subplot(n_rows, n_cols, ax_idx,
                                     projection=proj)

                var = f'{imt_var}_{stat}_{rp}yr'
                if var not in ds_maps:
                    ax.set_visible(False)
                    continue

                z = ds_maps[var].values

                # ── Cartopy base map ──────────────────────────────
                ax.set_extent([lon_min, lon_max,
                               lat_min, lat_max], crs=proj)
                ax.add_feature(cfeature.LAND,
                               facecolor='#F5F5F0', zorder=0)
                ax.add_feature(cfeature.OCEAN,
                               facecolor='#D6EAF8', zorder=0)
                ax.add_feature(cfeature.COASTLINE,
                               lw=0.6, edgecolor='#444', zorder=2)
                ax.add_feature(cfeature.BORDERS,
                               lw=0.3, linestyle=':', zorder=2)
                gl = ax.gridlines(draw_labels=True, lw=0.3,
                                  alpha=0.5, x_inline=False,
                                  y_inline=False)
                gl.top_labels   = False
                gl.right_labels = False
                gl.xlabel_style = {'size': 6}
                gl.ylabel_style = {'size': 6}

                # ── GMV scatter ───────────────────────────────────
                sc = ax.scatter(lon, lat, c=z, cmap='hot_r',
                                vmin=vmin, vmax=vmax,
                                s=8, transform=proj, zorder=3)
                plt.colorbar(sc, ax=ax, label='g (PGA)',
                             shrink=0.7, pad=0.08)

                # ── Selected sites overlay ────────────────────────
                for _, site in matched_sites.iterrows():
                    ax.plot(site['actual_lon'], site['actual_lat'],
                            marker='*', ms=12, color='cyan',
                            markeredgecolor='black',
                            markeredgewidth=0.5,
                            transform=proj, zorder=5)
                    ax.text(site['actual_lon'] + 0.05,
                            site['actual_lat'] + 0.05,
                            site['label'],
                            fontsize=6, color='black',
                            transform=proj, zorder=5,
                            bbox=dict(boxstyle='round,pad=0.2',
                                      facecolor='white',
                                      alpha=0.7, lw=0))

                # ── Titles and labels ─────────────────────────────
                if row == 0:
                    ax.set_title(stitle, fontsize=10,
                                 fontweight='bold', pad=8)
                if col == 0:
                    ax.text(-0.12, 0.5, f'{rp}yr RP',
                            transform=ax.transAxes,
                            fontsize=9, va='center',
                            rotation=90)

        plt.suptitle(f'Hazard Maps — {label} — calc #{calc_id}',
                     fontsize=13, fontweight='bold', y=1.01)
        plt.tight_layout()
        out = os.path.join(where2save, f'{output_prefix}_hazard_map_{imt_var}.png')
        plt.savefig(out, dpi=150, bbox_inches='tight')
        plt.close()
        print(f'  Saved: {out}')


# ══════════════════════════════════════════════════════════════════════
# PLOT 3 — UHS
# ══════════════════════════════════════════════════════════════════════

def plot_uhs(ds_uhs, matched_sites, calc_id, output_prefix, where2save):
    """
    1 figure total.
    Rows    = sites (N)
    Columns = PoEs  (max 5)
    """
    # Fix timedelta → float
    periods_raw = ds_uhs['period'].values
    if np.issubdtype(periods_raw.dtype, np.timedelta64):
        periods = periods_raw.astype('timedelta64[ns]').astype(float) / 1e9
    else:
        periods = periods_raw.astype(float)

    # Detect PoE labels
    poe_labels = sorted(set(
        re.search(r'poe_(\w+)_mean', v).group(1)
        for v in ds_uhs.data_vars
        if re.search(r'poe_(\w+)_mean', v)))

    def poe_label_to_title(poe_label):
        return poe_label.replace('pct', '% PoE in 50yr')

    n_sites = len(matched_sites)
    n_poes  = min(len(poe_labels), 5)   # cap at 5 columns
    poe_labels = poe_labels[:n_poes]

    fig, axes = plt.subplots(n_sites, n_poes,
                             figsize=(5 * n_poes, 4 * n_sites),
                             squeeze=False)

    for row, (_, site) in enumerate(matched_sites.iterrows()):
        sid   = site['site_id']
        label = site['label']

        for col, poe_label in enumerate(poe_labels):
            ax = axes[row, col]

            sa = {}
            for stat in ['mean', 'q10', 'q90']:
                var = f'poe_{poe_label}_{stat}'
                if var in ds_uhs:
                    sa[stat] = ds_uhs[var].sel(site=sid).values

            if 'q10' in sa and 'q90' in sa:
                ax.fill_between(periods, sa['q10'], sa['q90'],
                                alpha=0.25, color='#3498DB',
                                label='Q10–Q90')

            if 'mean' in sa:
                ax.plot(periods, sa['mean'],
                        color='#2C3E50', lw=2.0,
                        marker='o', ms=4, label='Mean')

            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=8)
            ax.set_xlim(left=0)
            ax.set_ylim(bottom=0)

            #ax.set_xscale('log')
            #ax.set_yscale('log')

            # Column headers — top row only
            if row == 0:
                ax.set_title(poe_label_to_title(poe_label),
                             fontsize=10, fontweight='bold')
                ax.legend(fontsize=7, loc='upper right')

            # Row labels — left column only
            if col == 0:
                ax.set_ylabel(f'{label}\nSa (g)', fontsize=9)

            # X label — bottom row only
            if row == n_sites - 1:
                ax.set_xlabel('Period T (s)', fontsize=9)

    plt.suptitle(f'Uniform Hazard Spectra — calc #{calc_id}',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    out = os.path.join(where2save, f'{output_prefix}_uhs.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {out}')


# ══════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Quick plots for OpenQuake Classical PSHA results.')
    parser.add_argument('--netcdf',        #required=True,
                        default='/home/ignatius.pranantyo/PSHA/Singapore/Simple_Deterministic/exercise__psha-based__SumatraSubduction__multiGMPEs/outputs/psha_summary_6.nc',
                        help='NetCDF file from post-processing_psha-calculation.py')
    parser.add_argument('--sites',         #required=True,
                        default='selected_sites.csv',
                        help='CSV with selected sites: site_label, lon, lat')
    parser.add_argument('--calc_id',       #required=True, 
                        type=int,
                        default=6,
                        help='Calculation ID (for plot titles)')
    parser.add_argument('--output_prefix', default=None,
                        help='Prefix for output PNG files')
    args = parser.parse_args()

    prefix = args.output_prefix or f'psha_{args.calc_id}'

    # ── Load data ─────────────────────────────────────────────────────
    print(f'Loading: {args.netcdf}')
    where2save = Path(args.netcdf).parent
    ds_curves = xr.open_dataset(args.netcdf, group='hazard_curves')
    ds_maps   = xr.open_dataset(args.netcdf, group='hazard_maps')
    ds_uhs    = xr.open_dataset(args.netcdf, group='uhs')

    # ── Detect structure ──────────────────────────────────────────────
    imts = detect_imts_from_ds(ds_curves)

    return_periods = sorted(set(
        int(re.search(r'_(\d+)yr$', v).group(1))
        for v in ds_maps.data_vars
        if re.search(r'_(\d+)yr$', v)))

    print(f'IMTs           : {imts}')
    print(f'Return periods : {return_periods} yr')

    # ── Selected sites ────────────────────────────────────────────────
    selected      = read_sites_csv(args.sites)
    matched_sites = find_nearest_sites(selected, ds_curves)

    # ── Plot 1 — Hazard curves ────────────────────────────────────────
    print('\nPlot 1 — Hazard curves...')
    plot_hazard_curves(ds_curves, matched_sites, imts,
                       return_periods, args.calc_id, prefix, where2save)

    # ── Plot 2 — Hazard maps ──────────────────────────────────────────
    print('\nPlot 2 — Hazard maps...')
    plot_hazard_maps(ds_maps, imts, return_periods,
                     matched_sites, args.calc_id, prefix, where2save)

    # ── Plot 3 — UHS ─────────────────────────────────────────────────
    print('\nPlot 3 — UHS...')
    plot_uhs(ds_uhs, matched_sites, args.calc_id, prefix, where2save)

    print('\nAll plots done!')
