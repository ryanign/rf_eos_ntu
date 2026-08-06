"""
Ryan Pranantyo
EOS, 6 August 2026

basic post-processing of Classical PSHA CSV output files, converted to a single NetCDF

useage:
    python post-processing_psha-calculation.py \
            --input_dir ./outputs/ \
            --calc_id 6  \
            --output psha_summary_6.nc

"""
import os
import re
import argparse
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from datetime import datetime

#-------------------------
# detect structure
#-------------------------
def read_oq_csv(path):
    """ Read OpenQuake CSV - skip comment line if present """
    with open(path) as f:
        first = f.readline()
    skip = 1 if first.startswith('#') else 0
    return pd.read_csv(path, skiprows=skip)

def get_site_coords(input_dir, calc_id):
    """ Read site coordinates from site_model_<calc_id>.csv """
    f = Path(input_dir) / f'site_model_{calc_id}.csv'
    df = read_oq_csv(f)
    # Normalise column names
    df.columns = [c.strip() for c in df.columns]
    # rename site_id to custom_site_id if needed
    if 'site_id' in df.columns and 'custom_site_id' not in df.columns:
        df = df.rename(columns={'site_id': 'custom_site_id'})
    return df[['custom_site_id', 'lon', 'lat', 'vs30']].sort_values('custom_site_id')

def detect_imts(input_dir, calc_id):
    """ Detect available IMTs from hazard curve filenames """
    pattern = re.compile(rf'hazard_curve-mean-(.+)_{calc_id}\.csv')
    imts = []
    for f in sorted(os.listdir(input_dir)):
        m = pattern.match(f)
        if m:
            imts.append(m.group(1))
    return imts

def detect_return_periods(input_dir, calc_id):
    """ Detect return periods from hazard_map filenames. """
    pattern = re.compile(rf'hazard_map-mean-(\d+)y_{calc_id}.csv')
    rps = []
    for f in sorted(os.listdir(input_dir)):
        m = pattern.match(f)
        if m:
            rps.append(int(m.group(1)))
    return sorted(rps)

def imt_to_varname(imt):
    """ convert IMT string to valid variable name """
    return (imt.replace('(','').replace(')','').replace('.','p'))

#-------------------------
# building hazard curves
#-------------------------
def build_hazard_curves(input_dir, calc_id, site_coords, imts):
    """
    Build xarray Dataset for hazard curves
    Dimensions: site, iml
    Variables: {IMT}_{stat} where stat = mean, q10, q50, q90
    """
    print(f'\n[hazard curves] Building ...')

    quantiles = {'q10' : '0.1', 'q50' : '0.5', 'q90' : '0.9'}
    stats     = {'mean': 'mean', **quantiles}

    sites = site_coords['custom_site_id'].values
    data_vars = {}

    # coordinates
    data_vars['lon']  = xr.DataArray(site_coords['lon'].values, dims='site')
    data_vars['lat']  = xr.DataArray(site_coords['lat'].values, dims='site')
    data_vars['vs30'] = xr.DataArray(site_coords['vs30'].values, dims='site')

    for imt in imts:
        imt_var = imt_to_varname(imt)
        iml_ref = None

        for stat_name, stat_key in stats.items():
            if stat_key == 'mean':
                fname = Path(input_dir) / f'hazard_curve-mean-{imt}_{calc_id}.csv'
            else:
                fname = Path(input_dir) / f'quantile_curve-{stat_key}-{imt}_{calc_id}.csv'

            if not fname.exists():
                print(f'  WARNING: {fname.name} not found, skipping ...')
                continue

            df = read_oq_csv(fname)
            if 'custom_site_id' in df.columns:
                df = df.sort_values('custom_site_id').reset_index(drop=True)

            # Extract IML columns (poe-*)
            poe_cols = [c for c in df.columns if c.startswith('poe-')]
            iml_vals = np.array([float(c.replace('poe-','')) for c in poe_cols])

            if iml_ref is None:
                iml_ref = iml_vals
                data_vars[f'{imt_var}_iml'] = xr.DataArray(
                        iml_vals, dims=f'{imt_var}_iml',
                        attrs={'imt': imt, 'units': 'g'})

            #PoE values: shape (n_sites, n_iml)
            poe_vals = df[poe_cols].values

            varname = f'{imt_var}_{stat_name}'
            data_vars[varname] = xr.DataArray(
                    poe_vals,
                    dims=['site', f'{imt_var}_iml'],
                    attrs={'imt': imt, 'statistic':stat_name, 'units': 'probability'})

            print(f'  {imt} {stat_name} : {poe_vals.shape}')

    ds = xr.Dataset(
            data_vars,
            coords={'site': sites},
            attrs={
                'group' : 'hazard_curves',
                'description' : 'Hazard curves - PoE vs IML per site',
                'imts' : ', '.join(imts),
                'statistics' : ['mean', 'q10', 'q50', 'q90'],
                'created' : datetime.now().isoformat(),
                })

    return ds

#-------------------------
# build hazard maps
#-------------------------
def build_hazard_maps(input_dir, calc_id, site_coords, imts, return_periods):
    """
    Build xarray Dataset for hazard maps
    """
    print(f'\n[hazard_maps] Building ...')

    quantiles = {'q10': '0.1', 'q50': '0.5', 'q90': '0.9'}
    stats     = {'mean': 'mean', **quantiles}

    sites     = site_coords['custom_site_id'].values
    data_vars = {}

    # Coordinates
    data_vars['lon']  = xr.DataArray(site_coords['lon'].values,  dims='site')
    data_vars['lat']  = xr.DataArray(site_coords['lat'].values,  dims='site')
    data_vars['vs30'] = xr.DataArray(site_coords['vs30'].values, dims='site')

    for stat_name, stat_key in stats.items():
        for rp in return_periods:
            if stat_key == 'mean':
                fname = Path(input_dir) / f'hazard_map-mean-{rp}y_{calc_id}.csv'
            else:
                fname = Path(input_dir) / f'quantile_map-{stat_key}-{rp}y_{calc_id}.csv'

            if not fname.exists():
                print(f'  WARNING: {fname.name} not found, skipping')
                continue

            df = read_oq_csv(fname)
            if 'custom_site_id' in df.columns:
                df = df.sort_values('custom_site_id').reset_index(drop=True)

            for imt in imts:
                imt_var = imt_to_varname(imt)

                # find column matching IMT
                # hazard_map columns: lon, lat, PGA, SA(1.0), ...
                imt_col = None
                for col in df.columns:
                    if col == imt or col.startswith(f'{imt}-'):
                        imt_col = col
                        break

                if imt_col is None:
                    print(f'  WARNING: {imt} column not found in {fname.name}')
                    continue

                vals = df[imt_col].values
                varname = f'{imt_var}_{stat_name}_{rp}yr'
                data_vars[varname] = xr.DataArray(
                        vals, dims='site',
                        attrs={
                            'imt' : imt,
                            'statistic' : stat_name,
                            'return_period': rp,
                            'units' : 'g',
                            })

                print(f'  {imt} {stat_name} {rp}yr: {vals.shape}')

    ds = xr.Dataset(
        data_vars,
        coords={'site': sites},
        attrs={
            'group'          : 'hazard_maps',
            'description'    : 'Hazard maps — IML at fixed return period per site',
            'imts'           : ', '.join(imts),
            'return_periods' : str(return_periods),
            'statistics'     : 'mean, q10, q50, q90',
            'created'        : datetime.now().isoformat(),
        })

    return ds

#-------------------------
# build UHS
#-------------------------
def build_uhs(input_dir, calc_id, site_coords):
    """
    Build xarray Dataset for Uniform Hazard Spectra
    """
    print(f'\n[uhs] Building ...')

    quantiles = {'q10': '0.1', 'q50': '0.5', 'q90': '0.9'}
    stats     = {'mean': 'mean', **quantiles}

    sites     = site_coords['custom_site_id'].values
    data_vars = {}

    # Coordinates
    data_vars['lon']  = xr.DataArray(site_coords['lon'].values,  dims='site')
    data_vars['lat']  = xr.DataArray(site_coords['lat'].values,  dims='site')
    data_vars['vs30'] = xr.DataArray(site_coords['vs30'].values, dims='site')

    period_ref = None

    for stat_name, stat_key in stats.items():
        if stat_key == 'mean':
            fname = Path(input_dir) / f'hazard_uhs-mean_{calc_id}.csv'
        else:
            fname = Path(input_dir) / f'quantile_uhs-{stat_key}_{calc_id}.csv'

        if not fname.exists():
            print(f'  WARNING: {fname.name} not found, skipping')
            continue

        df = read_oq_csv(fname)
        if 'custom_site_id' in df.columns:
            df = df.sort_values('custom_site_id').reset_index(drop=True)

        # UHS columns format: {poe}~{IMT}  e.g. 0.100000~PGA, 0.020000~SA(1.0)
        uhs_cols = [c for c in df.columns
                    if '~' in c and c not in ('lon','lat','custom_site_id')]

        # Detect unique PoEs
        poes = sorted(set(c.split('~')[0] for c in uhs_cols))

        for poe_str in poes:
            # poe label: 0.100000 -> 10pct, 0.020000 -> 2pct
            poe_val = float(poe_str)
            poe_label = f'{round(poe_val*100):g}pct'

            # Columns for this PoE
            cols_this_poe = [c for c in uhs_cols if c.startswith(poe_str)]

            # Extract periods from IMT names
            periods = []
            for col in cols_this_poe:
                imt = col.split('~')[1]
                if imt == 'PGA':
                    periods.append(0.0)
                else:
                    T = float(re.search(r'\(([\d.]+)\)', imt).group(1))
                    periods.append(T)

            sort_idx = np.argsort(periods)
            periods  = np.array(periods)[sort_idx]
            cols_sorted = [cols_this_poe[i] for i in sort_idx]

            if period_ref is None:
                period_ref = periods
                data_vars['period'] = xr.DataArray(
                    periods, dims='period',
                    attrs={'units': 'seconds', 'description': 'Spectral period (0=PGA)'})

            # Sa values: shape (n_sites, n_periods)
            sa_vals = df[cols_sorted].values

            varname = f'poe_{poe_label}_{stat_name}'
            data_vars[varname] = xr.DataArray(
                sa_vals,
                dims=['site', 'period'],
                attrs={
                    'poe'       : poe_val,
                    'statistic' : stat_name,
                    'units'     : 'g',
                })
            print(f'  PoE={poe_label} {stat_name}: {sa_vals.shape}')

    ds = xr.Dataset(
        data_vars,
        coords={'site': sites},
        attrs={
            'group'      : 'uhs',
            'description': 'Uniform Hazard Spectra — Sa at fixed PoE per site',
            'statistics' : 'mean, q10, q50, q90',
            'created'    : datetime.now().isoformat(),
        })

    return ds

#-------------------------
# MAIN
#-------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Build NetCDF from OpenQuake Classical PSHA CSV outputs, based on ver 3.25.1.')
    parser.add_argument(
            '--input_dir', 
            default='/home/ignatius.pranantyo/PSHA/Singapore/Simple_Deterministic/exercise__psha-based__SumatraSubduction__multiGMPEs/outputs',
            help='Directory containing OQ CSV exported files',
            )
    parser.add_argument(
            '--calc_id',
            type=int,
            default=6,
            help='OpenQuake calculation ID',
            )
    args = parser.parse_args()

    input_dir = args.input_dir
    calc_id   = args.calc_id
    output_nc = os.path.join(args.input_dir, f'psha_summary_{calc_id}.nc')

    print(f'='*60)
    print(f'Converting CSV files to NetCDF file ...')
    print(f'  input folder : {input_dir}')
    print(f'  calc_id      : {calc_id}'  )
    print(f'  output file  : {output_nc}')
    print(f'-'*60)

    # -- detect structure --
    site_coords    = get_site_coords(input_dir, calc_id)
    imts           = detect_imts(input_dir, calc_id)
    return_periods = detect_return_periods(input_dir, calc_id)
    
    print(f'\n')
    print(f'Sites : {len(site_coords)}')
    print(f'IMTs  : {imts}')
    print(f'RPs   : {return_periods} yr')

    # -- build groups
    ds_curves = build_hazard_curves(input_dir, calc_id, site_coords, imts)
    ds_maps   = build_hazard_maps(input_dir, calc_id, site_coords, imts, return_periods)
    ds_uhs    = build_uhs(input_dir, calc_id, site_coords)

    # -- save netcdf
    print(f'\n{"="*60}')
    print(f'Saving: {output_nc}')

    encoding = lambda ds: {v: {'zlib': True, 'complevel': 4} for v in ds.data_vars}

    groups = {
        'hazard_curves': ds_curves,
        'hazard_maps'  : ds_maps,
        'uhs'          : ds_uhs,
    }

    for i, (group, ds) in enumerate(groups.items()):
        mode = 'w' if i == 0 else 'a'
        ds.to_netcdf(output_nc, group=group, mode=mode, encoding=encoding(ds))
        print(f'  Written: /{group}/')

    print(f'\nDone! {output_nc}')
    print(f'  Groups: {list(groups.keys())}')

    # -- quick summary
    print(f'\nVariables per group:')
    for group, ds in groups.items():
        print(f'  /{group}/: {len(ds.data_vars)} variables')

    print(f'\n')
    print(f'-'*60)
    print(f' DONE : ..._{calc_id}.csv files can be removed if no longer needed')

    # -- file size summary
    print(f'  File size comparison:')
    
    # total csv size
    total_csv = 0
    csv_files = list(Path(input_dir).glob(f'*_{calc_id}.csv'))
    for f in csv_files:
        total_csv += f.stat().st_size
        print(f'   {f.name:55s} {f.stat().st_size/1024:8.1f} KB')

    # netcdf size
    nc_size = Path(output_nc).stat().st_size

    print(f'\n')
    print(f'   {"Total CSV":55s} {total_csv/1024/1024:8.2f} MB')
    print(f'   {"NetCDF (compressed)":55s} {nc_size/1024/1024:8.2f} MB')
    print(f'   {"Compression ratio":55s} {total_csv/nc_size:8.2f}x')



