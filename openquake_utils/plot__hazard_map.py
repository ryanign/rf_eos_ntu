"""
Ryan Pranantyo
EOS, 21 April 2026
Basic plotting to check hazard_map-mean_<ID>.csv resulted

usage:
    python plot__hazard_map.py --input_file path_to__hazard_map-mean_<ID>.csv

    OR

    python plot__hazard_map.py --input_file path_to__quantile_map-<IMT>_<ID>.csv

    IMT = investigation time
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

def classic_psha(args):
    input_file = Path(args.input_file)
    name = input_file.name[:-4]
    df = pd.read_csv(input_file, skiprows=1)
    variables = list(df.columns[2:])
    if 'depth' in variables:
        variables.remove('depth')
    ncols = int(np.ceil(np.sqrt(len(variables))))
    nrows = ncols

    where2save = Path(os.path.join(input_file.parent, 'figures'))
    where2save.mkdir(exist_ok=True)

    print(f'='*60)
    print(f' Basic plotting of hazard map mean values ')
    print(f'  Classical PSHA calculation ')
    print(f'  input_file = {name}')
    print(f'  variables  = {variables}')
    print(f'  domain     = {df.lon.min()}/{df.lon.max()}/{df.lat.min()}/{df.lat.max()}')
    print(f'='*60)

    fig = plt.figure(figsize=(10,10), constrained_layout=True)
    for ii, var in enumerate(variables):
        var_arr = var.split('-')
        title = f'{var_arr[0]} at {float(var_arr[1])*100}% PoE (g)'
        ax = fig.add_subplot(ncols, nrows, ii+1)
        ax.set_title(title, loc='left')
        val = ax.scatter(df['lon'], df['lat'], c=df[var],
                vmin=df[var].min(), vmax=df[var].max(),
                s=10, alpha=0.8
                )
        fig.colorbar(val, shrink=0.95, pad=0.01, orientation='horizontal')

    fout = os.path.join(where2save, f'{name}.png')
    fig.savefig(fout, dpi=300)
    plt.close()

    return df

if __name__ == '__main__':
   parser = argparse.ArgumentParser(
           description='Basic quick plotting all variables from hazard_map-mean_<ID>.csv',
           formatter_class=argparse.RawDescriptionHelpFormatter,
           epilog=__doc__,
           )
   parser.add_argument(
           '--input_file',
           default='/home/ignatius.pranantyo/PSHA/Singapore/Replication_DuPan2020/case02__BackgroundSourceClassicalPSHA/hasil_case02_a/hazard_map-mean_3.csv',
           help='Path to hazard_map-mean_<ID>.csv file',
           )
   parser.add_argument(
           '--calculation_type',
           default=1,
           help='1:Classical PSHA',
           )
   args = parser.parse_args()

   if args.calculation_type == 1:
       df = classic_psha(args)
