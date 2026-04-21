"""
Ryan Pranantyo
EOS, 21 April 2026
Basic plotting to check hazard_curve-mean_<IML>_<ID>.csv resulted

usage:
    python plot__hazard-curves_map.py --input_file <path_to__hazard_curve-mean-<IML>_<ID>.csv>

    OR

    python plot__hazard-curves_map.py --input_file <path_to__quantile_curve-<IMT>-<IML>_<ID>.csv'

    IML = intensity meassure level
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

    IML = name.split('-')[-1].split('_')[0]

    where2save = Path(os.path.join(input_file.parent, 'figures'))
    where2save.mkdir(exist_ok=True)

    print(f'='*60)
    print(f' Basic plotting of hazard curves map mean values ')
    print(f'  Classical PSHA calculation ')
    print(f'  input_file = {name}')
    print(f'  IML        = {IML}')
    print(f'  variables  = {variables}')
    print(f'  domain     = {df.lon.min()}/{df.lon.max()}/{df.lat.min()}/{df.lat.max()}')
    print(f'='*60)

    fig = plt.figure(figsize=(10,10), constrained_layout=True)
    for ii, var in enumerate(variables):
        var_arr = var.split('-')
        poe = float(var_arr[1])
        title = f'$>$ {poe} g'
        ax = fig.add_subplot(ncols, nrows, ii+1)
        ax.set_title(title, loc='left', fontsize=8, pad=0.01)
        val = ax.scatter(df['lon'], df['lat'], c=df[var],
                vmin=0.0, vmax=1.0,
                s=6, alpha=0.8
                )
        ax.tick_params(labelleft=False, labelbottom=False)
        #fig.colorbar(val, shrink=0.95, pad=0.01, orientation='horizontal')
        if ii == (ncols * 2) - 1:
            cax = ax.inset_axes([1.01,0., 0.05, 1.0])
            cb = fig.colorbar(val, cax=cax, orientation='vertical')
            cb.set_label(f'PoE')

    fout = os.path.join(where2save, f'{name}.png')
    fig.savefig(fout, dpi=300)
    plt.close()

    return df

if __name__ == '__main__':
   parser = argparse.ArgumentParser(
           description='Basic quick plotting all variables from hazard_curve-mean-<IMT>_<ID>.csv',
           formatter_class=argparse.RawDescriptionHelpFormatter,
           epilog=__doc__,
           )
   parser.add_argument(
           '--input_file',
           default='/home/ignatius.pranantyo/PSHA/Singapore/Replication_DuPan2020/case02__BackgroundSourceClassicalPSHA/hasil_case02_a/hazard_curve-mean-PGA_3.csv',
           help='Path to hazard_curve-mean-<IMT>_<ID>.csv file',
           )
   parser.add_argument(
           '--calculation_type',
           default=1,
           help='1:Classical PSHA',
           )
   args = parser.parse_args()

   if args.calculation_type == 1:
       df = classic_psha(args)
