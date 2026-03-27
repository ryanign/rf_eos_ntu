"""
Ryan Pranantyo
EOS, 25 March 2026

Description:
    A basic script to do quick analysis of:
        1)
        2)
        3)

Outputs:

Useage:
    python basic_analysis.py \\
        --input_file \\
        --file_type \\
        --option \\
        --output_dir \\
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import pyarrow.parquet as pq
from tqdm import tqdm

#-----------------------------------------------------------
# OPTION: 1a - basic
#-----------------------------------------------------------
def histogram_buildings_basic(input_file, output_dir, name):
    """
    It just to do a quick analysis on the number of buildings affected,
    flood depth intensity is not considered here
    """
    print(f'=== PLOTTING HISTOGRAM OF BUILDINGS AFFECTED ===')
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_parquet(input_file)
    cols = [c for c in df.columns if c not in ['FID', 'tile_name']]

    # number of buildings affected per-scenario
    # count in batches
    BATCH_SIZE = 50000
    n_buildings   = len(df)
    _n_scenarios_ = np.zeros(len(cols)).astype(int)
    for start in tqdm(range(0, n_buildings, BATCH_SIZE), desc='  Processing'):
        end = min(start + BATCH_SIZE, n_buildings)
        df_rows = df.iloc[start:end]
        _n_scenarios_ += (df_rows[cols]>0).sum(axis=0).values
    
    # simple histogram
    bins = np.linspace(0, 5000, 11)
    bins[0] = 1
    #bins = np.append(bins, np.max(_n_scenarios_)+5)
    hist,_ = np.histogram(_n_scenarios_, bins=bins)
    n_scenarios_ = (_n_scenarios_ != 0).sum()              ### number of scenarios affected buildings
    n_scenarios_gt5000 = (_n_scenarios_ > bins.max()).sum()

    #print(hist)
    fig = plt.figure(constrained_layout=True)
    ax = fig.add_subplot(1,1,1)
    ax.hist(_n_scenarios_, bins=bins)
    ax.set_xlabel('Num of buildings')
    ax.set_ylabel('Num of scenarios')
    
    text = f'Num of events : {len(cols):,}\nNum of events affected buildings : {n_scenarios_:,}\nNum of events affected more than 5,000 buildings : {n_scenarios_gt5000:,}'
    ax.text(4800, hist.max()-10, text, ha='right', va='top')

    fout = os.path.join(output_dir, f'n_buildings_affected__{name}.png')
    fig.savefig(fout, dpi=300)
    plt.close()

    print(f'  n_events                            : {len(cols):,}')
    print(f'  n_events caused inundation          : {n_scenarios_:,}')
    print(f'  n_events affected > 5,000 buildings : {n_scenarios_gt5000:,}')
    print(f'  histogram saved at {fout}')
    
    return df

#-----------------------------------------------------------
# OPTION: 1b - basic per Mw
#-----------------------------------------------------------
def histogram_buildings_basic_mw(input_file, output_dir, name):
    print(f'=== PLOTTING HISTOGRAM OF BUILDINGS AFFECTED PER MAGNITUDE ===')
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_parquet(input_file)
    cols = [c for c in df.columns if c not in ['FID', 'tile_name']]

    mw_cols = [c.split('__')[1].split('_')[1] for c in cols]
    mws     = np.unique(mw_cols)
    mw_arr  = np.array(mw_cols)

    # number of buildings affected per-mw
    # count per-mw
    BATCH_SIZE = 50000
    n_buildings = len(df)
    records = []
    for mw in tqdm(mws, desc='  Processing'):
        mask       = [c for c in cols if mw in c]
        df_m       = df[mask].values  # numpy array, lebih cepat
        n_affected = np.zeros(len(mask), dtype=np.int32)
    
        for start in range(0, n_buildings, BATCH_SIZE):
            end          = min(start + BATCH_SIZE, n_buildings)
            n_affected  += (df_m[start:end] > 0).sum(axis=0)
    
        records.append({
            'Mw'                : float(mw),
            'n_scenarios_per_mw': int((n_affected > 0).sum()),
            'n_scenarios_in_mw' : len(mask),
        })
    
    result = pd.DataFrame(records).sort_values('Mw').reset_index(drop=True)
    result['ratio_pct'] = result['n_scenarios_per_mw'] / result['n_scenarios_in_mw'] * 100.

    fig = plt.figure(constrained_layout=True)
    ax = fig.add_subplot(1,1,1)
    ax.set_title('Buildings affected per-Mw')
    ax.vlines(result['Mw'], ymin=0, ymax=result['n_scenarios_in_mw'], zorder=0, 
            linewidths=2,colors='black', label='num of events in catalogue')
    ax.vlines(result['Mw'], ymin=0, ymax=result['n_scenarios_per_mw'], zorder=1,
            linewidths=2,colors='red', label='num of events generated inundation')
    ax.legend(prop={'size':8})
    
    ax1 = ax.twinx()
    ax1.plot(result['Mw'], result['ratio_pct'], c='green', linestyle='--')

    ax.set_ylim(0,result['n_scenarios_in_mw'].max()+5)
    ax.set_xticks(result['Mw'])
    ax.set_ylabel('Num of scenarios')
    ax.set_xlabel('Mw')

    ax1.set_ylabel('% red/black')
    ax1.tick_params(color='green', labelcolor='green')
    ax1.yaxis.set_label_position('right')
    ax1.grid(linestyle=':', alpha=0.3, c='gray')
    ax1.set_ylim(0,100)

    fout = os.path.join(output_dir, f'n_buildings_affected_perMw__{name}.png')
    fig.savefig(fout, dpi=300)
    plt.close()

    print(f'  histogram saved at {fout}')

    return df

def main(args):
    file_type = args.file_type
    option    = args.option

    if file_type == 1:
        if option == '1a':
            df = histogram_buildings_basic(args.input_file, args.output_dir, args.name)
        elif option == '1b':
            df = histogram_buildings_basic_mw(args.input_file, args.output_dir, args.name)

    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '--input_file',
            default='/home/ignatius.pranantyo/Tsunamis/PTRA_SouthernJava//flood_matrix/flood_depth_matrix__Tile_2-10__NLSWE.parquet',
            help='Path to an input parquet file',
            )
    parser.add_argument(
            '--file_type',
            type=int,
            default=1,
            help='1: Flood depth matrix; 2: bla bla ;3:',
            )
    parser.add_argument(
            '--option',
            default='1b',
            help='1a: building affected basic; 1b: building affected per-Mw',
            )
    parser.add_argument(
            '--output_dir',
            type=str,
            default='/home/ignatius.pranantyo/Tsunamis/PTRA_SouthernJava/figs_analysis',
            help='Directory to save output files',
            )
    parser.add_argument(
            '--name',
            type=str,
            default='Tile_2-10__NLSWE',
            help='For file naming purpose',
            )
    args = parser.parse_args()


    df = main(args)





