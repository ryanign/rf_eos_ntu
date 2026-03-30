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
import gc
from scipy.stats import norm
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
        mask       = [c for c in cols if f'Mw_{mw}' in c]
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

#-----------------------------------------------------------
# OPTION: 1c - basic damage probability
#-----------------------------------------------------------
def damage_func(x, mu, sigma):
    """
    damage func in normal
    damage(x) = phi((x-mu)/sigma)
        phi   : cummulative standard-normal distribution function
        x     : flood depth
        mu    : mean coefficient
        sigma : standard deviation coef
    """
    return norm.cdf((x-mu) / sigma)

def damage_func_ln(x, mu, sigma):
    """
    damage func in log normal
    damage(x) = phi((np.log(x) - mu) / sigma)
        phi   : cummulative standard-normal distribution function
        x     : flood depth
        mu    : mean coefficient
        sigma : standard deviation coef
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        return norm.cdf((np.log(x) - mu) / sigma)

def damage_classification_syamsidik2023(px):
    """
    Classification of damage states based on Syamsidik et al. 2023, Table 4
    """
    ds = np.zeros(px.shape, dtype=np.int8)
    ds = np.where((px >= 0.01) & (px < 0.10), 1, ds)
    ds = np.where((px >= 0.10) & (px < 0.35), 2, ds)
    ds = np.where((px >= 0.35) & (px < 0.75), 3, ds)
    ds = np.where((px >= 0.75) & (px < 0.90), 4, ds)
    ds = np.where( px >= 0.90               , 5, ds)
    return ds

def histogram_damage_probability(input_file, output_dir, name):
    """
    A basic histogram of building damage ratio following Syamsidik et al. 2023.
        https://doi.org/10.1016/j.ijdrr.2023.103652

    Assumption:
        - All buildings are C1-La classl Reinforced Concrete 1 storey
          mu = 0.888; sigma = 0.839 (Table 3)
        - Classification of damage following Table 4
    """
    os.makedirs(output_dir, exist_ok=True)
    mu = 0.888
    sigma = 0.839

    print(f'='*60)
    print(f' BASIC ESTIMATION OF PROBABILITY DAMAGE LEVELS ')
    print(f' Assumptions:')
    print(f'  - All buildings are reinforced concrete 1 storey (C1-La) of Syamsidik et al. 2023')
    print(f'  - mu    = {mu}')
    print(f'  - sigma = {sigma}')
    print(f'='*60)
    print(f'\n  Loading data ...')
    df = pd.read_parquet(input_file)
    meta_cols = ['FID', 'tile_name']
    cols = [c for c in df.columns if c not in meta_cols]

    # number of buildings affected per-scenario
    # count in batches
    print(f'\n Calculating probability damage...')
    BATCH_SIZE = 25000
    n_buildings   = len(df)
    damage_arr = np.zeros((n_buildings, len(cols)))
    for start in tqdm(range(0, n_buildings, BATCH_SIZE), desc='  Processing:'):
        end = min(start + BATCH_SIZE, n_buildings)
        rows = df.iloc[start:end][cols].values / 100.
        damage_arr[start:end] = damage_func_ln(rows, mu, sigma)
    
    damage_df = df[meta_cols].copy()
    damage_df = pd.concat([
        damage_df.reset_index(drop=True), pd.DataFrame(damage_arr, columns=cols)],
        axis=1,
        )

    # classification
    print(f'\n  Damage states classification ...')
    damage_lev = np.zeros((n_buildings, len(cols)), dtype=np.int8)
    for start in tqdm(range(0, n_buildings, BATCH_SIZE), desc='  Processing:'):
        end = min(start + BATCH_SIZE, n_buildings)
        rows = damage_df.iloc[start:end][cols]
        damage_lev[start:end] = damage_classification_syamsidik2023(rows.values)
    damage_lev_df = df[meta_cols].copy()
    damage_lev_df = pd.concat([
        damage_lev_df.reset_index(drop=True), pd.DataFrame(damage_lev, columns=col)],
        axis=1,
        )

    # count per-Mw
    mw_cols = [c.split('__')[1].split('_')[1] for c in cols]
    mws     = np.unique(mw_cols)
    mw_arr  = np.array(mw_cols)
    BATCH_SIZE = 50000
    records = []

    fig = plt.figure(figsize=(10,20), constrained_layout=True)
    fid = 1
    for mw in tqdm(mws, desc='  Processing'):
        mask = [c for c in cols if f'Mw_{mw}' in c]
        df_m = damage_lev_df[mask]
        for ds in [1,2,3,4,5]:
            n_affected = np.zeros(len(mask), dtype=np.int32)
            for start in range(0, n_buildings, BATCH_SIZE):
                end = min(start + BATCH_SIZE, n_buildings)
                df_r = df_m.iloc[start:end]
                n_affected += (df_r==ds).sum(axis=0)

            records.append({
                'Mw'                 : float(mw),
                'DamageState'        : int(ds),
                'n_scenarios_in_mw'  : len(mask),
                'n_scenarios_per_mw' : int((n_affected >0).sum(axis=0)),
                })
            ax = fig.add_subplot(len(mws), 5, fid)
            ax.set_title(f'Mw {float(mw):.1f} -- DS {ds}', fontsize=6, pad=0.01)
            try:
                bins = np.linspace(1, np.max(n_affected), 5)
                ax.hist(n_affected, bins=bins)

            except ValueError:
                ax.text(1,1,'N.A', ha='center', va='center')
                ax.set_xlim(0.5,1.5)
                ax.set_ylim(0.5,1.5)
            ax.tick_params(labelsize=4, pad=0.01)
            ax.set_xlim(0,np.max(bins)+2)
            if '8.7' in mw:
                ax.set_xlabel('# buildings', fontsize=4)
            if ds == 1:
                ax.set_ylabel('# events', fontsize=4)
            fid += 1
    fout = os.path.join(output_dir, f'n_buildings_affected_DamageState_perMw__{name}.png')
    fig.savefig(fout, dpi=300)
    plt.close()

    result = pd.DataFrame(records).sort_values(['Mw', 'DamageState']).reset_index(drop=True)

    return damage_lev_df

def main(args):
    file_type = args.file_type
    option    = args.option

    if file_type == 1:
        if option == '1a':
            df = histogram_buildings_basic(args.input_file, args.output_dir, args.name)
        elif option == '1b':
            df = histogram_buildings_basic_mw(args.input_file, args.output_dir, args.name)
        elif option == '1c':
            df = histogram_damage_probability(args.input_file, args.output_dir, args.name)

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
            default='1c',
            help='1a: building affected basic; 1b: building affected per-Mw; 1c: damage probability',
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





