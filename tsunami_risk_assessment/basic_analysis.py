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
from tqdm import tqdm

def histogram_buildings(input_file, output_dir):
    print(f'=== PLOTTING HISTOGRAM OF BUILDINGS AFFECTED ===')
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_parquet(input_file)
    cols = df.columns[[c for c in df.columns!='FID']]

    # number of buildings affected per-scenario
    ##n_buildings = (df[cols]>0).sum(axis=0)
    # count in batches
    BATCH_SIZE = 50000
    n_buildings   = len(df)
    _n_scenarios_ = np.zeros(len(cols)).astype(int)
    for start in tqdm(range(0, n_buildings, BATCH_SIZE), desc='  Processing'):
        end = min(start + BATCH_SIZE, n_buildings)
        df_rows = df.iloc[start:end]
        _n_scenarios_ += (df_rows[cols]>0).sum(axis=0)
    
    # simple histogram
    bins = np.arange(10, 50000, 250)
    #hist,_ = np.histogram(n_buildings, bins=bins)

    #print(hist)
    fig = plt.figure(constrained_layout=True)
    ax = fig.add_subplot(1,1,1)
    ax.hist(_n_scenarios_, bins=bins)
    ax.set_xlabel('Num of buildings')
    ax.set_ylabel('Num of scenarios')
    
    fout = os.path.join(output_dir, f'n_buildings_affected.png')
    fig.savefig(fout, dpi=300)
    plt.close()

    print(f'  n_buildings : {n_buidlings:,}')
    print(f'  histogram saved at {fout}')
    
    return df

def main(args):
    file_type = args.file_type
    option    = args.option

    if file_type == 1:
        if option == '1a':
            df = histogram_buildings(args.input_file, args.output_dir)

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
            help='1: Building scenario count; 2: bla bla ;3:',
            )
    parser.add_argument(
            '--option',
            default='1a',
            help='1a: basic histogram of building impacted',
            )
    parser.add_argument(
            '--output_dir',
            default='/home/ignatius.pranantyo/Tsunamis/PTRA_SouthernJava/figs_analysis',
            help='Directory to save output files',
            )
    args = parser.parse_args()


    df = main(args)





