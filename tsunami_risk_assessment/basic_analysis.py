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


def histogram_buildings(input_file):
    print(f'=== PLOTTING HISTOGRAM OF BUILDINGS AFFECTED ===')
    df = pd.read_parquet(input_file)
    cols = df.columns[[c for c in df.columns!='FID']]

    # number of buildings affected per-scenario
    n_buildings = (df[cols]>0).sum(axis=0)
    
    # simple histogram
    bins = np.arange(10, 50000, 250)
    hist,_ = np.histogram(n_buildings, bins=bins)

    return df

def main(args):
    file_type = args.file_type
    option    = args.option

    if file_type == 1:
        if option == '1a':
            df = histogram_buildings(args.input_file)

    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '--input_file',
            default='/home/ryan/SJava_PTRA_analysis/flood_matrix/flood_depth_matrix__Tile_2-10__NLSWE.parquet',
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
            default='/home/ryan/SJava_PTRA_analysis/figs_analysis',
            help='Directory to save output files',
            )
    args = parser.parse_args()


    df = main(args)





