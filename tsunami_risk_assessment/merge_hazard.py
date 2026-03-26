"""
Ryan Pranantyo
EOS, 25 March 2026

Description:
    Merge hazard__<tile_name>.parquet files into one big matrix.

Output:

Useage:


"""
import os
import time
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm

#---------------------------
# LOADING TILES
#---------------------------
def load_tiles(input_files):
    print(f'\n=== LOADING FILES ===')

    dfs = []
    for f in tqdm(input_files, desc='  Loading'):
        tile = f.split('/')[-1].split('.')[0]
        df = pd.read_parquet(f)
        df.insert(1, 'tile_name', tile)
        dfs.append(df)
        print(f'  Loaded {tile} : {len(df):,} buildings, {len(df.columns)} columns')

    combined = pd.concat(dfs, ignore_index=True)

    n_total     = len(combined)
    n_unique    = combined['FID'].nunique()
    n_duplicate = n_total - n_unique
    print(f'\n  Combined shape  : {combined.shape}')
    print(f'  Unique FIDs     : {n_unique:,}')
    print(f'  Duplicated FIDs : {n_duplicate:,} (border buildings)')

    return combined

#---------------------------
# RESOLVING DUPLICATE FIDs
#---------------------------
def resolve_duplicates(combined):
    print(f'\n=== RESOLVING DUPLICATE FIDs ===')
    lambda_cols = [c for c in combined.columns if c.startswith('lambda_')]

    # use sum of lambda values as tiebreaker
    combined['_lambda_sum'] = combined[lambda_cols].fillna(0).sum(axis=1)
    combined_sorted = combined.sort_values(
            ['FID', '_lambda_sum'], 
            ascending=[True, False],
            )
    merged = combined_sorted.drop_duplicates(subset='FID', keep='first')
    merged = merged.drop(columns=['_lambda_sum']).reset_index(drop=True)

    assert merged['FID'].nunique() == len(merged), 'Duplicate FIDs remain after merge!'

    print(f'  Resolved shape : {merged.shape}')
    print(f'  Unique FIDs    : {merged["FID"].nunique():,}')

    # summary
    print(f'\n  Buildings per tile (after resolving duplicates):')
    tile_counts = merged['tile_name'].value_counts().sort_index()
    for tile, count in tile_counts.items():
        print(f'  {tile} : {count:,}')

    return merged

#---------------------------
# SAVING
#---------------------------
def save_output(merged, name, output_dir):
    print(f'\n=== SAVING ===')
    os.mkdirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f'hazard__all_tiles__{name}.parquet')
    merged.to_parquet(output_path, index=False, compression='snappy')
    
    size_mb = os.path.getsize(output_path) / 1024**2
    print(f'  Saved : {output_path} ({size_mb:.1f} MB)')
    print(f'  Shape : {merged.shape}')


#---------------------------
# MAIN
#---------------------------
def main(args):
    print(f'='*60)
    print(f' MERGING COMPUTED HAZARD DATAFRAMES ')
    print(f'='*60)
    print(f'  Input files : {args.input_files}')
    print(f'  Output dir  : {args.output_dir}')
    print(f'='*60)

    t0 = time.time()

    # Step 1: combine all input files
    combined = load_tiles(args.input_files)
    
    # Step 2: resolving duplicate FIDs
    merged = resolve_duplicates(combined)
    
    # Step 3: saving
    save_output(merged, args.name, args.output_dir)

    total = time.time() - t0
    print(f'\nTotal elapsed : {total:.1f}s ({total/60:.1f} min)')
    print(f'Done!')
    print(f'='*60)

    return merged

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '--input_files',
            type=str,
            nargs='+',
            default=[
                '/home/ryan/SJava_PTRA_analysis/flood_hazard/hazard__Tile_5__NLSWE.parquet',
                '/home/ryan/SJava_PTRA_analysis/flood_hazard/hazard__Tile_6__NLSWE.parquet',
                ],
            help='A list of hazard__<tile_name>.parquet files',
            )
    parser.add_argument(
            '--output_dir',
            default='/home/ryan/SJava_PTRA_analysis/flood_hazard',
            help='Path to save output',
            )
    parser.add_argument(
            '--name',
            default='NLSWE',
            help='For naming purpose only',
            )
    args = parser.parse_args()

    merged = main(args)
