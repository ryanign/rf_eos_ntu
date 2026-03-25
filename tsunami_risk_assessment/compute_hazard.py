"""
Ryan Pranantyo
EOS, 24 March 2026

Description:
    Compute flood depth hazard statistics from a flood depth matrix prepared using
    assign_values_to_buildings.py.
    It generates an intermediate output of building_scenario_count__<tile_name>.parquet. 
    Remove this file first if you need to analyse different flood depth thresholds.

Outputs:
    <inside --output_dir>
    1. building_scenario_count__<tile_name>.parquet
        FID | n_scenarios_h10 | n_scenarios_h20 | ... | n_scenarios_h3000 
        number of buildings flooded above 10 cm, 20 cm, ..., 3000 cm.

    2. hazard__<tile_name>.parquet
        FID | lambda_h10 | ... | var_h10 | ... | sigma_h10 | ... | H_500yr |
        - lambda_h10 = annual exceedance rate above 10 cm.
        - var_h10 = variance of lambda above 10 cm.
        - sigma_h10 = standard deviation of lambda above 10 cm.
        - H_500yr = flood depth hazard at 1/500 annual exceedance rate.

Usage:
    python compute_hazard.py \\
        --flood_matrix /path/to/flood_depth_matrix__<tile_name>.parquet \\ 
        --output_dir /path/to/save/output/files \\
        --thresholds 10 3010 10 \\ <start> <stop> <end> -> will generate a list of array
        --nyears 10000 \\ length of catalogue to consider in year
        --return_periods 500 2500 5000 10000 \\ list of return periods to consider for calculating hazard
        --n_jobs -2 \\ number of cpus to use

Note:
    Assisted by Claude (Anthropic, claude.ai) to make the script efficient and easier to follow
"""

import os
import time
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
from scipy.interpolate import interp1d
from pathlib import Path

#---------------------------------------------------------------------------
# Helper
#---------------------------------------------------------------------------
def plot_hazard_curves(lambda_df, sigma_df, df_count, thresholds, output_dir, name, nyears, n_top=15, n_random=5, seed=42):
    """
    Plot hazard curves with uncertainty for selected buildings.
        - n_top    : buildings with highest max lambda
        - n_random : random buildings with lambda > 0 (excluding top)
    Output:
        hazard_curves__<name>.png
    """
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    np.random.seed(seed)

    thresholds_arr = np.array(thresholds)

    # ── Pilih buildings ──────────────────────────────────────
    max_lambda  = lambda_df.max(axis=1)
    top_idx     = max_lambda.nlargest(n_top).index
    flooded_idx = max_lambda[max_lambda > 0].index.difference(top_idx)

    n_random_actual = min(n_random, len(flooded_idx))
    random_idx      = np.random.choice(flooded_idx, size=n_random_actual, replace=False)

    print(f"  Top {n_top} FIDs   : {df_count.loc[top_idx, 'FID'].values}")
    print(f"  Random {n_random_actual} FIDs : {df_count.loc[random_idx, 'FID'].values}")

    # ── Plot ─────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(11, 7))

    colors_top    = cm.Reds(np.linspace(0.4, 0.9, n_top))
    colors_random = cm.Blues(np.linspace(0.4, 0.9, n_random_actual))

    for i, idx in enumerate(top_idx):
        fid  = df_count.loc[idx, "FID"]
        lam  = lambda_df.loc[idx].values.astype(float)
        sig  = sigma_df.loc[idx].values.astype(float)
        mask = ~np.isnan(lam)
        ax.semilogy(thresholds_arr[mask], lam[mask],
                    color=colors_top[i], lw=1.5,
                    label=f"top (FID {fid})" if i == 0 else "_")
        #ax.fill_between(thresholds_arr[mask],
        #                np.maximum(lam[mask] - 1.96 * sig[mask], 1/nyears),
        #                lam[mask] + 1.96 * sig[mask],
        #                color=colors_top[i], alpha=0.15)

    for i, idx in enumerate(random_idx):
        fid  = df_count.loc[idx, "FID"]
        lam  = lambda_df.loc[idx].values.astype(float)
        sig  = sigma_df.loc[idx].values.astype(float)
        mask = ~np.isnan(lam)
        ax.semilogy(thresholds_arr[mask], lam[mask],
                    color=colors_random[i], lw=1.5, ls="--",
                    label=f"random (FID {fid})" if i == 0 else "_")
        #ax.fill_between(thresholds_arr[mask],
        #                np.maximum(lam[mask] - 1.96 * sig[mask], 1/nyears),
        #                lam[mask] + 1.96 * sig[mask],
        #                color=colors_random[i], alpha=0.15)

    # Reference return period lines
    for rp, ls in zip([100, 500, 1000, 2500], [":", "--", "-.", ":"]):
        ax.axhline(1/rp, color="gray", ls=ls, lw=0.8, alpha=0.6, label=f"1/{rp} yr")

    ax.set_ylim(1/10000,1/10)
    ax.set_xlabel("Flood Depth (cm)", fontsize=12)
    ax.set_ylabel("Annual Exceedance Rate λ", fontsize=12)
    ax.set_title(f"Flood Hazard Curves with Uncertainty — {name}\n"
                 f"Red = top {n_top} | Blue dashed = {n_random_actual} random flooded",
                 fontsize=12)
    ax.legend(loc="upper right", fontsize=7, ncol=2)
    ax.grid(True, which="both", alpha=0.3)
    plt.tight_layout()

    output_path = os.path.join(output_dir, f"hazard_curves__{name}.png")
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  Saved : {output_path}")

#---------------------------------------------------------------------------
# LOAD OR COMPUTE SCENARIO COUNT
#---------------------------------------------------------------------------
def _count_threshold(values, t):
    """Count number of scenarios exceeding threshold t per building."""
    return (values >= t).sum(axis=1)

def load_or_compute_scenario_count(flood_matrix_path, thresholds, output_dir, name, n_jobs):
    """
    Load building_scenario_count if it exists, otherwise compute from flood matrix.
    Saves intermediate result for future reuse.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    count_path = os.path.join(output_dir, f'building_scenario_count__{name}.parquet')

    if os.path.exists(count_path):
        print(f'  Found existing building_scenario_count, loading ...')
        df_count = pd.read_parquet(count_path)
        # load FIDs and n_scenarios from flood_matrix
        df_meta = pd.read_parquet(flood_matrix_path, columns=['FID'])
        n_scenarios = len(pd.read_parquet(flood_matrix_path).columns) - 1 # exclude FID
        print(f'  Loaded : {len(df_count):,} buildings, {len(thresholds)} thresholds')
        return df_count, n_scenarios

    print(f'  Loading flood depth matrix ...')

    df_fid      = pd.read_parquet(flood_matrix_path, columns = ['FID'])
    all_cols    = pd.read_parquet(flood_matrix_path).columns.tolist()
    event_cols  = [c for c in all_cols if c != 'FID']
    n_scenarios = len(event_cols)
    n_buildings = len(df_fid)

    print(f'  n_buildings : {n_buildings:,}')
    print(f'  n_scenarios : {n_scenarios:,}')

    counts = {f'n_scenarios_h{int(t)}': np.zeros(n_buildings, dtype=np.int32) for t in thresholds}

    # process in column batches
    COL_BATCH = 200
    for start in tqdm(range(0, n_scenarios, COL_BATCH), desc='  Processing batches'):
        end   = min(start + COL_BATCH, n_scenarios)
        cols  = event_cols[start:end]
        batch = pd.read_parquet(flood_matrix_path, columns=cols).values
        
        for t in thresholds:
            counts[f'n_scenarios_h{int(t)}'] += (batch >= t).sum(axis=1)
        del batch

    df_count = pd.DataFrame(counts)
    df_count.insert(0, 'FID', df_fid['FID'].values)

    df_count.to_parquet(count_path, index=False, compression='snappy')
    size_mb = os.path.getsize(count_path) / 1024**2
    print(f'  Saved intermediate : {count_path} ({size_mb:.1f} MB)')

    return df_count, n_scenarios

#---------------------------------------------------------------------------
# ANNUAL EXCEEDANCE RATE
#---------------------------------------------------------------------------
def compute_annual_exceedance_rate(df_count, thresholds, n_scenarios, nyears):
    """
    lambda(H>=h) = 1/nyears * sum(P(h))
    var(H>=h)    = 1/nyears^2 * sum(P(h) * [1 - P(h)])
    sigma(H>=h)  = sqrt(var)

    where P(h) = n_scenarios_exceeded / n_scenarios  per building
    """
    print(f'\n=== ANNUAL EXCEEDANCE RATE ===')
    annual_rate = n_scenarios / nyears
    lambda_dict = {}
    var_dict    = {}
    sigma_dict  = {}

    for t in tqdm(thresholds, desc='  Computing'):
        col = f'n_scenarios_h{int(t)}'
        Ph  = df_count[col].values / n_scenarios               # P(h) per building

        lambda_dict[f'lambda_h{int(t)}'] = Ph * annual_rate
        var = Ph * (1 - Ph) * (annual_rate ** 2)
        var_dict[f'var_h{int(t)}']     = var
        sigma_dict[f'sigma_f{int(t)}'] = np.sqrt(var)

    lambda_df = pd.DataFrame(lambda_dict)
    var_df    = pd.DataFrame(var_dict)
    sigma_df  = pd.DataFrame(sigma_dict)

    # conver 0 -> NaN (buildings never flooded)
    lambda_df = lambda_df.where(lambda_df != 0, np.nan)
    var_df    = var_df.where(var_df != 0, np.nan)
    sigma_df  = sigma_df.where(sigma_df != 0, np.nan)

    print(f'  Buildings with lambda > 0 : {lambda_df.iloc[:,0].notna().sum():,}')

    return lambda_df, var_df, sigma_df

#---------------------------------------------------------------------------
# COMPUTE FLOOD DEPTH AT RETURN PERIOD
#---------------------------------------------------------------------------
def _h_at_rate_single(lambda_vals, thresholds, yr_rate):
    """
    Interpolate flood depth at a given annual exceedance rate for one building.
    lambda_vals : array of lambda values per thresholds (may contain NaN)
    thresholds  : array of threshold values (cm)
    yr_rate     : target annual exceedance rate (1/return_period)
    """
    mask    = ~np.isnan(lambda_vals)
    x_clean = lambda_vals[mask]
    y_clean = thresholds[mask]

    if len(x_clean) == 0:
        return np.nan
    if yr_rate < x_clean.min():
        return y_clean[-1]
    if yr_rate > x_clean.max():
        return y_clean[0]

    func = interp1d(x_clean, y_clean, kind='linear')
    return float(func(yr_rate))

def _compute_h_batch(lambda_array, thresholds, yr_rate, start, end):
    """
    Process a batch of buildings for one return period
    """
    return np.array([
        _h_at_rate_single(lambda_array[i], thresholds, yr_rate)
        for i in range(start, end)
        ])

def compute_h_at_return_periods(lambda_df, thresholds, return_periods, n_jobs):
    """
    For each building and each return period, interpolate flood depth.
    Parallelized per return period
    """
    print(f'\n=== FLOOD DEPTH AT RETURN PERIODS ===')
    lambda_array = lambda_df.values
    thresholds   = np.array(thresholds)
    n_buildings  = len(lambda_df)
    BATCH_SIZE   = 10000                  # currently fixed

    rp_dict = {}
    for rp in tqdm(return_periods, desc='  Return periods'):
        yr_rate = 1.0 / rp
        # Parallel per batch of buildings
        batches = range(0, n_buildings, BATCH_SIZE)
        results = Parallel(n_jobs=n_jobs)(
                delayed(_compute_h_batch)(lambda_array, thresholds, yr_rate, start, min(start+BATCH_SIZE, n_buildings))
                for start in batches
                )
        rp_dict[f'H_{int(rp)}yr'] = np.concatenate(results)

    rp_df = pd.DataFrame(rp_dict)
    print(f'  Done. Shape : {rp_df.shape}')
    return rp_df

#---------------------------------------------------------------------------
# SAVE OUTPUT
#---------------------------------------------------------------------------
def save_output(df_count, lambda_df, var_df, sigma_df, rp_df, output_dir, name):
    print(f'\n=== SAVING ===')
    os.makedirs(output_dir, exist_ok=True)
    # combine all into one flat DataFrame
    hazard = pd.concat(
            [df_count[['FID']], lambda_df, var_df, sigma_df, rp_df],
            axis=1,
            )
    output_path = os.path.join(output_dir, f'hazard__{name}.parquet')
    hazard.to_parquet(output_path, index=False, compression='snappy')

    size_mb = os.path.getsize(output_path) / 1024**2
    print(f'  Saved   : {output_path} ({size_mb:.1f} MB)')
    print(f'  Shape   : {hazard.shape}')
    print(f'  Columns : FID + {len(lambda_df.columns)} lambda + '
          f'{len(var_df.columns)} var + {len(sigma_df.columns)} sigma + '
          f'{len(rp_df.columns)} H_at_rp')



def main(args):
   #print(args) 
   thresholds     = np.arange(args.thresholds[0], args.thresholds[1], args.thresholds[2])
   return_periods = np.array(args.return_periods)
   name           = args.tile_name

   print(f'='*60)
   print(f' FLOOD HAZARD COMPUTATION ')
   print(f'='*60)
   print(f'  Flood matrix   : {args.flood_matrix}')
   print(f'  Output dir     : {args.output_dir}')
   print(f'  Name           : {name}')
   print(f'  Thresholds     : {thresholds[0]:.0f} - {thresholds[-1]:.0f} cm, at interval {args.thresholds[2]} cm')
   print(f'  Nyears         : {args.nyears:,.0f}')
   print(f'  Return periods : {return_periods.tolist()}')
   print(f'  n_jobs         : {args.n_jobs}')
   print(f'='*60)

   t0 = time.time()

   # Step 1: load or compute scenario count
   print(f'\n=== BUILDING SCENARIO COUNT ===')
   df_count, n_scenarios = load_or_compute_scenario_count(
           args.flood_matrix, thresholds, args.output_dir, name, args.n_jobs,
           )
   print(f'  n_scenarios : {n_scenarios:,}')

   # Step 2: annual exceedance rate
   lambda_df, var_df, sigma_df = compute_annual_exceedance_rate(
           df_count, thresholds, n_scenarios, args.nyears,
           )

   # plot hazard curves, quick plot
   print(f'\n=== PLOT QUICK HAZARD CURVES ===')
   plot_hazard_curves(lambda_df, sigma_df, df_count, thresholds, args.output_dir, name, args.nyears)

   # Step 3: flood depth at return periods
   rp_df = compute_h_at_return_periods(
           lambda_df, thresholds, return_periods, args.n_jobs,
           )

   # Step 4: save
   save_output(df_count, lambda_df, var_df, sigma_df, rp_df, args.output_dir, name)

   total = time.time() - t0
   print(f'\nTotal elapsed : {total:.1f}s ({total/60:.1f} min)')
   print(f'DONE')
   print(f'='*60)


   return lambda_df, var_df, sigma_df, rp_df

if __name__ == '__main__':
   parser = argparse.ArgumentParser(
           description='Compute flood hazard from a flood depth matrix',
           formatter_class=argparse.RawDescriptionHelpFormatter,
           epilog=__doc__,
           )
   parser.add_argument(
           '--flood_matrix',
           default='/home/ryan/SJava_PTRA_analysis/flood_matrix/flood_depth_matrix__Tile_5__NLSWE.parquet',
           help='Path to flood_depth_matrix__<name>.parquet',
           )
   parser.add_argument(
           '--output_dir',
           default='/home/ryan/SJava_PTRA_analysis/flood_hazard',
           help='Directory to save output files',
           )
   parser.add_argument(
           '--tile_name',
           default='Tile_5__NLSWE',
           help='Tile identifier for naming',
           )
   parser.add_argument(
           '--thresholds',
           type=float,
           nargs=3,
           default=[10,3010,10],
           help='Flood depth thresholds in cm: start stop step',
           metavar=('START', 'STOP', 'STEP'),
           )
   parser.add_argument(
           '--nyears',
           type=int,
           default=10000,
           help='Length of catalogue in years',
           )
   parser.add_argument(
           '--return_periods',
           type=float,
           nargs='+',
           default=[500, 2500, 5000, 10000],
           help='Return periods in years',
           )
   parser.add_argument(
           '--n_jobs', type=int,
           default=2,
           help='Number of parallel jobs for joblib',
           )
   args = parser.parse_args()

   #lambda_df, var_df, sigma_df, rp_df = main(args)
