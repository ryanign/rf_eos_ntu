"""
Ryan Pranantyo
EOS, 24 March 2026

Description:
    bla bla

Output:
    -

Useage:


"""
import os
import time
import warnings
import argparse
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from rasterio.features import rasterize
from rasterio.transform import from_bounds
from tqdm import tqdm
from pathlib import Path

warnings.filterwarnings("ignore", message="Geometry is in a geographic CRS")

# -------------------------------------------------------------------------
# Step 1: loading data
# -------------------------------------------------------------------------
def load_data(buildings, flood_db):
    print(f'\n=== LOADING DATA ===')
    gdf = gpd.read_file(buildings)
    nc = xr.open_dataset(flood_db)
    print(f'  Building footprints : {len(gdf):,}')
    print(f'  Scenarios           : {nc.sizes["eventid"]:,}')
    print(f'  Grid size           : {nc.sizes["lat"]} x {nc.sizes["lon"]}')
    return gdf, nc

# -------------------------------------------------------------------------
# Step 2: rasterize building footprints onto flood grid
#  - building footprits are converted into raster for masking purpose
# -------------------------------------------------------------------------
def rasterize_buildings(gdf, lat_vals, lon_vals):
    print(f'\n=== RASTERIZE BUILDING FOOTPRINTS ===')
    nrows = len(lat_vals)
    ncols = len(lon_vals)

    transform = from_bounds(
            lon_vals.min(), lat_vals.min(),
            lon_vals.max(), lat_vals.max(),
            ncols, nrows,
            )

    shapes = [(geom, fid) for geom, fid in zip(gdf.geometry, gdf['FID'])]

    label_raster = rasterize(
            shapes=shapes,
            out_shape=(nrows,ncols),
            transform=transform,
            fill=0,
            dtype='int32',
            all_touched=False,     # only pixels whose centre falls inside polygon
            )

    print(f'  Unique buildings in raster : {len(np.unique(label_raster)) - 1:,}')
    print(f'  Pixels with building       : {(label_raster > 0).sum():,}')

    return label_raster

# -------------------------------------------------------------------------
# Step 3: centroid fallback for sub-pixel buildings
# -------------------------------------------------------------------------
def build_centroid_fallback(gdf, label_raster, lat_vals, lon_vals):
    print(f'\n=== CENTROID FALLBACK ===')
    nrows = len(lat_vals)
    ncols = len(lon_vals)

    covered_fids     = set(np.unique(label_raster)) - {0}
    mask_not_covered = ~gdf['FID'].isin(covered_fids)
    gdf_fallback     = gdf[mask_not_covered].copy()

    centroids          = gdf_fallback.geometry.centroid
    gdf_fallback['cx'] = centroids.x
    gdf_fallback['cy'] = centroids.y

    # convert lon/lat coordinates -> pixel cor/row indices
    gdf_fallback['col'] = (
            (gdf_fallback['cx'] - lon_vals.min()) / 
            (lon_vals.max() - lon_vals.min()) * (ncols - 1)
            ).round().astype(int).clip(0, ncols - 1)

    gdf_fallback['row'] = (
            (lat_vals.max() - gdf_fallback['cy']) /
            (lat_vals.max() - lat_vals.min()) * (nrows - 1)
            ).round().astype(int).clip(0, nrows - 1)

    outside = (
            (gdf_fallback['cx'] < lon_vals.min()) | (gdf_fallback['cx'] > lon_vals.max()) |
            (gdf_fallback['cy'] < lat_vals.min()) | (gdf_fallback['cy'] > lat_vals.max())
            )

    print(f'  Covered by raster   : {len(covered_fids):,}')
    print(f'  Centroid fallback   : {len(gdf_fallback):,}')
    print(f'  Outside grid extent : {outside.sum():,} (will be 0 for all scenarios)')

    return gdf_fallback, outside

# -------------------------------------------------------------------------
# Step 4: building a lookup table
# -------------------------------------------------------------------------
def build_lookup_table(gdf, label_raster, gdf_fallback, outside):
    print(f'\n=== BUILD LOOKUP TABLE ===')
    # from raster - buildings covering multiple pixes have multiple entries
    rows_r, cols_r = np.where(label_raster > 0)
    df_raster = pd.DataFrame({
            'FID'    : label_raster[rows_r, cols_r],
            'row'    : rows_r,
            'col'    : cols_r,
            'source' : 'raster',
            })

    # from centroid fallback - inside grid only
    gdf_inside = gdf_fallback[~outside].copy()
    df_centroid = pd.DataFrame({
        'FID'    : gdf_inside['FID'].values,
        'row'    : gdf_inside['row'].values,
        'col'    : gdf_inside['col'].values,
        'source' : 'centroid',
        })

    df_lookup = pd.concat([df_raster, df_centroid], ignore_index=True)

    print(f'  Raster entries         : {len(df_raster):,}')
    print(f'  Centroid entries       : {len(df_centroid):,}')
    print(f'  Total lookup entries   : {len(df_lookup):,}')
    print(f'  Buildings in lookup    : {df_lookup["FID"].nunique():,}')
    print(f'  Buildings outside grid : {len(gdf) - df_lookup["FID"].nunique():,}')

    return df_lookup

# -------------------------------------------------------------------------
# Step 5: scan non-zero scenarios
# -------------------------------------------------------------------------
def scan_nonzero_scenarios(nc):
    print(f'\n=== SCAN NON-ZERO SCENARIOS ===')
    n_events       = nc.sizes['eventid']
    nonzero_events = []
    scan_batch     = 50  # fixed smal batch for scanning

    for i in tqdm(range(0, n_events, scan_batch), desc='  Scanning'):
        end = min(i + scan_batch, n_events)
        batch = nc['flood_depth'].isel(eventid=slice(i, end)).values
        max_per_event = batch.max(axis=(1,2))
        nonzero_in_batch = np.where(max_per_event > 0)[0] + i
        nonzero_events.extend(nonzero_in_batch.tolist())

    nonzero_events = np.array(nonzero_events)
    print(f'  Non-zero scenarios  : {len(nonzero_events):,}')
    print(f'  Zero-only scenarios : {n_events - len(nonzero_events):,}')
    return nonzero_events

# -------------------------------------------------------------------------
# Step 6: extract flood depth values
# -------------------------------------------------------------------------
def extract_flood_depth(nc, df_lookup, all_fids, nonzero_events, batch_size):
    print(f'\n=== EXTRACT FLOOD DEPTH ===')
    n_buildings = len(all_fids)
    n_nonzero   = len(nonzero_events)
    n_events    = nc.sizes['eventid']

    # Build pixel lookup arrays
    lookup_rows    = df_lookup.groupby('FID')['row'].apply(np.array)
    lookup_cols    = df_lookup.groupby('FID')['col'].apply(np.array)
    fids_in_lookup = np.array(lookup_rows.index)
    fid_to_outrow  = {fid: i for i, fid in enumerate(all_fids)}
    out_indices    = np.array([fid_to_outrow[fid] for fid in fids_in_lookup])

    pixel_rows = np.array([lookup_rows[fid] for fid in fids_in_lookup], dtype=object)
    pixel_cols = np.array([lookup_cols[fid] for fid in fids_in_lookup], dtype=object)

    # Seperate single-pixel (vectorized) vs multi-pixel (small loop, take max)
    pixel_counts = np.array([len(r) for r in pixel_rows])
    single_mask  = pixel_counts == 1
    multi_mask   = pixel_counts > 1

    single_out_indices = out_indices[single_mask]
    single_rows = np.array([pixel_rows[i][0] for i in np.where(single_mask)[0]])
    single_cols = np.array([pixel_cols[i][0] for i in np.where(single_mask)[0]])
    multi_out_indices = out_indices[multi_mask]
    multi_pixel_rows  = pixel_rows[multi_mask]
    multi_pixel_cols  = pixel_cols[multi_mask]

    print(f'  Single-pixel buildings : {single_mask.sum():,}')
    print(f'  Multi-pixel buidlings  : {multi_mask.sum():,}')

    # Build sequential batches
    nonzero_set = set(nonzero_events.tolist())
    batches = []
    for batch_start in range(0, n_events, batch_size):
        batch_end = min(batch_start + batch_size, n_events)
        nonzero_in_batch = [i for i in range(batch_start, batch_end) if i in nonzero_set]
        if nonzero_in_batch:
            batches.append((batch_start, batch_end, nonzero_in_batch))
    print(f'  Total batches : {len(batches)}')

    # Result matrix: only non-zero scenarios (expanded to full n_events when saving)
    result = np.zeros((n_buildings, n_nonzero), dtype=np.int16)
    event_to_col = {ev: col for col, ev in enumerate(nonzero_events)}

    start_time = time.time()

    for batch_start, batch_end, nonzero_in_batch in tqdm(batches, desc='  Extracting'):
        # Load sequential slice: shape (batch_size, nrows, ncols)
        flood_batch = nc['flood_depth'].isel(eventid=slice(batch_start, batch_end)).values

        col_indices = np.array([event_to_col[i] for i in nonzero_in_batch])
        batch_rel   = np.array([i - batch_start for i in nonzero_in_batch])
        flood_nonzero = flood_batch[batch_rel]

        # Single-pixel buildings: fully vectorized
        result[np.ix_(single_out_indices, col_indices)] = \
                flood_nonzero[:, single_rows, single_cols].T 

        # Multi-pixel buildings: small loop, take max across pixels
        for b_idx, (p_rows, p_cols) in enumerate(zip(multi_pixel_rows, multi_pixel_cols)):
            vals = flood_nonzero[:, p_rows, p_cols]     # (n_nonzero_in_batch, n_pixels)
            result[multi_out_indices[b_idx], col_indices] = vals.max(axis=1)

    elapsed = time.time() - start_time
    print(f'  Done in {elapsed:.1f}s ({elapsed/60:.1f} min')
    print(f'  Non-zero values : {(result > 0).sum():,}')

    return result
    
# -------------------------------------------------------------------------
# Step 7: saving
# -------------------------------------------------------------------------
def save_output(result, nonzero_events, all_fids, nc, output_dir, tile_name):
    print(f'\n=== SAVING OUTPUT ===')
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    n_buildings  = len(all_fids)
    n_events     = nc.sizes['eventid']
    all_eventids = nc['eventid'].values

    # Expand result (n_buildings x n_nonzero) -> full matrix (n_builds x n_events)
    flood_matrix = pd.DataFrame(
            np.zeros((n_buildings, n_events), dtype=np.int16),
            columns=all_eventids,
            )
    flood_matrix.iloc[:, nonzero_events] = result
    flood_matrix.insert(0, 'FID', all_fids)

    output_path = os.path.join(output_dir, f'flood_depth_matrix__{tile_name}.parquet')
    flood_matrix.to_parquet(output_path, index=False, compression='snappy')

    size_mb = os.path.getsize(output_path) / 1024**2
    print(f'  Saved : {output_path} ({size_mb:.1f} MB)')
    print(f'  Shape : {flood_matrix.shape}')
    
    return flood_matrix


def main(args):
    print(f'='*60)
    print(f' FLOOD DEPTH MATRIX EXTRACTION ')
    print(f'='*60)
    print(f' -Buildings  : {args.buildings}')
    print(f' -Flood DB   : {args.flood_db}')
    print(f' -Output dir : {args.output_dir}')
    print(f' -Tile       : {args.tile_name}')
    print(f' -Batch size : {args.batch_size}')
    print(f'='*60)

    t0 = time.time()

    # Step 1: load data
    gdf, nc = load_data(args.buildings, args.flood_db)
    lon_vals = nc['lon'].values
    lat_vals = nc['lat'].values

    # Step 2: spatial lookup
    label_raster = rasterize_buildings(gdf, lat_vals, lon_vals)
    
    # Step 3: 
    gdf_fallback, outside = build_centroid_fallback(gdf, label_raster, lat_vals, lon_vals)

    # Step 4: build lookup table
    df_lookup = build_lookup_table(gdf, label_raster, gdf_fallback, outside)

    # Step 5: scan non-zeros scenarios
    nonzero_events = scan_nonzero_scenarios(nc)

    # Step 6: extract
    all_fids = gdf['FID'].values
    result   = extract_flood_depth(nc, df_lookup, all_fids, nonzero_events, args.batch_size)

    # Step 7: save output
    flood_matrix = save_output(result, nonzero_events, all_fids, nc, args.output_dir, args.tile_name)

    total = time.time() - t0
    print(f'\n')
    print(f'='*60)
    print(f'Total elapsed : {total:.1f}s ({total/60:.1f} min)')
    print(f'DONE!')
    print(f'='*60)


    return gdf, nc, flood_matrix

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Extract flood depth matrix for building footprints from a NetCDF flood database.',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=__doc__,
            )
    parser.add_argument(
            '--buildings', 
            #required=True,
            default='/home/ryan/Microsoft__GlobalML__Indonesia/extract/within_sfincs_active_cells/tile_5__extracted.shp',
            help='Path to building footprints shapefile (.shp)'
            )
    parser.add_argument(
            '--flood_db',
            #required=True,
            default='/home/ryan/SJava_PTHA_analysis/simulations_output/CATALOGUE-010000-years/onshore_database/flood_depth__tile_5__CATALOGUE__10000-YEARS__NONLINEAR.nc',
            help='Path to flood depth NetCDF file (.nc), in this case from Pranantyo et al., 2026'
            )
    parser.add_argument(
            '--output_dir',
            #required=True,
            default='/home/ryan/SJava_PTRA_analysis/flood_matrix',
            help='Directory to save output files',
            )
    parser.add_argument(
            '--tile_name',
            #required=True,
            default='Tile_5__NLSWE',
            help='Tile identified used in output filename',
            )
    parser.add_argument(
            '--batch_size',
            type=int,
            default=100,
            help='Scenarios loaded per batch. Higher = faster but more RAM.',
            )
    args = parser.parse_args()

    gdf, nc, flodd_matrix = main(args)

