"""
Ryan Pranantyo
EOS, 21 May 2026

Description:
    Script to estimate average shear wave velocity at the upper 30 m depth (Vs30) following topographic slope
    approximation of Allen and Wald (2007, 2009). There are two options: (i) active tectonic (ii) stable continent.
    I use Table 1 from Allen and Wald (2009) as the basic approximation. The table shows a range of Vs30 (m/s).

    The main input as today is a slope angle with m/m unit in a GeoTiff format. I used QGIS to generate this file.

Usage:
    python Vs30_estimate_topographic_slope 
            --slope_geotiff_file <path to a slope geotiff file in m/m unit>
            --region ['active', 'stable']
            --output_dir <path to save outputs>
            --convert_projection [True, False]
            --projection 'epsg:4326' <if convert_projection True>

"""
import os
import numpy as np
import pandas as pd
import xarray as xr
import rioxarray
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import argparse
from pathlib import Path
import gc
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ─────────────────────────────────────────────────────────────────────────────
# Lookup table — Allen & Wald (2007) Table 1
# Format: (slope_min, slope_max, vs30_min, vs30_max, nehrp_class)
# slope in m/m, vs30 in m/s
# ─────────────────────────────────────────────────────────────────────────────

LOOKUP = {
    'active': [
        # slope_min,    slope_max,  vs30_low, vs30_high, nehrp
        (0.0,           3e-4,       100,       180,       'E'),
        (3e-4,          3.5e-3,     180,       240,       'E'),
        (3.5e-3,        0.010,      240,       300,       'D'),
        (0.010,         0.024,      300,       360,       'D'),
        (0.024,         0.08,       360,       490,       'D'),
        (0.08,          0.14,       490,       620,       'C'),
        (0.14,          0.20,       620,       760,       'C'),
        (0.20,          9999,       760,      1600,       'B'),
    ],
    'stable': [
        # slope_min,    slope_max,  vs30_low, vs30_high, nehrp
        (0.0,           1e-4,       100,       180,       'E'),
        (1e-4,          4.5e-3,     180,       240,       'E'),
        (4.5e-3,        8.5e-3,     240,       300,       'D'),
        (8.5e-3,        0.013,      300,       360,       'D'),
        (0.013,         0.022,      360,       490,       'D'),
        (0.022,         0.03,       490,       620,       'C'),
        (0.03,          0.04,       620,       760,       'C'),
        (0.04,          9999,       760,      1600,       'B'),
    ],
}

NEHRP_CODE = {'E': 1, 'D': 2, 'C': 3, 'B': 4}
NEHRP_NAME = {1: 'E', 2: 'D', 3: 'C', 4: 'B'}

def load_slope(slope_file, convert_projection, projection):
    """
    Load slope GeoTIFF, optionally reproject.
    """
    print(f'\n=== Loading slope file ===')
    slope_da = rioxarray.open_rasterio(slope_file, masked=True)
    print(f'  CRS        : {slope_da.rio.crs}')
    print(f'  Shape      : {slope_da.shape}')
    print(f'  Resolution : {slope_da.rio.resolution()}')
    print(f'  Bounds     : {slope_da.rio.bounds()}')

    if convert_projection:
        target_crs = projection.upper() if not projection.startswith('epsg') \
                else projection
        if str(slope_da.rio.crs).lower() != target_crs.lower():
            print(f'    Reprojecting to {target_crs} ...')
            slope_da = slope_da.rio.reproject(target_crs)
            print(f'  New shape  : {slope_da.shape}')
        else:
            print(f'  Already in {target_crs}, skipping reproject.')

    # Extract 2D array (band 0), ensure positive values
    slope_arr = slope_da.values[0].astype(np.float32)
    slope_arr = np.abs(slope_arr)
    slope_arr = np.where(slope_arr == 0, np.nan, slope_arr)

    # Mask nodata
    nodata = slope_da.rio.nodata
    if nodata is not None:
        slope_arr[slope_arr == nodata] = np.nan

    valid = np.isfinite(slope_arr)
    print(f'  Slope range: {np.nanmin(slope_arr):.6f} to {np.nanmax(slope_arr):.6f} m/m')
    print(f'  Valid pixel: {valid.sum():,} of {slope_arr.size:,}')
    
    return slope_da, slope_arr

def slope_to_vs30(slope_arr, region):
    """
    CORE FUNCTION to convert slope to vs30 and NEHRP site class
    """
    lookup = LOOKUP[region]

    vs30  = np.full(slope_arr.shape, np.nan, dtype=np.float32)
    nehrp = np.full(slope_arr.shape, np.nan, dtype=np.float32)

    for (s_min, s_max, v_low, v_high, nehrp_class) in lookup:
        mask = (slope_arr >= s_min) & (slope_arr < s_max)
        if not mask.any():
            continue

        if s_max >= 9999:
            # above max threshold, fixed Vs30
            vs30[mask] = v_high
        else:
            # linear interpolation within range
            t          = (slope_arr[mask] - s_min) / (s_max - s_min)
            vs30[mask] = v_low + t * (v_high - v_low)

        nehrp[mask] = NEHRP_CODE[nehrp_class]

    return vs30, nehrp

def save_outputs(vs30, nehrp, slope_da, output_dir, stem, region):
    print(f'\n=== SAVING Vs30 GeoTIFF ===')
    os.makedirs(output_dir, exist_ok=True)
    vs30_da = slope_da.copy(data=vs30[np.newaxis, :, :].astype(np.float32))
    vs30_da.attrs['long_name'] = 'Vs30 (m/s)'
    vs30_da.attrs['units']     = 'm/s'
    vs30_da.attrs['source']    = 'Allen & Wald (2009) Table 1'
    vs30_da.attrs['region']    = region

    vs30_path = os.path.join(output_dir, f'{stem}__vs30__{region}.tif')
    vs30_da.rio.to_raster(vs30_path, dtype='float32')
    print(f'  Saved Vs30  : {vs30_path}')

    print(f'\n=== SAVING NEHRP Site Class GeoTIFF ===')
    nehrp_da = slope_da.copy(data=nehrp[np.newaxis, :, :].astype(np.float32))
    nehrp_da.attrs['long_name'] = 'NEHRP site class (1=E, 2=D, 3=C, 4=B)'
    nehrp_da.attrs['source']    = 'Allen & Wald (2009) Table 1'
    nehrp_da.attrs['region']    = region

    nehrp_path = os.path.join(output_dir, f'{stem}__nehrp__{region}.tif')
    nehrp_da.rio.to_raster(nehrp_path, dtype='float32')
    print(f'  Saved NEHRP : {nehrp_path}')


    print(f'\n=== QUICK PLOTTING ===')
    xmin, ymin, xmax, ymax = slope_da.rio.bounds()
    fig = plt.figure(figsize=(8,3.8), constrained_layout=True)
    for ii, var in enumerate([vs30, nehrp]):
        ax = fig.add_subplot(1,2,ii+1, projection=ccrs.PlateCarree())
        #ax.coastlines(linewidth=1, zorder=10, resolution='10m')

        if ii == 0:
            ax.set_title(f'Vs30 Estimation - {region.capitalize()}', fontsize=8)
            im = ax.imshow(
                    vs30, origin='upper',
                    extent=[xmin, xmax, ymin, ymax],
                    cmap='RdYlGn', vmin=100, vmax=800,)
            cb = fig.colorbar(im, shrink=0.8, orientation='horizontal', extend='both',)
            cb.set_label('Vs30 (m/s)', fontsize=8)

        else:
            ax.set_title(f'NEHRP Classifcation', fontsize=8)
            cmap_nehrp = mcolors.ListedColormap(['#d73027', '#fc8d59', '#91cf60', '#1a9850'])
            bounds     = [0.5, 1.5, 2.5, 3.5, 4.5]
            norm_nehrp = mcolors.BoundaryNorm(bounds, cmap_nehrp.N)
            im = ax.imshow(
                    nehrp, origin='upper',
                    extent=[xmin, xmax, ymin, ymax],
                    cmap=cmap_nehrp, norm=norm_nehrp,)
            cb = fig.colorbar(im, shrink=0.8, orientation='horizontal', ticks=[1,2,3,4],)
            cb.set_ticklabels(['E', 'D', 'C', 'B'])
            cb.set_label('NEHRP Site Class', fontsize=8)

    png_path = os.path.join(output_dir, f'{stem}__vs30__{region}.png')
    fig.savefig(png_path, dpi=300)
    plt.close()
    print(f'  Saved PNG : {png_path}')


def main(args):
    region = args.region.lower()
    if region not in ('active', 'stable'):
        raise ValueError(f'--region must be "active" or "stable", got "{args.region}"')

    stem       = Path(args.slope_gtiff_file).stem
    output_dir = args.output_dir if args.output_dir \
            else str(Path(args.slope_gtiff_file).parent)

    print(f'='*60)
    print(f' Vs30 ESTIMATION FROM TERRAIN SLOPE')
    print(f' Allen & Wald (2009) Table 1')
    print(f'='*60)
    print(f'  Input  : {args.slope_gtiff_file}')
    print(f'  Region : {region}')
    print(f'  Output : {output_dir}')
    if args.convert_projection:
        print(f'    - Reproject to : {args.projection}')
    print(f'='*60)


    ### Load slope tiff file
    slope_da, slope_arr = load_slope(
            args.slope_gtiff_file,
            args.convert_projection,
            args.projection,
            )

    ### Classify
    print(f'\n=== CLASSIFYING ===')
    vs30, nehrp = slope_to_vs30(slope_arr, region)
    print(f'  Done.')
    print(f'  Vs30 shape: {vs30.shape}')

    ### saving
    save_outputs(vs30, nehrp, slope_da, output_dir, stem, region)

    return slope_da, vs30, nehrp



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description = 'Vs30 estimation based on terrain slope of Allen & Wald (2007, 2009)',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=__doc__,
            )
    parser.add_argument(
            '--slope_gtiff_file',
            default = '/home/ryan/OneDrive_NTU_Projects/202602XX__SHA_Singapore/raw_data/DEM_Singapore/slope__fathomdem__singapore__utm48n.tif',
            type = str,
            help = 'Slope angle in m/m unit in a GeoTIFF file format',
            )
    parser.add_argument(
            '--region',
            type = str,
            default = 'stable',
            choices = ['active', 'stable'],
            help = 'active: Active tectonic | stable: Stable continent (default: stable)',
            )
    parser.add_argument(
            '--output_dir',
            type = str,
            default = '',
            help = 'Output directory (default: same as input file)',
            )
    parser.add_argument(
            '--convert_projection',
            type = bool,
            default = True,
            help = 'True or False to convert projection',
            )
    parser.add_argument(
            '--projection',
            type = str,
            default = 'epsg:4326',
            help = 'if convert_projection is True, an epsg projection code should be added',
            )
    args = parser.parse_args()

    slope_da, vs30, nehrp = main(args)
