"""
Ryan Pranantyo
EOS, February 2025

this script is just for visualisastion of stochastic slip models generated
-- later will be added to generate one tiff or grd file per model for explicit modelling
"""
import os, sys
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import xarray as xr
import rioxarray
import argparse
from pathlib import Path

def export2geotiff(da, disp, fin, nn, outdir):
    """ function to export displacement to a tiff file"""
    print("    >> converting to geotiff ...")
    da.data = disp
    outdir = Path(os.path.join(outdir, "gtiff_temp"))
    outdir.mkdir(exist_ok = True)
    outf = os.path.join(outdir, f"{fin.name[:-4]}_{nn}.tif")
    da.rio.to_raster(outf)
    return outf

def convert2grd(gtiff, outdir, bbox):
    """ function to convert displacement file for a tiff to grd file
        this is for an input file for JAGURS """
    ### convert gtiff to a temp grd file
    tmpfile = f"{gtiff[:-4]}_tmp.grd"
    cmd = f"gmt grdconvert {gtiff} -G{tmpfile}"
    os.system(cmd)

    ### make sure the gtiff file follow jagurs bbox
    tmpfile2 = f"{tmpfile[:-4]}-2.grd"
    cmd = f"gmt grdcut {tmpfile} -R{bbox[0]-1}/{bbox[1]+1}/{bbox[2]-1}/{bbox[3]+1} -N0 -G{tmpfile2}"
    os.system(cmd)
    print(f"      <<{cmd}>>")

    ### prepare the grdfile
    outfile = f"{gtiff[:-4]}.grd"
    cmd = f"gmt grdsample {tmpfile2} -R{bbox[0]}/{bbox[1]}/{bbox[2]}/{bbox[3]} -I{bbox[4]}/{bbox[4]} -T -G{outfile}=10"
    os.system(cmd)
    print(f"      <<{cmd}>>")

    ### cleaning up
    cmd = f"mv -f {outfile} {outdir}"
    os.system(cmd)
    cmd = f"rm -f {tmpfile} {tmpfile2}"
    os.system(cmd)
    return

def main_plotting(args):
    print(args)
    fin = Path(args.SFFM_model)
    outdir = fin.parent
    df = pd.read_csv(fin)

    num_samples = args.Nsamples

    data_template = rioxarray.open_rasterio(df['unit_source_filename'][0])
    #fig = plt.figure(figsize=(10,10), constrained_layout=True)
    for nn in range(num_samples):
        cols = ['unit_source_index', 'unit_source_filename', f'unit_source_slip__{nn}']
        dfn = df[cols]
        dfn = dfn.loc[(dfn[f'unit_source_slip__{nn}'] != 0)]
        print(dfn)

        fig = plt.figure(figsize=(7,5), constrained_layout = True)
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
        ax.set_title(f'random slip {nn} : ')
        ax.coastlines()
        disp = np.zeros_like(data_template.data)
        for flt in dfn.index:
            #print(flt)
            ds = rioxarray.open_rasterio(dfn['unit_source_filename'][flt])
            disp = disp + ds.data * dfn[f'unit_source_slip__{nn}'][flt]
            

        if args.export2grd:
            print("  converting displacement to a grd file ...")
            gtiff = export2geotiff(data_template, disp, fin, nn, args.displacement_outdir)
            convert2grd(gtiff, args.displacement_outdir, args.grd_bbox)

        if args.plotting_source:
            fig = plt.figure(figsize=(7,5), constrained_layout = True)
            ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
            ax.coastlines()
            source = ax.imshow(disp[0][::-1], origin='lower',
                     extent=[ds.x.min(), ds.x.max(), ds.y.min(), ds.y.max()],
                     vmin = -2, vmax = 2,
                     cmap = 'RdBu_r')
            if args.bbox_plot is not None:
                bbox = args.bbox_plot
                ax.set_xlim(bbox[0], bbox[1])
                ax.set_ylim(bbox[2], bbox[3])
            
            fig.colorbar(source, orientation='vertical', shrink=0.7, extend='both')
            fout = os.path.join(outdir, f"{fin.name[:-4]}__{nn}.png")
            fig.savefig(fout)
            plt.close(fig)







if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--SFFM_model", type=str,
            default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/OUTPUTS/stochastic_slips__sumatera_jawa__slab2__edited/stochastic_sources__Mw_7.900000__Lon_104.104800__Lat_-7.122200__table_simplified.csv",
            help = "a list of stochastic finite fault model generated from automate_stochastic_generation.py")
    parser.add_argument("--Nsamples" , type=int,
            default = 25,
            help = "Number of visulisation to be done, max is number of stochastic samples per file")
    parser.add_argument("--bbox_plot", type=float, nargs="+", default=None,
            help = "bbox for plotting [xmin, xmax, ymin, ymax]")
    parser.add_argument("--export2grd", type=bool, default=True, #False,
            help = "convert to a grd file?")
    parser.add_argument("--grd_bbox", type=float, nargs="+", default=[86,122,-17,18,0.0135],
            help = "disp grid file for JAGURS: xmin, xmax, ymin, ymax, resolution in degrees")
    parser.add_argument("--displacement_outdir", type=str, 
            default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/OUTPUTS/displacement_ready2use",
            help = "where to save displacement grd file for JAGURS")
    parser.add_argument("--plotting_source", type=bool, default=True)

    args = parser.parse_args()
    main_plotting(args)
