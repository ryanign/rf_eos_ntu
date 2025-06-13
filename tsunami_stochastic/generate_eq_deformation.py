"""
Ryan Pranantyo
EOS, 3 June 2025

A master script to combine unit sources tif file into one single deformation file in GeoTIFF.
Execute this once you have an events catalogue list
"""
import os, sys
import numpy as np
import pandas as pd
import xarray as xr
import rioxarray
import argparse
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pathlib import Path
from joblib import Parallel, delayed


def combine_disp(df, ii, disp_path, convert2grd=False, plotting_source=False, grd_bbox=None):
    if df.loc[ii, 'EVENTID'] != np.zeros:
        """ 
        giving EVENTID if not yet done 
        EVENTID will be the base name for everything
        """
        sffm_fname = Path(df["SFFM_filename"][ii]).name[:-4]
        sffm_fname_l = sffm_fname.split("__")[2:5]
        EVENTID_l = ["SFFM"] + sffm_fname_l + [df["SFFM_src_id"][ii]]
        EVENTID = "__".join(EVENTID_l)
    else:
        EVENTID = df["EVENTID"][ii]

    print(EVENTID)

    sffm_src_id = df["SFFM_src_id"][ii]
    cols2take = [sffm_src_id, "unit_source_filename"]

    sffm_df = pd.read_csv(df["SFFM_filename"][ii])
    sffm_df = sffm_df[cols2take].iloc[:-6]
    sffm_df[sffm_src_id] = sffm_df[sffm_src_id].astype(float)
    sffm_df = sffm_df[sffm_df[sffm_src_id] != 0.0]

    ### start to combine!
    disp_template = rioxarray.open_rasterio(sffm_df.iloc[0]["unit_source_filename"])
    disp = np.zeros_like(disp_template.data)

    for flt in sffm_df.index:
        ds = rioxarray.open_rasterio(sffm_df["unit_source_filename"][flt])
        disp = disp + ds.data * sffm_df[sffm_src_id][flt]
    
    if plotting_source:
        disp_path_f = Path(os.path.join(disp_path, "figures"))
        disp_path_f.mkdir(exist_ok = True)

        fig = plt.figure(figsize = (7,5), constrained_layout = True)
        ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
        ax.coastlines()
        source = ax.imshow(disp[0][::-1], origin = "lower",
                extent = [ds.x.min(), ds.x.max(), ds.y.min(), ds.y.max()],
                vmin = -disp.max(), vmax = disp.max(),
                cmap = "RdBu_r")
        cb = fig.colorbar(source, orientation = "vertical", shrink = 0.7, pad = 0,
                label = "displacement (m)")
        fout = os.path.join(disp_path_f, f"{EVENTID}.png")
        fig.savefig(fout)
        plt.close()

    if convert2grd:
        print("exporting to grd for JAGURS")
        gtiff_f = export2geotiff(disp_template, disp, EVENTID, disp_path)
        #print(gtiff_f, disp_path, grd_bbox)
        convert2grd_func(gtiff_f, disp_path, grd_bbox)
        ###to grd
    return disp

def convert2grd_func(gtiff, outdir, bbox):
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

def export2geotiff(da, disp, fname, outdir):
    """exporting to a geotiff file following da as the template"""
    da.data = disp
    outdir = Path(os.path.join(outdir, "Gtiff"))
    outdir.mkdir(exist_ok = True)
    outf = os.path.join(outdir, f"{fname}.tif")
    da.rio.to_raster(outf)
    return outf

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
            "--EVENT_CATALOGUE", 
            type = str,
            default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/input_files__SouthernJava/earthquake_catalogue__region/20250526__cat_6.5-8.7_100k__CATALOGUE-1.dat__EVENT_LIST__Mw6.95+.csv",
            help = "event catalogue file after executing ./source_catalogue_analysis/select_events_catalogue2use.py")
    parser.add_argument(
            "--DISP_PATH",
            type = str,
            default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/input_files__SouthernJava/earthquake_displacements/disp__Catalogue--1",
            help = "where to save the displacement file")
    parser.add_argument(
            "--DISP_COMPILATION_PATH",
            type = str,
            default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/input_files__SouthernJava/earthquake_displacements/disp__collections",
            help = "collection of displacement files have been generated by other event catalogue, they are just a symbolic link file")
    parser.add_argument(
            "--convert2grd",
            type = bool,
            default = True,
            help = "convert displacement file to a grd file for JAGURS, bbox is required")
    parser.add_argument(
            "--grd_bbox",
            type = float,
            nargs = "+",
            default = [102.0000, 117.0255, -12.5000, -3.7925, 0.0135],
            help = "disp grid file for JAGURS: xmin, xmax, ymin, ymax, resolution in degrees, MANDATORY if convert2grd is True")
    parser.add_argument(
            "--plotting_source",
            type = bool,
            default = True)
    parser.add_argument(
            "--ncpus",
            type = int,
            default = 1,
            help = "number of cpus to used")
    args = parser.parse_args()


    ### main program ###
    disp_path = Path(args.DISP_PATH)
    disp_path.mkdir(exist_ok = True)

    df = pd.read_csv(args.EVENT_CATALOGUE)
    ### giving EVENTID and path to the displacement
    for ii in df.index:
        sffm_fname = Path(df["SFFM_filename"][ii]).name[:-4]
        sffm_fname_l = sffm_fname.split("__")[2:5]
        EVENTID_l = ["SFFM"] + sffm_fname_l + [df["SFFM_src_id"][ii]]
        EVENTID = "__".join(EVENTID_l)
        df.loc[ii, "EVENTID"] = EVENTID
        df.loc[ii, "disp_path"] = os.path.join(disp_path, EVENTID)
    ### resaving the catalogue
    cat_fout = f"{args.EVENT_CATALOGUE[:-4]}__updated.csv"
    df.to_csv(cat_fout)

    ### check if there displacement files have been generated from other catalogues
    index2drop = list()
    for ii in df.index:
        EVENTID = df["EVENTID"][ii]
        disp2check_f = os.path.join(args.DISP_COMPILATION_PATH, f"{EVENTID}.grd")
        if os.path.exists(disp2check_f):
            print(f"drop index {ii} for the generation")
            index2drop.append(ii)
    df = df.drop(index = index2drop)


    plotting_source = args.plotting_source
    convert2grd = args.convert2grd
    if convert2grd:
        grd_bbox = args.grd_bbox
        if grd_bbox is None:
            print(f"you need to define grd_bbox!")
            sys.exit()
        
    
    #sys.exit()    
    ### ncpus
    if len(df) > 0:
        ncpus = args.ncpus
        Parallel(n_jobs = ncpus)(delayed(combine_disp)(df, ii, disp_path, convert2grd, plotting_source, grd_bbox) for ii in df.index)




                
