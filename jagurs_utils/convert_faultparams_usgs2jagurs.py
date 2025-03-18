"""
Ryan Pranantyo
EOS, March 2025

A basic script to convert finite fault inversion solution 
downloaded from USGS to a standard format for JAGURS
> we use .param file
> later, .fsp file can be used as well
"""

import os, sys
import numpy as np
import pandas as pd
import argparse
from pathlib import Path

def write_jagurs_fault(df, infile):
    """ 
    JAGURS needs
    lat(deg) lon(deg) depth(km) length(km) width(km) dip(deg) strike(deg) rake(deg) slip_amp(m)
    lat, lon, and depth = top-right position
    """
    parent_dir = infile.parent
    fname = infile.name[:-len(infile.suffix)]

    df["slip"] = df["slip"] / 100.
    df_jagurs = df[["lat_tc", "lon_tc", "depth_tc", "length", "width", "dip", "strike", "rake", "slip"]]
    fout = os.path.join(parent_dir, f"jagurs__{fname}.flt")

    #ff = open(fout, "w")
    #ff.write("!lat lon depth length width dip strike rake slip\n")
    #for idx in df_jagurs.index:
    #    ff.write(f"{df['lat_tc'][idx]:.5f}\n")
    #ff.close()
    df_jagurs = df_jagurs.round(decimals=5).astype("float32")
    df_jagurs.to_csv(fout, header=None, index=None, sep=" ")
    
    return df_jagurs

def find_solution_start_from_row(infile):
    ff = open(infile, "r")
    lines = ff.readlines()
    text = "#Lat. Lon. depth slip rake strike dip t_rup t_ris t_fal mo"
    for row, line in enumerate(lines):
        if text in line:
            return row

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--USGS_ffi_solution", type = str,
            default = "./examples/2006Java__basic_inversion.param")
    parser.add_argument("--dx_dy", type = float, nargs = "+",
            default = [20.00, 11.00],
            help = "dx = along strike length, dy = downdip width")
    args = parser.parse_args()

    infile = Path(args.USGS_ffi_solution)
    dx,dy = args.dx_dy[0], args.dx_dy[1]

    row = find_solution_start_from_row(infile)
    df = pd.read_table(infile, skiprows=row, sep='\\s+')
    df = df.rename(columns = {"#Lat." : "lat", "Lon."  : "lon"})
    
    df["length"] = dx
    df["width"]  = dy
    df["depth_tc"] = df.depth - (np.sin(np.radians(df.dip)) * (dy / 2))

    dist_dx = np.sin(np.radians(360 - df.strike)) * dx/2 / 111.11
    dist_dy = np.cos(np.radians(360 - df.strike)) * dy/2 / 111.11

    df["lon_tc"] = df["lon"] + dist_dx
    df["lat_tc"] = df["lat"] - dist_dy

    df_jagurs = write_jagurs_fault(df, infile)

    print(df_jagurs)




