"""
Ryan Pranantyo
EOS, 12 June 2025

a script to convert finite-fault parameters from http://equake-rc.info/srcmod/ to a JAGURS fault format
"""
import os, sys
import numpy as np
import pandas as pd
import argparse
from pathlib import Path
from pyproj import Geod

#def fsp2fault():
#    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fsp_file", 
                        type = str,
                        default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/Benchmarking__Historical/2006Java/s2006SOUTHE01JIxx.fsp",
                        help = "an .fsp file from srcmod format")
    parser.add_argument("--dx_dz",
                        type = float,
                        nargs = "+",
                        default = [15.00, 11.00],
                        help = "subfault length along strike (km) and subfault width along downdip (km), check inside .fsp file")
    parser.add_argument("--strike_dip",
                        type = float,
                        nargs = "+",
                        default = [288.940, 10.350],
                        help = "strike in degree and dip angle in degree, check inside .fsp file")
    parser.add_argument("--first_row",
                        type = int,
                        default = 47,
                        help = "to help finding the first row of the finite fault param")
    args = parser.parse_args()

    infile = Path(args.fsp_file)
    dx , dz = args.dx_dz[0] , args.dx_dz[1] 
    strike , dip = args.strike_dip[0], args.strike_dip[1]


    ### load the data
    df = pd.read_table(infile, skiprows = args.first_row-1, header = None, sep = "\\s+")
    df = df.rename(
            columns = {
                0 : "lat",
                1 : "lon",
                2 : "X==EW",
                3 : "Y==NS",
                4 : "z",
                5 : "slip",
                6 : "rake"}
            )
    df["strike"] = strike
    df["dip"] = dip
    df["length"] = dx
    df["width"] = dz

    ### start to translate lon-lat-z from center-top to right-top corner
    g = Geod(ellps="WGS84")
    for ii in df.index:
        xrt, yrt, _ = g.fwd(df["lon"][ii], df["lat"][ii], df["strike"][ii] - 180, 0.5 * df["length"][ii] * 1000)
        df.loc[ii, "lon_tr"] = xrt
        df.loc[ii, "lat_tr"] = yrt


    ### fsp for jagurs format
    df_jgrs = df[["lat_tr", "lon_tr", "z", "length", "width", "dip", "strike", "rake", "slip"]]
    df_jgrs = df_jgrs.round(decimals = 5).astype("float32")

    ### exporting 
    fname = infile.name[:-4] + "__jagurs.flt"
    fout = os.path.join(infile.parent, fname)
    df_jgrs.to_csv(fout, header = None, index = None, sep = " ")

    print(df_jgrs)
    














