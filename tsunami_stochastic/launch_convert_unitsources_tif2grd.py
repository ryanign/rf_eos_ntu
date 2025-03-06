"""
Ryan Pranantyo
EOS, March 2025

this is to launch unit source displacement file from tif to grd
"""
import os, sys
import pandas as pd
from joblib import Parallel, delayed

def launch(script, gtiff, bbox, outdir):
    cmd = f"python -W ignore {script} --gtiff_unitsource {gtiff} --grid_bbox {bbox} --displacement_outdir {outdir}"
    print(cmd)
    os.system(cmd)
    return cmd

if __name__ == "__main__":
    gtiff_path = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/OUTPUTS/Unit_source_data/sumatera_jawa__slab2__edited"
    grid_bbox = [102.0000, 117.0660, -12.5000, -3.7520, 0.027]

    main_outdir = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/jagurs_runs/uji_coba/unit_sources_gf/unit_sources__grid_2700m"
    
    script = "convert_unitsources_tif2grd.py"
    
    unit_src_df = pd.read_csv("/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/OUTPUTS/unit_source_grid/unit_sources_torun_testing_jawa.csv")
    
    bbox = f"{grid_bbox[0]} {grid_bbox[1]} {grid_bbox[2]} {grid_bbox[3]} {grid_bbox[4]}"

    gtiff_file = []
    unit_outdir = []
    for ii in unit_src_df.index:
        downdip = unit_src_df['dwndp_n'][ii]
        alongstrike = unit_src_df['alngst_'][ii]
        infile = os.path.join(gtiff_path, f"sumatera_jawa__slab2__edited_{downdip}_{alongstrike}.tif")
        gtiff_file.append(infile)
        outfile = os.path.join(main_outdir, f"unit_source__{downdip}_{alongstrike}")
        unit_outdir.append(outfile)

    unit_src_df['gtiff_file'] = gtiff_file
    unit_src_df['grd_outdir'] = unit_outdir

    ###launch
    ii = 0
    cmd = launch(script, unit_src_df["gtiff_file"][ii], bbox, unit_src_df["grd_outdir"][ii])

    Parallel(n_jobs = 8)(delayed(launch)(script,unit_src_df["gtiff_file"][ii], bbox, unit_src_df["grd_outdir"][ii]) for ii in unit_src_df.index)




