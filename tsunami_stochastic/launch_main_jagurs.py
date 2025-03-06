"""
Ryan Pranantyo
EOS, March 2025

script to automate launch many jagurs simulations 
"""
import os, sys
import pandas as pd
from pathlib import Path
from joblib import Parallel, delayed

def launch(script, disp, jagurs_template, outdir_path, job_name, ncpus=4, ngrids=2):

    cmd = f"python -W ignore {script} --displacement_grdfile {disp} --jagurs_template {jagurs_template} --where_to_run {outdir_path} --job_name {job_name} --ncpus {ncpus}  --num_of_grids {ngrids}"
    print(cmd)
    os.system(cmd)
    return cmd

if __name__ == "__main__":
    script = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/scripts/main_jagurs.py"
    jagurs_template = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/jagurs_runs/uji_coba/template_20250228_jawa_2700m_to_0900m"
    main_displacement_path = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/jagurs_runs/uji_coba/unit_sources_gf/unit_sources__grid_2700m"
    main_outdir_path = Path("/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/jagurs_runs/uji_coba/unit_sources_gf__tsunami/tsunami_gf__2700m")
    main_outdir_path.mkdir(exist_ok = True)

    ### main ###
    unit_src_df = pd.read_csv("/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/OUTPUTS/unit_source_grid/unit_sources_torun_testing_jawa.csv")
    disp_file = []
    outdir_paths = []
    job_names = []
    for ii in unit_src_df.index:
        downdip = unit_src_df['dwndp_n'][ii]
        alongstrike = unit_src_df['alngst_'][ii]
        infile = os.path.join(main_displacement_path, f"unit_source__{downdip}_{alongstrike}", f"sumatera_jawa__slab2__edited_{downdip}_{alongstrike}.grd")
        disp_file.append(infile)
        outfile = os.path.join(main_outdir_path, f"gf__{downdip}_{alongstrike}")
        outdir_paths.append(outfile)
        job_names.append(f"gf_{downdip}_{alongstrike}")

    unit_src_df["disp_file"] = disp_file
    unit_src_df["outdir_path"] = outdir_paths
    unit_src_df["job_name"] = job_names

    #ii = 0
    for ii in range(100, len(unit_src_df)):
        cmd = launch(script, unit_src_df["disp_file"][ii], jagurs_template, unit_src_df["outdir_path"][ii], unit_src_df["job_name"][ii], 4, 2)
    
