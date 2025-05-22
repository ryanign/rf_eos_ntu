"""
Ryan Pranantyo
EOS, 28 February 2025

> main script to launch any type of simulations
  >> template_folder
  >> a displacement file, configuration of the disp grd file should as same as the dem grd file

"""
import os, sys
import numpy as np
import pandas as pd
import argparse
from pathlib import Path

def write_pbs(args):
    pbs = open("jagurs.pbs", "w")
    pbs.write(f"#PBS -N {args.job_name}\n")
    pbs.write(f"#PBS -P {args.project_code}\n")
    pbs.write(f"#PBS -q {args.qtype}\n")
    pbs.write(f"#PBS -l walltime={args.walltime:02d}:00:00\n")
    pbs.write(f"#PBS -l select=1:ncpus={args.ncpus}\n")
    pbs.write(f"cd $PBS_O_WORKDIR\n")
    pbs.write(f"module load openmpi/4.1.1\n")
    pbs.write(f"export LD_LIBRARY_PATH='LD_LIBRARY_PATH':/home/ignatius.pranantyo/apps/mylibs/proj-4.9.3/lib:/home/ignatius.pranantyo/apps/mylibs/netcdf-fortran-4.5.3/lib\n")
    pbs.write(f"time ./jagurs_serial_ncdio par=tsun.par\n")
    ### to compress output files
    pbs.write(f"source /home/ignatius.pranantyo/.bashrc\n")
    pbs.write(f"~/apps/miniconda3/bin/activate\n")
    pbs.write(f"conda activate tsunamis_py39\n")
    for nn in range(args.num_of_grids):
        #pbs.write(f"nccopy -d2 SD{nn:02d}.nc SD{nn:02d}_tmp.nc\n")
        #pbs.write(f"mv SD{nn:02d}_tmp.nc SD{nn:02d}.nc\n")
        pbs.write(f"python utils_compress_jagurs_nc.py --jagurs_nc SD{nn:02d}.nc\n")
    pbs.close()

def launch(args):
    #print(args)
    ### paths preparation
    curdir = os.getcwd()
    where_to_run = Path(args.where_to_run)
    where_to_run.mkdir(exist_ok = True)
    os.system(f"rm -rf {where_to_run}/*")
    print(f"Preparing scenario = {where_to_run}")

    ### copy jagurs template folder to where_to_run
    jagurs_template = args.jagurs_template
    cmd = f"cp -rP {jagurs_template}/* {where_to_run}"
    os.system(cmd)

    ### copy symbolic link of the displacement file to where_to_run
    os.chdir(where_to_run)
    disp_file = args.displacement_grdfile
    cmd = f"ln -s {disp_file} disp__SD00_____.grd"
    os.system(cmd)

    ### modify jagurs.pbs script
    write_pbs(args)

    ### launch job
    os.system("qsub jagurs.pbs")

    ### go back to script directory
    os.chdir(curdir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--displacement_grdfile", type=str,
                        default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/jagurs_runs/uji_coba/displacements_20250228/scenario_00001/disp__SD00____.grd",
                        help = "Displacement grd file, can be a unit source or a heterogenous slip model.")
    parser.add_argument("--jagurs_template", type=str,
                        default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/jagurs_runs/uji_coba/template_20250228",
                        help = "JAGURS template directory contains all basic input grid files")
    parser.add_argument("--where_to_run", type=str,
                        default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/jagurs_runs/uji_coba/runs_20250228/runs_00001",
                        help = "Path to run simulation")
    parser.add_argument("--job_name", type=str, default = "jagurs00001",
                        help = "just for job naming")
    parser.add_argument("--project_code", type=str, default = "eos_luca.dalzilio",
                        help = "Project code name on WildFly")
    parser.add_argument("--qtype", type=str, default = "qintel_wfly",
                        help = "queue type in Wildfly: qintel_wfly, qamd_wfly, dev")
    parser.add_argument("--walltime", type=int, default = 8,
                        help = "computation time request")
    parser.add_argument("--ncpus", type=int, default = 4,
                        help = "number of cpus requested to run this job")
    parser.add_argument("--num_of_grids", type=int, default = 2,
                        help = "nested grids used in gridfile.dat")

    args = parser.parse_args()

    launch(args)
