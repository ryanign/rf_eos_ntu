"""
Ryan Pranantyo
EOS, March 2025

a basic script to generate tsunami offshore footprint,
could be used to generate a new timeseries as well.
"""
import os, sys
import numpy as np
import pandas as pd
import rioxarray
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import argparse
import glob

from pathlib import Path
from joblib import Parallel, delayed

def offshore_footprint(args):
    print(args)
    
    fin = Path(args.SFFM_model)
    model_name = fin.name[:-4]
    df = pd.read_csv(fin)
    outdir = Path(os.path.join(args.where_to_save, f"footprints_{model_name}"))
    outdir.mkdir(exist_ok = True)

    num_samples = args.Nsamples

    ### find a template SD01.nc
    temp_folder = glob.glob(os.path.join(args.tsunami_gf_path, "gf_*"))[0]
    dtemp = xr.open_dataset(os.path.join(temp_folder, "SD01.nc"))
    if "step" in dtemp: dtemp = dtemp.drops_vars(["step"])
    if "max_velocity" in dtemp: dtemp = dtemp.drops_vars(["max_velocity"])

    for nn in range(num_samples):
        cols = ['unit_source_index', 'unit_source_filename', f'unit_source_slip__{nn}']
        dfn = df[cols]
        dfn = dfn.loc[(dfn[f'unit_source_slip__{nn}'] != 0)]
        print(dfn)
        
        
        ### read slip needed from each unit source
        ### do G*m = d
        initial_disp = np.zeros_like(dtemp.initial_displacement)
        max_height = np.zeros_like(dtemp.max_height)
        wave_height = np.zeros_like(dtemp.wave_height)
        for idx in dfn.index:
            unit_src = Path(dfn['unit_source_filename'][idx]).name[:-4].split("_")
            downdip = int(unit_src[-2])
            alongstrike = int(unit_src[-1])
            slip = dfn[f"unit_source_slip__{nn}"][idx].astype("float16")
            print(downdip, alongstrike, slip)
            
            ### there is likely that the tsunami greens function has not been generated
            ### read tsunami greens function footprint
            gf_path = os.path.join(args.tsunami_gf_path, f"gf__{downdip}_{alongstrike}")
            if os.path.exists(gf_path):
                ### here I read only result from highest grid
                gf_01 = xr.open_dataset(os.path.join(gf_path, "SD01.nc"))
                if "step" in gf_01: gf_01 = gf_01.drop_vars(["step"])
                if "max_velocity" in gf_01: gf_01 = gf_01.drop_vars(["max_velocity"])
                gf_01 = gf_01.fillna(0)
                gf_01 = gf_01 * slip                                                #### G(i) * m(i)
                initial_disp += gf_01.initial_displacement
                max_height += gf_01.max_height
                wave_height =+ gf_01.wave_height
            
                gf_01.close()
                status = "generate_footprint"
            else:
                print(f" >> INCOMPLETE UNIT SOURCE GF FOR SAMPLE NUMBER {nn}")
                print(f" >> ERROR = TSUNAMI GF FOR DOWNDIP {downdip} ALONG STIRKE {alongstrike} HAS NOT BEEN GENERATED!")
                status = "do_not_generate_footprint"

        if status == "generate_footprint":
            dcp = dtemp.copy()  
            dcp.initial_displacement.data = initial_disp
            dcp.max_height.data = max_height
            dcp.wave_height.data = wave_height
        
            ### check result throug figure
            fig = plt.figure(figsize=(10,5), constrained_layout = True)
            ax1 = fig.add_subplot(1,2,1, projection=ccrs.PlateCarree())
            ax1.coastlines()
            disp = ax1.imshow(initial_disp[0], 
                    extent = [dtemp.lon.min(), dtemp.lon.max(), dtemp.lat.min(), dtemp.lat.max()],
                    vmin = -1, vmax = 1, cmap = "RdBu_r",
                    origin="lower")
            fig.colorbar(disp, orientation="horizontal", extend="both")
            ax2 = fig.add_subplot(1,2,2, projection=ccrs.PlateCarree())
            ax2.coastlines()
            fpmax = ax2.imshow(max_height,
                    extent = [dtemp.lon.min(), dtemp.lon.max(), dtemp.lat.min(), dtemp.lat.max()],
                    vmin = 0, vmax = 2, cmap = "hot_r",
                    origin="lower")
            fig.colorbar(fpmax, orientation="horizontal", extend="max")
            fout = os.path.join(outdir, f"footprints__sample_num__{nn}.png")
            fig.savefig(fout)
            plt.close()
        
            ### export a new netcdf, aggregation
            fout = os.path.join(outdir, f"footprints__sample_num__{nn}.nc")
            dcp.to_netcdf(fout)
            dcp.close()

    dtemp.close()



    return df, dfn, dcp

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--SFFM_model", type=str,
            default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/OUTPUTS/stochastic_slips__sumatera_jawa__slab2__edited/stochastic_sources__Mw_7.900000__Lon_104.104800__Lat_-7.122200__table_simplified.csv",
            help = "a list of stochastic finite fault model generated from automate_stochastic_generation.py")
    parser.add_argument("--Nsamples", type=int,
            default = 25,
            help = "number of footprints sample to be generated, max is number of stochastic samples per-SFFM model file")
    parser.add_argument("--tsunami_gf_path", type=str,
            default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/jagurs_runs/uji_coba/unit_sources_gf__tsunami/tsunami_gf__2700m",
            help = "path to tsunami Green's Function to be used")
    parser.add_argument("--where_to_save", type=str,
            default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/jagurs_runs/uji_coba/tsunami_footprints/sample_jawa__2700m_to_0900m")
    args = parser.parse_args()

    main_outdir = Path(args.where_to_save)
    main_outdir.mkdir(exist_ok = True)

    df, dfn, dcp = offshore_footprint(args)
