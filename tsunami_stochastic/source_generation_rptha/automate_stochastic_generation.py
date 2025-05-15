"""
Ryan Pranantyo
EOS, February 2025

This script is used to automate the proces of generating stochastic slip model
    main script is combine_tsunami_sources.R
    what to be automated :
        - desired_Mw = from Mw_min to Max_max with 0.1 bin
        - target_location = following centroids given (the finite fault model)
        - number_of_sffm = 100 generate N random samples
        - stochastic_slip_events_table output file name

"""
import os, sys
import numpy as np
import pandas as pd
import argparse
from joblib import Parallel, delayed
from pathlib import Path

def write_combine_source_script(args, Mw, target_lon, target_lat, num_of_sffm, scaling, where_to_save):
    print(Mw, target_lon, target_lat, num_of_sffm, scaling, where_to_save)

    final_destination = Path(os.path.join(where_to_save, f"stochastic_slips__{args.sourcename}"))
    final_destination.mkdir(exist_ok = True)

    raster_grid_list = os.path.join(final_destination, f"unit_source_grid_raster_filename_index__{args.sourcename}.csv")
    fin_name = f"stochastic_sources__Mw_{Mw:3f}__Lon_{target_lon:4f}__Lat_{target_lat:4f}"
    stochastic_table_list = os.path.join(final_destination, f"{fin_name}__table.csv")

    fin = open(f"{fin_name}.R", 'w')
    fin.write(f"library(rptha)\n")
    fin.write(f"unit_source_dirname = '{args.unit_source_dirname}'\n")
    fin.write(f"sourcename = '{args.sourcename}'\n")
    fin.write(f"all_discretized_source_RDS = '{args.all_discretized_source_RDS}'\n")
    fin.write(f"desired_Mw = {Mw}\n")
    fin.write(f"target_location = c({target_lon},{target_lat})\n")
    fin.write(f"number_of_sffm = {num_of_sffm}\n")
    fin.write(f"### end of input ###\n")
    fin.write(f"discretized_source = readRDS(all_discretized_source_RDS)[[sourcename]]\n")
    fin.write(f"discretized_source_statistics = discretized_source_approximate_summary_statistics(discretized_source)\n")
    fin.write(f"\n")
    fin.write( "unit_source_raster_files = paste(unit_source_dirname, '/', sourcename, '/',\n")
    fin.write( "    sourcename, '_', discretized_source_statistics$downdip_number, '_',\n")
    fin.write( "    discretized_source_statistics$alongstrike_number, '.tif', sep='')\n")
    fin.write( "\n")
    fin.write( "if(!(all(file.exists(unit_source_raster_files)))){stop('Could not find some unit source raster files')}\n")
    fin.write( "\n")
    fin.write(f"unit_source_rasters = lapply(as.list(unit_source_raster_files), f<-function(x) raster(x))\n")
    fin.write(f"stochastic_slip_events = sffm_make_events_on_discretized_source(\n")
    fin.write(f"        discretized_source_statistics = discretized_source_statistics,\n")
    fin.write(f"        target_location = target_location,\n")
    fin.write(f"        target_event_mw = desired_Mw,\n")
    fin.write(f"        num_events = number_of_sffm,\n")
    fin.write(f"        uniform_slip = FALSE,\n")
    fin.write(f"        expand_length_if_width_limited = 'random',\n")
    fin.write(f"        use_deterministic_LWkc=FALSE,\n")
    fin.write(f"        clip_random_parameters_at_2sd = TRUE,\n")
    fin.write(f"        relation = '{scaling}', \n")
    fin.write(f"        peak_slip_location_near_centre=FALSE)\n")
    fin.write(f"stochastic_slip_events_table = sffm_events_to_table(stochastic_slip_events)\n")
    fin.write(f"### Saving output \n")
    fin.write(f"write.csv(unit_source_raster_files, '{raster_grid_list}', row.names=TRUE)\n")
    fin.write(f"write.csv(stochastic_slip_events_table, '{stochastic_table_list}', row.names=TRUE)\n")
    fin.close()

    return f"{fin_name}", final_destination

def launch_combine_source_script(args, df_comb, Nsamples, scaling, where_to_save, ii):
    """ to launch """
    fin_name, final_destination = write_combine_source_script(
            args, df_comb['Mw'][ii], df_comb['target_lon'][ii], df_comb['target_lat'][ii], Nsamples, scaling, where_to_save)
    ### execute!
    cmd = f"Rscript {fin_name}.R"
    print(f"  {cmd}")
    os.system(cmd)

    ### cleaning up the stochastic slip table
    stoch_f = os.path.join(final_destination, f"{fin_name}__table.csv")
    df_stoch = pd.read_csv(stoch_f)
    df_stoch = df_stoch.rename(columns={'Unnamed: 0' : 'event_id'})
    df_stoch['sourcename'] = args.sourcename
    df_stoch.to_csv(stoch_f, index = False)

    if args.clean_rscripts:
        print(f"  removing all R (temporary) script to generte stochastic slip models ...")
        cmd = f"rm -f {fin_name}.R"
        os.system(cmd)
    
    print(f">>> DONE FOR {ii} <<<")
    
    return df_stoch

def combine_unit_source_and_sffm(ii, df_comb, unit_df, args):
    sffm_in = os.path.join(args.where_to_save, f"stochastic_slips__{args.sourcename}",
            f"stochastic_sources__Mw_{df_comb.Mw[ii]:3f}__Lon_{df_comb.target_lon[ii]:4f}__Lat_{df_comb.target_lat[ii]:4f}__table.csv")
    sffm_df = pd.read_csv(sffm_in)

    df = pd.DataFrame()
    for src in sffm_df.index:
        print(f"   decomposing event_id {src+1} ...")
        flt_index_str = sffm_df.event_index_string[src]
        flt_slip_str  = sffm_df.event_slip_string[src]
        flt_index = list(map(int, flt_index_str.split('-')[:-1]))
        flt_slip  = list(map(float, flt_slip_str.split('_')[:-1]))
        tmp_df = pd.DataFrame(data = {'unit_source_index' : flt_index,
                                      f'unit_source_slip__{src}'  : flt_slip})
        tmp_df2 = pd.merge(unit_df, tmp_df, on='unit_source_index', how='left')
        df = pd.concat([df, tmp_df2], axis=1)
    df = df.loc[:,~df.columns.duplicated()].copy()
    df = df.fillna(0)
    print(df)

    sffm_out = os.path.join(args.where_to_save, f"stochastic_slips__{args.sourcename}",
            f"stochastic_sources__Mw_{df_comb.Mw[ii]:3f}__Lon_{df_comb.target_lon[ii]:4f}__Lat_{df_comb.target_lat[ii]:4f}__table_simplified.csv")
    df.to_csv(sffm_out, index=False)

    return df

def main(args):
    if args.Mw_bins is None:
        Mw_series = [args.Mw_target]
    else:
        Mw_min = args.Mw_bins[0]
        Mw_max = args.Mw_bins[1]
        Mw_bin = args.Mw_bins[2]
        Mw_series = np.arange(Mw_min , Mw_max + Mw_bin, Mw_bin)
    
    print(f"Going to generate Mw = {Mw_series}")

    Nsamples = args.NumSamples
    centroids = pd.read_csv(args.epicentre_list)
    scaling = args.scaling_relationship
    where_to_save = Path(args.where_to_save)
    where_to_save.mkdir(exist_ok = True)

    ### just to create temp dataframe
    Mw_list, lons, lats = [], [], []
    for Mw in Mw_series:                ### Mw_series
        for jj in centroids.index:       #[0, 1]:           ### centroids index
            Mw_list.append(Mw)
            lons.append(centroids.lon[jj])
            lats.append(centroids.lat[jj])
    df_comb = pd.DataFrame(
            data = {'Mw' : Mw_list,
                    'target_lon' : lons,
                    'target_lat' : lats})
   
    print(df_comb)

    ### launch combine_source_script.R
    #ii = 0
    #df_stoch = launch_combine_source_script(args, df_comb, Nsamples, scaling, where_to_save, ii)
    ncpus = args.ncpus
    Parallel(n_jobs = args.ncpus)(delayed(launch_combine_source_script)(
        args, df_comb, Nsamples, scaling, where_to_save, ii) for ii in df_comb.index)

    ### cleaning up unit_source_grid_raster_filename
    unit_f = os.path.join(where_to_save, f"stochastic_slips__{args.sourcename}", 
            f"unit_source_grid_raster_filename_index__{args.sourcename}.csv")
    unit_df = pd.read_csv(unit_f)
    unit_df = unit_df.rename(columns = {'Unnamed: 0' : 'unit_source_index',
                                        'x'          : 'unit_source_filename'})
    unit_df.to_csv(unit_f, index=False)

    
    if args.combine_sffm_and_unit_source == True:
        """ WE CAN REMOVE THIS LATER """
        ### combine unit_df and df_stoch (stochastic slip model into one DF ###
        Parallel(n_jobs = args.ncpus)(delayed(combine_unit_source_and_sffm)(ii, df_comb, unit_df, args) for ii in df_comb.index)

    ### cleaning up working folder
    print(f"cleaning up working folder ...")
    if args.clean_rscripts:
        print(f"  removing all R (temporary) scripts to generte stochastic slip models ...")
        cmd = f"rm -f stochastic_sources__Mw_*__Lon_*__Lat_*.R"
        os.system(cmd)
    if args.clean_raw_sffm_tables:
        print(f"  removing raw SFFM tables: {where_to_save}/stochastic_slips__{args.sourcename}/stochastic_sources__Mw*__table.csv")
        cmd = f"rm -f {where_to_save}/stochastic_slips__{args.sourcename}/stochastic_sources__Mw*__Lon*__Lat*__table.csv"
        os.system(cmd)

    ### DONE ###
    print("MAIN PROGRAM TO GENERATE STOCHASTIC HETEROGENOUS SLIP MODEL IS DONE")
    print("  execute the next step to produce a raster file of the source")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--Mw_bins", type=float, nargs="+", default=None, #[7.0, 9.0, 0.1],
            help = "[Mw_min, Mw_max, bin], if want to generate a series of Mw")
    parser.add_argument("--Mw_target", type=float, default=7.3,
            help = "target Mw to be generated")
    parser.add_argument("--NumSamples", type=int, default=100,
            help = "Number of random samples on every epicentre")
    parser.add_argument("--epicentre_list", type=str, default="/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/input_files/centroids__sumatera_jawa__slab2__20kmX20km.csv",
            help = "list of epicentres to generate stochastic slip models")
    parser.add_argument("--unit_source_dirname", type=str,
            default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/OUTPUTS/Unit_source_data",
            help = "path to unit sources generated from the previous step")
    parser.add_argument("--sourcename", type=str, default="sumatera_jawa__slab2__edited",
            help = "name of the source zone")
    parser.add_argument("--all_discretized_source_RDS", type=str,
            default = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/OUTPUTS/all_discretized_sources.RDS",
            help = "path to all_discretized_sources.RDS generated from the the previous step")
    parser.add_argument("--scaling_relationship", type=str, default="Strasser",
            help = "scalling relationship between MW, rupture area, and slip. Read main script for more options")
    parser.add_argument("--where_to_save", type=str, default = "../OUTPUTS/",
            help = "where to save the stochastic tables")
    parser.add_argument("--ncpus", type=int, default = 2,
            help = "number of cpus to use")
    parser.add_argument("--clean_rscripts", type=bool, default=True,
            help = "delete Rscripts thaat used to generate SFFM at the end of the script")
    parser.add_argument("--clean_raw_sffm_tables", type=bool, default=True,
            help = "delete ram SFFM tables at the end of the script")
    parser.add_argument("--combine_sffm_and_unit_source", type=bool, default=False,
            help = "combine unit_source filename and SFFM table to make it simple")
    args = parser.parse_args()

    main(args)

