"""
Ryan Pranantyo
EOS, May 2025

Giuseppe provides a list of events catalogue from his 2D ETAS method. The catalogue contains:
    - time, mag, lon, lat (epicentre only)

There are two options:
    (1) generate random stochastic slip models following the catalogue from Giuseppe, 1 event = N realisations
    (2) I have generated some random slip models by giving Mw and epicentre. The epicentre given is centroid of 
        subfaults used and Mw is in 0.1 increament. However, there is a chance that we will update the fault grids
        and redo generating the unit sources.

    - Option (1) would be the easier one. But when we need to update the catalogue for several reasons,
      I would need to regenerate the random slip model.
    - Option (2) I can generate as many random slip as I could, then select the Mw and epicentre from the database.
      At the same time, I can also select which random slips that are realistic.
      The main drawback of option (2) is if I had to regenerate the fault grids, then we need to redo the unit
      sources generation and redo the random slip generation.

    22 May 2025:
    I AM GOING TO USE OPTION (2)
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmcrameri.cm as cm
from pathlib import Path

### FUNCTIONS
R_EARTH = 6371.0 # Earth radius in km

def haversine_dist(lon1, lat1, lon2, lat2):
    """Calculate haversine distance between two (lon, lat) points."""
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    distance = R_EARTH * c
    return distance


###############################################################################

### MAIN CODE
# list of 2D ETAS Catalogue
etas_catalogue_f = Path("/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/input_files__SouthernJava/earthquake_catalogue__region/20250520__cat_6.5-8_100k.dat")
etas_df = pd.read_fwf(etas_catalogue_f, header = None)
etas_df = etas_df.rename(columns = {0 : 'TIME', 1 : 'Mw', 2 : 'LON', 3 : 'LAT', 4 : 'NN'})

# list of centroids used to generate random slip models
random_slip_epicentre_f = Path("/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/input_files/SLAB2__Jawa__Centroids__reduced.csv")
epi_df = pd.read_csv(random_slip_epicentre_f)


# list of SFFM tables ready2use
SFFM_path = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/OUTPUTS__Slab2__Jawa/SourceCombinations__20250516/stochastic_slips__SLAB2__Jawa/SFFM_tables__collection/sffm_read2use"
SFFM_f_fmt = "filtered__N5280__stochastic_sources__Mw_{0:.6f}__table.csv"

### find the closest epi_df lon, lat to etas_df
etas_df["closest_lon"] = np.zeros
etas_df["closest_lat"] = np.zeros
etas_df["target_Mw"] = etas_df["Mw"].round(decimals=1)
etas_df["EVENTID"] = np.zeros
etas_df["SFFM_filename"] = np.zeros
etas_df["SFFM_src_id"] = np.zeros

### just for testing
etas_df = etas_df[etas_df['Mw'] > 6.95]





print(etas_df)

### assign SFFM to use
print(f"Looking for SFFM realisation to use")
for Mw in np.sort(etas_df["target_Mw"].unique()):
    print(f" target Mw = {Mw}")
    # temp filtering
    etas_Mw = etas_df[etas_df["target_Mw"] == Mw].copy()
    etas_Mw["SFFM_filename"] = np.zeros
    etas_Mw["SFFM_src_id"] = np.zeros
    
    SFFM_f = os.path.join(SFFM_path, SFFM_f_fmt.format(Mw))
    df = pd.read_csv(SFFM_f, low_memory = False)
    #sffm_df = df.iloc[:-6]
    #info_df = df.iloc[-6:]
    
    for ii in etas_Mw.index:
        ### to find the closest lon and lat from SFFM generated
        target_lon = etas_df["LON"][ii]
        target_lat = etas_df["LAT"][ii]
        distances = haversine_dist(target_lon, target_lat, epi_df["lon"], epi_df["lat"])
        min_idx = np.argmin(distances)
        etas_Mw.loc[ii, "closest_lon"] = epi_df["lon"][min_idx]
        etas_Mw.loc[ii, "closest_lat"] = epi_df["lat"][min_idx]
        #etas_df.loc[ii, "target_Mw"] = np.round(etas_df['Mw'][ii], decimals = 1)

    
        ### start to find SFFM to use
        target_x = etas_Mw["closest_lon"][ii]
        target_y = etas_Mw["closest_lat"][ii]

        sffm_sample_id = df.columns[np.isin(df.values, [str(target_x), str(target_y)]).any(0)].tolist()
        #print(sffm_sample_id)

        if len(sffm_sample_id) >= 1:
            ##THERE IS A CHANCE THAT I WILL TAKE THE SAME sffm_sample_id,
            ##mainly because of the target_Mw, target_x, and target_y are same values
            ##to avoid this, I remove sffm_sample_id[0] after it has been assigned to the first
            ##   duplicate event
            ##
            ##Let's say:
            ## - event_1 and event_2 have same target_Mw, target_x, and target_y; and
            ## - sffm_sample_id = ['sffm_89', 'sffm_111', 'sffm_345']
            ##Hence:
            ## - event_1 will get 'sffm_89' -- then 'sffm_89' will be removed from df.columns
            ## - event_2 will get 'sffm_111'
            ##
            ##BUT, if there are three events but we only have two sffm_sample_id:
            ## the 3rd event will not be assigned sffm_model. Therefore, we would need another set of SFFM 
            ##
            ##I think, the option is to generate another set of SFFM then merge with the last set seperately

            etas_Mw.loc[ii, "SFFM_filename"] = SFFM_f
            etas_Mw.loc[ii, "SFFM_src_id"  ] = sffm_sample_id[0]
            df = df.drop(columns=[sffm_sample_id[0]])

            ### giving EVENTID
            ntime = f"{etas_Mw['TIME'][ii].astype(int):010d}"
            nMw_tmp = etas_Mw["target_Mw"][ii].astype(str).split('.')
            nMw = '_'.join(nMw_tmp)
            nLon_tmp = etas_Mw["closest_lon"][ii].astype(str).split('.')
            nLon = '_'.join(nLon_tmp)
            nLat_tmp = etas_Mw["closest_lat"][ii].astype(str).split('.')
            nLat = '_'.join(nLat_tmp)
            sample_id = sffm_sample_id[0].split('__')[-1]
            EVENTID = f"{ntime}__Mw_{nMw}__Sample_{sample_id}__Lon_{nLon}__Lat_{nLat}"
            etas_Mw.loc[ii, "EVENTID"] = EVENTID

    ### remove np.zeros SFFM_filename and SFFM_src_id
    etas_Mw = etas_Mw[etas_Mw["SFFM_src_id"] != np.zeros]
    ### remove duplications at SFFM_filename and SFFM_src_id
    ### means, I need more realisation
    etas_Mw = etas_Mw.drop_duplicates(subset = ["SFFM_filename", "SFFM_src_id"], keep="first")

    etas_df.loc[etas_Mw.index, "SFFM_filename"] = etas_Mw.loc[etas_Mw.index, "SFFM_filename"]
    etas_df.loc[etas_Mw.index, "SFFM_src_id"] = etas_Mw.loc[etas_Mw.index, "SFFM_src_id"]


print(f"\n===============================")
print(f"NEED MORE REALISASTION SFFM ...")
print(etas_df[etas_df["SFFM_filename"] == np.zeros])
print(f"===============================\n")

print(f"\n===============================")
print(f"READY TO USE ...")
print(etas_df[etas_df["SFFM_filename"] != np.zeros])
print(f"===============================\n")


###
# RESERVED FOR SAVING THE OUTPUT
###
sys.exit()
