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

def find_closest_coordinates(df, ref_df, order=0):
    ### find the closest lon and lat from SFFM database
    ### df = etas_df, catalogue to be filled up
    ### ref_def = epi_df, database of epicentres used to generate SFFM
    ### order = 0, 1, 2, there are chances that we're running out SFFM because of the use of the closest epicentre from database,
    ###  in the case, we will relook at the second or the third closest database.
    for ii in df.index:
        target_lon = df["LON"][ii]
        target_lat = df["LAT"][ii]
        distances = haversine_dist(target_lon, target_lat, ref_df["lon"], ref_df["lat"])
        min_idx = np.argpartition(distances, order)[order]
        df.loc[ii, "closest_lon"] = ref_df["lon"][min_idx]
        df.loc[ii, "closest_lat"] = ref_df["lat"][min_idx]
    return df

def find_sffm_to_use(SFFM_path, SFFM_f_fmt, df, bg_df, Mw = 7.3):
    ##FIRST NOTE -----
    ## df = etas_df to be looked up for the SFFM filename and SFFM src id
    ## bg_df = 'background df' to check whether SFFM filename and SFFM src id has been used previously
    ##   if it has been used, we will not use it
    ## for the first iteration, df = bg_df
    
    ##SECOND NOTE -----
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

    print(Mw)
    etas_Mw = df[df["target_Mw"] == Mw].copy()
    unique_coords = set(zip(etas_Mw['closest_lon'], etas_Mw['closest_lat']))
    for jj, pts in enumerate(unique_coords):
        etas_xx_yy = etas_Mw[etas_Mw['closest_lon'] == pts[0]]
        etas_xx_yy = etas_Mw[etas_Mw['closest_lat'] == pts[1]]

        ### find SFFM file:
        sffm_f = os.path.join(SFFM_path, SFFM_f_fmt.format(Mw, pts[0], pts[1]))
        if os.path.exists(sffm_f) == True:
            sffm_df = pd.read_csv(sffm_f)
            sffm_samples_id = sffm_df.columns[1:-1]
            
            ### check if it has been used by other events
            bg_f = bg_df[bg_df['SFFM_filename'] == sffm_f]
            if len(bg_f) > 0:
                bg_sample = bg_f['SFFM_src_id']
                bg_sample = bg_sample[bg_sample != np.zeros]
                sffm_samples_id = sffm_samples_id.drop(bg_sample)

            if len(sffm_samples_id) > 0:
                for kk in etas_xx_yy.index:
                    df.loc[kk, "SFFM_filename"] = sffm_f
                    if len(sffm_samples_id) > 0:
                        df.loc[kk, "SFFM_src_id"] = sffm_samples_id[0]
                        sffm_samples_id = sffm_samples_id.drop(sffm_samples_id[0])
                    #else:
                    #    df.loc[jj, "SFFM_src_id"] = "need more realisastion"
    return df

###############################################################################

### MAIN CODE
# list of 2D ETAS Catalogue
#etas_catalogue_f = Path("/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/input_files__SouthernJava/earthquake_catalogue__region/20250526__cat_6.5-8_100k.dat")
etas_catalogue_f = Path("/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/input_files__SouthernJava/earthquake_catalogue__region/20250603__cat_6.5-8.7_100k__7-CATALOGUES.dat")
etas_df = pd.read_fwf(etas_catalogue_f, header = None)
etas_df = etas_df.rename(columns = {0 : 'TIME', 1 : 'Mw', 2 : 'LON', 3 : 'LAT', 4 : 'NN'})

# list of centroids used to generate random slip models
#random_slip_epicentre_f = Path("/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/input_files/SLAB2__Jawa__Centroids__reduced.csv")
random_slip_epicentre_f = Path("/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/input_files/SLAB2__Jawa__Centroids__reduced__20250526__used.csv")
epi_df = pd.read_csv(random_slip_epicentre_f)


# list of SFFM tables ready2use
#SFFM_path = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/OUTPUTS__Slab2__Jawa/SourceCombinations__20250516/stochastic_slips__SLAB2__Jawa/SFFM_tables__collection/sffm_read2use"
#SFFM_f_fmt = "filtered__N5280__stochastic_sources__Mw_{0:.6f}__table.csv"

SFFM_path = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/OUTPUTS__Slab2__Jawa/SourceCombinations__20250523/SFFM_realistic__SLAB2__Jawa"
SFFM_f_fmt = "realistic__stochastic_sources__Mw_{0:.6f}__Lon_{1:.6f}__Lat_{2:.6f}__table.csv"

### find the closest epi_df lon, lat to etas_df
etas_df["closest_lon"] = np.zeros
etas_df["closest_lat"] = np.zeros
etas_df["target_Mw"] = etas_df["Mw"].round(decimals=1)
etas_df["EVENTID"] = np.zeros
etas_df["SFFM_filename"] = np.zeros
etas_df["SFFM_src_id"] = np.zeros

### just for testing
Mw_threshold = 6.95
etas_df = etas_df[etas_df['Mw'] > Mw_threshold]
print(etas_df)

"""
I will group the mapping based on (i) target_Mw and (ii) closest_lon and closest_lat
"""

### initiate the first round of look up the SFFM to use
print(f"\n first round ...")
etas_df = find_closest_coordinates(etas_df, epi_df, order=0)

for Mw in np.sort(etas_df['target_Mw'].unique()):
    etas_df = find_sffm_to_use(SFFM_path, SFFM_f_fmt, etas_df, etas_df, Mw)

### if need second round
if len(etas_df[etas_df['SFFM_src_id'] == np.zeros]) > 0:
    print(f"\n second round ...")
    etas2nd_df = find_closest_coordinates(
                    etas_df[etas_df['SFFM_src_id'] == np.zeros],
                    epi_df,
                    1)
    etas_df = etas_df[etas_df['SFFM_src_id'] != np.zeros]
    bg_df = etas_df
    for Mw in np.sort(etas2nd_df['target_Mw'].unique()):
        etas2nd_df = find_sffm_to_use(SFFM_path, SFFM_f_fmt, etas2nd_df, bg_df, Mw)
    #etas_df = etas_df[etas_df['SFFM_src_id'] != np.zeros]
    round_2nd = True
else:
    round_2nd = False


### FROM BELOW IS NOT THAT EFFECTIVE!
### MIGHT NEED TO CHANGE!


### if need 3rd round
if len(etas2nd_df[etas2nd_df['SFFM_src_id'] == np.zeros]) > 0:
    print(f"\n third round ...")
    etas3rd_df = find_closest_coordinates(
                    etas2nd_df[etas2nd_df['SFFM_src_id'] == np.zeros],
                    epi_df,
                    2)
    etas2nd_df = etas2nd_df[etas2nd_df['SFFM_src_id'] != np.zeros]
    bg_df = pd.concat([etas_df, etas2nd_df])
    for Mw in np.sort(etas3rd_df['target_Mw'].unique()):
        etas3rd_df = find_sffm_to_use(SFFM_path, SFFM_f_fmt, etas3rd_df, bg_df, Mw)
    #etas2nd_df = etas2nd_df[etas2nd_df['SFFM_src_id'] != np.zeros]
    round_3rd = True
else:
    round_3rd = False

### if need 4th round
if round_3rd == True and len(etas3rd_df[etas3rd_df['SFFM_src_id'] == np.zeros]) > 0:
    print(f"\n fourth round ...")
    etas4th_df = find_closest_coordinates(
                    etas3rd_df[etas3rd_df['SFFM_src_id'] == np.zeros],
                    epi_df,
                    3)
    etas3rd_df = etas3rd_df[etas3rd_df['SFFM_src_id'] != np.zeros]
    bg_df = pd.concat([etas_df, etas2nd_df, etas3rd_df])
    for Mw in np.sort(etas4th_df['target_Mw'].unique()):
        etas4th_df = find_sffm_to_use(SFFM_path, SFFM_f_fmt, etas4th_df, bg_df, Mw)
    #etas3rd_df = etas3rd_df[etas3rd_df['SFFM_src_id'] != np.zeros]
    round_4th = True
else:
    round_4th = False

### if need 5th round
if round_4th == True and len(etas4th_df[etas4th_df['SFFM_src_id'] == np.zeros]) > 0:
    print(f"\n fifth round ...")
    etas5th_df = find_closest_coordinates(
                    etas4th_df[etas4th_df['SFFM_src_id'] == np.zeros],
                    epi_df,
                    4)
    etas4th_df = etas4th_df[etas4th_df['SFFM_src_id'] != np.zeros]
    bg_df = pd.concat([etas_df, etas2nd_df, etas3rd_df, etas4th_df])
    for Mw in np.sort(etas5th_df['target_Mw'].unique()):
        etas5th_df = find_sffm_to_use(SFFM_path, SFFM_f_fmt, etas5th_df, bg_df, Mw)
    #etas4th_df = etas4th_df[etas4th_df['SFFM_src_id'] != np.zeros]
    round_5th = True
else:
    round_5th = False

### if need 6th round
if round_5th == True and len(etas5th_df[etas5th_df['SFFM_src_id'] == np.zeros]) > 0:
    print(f"\n sixth round ...")
    etas6th_df = find_closest_coordinates(
                    etas5th_df[etas5th_df['SFFM_src_id'] == np.zeros],
                    epi_df,
                    5)
    etas5th_df = etas5th_df[etas5th_df['SFFM_src_id'] != np.zeros]
    bg_df = pd.concat([etas_df, etas2nd_df, etas3rd_df, etas4th_df, etas5th_df])
    for Mw in np.sort(etas6th_df['target_Mw'].unique()):
        etas6th_df = find_sffm_to_use(SFFM_path, SFFM_f_fmt, etas6th_df, bg_df, Mw)
    #etas5th_df = etas5th_df[etas5th_df['SFFM_src_id'] != np.zeros]
    round_6th = True
else:
    round_6th = False


### check after round 6th if we still have empty SFFM_filename and SFFM_src_id
### potentially they come from small Mw events
if round_6th == True:
    check_df = etas6th_df[etas6th_df['SFFM_src_id'] != np.zeros]
    if len(check_df) == 0:
        print(f'events below will be removed as we would need more realisastion events and it is tricky!')
        print(etas6th_df)
        etas6th_df = etas6th_df[etas6th_df['SFFM_src_id'] != np.zeros]


### merging all back into one DF
collect_df = etas_df[etas_df['SFFM_src_id'] != np.zeros]
if round_2nd == True and round_3rd == False and round_4th == False and round_5th == False and round_6th == False:
    collect_df = pd.concat([collect_df, etas2nd_df])
elif round_2nd == True and round_3rd == True and round_4th == False and round_5th == False and round_6th == False:
    collect_df = pd.concat([collect_df, etas2nd_df, etas3rd_df])
elif round_2nd == True and round_3rd == True and round_4th == True and round_5th == False and round_6th == False:
    collect_df = pd.concat([collect_df, etas2nd_df, etas3rd_df, etas4th_df])
elif round_2nd == True and round_3rd == True and round_4th == True and round_5th == True and round_6th == False:
    collect_df = pd.concat([collect_df, etas2nd_df, etas3rd_df, etas4th_df, etas5th_df])
elif round_2nd == True and round_3rd == True and round_4th == True and round_5th == True and round_6th == True:
    collect_df = pd.concat([collect_df, etas2nd_df, etas3rd_df, etas4th_df, etas5th_df, etas6th_df])

### cleaning up if there are still empty SFFM_filename and SFFM_src_id
collect_df = collect_df[collect_df['SFFM_src_id'] != np.zeros]

### saving
where_to_save = etas_catalogue_f.parent
fname = f'{etas_catalogue_f.name}__EVENT_LIST__Mw{Mw_threshold}+.csv'
fout = os.path.join(where_to_save, fname)
collect_df.to_csv(fout, index = False)

print(fout)
print(collect_df)























sys.exit()
for Mw in np.sort(etas_df['target_Mw'].unique()):
    print(Mw)
    etas_Mw = etas_df[etas_df["target_Mw"] == Mw].copy()

    unique_coords = set(zip(etas_Mw['closest_lon'], etas_Mw['closest_lat']))
    for jj, pts in enumerate(unique_coords):
        #print(xx, yy)
        etas_xx_yy = etas_Mw[etas_Mw['closest_lon'] == pts[0]]
        etas_xx_yy = etas_Mw[etas_Mw['closest_lat'] == pts[1]]
        #print(etas_xx_yy)

        ### find SFFM file:
        sffm_f = os.path.join(SFFM_path, SFFM_f_fmt.format(Mw, pts[0], pts[1]))
        if os.path.exists(sffm_f) == False:
            print('not available')
            etas_df.loc[jj, "SFFM_filename"] = np.zeros
            etas_df.loc[jj, "SFFM_src_id"] = np.zeros
        else:
            sffm_df = pd.read_csv(sffm_f)
            sffm_samples_id = sffm_df.columns[1:-1]

            for jj in etas_xx_yy.index:
                etas_df.loc[jj, "SFFM_filename"] = sffm_f
                if len(sffm_samples_id) > 0:
                    etas_df.loc[jj, "SFFM_src_id"] = sffm_samples_id[0]
                    sffm_samples_id = sffm_samples_id.drop(sffm_samples_id[0])
                else:
                    etas_df.loc[jj, "SFFM_src_id"] = "need more realisastion"

print(etas_df)




sys.exit()
sys.exit()


sys.exit()
for ii in etas_df.index:
    ### find the closest lon and lat from SFFM database
    target_lon = etas_df["LON"][ii]
    target_lat = etas_df["LAT"][ii]
    distances = haversine_dist(target_lon, target_lat, epi_df["lon"], epi_df["lat"])
    min_idx = np.argmin(distances)
    etas_df.loc[ii, "closest_lon"] = epi_df["lon"][min_idx]
    etas_df.loc[ii, "closest_lat"] = epi_df["lat"][min_idx]

    ### find SFFM file
    target_Mw = etas_df["target_Mw"][ii]

    # check if available
    if os.path.exists(sffm_f):
        print(ii, f"available")
        etas_df.loc[ii, "SFFM_filename"] = sffm_f

        sffm_df = pd.read_csv(sffm_f)
        slip_samples = slip_df.columns[1:-1].to_list()

sys.exit()


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
