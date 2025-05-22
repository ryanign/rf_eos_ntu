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

### find the closest epi_df lon, lat to etas_df
etas_df["closest_lon"] = np.zeros
etas_df["closest_lat"] = np.zeros
etas_df["target_Mw"] = np.zeros
for ii in etas_df.index:
    target_x = etas_df["LON"][ii]
    target_y = etas_df["LAT"][ii]

    distances = haversine_dist(target_x, target_y, epi_df["lon"], epi_df["lat"])
    min_idx = np.argmin(distances)

    etas_df.loc[ii, "closest_lon"] = epi_df["lon"][min_idx]
    etas_df.loc[ii, "closest_lat"] = epi_df["lat"][min_idx]
    etas_df.loc[ii, "target_Mw"] = np.round(etas_df['Mw'][ii], decimals = 1)

print(etas_df)

