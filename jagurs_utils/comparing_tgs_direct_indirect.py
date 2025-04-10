"""
Ryan Pranantyo
EOS, March 2025

a script to compare time series between 
tgs (direct output from JAGURS) and the one extracted from *nc file

this is to confirm that this script: /home/ignatius.pranantyo/apps/rf_eos_ntu/tsunami_stochastic/main_extract_elev_timeseries.py extract the correct values
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

virtual_gauges = [1, 18, 84, 138]

### direct output from JAGURS
#tgs_f1_path = "/home/ignatius.pranantyo/Tsunamis/SensitivityTests__MiodelResolutions/2006_SouthernJava/Java_lowRes__LSWE_to_LSWE"i
tgs_f1_path = "/home/ignatius.pranantyo/Tsunamis/SensitivityTests__ModelResolutions/2006_SouthernJava/Java_lowRes__LSWE_to_LSWE"
params_dict = {"step=" : "", 
               "t=" : "," , 
               "hz=" : "," , 
               "fx=" : "," , 
               "fy=" : ","}

df_f1 = pd.DataFrame()
for tgs in virtual_gauges:
    ff = Path(os.path.join(tgs_f1_path, f"tgs_station.{tgs:06d}"))
    fname = ff.name
    with open(ff, "r") as file:
        data = file.read()
        for i,j in params_dict.items():
            data = data.replace(i,j)
    fout = os.path.join(tgs_f1_path, f"{fname}.csv")
    with open(fout, "w") as file:
        file.write(data)

    df = pd.read_csv(fout, skiprows=1, header=None)
    df = df.rename(columns={0:'step', 1:'time', 2:tgs, 3:'fx', 4:'fy'})
    df = df[['time', tgs]]
    df_f1 = pd.concat([df_f1, df], axis=1)

df_f1 = df_f1.loc[:,~df_f1.columns.duplicated()].copy()


### indirect output extracted from SD01.nc file
tgs_f2 = "/home/ignatius.pranantyo/Tsunamis/SensitivityTests__ModelResolutions/2006_SouthernJava/Java_lowRes__LSWE_to_LSWE/extracted__SD01/timeseries__at_points__open_bc_line__pangandaran.csv"


###
df_f2 = pd.read_csv(tgs_f2)
df_f2 = df_f2[df_fd["time", "pts_1", "pts_18", "pts_84", "pts_138"]]


### start to compare
fig = plt.figure(figsize = (8,8), constrained_layout = True)
for ii, tgs in enumerate(virtual_gauges):
    ax = fig.add_subplot(2,2, ii+1)
    ax.set_title(tgs)
    ax.plot(df_f1[tgs], c='red', label="direct")
    ax.plot(df_f2[f"pts_{tgs}"][1:]/10000, c='blue', label="indirect")
ax.legend()
fig.savefig("comparison.png", dpi=300)
plt.close()

