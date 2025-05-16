"""
Ryan Pranantyo
EOS, May 2025

just a basic script for plotting GR from a given a- & b- values
also to count N of events to build X-years of catalogue
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

a_val = 5.99
b_val = 1.15
Mw_bins = np.arange(7.0, 8.8, 0.1)

###  1-year
N_1_year = 10**( a_val - b_val * Mw_bins) 


### 100,000 years
N_100k_year = N_1_year * 100000

fig = plt.figure(constrained_layout = True)
ax = fig.add_subplot(1,1,1)
ax.set_title(f'Number of events from: a={a_val} b={b_val}')

df = pd.DataFrame()
df['Mw'] = Mw_bins
for year in [1, 10, 50, 100, 1000, 5000, 10000, 50000, 100000]:
    df[f'N_{year}_year'] = (N_1_year * year).astype(int)
    ax.plot(df['Mw'], df[f'N_{year}_year'], label = year)
ax.legend()
ax.set_xticks(Mw_bins)
ax.grid()
ax.set_xlabel('Mw')
ax.set_ylabel('Num of events (log)')
ax.set_yscale('log')
#ax.set_xlim(8.0, 8.8)

fig.savefig(f'GR_curve_basic.png')
plt.close()
