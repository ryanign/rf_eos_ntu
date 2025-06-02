import os, sys
import numpy as np
import pandas as pd
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--Mw", type = float, default = 8.3)
args = parser.parse_args()

Mw = args.Mw

SFFM_files = glob.glob(f'/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/OUTPUTS__Slab2__Jawa/SourceCombinations__20250523/SFFM_realistic__SLAB2__Jawa/realistic__stochastic_sources__Mw_{Mw:.6f}__*__table.csv')
SFFM_files.sort()

data = list()

Nrealization = 0
for ii, sffm_f in enumerate(SFFM_files):
    print(f'{ii+1} of {len(SFFM_files)}')
    print(f'  >> number of realistic SFFM = {Nrealization}')
    sffm_df = pd.read_csv(sffm_f)
    unit_slips = sffm_df.columns[1:-1]
    Nrealization += len(unit_slips)

    slip_data = sffm_df[unit_slips]
    data.append(slip_data)
    #data = np.hstack([data, slip_data])

data = np.hstack(data)

colname = list()
for ii in range(Nrealization):
    colname.append(f'unit_source_slip__{ii}')

df = pd.DataFrame(
                data = data,
                columns = colname
                )
df['unit_source_index'] = sffm_df['unit_source_index']
df['unit_source_filename'] = sffm_df['unit_source_filename']

cols_reoder = ['unit_source_index'] + colname + ['unit_source_filename']
df = df[cols_reoder]

fout = os.path.join('/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/OUTPUTS__Slab2__Jawa/SourceCombinations__20250523/merge__SFFM_realistic__SLAB2__Jawa', f'merged__realistic_SFFM__Mw{Mw:.6f}.csv')
df.to_csv(fout, index = False)

print(df)
