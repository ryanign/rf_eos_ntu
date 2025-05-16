"""
Ryan Pranantyo
EOS, May 2025

a quick script to extract sfincs results
"""
import os, sys
import pandas as pd
from joblib import Parallel, delayed
from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--correction_by", type=str,
        default="-0.0")
parser.add_argument("--resolution_m", type=int,
        default=90)
args = parser.parse_args()

events_catalogue = '/home/ignatius.pranantyo/Tsunamis/AOGS2025__StressCaseScenarios/EventsCatalogue__Sample__AOGS2025__20250513.csv'
events_df = pd.read_csv(events_catalogue)

sfincs_script = '/home/ignatius.pranantyo/apps/rf_eos_ntu/sfincs_utils/sfincs_check_result.py'

#sfincs_target__path = '/home/ignatius.pranantyo/Tsunamis/AOGS2025__StressCaseScenarios/SFINCS_onshore_simulations/SIMULATIONS__AOGS2025/'
correction = args.correction_by
resolution2extract = args.resolution_m

sfincs_target__path = Path(f'/home/ignatius.pranantyo/Tsunamis/AOGS2025__StressCaseScenarios/SFINCS_onshore_simulations/SIMULATIONS__AOGS2025/DeltaDTM__corrected_by_{correction}m__resolution_{resolution2extract}m')


### start to launch simulations
def extract(script, sfincs_model):
    print(f'  goint to extract {sfincs_model}')
    cmd = f'python -W ignore {script} --sfincs_model {sfincs_model} --steps 0'
    #print(cmd)
    os.system(cmd)

for ii in events_df.index:
    event_id = events_df['EVENT_ID'][ii]
    event_name = events_df['EVENT_NAME'][ii]
    scenario = f"{event_id:08d}__{event_name}"
    print(scenario)

    openbc__path = events_df['openbc_timeseries_path'][ii]
    tiles2run = pd.read_csv(events_df['tiles2run'][ii], header = None)[0].to_numpy()

    sfincs_models = []
    for tile in tiles2run:
        sfincs_models.append(os.path.join(sfincs_target__path, scenario, f'{tile}__{scenario}__resolution__{resolution2extract}_m'))

    Parallel(n_jobs = 3)(delayed(extract)(sfincs_script, sfincs_model) for sfincs_model in sfincs_models)


