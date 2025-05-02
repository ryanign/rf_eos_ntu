"""
Ryan Pranantyo
EOS, April 2025

to launch big SFINCS simulations!

cmd = 
    python -W ignore sfincs_build_event_model.py --sfincs_template < > --target_dir < > --event_name < > --domain_name < > --bc_waterlevel < > --bc_points_csv < > --scale_ratio 1
"""
import os, sys
import numpy as np
import pandas as pd
from pathlib import Path

events_catalogue = '../EventsCatalogue__Jawa__TESTING__20250502.csv'
events_df = pd.read_csv(events_catalogue)

sfincs_script = '/home/ignatius.pranantyo/apps/rf_eos_ntu/sfincs_utils/sfincs_build_event_model.py'

bc_points_csv_f = Path('/home/ignatius.pranantyo/Tsunamis/AOGS2025__StressCaseScenarios/SFINCS_config/vectors/domains_tile/points_open_bc_line_jawa.csv')
bc_points_csv__path = bc_points_csv_f.parent

sfincs_templates__path = '/home/ignatius.pranantyo/Tsunamis/AOGS2025__StressCaseScenarios/SFINCS_onshore_templates/DeltaDTM__correctedby-0.5m__update'
resolution2run = 90

sfincs_target__path = '/home/ignatius.pranantyo/Tsunamis/AOGS2025__StressCaseScenarios/SFINCS_onshore_simulations/SIMULATIONS'

### start to launch simulations
for ii in events_df.index:
    event_id = events_df['EVENT_ID'][ii]
    event_name = events_df['EVENT_NAME'][ii]
    scenario = f"{event_id:08d}__{event_name}"
    print(scenario)

    openbc__path = events_df['openbc_timeseries_path'][ii]
    tiles2run = pd.read_csv(events_df['tiles2run'][ii], header = None)[0].to_numpy()

    for tile in tiles2run:
        print(f'  going to launch sfincs for {tile}')
        
        tile_id = tile.split('_')[-1] 
        sfincs_template = os.path.join(sfincs_templates__path, f'domain_id_{tile_id}', f'resolution__{resolution2run}_m')

        where2run_sfincs = Path(os.path.join(sfincs_target__path, scenario))
        where2run_sfincs.mkdir(exist_ok = True)
        #where2run_sfincs = Path(os.path.join(where2run_sfincs, tile))
        
        bc_points_csv = os.path.join(bc_points_csv__path, f'{bc_points_csv_f.name[:-4]}__{tile}.csv')

        bc_waterlevel = os.path.join(openbc__path, f'sfincs-open-bc__at__{tile}', 'extracted__SD01', f'timeseries__at_{bc_points_csv_f.name[:-4]}__{tile}.csv')

        ### run
        cmd = f'python -W ignore {sfincs_script} --sfincs_template {sfincs_template} --target_dir {where2run_sfincs} --event_name {scenario} --domain_name {tile} --bc_waterlevel {bc_waterlevel} --bc_points_csv {bc_points_csv} --scale_ratio 1'
        print(cmd)
        os.system(cmd)




