"""
Ryan Pranantyo
EOS, May 2025

to launch events selection
"""
import os, sys
import glob
from joblib import Parallel, delayed
from joblib.externals.loky import get_reusable_executor

ncpus = 6
script = 'select_realistic_events.py'

#SFFM_files = ['/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/OUTPUTS__Slab2__Jawa/SourceCombinations__20250523/stochastic_slips__SLAB2__Jawa/stochastic_sources__Mw_8.300000__Lon_107.687250__Lat_-8.436890__table.csv',
#        '/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/OUTPUTS__Slab2__Jawa/SourceCombinations__20250523/stochastic_slips__SLAB2__Jawa/stochastic_sources__Mw_7.900000__Lon_110.214210__Lat_-9.425500__table.csv',
#        '/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/OUTPUTS__Slab2__Jawa/SourceCombinations__20250523/stochastic_slips__SLAB2__Jawa/stochastic_sources__Mw_8.400000__Lon_104.659960__Lat_-5.777790__table.csv',
#        '/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/OUTPUTS__Slab2__Jawa/SourceCombinations__20250523/stochastic_slips__SLAB2__Jawa/stochastic_sources__Mw_8.700000__Lon_118.918960__Lat_-9.669250__table.csv']

SFFM_files = glob.glob(os.path.join('/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/OUTPUTS__Slab2__Jawa/SourceCombinations__20250523/stochastic_slips__SLAB2__Jawa', 'stochastic_sources__*table.csv'))

def launch(script, SFFM_file, sourcename, grid_source, plot=True):
    cmd = f"time python -W ignore {script} --SFFM_file {SFFM_file} --sourcename {sourcename} --grid_source {grid_source} --SFFM_plot {plot}"
    print(cmd)
    os.system(cmd)

sourcename = "SLAB2__Jawa"
grid_source = "/home/ignatius.pranantyo/Tsunamis/Stochastic__Sumatera_Java/PUSGEN2017__Segmentatations/OUTPUTS__Slab2__Jawa/unit_source_grid/SLAB2__Jawa.shp"
plot = False #True

Parallel(n_jobs = ncpus)(delayed(launch)(script, infile, sourcename, grid_source, plot) for infile in SFFM_files)








get_reusable_executor().shutdown(wait=True)
