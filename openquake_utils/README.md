# UTILITIES TOOLS TO PREPARE INPUT AND POST-PROCESS OUTPUT FILES FOR THE OPENQUAKE ENGINE

## Introduction
This sub-repo is collection of my scripts to prepare input and post-process output files for the OpenQuake Engine.

## How to install
### On local computer
1. I installed OQ engine version 3.25 following [on Linux](https://docs.openquake.org/oq-engine/manual/latest/getting-started/installation-instructions/universal.html#user-installation)
2. You will have a dedicated openquake environment, it is installed at '~/openquake/'

### On WildFly (NTU's HPC)
1. I created an openquake environment using Anaconda.
2. Installed OQ engine version 3.25 following [on Linux](https://docs.openquake.org/oq-engine/manual/latest/getting-started/installation-instructions/universal.html#user-installation)
3. You will have a dedicated openquake evnironment, it is installed at '~/openquake/' and different than the openquake env under Anaconda. 

>[!NOTE]
> The openquake environment under Anaconda was only used to tackle sqlite3 library which is not available on WildFly. If there is no issue with sqlite3, you can follow installation process On local computer.

## Basic steps
1. Prepare all input files required in one folder. Read the documentation for the details (too many information).
2. Activate openquake environment.
3. Inside the working folder, do 'oq engine --run job.ini --log_file run.log'
4. Once finished, check simulation ID of the simulation inside run.log (at the very last line). This indicate simulation output filename inside '~/oqdata/'
5. Export results to a bunch of csv files using 'oq engine --export-outputs ID TARGET_FOLDER'
6. Post-analysis.

### Typical output files
Inside my job.ini file, I set up:
'''
intensity_measure_types_and_levels = {
  "PGA": logscale(0.005, 2.13, 45),
  "SA(1.0)": logscale(0.005, 2.13, 45)
  
quantiles = 0.25 0.50 0.95
hazard_maps = true
uniform_hazard_spectra = true
poes = 0.1 0.02
'''
