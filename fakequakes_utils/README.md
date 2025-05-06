Requirements:
- clone MudPy from https://github.com/UO-Geophysics/MudPy.git
- I created a new conda env called mudpy
-- in mudpy, I installed obspy, mpich, mpi4py, geod, pyproj
- go to src/fk
-- make clean
-- make
- update your .bashrc file:
-- export PATH="/path_to_MudPy/MudPy/src/fk:$PATH"
-- export PYTHONPATH=/path_to_MudPy/MudPy/src/python:$PYTHONPATH"
-- export MUD="/path_to_MudPy/MudPy"

To execute:

Utilities file:



