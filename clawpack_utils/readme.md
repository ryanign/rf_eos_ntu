# My utility scripts for Clawpack

Utilities scripts to perform pre- and post-simulations using Clawpack, especially using the geoclaw library.

Ryan Pranantyo
EOS, 29 September 2026

## Instalation procedures on WildFly
```
# load required modules
module load python/3.11.5
module load gnu/gcc-12.3
module load openmpi/4.1.5-gcc12.3.0

# create an virtualenv
python3 -m venv ~/apps/clawpack-env
source ~/apps/clawpack-env/bin/activate

# install basic libraries
pip install numpy
pip install matplotlib

# install clawpack
pip install meson-python ninja
pip install --src=$HOME/apps/clawpack_src --no-build-isolation -e git++https://github.com/clawpack/clawpack.git@v5.14.0#egg=clawpack
```

once finished without error, add below to activate environment and libraries required automatically when activating 'clawpack-env'

```
cat >> ~/apps/clawpack-env/bin/activate << 'EOF'

# Clawpack environment
module load python/3.11.5
module load gnu/gcc-12.3
module load openmpi/4.1.5-gcc12.3.0
export CLAW=$HOME/apps/clawpack_src/clawpack
export FC=gfortran
EOF
```

Test the installation

```
source ~/apps/clawpack-env/bin/activate
cd $CLAW/geoclaw/examples/tsunami/chile2010
python maketopo.py
make .output
make plots
```

Hope, there is no error!
