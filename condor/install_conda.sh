if [[ -z "$CONDA_PREFIX" ]]; then
    echo "Must provide CONDA_PREFIX in environment" 1>&2
    exit 1
fi

export DATA=/data/xchem-fragalysis
rm -r $CONDA_PREFIX
bash $DATA/shared/Miniconda3-latest-Linux-x86_64.sh -p $CONDA_PREFIX -b

source $CONDA_PREFIX/etc/profile.d/conda.sh
export PIP_NO_CACHE_DIR=1
export PIP_NO_USER=1
export PYTHONUSERBASE=$CONDA_PREFIX
conda activate base
conda update -n base -y -c defaults conda

# Jupyter stuff
conda install -y -n base -c conda-forge distro nodejs sqlite jupyterlab jupyter_http_over_ws nb_conda_kernels
conda update -y -c conda-forge nodejs   # peace among worlds

# install whatever you want here
pip install -q pandas plotly seaborn
# pip install fragmenstein nglview pebble pandarallel pandera pyrosetta_help
# pip install $DATA/shared/pyrosetta-2023.27+release.e3ce6ea9faf-cp311-cp311-linux_x86_64.whl

conda clean --all -y