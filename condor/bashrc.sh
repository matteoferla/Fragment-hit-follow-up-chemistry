# ========================
# ``source /data/xchem-fragalysis/shared/bashrc.sh``
# Sets presents...
# tries to source $HOME/.bashrc
# or $HOME2/.bashrc;
# or a fallback
# ========================

# Matteo's .bashrc has
# export CONDA_PREFIX=${CONDA_PREFIX:-$HOME2/neoconda}
# source $CONDA_PREFIX/etc/profile.d/conda.sh
# conda activate
# export CONDA_ENVS_PATH=$CONDA_ENVS_PATH:$DATA/mferla/.conda/envs:$DATA/sanchezg/app/miniconda3_2/envs:$DATA/mferla/rocky-conda/envs
# export JUPYTER_CONFIG_DIR=$DATA/mferla/rocky-conda/jupyter # for jupyter
# export CONDA_ALWAYS_YES=yes;

# universal fixes
export PS1="[\u@\h \W]\$"
export LANG=en_GB.UTF-8
export DATA=/data/xchem-fragalysis
export HOST=${HOST:-$(hostname)}
export USER=${USER:-$(users)}
export HOME=${HOME:-$_CONDOR_SCRATCH_DIR}
export SHELL=/bin/bash
source /etc/os-release;
export PIP_NO_CACHE_DIR=1
export PIP_NO_USER=1
export NUMEXPR_MAX_THREADS=$(lscpu -p=CPU | tail -n 1 | xargs)
# Jedi cache
mkdir -p $HOME2/.cache
export XDG_CACHE_HOME=$HOME2/.cache

# -------------------------------------

if [ -f $HOME/.bashrc ]
then
    source $HOME/.bashrc;
elif [ -f $HOME2/.bashrc ]
then
    source $HOME2/.bashrc;
else
    echo "No .bashrc found in $HOME or $HOME2" 1>&2;
    exit 1;
fi

export PYTHONUSERBASE=$CONDA_PREFIX;

#export JUPYTER_CONFIG_DIR=${JUPYTER_CONFIG_DIR:-$DATA/mferla/rocky-conda/jupyter}

after_install() {
    # multi-user conda
    conda clean --all -y 2>&1 > /dev/null;
    chmod -r -f a+rwX $CONDA_PREFIX 2>&1 > /dev/null;
}
sleep 1;
