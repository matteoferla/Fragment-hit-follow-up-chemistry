#!/bin/bash

# ========================
# Run $JOB_INNER_SCRIPT on a singularity container
# if JOB_INIT_SCRIPT is specified it will run that first.
# for example a ssh connection script!
# ========================

export HOST=${HOST:-$(hostname)}
export USER=${USER:-$(users)}
export HOME=${HOME:-$_CONDOR_SCRATCH_DIR}
source /etc/os-release;
echo "Running script ${0} as $USER in $HOST which runs $PRETTY_NAME."
# ---------------------------------------------------------------

export APPTAINER_CONTAINER=${APPTAINER_CONTAINER:-/data/xchem-fragalysis/shared/singularity/cuda113ubuntu20cudnn8_latest.sif}
export APPTAINER_BIND='/data/xchem-fragalysis/:/data/xchem-fragalysis/'
export APPTAINER_HOSTNAME=${APPTAINER_HOSTNAME:-singularity}
echo "Singularity $APPTAINER_CONTAINER with target script: $JOB_INNER_SCRIPT"
export APPTAINERENV_HOME2=$HOME2
export APPTAINERENV_DATA=$DATA
export APPTAINERENV_CONDA_ENVS_PATH=${APPTAINERENV_CONDA_ENVS_PATH:-$CONDA_ENVS_PATH}
export APPTAINERENV_JUPYTER_CONFIG_DIR=${APPTAINERENV_JUPYTER_CONFIG_DIR:-$JUPYTER_CONFIG_DIR}
export APPTAINER_WORKDIR=${APPTAINER_WORKDIR:-/tmp}
#export APPTAINER_WRITABLE_TMPFS=${APPTAINER_WRITABLE_TMPFS:-true}
export SSH_FORWARD_PORT=$JOB_PORT;
export APPTAINERENV_JUPYTER_PORT=${APPTAINERENV_JUPYTER_PORT:-$JOB_PORT}
# ---------------------------------------------------------------

if [ -n "$JOB_INIT_SCRIPT" ]; then
    bash $JOB_INIT_SCRIPT &
fi

# ---------------------------------------------------------------

echo 'Running singularity ...'
if ls /dev | grep -q '^nvidia'; then
    echo "NVIDIA GPU is present."
    /usr/bin/singularity exec --nv --writable-tmpfs $APPTAINER_CONTAINER /bin/bash $JOB_INNER_SCRIPT;
else
    echo "No NVIDIA GPU found."
    /usr/bin/singularity exec --writable-tmpfs $APPTAINER_CONTAINER /bin/bash $JOB_INNER_SCRIPT;
fi


echo 'DIED!' 1>&2;
