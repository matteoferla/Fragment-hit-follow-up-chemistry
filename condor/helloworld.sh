#!/bin/bash

export HOST=${HOST:-$(hostname)}
export USER=${USER:-$(users)}
export HOME=${HOME:-$_CONDOR_SCRATCH_DIR}
source /etc/os-release;


echo "************************"
echo "HELLO WORLD!"
echo "************************"
echo "Greet from script ${0} as $USER in $HOST which runs $PRETTY_NAME"


echo "ls $PWD"
ls $PWD

export HOME=${HOME:-$_CONDOR_SCRATCH_DIR}

echo "ls $HOME"
ls $HOME

echo "printenv"

printenv