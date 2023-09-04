export HOST=${HOST:-$(hostname)}
export USER=${USER:-$(users)}
export HOME=${HOME:-$_CONDOR_SCRATCH_DIR}
source /etc/os-release;

if [ -n "$JUPYTER_PORT" ]; then
    echo "$JUPYTER_PORT set"
elif [ -n "$JOB_PORT" ]; then
    export JUPYTER_PORT=$JOB_PORT
elif [ -n "$SSH_FORWARD_PORT" ]; then
    export JUPYTER_PORT=$SSH_FORWARD_PORT
else
    raise error "Your JUPYTER_PORT is not specified"
fi

if [ -z "$JUPYTER_CONFIG_DIR" ]; then
    raise error "Your JUPYTER_CONFIG_DIR is not specified either"
fi



echo "************************"
echo "HELLO JUPYTER!"
echo "************************"
echo "Greet from Jupyter lab script ${0} as $USER in $HOST which runs $PRETTY_NAME on $JUPYTER_PORT with settings from $JUPYTER_CONFIG_DIR"

source /data/xchem-fragalysis/shared/bashrc.sh;
# conda activate

# First time? Remember to set:
# jupyter notebook --generate-config
# yes invasion | jupyter server password

# port is JUPYTER_PORT
while true
do
jupyter lab --ip="0.0.0.0" --no-browser
done
