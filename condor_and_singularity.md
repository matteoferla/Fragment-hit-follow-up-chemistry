## Diamond / STFC cluster job submission

> Details of the cluster are not posted here. Ask for them

The cluster accepts k8s, condor and slurm jobs.
The former is limited for a regular person, while the HTCondor is very flexible and has a lot to offer.

## HTCondor

The CPU nodes and the GPU nodes are accessible via different submission nodes.
The submission nodes and the compute nodes have a data folder mounted at `/data/xchem-fragalysis` (See details files).
The folder `/data/xchem-fragalysis/shared` has a few files that may be useful,
that rely heavily on environment variables.

### target_script.condor

The submission script [target_script.condor](condor/target_script.condor) runs `$JOB_NODE_SCRIPT` 
within initial dir `$HOME2`.

Envs used:
* $HOME2 the fake home, _e.g._ /data/xchem-fragalysis/mferla
* $JOB_SCRIPT the bash script to run —which has to be on the compute node (i.e. data folder).

This script and a few other shared here use a custom notation `$JOB_*` to indicate a variable that is set for a job.

Say the data folder did not exist, the $HOME2 would need to be set to `/tmp`
and extra arguments for transfer of files would need to be added.

To specify a particular machine 
add `-a 'Requirements=(machine == "orpheus-worker-gpu-13.novalocal")'` as a cmd arg
to specify a particular node as I have been unable to do so via environment variables.

nice scripts to use for `$JOB_SCRIPT`:

* helloworld.sh
* singularity.sh

### singularity.sh

To run a calculation I create a job.sh in the folder of interest, used with the file below with singularity.sh.
The [singularity.sh](condor/singularity.sh) script is a wrapper for the singularity command,
which is basically

```bash
/usr/bin/singularity exec --writable-tmpfs  --bind $APPTAINER_BIND $APPTAINER_CONTAINER /bin/bash $JOB_INNER_SCRIPT;
```

but with a few extra envs set if not already set, both config variables 
(`$APPTAINER_CONTAINER`, `$APPTAINER_BIND`, `$APPTAINER_HOSTNAME`) and a few environment variables visible within
the container (`$APPTAINERENV_HOME2`, `$APPTAINERENV_DATA`, `$APPTAINERENV_CONDA_ENVS_PATH`,
`$APPTAINERENV_JUPYTER_CONFIG_DIR`, `$APPTAINERENV_JUPYTER_PORT` (based on `$JOB_PORT`)
—any `$APPTAINERENV_*` becomes a variable in the container).

Before the singularity container is run, the `$JOB_INIT_SCRIPT` is called (bash on background) if present.
Say, to set up a reverse port forward ssh connection to a remote machine.

In `/data/xchem-fragalysis/shared/singularity` folder are two sif files:

* rockyplus.sif
* cuda113ubuntu20plus.sif

The definition files used to create them can be found [here](condor).

Example of a calculation:

```bash
```bash
export JOB_SCRIPT=/data/xchem-fragalysis/shared/singularity.sh
export JOB_INNER_SCRIPT=/data/xchem-fragalysis/mferla/target-name/job.sh
export APPTAINER_CONTAINER=/data/xchem-fragalysis/shared/singularity/rockyplus.sif
condor_submit /data/xchem-fragalysis/shared/target_script.condor -a 'Requirements=(machine == "orpheus-worker-67.novalocal")';
```
where `target-name/job.sh` is something like:

```bash
# how this job was actually run:

export HOST=${HOST:-$(hostname)}
export USER=${USER:-$(users)}
export HOME=${HOME:-$_CONDOR_SCRATCH_DIR}
export N_CORES=$(cat /proc/cpuinfo | grep processor | wc -l);
source /etc/os-release;
echo "Running script ${0} as $USER in $HOST which runs $PRETTY_NAME"
source /data/xchem-fragalysis/mferla/activate_conda.sh;

cd /data/xchem-fragalysis/mferla/target-name;
nice -19 python fragmenstein_merge_sw_place.mod.py \
--n_cores $(($N_CORES - 1)) \
--template reference.pdb \
--suffix _iter3 \
--hits new.sdf \
--second_hits fragmented.sdf \
--weights weights.json \
--sw_databases REAL-Database-22Q1.smi.anon \
--weights weights.json \
--combination_size 2 \
--timeout 600;
```

### Set up conda

There's a conda installer in shared (if not `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh`),
and a script [install_conda.sh](condor/install_conda.sh) to copy and embellish.

#### Automatic

```bash
export JOB_SCRIPT=/data/xchem-fragalysis/shared/singularity.sh
export APPTAINER_CONTAINER=/data/xchem-fragalysis/shared/singularity/rockyplus.sif
export JOB_INNER_SCRIPT=/data/xchem-fragalysis/shared/install_conda.sh
export APPTAINERENV_CONDA_PREFIX='/data/xchem-fragalysis/👾👾/conda' # or whatever
condor_submit /data/xchem-fragalysis/shared/target_script.condor
```
Conda when activated adds the variable `CONDA_PREFIX`, which is the folder where conda lives.
Herein, I am using this ahead of time.
#### Interactive

```bash
condor_submit -interactive initialdir=$HOME2"
export APPTAINER_BIND='/data/xchem-fragalysis/:/data/xchem-fragalysis/';
export APPTAINER_CONTAINER='/data/xchem-fragalysis/shared/singularity/rockyplus.sif';
export APPTAINER_HOSTNAME=rockylinux;
export APPTAINERENV_PIP_NO_CACHE_DIR=1
export APPTAINERENV_PIP_NO_USER=1
export APPTAINERENV_PYTHONUSERBASE=$CONDA_PREFIX
/usr/bin/singularity shell --writable-tmpfs $APPTAINER_CONTAINER;
```
And copy-paste to flavour. The `$PIP_NO_CACHE_DIR`, `$PIP_NO_USER`, `$PYTHONUSERBASE` are needed to avoid errors.
As there is little or no home folder on the singularity container (16MB or 64GB). To check this, `df -h`.
These declared within are:

```bash
```bash
export PIP_NO_CACHE_DIR=1
export PIP_NO_USER=1
export PYTHONUSERBASE=$CONDA_PREFIX
```

Also, in a Jupyter notebook, the autocompletion, Jedi, will try to write to the home folder, stop it!

```bash
mkdir -p $HOME2/.cache
export XDG_CACHE_HOME=$HOME2/.cache
```

#### Other

Note that conda installed interactively gunks up your `.bashrc` file.
To avoid this, simply source conda.sh in etc/profile.d after installing conda:

```bash
source $CONDA_PREFIX/etc/profile.d/conda.sh
```

If you installed Mamba at `$CONDA_PREFIX`, then it's

```bash

source $CONDA_PREFIX/etc/profile.d/conda.sh
source $CONDA_PREFIX/etc/profile.d/mamba.sh
mamba activate
```

## Footnote

### Fluff
Fluff in a script that is useful for debug:
(The variable `$PRETTY_NAME` comes from `/etc/os-release` and is the distro name)

```bash
export HOST=${HOST:-$(hostname)}
export USER=${USER:-$(users)}
export HOME=${HOME:-$_CONDOR_SCRATCH_DIR}
source /etc/os-release;
echo "Running script ${0} as $USER in $HOST which runs $PRETTY_NAME."
```

Also for sanity

```bash
export PS1="[\u@\h \W]\$"
export LANG=en_GB.UTF-8
```

The former sets the prompt to something more useful: `[foouser@foohost foofolder]$`.