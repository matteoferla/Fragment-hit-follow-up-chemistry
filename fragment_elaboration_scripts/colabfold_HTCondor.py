"""
This is a script to make it easy run ColabFold (AlphaFold2 + MMSeqs2) on HTCondor.

.. code-block:: python

    import alphacondor
    print(  list(alphacondor.get_available('cuda'))  )

    barnase: str = alphacondor.retrieve_uniprot_sequence('P00648')
    barnstar: str = alphacondor.retrieve_uniprot_sequence('P11540')
    alphacondor.submit_prediction(name ='barnase-barnstar',
                              sequence = f'{barnase}:{barnstar}',
                              n_cpu=24,
                              n_gpu=1,
                              recycle=16,
                             )


The script might be out of date.
EDIT. Last run I did:

.. code-block:: bash

    colabfold_batch $HOME2/ASAP/A71/DENND5B_RAB6A_953a9.a3m $HOME2/ASAP/A71/model \
            --amber --num-recycle 16 --num-models 5 --pair-mode 'unpaired_paired' --msa-mode 'mmseqs2_uniref' \
            --rank auto --data $HOME2/.cache/colabfold

So the values for the keys do need to be updated.
Oddly, it was not making a A3M file, so I run on Colab, downloaded the A3M file, killed it, and running in the cluster.

For mega large proteins I had used these but I am not sure which is the one that makes it work:

.. code-block:: bash

    unset TF_FORCE_UNIFIED_MEMORY;
    unset XLA_PYTHON_CLIENT_MEM_FRACTION;

    export TF_FORCE_GPU_ALLOW_GROWTH=true
    export TF_FORCE_UNIFIED_MEMORY=1
    export XLA_PYTHON_CLIENT_PREALLOCATE=true
    export XLA_PYTHON_CLIENT_MEM_FRACTION=0.5
    export XLA_PYTHON_CLIENT_ALLOCATOR=platform
    export XLA_FLAGS='--xla_gpu_force_compilation_parallelism=1'
"""


import functools
import os
import re, json
import subprocess
from types import FunctionType
from typing import Dict, Any, Optional, List, Tuple
import csv, os, re
import subprocess
import pandas as pd

import requests

root_folder = os.getcwd()
logs_folder = os.path.join(root_folder, 'logs')
protein_folder = os.path.join(root_folder, 'protein_modelling')
# os.getcwd()
results_folder = os.path.join(protein_folder, 'results')
temp_folder = os.path.join(protein_folder, 'temp')


def run_command(cmd, json_output: bool = True):
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    out, err = p.communicate()
    return json.loads(out.decode()) if json_output else out.decode()


def get_available(subset='') -> List[str]:
    machines = get_status_table()
    return machines.loc[(machines.State == 'Unclaimed') & \
                        (~machines.duplicated('Name')) & \
                        (machines.Machine.str.contains(subset))
                        ] \
        .Machine.to_list()


def get_status() -> List[Dict[str, Any]]:
    return run_command('condor_status -json', json_output=True)


def get_status_table() -> pd.DataFrame:
    """
    Get a pandas DataFrame of the machines present as given by ``condor_status``.
    These extra columns are added:

    * The column ``TotalMemory_GB`` is the ``TotalMemory`` but in Gigabytes...
    * The column ``is_compute`` is if its a compute node
    * The column ``number`` is its number in the group
    * The column ``group`` is the group

    This function has been patched with the extra attributes ``specs_cols``,
    which can be used to see the specs, and ``grouping_cols``, which are the aforemention grouping columns.
    """
    machines = pd.DataFrame(get_status())
    machines['TotalMemory_GB'] = machines['TotalMemory'] / 2 ** 10
    machines['is_compute'] = machines.Machine.str.contains('pulsar-exec-node')
    machines['number'] = machines.Machine.str.extract(r'.*(\d+)').astype(int)
    machines['group'] = machines.Machine.str.replace(r'exec-node-(\d+)', r'exec-node-unlabelled-\1', regex=True) \
        .str.extract(r'exec-node-?([\D]+)\d+')[0] \
        .str.replace(r'^-', '', regex=True) \
        .str.replace(r'-$', '', regex=True)
    addresses = pd.DataFrame(machines['AddressV1'].apply(split_address).to_list())
    machines = pd.concat([machines, addresses], axis='columns')
    machines['TotalCpus'] = machines.TotalCpus.astype(int)
    return machines


setattr(get_status_table, 'grouping_cols', ['Machine', 'is_compute', 'number', 'group'])

setattr(get_status_table, 'specs_cols',
        ['Machine', 'TotalCpus', 'TotalMemory_GB', 'TotalDisk', 'CUDACapability', 'CUDADeviceName']
        )


def split_address(address):
    """
    This simply splits up `AddressV1`.
    """
    info = {}
    for term in re.findall('p="(.*?)"', address):
        m = re.search(r'"' + term + r'"; (.*?); \]', address).group(1).split('; ')
        _ = [mm.split('=') for mm in m]
        info = {**info, **{f'{term}_{k}': v for k, v in _}}
    return info


def retrieve_uniprot_data(uniprot: str) -> dict:
    response = requests.get(f'https://rest.uniprot.org/uniprotkb/{uniprot}.json')
    response.raise_for_status()
    return response.json()


def retrieve_uniprot_sequence(uniprot: str) -> str:
    response = requests.get(f'https://rest.uniprot.org/uniprotkb/{uniprot}.fasta')
    response.raise_for_status()
    return ''.join(response.text.split('\n')[1:])


def submit_prediction(name: str, sequence: str, use_existant=True, **options) -> str:
    folder = os.path.join(results_folder, name)
    if not os.path.exists(folder):
        os.mkdir(folder)
    if use_existant and os.path.exists(os.path.join(folder, f'{name}.a3m')):
        filepath = os.path.join(folder, f'{name}.a3m')
    else:
        # write input file
        csv_filename = f'{name}.csv'
        filepath = os.path.join(folder, csv_filename)
        with open(filepath, 'w') as fh:
            w = csv.writer(fh)
            w.writerow(['id', 'sequence'])
            w.writerow([name, sequence])
    # create colabfold shell script
    # recycle, models, multimer
    colab_options = filter_options(options, create_colabfold_script)
    colab_options['multimer'] = colab_options.get('multimer', ':' in sequence)
    script_path: str = create_colabfold_script(filepath, folder, **colab_options)
    # machine, n_cpu, n_mem, n_gpu
    condor_options = filter_options(options, create_condor)
    condor_path = create_condor(script_path, **condor_options)
    out = run_command(f'condor_submit {condor_path}', json_output=False)
    return re.search(r'\s(\d+)\.$', out.strip()).group(1)


def filter_options(options: Dict[str, Any], fun: FunctionType) -> Dict[str, Any]:
    keys = set(options.keys()).intersection(fun.__annotations__.keys())
    return {k: options[k] for k in keys if k in options}


def create_colabfold_script(target_path: str,
                            out_folder: str,
                            recycle: int = 20, models: int = 5,
                            multimer: bool = True,
                            msa_mode: str = 'mmseqs2_uniref',
                            cpu: bool = False) -> str:
    assert os.path.exists(target_path), f'no {target_path}'
    name = os.path.splitext(os.path.split(target_path)[-1])[0]
    cmd = f'''

    conda activate af2023
    export PATH='/data/xchem-fragalysis/mferla/.conda/envs/af2023/bin/:$PATH'
    echo $CONDA_PREFIX
    export TF_FORCE_GPU_ALLOW_GROWTH=true
    export TF_FORCE_UNIFIED_MEMORY=1
    export XLA_PYTHON_CLIENT_PREALLOCATE=true
    export XLA_PYTHON_CLIENT_MEM_FRACTION=0.5
    export XLA_PYTHON_CLIENT_ALLOCATOR=platform
    export XLA_FLAGS='--xla_gpu_force_compilation_parallelism=1';
    colabfold_batch {'--cpu' if cpu else ''} --amber --num-recycle {recycle} --num-models {models} \
    --pair-mode 'unpaired_paired' \
    --msa-mode '{msa_mode}' \
    --data {protein_folder}/params \
    --model-type {'alphafold2_multimer_v3' if multimer else 'alphafold2'} \
    --rank {'multimer' if multimer else 'ptmscore'}\
    {target_path} {out_folder}
    '''.replace('\n    ', '\n')
    return create_script(cmd=cmd, name=name, out_folder=out_folder)


model_types = ['auto', 'alphafold2', 'alphafold2_ptm', 'alphafold2_multimer_v1', 'alphafold2_multimer_v2',
               'alphafold2_multimer_v3']


def create_script(cmd: str, name: str = 'job', out_folder='temp') -> str:
    """
    This function creates a script that requires one argument when run, the logfile prefix
    """
    cmd = f'''
    cd {root_folder}
    echo {name} >> $1'.log'
    top -b -u mferla -d 300.0 >> $1'.usage' &
    echo '###################################'
    hostname
    grep -c ^processor /proc/cpuinfo
    nvidia-smi
    nvcc --version
    echo 'cuda version:'$(cat /usr/local/cuda/version.txt)
    echo 'CUDA_VISIBLE_DEVICES:'$CUDA_VISIBLE_DEVICES
    echo $1
    echo '###################################'


    source {root_folder}/.bashrc
    echo {name}

    {cmd}
    curl -X POST -H 'Content-type: application/json' --data '{{"text":"job complete: {name}"}}' $SLACK_KEY

    0
    '''.replace('\n    ', '\n')
    filepath = os.path.join(out_folder, name + '.sh')
    with open(filepath, 'w') as fh:
        fh.write(cmd)
    return filepath


def create_condor(filepath: str,
                  machine: Optional[int] = None,
                  n_cpu: Optional[int] = None,
                  n_mem: Optional[int] = None,
                  n_gpu: Optional[int] = None):
    outpath = filepath.replace('.sh', '.condor')
    out_folder = os.path.dirname(outpath)
    with open(outpath, 'w') as fh:
        fh.write('Executable      = /usr/bin/bash\n')
        fh.write(f'Arguments       = {filepath} {out_folder}/$(Cluster).$(Process)\n')
        fh.write('Universe        = vanilla\n')
        fh.write(f'Output          = {out_folder}/$(Cluster).$(Process).out\n')
        fh.write(f'Error           = {out_folder}/$(Cluster).$(Process).err\n')
        fh.write(f'Log             = {out_folder}/$(Cluster).$(Process).log\n')
        fh.write('notify_user = matteo.ferla@stats.ox.ac.uk\n')
        fh.write('notification = Complete\n')
        fh.write(f'request_cpus = {n_cpu if isinstance(n_cpu, int) else "TotalSlotCpus"}\n')
        if n_gpu:
            fh.write(f'request_gpus = {n_gpu if isinstance(n_gpu, int) else "TotalGpus"}\n')
        # fh.write(f'request_memory = {n_mem if n_mem else "TotalSlotMemory"}\n')
        fh.write('+RequiresWholeMachine = True\n')
        if machine:
            fh.write(f'requirements = (TARGET.Machine == "{machine}")\n')
        fh.write('Queue\n')
    return outpath


def get_logs(job_id) -> Dict[str, str]:
    """
    This is no longer needed...
    Returns a dictionary of the logs in ``alphacondor.logs_folder``,
    for the ``job_id`` ID.

    Keys will be 'out', 'err', 'log', 'usage'
    """
    log_filenames = *filter(functools.partial(re.match, job_id), os.listdir(logs_folder)),
    logs = {}
    for fn in log_filenames:
        with open(os.path.join(logs_folder, fn)) as fh:
            ext = os.path.splitext(fn)[1][1:]
            logs[ext] = fh.read()
    return logs
