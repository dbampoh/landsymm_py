import os
from math import log10

# Setup
# -----
# Node list
NODES = [55, 64, ]
NODE_IO = 55
PARTITION_IO = 'ivyshort'
# Walltime
WT_DAYS = 0
WT_HOURS = 12
WT_MINUTES = 0
WT_SECONDS = 0
# Job name
JOBNAME = 'lctr'
# Gridlist file
FNAME_GRIDLIST = 'lu_1700_luh2_aggregate_sum2x2_midpoint_nourban_gcp2022_2022_08_10.txt.NOHEADER'
# TEST? True WILL generate the scripts but WON'T submit jobs
TEST = False

# Keal HPC information
# --------------------
NODES_IN_PARTITION = {
                  'ivy':[55,56,57,59,61,62,64,65,66,67,68,69,71,72,73],
                  'rome':[84,85,86,87],
                  'cclake':[34,35,36,37],
                 }

CPUS_PER_NODE = {
                  'ivy':40,
                  'rome':64,
                  'cclake':80,
                 }


def get_nlines(fname_gridlist):

    with open(fname_gridlist) as f:
        nlines = sum(1 for line in f)

    return nlines


def get_walltime(days, hours, minutes, seconds):

    wt = str(days) + '-' \
        + str(hours).rjust(2,'0') + ":" \
        + str(minutes).rjust(2,'0') + ":" \
        + str(seconds).rjust(2,'0')

    return wt


def get_partition():

    partitions = []
    for node in NODES:
        partitions.extend(k for k in NODES_IN_PARTITION if node in NODES_IN_PARTITION[k])
    partitions = set(partitions)

    if len(partitions) > 1:
        exit('Nodes belong to different partitions.')

    return partitions.pop()


def get_nodelist(ntasks, partition):

    end = ntasks//CPUS_PER_NODE[partition] + 1
    nodelist = f"kea[{','.join([str(node) for node in NODES[:end]])}]"

    return nodelist


def get_ntasks(partition, nlines_gridlist):

    ntasks = len(NODES)*CPUS_PER_NODE[partition]
    ntasks = min(ntasks, nlines_gridlist)

    return ntasks


def array_py2sh(array_py, array_name):

    if len(array_py) == 0:
        return array_name + '=( )'

    array_sh = array_name + '=('
    last = array_py[-1]
    for elem in array_py:
        array_sh += str(elem) + ' '
    array_sh = array_sh[:-1] + ')'

    return array_sh


def get_starts_nlines(ntasks, nlines_gridlist):

    nlines_per_task = nlines_gridlist//ntasks
    nlines_left = nlines_gridlist%ntasks
    nlines = [nlines_per_task + int(i<nlines_left) for i in range(ntasks)]
    ends = [ sum(nlines[:i+1]) for i in range(ntasks) ]
    starts = [ i - j for i, j in zip(ends, nlines) ]

    nlines = array_py2sh(nlines, 'nlines')
    starts = array_py2sh(starts, 'starts')

    return starts, nlines


def generate_mktable_script(partition, nodelist, ntasks, starts, nlines):

    walltime = get_walltime(WT_DAYS, WT_HOURS, WT_MINUTES, WT_SECONDS)

    # Header
    lines = [
    '#!/bin/bash',
    '#SBATCH --partition=' + partition,
    '#SBATCH --nodelist=' + nodelist,
    '#SBATCH --ntasks=' + str(ntasks),
    '#SBATCH --output=' + 'log_o.txt',
    '#SBATCH --error=' + 'log_e.txt',
    '#SBATCH --time=' + walltime,
    '#SBATCH --job-name=' + JOBNAME,
    '#SBATCH --exclusive',
    '',

    # Load python 3.7
    'module load app/python/anaconda-3.7',
    '',

    # Define starts and nlines arrays
    starts,
    nlines,
    '',

    # Launch jobs
    'for r in {' + '0'*int(log10(ntasks)+1) + '..' + str(ntasks-1) + '}',
    'do',
    #'  echo "python hildap_tables.py ${starts[10#$r]} ${nlines[10#$r]} $r > out$r.txt 2>&1 &"',
    '  python hildap_tables.py ${starts[10#$r]} ${nlines[10#$r]} $r > out$r.txt 2>&1 &',
    'done',
    'wait',
    ]

    #print('\n'.join(lines))
    with open('submit_mktable.sh', 'w') as f:
        f.write('\n'.join(lines))


def generate_append_script(ntasks):
    
    nzeros = int(log10(ntasks))
    lines = [
    '#!/bin/bash',
    '#SBATCH --partition=' + PARTITION_IO,
    '#SBATCH --nodelist=kea['+str(NODE_IO)+']',
    '#SBATCH --ntasks=1',
    '#SBATCH --output=' + 'log_ao.txt',
    '#SBATCH --error=' + 'log_ae.txt',
    '#SBATCH --time=0-04:00:00',
    '#SBATCH --job-name=' + JOBNAME + '_app',
    '#SBATCH --exclusive',
    '',

    #'echo "cat hildaplus_gross_1900_2019.txt.00 > hildaplus_gross_1900_2019.txt"',
    'mkdir -p outfiles',
    'rm -f outfiles/*',
    'cat hildaplus_gross_1901_2020.txt.00 > hildaplus_gross_1901_2020.txt',
    'cat hildaplus_netfrac_1901_2020.txt.00 > hildaplus_netfrac_1901_2020.txt',
    'cat hildaplus_netfrac_check_1901_2020.txt.00 > hildaplus_netfrac_check_1901_2020.txt',
    'cat hildaplus_forestfrac_1901_2020.txt.00 > hildaplus_forestfrac_1901_2020.txt',
    'cat hildaplus_forestfrac_check_1901_2020.txt.00 > hildaplus_forestfrac_check_1901_2020.txt',
    'for r in {' + '0'*nzeros + '1..' + str(ntasks-1) + '}',
    'do',
    #'  echo "cat hildaplus_gross_1900_2019.txt.$r | sed 1d >> hildaplus_gross_1900_2019.txt"',
    '  cat hildaplus_gross_1901_2020.txt.$r | sed 1d >> hildaplus_gross_1901_2020.txt',
    '  cat hildaplus_netfrac_1901_2020.txt.$r | sed 1d >> hildaplus_netfrac_1901_2020.txt',
    '  cat hildaplus_netfrac_check_1901_2020.txt.$r | sed 1d >> hildaplus_netfrac_check_1901_2020.txt',
    '  cat hildaplus_forestfrac_1901_2020.txt.$r | sed 1d >> hildaplus_forestfrac_1901_2020.txt',
    '  cat hildaplus_forestfrac_check_1901_2020.txt.$r | sed 1d >> hildaplus_forestfrac_check_1901_2020.txt',
    'done',
    f'mv hildaplus_gross_1901_2020.txt.{"?"*(nzeros+1)} outfiles',
    f'mv hildaplus_netfrac_1901_2020.txt.{"?"*(nzeros+1)} outfiles',
    f'mv hildaplus_netfrac_check_1901_2020.txt.{"?"*(nzeros+1)} outfiles',
    f'mv hildaplus_forestfrac_1901_2020.txt.{"?"*(nzeros+1)} outfiles',
    f'mv hildaplus_forestfrac_check_1901_2020.txt.{"?"*(nzeros+1)} outfiles',
    f'mv out{"?"*(nzeros+1)}.txt outfiles',
    ]

    #print('\n'.join(lines))
    with open('submit_append.sh', 'w') as f:
        f.write('\n'.join(lines))


def generate_submit_script():

    lines = [
    "SLURM_JOBID=$(sbatch submit_mktable.sh)",
    "HPTR_ID=$(echo ${SLURM_JOBID} | sed 's,Submitted batch job ,,g')",
    "sbatch --dependency=afterok:${HPTR_ID} submit_append.sh",
    ]

    with open('submit.sh', 'w') as f:
        f.write('\n'.join(lines))


def main():

    partition = get_partition()
    nlines_gridlist = get_nlines(FNAME_GRIDLIST)
    ntasks = get_ntasks(partition, nlines_gridlist)
    starts, nlines = get_starts_nlines(ntasks, nlines_gridlist)
    nodelist = get_nodelist(ntasks, partition)

    generate_mktable_script(partition, nodelist, ntasks, starts, nlines)
    generate_append_script(ntasks)
    generate_submit_script()
    
    if not TEST:
        os.system('chmod +x submit_mktable.sh')
        os.system('chmod +x submit_append.sh')
        os.system('chmod +x submit.sh')
        os.system('./submit.sh')


if __name__ == '__main__':

    main()
    
