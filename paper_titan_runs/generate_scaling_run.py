# Generates directories, input files and batch scripts for a cooldown, i.e. a series of dca/analysis
# runs with gradually decreased temperature.
# By default, this script and the provided template input files (input_sp.json.in and
# input_sp.json.in) are configured for DCA(+) calculations of the 2D single-band Hubbard model with
# on-site Coulomb interaction U and fixed density d.
#
# Usage: 1. Configure the EDIT block.
#        2. Execute the script
#
# Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)

import os
from sys import argv

################################################################################
##################################### EDIT #####################################
################################################################################

# Numbers of compute nodes
times = ["2:00"]
node_counts = [1]

# Measurements
measurements_per_node = False # Set to False for strong scaling.

# System: options: "titan" "summit"
system= "titan"

# Code version: options: "old" "new"
version= "new"

################################################################################
################################## UNTIL HERE ##################################
################################################################################

################################################################################
############################ System Specific settings ##########################
################################################################################

########## Try to optimize parameters depending on system and version #########
walkers = 5
accumulators = 1
shared = False
smt = 1

if len(argv) >= 3:
    walkers = int(argv[1])
    accumulators = int(argv[2])
    shared = walkers == accumulators
    threads = walkers + accumulators if not shared else walkers
    if threads > 14 : smt = 4
    elif threads > 7 : smt = 2
    else : smt = 1


if len(argv) >= 4 : version = argv[3]

if len(argv) < 3 : 
    print("Usage: ./generate.py <walkers> <accumulators> [version]")

threads = 8
smt = 4

measurements = 1000 #if version == "old" else 10000          # Something like 10 to 20M for strong scaling.

################################################################################

# Inverse temperature
beta = 50

# DCA cluster (see cluster_definitions.txt)
size_x = 6
size_y = 6

################################################################################
################################################################################
################################################################################

def boolToStr(boolean):
    return "true"  if boolean else "false"

# Prepares the input file for the given temperature.
def prepare_input_file(filename, outname, nodes):
    file = open(filename, 'r')
    text = file.read()
    file.close()

    text = text.replace('BETA_CURRENT_TEMP', str(beta))
    text = text.replace('NODES', str(nodes))
    text = text.replace('SIZE_X', str(size_x))
    text = text.replace('SIZE_Y', str(size_y))

    tot_meas = measurements * nodes if measurements_per_node else measurements
    text = text.replace('MEASUREMENTS', str(int(tot_meas)))

    text = text.replace('ACCUMULATORS', str(accumulators))
    text = text.replace('WALKERS', str(walkers))
    text = text.replace('SHARED', boolToStr(shared))
    text = text.replace('THREADS', str(threads))

    file = open(outname, 'w')
    file.write(text)
    file.close()

################################################################################

batch_tmpl = system + '_scaling.job.in'
tag="w" + str(walkers) + "_a" + str(accumulators)

print('Generating directories and input files:')

run_type = "weak" if measurements_per_node else "strong"
parent_dir = "run_dca_" + system + "_" + version + "_" + run_type

os.system('mkdir -p ' + parent_dir)

for n_ind, n in enumerate(node_counts):
    print('N = ' + str(n))

    # Create directory.
    dir_str = parent_dir 

    # Generate the input file.
    input = 'input_' +  tag + '.json'
    prepare_input_file('input.json.in', dir_str+ "/" + input, n)

    # Add job.
    # Generate the dca batch script.
    batch_name_dca = parent_dir + "/job_" + tag +  ".job"
    print('\nGenerating the dca batch script: ' + batch_name_dca)
    file = open(batch_tmpl, 'r')
    text_dca = file.read()
    file.close()
    executable =  '../main_dca_' + version
    text_dca = text_dca.replace('EXEC', executable)
    text_dca = text_dca.replace('INPUT', input)
    text_dca = text_dca.replace('OUTPUT', 'output_' + tag + ".txt")
    text_dca = text_dca.replace('DIRECTORY', '')

    text_dca = text_dca.replace('TIME', times[n_ind])
    text_dca = text_dca.replace('SMT', str(smt))

    text_dca = text_dca.replace('NODES', str(n))
    processes_per_node = 6 if system == "summit" else 1
    text_dca = text_dca.replace('PROCESSES', str(processes_per_node * n))
    text_dca = text_dca.replace('JOBID', tag)

    file = open(batch_name_dca, 'w')
    file.write(text_dca)
    file.close()
    os.chmod(batch_name_dca, 0o755)


    submission_exec = "bsub" if system == "summit" else "qsub"
    #submission_script += submission_exec + " " + "dca_nodes_" + str(n) +  ".job\n"



