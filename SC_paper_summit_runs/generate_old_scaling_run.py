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

################################################################################
##################################### EDIT #####################################
################################################################################

# Numbers of compute nodes
times = ["2:00", "18:00"]
node_counts = [10, 4608] # run with 10 nodes is for testing.

# Measurements
measurements_per_node = False # Set to False for strong scaling.
measurements = 80e6           # Something like 10 to 20M for strong scaling.

# System: options: "titan" "summit"
system= "summit"

# Code version: options: "old" "new"
version= "old"

################################################################################
################################## UNTIL HERE ##################################
################################################################################

################################################################################
############################ System Specific settings ##########################
################################################################################

########## Try to optimize parameters depending on system and version #########
walkers = 0
accumulators = 0
shared = False
smt = 1
threads = 7 if system == "summit" else 8


if system == "summit" and version == "new" :
    walkers = accumulators = 7
    shared = True

elif system == "summit" and version == "old" :
    walkers = 1
    accumulators = 13
    smt = 2

elif system == "titan" and version == "new" :
    walkers = 5
    accumulators = 1

elif system == "titan" and version == "old" :
    walkers = 3
    accumulators = 11

else : raise "Invalid options."

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

print('Generating directories and input files:')

run_type = "weak" if measurements_per_node else "strong"
parent_dir = "run_dca_" + system + "_" + version + "_" + run_type
submission_script = "#!/bin/bash\n\n"

for n_ind, n in enumerate(node_counts):
    print('N = ' + str(n))

    # Create directory.
    dir_str = parent_dir + '/nodes_' + str(n)
    os.system('mkdir -p ' + dir_str)

    # Generate the input file.
    input = dir_str + '/input.json'
    prepare_input_file('input.json.in', input, n)

    # Add job.
    # Generate the dca batch script.
    batch_name_dca = parent_dir + "/dca_nodes_" + str(n) +  ".job"
    print('\nGenerating the dca batch script: ' + batch_name_dca)
    file = open(batch_tmpl, 'r')
    text_dca = file.read()
    file.close()
    executable =  '../main_dca_' + version
    text_dca = text_dca.replace('EXEC', executable)
    text_dca = text_dca.replace('INPUT', './nodes_' + str(n) + "/input.json")
    text_dca = text_dca.replace('OUTPUT', './nodes_' + str(n) + '/output')

    text_dca = text_dca.replace('TIME', times[n_ind])
    text_dca = text_dca.replace('SMT', str(smt))

    text_dca = text_dca.replace('NODES', str(n))
    processes_per_node = 6 if system == "summit" else 1
    text_dca = text_dca.replace('PROCESSES', str(processes_per_node * n))

    file = open(batch_name_dca, 'w')
    file.write(text_dca)
    file.close()
    os.chmod(batch_name_dca, 0o755)


    submission_exec = "bsub" if system == "summit" else "qsub"
    submission_script += submission_exec + " " + "dca_nodes_" + str(n) +  ".job\n"


file = open(parent_dir + "/submit_all.bash", "w")
file.write(submission_script)
file.close()
