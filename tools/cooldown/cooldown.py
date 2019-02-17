# Generates directories, input files and batch scripts for a cooldown, i.e. a series of dca/analysis
# runs with gradually decreased temperature.
# By default, this script and the provided template input files (input_sp.json.in and
# input_sp.json.in) are configured for DCA(+) calculations of the 2D single-band Hubbard model with
# on-site Coulomb interaction U and fixed density d.
#
# Usage: 1. Configure the EDIT block.
#        2. Execute the script: python cooldown.py
#
# See https://github.com/CompFUSE/DCA/wiki/Running for more details on how to use this script and
# how to run a DCA(+) calculation.
#
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)

import os

################################################################################
##################################### EDIT #####################################
################################################################################

# Numbers of compute nodes
nodes = 1

# DCA cluster (see cluster_definitions.txt)
size_x = 6
size_y = 6


# DCA(+) iterations
iters_ht = 10  # first/highest temperature
iters_last = 10  # last/lowest temperature
iters = 5     # other temperatures

# Inverse temperatures for cooldown
betas = [1, 2, 5, 10, 20, 40, 50]
# Starting temperature for computing two-particle quantities and running the analysis.
beta_analysis = 50

#Measurements
measurements_per_node = False
measurements = 1000000

################################################################################
################################## UNTIL HERE ##################################
################################################################################

############################ System Specific settings ##########################
# Command (including options) to run (MPI) programs, e.g. 'jsrun ...'.
# Used to execute main_dca.
batch_tmpl = 'summit/summit.job'
run_command_dca = 'jsrun -n PROCESSES -a 1 -g 1 -c 7 -b rs'
# The run command to execute main_analysis.
#run_command_analysis = run_command_dca
################################################################################


# Formats the inverse temperature beta.
def beta_format(x):
    return ('%.6g' % x)

# Prepares the input file for the given temperature.
def prepare_input_file(filename, b_ind, type):
    file = open(filename, 'r')
    text = file.read()
    file.close()

    text = text.replace('BETA_CURRENT_TEMP', str(betas[b_ind]))
    text = text.replace('SIZE_X', str(size_x))
    text = text.replace('SIZE_Y', str(size_y))

    if measurements_per_node : meaasurements *= nodes
    text = text.replace('MEASUREMENTS', str(measurements))

    # For the first temperature set the initial self-energy to zero.
    if (b_ind == 0):
        text = text.replace('./beta_BETA_PREVIOUS_TEMP/dca_sp.hdf5', 'zero')
        text = text.replace('./beta_BETA_PREVIOUS_TEMP/configuration', '')
    else:
        previous_beta = betas[b_ind-1] if type == 'sp' else betas[b_ind]
        text = text.replace('BETA_PREVIOUS_TEMP', str(previous_beta))

    if (b_ind == 0):
        text = text.replace('ITERS', str(iters_ht))
    elif (b_ind == len(betas)-1):
        text = text.replace('ITERS', str(iters_last))
    else:
        text = text.replace('ITERS', str(iters))

    file = open(filename, 'w')
    file.write(text)
    file.close()

################################################################################

batch_str_dca_sp = ''
batch_str_dca_tp = ''
#batch_str_analysis = ''

print('Generating directories and input files:')

for b_ind, beta in enumerate(betas):
    print('BETA = ' + str(beta))

    # Create directory.
    dir_str = './beta_' + str(beta)
    cmd = 'mkdir -p ' + dir_str
    os.system(cmd)
    cmd = 'mkdir -p ' + dir_str + "/configuration"
    os.system(cmd)

    input_sp = dir_str + '/input_sp.json'
    input_tp = dir_str + '/input_tp.json'

    data_dca_sp   = dir_str + '/dca_sp.hdf5'
    data_dca_tp   = dir_str + '/dca_tp.hdf5'
    #data_analysis = dir_str + '/analysis.hdf5'

    # dca sp
    # Generate the sp input file.
    cmd = 'cp ./input_sp.json.in ' + input_sp
    os.system(cmd)
    prepare_input_file(input_sp, b_ind, 'sp')

    # Add job.
    batch_str_dca_sp += run_command_dca + ' ./main_dca ' + input_sp + '\n'

    if (beta >= beta_analysis):
        # dca tp
        # Generate the tp input file.
        cmd = 'cp ./input_tp.json.in ' + input_tp
        os.system(cmd)
        prepare_input_file(input_tp, b_ind, 'tp')

        # Add job.
        batch_str_dca_tp += run_command_dca + ' ./main_dca ' + input_tp + '\n'

        # analysis
        # Add job.
        #batch_str_analysis = batch_str_analysis + run_command_analysis + ' ./main_analysis ' + input_tp + '\n'


# Get filename extension of batch script.
_, extension = os.path.splitext(batch_tmpl)

# Generate the dca batch script.
for type in ['sp', 'tp']:
    batch_name_dca = 'dca_' + type +'.job'
    print('\nGenerating the dca batch script: ' + batch_name_dca)
    file = open(batch_tmpl, 'r')
    text_dca = file.read()
    file.close()

    text_dca = text_dca.replace('JOBS', batch_str_dca_sp if type == 'sp' else batch_str_dca_tp)

    text_dca = text_dca.replace('APPLICATION', 'dca')
    text_dca = text_dca.replace('NODES', str(nodes))
    text_dca = text_dca.replace('PROCESSES', str(6 * nodes))

    file = open(batch_name_dca, 'w')
    file.write(text_dca)
    file.close()

# Generate the analysis batch script.
#batch_name_analysis = 'job.analysis_U=' + str(U) + '_d=' + str(d) + '_Nc=' + str(Nc) + extension
#print('Generating the analysis batch script: ' + batch_name_analysis)
#file = open(batch_tmpl, 'r')
#text_analysis = file.read()
#file.close()

# text_analysis = text_analysis.replace('APPLICATION', 'analysis')
# text_analysis = text_analysis.replace('HUBBARDU', str(U))
# text_analysis = text_analysis.replace('DENS', str(d))
# text_analysis = text_analysis.replace('SIZE', str(Nc))
# text_analysis = text_analysis.replace('JOBS', batch_str_analysis)

# file = open(batch_name_analysis, 'w')
# file.write(text_analysis)
# file.close()

# Make batch scripts executable.
os.chmod(batch_name_dca, 0o755)
# os.chmod(batch_name_analysis, 0o755)
