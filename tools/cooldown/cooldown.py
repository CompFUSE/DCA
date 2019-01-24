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

import os
import sys

################################################################################
##################################### EDIT #####################################
################################################################################

# Filename (including path) of the batch script template.
# The file should contain the placeholder JOBS.
# Other optional placeholders: APPLICATION, HUBBARDU, DENS, SIZE.
batch_tmpl = 'EDIT'

# Command (including options) to run (MPI) programs, e.g. 'mpirun -n 8'.
# Used to execute main_dca.
run_command_dca = 'EDIT'
# The run command to execute main_analysis.
run_command_analysis = run_command_dca

# Simulation parameters
dca_plus = False  # True|False
U = 6             # on-site Coulomb interaction
d = 0.95          # density

# DCA cluster (see cluster_definitions.txt)
Nc = 4 # Cluster size
cluster_vec_1 = [2, 0]
cluster_vec_2 = [0, 2]

# DCA(+) iterations
iters_ht = 8  # first/highest temperature
iters = 6     # other temperatures

# Temperatures for cooldown
temps = [1, 0.75, 0.5, 0.25, 0.125, 0.1, 0.09, 0.08, 0.07]
# Starting temperature for computing two-particle quantities and running the analysis.
T_analysis = 0.1

################################################################################
################################## UNTIL HERE ##################################
################################################################################

# Formats the inverse temperature beta.
def beta_format(x):
    return ('%.6g' % x)

# Prepares the input file for the given temperature.
def prepare_input_file(filename, T_ind):
    file = open(filename, 'r')
    text = file.read()
    file.close()

    text = text.replace('CURRENT_TEMP', str(temps[T_ind]))
    text = text.replace('BETA', str(beta_format(1./temps[T_ind])))
    text = text.replace('DENS', str(d))
    text = text.replace('HUBBARDU', str(U))

    # For the first temperature set the initial self-energy to zero.
    if (T_ind == 0):
        text = text.replace('./T=PREVIOUS_TEMP/dca_sp.hdf5', 'zero')
    else:
        text = text.replace('PREVIOUS_TEMP', str(temps[T_ind-1]))

    if (T_ind == 0):
        text = text.replace('ITERS', str(iters_ht))
    else:
        text = text.replace('ITERS', str(iters))

    # Use lower() since Python's booleans start with a capital letter, while JSON's booleans don't.
    text = text.replace('DO_DCA_PLUS', str(dca_plus).lower())

    text = text.replace('VEC1', str(cluster_vec_1))
    text = text.replace('VEC2', str(cluster_vec_2))

    file = open(filename, 'w')
    file.write(text)
    file.close()

################################################################################

batch_str_dca = ''
batch_str_analysis = ''

print('Generating directories and input files:')

for T_ind, T in enumerate(temps):
    print('T = ' + str(T))

    # Create directory.
    dir_str = './T=' + str(T)
    cmd = 'mkdir -p ' + dir_str
    os.system(cmd)
    cmd = 'mkdir -p ' + dir_str + "/configuration"
    os.system(cmd)

    input_sp = dir_str + '/input_sp.json'
    input_tp = dir_str + '/input_tp.json'

    data_dca_sp   = dir_str + '/dca_sp.hdf5'
    data_dca_tp   = dir_str + '/dca_tp.hdf5'
    data_analysis = dir_str + '/analysis.hdf5'

    # dca sp
    # Generate the sp input file.
    cmd = 'cp ./input_sp.json.in ' + input_sp
    os.system(cmd)
    prepare_input_file(input_sp, T_ind)

    # Add job.
    batch_str_dca = batch_str_dca + run_command_dca + ' ./main_dca ' + input_sp + '\n'

    if (T <= T_analysis):
        # dca tp
        # Generate the tp input file.
        cmd = 'cp ./input_tp.json.in ' + input_tp
        os.system(cmd)
        prepare_input_file(input_tp, T_ind)

        # Add job.
        batch_str_dca = batch_str_dca + run_command_dca + ' ./main_dca ' + input_tp + '\n'

        # analysis
        # Add job.
        batch_str_analysis = batch_str_analysis + run_command_analysis + ' ./main_analysis ' + input_tp + '\n'


# Get filename extension of batch script.
_, extension = os.path.splitext(batch_tmpl)

# Generate the dca batch script.
batch_name_dca = 'job.dca_U=' + str(U) + '_d=' + str(d) + '_Nc=' + str(Nc) + extension
print('\nGenerating the dca batch script: ' + batch_name_dca)
file = open(batch_tmpl, 'r')
text_dca = file.read()
file.close()

text_dca = text_dca.replace('APPLICATION', 'dca')
text_dca = text_dca.replace('HUBBARDU', str(U))
text_dca = text_dca.replace('DENS', str(d))
text_dca = text_dca.replace('SIZE', str(Nc))
text_dca = text_dca.replace('JOBS', batch_str_dca)

file = open(batch_name_dca, 'w')
file.write(text_dca)
file.close()

# Generate the analysis batch script.
batch_name_analysis = 'job.analysis_U=' + str(U) + '_d=' + str(d) + '_Nc=' + str(Nc) + extension
print('Generating the analysis batch script: ' + batch_name_analysis)
file = open(batch_tmpl, 'r')
text_analysis = file.read()
file.close()

text_analysis = text_analysis.replace('APPLICATION', 'analysis')
text_analysis = text_analysis.replace('HUBBARDU', str(U))
text_analysis = text_analysis.replace('DENS', str(d))
text_analysis = text_analysis.replace('SIZE', str(Nc))
text_analysis = text_analysis.replace('JOBS', batch_str_analysis)

file = open(batch_name_analysis, 'w')
file.write(text_analysis)
file.close()

# Make batch scripts executable.
os.chmod(batch_name_dca, 0o755)
os.chmod(batch_name_analysis, 0o755)
