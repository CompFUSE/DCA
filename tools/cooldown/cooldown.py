# Creates directories, input files and batch scripts for a cooldown, i.e. a series of dca/analysis
# runs with gradually decreased temperature.
# By default this script and the provided template input files (input_sp.json.in and
# input_sp.json.in) are configured for the single-band Hubbard model with on-site Coulomb
# repulsion U and simulations in the grand-canonical ensemble with fixed density d.
#
# See https://github.com/eth-cscs/dca_ethz/wiki/Running for details on how to use this script.
#
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)

import os

################################################################################
##################################### EDIT #####################################
################################################################################

batch_tmpl = "EDIT"                          # Batch script template
_, extension = os.path.splitext(batch_tmpl)  # Get filename extension of batch script.

application = "dca"  # "dca" | "analysis"

# Resources
# Used to substitute placeholders in the batch script template and the run command.
nodes = 1              # Replaces NODES.
walltime = "01:00:00"  # Format "hh:mm:ss". Replaces WALLTIME.
# These options can be used by the 'run_command'.
tasks_per_node = 1
tasks = 1
threads = 1

# The command (including options) to run (MPI) programs.
# The batch script will use it as
# run_command ./main_dca ./T=0.1/input_sp.json
run_command = "EDIT"

# Parameters
dca_plus = False
U = 8     # on-site Coulomb repulsion
d = 0.95  # density

# DCA(+) iterations
iters_ht = 8  # first/highest temperature
iters = 3     # other temperatures

# DCA cluster (see cluster_definitions.txt)
Nc = 4
cluster_vec_1 = [2, 0]
cluster_vec_2 = [0, 2]

# Temperatures for cooldown
temps = [1, 0.75, 0.5, 0.25, 0.2, 0.15, 0.125, 0.1, 0.09, 0.085, 0.08, 0.075, 0.07]
# Starting temperature for computing two-particle quantities and running the analysis.
T_analysis = 0.1 

################################################################################
################################## UNTIL HERE ##################################
################################################################################

# Formats the inverse temperature beta.
def beta_format(x):
    return ("%.6g" % x)

# Prepares the input file for the given temperature.
def prepare_input_file(filename, T_ind):
    file = open(filename, "r")
    text = file.read()
    file.close()

    text = text.replace("CURRENT_TEMP", str(temps[T_ind]))
    text = text.replace("BETA", str(beta_format(1./temps[T_ind])))
    text = text.replace("DENS", str(d))
    text = text.replace("HUBBARDU", str(U))

    # For the first temperature set the initial self-energy to zero.
    if (T_ind == 0):
        text = text.replace("./T=PREVIOUS_TEMP/dca_sp.hdf5", "zero")
    else:
        text = text.replace("PREVIOUS_TEMP", str(temps[T_ind-1]))
    
    if (T_ind == 0):
        text = text.replace("ITERS", str(iters_ht))
    else:
        text = text.replace("ITERS", str(iters))

    # Use lower() since Python's booleans start with a capital letter, while JSON's booleans don't.
    text = text.replace("DO_DCA_PLUS", str(dca_plus).lower())
            
    text = text.replace("VEC1", str(cluster_vec_1))
    text = text.replace("VEC2", str(cluster_vec_2))

    file = open(filename, "w")
    file.write(text)
    file.close()


batch_str = ""

for T_ind, T in enumerate(temps):
    print "T =", T

    # Create directory if it doesn't exist yet.
    if (not os.path.exists("./T=" + str(T))):
        cmd = "mkdir T=" + str(T)
        os.system(cmd)

    dir_str = "./T=" + str(T)

    input_sp = dir_str + "/input_sp.json"
    input_tp = dir_str + "/input_tp.json"

    data_dca_sp   = dir_str + "/dca_sp.hdf5"
    data_dca_tp   = dir_str + "/dca_tp.hdf5"
    data_analysis = dir_str + "/analysis.hdf5"

    # dca sp
    if (application == "dca" and not os.path.exists(data_dca_sp)):
        # Generate the input file.
        cmd = "cp ./input_sp.json.in " + input_sp
        os.system(cmd)
        prepare_input_file(input_sp, T_ind)

        # Add job.
        batch_str = batch_str + run_command + " ./main_dca " + input_sp + "\n"

    # dca tp
    if (application == "dca" and T<=T_analysis and not os.path.exists(data_dca_tp)):
        # Generate the input file.
        cmd = "cp ./input_tp.json.in " + input_tp
        os.system(cmd)
        prepare_input_file(input_tp, T_ind)

        # Add job.
        batch_str = batch_str + run_command + " ./main_dca " + input_tp + "\n"

    # analysis
    if (application == "analysis" and os.path.exists(data_dca_tp) and not os.path.exists(data_analysis)):
        # Add job.
        batch_str = batch_str + run_command + " ./main_analysis " + input_tp + "\n"

# Generate the batch script.
file = open(batch_tmpl, "r")
text = file.read()
file.close()

text = text.replace("APPLICATION", application)
text = text.replace("NODES", str(nodes))
text = text.replace("WALLTIME", walltime)
text = text.replace("HUBBARDU", str(U))
text = text.replace("DENS", str(d))
text = text.replace("SIZE", str(Nc))
text = text.replace("JOBS", batch_str)

batch_name = "job." + application + "_U=" + str(U) + "_d=" + str(d) + "_Nc=" + str(Nc) + extension
file = open(batch_name, "w")
file.write(text)
file.close()
