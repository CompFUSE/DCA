import os
import math

section = "dca"  # "dca"|"analysis"
algorithm = "DCA"  # "DCA+"|"DCA"

# Parameters
U = 7
V = 0  # assume V = V-prime now
d = 0.942

period = 2

temps = [1, 0.8, 0.7, 0.6, 0.5, 0.33, 0.25, 0.2, 0.167, 0.15, 0.125]
T_start_analysis = 0.05  # Starting temperature for computing two-particle quantities and doing the analysis.

Nc = 16  # See cluster_definitions
vecx = [2, 4]
vecy = [4, 0]

channel = "PARTICLE_PARTICLE_UP_DOWN"
qvec = [0, 0]

sigma_cutoff = 0.5
radius_cutoff = 1.5


def my_format(x):
    return ("%.6g" % x)

def prepare_input_file(filename, T_ind):
    file = open(filename, "r")
    text = file.read()
    file.close()

    text = text.replace("ALGORITHM", algorithm)
    
    if (T_ind == 0):
        text = text.replace("./T=PREVIOUS_TEMP/data."+algorithm+"_sp.hdf5", "zero")
    else:
        text = text.replace("PREVIOUS_TEMP", str(temps[T_ind-1]))

    text = text.replace("TEMP", str(temps[T_ind]))
        
    text = text.replace("BETA", str(my_format(1./temps[T_ind])))
    text = text.replace("DENS", str(d))
    text = text.replace("UVAL", str(U))
    text = text.replace("VVAL", str(V))
        
    if (temps[T_ind] >= 1.0):
        text = text.replace("ITERS", "8")
    else:
        text = text.replace("ITERS", "4")

    text = text.replace("VECX", str(vecx))
    text = text.replace("VECY", str(vecy))

    text = text.replace("PERIOD", str(period))

    text = text.replace("WARM", "20")
    text = text.replace("SWEEP", "1")
    if (temps[T_ind] >= 0.06):
        text = text.replace("MEAS", "100")
    else:
        text = text.replace("MEAS", "100")

    text = text.replace("SEED", "985456376")

    text = text.replace("CHANNEL", channel)
    text = text.replace("QVEC", str(qvec))
    text = text.replace("SIGMACUTOFF", str(sigma_cutoff))
    text = text.replace("RADIUSCUTOFF", str(radius_cutoff))

    file = open(filename, "w")
    file.write(text)
    file.close()


batch_str = ""

for T_ind, T in enumerate(temps):
    print my_format(T)

    if (not os.path.exists("./T=" + str(T))):
        cmd = "mkdir T=" + str(temps[T_ind])
        os.system(cmd)

    dir_str = "./T=" + str(T)

    input_sp = dir_str + "/input.sp.json"
    input_tp = dir_str + "/input.tp.json"

    data_dca_sp   = dir_str + "/data."+algorithm+"_sp.hdf5"
    data_dca_tp   = dir_str + "/data."+algorithm+"_tp.hdf5"
    data_analysis = dir_str + "/data.BSE.hdf5"

    # dca sp
    if (section == "dca" and not os.path.exists(data_dca_sp)):
        cmd = "cp ./input.sp.json.in " + input_sp
        os.system(cmd)
        prepare_input_file(input_sp, T_ind)

        batch_str = batch_str + "aprun -B ./main_dca " + input_sp + "\n"

    # dca tp
    if (section == "dca" and T<=T_start_analysis and not os.path.exists(data_dca_tp)):
        cmd = "cp ./input.tp.json.in " + input_tp
        os.system(cmd)
        prepare_input_file(input_tp, T_ind)

        batch_str = batch_str + "aprun -B ./main_dca " + input_tp + "\n"

    # analysis
    if (section == "analysis" and os.path.exists(data_dca_tp) and not os.path.exists(data_analysis)):
        batch_str = batch_str + "aprun -B ./main_analysis " + input_tp + "\n"

if (section == "dca"):
    file = open("job.dca.slm.in", "r")
    text = file.read()
    file.close()
    
    batch_script_name = "job.U="+str(U)+"_d="+str(d)+"_Nc="+str(Nc)+".dca_"+algorithm+".slm"
        
elif (section == "analysis"):
    file = open("job.analysis.slm.in", "r")
    text = file.read()
    file.close()
    
    batch_script_name = "job.U="+str(U)+"_d="+str(d)+"_Nc="+str(Nc)+".analysis_"+algorithm+".slm"
    
text = text.replace("DENS", str(d))
text = text.replace("HUBBARDU", str(U))
text = text.replace("VVAL"    , str(V))
text = text.replace("SIZE", str(Nc))
text = text.replace("ALGORITHM", algorithm)
text = text.replace("JOBS", batch_str)

file = open(batch_script_name, "w")
file.write(text)
file.close()
