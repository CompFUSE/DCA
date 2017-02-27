import shutil
import os
import sys
import time
import csv

import h5py
import math

Nc    = ["16"]
temp  = [0.05, 0.045, 0.04, 0.035, 0.03, 0.025, 0.02]
hubbard_u = ["4"]
dens  =["0.900"]
ch   = "PARTICLE_PARTICLE_UP_DOWN"

index  = ["16"]
sigma  = ["0.5"]
radius = ["1.5"]

def my_format(x):
    return ("%0.6f" % x).rjust(10)

for k in range(0, len(hubbard_u)):
    for l in range(0, len(dens)):
        for j in range(0, len(Nc)):
            print("\t-----------------------------------------")
            print("\t Nc = " + Nc[j] + ", U = " + hubbard_u[k] + ", d = " + dens[l])
          #  print("\t-----------------------------------------")
            for t in range(0, len(temp)):

                text = ""

                for ia in range(0,len(index)):
                    for ib in range(0,len(sigma)):
                        for ic in range(0,len(radius)):

                            data_file_BSE = "./T="+str(temp[t])+"/data.BSE.hdf5"

               	            if(os.path.exists(data_file_BSE)):

                              #  print(data_file_BSE)

                                data = h5py.File(data_file_BSE,'r')
                                data_domains = data["domains"]

                                T_val = 1./data["parameters"]["physics-parameters"]["beta"][0]

                                for key_name in data.keys():

                                    if(key_name == "functions" or
                                       key_name == "DCA-loop-functions" or
                                       key_name == "CPE-functions" or
                                       key_name == "spectral-functions" or
                                       key_name == "analysis-functions" or
                                       key_name == "CT-AUX-SOLVER-functions"):

                                        data_functions = data[key_name]

                                        for function_name in data_functions.keys():

                                            #function_name = data[key_name].keys()[l_1]                                                                                                                                                                          
                                            data_function = data[key_name][function_name]

                                            #print '\tread '+ function_name + " ... ? "                                                                                                                                                                                                                                                                                                                                                  
                                            if(function_name=="leading-eigenvalues"):

                                                if(False):

                                                    out_str = "" +my_format(i) + "\t" + my_format(T_val)

                                                    for l0_tmp in range(0, 2):
                                                        out_str = out_str + "\t" + my_format(data_function["data"][l0_tmp,0]) #+ " ( "                                                                                                                                                                                    #print data[key_name]["leading-symmetry-decomposition"]["data"].shape
                                                        for l1_tmp in range(0, 3):
                                                            out_str = out_str + my_format(data[key_name]["leading-symmetry-decomposition"]["data"][16, l1_tmp, 0, 0, l0_tmp, 0]) + " , "

                                                    text = text + out_str + "\n"

                                                else:
                                                    tmp=[-1]
                                        
                                                    for l0_tmp in range(0, 10):
                                                        d_wave_proj = data[key_name]["leading-symmetry-decomposition"]["data"][16, 2, 0, 0, l0_tmp, 0]

                                                        if(abs(abs(d_wave_proj)-1)<0.1):
                                                   #     if(abs(d_wave_proj-1)<0.1):
                                                            tmp.append(data_function["data"][l0_tmp,0])

                                                    lam_max = max(tmp)

                                                    if(lam_max>0):
                                                        out_str = "\t" + my_format(T_val) + "\t" + my_format(lam_max)

                                                        print(out_str)

                                                        text = text + out_str + "\n"

#            if(True and not text == ""):

#                file = open("./data/data_d="+str(dens[l])+"_Nc="+str(Nc[j])+"_U="+str(hubbard_u[k])+"_r="+radius+"_s="+sigma+"_i="+index+".txt", "w")
#                file.write(text)
#                file.close()

