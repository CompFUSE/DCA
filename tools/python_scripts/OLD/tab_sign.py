import commands
import shutil
import os
import sys
import time
import csv

import h5py
import math

Nc    = ["24"]
#Nc    = ["08", "10", "12", "14", "16A", "16B", "20", "24"]
temp  = [0.2, 0.167, 0.125, 0.1 ,0.08, 0.06]# 0.08, 0.06, 0.055, 0.05, 0.045, 0.04, 0.035, 0.03, 0.025, 0.02]
hubbard_u = ["7"]
#dens      = ["1.000", "0.975", "0.95", "0.925", "0.900"]#, "0.875", "0.850", "0.825", "0.800"]
dens      = ["0.900"]

def my_format(x):
    return ('%.6g' % x)

for k in range(0, len(hubbard_u)):
    for p in range(0, len(dens)):
        for j in range(0, len(Nc)):

            print("\t Nc = " + Nc[j] + ", U = " + hubbard_u[k] + ", d = " + dens[p] )
            for t in range(0, len(temp)):

                for l in range(0, 1):        

                    filename = "./T_"+str(temp[t])+"/data_d="+str(dens[p])+"_Nc="+str(Nc[j])+"_U="+str(hubbard_u[k])+"_T="+str(temp[t])+"_DCA+.hdf5"
        
                    if(os.path.exists(filename)):
    
                        data = h5py.File(filename,'r')

                        data_domains = data["domains"]

                        T_val = 1./data["parameters"]["physics-parameters"]["beta"][0]

                        for i in range(0, len(data.keys())):

                            key_name = data.keys()[i]

                            if(key_name == "functions" or
                               key_name=="DCA-loop-functions" or
                               key_name=="CPE-functions" or
                               key_name=="spectral-functions" or
                               key_name=="CT-AUX-SOLVER-functions" or
                               key_name=="analysis-functions"):
                    
                                data_functions = data[key_name]

                                for q in range(0, len(data[key_name].keys())):
                    
                                    function_name = data[key_name].keys()[q]
                    
                                    data_function = data[key_name][function_name]
                    
                                    #print '\tplot '+ function_name + " ... ? "
                
                                    if(function_name=="sign"):
                                        print "\t" + my_format(T_val) + "\t" + my_format(data_function["data"][0])

