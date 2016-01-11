import commands
import os
import sys
import time

import h5py

figext = '1'

if(figext == '1'):
    FIGURE_EXTENSION = ".png"

if(figext == '2'):
    FIGURE_EXTENSION = ".eps"

if(figext == '3'):
    FIGURE_EXTENSION = ".pdf"

def my_format(x):
    return ("%0.6f" % x).rjust(10)

Nc    = "16B"
temp  = "0.06"
hubbard_u = "7"
dens  ="1.000"
ch   = "PARTICLE_PARTICLE_SUPERCONDUCTING"

index  = "16"
sigma  = "0.5"
radius = "1.5"

for l in [1]:

    file_name = "./T_"+temp+"/data_d="+dens+"_Nc="+Nc+"_U="+hubbard_u+"_T="+temp+"_BSE_"+ch+"_tpcluster144.hdf5"
    print(file_name)

    data = h5py.File(file_name,'r')

    data_domains = data["domains"]

    for i in range(0, len(data.keys())):

        key_name = data.keys()[i]
        print key_name

        if(key_name == "functions"          or 
           key_name == "DCA-loop-functions" or 
           key_name == "CPE-functions"      or 
           key_name == "spectral-functions" or
           key_name == "analysis-functions" or
           key_name == "domains"            or
           key_name == "parameters" or
           key_name == "CT-AUX-SOLVER-functions" or
           key_name == "SS-HYB-SOLVER-functions" or
           key_name == "A-versus-G-functions"):
            
            data_functions = data[key_name]

            for j in range(0, len(data[key_name].keys())):
        
                function_name = data[key_name].keys()[j]

                data_function = data[key_name][function_name]

                print '\treading  '+ function_name + " ... "

                if(key_name == "analysis-functions"):
                    print data_function["data"]

                if(function_name == "leading-eigenvalues"):

                    print("\t-----------------------------------------")

                    print "\t " + function_name

                    print data_function["data"]

                    for pp in range(0, 10):
               
                        value1 = data_function["data"][pp, 0]
                        value2 = data_function["data"][pp, 1]
 
                        out = str(value1) + "\t" + str(value2) 
                        print(out)
                    print("\t-----------------------------------------")

                if(function_name == "leading-symmetry-decomposition"):
                    print("\t-----------------------------------------")
                    print "\t " + function_name

                    print data_function["data"]
                    print("\tdata[16,2,0,0,:,0]:")
        
                    for pp in range(0, 10):

                        value1 = data_function["data"][16, 2, 0, 0, pp, 0]
                        out = str(value1)
                        print(out)

                    print("\t-----------------------------------------")

                if(function_name == "standard_deviation"):
                    print("\t-----------------------------------------")
                    print "\t " + function_name

                    print data_function["data"]

                    print("\t Sigma_stddev for the last DCA+ iteration:")
                    for pp in range(0, 16):

                        value1 = data_function["data"][2, pp, 0, 0]
                        value2 = data_function["data"][2, pp, 1, 0]
                        out = "k" + str(pp) + "\t" + my_format(value1) + "\t" + my_format(value2)
                        print(out)
                    print("\t-----------------------------------------")

                if(function_name == "Greens-k-lattice"):
                    print("\t-----------------------------------------")
                    print "\t " + function_name

                    print data_function["data"]

                    for pp in range(0, 16):

                        value1 = data_function["data"][pp, 0, 0, 0, 0, 0]
                        value2 = data_function["data"][pp, 0, 0, 1, 0, 0]
                        value3 = data_function["data"][pp, 1, 0, 0, 0, 0]
                        value4 = data_function["data"][pp, 1, 0, 1, 0, 0]
                        value5 = data_function["data"][pp, 0, 0, 0, 0, 1]
                        value6 = data_function["data"][pp, 0, 0, 1, 0, 1]
                        value7 = data_function["data"][pp, 1, 0, 0, 0, 1]
                        value8 = data_function["data"][pp, 1, 0, 1, 0, 1]

                        out = str(value1) + "\t" + str(value2) + "\t" + str(value3) + "\t" + str(value4) + "\t" + str(value5) + "\t" + str(value6) + "\t" + str(value7) + "\t" + str(value8)
                       # print(out)
                    print("\t-----------------------------------------")

                if(function_name == "L2_Sigma_difference"):
                    print("\t-----------------------------------------")

                    print "\t " + function_name

                    print data_function["data"]

                    for pp in range(0, 3):

                        value1 = data_function["data"][pp]

                        print(str(value1))

                if(function_name == "sign"):
                    print("\t-----------------------------------------")

                    print "\t " + function_name

                    print data_function["data"]

                    for pp in range(0, 3):

                        value1 = data_function["data"][pp]

                        print(str(value1))

