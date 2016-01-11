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

N_iter = 1
Nc = 16

Ncc    = "16B"
temp  = "0.06"
hubbard_u = "7"
dens  ="1.000"

for l in [1]:

    file_name = "./T_"+temp+"/data_d="+dens+"_Nc="+Ncc+"_U="+hubbard_u+"_T="+temp+"_DCA+_PARTICLE_PARTICLE_SUPERCONDUCTING.hdf5"

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

                if(key_name == "functions"):
                   print data_function["data"]

                if(function_name == "Self_Energy"):
                    print("\t-----------------------------------------")

                    print "\t " + function_name

                    print data_function["data"]

                    for pp in range(0, N_iter):

                        value1 = data_function["data"][pp]

#                        print(str(value1))


                if(function_name == "Self_Energy-stddev"):

                    print("\t-----------------------------------------")

                    print "\t " + function_name

                    print data_function["data"]

                    for pp in range(0, Nc):
               
                        value1 = data_function["data"][20, pp, 0, 0, 0, 0, 0]
                        value2 = data_function["data"][20, pp, 1, 0, 0, 0, 0]
                        value3 = data_function["data"][20, pp, 0, 0, 1, 0, 0]
                        value4 = data_function["data"][20, pp, 1, 0, 1, 0, 0]
                        value5 = data_function["data"][20, pp, 0, 0, 0, 0, 1]
                        value6 = data_function["data"][20, pp, 1, 0, 0, 0, 1]
                        value7 = data_function["data"][20, pp, 0, 0, 1, 0, 1]
                        value8 = data_function["data"][20, pp, 1, 0, 1, 0, 1]  
 
                        out = str(value1) + "\t" + str(value2) + "\t" + str(value3) + "\t" + str(value4) + "\t" + str(value5) + "\t" + str(value6) + "\t" + str(value7) + "\t" + str(value8)
                        print(out)
                    print("\t-----------------------------------------")

                if(function_name == "Sigma_zero_moment"):
                    print("\t-----------------------------------------")
                    print "\t " + function_name

                    print data_function["data"]

                    out = "Iter       0"
                    print(out)

                    for pp in range(0, Nc):
                        out = "k" + str(pp)
                        for iter in range(0, N_iter):
                            value = data_function["data"][iter, pp, 0, 0]
                            out = out + "\t" + my_format(value)
                        print(out)

                    print("\t-----------------------------------------")


                if(function_name == "standard_deviation"):
                    print("\t-----------------------------------------")
                    print "\t " + function_name

                    print data_function["data"]

                    print("\t Sigma_stddev for the last DCA+ iteration:")
                    for pp in range(0, Nc):

                        value1 = data_function["data"][N_iter-1, pp, 0, 0]
                        value2 = data_function["data"][N_iter-1, pp, 1, 0]
                        out = "k" + str(pp) + "\t" + my_format(value1) #+ "\t" + my_format(value2)
                        print(out)
                    print("\t-----------------------------------------")

                if(function_name == "Greens-k-lattice"):
                    print("\t-----------------------------------------")
                    print "\t " + function_name

                    print data_function["data"]

                    for pp in range(0, Nc):

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

                    for pp in range(0, N_iter):

                        value1 = data_function["data"][pp]

                        print(str(value1))

                if(function_name == "sign"):
                    print("\t-----------------------------------------")

                    print "\t " + function_name

                    print data_function["data"]

                    for pp in range(0, N_iter):

                        value1 = data_function["data"][pp]

                        print(str(value1))

