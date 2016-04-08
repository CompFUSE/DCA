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

N_iter = 0
Nc = 16

Ncc    = "16"
temp  = "0.125"
hubbard_u = "7"
dens  ="0.942"

for l in [1]:

    file_name = "./T=0.125/data.DCA_sp.hdf5"
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

                if(key_name == "functions"):
                   print data_function["data"]

                if(function_name == "Self_Energy"):
                    print("\t-----------------------------------------")

                    print "\t " + function_name

                    print data_function["data"]

                    k_dmn    = data["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"]["data"]
                    shape = k_dmn.shape
                    Kx = k_dmn[:,0]
                    Ky = k_dmn[:,1]

               #     for iK in range(0,shape[0]):
               #         print iK, "  kx = ",Kx[iK],"ky = ",Ky[iK]

                   # for pp in range(240, 272):

                   #     value1 = data_function["data"][pp,2,0,0,0,0]
                   #     value2 = data_function["data"][pp,8,0,0,0,0]
                   #     a = str(data["domains"]["frequency-domain"]["elements"][pp])+"  "+str(value1)+"  "+str(value2)
                   #     print(a)

                if(function_name == "Self_Energy-stddev"):

                    print("\t-----------------------------------------")

                    print "\t " + function_name

                    print data_function["data"]

                    # sigma(K,pi*T)
                    for pp in range(0, Nc):

                        value = data_function["data"][256, pp, 0, 0, 0, 0, 1]

                        print Kx[pp], "  ", Ky[pp], "  ", data["functions"]["Self_Energy"]["data"][256, pp, 0, 0, 0, 0, 1], "  ", value
                    print("\t-----------------------------------------")

                    # sigma(K=(pi,0),iwn)
                    for pp in range(240, 272):

                        value = data_function["data"][pp, 8, 0, 0, 0, 0, 1]

                        print data["domains"]["frequency-domain"]["elements"][pp], "  ", data["functions"]["Self_Energy"]["data"][pp, 8, 0, 0, 0, 0, 1], "  ", value
                    print("\t-----------------------------------------")

                if(function_name == "Sigma_zero_moment"):
                    print("\t-----------------------------------------")
                    print "\t " + function_name

                    print data_function["data"]

                    out = "Iter       0          1          2         3"
                   # print(out)

                    for pp in range(0, Nc):
                        out = "k" + str(pp)
                        for iter in range(0, N_iter):
                            value = data_function["data"][iter, pp, 0, 0]
                            out = out + "\t" + my_format(value)
                    #    print(out)

                    print("\n Further analyze the Sigma difference between iterations:")

                    out = "Diff_Iter       01          12          23"
                    print(out)

                    for pp in range(0, Nc):
                        out = "k" + str(pp)
                        for iter in range(0, N_iter-1):
                            diff = abs( data_function["data"][iter+1, pp, 0, 0] - data_function["data"][iter, pp, 0, 0])
                            out = out + "\t" + my_format(diff)
                    #    print(out)

                    print("\t-----------------------------------------")


                if(function_name == "standard_deviation"):
                    print("\t-----------------------------------------")
                    print "\t " + function_name

                    print data_function["data"]

                    print("\t Sigma_stddev for the last DCA+ iteration:")

                    out = "Iter       0          1          2         3"
                    print(out)

                    for pp in range(0, Nc):
                        out = "k" + str(pp)
                        for iter in range(0, N_iter):
                            value = data_function["data"][iter, pp, 0, 0]
                            out = out + "\t" + my_format(value)
                   #     print(out)

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

                if(function_name == "Monte-Carlo-Integration"):
                    print("\t-----------------------------------------")

                    print data.keys()
                    print key_name

