import commands
import os
import sys
import time

import matplotlib
matplotlib.use('Agg')

import shutil
import os

import h5py

from pylab import *

import numpy as np
import scipy
from scipy.interpolate import Rbf
import matplotlib.pyplot as plt
from matplotlib import cm

import json
import pylab
import numpy

from pylab import *
from scipy import *

from scipy          import optimize
from scipy.optimize import curve_fit

figext = '3'

if(figext == '1'):
    FIGURE_EXTENSION = ".png"

if(figext == '2'):
    FIGURE_EXTENSION = ".eps"

if(figext == '3'):
    FIGURE_EXTENSION = ".pdf"

def my_format(x):
    return ('%.3g' % x)

def my_color(i, j):
    return cm.jet((i+0.)/(j+0.))

def configure_figure(x_label, y_label):

    ax = gca()

    ax.aspect = (sqrt(5)-1.0)/2.0
    
    ax.tick_params(direction='in', length=8, width=1.5, pad=10, labelsize=16)

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    xlabel(x_label, fontsize=30, labelpad=12)
    ylabel(y_label, fontsize=24, labelpad=12)
    
    subplots_adjust(left=0.15, bottom=0.125, right=0.975, top=0.955) 

    xlim(0, 1.1)
    
    legend(prop={'size':16})

def append_data(params, domains, functions):

    beta = params["physics-parameters"]["beta"][0]

    T_list.append(1./beta)

    G_k_list.append(functions["G_k_beta_over_2"        ]["data"][:])
    A_w_list.append(functions["integrated_A_k_div_cosh"]["data"][:])


U_str = "6"

pict_name = "G_k_versus_A_w_U="+U_str

file_name_tmpl = "./run_I/data_d=1.000_Nc=128_U="+U_str+"_CPE.hdf5" 

list = [0,1,2,3,4,5,6,7]

T_list = []

G_k_list = []
A_w_list = []

for l in list:

    n_ind = l
    N_tot = len(list)

    file_name = file_name_tmpl.replace("I", str(l))

    print "reading " + file_name + " ..."

    data = h5py.File(file_name,'r')

    data_domains = data["domains"]
    data_params  = data["parameters"]

    for i in range(0, len(data.keys())):

        key_name = data.keys()[i]
        #print key_name

        if(key_name == "functions"          or 
           key_name == "DCA-loop-functions" or 
           key_name == "CPE-functions"      or 
           key_name == "spectral-functions" or
           key_name == "analysis-functions" or
           key_name == "CT-AUX-SOLVER-functions" or
           key_name == "SS-HYB-SOLVER-functions" or
           key_name == "A-versus-G-functions"):
            
            data_functions = data[key_name]

            for j in range(0, len(data[key_name].keys())):
        
                function_name = data[key_name].keys()[j]

                data_function = data[key_name][function_name]

                if(function_name == "G_k_beta_over_2"):
                    append_data(data_params, data_domains, data_functions)


print T_list

G_k_list = transpose(G_k_list)
A_w_list = transpose(A_w_list)

figure(num=None)

for l in range(0,128):

    print T_list

    plot(T_list, G_k_list[l], "ko-")
    plot(T_list, A_w_list[l], "r+")

configure_figure(r"$\omega$", r"$\mathcal{A}(\omega)$")
               
savefig('pic_' + pict_name + FIGURE_EXTENSION)
