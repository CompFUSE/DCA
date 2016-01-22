
import commands
import os
import sys
import time

import matplotlib
matplotlib.use('Agg')

import shutil
import os

from pylab import *

import json
import csv

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


filename = raw_input('Please enter input-filename : ... \n\t 1 = data.json \n\t 2 = data_analysis.json \n\t 3 = data_susceptibilities.json \n\t 4 = data_spectrum.json \n\n   ')

if(filename == '1'):
    filename = 'data.json'

if(filename == '2'):
    filename = 'data_analysis.json'

if(filename == '3'):
    filename = 'data_susceptibilities.json'

if(filename == '4'):
    filename = 'data_spectrum.json'

figext = raw_input('Please enter figure-extension : ... \n\t 1 = png \n\t 2 = eps \n\t 3 = pdf \n\n   ')

if(figext == '1'):
    FIGURE_EXTENSION = ".png"

if(figext == '2'):
    FIGURE_EXTENSION = ".eps"

if(figext == '3'):
    FIGURE_EXTENSION = ".pdf"

left   = 0.125 # the left side of the subplots of the figure
right  = 0.9   # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top    = 0.9   # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.8   # the amount of height reserved for white space between subplots

F = figure(1)
#DefaultSize = F.get_size_inches()
#F.set_size_inches( (DefaultSize[0]*4, DefaultSize[1]*3) )


def my_colormap(x):
    return cm.jet(x)

def my_format(x):
    return ('%.3g' % x)

def read_function_name(line):

    function_name_stri = '{' + line[0:len(line)-2] + '}'
    function_name_json = json.loads(function_name_stri)
    function_name = function_name_json['name']
    return function_name

def read_function_matrix(file):
    
    stringmatrix = '{'

    for line in file:
        if "data" in line:
            stringmatrix = stringmatrix+line

            for line in file:
                if "}" in line:
                    break
                else:
                    stringmatrix = stringmatrix+line
            break        

    stringmatrix = stringmatrix+'}'
    
    return stringmatrix

def plot_contourf_hot(min, max, x, y, z, function_name, pic_name, pic_title):

    import numpy as np
    from scipy.interpolate import griddata
    import matplotlib.pyplot as plt
    import numpy.ma as ma
    from numpy.random import uniform

    # define grid.
    xi = np.linspace(min,max,200)
    yi = np.linspace(min,max,200)

    # grid the data.
    zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='linear')
    #zi = griddata((x,y),z,xi,yi)

    figure(num=None)    
    grid(True)

# contour the gridded data, plotting dots at the randomly spaced data points.
    CS = plt.contour(xi,yi,zi,5,linewidths=0.5,colors='r')
    CS = plt.contourf(xi,yi,zi,100,cmap=plt.cm.hot)

    plt.colorbar() # draw colorbar

# plot data points.
    plt.scatter(x,y,marker='o',c='k',s=5)
    plt.xlim(min,max)
    plt.ylim(min,max)
    plt.title(pic_title)

    xticks([min, min/2, 0, max/2, max], [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])
    yticks([min, min/2, 0, max/2, max], [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])

    savefig('./'+filename+'/'+function_name+'/'+pic_name)

def plot_contourf_jet(min, max, x, y, z, function_name, pic_name, pic_title):
    import numpy as np
    from scipy.interpolate import griddata
    import matplotlib.pyplot as plt
    import numpy.ma as ma
    from numpy.random import uniform

    # define grid.
    xi = np.linspace(min,max,200)
    yi = np.linspace(min,max,200)

    # grid the data.
    zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
    #zi = griddata((x,y),z,xi,yi)

    figure(num=None)    

    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

    CS = plt.contour(xi,yi,zi,10,linewidths=1,colors='k')
    CS = plt.contourf(xi,yi,zi,100,cmap=plt.cm.jet)
    plt.colorbar() # draw colorbar

    plt.xlim(min,max)
    plt.ylim(min,max)

    xticks([min, min/2, 0, max/2, max], [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])
    yticks([min, min/2, 0, max/2, max], [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])

    plt.title(pic_title)

    savefig('./'+filename+'/'+function_name+'/'+pic_name)


def plot_contourf_jet_nearest(min, max, x, y, z, function_name, pic_name, pic_title):
    import numpy as np
    from scipy.interpolate import griddata
    import matplotlib.pyplot as plt

    # define grid.
    xi = np.linspace(min,max,64)
    yi = np.linspace(min,max,64)


    #zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
    #zi = griddata((x,y),z,xi,yi)

    zi = []

    for xind in range(0, len(xi)):

        ztmp = []
        
        for yind in range(0, len(yi)):

            d   = []
            tmp = []

            for l in range(0, len(x)):
                d.append( (xi[xind]-x[l])**2 + (yi[yind]-y[l])**2 )
                tmp.append(z[l])

            f = transpose([d, tmp])
            f = f[f[:,0].argsort(),]
            f = transpose(f)

            ztmp.append(f[1][0])
            
        zi.append(ztmp)
            
    figure(num=None)    

    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

    #CS = plt.contour(xi,yi,zi,10,linewidths=1,colors='k')
    CS = plt.contourf(xi,yi,zi,100,cmap=plt.cm.jet)
    plt.colorbar() # draw colorbar

    plt.xlim(min,max)
    plt.ylim(min,max)

    xticks([min, min/2, 0, max/2, max], [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])
    yticks([min, min/2, 0, max/2, max], [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])

    plt.title(pic_title)

    savefig('./'+filename+'/'+function_name+'/'+pic_name)


def plot_pcolor(min, max, x, y, z, function_name, pic_name, pic_title):

    import numpy as np
    import scipy 
    from scipy.interpolate import Rbf
    import matplotlib.pyplot as plt
    from matplotlib import cm

    ti = np.linspace(min, max, 100)
    XI, YI = np.meshgrid(ti, ti)
    rbf = Rbf(x, y, z, epsilon=2)
    ZI = rbf(XI, YI)

    figure(num=None)    
    grid(True)
    
    n = plt.normalize(min, max)
    plt.subplot(1, 1, 1)
    plt.pcolor(XI, YI, ZI, cmap=cm.jet)
    plt.colorbar()
    plt.scatter(x, y, 50, 'k')
    plt.title(pic_title)
    plt.xlim(min, max)
    plt.ylim(min, max)

    savefig('./'+filename+'/'+function_name+'/'+pic_name)

def plot_spectrum(domains, function_name, function_matrix):

    r_freq_dmn = domains["domains"]["frequency-domain-real-axis"]["elements"]

    if os.path.exists('./'+filename+'/' + function_name):
        shutil.rmtree('./'+filename+'/' + function_name)
    os.makedirs('./'+filename+'/' + function_name)

#     bands = max(function_matrix[:,0]) + 1;
    times = max(function_matrix[:,0]) + 1;

    f_real = []
#     f_imag = []

#     for i in range(0,bands):
#         for j in range(0,bands):
    for l in range(0,times): 
        f_real.append(function_matrix[int(l),1])
#         f_imag.append(function_matrix[int(l),2])
                     
    figure(num=None)#, figsize=(10, 10), dpi=200)   
                    
    grid(True)
    plot(r_freq_dmn, f_real, 'k')

    #xlim(-5, 5)
    
    savefig('./'+filename+'/' + function_name + '/'+function_name+FIGURE_EXTENSION)
    
    f_real = [] 
    f_imag = [] 


def plot_real_axis(domains, function_name, function_matrix):

    r_freq_dmn = domains["domains"]["frequency-domain-real-axis"]["elements"]#domains["Pade_real_frequency_domain"]

    if os.path.exists('./'+filename+'/' + function_name):
        shutil.rmtree('./'+filename+'/' + function_name)
    os.makedirs('./'+filename+'/' + function_name)

    bands = max(function_matrix[:,0]) + 1;
    sites = max(function_matrix[:,4]) + 1;
    freqs = max(function_matrix[:,5]) + 1;

    f_real = []
    f_imag = []

    for i in range(0,bands):
        for j in range(0,bands):
            for r in range(0,sites):
                for l in range(0,freqs): 
                    f_real.append(function_matrix[int((i+bands)+(j+bands)*2*bands+r*4*bands*bands+l*4*bands*bands*sites),6])
                    f_imag.append(function_matrix[int((i+bands)+(j+bands)*2*bands+r*4*bands*bands+l*4*bands*bands*sites),7])
                    
                figure(num=None)#, figsize=(10, 10), dpi=200)   
                    
                grid(True)
                plot(r_freq_dmn, f_real, 'k')
                plot(r_freq_dmn, f_imag, 'r')
                legend(('real', 'imag'), 'upper right', shadow=True, fancybox=True)
                
                savefig('./'+filename+'/' + function_name + '/'+function_name+'_{' + str(r) + '_'+ str(i) + '_' + str(j) + '}'+FIGURE_EXTENSION)
                
                f_real = [] 
                f_imag = [] 

    answer = 'n'#raw_input('plot f(pi/T)? (y/n)')

    if(answer == 'y'):

        import numpy as np
        import scipy 
        from scipy.interpolate import Rbf
        import matplotlib.pyplot as plt
        from matplotlib import cm

        for i in range(0,bands):
            for j in range(0,bands):
                
                figure(num=None)
                 
                x = []
                y = []
                
                vec_re = []
                vec_im = []
                
                for k_0 in range(0,sites):
                    for l0 in [-1, 0, 1]:
                        for l1 in [-1, 0, 1]:
                            vec_re.append(function_matrix[int((i+bands)+(j+bands)*2*bands+r*4*bands*bands+(freqs/2.+1)*4*bands*bands*sites),6])
                            vec_im.append(function_matrix[int((i+bands)+(j+bands)*2*bands+r*4*bands*bands+(freqs/2.+1)*4*bands*bands*sites),7])
                            
                            K0 = domains["cluster"]["BASIS"]["K"]["K_1"][0];
                            K1 = domains["cluster"]["BASIS"]["K"]["K_2"][0];
                            
                            x.append(domains["cluster"]["FULL-CLUSTER"]["K"][k_0][0] + l0*K0[0]+l1*K1[0])
                            y.append(domains["cluster"]["FULL-CLUSTER"]["K"][k_0][1] + l0*K0[1]+l1*K1[1])

                plot_pcolor(-3.141592, 3.141592, x, y, vec_re, function_name, 'w=0_Re' + str(i) + '_' + str(i) + FIGURE_EXTENSION, 'Re[Sigma]')
                plot_pcolor(-3.141592, 3.141592, x, y, vec_im, function_name, 'w=0_Im' + str(i) + '_' + str(i) + FIGURE_EXTENSION, 'Im[Sigma]')



def plot_fine_real_axis(domains, function_name, function_matrix):

#     r_freq_dmn = domains["Pade_real_frequency_domain"]
#     i_freq_dmn = domains["Pade_imag_frequency_domain"]

    r_freq_dmn = domains["domains"]["frequency-domain-real-axis"]["elements"]
        
#    print r_freq_dmn

    if os.path.exists('./'+filename+'/' + function_name):
        shutil.rmtree('./'+filename+'/' + function_name)
    os.makedirs('./'+filename+'/' + function_name)
        
    bands = max(function_matrix[:,0]) + 1;
    sites = max(function_matrix[:,4]) + 1;
    times = max(function_matrix[:,5]) + 1;

    f_real = []
    f_imag = []

    for i in range(0,bands):
        for j in range(0,bands):
            for r in range(0,sites):
                for l in range(0,times): 
                    f_real.append(function_matrix[int(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites),6])
                    #f_imag.append(function_matrix[int(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites),7])
                    
                    
                figure(num=None)#, figsize=(10, 10), dpi=200)   
                    
                grid(True)
                plot(r_freq_dmn, f_real, 'k')
                #plot(f_real, 'k')
                
                savefig('./'+filename+'/' + function_name + '/' + function_name + '_{'  + str(r) + '_' + str(i) + '_' + str(j) + '}'+FIGURE_EXTENSION)
                
                f_real = [] 
                f_imag = [] 

def plot_feynamd_expansion_order(domains, function_name, expansion_order):
       if not os.path.exists('./'+filename+'/various'):
            os.makedirs('./'+filename+'/various')

       max_k = max(expansion_order[:,0]) + 1;

       k   = (expansion_order[:,0]).flatten();
       P_k = (expansion_order[:,1]).flatten();
        
       f = figure(num=None)#, figsize=(10, 10), dpi=200)

       subplots_adjust(left, bottom, right, top, wspace, hspace)

       grid(True)

       xlabel('k')
       ylabel('P(k)')
        
       title('diagram expansion order')

       f.text(0.5,0.75,'<k>= ' + str(average(k, 1, P_k)[0,0])[0:5])

       mean = average(k, 1, P_k)[0,0]
       plot(expansion_order[0:int(float(2*mean)),0], expansion_order[0:int(float(2*mean)),1]/sum(P_k)) 
       savefig('./'+filename+'/various/' + function_name+FIGURE_EXTENSION)

def plot_DCA_iteration_function(domains, function_name, function_matrix):
    if not os.path.exists('./'+filename+'/various'):
        os.makedirs('./'+filename+'/various')

    f = figure(num=None)
    grid(True)
    xlabel('DCA-iteration')
    #plot(domains['DCA_iteration_domain'], function_matrix[:,1]) 
    plot(function_matrix[:,1]) 
    savefig('./'+filename+'/various/' + function_name+FIGURE_EXTENSION)

def plot_DCA_iteration_function_log(domains, function_name, function_matrix):
    if not os.path.exists('./'+filename+'/various'):
        os.makedirs('./'+filename+'/various')

    f = figure(num=None)
    grid(True)
    xlabel('DCA-iteration')
    #semilogy(domains['DCA_iteration_domain'], function_matrix[:,1]) 
    semilogy(function_matrix[:,1]) 
    savefig('./'+filename+'/various/' + function_name+FIGURE_EXTENSION)


def plot_Sigma_zero_moment(domains, function_name, function_matrix):
    if not os.path.exists('./'+filename+'/various'):
        os.makedirs('./'+filename+'/various')

    bands   = max(function_matrix[:,0]) + 1;
    spins   = max(function_matrix[:,1]) + 1;
    sites   = max(function_matrix[:,2]) + 1;
    Nb_itrs = max(function_matrix[:,3]) + 1;

    f = figure(num=None)
    grid(True)
    xlabel('DCA-iteration')

    function_matrix_real  = []

    for i in range(0,bands):
        for s in range(0,spins):
            for r in range(0,sites):
                for itrs in range(0,Nb_itrs):
                    function_matrix_real.append(function_matrix[int(i+s*bands+r*spins*bands+itrs*spins*bands*sites),4])
                plot(domains['DCA_iteration_domain'], function_matrix_real) 
                function_matrix_real  = []
                
    savefig('./'+filename+'/various/' + function_name+FIGURE_EXTENSION)

def plot_densities(domains, function_name, function_matrix):
    if not os.path.exists('./'+filename+'/various'):
        os.makedirs('./'+filename+'/various')

    
    bands   = max(function_matrix[:,0]) + 1;
    spins   = max(function_matrix[:,1]) + 1;
    Nb_itrs = max(function_matrix[:,2]) + 1;

    f = figure(num=None)
    grid(True)
    xlabel('DCA-iteration')

    function_matrix_real  = []
    
    for i in range(0,bands):
        for s in range(0,spins):
            for itrs in range(0,Nb_itrs):
                function_matrix_real.append(function_matrix[int(i+s*bands+itrs*spins*bands),3])
            plot(domains['DCA_iteration_domain'], function_matrix_real) 
            function_matrix_real  = []
                
    savefig('./'+filename+'/various/' + function_name+FIGURE_EXTENSION)

def plot_Sigma_r_decay(domains, function_name, function_matrix):

    bands = max(function_matrix[:,0]) + 1;
    sites = max(function_matrix[:,4]) + 1;

    dist=[]
    for k_0 in range(0,sites):

        z=[]
        for l0 in [-1, 0, 1]:
            for l1 in [-1, 0, 1]:

                r = domains["domains"]["LATTICE_SP"]["REAL_SPACE"]["elements"][k_0]
   
                R0 = domains["domains"]["LATTICE_SP"]["REAL_SPACE"]["super-basis"][0]
                R1 = domains["domains"]["LATTICE_SP"]["REAL_SPACE"]["super-basis"][1]

#                 r  = domains["HOST"]["cluster"]["FULL-CLUSTER"]["R"][k_0];
#                 R0 = domains["HOST"]["cluster"]["SUPER-BASIS"]["R"]["R_0"][0];
#                 R1 = domains["HOST"]["cluster"]["SUPER-BASIS"]["R"]["R_1"][0];
                
                x=r[0] + l0*R0[0]+l1*R1[0]
                y=r[1] + l0*R0[1]+l1*R1[1]

                if(x==y):
                    z.append( sqrt(x*x+y*y))

        if(len(z)>0):
            z = sort(z)
            dist.append(z[0])


    #print dist

    function_matrix_real = []
    function_matrix_imag = []
    function_matrix_abs  = []

    figure(num=None)
    grid(True)

    Sigma_re = []
    Sigma_im = []
    Sigma_d  = []            
    Sigma_t = []                
    
    for i in range(0,1):
        j=i
        for r in range(0,sites):

            #r_v  = domains["HOST"]["cluster"]["FULL-CLUSTER"]["R"][r];
            r_v = domains["domains"]["LATTICE_SP"]["REAL_SPACE"]["elements"][r]

            if(r_v[0]==r_v[1]):
                function_matrix_real .append(function_matrix[int(i+j*2*bands+r*4*bands*bands),5])
                function_matrix_imag .append(function_matrix[int(i+j*2*bands+r*4*bands*bands),6])
                function_matrix_abs  .append( sqrt(function_matrix[int(i+j*2*bands+r*4*bands*bands),5]**2 + function_matrix[int(i+j*2*bands+r*4*bands*bands),6]**2) )

#         print function_matrix_abs
#         print transpose([dist,  function_matrix_real,  function_matrix_imag,  function_matrix_abs])

        f = transpose([dist, function_matrix_real])
        f = f[f[:,0].argsort(),]
        f = transpose(f)

        Sigma_d = f[0]
        Sigma_re = f[1]
                
        subplot(221)
        plot(f[0], f[1], 'ok')
        subplot(222)
        semilogy(f[0], abs(f[1]), 'ko')

        f = transpose([dist, function_matrix_imag])
        f = f[f[:,0].argsort(),]
        f = transpose(f)

        subplot(221)
        plot(f[0], f[1], 'or')
        subplot(222)
        semilogy(f[0], abs(f[1]), 'ro')
        
        Sigma_im = f[1]        
                
        f = transpose([dist, function_matrix_abs])
        f = f[f[:,0].argsort(),]
        f = transpose(f)

        Sigma_t = f[1]

        subplot(221)
        plot(f[0],  f[1], 'b--')
        plot(f[0], -f[1], 'b--')
        
        subplot(222)
        semilogy(f[0], f[1], 'b--')


    np.set_printoptions(suppress=False)
    arr = transpose(np.array([Sigma_d, Sigma_re, Sigma_im, Sigma_t]))    
    print arr.tolist()        

    dist=[]
    function_matrix_real = []
    function_matrix_imag = []
    function_matrix_abs  = []

    for i in range(0,1):
        j=i

        for r in range(0,sites):

            x=[]
            y=[]
            z=[]
            w=[]
            for l0 in [-1, 0, 1]:
                for l1 in [-1, 0, 1]:
                
                    r_v = domains["domains"]["LATTICE_SP"]["REAL_SPACE"]["elements"][r]
   
                    R0 = domains["domains"]["LATTICE_SP"]["REAL_SPACE"]["super-basis"][0]
                    R1 = domains["domains"]["LATTICE_SP"]["REAL_SPACE"]["super-basis"][1]

#                     r_v = domains["HOST"]["cluster"]["FULL-CLUSTER"]["R"][r];
#                     R0  = domains["HOST"]["cluster"]["SUPER-BASIS"]["R"]["R_0"][0];
#                     R1  = domains["HOST"]["cluster"]["SUPER-BASIS"]["R"]["R_1"][0];

                    r0 = r_v[0] + l0*R0[0]+l1*R1[0]
                    r1 = r_v[1] + l0*R0[1]+l1*R1[1]

                    if(abs(r0) < 0.00001):
                        
                        x.append(function_matrix[int(i+j*2*bands+r*4*bands*bands),5])
                        y.append(function_matrix[int(i+j*2*bands+r*4*bands*bands),6])
                        z.append( sqrt(function_matrix[int(i+j*2*bands+r*4*bands*bands),5]**2 + function_matrix[int(i+j*2*bands+r*4*bands*bands),6]**2) )
                        w.append( sqrt(r0*r0+r1*r1))

#                         function_matrix_real .append(function_matrix[int(i+j*2*bands+r*4*bands*bands),5])
#                         function_matrix_imag .append(function_matrix[int(i+j*2*bands+r*4*bands*bands),6])
#                         function_matrix_abs  .append( sqrt(function_matrix[int(i+j*2*bands+r*4*bands*bands),5]**2 + function_matrix[int(i+j*2*bands+r*4*bands*bands),6]**2) )

            if(len(w)>0):
                w = sort(w)
                dist.append(w[0])

                function_matrix_real.append(x[0])
                function_matrix_imag.append(y[0])
                function_matrix_abs .append(z[0])
        
        f = transpose([dist, function_matrix_real])
        f = f[f[:,0].argsort(),]
        f = transpose(f)
        
        subplot(223)
        plot(f[0], f[1], 'ok')
        subplot(224)
        semilogy(f[0], abs(f[1]), 'ko')


        f = transpose([dist, function_matrix_imag])
        f = f[f[:,0].argsort(),]
        f = transpose(f)

        subplot(223)
        plot(f[0], f[1], 'or')
        subplot(224)
        semilogy(f[0], abs(f[1]), 'ro')

        
        f = transpose([dist, function_matrix_abs])
        f = f[f[:,0].argsort(),]
        f = transpose(f)
        
        subplot(223)
        plot(f[0],  f[1], 'b--')
        plot(f[0], -f[1], 'b--')
        
        subplot(224)
        semilogy(f[0], f[1], 'b--')


        savefig('./'+filename+'/decay'+FIGURE_EXTENSION)
        savefig('./'+filename+'/various/decay'+FIGURE_EXTENSION)
                
        function_matrix_real  = []
        function_matrix_imag  = []

def plot_phi_r(domains, function_name, function_matrix):
    
    freq_dmn = domains["frequency_domain"]
        
    if not os.path.exists('./'+filename+'/' + function_name):
        os.makedirs('./'+filename+'/' + function_name)

    bands = max(function_matrix[:,0]) + 1;
    sites = max(function_matrix[:,2]) + 1;
    freqs = max(function_matrix[:,3]) + 1;
    
    dist=[]
    for k_0 in range(0,sites):
        r = domains["PCM"]["cluster"]["FULL-CLUSTER"]["R"][k_0];

        R0 = domains["PCM"]["cluster"]["SUPER-BASIS"]["R"]["R_1"][0];
        R1 = domains["PCM"]["cluster"]["SUPER-BASIS"]["R"]["R_2"][0];

        z=[]
        for l0 in [-1, 0, 1]:
            for l1 in [-1, 0, 1]:
                x=r[0] + l0*R0[0]+l1*R1[0]
                y=r[1] + l0*R0[1]+l1*R1[1]
                z.append( sqrt(x*x+y*y))

        z = sort(z)
        dist.append(z[0])

    print dist

    function_matrix_real  = []
    function_matrix_imag  = []

    figure(num=None)
    grid(False)

    for i in range(0,bands):
        j=i
        for l in range(0,10):

            freq_ind = freqs/2+l

            for r in range(0,sites):
                function_matrix_real .append(function_matrix[int(i+j*bands+r*bands*bands+freq_ind*bands*bands*sites),4])
                function_matrix_imag .append(function_matrix[int(i+j*bands+r*bands*bands+freq_ind*bands*bands*sites),5])

            f = transpose([dist, function_matrix_real])
            f = f[f[:,0].argsort(),]
            f = transpose(f)
            
            plot(f[0], f[1], color=my_colormap(float(l)/10.), linestyle='none', marker='o', ms=12.5, label=str(r'$\varpi_{'+str(l+1)+'}$'))

#             f = transpose([dist, function_matrix_imag])
#             f = f[f[:,0].argsort(),]
#             f = transpose(f)
#             plot(f[0], f[1], 'or', label='Im')
                
            function_matrix_real  = []
            function_matrix_imag  = []

    legend(shadow=True, fancybox=True)

    subplots_adjust(left=0.20, right=0.95, top=0.95, bottom=0.125)
    
    xlabel(r'$\vert \vec{r} \vert$', fontsize=25)
    ylabel(r'$ \sum_{\lambda}  \frac{\lambda}{\lambda_0} \phi^{\lambda}(r,\varpi_i) $', fontsize=25)

    savefig('./'+filename+'/' + function_name + '/decay'+FIGURE_EXTENSION)

def plot_Sigma_r_HOST(domains, function_name, function_matrix):
    
    if not os.path.exists('./'+filename+'/' + function_name):
        os.makedirs('./'+filename+'/' + function_name)

    bands = max(function_matrix[:,0]) + 1;
    sites = max(function_matrix[:,4]) + 1;
    
    function_matrix_real  = []
    function_matrix_imag  = []

    x=[]
    y=[]
    for k_0 in range(0,sites):
        r = domains["HOST"]["cluster"]["FULL-CLUSTER"]["R"][k_0];

        R0 = domains["HOST"]["cluster"]["SUPER-BASIS"]["R"]["R_1"][0];
        R1 = domains["HOST"]["cluster"]["SUPER-BASIS"]["R"]["R_2"][0];

        x.append(r[0])
        y.append(r[1])

        function_matrix_real.append(log10(abs(function_matrix[int(k_0*4*bands*bands),5])))
        function_matrix_imag.append(log10(abs(function_matrix[int(k_0*4*bands*bands),6])))

#    plot_contourf(-10, 10, x, y, function_matrix_real, function_name, "Sigma_r_RE", "Sigma_r_RE")
#    plot_contourf(-10, 10, x, y, function_matrix_imag, function_name, "Sigma_r_IM", "Sigma_r_IM")
        

    

    figure()
    grid(False)

    subplot(121)
                                  #D = -min(function_matrix_real)
    scatter(x,y,c=function_matrix_real, s=35)#, vmin=-D, vmax=D)
    colorbar()

    subplot(122)

    #D = -min(function_matrix_imag)
    scatter(x,y,c=function_matrix_imag, s=35)#, vmin=D, vmax=D)
    colorbar()

    savefig('./'+filename+'/' + function_name + '/decay'+FIGURE_EXTENSION)

    plot_contourf_hot(0, 16, x, y, function_matrix_real, function_name, "Sigma_r_RE", "Sigma_r_RE")
    plot_contourf_hot(0, 16, x, y, function_matrix_imag, function_name, "Sigma_r_IM", "Sigma_r_IM")

    function_matrix_real  = []
    function_matrix_imag  = []

    x=[]
    y=[]

def plot_Sigma_k_HOST(domains, function_name, function_matrix):
    
    if not os.path.exists('./'+filename+'/' + function_name):
        os.makedirs('./'+filename+'/' + function_name)

    bands = max(function_matrix[:,0]) + 1;
    sites = max(function_matrix[:,4]) + 1;
    
    vec_re = []
    vec_im = []

    x=[]
    y=[]
    for k_0 in range(0,sites):
        for l0 in [-1, 0, 1]:
            for l1 in [-1, 0, 1]:                         
                kx = domains["domains"]["LATTICE_SP"]["MOMENTUM_SPACE"]["elements"][k_0][0];
                ky = domains["domains"]["LATTICE_SP"]["MOMENTUM_SPACE"]["elements"][k_0][1];
   
                K0 = domains["domains"]["LATTICE_SP"]["MOMENTUM_SPACE"]["super-basis"][0]
                K1 = domains["domains"]["LATTICE_SP"]["MOMENTUM_SPACE"]["super-basis"][1]
            
                x.append(kx + l0*K0[0]+l1*K1[0])
                y.append(ky + l0*K0[1]+l1*K1[1])

                vec_re.append(function_matrix[int(k_0*4*bands*bands),5])
                vec_im.append(function_matrix[int(k_0*4*bands*bands),6])

    print transpose([x, y, vec_re])
    print transpose([x, y, vec_im])

    plot_contourf_jet(-3.141592, 3.141592, x, y, vec_re, function_name, function_name+'_Re'+FIGURE_EXTENSION, r'Re[$\Sigma$]')
    plot_contourf_jet(-3.141592, 3.141592, x, y, vec_im, function_name, function_name+'_Im'+FIGURE_EXTENSION, r'Im[$\Sigma$]')

    #plot_pcolor(-3.141592, 3.141592, x, y, vec_re, function_name, 'Sigma_Re.pdf', 'Re[Sigma]')
    #plot_pcolor(-3.141592, 3.141592, x, y, vec_im, function_name, 'Sigma_Im.pdf', 'Im[Sigma]')

def plot_wn_function(domains, function_name, function_matrix):

    freq_dmn = domains["domains"]["frequency-domain"]["elements"]
        
    if os.path.exists('./'+filename+'/' + function_name):
        shutil.rmtree('./'+filename+'/' + function_name)
    os.makedirs('./'+filename+'/' + function_name)

    bands = max(function_matrix[:,0]) + 1;
    sites = max(function_matrix[:,4]) + 1;
    freqs = max(function_matrix[:,5]) + 1;

    dmn                   = []    
    function_matrix_real  = []
    function_matrix_imag  = []

    function_matrix_real_up  = []
    function_matrix_imag_up  = []
    function_matrix_real_dn  = []
    function_matrix_imag_dn  = []
    
    for r in range(0,sites):
        figure(num=None)
        for i in range(0,bands):
            for j in range(0,bands):

                for l in range(int(freqs/2-64), int(freqs/2+min(64,freqs/2))):
                #for l in range(int(0), int(freqs)):
                    dmn.append(freq_dmn[l])
                    function_matrix_real_up .append(function_matrix[int(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites),6])
                    function_matrix_imag_up .append(function_matrix[int(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites),7])
                    function_matrix_real_dn .append(function_matrix[int(i+bands+(j+bands)*2*bands+r*4*bands*bands+l*4*bands*bands*sites),6])
                    function_matrix_imag_dn .append(function_matrix[int(i+bands+(j+bands)*2*bands+r*4*bands*bands+l*4*bands*bands*sites),7])
                
                subplot(bands, bands, i+1+bands*j)

                #grid(True)
                plot(dmn, function_matrix_real_up, 'k.-', label="Re[ "+function_name+" ]")
                plot(dmn, function_matrix_imag_up, 'r.-', label="Re[ "+function_name+" ]")
#                 plot(dmn, function_matrix_real_dn, 'g.-')
#                 plot(dmn, function_matrix_imag_dn, 'b.-')

                #legend(('real', 'imag'), 'upper right', shadow=True, fancybox=True)
                
                savefig('./'+filename+'/' + function_name + '/' + function_name + '_band_overview_'+ str(r) + FIGURE_EXTENSION)
                
                dmn                   = []
                function_matrix_real  = []
                function_matrix_imag  = []
                function_matrix_real_up  = []
                function_matrix_imag_up  = []
                function_matrix_real_dn  = []
                function_matrix_imag_dn  = []

    figure(num=None)
    for i in range(0,bands):
        j=i
        #for j in range(0,bands):
        for r in range(0,sites):
            for l in range(int(freqs/2), int(freqs/2+min(32,freqs/2))):
                dmn.append(freq_dmn[l])
                function_matrix_real .append(function_matrix[int(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites),6])
                function_matrix_imag .append(function_matrix[int(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites),7])

            grid(True)
            plot(dmn, function_matrix_real, 'k-')
            plot(dmn, function_matrix_imag, 'r-')
            legend(('real', 'imag'), 'upper right', shadow=True, fancybox=True)
                
            dmn                   = []
            function_matrix_real  = []
            function_matrix_imag  = []
    
    savefig('./'+filename+'/' + function_name + '/' + 'overview'+FIGURE_EXTENSION)

    if(False):

        import numpy as np
        import scipy 
        from scipy.interpolate import Rbf
        import matplotlib.pyplot as plt
        from matplotlib import cm

        for i_0 in range(0,bands):
            for i_1 in range(0,bands):
                
                figure(num=None)
                 
                x = []
                y = []
                
                vec_re = []
                vec_im = []
                
                for k_0 in range(0,sites):
                    for l0 in [-1, 0, 1]:
                        for l1 in [-1, 0, 1]:
                            vec_re.append(function_matrix[int(i_0 + i_1*2*bands + k_0*4*bands*bands + freqs/2*4*bands*bands*sites),6])
                            vec_im.append(function_matrix[int(i_0 + i_1*2*bands + k_0*4*bands*bands + freqs/2*4*bands*bands*sites),7])
                            
                            K0 = domains["cluster"]["BASIS"]["K"]["K_1"][0];
                            K1 = domains["cluster"]["BASIS"]["K"]["K_2"][0];
                            
                            x.append(domains["cluster"]["FULL-CLUSTER"]["K"][k_0][0] + l0*K0[0]+l1*K1[0])
                            y.append(domains["cluster"]["FULL-CLUSTER"]["K"][k_0][1] + l0*K0[1]+l1*K1[1])

        #plot_pcolor(-3.141592, 3.141592, x, y, vec_re, function_name, 'Sigma_Re' + str(i_0) + '_' + str(i_1) + FIGURE_EXTENSION, 'Re[Sigma]')
        #plot_pcolor(-3.141592, 3.141592, x, y, vec_im, function_name, 'Sigma_Im' + str(i_0) + '_' + str(i_1) + FIGURE_EXTENSION, 'Im[Sigma]')

        plot_contourf_jet_nearest(-3.141592, 3.141592, x, y, vec_re, function_name, 'Sigma_Re_nearest' + FIGURE_EXTENSION, r'Re[$\Sigma_{QMC}$]')
        plot_contourf_jet_nearest(-3.141592, 3.141592, x, y, vec_im, function_name, 'Sigma_Im_nearest' + FIGURE_EXTENSION, r'Im[$\Sigma_{QMC}$]')

def plot_expansion(domains, function_name, function_matrix):

    bands = max(function_matrix[:,0]) + 1;
    expan = max(function_matrix[:,4]) + 1;
    itera = max(function_matrix[:,5]) + 1;

    print bands
    print expan
    print itera

    function_matrix_real  = []
    function_matrix_imag  = []
    
    figure(figsize=(10,4))

    for it in range(0,itera):
        for i in range(0,1):
            j=i
               
            for l in range(0, expan):
                function_matrix_real .append(function_matrix[int(i+j*2*bands+l*4*bands*bands+it*4*bands*bands*expan),6])
                function_matrix_imag .append(function_matrix[int(i+j*2*bands+l*4*bands*bands+it*4*bands*bands*expan),7])

            z=[]
            for l in range(0, expan):
                z.append(sqrt(function_matrix_real[l]**2 + function_matrix_imag[l]**2))

            print function_matrix_real

            subplot(1,3,1)
            plot(function_matrix_real, '.-', label = "DCA-it : " + str(it))
            
            subplot(1,3,2)
            plot(function_matrix_imag, '.-', label = "DCA-it : " + str(it))

            subplot(1,3,3)
            plot(z, '.-', label = "DCA-it : " + str(it))
     
            function_matrix_real  = []
            function_matrix_imag  = []
           
    legend()
    savefig('./'+filename+'/expansion' + FIGURE_EXTENSION)
     


def plot_local_chi_quantities(domains,function_name, function_matrix):

    freq_dmn = domains["extended_bosonic_vertex_frequency_domain"]
        
    if os.path.exists('./'+filename+'/' + function_name):
        shutil.rmtree('./'+filename+'/' + function_name)
    os.makedirs('./'+filename+'/' + function_name)

    bands = max(function_matrix[:,0]) + 1;
    sites = max(function_matrix[:,4]) + 1;
    freqs = max(function_matrix[:,6]) + 1;

    function_matrix_real  = []
    function_matrix_imag  = []
    
    for b1 in range(0,bands):
        for b2 in range(0,bands):
            for b3 in range(0,bands):
                for b4 in range(0,bands):
                
                    for q in range(0,sites):
                        print q

                        for k in range(0,sites):
                        
                            for m in range(0,freqs):
                                function_matrix_real .append(function_matrix[int(b1
                                                                                 +b2*bands
                                                                                 +b3*bands*bands
                                                                                 +b4*bands*bands*bands

                                                                                 +k *bands*bands*bands*bands
                                                                                 +q *bands*bands*bands*bands*sites
                                                                                 +m *bands*bands*bands*bands*sites*sites),7])
                                function_matrix_imag .append(function_matrix[int(b1
                                                                                 +b2*bands
                                                                                 +b3*bands*bands
                                                                                 +b4*bands*bands*bands

                                                                                 +k *bands*bands*bands*bands
                                                                                 +q *bands*bands*bands*bands*sites
                                                                                 +m *bands*bands*bands*bands*sites*sites),8])
                
                            figure(num=None)#, figsize=(10, 10), dpi=200)


                            print function_matrix_real

                            grid(True)
                            plot(freq_dmn, function_matrix_real, 'k')
                            plot(freq_dmn, function_matrix_imag, 'r')
                            plot(freq_dmn, function_matrix_real, 'ok')
                            plot(freq_dmn, function_matrix_imag, 'or')
                            #legend(('real', 'imag'), 'upper right', shadow=True, fancybox=True)

                            savefig('./'+filename+'/' + function_name + '/' + function_name 
                                    + '_'+ str(b1) 
                                    + '_'+ str(b2) 
                                    + '_'+ str(b3) 
                                    + '_'+ str(b4) 
                                    + '_'+ str(k) + '_'+ str(q) + FIGURE_EXTENSION)
                
                            function_matrix_real  = []
                            function_matrix_imag  = []

def plot_local_susceptibilities(domains,function_name, function_matrix):

    freq_dmn = domains["domains"]["vertex-frequency-domain (EXTENDED_BOSONIC)"]["elements"]
        
    if os.path.exists('./'+filename+'/' + function_name):
        shutil.rmtree('./'+filename+'/' + function_name)
    os.makedirs('./'+filename+'/' + function_name)

    bands = max(function_matrix[:,0]) + 1;
    sites = max(function_matrix[:,4]) + 1;
    freqs = max(function_matrix[:,5]) + 1;

    function_matrix_real  = []
    function_matrix_imag  = []
    
    for b1 in range(0,bands):
        for b2 in range(0,bands):
            for b3 in range(0,bands):
                for b4 in range(0,bands):

                    for r in range(0,sites):
                        for m in range(0,freqs):
                            function_matrix_real .append(function_matrix[int(b1
                                                                             +b2*bands
                                                                             +b3*bands*bands
                                                                             +b4*bands*bands*bands
                                                                             +r *bands*bands*bands*bands
                                                                             +m *bands*bands*bands*bands*sites),6])
                            function_matrix_imag .append(function_matrix[int(b1
                                                                             +b2*bands
                                                                             +b3*bands*bands
                                                                             +b4*bands*bands*bands
                                                                             +r *bands*bands*bands*bands
                                                                             +m *bands*bands*bands*bands*sites),7])
                
                        figure(num=None)#, figsize=(10, 10), dpi=200)

                        grid(True)
                        plot(freq_dmn, function_matrix_real, 'k')
                        plot(freq_dmn, function_matrix_imag, 'r')
                        #plot(freq_dmn, function_matrix_real, 'ok')
                        #plot(freq_dmn, function_matrix_imag, 'or')
                        #legend(('real', 'imag'), 'upper right', shadow=True, fancybox=True)
                        
                        savefig('./'+filename+'/' + function_name + '/' + function_name 
                                + '_'+ str(b1)
                                + '_'+ str(b2)
                                + '_'+ str(b3)
                                + '_'+ str(b4)
                                + '_'+ str(r) + FIGURE_EXTENSION)
                
                        function_matrix_real  = []
                        function_matrix_imag  = []


#     import numpy as np
#     import scipy 
#     from scipy.interpolate import Rbf
#     import matplotlib.pyplot as plt
#     from matplotlib import cm

#     figure(num=None)
    
#     x = []
#     y = []
    
#     vec_re = []

#     print int(freqs)
#     print int((freqs-1)/2)

#     for k_0 in range(0,sites):
#         print function_matrix[int(0
#                                   +0*bands
#                                   +0*bands*bands
#                                   +0*bands*bands*bands
#                                   +k_0 *bands*bands*bands*bands
#                                   +int((freqs-1)/2)*bands*bands*bands*bands*sites),6]

#         for l0 in [-2, -1, 0, 1, 2]:
#             for l1 in [-2, -1, 0, 1, 2]:
#                 vec_re.append(function_matrix[int(0
#                                                   +0*bands
#                                                   +0*bands*bands
#                                                   +0*bands*bands*bands
#                                                   +k_0 *bands*bands*bands*bands
#                                                   +int((freqs-1)/2)*bands*bands*bands*bands*sites),6])
                
#                 K0 = domains["cluster"]["BASIS"]["K"]["K_1"][0];
#                 K1 = domains["cluster"]["BASIS"]["K"]["K_2"][0];

#                 x.append(domains["cluster"]["FULL-CLUSTER"]["K"][k_0][0] + l0*K0[0]+l1*K1[0])
#                 y.append(domains["cluster"]["FULL-CLUSTER"]["K"][k_0][1] + l0*K0[1]+l1*K1[1])

                
# #     print x
# #     print y

#     plot_pcolor(-3.141592, 3.141592, x, y, vec_re, function_name, 'Re_real_local_chi_' + FIGURE_EXTENSION, 'local $\chi$ Re[$\phi(\pi\:T)$]')

def plot_k_cut_spectrum(domains, function_name, function_matrix):

    freq_dmn = domains["Pade_real_frequency_domain"]

    import numpy as np
    import scipy 
    from scipy.interpolate import Rbf
    import matplotlib.pyplot as plt
    from matplotlib import cm

    

    if os.path.exists('./'+filename+'/' + function_name):
        shutil.rmtree('./'+filename+'/' + function_name)
    os.makedirs('./'+filename+'/' + function_name)
    
    bands = 1#max(function_matrix[:,0]) + 1;
    sites = max(function_matrix[:,4]) + 1;
    freqs = max(function_matrix[:,5]) + 1;
    print sites
    print freqs

    x = []
    y = []
    z = []

    for i in range(0,bands):
        j=i
        for j in range(0,bands):

            for k in range(0,sites):
                for w in range(0,freqs):
                    x.append(k)
                    y.append(freq_dmn[w])
                    z.append(-function_matrix[int(i+j*2*bands+k*4*bands*bands+w*4*bands*bands*sites),7])

    xi = np.linspace(1., 192., 192)
    yi = np.linspace(-10., 10., 2048)
    zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='linear')
    
    d = np.linspace(0., 2., 100)

    # contour the gridded data, plotting dots at the randomly spaced data points.
    figure(1)
    plt.contourf(xi,yi,zi,d, cmap=plt.cm.jet)
    plt.colorbar() 

    savefig('./'+filename+'/' + function_name + '/spectrum.pdf')

#     from mpl_toolkits.mplot3d.axes3d import Axes3D
#     import matplotlib.pyplot as plt

# # imports specific to the plots in this example
#     import numpy as np
#     from matplotlib import cm
#     from mpl_toolkits.mplot3d.axes3d import get_test_data

#     if(False):
# # Twice as wide as it is tall.
#         fig = figure(2)
    
# #---- First subplot
#         ax = fig.add_subplot(1, 1, 1, projection='3d')
#         surf = ax.plot_surface(xi, yi, zi, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=False)
#         ax.set_zlim3d(-1.01, 1.01)
#         fig.colorbar(surf, shrink=0.5, aspect=10)
        
#         savefig('./'+filename+'/' + function_name + '/spectrum_3D_A.pdf')

# #---- Second subplot
#     if(False):
#         fig = figure(3)
#         ax = fig.add_subplot(1, 1, 1, projection='3d')
#         ax.plot_wireframe(xi, yi, zi, rstride=10, cstride=10)
#         savefig('./'+filename+'/' + function_name + '/spectrum_3D_B.pdf')


def plot_tau_function(domains, function_name, function_matrix):

    time_dmn = domains["domains"]["time-domain"]["elements"]
        
    if os.path.exists('./'+filename+'/' + function_name):
        shutil.rmtree('./'+filename+'/' + function_name)
    os.makedirs('./'+filename+'/' + function_name)
    
    bands = max(function_matrix[:,0]) + 1;
    sites = max(function_matrix[:,4]) + 1;
    times = max(function_matrix[:,5]) + 1;

    f = []

    for i in range(0,bands):
        #j=i
        for j in range(0,bands):
            for r in range(0,sites):
                for l in range(0,times):
                    f.append(function_matrix[int(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites),6])
                    
                figure(num=None)#, figsize=(10, 10), dpi=200)
                    
                subplots_adjust(left, bottom, right, top,
                                wspace, hspace)
                
                grid(True)
                plot(time_dmn, f) 
                
                savefig('./'+filename+'/' + function_name + '/' + function_name + '_' + str(r) + '_'+ str(i) + '_' + str(j)+FIGURE_EXTENSION)
                
                f = [] 

def plot_phi_wn_quantities(domains, function_name, function_matrix):

    freq_dmn = domains["frequency_domain"]
        
    if os.path.exists('./'+filename+'/' + function_name):
        shutil.rmtree('./'+filename+'/' + function_name)
    os.makedirs('./'+filename+'/' + function_name)
    
    phis  = max(function_matrix[:,0]) + 1;
    bands = max(function_matrix[:,1]) + 1;
    sites = max(function_matrix[:,5]) + 1;
    freqs = max(function_matrix[:,6]) + 1;

    f_real = []
    f_imag = []

    for phi in range(0,phis):
        for i in range(0,bands):
        #j=i
            for j in range(0,bands):
                for r in range(0,sites):
                    for l in range(0,freqs):
                        f_real.append(function_matrix[int(phi + phis*(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites)),7])
                        f_imag.append(function_matrix[int(phi + phis*(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites)),8])

                    figure(num=None)#, figsize=(10, 10), dpi=200)
                    
                    subplots_adjust(left, bottom, right, top, wspace, hspace)
                    
                    grid(True)
                    plot(freq_dmn, f_real, 'k') 
                    plot(freq_dmn, f_imag, 'r') 
                    legend(('real', 'imag'), 'upper right', shadow=True, fancybox=True)
                                   
                    savefig('./'+filename+'/' + function_name + '/' + function_name + '_' + str(phi) + '_' + str(r) + '_'+ str(i) + '_' + str(j)+FIGURE_EXTENSION)
                    
                    f_real = [] 
                    f_imag = [] 

def plot_phi_t_quantities(domains, function_name, function_matrix):

    time_dmn = domains["time_domain"]
        
    if os.path.exists('./'+filename+'/' + function_name):
        shutil.rmtree('./'+filename+'/' + function_name)
    os.makedirs('./'+filename+'/' + function_name)
    
    phis  = max(function_matrix[:,0]) + 1;
    bands = max(function_matrix[:,1]) + 1;
    sites = max(function_matrix[:,5]) + 1;
    times = max(function_matrix[:,6]) + 1;

    f_real = []
    f_imag = []

    for phi in range(0,phis):
        for i in range(0,bands):
        #j=i
            for j in range(0,bands):
                for r in range(0,sites):
                    for l in range(0,times):
                        f_real.append(function_matrix[int(phi + phis*(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites)),7])
                        f_imag.append(function_matrix[int(phi + phis*(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites)),8])

                    figure(num=None)#, figsize=(10, 10), dpi=200)
                    
                    subplots_adjust(left, bottom, right, top, wspace, hspace)
                    
                    grid(True)
                    plot(time_dmn, f_real, 'k') 
                    #plot(time_dmn, f_imag, 'r') 
                    #legend(('real', 'imag'), 'upper right', shadow=True, fancybox=True)
                                   
                    savefig('./'+filename+'/' + function_name + '/' + function_name + '_' + str(phi) + '_' + str(r) + '_'+ str(i) + '_' + str(j)+FIGURE_EXTENSION)
                    
                    f_real = [] 
                    f_imag = [] 

def plot_imaginary_axis(file, line, domains, datamatrix):

    r_freq_dmn = domains["Pade_real_frequency_domain"]
    i_freq_dmn = domains["Pade_imag_frequency_domain"]

    dmn = []

    a_wn = datamatrix
    S_wn = [[]]

    index = 0
    for line in file:

        if "\"S_wn\"" in line:
            stringmatrix = read_function_matrix(file)
            functiondata = json.loads(stringmatrix)
            S_wn = matrix(functiondata)
            index = index + 1

        if(index == 1):
            break

    if os.path.exists('./'+filename+'/a_wn_versus_S_wn/'):
        shutil.rmtree('./'+filename+'/a_wn_versus_S_wn/')
    os.makedirs('./'+filename+'/a_wn_versus_S_wn/')

    bands = max(a_wn[:,0]) + 1;
    sites = max(a_wn[:,4]) + 1;
    times = min(50, max(a_wn[:,5]) + 1);

    a_wn_real = []
    a_wn_imag = []
    S_wn_real = []
    S_wn_imag = []
    S_dmn_new = []
       
    for i in range(0,bands):
        for j in range(0,bands):
            for r in range(0,sites):
                for l in range(0,times): 
                    dmn.append(i_freq_dmn[l])
                    a_wn_real.append(a_wn[int(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites),6])
                    a_wn_imag.append(a_wn[int(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites),7])
                    S_wn_real.append(S_wn[int(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites),6])
                    S_wn_imag.append(S_wn[int(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites),7])
                    
                figure(num=None)
                    
                grid(True)
                plot(dmn, S_wn_real, 'k')
                plot(dmn, S_wn_imag, 'r')
                plot(dmn, a_wn_real, 'k+')
                plot(dmn, a_wn_imag, 'r+')
                legend(('Re[S]', 'Im[S]','Re[a]', 'Im[a]'), 'upper right', shadow=True, fancybox=True)
                
                savefig('./'+filename+'/a_wn_versus_S_wn/a_wn_versus_S_wn_{'  + str(r) + '_' + str(i) + '_' + str(j) + '}'+FIGURE_EXTENSION)
                
                dmn       = []
                a_wn_real = []
                a_wn_imag = []
                S_wn_real = []
                S_wn_imag = []
                S_dmn_new = []


def plot_imaginary_axis_new(file, line, domains, datamatrix):

#     r_freq_dmn = domains["Pade_real_frequency_domain"]
    i_freq_dmn = freq_dmn = domains["domains"]["frequency-domain"]["elements"]#domains["Pade_imag_frequency_domain"]

    dmn = []

    f_approx = datamatrix
    f_source = [[]]

    index = 0
    for line in file:

        if "\"f-measured\"" in line:
            stringmatrix = read_function_matrix(file)
            functiondata = json.loads(stringmatrix)
            f_source = matrix(functiondata["data"])
            index = index + 1

        if(index == 1):
            break

    if os.path.exists('./'+filename+'/f_approx_versus_f_source/'):
        shutil.rmtree('./'+filename+'/f_approx_versus_f_source/')
    os.makedirs('./'+filename+'/f_approx_versus_f_source/')

    bands = max(f_approx[:,0]) + 1;
    sites = max(f_approx[:,4]) + 1;
    times = min(50, max(f_approx[:,5]) + 1);

    f_approx_real = []
    f_approx_imag = []
    f_source_real = []
    f_source_imag = []
    S_dmn_new = []
       
    for i in range(0,bands):
        for j in range(0,bands):
            for r in range(0,sites):
                for l in range(0,times): 
                    dmn.append(i_freq_dmn[l+len(i_freq_dmn)/2])
                    f_approx_real.append(f_approx[int(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites),6])
                    f_approx_imag.append(f_approx[int(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites),7])
                    f_source_real.append(f_source[int(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites),6])
                    f_source_imag.append(f_source[int(i+j*2*bands+r*4*bands*bands+l*4*bands*bands*sites),7])
                    
                figure(num=None)
                    
                grid(True)
                plot(dmn, f_source_real, 'k.-')
                plot(dmn, f_source_imag, 'r.-')
                plot(dmn, f_approx_real, 'k+')
                plot(dmn, f_approx_imag, 'r+')
#                 legend(('Re[f]', 'Im[f]','Re[F] CPE', 'Im[F] CPE'), 'upper right', shadow=True, fancybox=True)
                
                savefig('./'+filename+'/f_approx_versus_f_source/f_approx_versus_f_source_{'  + str(r) + '_' + str(i) + '_' + str(j) + '}'+FIGURE_EXTENSION)
                
                dmn       = []
                f_approx_real = []
                f_approx_imag = []
                f_source_real = []
                f_source_imag = []
                S_dmn_new = []


def plot_band_structure(domains, function_name, function_matrix):

    if not os.path.exists('./'+filename+'/various'):
        os.makedirs('./'+filename+'/various')
    
    bands = max(function_matrix[:,0]) + 1;
    spin  = max(function_matrix[:,1]) + 1;
    sites = max(function_matrix[:,2]) + 1;

    figure(num=None)
    
    for i in range(0,bands):
        f = []
        for r in range(0,sites):
            f.append(function_matrix[int(i+0*bands + r*bands*spin),3])

        grid(True)
        plot(f, 'k') 

    ylabel('Energy (eV)')
    xticks ([0, 100, 200, 300],
            ['$\Gamma$', 'K', 'X', '$\Gamma$'])

    savefig('./'+filename+'/various/' + function_name + FIGURE_EXTENSION)
    
    f = [] 

def plot_S_band_structure(file, line, function_name, function_matrix):

    if not os.path.exists('./'+filename+'/various'):
        os.makedirs('./'+filename+'/various')
    
    bands = max(function_matrix[:,0]) + 1;
    spin  = max(function_matrix[:,1]) + 1;
    sites = max(function_matrix[:,2]) + 1;

    Sigma_qmc = function_matrix

    Sigma_c = [[]]
    Sigma_l = [[]]

    Sigma_interp = [[]]
    Sigma_coarse = [[]]
    
    index = 0
    for line in file:
        
        if "\"Sigma-cluster-band-structure\"" in line:
            print "\"Sigma-cluster-band-structure\""

            stringmatrix = read_function_matrix(file)
            functiondata = json.loads(stringmatrix)
            Sigma_c = matrix(functiondata["data"])
            index = index + 1

        if(index == 1):
            break

    index = 0
    for line in file:
        
        if "\"Sigma-lattice-band-structure\"" in line:
            print "\"Sigma-lattice-band-structure\""

            stringmatrix = read_function_matrix(file)
            functiondata = json.loads(stringmatrix)
            Sigma_l = matrix(functiondata["data"])
            index = index + 1

        if(index == 1):
            break

    index = 0
    for line in file:
        
        if "\"Sigma-band-structure-interpolated\"" in line:
            print "\"Sigma-band-structure-interpolated\""

            stringmatrix = read_function_matrix(file)
            functiondata = json.loads(stringmatrix)
            Sigma_interp = matrix(functiondata["data"])
            index = index + 1

        if(index == 1):
            break

    index = 0
    for line in file:
        
        if "\"Sigma-band-structure-coarsegrained\"" in line:
            print "\"Sigma-band-structure-coarsegrained\""
            
            stringmatrix = read_function_matrix(file)
            functiondata = json.loads(stringmatrix)
            Sigma_coarse = matrix(functiondata["data"])
            index = index + 1

        if(index == 1):
            break

    for i in range(0,bands):

        Sqmc_re = []
        Sqmc_im = []

        Sc_re = []
        Sc_im = []

        Sl_re = []
        Sl_im = []

        S_interp_re = []
        S_interp_im = []

        S_coarse_re = []
        S_coarse_im = []

        for r in range(0,sites):
            Sqmc_re.append(Sigma_qmc[int(i+0*bands + r*bands*spin),3])
            Sqmc_im.append(Sigma_qmc[int(i+0*bands + r*bands*spin),4])

            Sc_re.append(Sigma_c[int(i+0*bands + r*bands*spin),3])
            Sc_im.append(Sigma_c[int(i+0*bands + r*bands*spin),4])

            Sl_re.append(Sigma_l[int(i+0*bands + r*bands*spin),3])
            Sl_im.append(Sigma_l[int(i+0*bands + r*bands*spin),4])

            S_interp_re.append(Sigma_interp[int(i+0*bands + r*bands*spin),3])
            S_interp_im.append(Sigma_interp[int(i+0*bands + r*bands*spin),4])

            S_coarse_re.append(Sigma_coarse[int(i+0*bands + r*bands*spin),3])
            S_coarse_im.append(Sigma_coarse[int(i+0*bands + r*bands*spin),4])

    
        print max(Sl_im)
        mx = 1.1*max([max(Sl_re), max(Sl_im), 0])
        mn = 1.1*min([min(Sl_re), min(Sl_im)])

        if(False):
            #figure(1, figsize=(4,3))
            figure(num=None, figsize=(8,10))

            grid(True)
            
            subplot(311)
            plot(S_interp_re, 'k-', label=r"Re[$\Sigma_{cl}(\vec{k})$]")
            plot(S_interp_im, 'r-', label=r"Im[$\Sigma_{cl}(\vec{k})$]")
            
            plot(S_coarse_re, 'k.', label=r"$\int d\vec{k'} \phi_{k}(k')$"+"Re["+r"$\Sigma_{l}(k')$"+"]", mew=0, markevery=8)
            plot(S_coarse_im, 'r.', label=r"$\int d\vec{k'} \phi_{k}(k')$"+"Im["+r"$\Sigma_{l}(k')$"+"]", mew=0, markevery=8)
            
            plot(Sl_re, 'k--', label=r"Re[$\Sigma_l(\vec{k})$]") 
            plot(Sl_im, 'r--', label=r"Im[$\Sigma_l(\vec{k})$]") 
            
            ylim(mn, mx)
            xlim(0, 404)
            xticks ([0, 101, 202, 303, 404], [r'$(0,0)$', r'$(\pi,\pi)$', r'$(\pi,0)$', r'$(0,\pi)$', r'$(0,0)$'])
            
            legend(loc="upper right")#, bbox_to_anchor=[1.6, 1.0])
            
            subplot(312)
            plot(Sqmc_re    , 'k-', label=r"Re[$\Sigma_{QMC}$]")
            plot(Sqmc_im    , 'r-', label=r"Im[$\Sigma_{QMC}$]")
            plot(S_interp_re, 'k-', label=r"Re[$\Sigma_{cl}$]")
            plot(S_interp_im, 'r-', label=r"Im[$\Sigma_{cl}$]")
            
            ylim(mn, mx)
            xlim(0, 404)
            xticks ([0, 101, 202, 303, 404], [r'$(0,0)$', r'$(\pi,\pi)$', r'$(\pi,0)$', r'$(0,\pi)$', r'$(0,0)$'])
            
            legend(loc="upper right")#, bbox_to_anchor=[1.4, 1.0])
            
            subplot(313)
            plot(Sqmc_re, 'k-', label=r"Re[$\Sigma_{QMC}$]")
            plot(Sqmc_im, 'r-', label=r"Im[$\Sigma_{QMC}$]")

            plot(Sc_re, 'k.', label=r"Re[$\Sigma_{cg}$]", mew=0, markevery=8) 
            plot(Sc_im, 'r.', label=r"Im[$\Sigma_{cg}$]", mew=0, markevery=8)  
            
            ylim(mn, mx)
            xlim(0, 404)
            xticks ([0, 101, 202, 303, 404], [r'$(0,0)$', r'$(\pi,\pi)$', r'$(\pi,0)$', r'$(0,\pi)$', r'$(0,0)$'])
            
            legend(loc="upper right")#, bbox_to_anchor=[1.4, 1.0])
            
            #subplots_adjust(left=0.1, right=0.65, top=0.95, bottom=0.05)

            if(False):
                data = transpose([Sqmc_re, Sqmc_im, Sc_re, Sc_im, Sl_re, Sl_im,S_interp_re, S_interp_im, S_coarse_re, S_coarse_im])
                writer = csv.writer(open("S_cut.txt", "w"), delimiter="\t")
                writer.writerows(data)

            savefig('./'+filename+'/' + function_name + FIGURE_EXTENSION)
            savefig('./'+filename+'/various/' + function_name + FIGURE_EXTENSION)

        if(True):
            figure(1)

            plot(Sqmc_re, 'k-', label=r"Re[$\Sigma_{\vec{K}}$]")
            plot(Sc_re, 'k.', label=r"Re[$\bar{\Sigma}_{\vec{K}}$]", mew=0, markevery=8)
            plot(Sl_re, 'k--', label=r"Re[$\Sigma(\vec{k})$]") 

            plot(Sqmc_im, 'r-', label=r"Im[$\Sigma_{\vec{K}}$]")
            plot(Sc_im, 'r.', label=r"Im[$\bar{\Sigma}_{\vec{K}}$]", mew=0, markevery=8)  
            plot(Sl_im, 'r--', label=r"Im[$\Sigma(\vec{k})$]") 

            #ylim(mn, mx)

            #ylim(-2, 1)
            xlim(0, 399)

            yticks(fontsize=16)
            #xticks([0, 101, 202, 303, 404], [r'$(0,0)$', r'$(\pi,\pi)$', r'$(\pi,0)$', r'$(0,\pi)$', r'$(0,0)$'], fontsize=16, rotation=45)
            xticks([0, 100, 200, 300, 400], [r'$(0,0)$', r'$(\pi,\pi)$', r'$(\pi,0)$', r'$(0,\pi)$', r'$(0,0)$'], fontsize=16, rotation=45)           
            
            subplots_adjust(top=0.95, bottom=0.1, left=0.10, right=0.75, wspace=0.75)

            #legend(loc="upper left", prop={'size':15}, ncol=2, columnspacing=0.2, bbox_to_anchor=[0.0,-0.1])
            legend(loc=(1, 0.49), prop={'size':18}, ncol=1, columnspacing=0.2)
            #legend(loc='center right', bbox_to_anchor=(0.0, 0.0), ncol=2, fancybox=True, shadow=True,prop={'size':15})
            
            savefig('./'+filename+'/Sigma_overview' + FIGURE_EXTENSION)
    

def plot_cut(domains, function_name, function_matrix):

    if not os.path.exists('./'+filename+'/various'):
        os.makedirs('./'+filename+'/various')
    
    k_points = max(function_matrix[:,0]) + 1;

    f_re = []
    f_im = []

    for k in range(0,k_points):
        f_re.append(function_matrix[int(k),1])
        f_im.append(function_matrix[int(k),2])

    figure(figsize=(5,5))

    grid(True)

#     if(function_name == 'S_k_PCM_w_cut'):
#         ylim([-3,1])
    
    if( function_name == 'G0_k_PCM_w_cluster_excluded_cut'):
        ylim([-2,2])

    if(function_name == 'S_k_DCA_w_cut'):
        ylim([-12,12])
    
    plot(f_re, 'k') 
    plot(f_im, 'r') 

    legend(('Re[f] ', 'Im[f] '), 'upper right', shadow=True, fancybox=True)

    ylabel('Energy (eV)')
    xticks ([0, (k_points-1)/3, 2*(k_points-1)/3-1, k_points-1],
            ['[0,0]', '[pi,pi]', '[0,pi]', '[0,0]'])

    savefig('./'+filename+'/various/' + function_name + FIGURE_EXTENSION)
    
    f_re = [] 
    f_im = [] 

def plot_numerical_error(domains, function_name, function_matrix):

    if not os.path.exists('./'+filename+'/various'):
        os.makedirs('./'+filename+'/various')

    bins = domains["numerical_error_domain"];

    average=0
    tot    =0

    error = []
    for l in range(0,len(bins)):
        error.append(function_matrix[int(l),1]+0.1)

        tot     = tot     + function_matrix[int(l),1]
        average = average + bins[l]*function_matrix[int(l),1]

    average = average/tot

    print average

    figure(num=None)

    loglog(bins, error, "ko--")
    
    text(0.000001,0.1*max(error),r'$\epsilon$ = ' + str(my_format(average)))
        
    ylabel('numerical error')
    xlabel(r'$\epsilon$')

    ylim(0.9, 1.1*max(error))

    savefig('./'+filename+'/various/' + function_name + FIGURE_EXTENSION)

def plot_G4_k_k_w_w(function_name, function_matrix):

    if os.path.exists('./'+filename+'/' + function_name):
        shutil.rmtree('./'+filename+'/' + function_name)
    os.makedirs('./'+filename+'/' + function_name)
    
    bands = max(function_matrix[:,0]) + 1;
    sites = max(function_matrix[:,4]) + 1;
    freqs = max(function_matrix[:,6]) + 1;

    print bands
    print sites
    print freqs

    matrixdim = int(freqs)
    print matrixdim

    dim = [];
    dim.append(matrixdim)
    dim.append(matrixdim)
    print dim

    x = []
    y = []

    G4_re = []
    G4_im = []
    G4_lg = []

#     G4_re = zeros([matrixdim, matrixdim]) 
#     G4_im = zeros([matrixdim, matrixdim]) 

#    figure(num=None)
    
    for i_0 in range(0,bands):
        for i_1 in range(0,bands):
            for j_0 in range(0,bands):
                for j_1 in range(0,bands):
                    for k_0 in range(0,sites):
                        for k_1 in range(0,sites):
                            #if((k_0 == 0) and (k_1 == 0)):
                                for w_0 in range(0,freqs):
                                    w_1 = w_0
                                    #for w_1 in range(0,freqs):
                                        #x.append(domains["compact_vertex_frequency_domain"][w_0])
                                        #y.append(domains["compact_vertex_frequency_domain"][w_1])
                                    
                                    re = function_matrix[int(i_0 + bands*i_1 + bands*bands*j_0 + bands*bands*bands*j_1 + bands*bands*bands*bands*k_0 + bands*bands*bands*bands*sites*k_1 + bands*bands*bands*bands*sites*sites*w_0 + bands*bands*bands*bands*sites*sites*freqs*w_1),8]
                                    im = function_matrix[int(i_0 + bands*i_1 + bands*bands*j_0 + bands*bands*bands*j_1 + bands*bands*bands*bands*k_0 + bands*bands*bands*bands*sites*k_1 + bands*bands*bands*bands*sites*sites*w_0 + bands*bands*bands*bands*sites*sites*freqs*w_1),9]
                                        
                                    G4_re.append(re)
                                    G4_im.append(im)

                                        #G4_lg.append(log10(sqrt(re*re+im*im)))
                                
                                #plot_contourf(x[0], -x[0], x, y, G4_re, function_name, 'Re_G4_k1='+str(k_0)+'_k2='+str(k_1)+FIGURE_EXTENSION, 'Re[$G_4$]')
                                #plot_contourf(x[0], -x[0], x, y, G4_im, function_name, 'IM_G4_k1='+str(k_0)+'_k2='+str(k_1)+FIGURE_EXTENSION, 'Im[$G_4$]')
                                #plot_contourf(x[0], -x[0], x, y, G4_lg, function_name, 'Log10_G4'+'_i1='+str(i_0)+'_i2='+str(i_1)+'_j1='+str(j_0)+'_j2='+str(j_1)+'_k1='+str(k_0)+'_k2='+str(k_1)+FIGURE_EXTENSION, '$\log_{10}(|G_4|)$')

                                #x     = []
                                #y     = []
                                
                                figure(num=None)
                    
                                grid(True)
                                plot(G4_re, 'k')
                                plot(G4_im, 'r')
                                legend(('Re[G4]', 'Im[G4]'), 'upper right', shadow=True, fancybox=True)
                                    
                                savefig('./'+filename+ '/' + function_name + '/' + 'G4'+'_i1='+str(i_0)+'_i2='+str(i_1)+'_j1='+str(j_0)+'_j2='+str(j_1)+'_k1='+str(k_0)+'_k2='+str(k_1)+FIGURE_EXTENSION)

                                G4_re = []
                                G4_im = []
                                G4_lg = []
                            
  
def plot_leading_eigenvalues(domains,function_name, function_matrix):
    
    if( not (os.path.exists('./'+filename+'/various'))):
        os.makedirs('./'+filename+'/various')

    nb_eigenvals = max(function_matrix[:,0]) + 1;

    figure(num=None)

    val_re = []
    val_im = []

    for i in range(0,nb_eigenvals):
        val_re.append(function_matrix[int(i),1])
        val_im.append(function_matrix[int(i),2])

    grid(True)
    plot(val_re, 'k')
    plot(val_im, 'r')

    savefig('./'+filename+'/various/' + function_name + FIGURE_EXTENSION)

def plot_leading_eigenvectors(domains,function_name, function_matrix):

    import numpy as np
    import scipy 
    from scipy.interpolate import Rbf
    import matplotlib.pyplot as plt
    from matplotlib import cm
    
    if os.path.exists('./'+filename+'/' + function_name):
        shutil.rmtree('./'+filename+'/' + function_name)
    os.makedirs('./'+filename+'/' + function_name)

#    if os.path.exists('./'+filename+'/wave_function/'):
#shutil.rmtree('./'+filename+'/wave_function/')
#    os.makedirs('./'+filename+'/wave_function/')

    #freq_dmn = domains["compact_vertex_frequency_domain"]

    freq_dmn = domains["domains"]["vertex-frequency-domain (COMPACT)"]["elements"]

    Number = max(function_matrix[:,0]) + 1;    
    bands = max(function_matrix[:,1]) + 1;
    sites = max(function_matrix[:,3]) + 1;
    freqs = max(function_matrix[:,4]) + 1;

    host_cluster=False

    print domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"]

    if( len(domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"])!=sites):
        host_cluster=True

    for Nb in range(0, min(3,Number)):
                    
        for i_0 in range(0,bands):
            i_1=i_0

            figure(num=None)
            grid(True)

            for  k_0 in range(0,sites):
                   
                vec_re = []
                vec_im = []
                    
                for w_0 in range(0,freqs):
                    vec_re.append(function_matrix[int(Nb + Number*(i_0 + i_1*bands + k_0*bands*bands + w_0*bands*bands*sites)),5])
                    vec_im.append(function_matrix[int(Nb + Number*(i_0 + i_1*bands + k_0*bands*bands + w_0*bands*bands*sites)),6])
                        
                subplot(1,2,1)
                plot(freq_dmn, vec_re, '-')#, label=str(my_format(domains["cluster"]["FULL-CLUSTER"]["K"][k_0][0])+","+my_format(domains["cluster"]["FULL-CLUSTER"]["K"][k_0][1])))
                    #plot(freq_dmn, vec_im, label=str(my_format(domains["cluster"]["FULL-CLUSTER"]["K"][k_0][0])+","+my_format(domains["cluster"]["FULL-CLUSTER"]["K"][k_0][1])))

                subplot(1,2,2)
                plot(freq_dmn, vec_im, '-')

                legend(loc="upper right")

                ylabel('$\phi(k, \omega_n)$',fontsize=16)
                xlabel('$\omega_n$',fontsize=16)
                                    
            savefig('./'+filename+'/'+ function_name +'/eigenvectors_' + str(Nb) + FIGURE_EXTENSION)

    #savefig('./'+filename+'/' + function_name +'/' + function_name + '_' + str(Nb) + '_' + str(i_0) + '_' + str(i_1) + '_' + str(k_0) + FIGURE_EXTENSION)

        freq_ind = freqs/2
        for i_0 in range(0,bands):
            for i_1 in range(0,bands):

                i = 0.
                x = []
            
                vec_re = []
                vec_im = []
                

                figure(num=None)

                kxvec=[]
                Gapvec=[]
                
                for k_0 in range(0,sites):
                        
                    freq_ind = freqs/2
                        
                    if(host_cluster):
                        K0 = domains["domains"]["LATTICE_TP"]["MOMENTUM_SPACE"]["elements"][k_0][0];
                        K1 = domains["domains"]["LATTICE_TP"]["MOMENTUM_SPACE"]["elements"][k_0][1];
                    else:
                        K0 = domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"][k_0][0];
                        K1 = domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"][k_0][1];
                                               
                    if(abs(K0-3.141592)<0.001 and abs(K1-0)<0.001):
                        freq_ind = freqs/2
                        maxvec = function_matrix[int(Nb + Number*(i_0 + i_1*bands + k_0*bands*bands + freq_ind*bands*bands*sites)),5]
                
                for k_0 in range(0,sites):

                    if(host_cluster):
                        K0 = domains["domains"]["LATTICE_TP"]["MOMENTUM_SPACE"]["elements"][k_0][0];
                        K1 = domains["domains"]["LATTICE_TP"]["MOMENTUM_SPACE"]["elements"][k_0][1];
                    else:
                        K0 = domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"][k_0][0];
                        K1 = domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"][k_0][1];
                                   
                    if(abs(K0+K1-3.141592)<0.003):
                        x.append(i)
                        i = i+1.

                        freq_ind = freqs/2
                        
                        vec_re.append(function_matrix[int(Nb + Number*(i_0 + i_1*bands + k_0*bands*bands + freq_ind*bands*bands*sites)),5]/maxvec)
                        vec_im.append(function_matrix[int(Nb + Number*(i_0 + i_1*bands + k_0*bands*bands + freq_ind*bands*bands*sites)),6]/maxvec)
                    
                        vec_w_re = []
                        vec_w_im = []
                        
                        for w_0 in range(0,freqs):
                            vec_w_re.append(function_matrix[int(Nb + Number*(i_0 + i_1*bands + k_0*bands*bands + w_0*bands*bands*sites)),5]/maxvec)
                            vec_w_im.append(function_matrix[int(Nb + Number*(i_0 + i_1*bands + k_0*bands*bands + w_0*bands*bands*sites)),6]/maxvec)
                        
                        kxvec .append(K0)
                        Gapvec.append(vec_w_re)

                gapmin=0
                gapmax=0
                for l0 in range(0, len(kxvec)):                
                    
                    if(gapmin>min(Gapvec[l0])):
                        gapmin=min(Gapvec[l0])

                    if(gapmax<max(Gapvec[l0])):
                        gapmax=max(Gapvec[l0])
                    
                    plot(freq_dmn, Gapvec[l0], 'o-', color=plt.cm.jet((l0+0.0)/((len(kxvec)-1.))), label=r'$\vec{k}$ = ('+str(my_format(kxvec[l0]))+', '+str(my_format(3.141592-kxvec[l0]))+')')

                gapdelta   = (gapmax-gapmin)/2                      
                gapaverage = (gapmax+gapmin)/2                      

                ylim(gapaverage-1.1*gapdelta, gapaverage+1.1*gapdelta)            
                 
                #ylim(0.5, 1.1)
                        
                for l0 in range(0, len(x)):            
                    x[l0] = (x[l0])/(len(x)-1)            
                                                    
                legend(loc="upper right")
            
                xticks(fontsize=16)
                yticks(fontsize=16)
                ylabel(r'$\Phi(\vec{k}, \omega_n)$',fontsize=25, labelpad=10)
                xlabel(r'$\omega_n$'               ,fontsize=25, labelpad=10)
            
                
                a = axes([0.225, 0.7, .225, .225], axisbg='w')
                plot(x, vec_re, 'ko-')
                title(r'$\Phi(\vec{k}, \pi\,T)$')
                xticks( [0, 0.25, 0.5, 0.75, 1.], ('($\pi$,$0$)', '($\\frac{\pi}{4}$, $\\frac{3\pi}{4}$)', '($\\frac{\pi}{2}$, $\\frac{\pi}{2}$)', '($\\frac{3\pi}{4}$, $\\frac{\pi}{4}$)', '($\pi$, $0$)'), rotation=45 )
                setp(a, xlim=(0,1.))

                subplots_adjust(left=0.15, right=0.985, top=0.985, bottom=0.115)
                savefig('./'+filename+'/' + function_name +'/Gap_function_' + str(Nb) + '_' + str(i_0) + '_' + str(i_1) + FIGURE_EXTENSION)

        if(True):
            import numpy as np
            import scipy 
            from scipy.interpolate import Rbf
            import matplotlib.pyplot as plt
            from matplotlib import cm

            freq_ind = freqs/2

            for i_0 in range(0,bands):
                for i_1 in range(0,bands):

                    figure(num=None)
                
                    x = []
                    y = []
                
                    vec_re = []
                    vec_im = []

                    for k_0 in range(0,sites):
                        for l0 in [-1, 0, 1]:
                            for l1 in [-1, 0, 1]:
                                vec_re.append(function_matrix[int(Nb + Number*(i_0 + i_1*bands + k_0*bands*bands + freq_ind*bands*bands*sites)),5])
                                vec_im.append(function_matrix[int(Nb + Number*(i_0 + i_1*bands + k_0*bands*bands + freq_ind*bands*bands*sites)),6])

                                if(host_cluster):
                                    K0 = domains["domains"]["LATTICE_TP"]["MOMENTUM_SPACE"]["super-basis"][0];
                                    K1 = domains["domains"]["LATTICE_TP"]["MOMENTUM_SPACE"]["super-basis"][1];

                                    x.append(domains["domains"]["LATTICE_TP"]["MOMENTUM_SPACE"]["elements"][k_0][0] + l0*K0[0]+l1*K1[0])
                                    y.append(domains["domains"]["LATTICE_TP"]["MOMENTUM_SPACE"]["elements"][k_0][1] + l0*K0[1]+l1*K1[1])
                                else:
                                    K0 = domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["super-basis"][0];
                                    K1 = domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["super-basis"][1];

                                    x.append(domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"][k_0][0] + l0*K0[0]+l1*K1[0])
                                    y.append(domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"][k_0][1] + l0*K0[1]+l1*K1[1])

                    xticks( arange(3), ('$-\pi$', '$0$', '$\pi$') )
                    yticks( arange(3), ('$-\pi$', '$0$', '$\pi$') )
                
                    plot_pcolor(-3.141592, 3.141592, x, y, vec_re, function_name, 'Re_wavefunction_' + str(Nb) + '_' + str(i_0) + '_' + str(i_1) + FIGURE_EXTENSION, 'leading eigenvector Re[$\phi(\pi/T)$]')
                    plot_pcolor(-3.141592, 3.141592, x, y, vec_im, function_name, 'Im_wavefunction_' + str(Nb) + '_' + str(i_0) + '_' + str(i_1) + FIGURE_EXTENSION, 'leading eigenvector Im[$\phi(\pi/T)$)]')

def fit_Lorentz_curve(phi_x, phi_y):
    
    print "\n\t\ starting to fit \n\n"

    fitfunc = lambda p, x: p[0]+p[1]/(pow(p[2],2)+pow(x,2)) # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function

    Xlist = np.array(phi_x)
    Ylist = np.array(phi_y)

    p0 = [0., 1., 1.] # Initial guess for the parameters                                                                                                                   
    p1, success = optimize.leastsq(errfunc, p0[:], args=(Xlist, Ylist))
    
    print p1[0]
    print p1[1]
    print p1[2]
    
    Xfine = linspace(min(phi_x), max(phi_x), 100)
    Yfine = fitfunc(p1, Xfine)
    
    return [p1[0], p1[1], p1[2], Xfine, Yfine]


def fit_leading_eigenvectors(domains,function_name, function_matrix):
    
    print "\n\t\ starting to plot \n\n"

    if os.path.exists('./'+filename+'/' + function_name):
        shutil.rmtree('./'+filename+'/' + function_name)
    os.makedirs('./'+filename+'/' + function_name)

    freq_dmn = domains["domains"]["vertex-frequency-domain (COMPACT)"]["elements"]

    Number = max(function_matrix[:,0]) + 1;    
    bands = max(function_matrix[:,1]) + 1;
    sites = max(function_matrix[:,3]) + 1;
    freqs = max(function_matrix[:,4]) + 1;

    host_cluster=False

    print domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"]

    if( len(domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"])!=sites):
        host_cluster=True

    for Nb in range(0, min(3,Number)):
                    
        freq_ind = freqs/2
        for i_0 in range(0,bands):
            for i_1 in range(0,bands):

                i = 0.
                x = []
            
                vec_re = []
                vec_im = []
                

                figure(num=None)

                kxvec      = []
                Gapvec     = []
                fullgapvec = []

                for k_0 in range(0,sites):
                        
                    freq_ind = freqs/2
                        
                    if(host_cluster):
                        K0 = domains["domains"]["LATTICE_TP"]["MOMENTUM_SPACE"]["elements"][k_0][0];
                        K1 = domains["domains"]["LATTICE_TP"]["MOMENTUM_SPACE"]["elements"][k_0][1];
                    else:
                        K0 = domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"][k_0][0];
                        K1 = domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"][k_0][1];
                                               
                    if(abs(K0-3.141592)<0.001 and abs(K1-0)<0.001):
                        freq_ind = freqs/2
                        maxvec = max(abs(function_matrix[int(Nb + Number*(i_0 + i_1*bands + k_0*bands*bands + freq_ind*bands*bands*sites)),5]),
                                     abs(function_matrix[int(Nb + Number*(i_0 + i_1*bands + k_0*bands*bands + freq_ind*bands*bands*sites)),6]))

                for k_0 in range(0,sites):

                    if(host_cluster):
                        K0 = domains["domains"]["LATTICE_TP"]["MOMENTUM_SPACE"]["elements"][k_0][0];
                        K1 = domains["domains"]["LATTICE_TP"]["MOMENTUM_SPACE"]["elements"][k_0][1];
                    else:
                        K0 = domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"][k_0][0];
                        K1 = domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"][k_0][1];
                                   
                    if(abs(K0+K1-3.141592)<0.003):
                        x.append(i)
                        i = i+1.

                        freq_ind = freqs/2
                        
                        vec_re.append(function_matrix[int(Nb + Number*(i_0 + i_1*bands + k_0*bands*bands + freq_ind*bands*bands*sites)),5]/maxvec)
                        vec_im.append(function_matrix[int(Nb + Number*(i_0 + i_1*bands + k_0*bands*bands + freq_ind*bands*bands*sites)),6]/maxvec)
                    
                        vec_w_re = []
                        vec_w_im = []
                        
                        for w_0 in range(0,freqs):
                            vec_w_re.append(function_matrix[int(Nb + Number*(i_0 + i_1*bands + k_0*bands*bands + w_0*bands*bands*sites)),5]/maxvec)
                            vec_w_im.append(function_matrix[int(Nb + Number*(i_0 + i_1*bands + k_0*bands*bands + w_0*bands*bands*sites)),6]/maxvec)
                        
                        kxvec .append(K0)
                        Gapvec.append(vec_w_re)

                gapmin=0
                gapmax=0
                for l0 in range(0, len(kxvec)):                
                    
                    if(gapmin>min(Gapvec[l0])):
                        gapmin=min(Gapvec[l0])

                    if(gapmax<max(Gapvec[l0])):
                        gapmax=max(Gapvec[l0])
                    
                    tmpx = freq_dmn
                    tmpy = Gapvec[l0]

                    [a,b,c, Xfine, Yfine] = fit_Lorentz_curve(tmpx, tmpy)

                    plot(freq_dmn, Gapvec[l0], 'o', color=plt.cm.jet((l0+0.0)/((len(kxvec)-1.))), label=r'$\vec{k}$ = ('+str(my_format(kxvec[l0]))+', '+str(my_format(3.141592-kxvec[l0]))+')')                
                    plot(Xfine, Yfine, '--', color=plt.cm.jet((l0+0.0)/((len(kxvec)-1.))), label=r"$\gamma=$"+my_format(c))



                gapdelta   = (gapmax-gapmin)/2                      
                gapaverage = (gapmax+gapmin)/2                      

                ylim(gapaverage-1.1*gapdelta, gapaverage+1.1*gapdelta)            
                 
                #ylim(0.5, 1.1)
                        
                for l0 in range(0, len(x)):            
                    x[l0] = (x[l0])/(len(x)-1)            
                                                    
                legend(loc="upper right", ncol=2, prop={'size':6})
            
                xticks(fontsize=16)
                yticks(fontsize=16)
                ylabel(r'$\Phi(\vec{k}, \omega_n)$',fontsize=25, labelpad=10)
                xlabel(r'$\omega_n$'               ,fontsize=25, labelpad=10)
                            
                a = axes([0.225, 0.7, .225, .225], axisbg='w')
                plot(x, vec_re, 'ko-')
                title(r'$\Phi(\vec{k}, \pi\,T)$')
                xticks( [0, 0.25, 0.5, 0.75, 1.], ('($\pi$,$0$)', '($\\frac{\pi}{4}$, $\\frac{3\pi}{4}$)', '($\\frac{\pi}{2}$, $\\frac{\pi}{2}$)', '($\\frac{3\pi}{4}$, $\\frac{\pi}{4}$)', '($\pi$, $0$)'), rotation=45 )
                setp(a, xlim=(0,1.))

                subplots_adjust(left=0.15, right=0.985, top=0.985, bottom=0.115)
                savefig('./'+filename+'/' + function_name +'/Gap_function_' + str(Nb) + '_' + str(i_0) + '_' + str(i_1) + FIGURE_EXTENSION)

def plot_chi_0(domains,function_name, function_matrix):
    
    if not os.path.exists('./'+filename+'/'+'various/'):
        os.makedirs('./'+filename+'/'+'various/')

    freq_dmn = domains["domains"]["vertex-frequency-domain (COMPACT)"]["elements"]

    bands = max(function_matrix[:,0]) + 1;
    sites = max(function_matrix[:,2]) + 1;
    freqs = max(function_matrix[:,3]) + 1;
               
    host_cluster=False
    if( len(domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"])!=sites):
        host_cluster=True
     
    for i_0 in range(0,bands):
        i_1=i_0

        figure(num=None)
        grid(True)

        for  k_0 in range(0,sites):
                   
            vec_re = []
            vec_im = []
                        
            for w_0 in range(0,freqs):
                vec_re.append(function_matrix[int(i_0 + i_1*bands + k_0*bands*bands + w_0*bands*bands*sites),4])
                vec_im.append(function_matrix[int(i_0 + i_1*bands + k_0*bands*bands + w_0*bands*bands*sites),5])
                
#             K0 = domains["HOST"]["cluster"]["FULL-CLUSTER"]["K"][k_0][0]
#             K1 = domains["HOST"]["cluster"]["FULL-CLUSTER"]["K"][k_0][1]
#             if( abs(K0+K1-3.141592)<0.001):
#                 plot(freq_dmn, vec_re, '-', label = str(K0)+" ; "+str(K1))
#                 legend()
#                 ylabel('$\phi(k, \omega_n)$',fontsize=25)
#                 xlabel('$\omega_n$',fontsize=25)
#             #legend(('k=[$\pi$,$0$]','k=[$\\frac{3\pi}{4}$, $\\frac{\pi}{4}$]','k=[$\\frac{\pi}{2}$, $\\frac{\pi}{2}$]', 'k=[$\\frac{\pi}{4}$, $\\frac{3\pi}{4}$]','k=[$0$,$\pi$]'), 'upper right', shadow=True, fancybox=True)
                
            plot(freq_dmn, vec_re)
            savefig('./'+filename+'/'+'various/'+ function_name  + FIGURE_EXTENSION)

        if(True):
            import numpy as np
            import scipy 
            from scipy.interpolate import Rbf
            import matplotlib.pyplot as plt
            from matplotlib import cm

            freq_ind = freqs/2

            for i_0 in range(0,bands):
                for i_1 in range(0,bands):

                    figure(num=None)
                
                    x = []
                    y = []
                
                    vec_re = []

                    for k_0 in range(0,sites):
                        for l0 in [-1, 0, 1]:
                            for l1 in [-1, 0, 1]:
                                vec_re.append(function_matrix[int(i_0 + i_1*bands + k_0*bands*bands + freq_ind*bands*bands*sites),4])
                                vec_im.append(function_matrix[int(i_0 + i_1*bands + k_0*bands*bands + freq_ind*bands*bands*sites),5])

                                if(host_cluster):
                                    K0 = domains["domains"]["LATTICE_TP"]["MOMENTUM_SPACE"]["super-basis"][0];
                                    K1 = domains["domains"]["LATTICE_TP"]["MOMENTUM_SPACE"]["super-basis"][1];

                                    x.append(domains["domains"]["LATTICE_TP"]["MOMENTUM_SPACE"]["elements"][k_0][0] + l0*K0[0]+l1*K1[0])
                                    y.append(domains["domains"]["LATTICE_TP"]["MOMENTUM_SPACE"]["elements"][k_0][1] + l0*K0[1]+l1*K1[1])
                                else:
                                    K0 = domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["super-basis"][0];
                                    K1 = domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["super-basis"][1];

                                    x.append(domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"][k_0][0] + l0*K0[0]+l1*K1[0])
                                    y.append(domains["domains"]["CLUSTER"]["MOMENTUM_SPACE"]["elements"][k_0][1] + l0*K0[1]+l1*K1[1])

                    xticks( arange(3), ('$-\pi$', '$0$', '$\pi$') )
                    yticks( arange(3), ('$-\pi$', '$0$', '$\pi$') )
                
                    plot_pcolor(-3.141592, 3.141592, x, y, vec_re, 'various', function_name + "_k_re" + FIGURE_EXTENSION, 'chi_0')
                    #plot_pcolor(-3.141592, 3.141592, x, y, vec_im, 'various', function_name + "_k_im" + FIGURE_EXTENSION, 'chi_0')


def plot_wave_function(domains,function_name, function_matrix):
    
    import numpy as np
    import scipy 
    from scipy.interpolate import Rbf
    import matplotlib.pyplot as plt
    from matplotlib import cm

    if os.path.exists('./'+filename+'/wave_function/'):
        shutil.rmtree('./'+filename+'/wave_function/')
    os.makedirs('./'+filename+'/wave_function/')

    bands = max(function_matrix[:,0]) + 1;
    sites = max(function_matrix[:,2]) + 1;

    for i_0 in range(0,bands):
        for i_1 in range(0,bands):

            figure(num=None)
                
            x = []
            y = []

            vec_re = []
            vec_im = []

            for k_0 in range(0,sites):
                for l0 in [-2, -1, 0, 1, 2]:
                    for l1 in [-2, -1, 0, 1, 2]:
                        vec_re.append(function_matrix[int(i_0 + i_1*bands + k_0*bands*bands),3])
                        vec_im.append(function_matrix[int(i_0 + i_1*bands + k_0*bands*bands),4])

                        K0 = domains["cluster"]["BASIS"]["K"]["K_1"][0];
                        K1 = domains["cluster"]["BASIS"]["K"]["K_2"][0];

                        x.append(domains["cluster"]["FULL-CLUSTER"]["K"][k_0][0] + l0*K0[0]+l1*K1[0])
                        y.append(domains["cluster"]["FULL-CLUSTER"]["K"][k_0][1] + l0*K0[1]+l1*K1[1])

#             print x
#             print y

            plot_pcolor(-3.141592, 3.141592, x, y, vec_re, 'wave_function', 'Re_wavefunction_' + str(i_0) + '_' + str(i_1) + FIGURE_EXTENSION, 'leading eigenvector Re[$\phi(\pi\:T)$]')
            plot_pcolor(-3.141592, 3.141592, x, y, vec_im, 'wave_function', 'Im_wavefunction_' + str(i_0) + '_' + str(i_1) + FIGURE_EXTENSION, 'leading eigenvector Im[$\phi(\pi\:T)$]')

def plot_f_k_function(domains,function_name, function_matrix):    

    import numpy as np
    import scipy 
    from scipy.interpolate import Rbf
    import matplotlib.pyplot as plt
    from matplotlib import cm

    if not os.path.exists('./'+filename+'/various'):
        os.makedirs('./'+filename+'/various')

    bands = max(function_matrix[:,0]) + 1;
    sites = max(function_matrix[:,2]) + 1;
    DCAit = max(function_matrix[:,3]) + 1;

    for i_0 in range(0,bands):

        figure(num=None)
                
        x = []
        y = []

        vec_re = []
        vec_im = []

        for k_0 in range(0,sites):
            for l0 in [-2, -1, 0, 1, 2]:
                for l1 in [-2, -1, 0, 1, 2]:
                    vec_re.append(function_matrix[int(i_0 + 0*bands + k_0*2*bands + (DCAit-1)*2*bands*sites),4])
                    vec_im.append(function_matrix[int(i_0 + 0*bands + k_0*2*bands + (DCAit-1)*2*bands*sites),5])

                    K0 = domains["cluster"]["BASIS"]["K"]["K_0"][0];
                    K1 = domains["cluster"]["BASIS"]["K"]["K_1"][0];

                    x.append(domains["cluster"]["FULL-CLUSTER"]["K"][k_0][0] + l0*K0[0]+l1*K1[0])
                    y.append(domains["cluster"]["FULL-CLUSTER"]["K"][k_0][1] + l0*K0[1]+l1*K1[1])

#             print x
#             print y

        plot_pcolor(-3.141592, 3.141592, x, y, vec_re, 'various', function_name + '_' + str(i_0) + FIGURE_EXTENSION, 'A(k)')
#         plot_pcolor(-3.141592, 3.141592, x, y, vec_im, './'+filename+'/various', 'Im_wavefunction_' + str(i_0) + '_' + str(i_1) + FIGURE_EXTENSION, 'leading eigenvector Im[$\phi(\pi\:T)$]')

def plot_f_k_HOST_function(domains,function_name, function_matrix, real_or_imag):    

    import numpy as np
    import scipy 
    from scipy.interpolate import Rbf
    import matplotlib.pyplot as plt
    from matplotlib import cm

    if not os.path.exists('./'+filename+'/various'):
        os.makedirs('./'+filename+'/various')

    bands = max(function_matrix[:,0]) + 1;
    sites = max(function_matrix[:,4]) + 1;

    for i_0 in range(0,1):

        figure(num=None)
                
        x = []
        y = []

        vec = []

        for k_0 in range(0,sites):
            for l0 in [-1,0,1]:
                for l1 in [-1,0,1]:
                    if(real_or_imag=='real'):
                        vec.append( function_matrix[int(k_0*2*bands*2*bands),5])
                    else:
                        vec.append(-function_matrix[int(k_0*2*bands*2*bands),6])

                    K0 = domains["HOST"]["cluster"]["BASIS"]["K"]["K_0"][0];
                    K1 = domains["HOST"]["cluster"]["BASIS"]["K"]["K_1"][0];

                    x.append(domains["HOST"]["cluster"]["FULL-CLUSTER"]["K"][k_0][0] + l0*K0[0]+l1*K1[0])
                    y.append(domains["HOST"]["cluster"]["FULL-CLUSTER"]["K"][k_0][1] + l0*K0[1]+l1*K1[1])

        plot_contourf_hot(-3.141592, 3.141592, x, y, vec, '/various', function_name + FIGURE_EXTENSION, function_name)


def plot_n_k_HOST_function(domains,function_name, function_matrix):    

    import numpy as np
    import scipy 
    from scipy.interpolate import Rbf
    import matplotlib.pyplot as plt
    from matplotlib import cm

    if not os.path.exists('./'+filename+'/various'):
        os.makedirs('./'+filename+'/various')

    bands = max(function_matrix[:,0]) + 1;
    sites = max(function_matrix[:,2]) + 1;

    for i_0 in range(0,1):

        figure(num=None)
                
        x = []
        y = []

        vec_re = []
        vec_im = []

        for k_0 in range(0,sites):
            for l0 in [-1,0,1]:
                for l1 in [-1,0,1]:
                    vec_re.append( function_matrix[int(k_0*2*bands),3])
                    vec_im.append(-function_matrix[int(k_0*2*bands),4])

                    K0 = domains["HOST"]["cluster"]["BASIS"]["K"]["K_0"][0];
                    K1 = domains["HOST"]["cluster"]["BASIS"]["K"]["K_1"][0];

                    x.append(domains["HOST"]["cluster"]["FULL-CLUSTER"]["K"][k_0][0] + l0*K0[0]+l1*K1[0])
                    y.append(domains["HOST"]["cluster"]["FULL-CLUSTER"]["K"][k_0][1] + l0*K0[1]+l1*K1[1])

        plot_contourf_hot(-3.141592, 3.141592, x, y, vec_re, '/various', function_name + FIGURE_EXTENSION, function_name)


print '\n \t start reading ' + filename  + ' \n'

file = open(filename, 'r')

print filename

filename = filename.rstrip('.json')

print filename

function_names = []

domains_data = ''

for line in file:

    if not ( ("\"functions\"" in line) or ("\"spectral-functions\"" in line) or ("\"CPE-functions\"" in line) ):
        domains_data = domains_data + line
    else:
        #domains_data_stri = domains_data[0:len(domains_data)-2] + '}'
        domains_data_stri = domains_data[0:len(domains_data)-2] + "\"tmp\" : 1 }"

        #print domains_data_stri

        domains = json.loads(domains_data_stri)

        for line in file:

            if "\"name\" : " in line:

                function_name = read_function_name(line)

                yes_or_no = raw_input('\n\t do you want to plot ' + function_name + '? (y/n) : ')

                if(yes_or_no == 'y'):

                    datamatrix = [[]]

                    stringmatrix = read_function_matrix(file)
                    functiondata = json.loads(stringmatrix)
                    datamatrix   = matrix(functiondata["data"])

                    if(function_name == 'sigma_lambda'):
                        plot_expansion(domains, function_name, datamatrix)

                    if(function_name == 'S_k_PCM_w_cut'
                       or function_name == 'alpha_k_cut'
                       or function_name == 'S_k_DCA_w_cut'
                       or function_name == 'G0_k_PCM_w_cluster_excluded_cut'):
                        plot_cut(domains, function_name, datamatrix)

                    if(function_name == 'Self_Energy'
                       or function_name == 'Self_Energy_r_w'
                       or function_name == 'cluster_greens_function_G_k_w'
                       or function_name == 'cluster_greens_function_G_r_w'
                       or function_name == 'free_cluster_greens_function_G0_k_w'
                       or function_name == 'free_cluster_greens_function_G0_r_w'
                       or function_name == 'cluster_excluded_greens_function_G0_k_w'
                       or function_name == 'cluster_excluded_greens_function_G0_r_w'
                       or function_name == 'M_r_w'
                       or function_name == 'M_r_w_squared'
                       or function_name == 'M_r_w_stddev'
                       or function_name == 'M_k_w_function'
                       or function_name == 'cluster_hybridization_F_k_w'
                       or function_name == 'cluster_greens_function_G_k_w_approximant'
                       or function_name == 'G_k_w_measured'
                       or function_name == 'cluster_greens_function_G_k_w_measured'
                       or function_name == 'G0_k_PCM_w_cluster_excluded'
                       or function_name == 'G0_k_PCM_w'
                       or function_name == 'S_k_PCM_w'
                       #or function_name == 'S_r_PCM_w'
                       or function_name == 'G_k_PCM_w'
                       or function_name == 'S_k_DCA_w'
                       or function_name == 'G_k_DCA_w'
                       or function_name == 'G0_r_PCM_w_cluster_excluded'
                       or function_name == 'G0_k_PCM_w_cluster_excluded'
                       or function_name == 'M_k_PCM_w'
                       or function_name == 'M_r_PCM_w'
                       or function_name == 'M_k_DCA_w'
                       or function_name == 'Sigma-1st-order'
                       or function_name == 'Sigma-2nd-order'
                       or function_name == 'perturbation-Sigma'
                       or function_name == 'Self-Energy-cluster'
                       or function_name == 'Self-energy-lattice'):
                        plot_wn_function(domains, function_name, datamatrix)

                    if(function_name == 'S_r_PCM'
                       or function_name == 'alpha_r_DCA'):
                        plot_Sigma_r_wn_function(domains, function_name, datamatrix)

                    if(function_name == 'cluster_greens_function_G_k_t'
                       or function_name == 'cluster_greens_function_G_r_t'
                       or function_name == 'free_cluster_greens_function_G0_k_t'
                       or function_name == 'free_cluster_greens_function_G0_r_t'
                       or function_name == 'cluster_excluded_greens_function_G0_k_t'
                       or function_name == 'cluster_excluded_greens_function_G0_r_t'
                       or function_name == 'cluster_hybridization_F_k_t'
                       or function_name == 'cluster_hybridization_F_r_t'
                       or function_name == 'G0_r_PCM_t_cluster_excluded'
                       or function_name == 'G0_k_PCM_t'
                       or function_name == 'M_r_t'
                       or function_name == 'K_r_t'
                       or function_name == 'G_r_t_measured'
                       or function_name == 'G_r_t_stddev'):
                        plot_tau_function(domains, function_name, datamatrix)

                    if(function_name == '<k>'):
                        plot_feynamd_expansion_order(domains, function_name, datamatrix)

                    if(function_name == 'sign'
                       or function_name == 'density'
                       or function_name == 'chemical-potential'
                       or function_name == 'expansion_order'
                       or function_name == 'updates'
                       or function_name == 'successfull_updates'
                       or function_name == 'unsuccessfull_updates'):
                        plot_DCA_iteration_function(domains, function_name, datamatrix)

                    if(function_name == 'L2_Sigma_difference'):
                        plot_DCA_iteration_function_log(domains, function_name, datamatrix)

                    if(function_name == 'orbital_occupancies'):
                        plot_densities(domains, function_name, datamatrix)

                    if(function_name == 'spectral_density' 
                       or function_name == 'free_spectral_density'
                       or function_name == 'Im_G_w' 
                       or function_name == 'Im_G0_w'):
                        plot_spectrum(domains, function_name, datamatrix)

                    if(function_name == 'cluster_greens_function_G_k_w_real_axis' 
                       or function_name == 'cluster_greens_function_G0_k_w_real_axis'
                       or function_name == 'cluster_greens_function_G_r_w_real_axis' 
                       or function_name == 'cluster_greens_function_G0_r_w_real_axis'
                       or function_name == 'Self_energy_real_axis'):
                        plot_real_axis(domains, function_name, datamatrix)

                    if(function_name == 'a_x' 
                       or function_name == 'alpha'):
                        plot_fine_real_axis(domains, function_name, datamatrix)

                    if(function_name == 'a_wn'):
                        plot_imaginary_axis(file, line, domains, datamatrix)

                    if(function_name == 'f-approx'):
                        plot_imaginary_axis_new(file, line, domains, datamatrix)

                    if(function_name == 'band-structure'):
                        plot_band_structure(domains, function_name, datamatrix)

                    if(function_name == 'Sigma-band-structure'):
                        plot_S_band_structure(file, line, function_name, datamatrix)

                    if(function_name == 'G_k_dmn_t_w'):
                        plot_k_cut_spectrum(domains, function_name, datamatrix)

#                     if(function_name == 'S_k_PCM_w_cut'
#                        or function_name == 'alpha_k_cut'
#                        or function_name == 'Gamma_k_PCM_w_cut'
#                        or function_name == 'G0_k_PCM_w_cluster_excluded_cut'):
#                          plot_cut(domains, function_name, datamatrix)

                    if(function_name == 'G4_k_k_w_w'):
                        plot_G4_k_k_w_w(function_name, datamatrix)

                    if(function_name == 'Sigma_zero_moment'
#                        or function_name == 'A_k'
#                        or function_name == 'n_k'
                       ) :
                        plot_Sigma_zero_moment(domains,function_name, datamatrix)

#                     if(function_name == 'A_k'
#                        or function_name == 'n_k'):
#                          plot_f_k_function(domains,function_name, datamatrix)

                    if( function_name == 'Fermi_surface'
                        or function_name == 'DOS_surface'):
                        plot_f_k_HOST_function(domains,function_name, datamatrix, 'imag')

                    if(function_name == 'Z_k'):
                        plot_f_k_HOST_function(domains,function_name, datamatrix, 'real')

                    if(function_name == 'n_k'
                       or function_name == 'grad_n_k'):
                        plot_n_k_HOST_function(domains,function_name, datamatrix)

                    if(function_name == 'leading_eigenvector'
                       or function_name == 'Gamma_diagonal'):
                        plot_leading_eigenvector(domains,function_name, datamatrix)

                    if(function_name == 'leading_eigenvectors'
                       or function_name == 'leading_U_K'
                       or function_name == 'leading_Vt_K'
                       or function_name == 'leading_U_k_interpolated'
                       or function_name == 'leading_Vt_k_interpolated'
                       or function_name == 'leading_VR_vectors_dca'
                       or function_name == 'leading_VR_vectors_pcm'
                       or function_name == 'leading_VR_vectors_dca_inverse'):
                        #plot_leading_eigenvectors(domains,function_name, datamatrix)
                        fit_leading_eigenvectors(domains,function_name, datamatrix)

                    if(function_name == 'phi_r'):
                        plot_phi_r(domains,function_name, datamatrix)
                    
                    if(function_name == 'leading_alpha'
                        or function_name == 'leading_eigenvalues'):
                       plot_leading_eigenvalues(domains,function_name, datamatrix)
                                               
                    if(function_name == 'wave_function'):
                        plot_wave_function(domains,function_name, datamatrix)

                    if(function_name == 'chi_0_function'
                       or function_name == 'full_chi_0_function'
                       or function_name == 'chi_0'
                       or function_name == 'G4_0'):
                        plot_chi_0(domains,function_name, datamatrix)
 
                    if(function_name == 'local_chi_ph_magnetic'
                       or function_name == 'local_chi_ph_charge'
                       or function_name == 'local_chi_pp'
                       or function_name == 'local_chi_0_ph'
                       or function_name == 'ph-bubble'
                       or function_name == 'pp-bubble'
                       or function_name == 'U_eff'):
                        plot_local_susceptibilities(domains,function_name, datamatrix)

                    if(function_name == 'G0_k_w_times_G0_k_plus_q_w_plus_W'
                       or function_name == 'G0_k_w_times_G0_q_min_k_W_min_w'):
                        plot_local_chi_quantities(domains,function_name, datamatrix)

                    if(function_name == 'phi_dependent_G_k_PCM_w'
                       or function_name == 'phi_dependent_G0_k_PCM_w'
                       or function_name == 'phi_dependent_G0_k_PCM_t'
                       or function_name == 'phi_dependent_S_k_PCM_w'
                       or function_name == 'phi_dependent_M_r_PCM_w'
                       or function_name == 'phi_dependent_M_k_PCM_w'
                       or function_name == 'phi_dependent_G0_k_PCM_w_cluster_excluded'):
                        plot_phi_wn_quantities(domains,function_name, datamatrix)

                    if(function_name == 'phi_dependent_G0_r_PCM_t'
                       or function_name == 'phi_dependent_G0_r_PCM_t_cluster_excluded'):
                        plot_phi_t_quantities(domains,function_name, datamatrix)

                    if(function_name=="Sigma-r-lattice"):
                        plot_Sigma_r_decay(domains, function_name, datamatrix)

                    if(function_name=="Sigma-k-lattice"
                       or function_name=="Greens-k-lattice"):
                        plot_Sigma_k_HOST(domains, function_name, datamatrix)

                    if(function_name=="numerical-error"):
                        plot_numerical_error(domains, function_name, datamatrix)

                    datamatrix = [[]]
                else:
                    if os.path.exists('./'+filename+'/' + function_name):
                        shutil.rmtree('./'+filename+'/' + function_name)

