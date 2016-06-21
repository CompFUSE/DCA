import commands
import shutil
import os
import sys
import time

if os.path.exists('../bin'):
    shutil.rmtree('../bin')
os.makedirs('../bin')

if not os.path.exists('../../version.txt'):
    print "\n\t--> print version.txt"
    cmd = "cd ../..; git log -n 1 > version.txt"
    os.system(cmd)

version_text=""

for line in open("../../version.txt", "r"):
    print line
    line2 = line.strip()
    version_text = version_text + "\n\t\"" + line2 + "\\n\""

print version_text

###############################
###                         ###
###  write include_files.h  ###
###                         ###
###############################

file = open("include_files.tmpl", "r")
text = file.read()
file.close()

#
#  concurrency
#

concurrency = raw_input('Please enter execution-mode : \n\t(1) SERIAL \n\t(2) MPI \n\t(3) MPI_FT \n\t(4) OPENMPI_FT \n\n')

if(concurrency == '1'):    
    text = text.replace("//#include \"include_concurrency_serial.h\""     , "#include \"include_concurrency_serial.h\"")

if(concurrency == '2'):    
    text = text.replace("//#include \"include_concurrency_mpi.h\""     , "#include \"include_concurrency_mpi.h\"")  

if(concurrency == '3'):    
    text = text.replace("//#include \"include_concurrency_mpi_ft.h\""     , "#include \"include_concurrency_mpi_ft.h\"")  

if(concurrency == '4'):    
    text = text.replace("//#include \"include_concurrency_ompi_ft.h\""     , "#include \"include_concurrency_ompi_ft.h\"")  

file = open("../bin/include_files.h", "w")
file.write(text) 
file.close() 


#####################################
###                               ###
###  write compiler_directives.h  ###
###                               ###
#####################################

file = open("compiler_directives.tmpl", "r")
text = file.read()
file.close()
#
#  profiler
#

profiler = raw_input('Please enter profiling-mod : \n\t(1) NONE \n\t(2) COUNTING \n\t(3) PAPI \n\n')

if(profiler == '1'):    
    text = text.replace("//#define NO_PROFILING"       , "#define NO_PROFILING")

if(profiler == '2'):    
    text = text.replace("//#define COUNTING_PROFILING" , "#define COUNTING_PROFILING")  

if(profiler == '3'):    
    text = text.replace("//#define DCA_PAPI_PROFILING", "#define DCA_PAPI_PROFILING")  

#
#  Random number generator
#

rng = raw_input('Please enter random-number generator : \n\t(1) Numerical recipes \n\t(2) SPRNG \n\n')

if(rng == '1'):    
    text = text.replace("//#define RNG_NUMERICAL_RECIPES", "#define RNG_NUMERICAL_RECIPES")

if(rng == '2'):    
    text = text.replace("//#define RNG_SPRNG"             , "#define RNG_SPRNG")  

#
#  assert
#

rng = raw_input('do you wish to assert? (y/n) \n\n\n')

if(rng == 'y'):    
    text = text.replace("#define NDEBUG    // --> comment out if you want to assert !! ", "//#define NDEBUG ")

if(rng == 'n'):    
    text = text.replace("#define NDEBUG    // --> comment out if you want to assert !! ", "#define NDEBUG ")

#
#  measurements precision
#

precision = raw_input('do you wish to measure in single precision? (y/n) \n\n\n')

if(precision == 'y'):    
    text = text.replace("//#define SINGLE_PRECISION_MEASUREMENTS", "#define SINGLE_PRECISION_MEASUREMENTS")

#
#  CUDA
#

gpu = raw_input('do you want GPU support? (y/n) \n\n\n')

if(gpu == 'y'):    
    text = text.replace("//#define USE_GPU", "#define USE_GPU")

#
#  std deviation + error
#

stddev = raw_input('do you want l1/l2/linf norm of errors? (y/n) \n\n\n')

if(stddev == 'y'):    
    text = text.replace("//#define MEASURE_ERROR_BARS", "#define MEASURE_ERROR_BARS")


file = open("../bin/compiler_directives.h", "w")
file.write(text) 
file.close() 

###############################
###                         ###
###  write main.cpp         ###
###                         ###
###############################

file = open("main.tmpl", "r")
text = file.read()
file.close()

if(concurrency == '1'):    
    text = text.replace("MPI_LIB_TYPE"     , "SERIAL_LIBRARY")

if(concurrency == '2'):    
    text = text.replace("MPI_LIB_TYPE"     , "MPI_LIBRARY")

if(concurrency == '3'):    
    text = text.replace("MPI_LIB_TYPE"     , "MPI_FT_LIBRARY")

if(concurrency == '4'):    
    text = text.replace("MPI_LIB_TYPE"     , "OPENMPI_FT_LIBRARY")

MC_SOLVER = raw_input('Please enter Monte-Carlo algorithm : \n\t(1) continuous-time auxilary field \n\t(2) single-site Hybridization \n\t(3) Hybridization \n\t(4) projected-cluster method \n\n')

if(MC_SOLVER == '1'):    
    text = text.replace("MC_ALG_TYPE"    , "CT_AUX")
    text = text.replace("MC_ACC_TYPE"    , "MATSUBARA")

if(MC_SOLVER == '2'):    
    text = text.replace("MC_ALG_TYPE"    , "HYBRIDIZATION")

    MC_ACCUMULATOR = raw_input('Please enter accumulator method : \n\t(1) matsubara \n\t(2) legendre \n\n')
    if(MC_ACCUMULATOR == '2'):
        text = text.replace("MC_ACC_TYPE"    , "LEGENDRE")
    else:
        text = text.replace("MC_ACC_TYPE"    , "MATSUBARA")

if(MC_SOLVER == '3'):    
    text = text.replace("MC_ALG_TYPE"    , "HYBRIDIZATION_FULL")
    text = text.replace("MC_ACC_TYPE"    , "MATSUBARA")

if(MC_SOLVER == '4'):    
    text = text.replace("MC_ALG_TYPE"    , "PCM")
    text = text.replace("MC_ACC_TYPE"    , "MATSUBARA")
    
    SELF_ENERGY_METHOD = raw_input('Please enter Self-energy treatment : \n\t(1) coarsegrain \n\t(2) interpolation \n\n')

    file = open("../bin/include_files.h", "r")
    text2 = file.read()
    file.close()

    text2 = text2.replace("#define USE_REDUCED_VERTEX_FUNCTION", "//#define USE_REDUCED_VERTEX_FUNCTION")

    if(SELF_ENERGY_METHOD == '1'):
        text2 = text2.replace("//#define COARSEGRAIN_SELF_ENERGY","#define COARSEGRAIN_SELF_ENERGY")

    file = open("../bin/include_files.h", "w")
    file.write(text2) 
    file.close() 

file = open("../bin/include_files.h", "r")
text3 = file.read()
file.close()

if  (MC_SOLVER == '1'):
    text3 = text3.replace("//#include \"include_CT_AUX.h\"","#include \"include_CT_AUX.h\"")
elif(MC_SOLVER == '2'):
    text3 = text3.replace("//#include \"include_single_site_Hybridization.h\"","#include \"include_single_site_Hybridization.h\"")
elif(MC_SOLVER == '3'):
    text3 = text3.replace("//#include \"include_hybridization_full.h\"","#include \"include_hybridization_full.h\"")
elif(MC_SOLVER == '4'):
    text3 = text3.replace("//#include \"include_CT_AUX.h\"","#include \"include_CT_AUX.h\"")
    text3 = text3.replace("//#include \"include_PCM.h\"","#include \"include_PCM.h\"")

file = open("../bin/include_files.h", "w")
file.write(text3)
file.close()

text = text.replace("\"DEFAULT\"", version_text)

file = open("../bin/main.cpp", "w")
file.write(text) 
file.close()

#
#  write analysis.cpp
#

file = open("analysis.tmpl", "r")
text = file.read()
file.close()

if(MC_SOLVER == '1'):    
    ANALYSIS_METHOD = raw_input('do you want an interpolation-schem? (y/n)\n\n')

    if(ANALYSIS_METHOD == 'y'):
        text = text.replace("MC_ALG_TYPE"    , "ANALYSIS_INTERPOLATION")
    else:
        text = text.replace("MC_ALG_TYPE"    , "CT_AUX")

if(MC_SOLVER == '2'):    
    text = text.replace("MC_ALG_TYPE"    , "HYBRIDIZATION")

if(MC_SOLVER == '3'):    
    text = text.replace("MC_ALG_TYPE"    , "HYBRIDIZATION_FULL")

if(MC_SOLVER == '4'):    
    text = text.replace("MC_ALG_TYPE"    , "PCM")

text = text.replace("\"DEFAULT\"", version_text)

file = open("../bin/analysis.cpp", "w")
file.write(text) 
file.close()

####################################
###                              ###
###   write type_definitions.h   ###
###                              ###
####################################

cmd = "cp ./type_definitions.h ../bin/type_definitions.h"
os.system(cmd)

####################################
###                              ###
###   write MakeFile             ###
###                              ###
####################################

machine = raw_input('Please enter machine : \n\t(1) MAC \n\t(2) MAC SNOW-LEOPARD \n\t(3) JAGUAR \n\t(4) LINUX \n\t(5) ROSA \n\t(6) EIGER \n\t(7) TODI \n\n')

if(machine == '1' or machine == '2'):

    file = open("Makefile.tmpl", "r")
    text = file.read()
    file.close()

    if(concurrency == '2'):
        if(machine == '1'):
            text = text.replace("COMPILER_C++"   , "../../../mpich2-1.2.1p1_gcc45/bin/mpic++")
        if(machine == '2'):
            text = text.replace("COMPILER_C++"   , "mpicxx")
    else:
        text = text.replace("COMPILER_C++"   , "g++")

    text = text.replace("COMPILER_C"     , "gcc")
    text = text.replace("FLAGOPTIONS"    , "-O3 -Wall -Werror -funroll-loops -ffast-math -msse4 -finline-functions")
    text = text.replace("LIBRARIES"      , "-llapack -lblas -lm ../../levmar-2.5/liblevmar.a ../../nfft-3.1.3/build/lib/libnfft3.a ../../fftw-3.2.2/build/lib/libfftw3.a " )
    text = text.replace("MAKE_FFT_THREADS", " -j 4 ")
    text = text.replace("OPTIONS_FFTW"   , "--prefix=`pwd`/build --disable-fortran")
    text = text.replace("OPTIONS_NFFT"   , "--prefix=`pwd`/build --with-fftw3=`pwd`/../fftw-3.2.2/build --without-apple-gcc-arch")
    text = text.replace("-SRC_DIR_SPRNG" , "")
    text = text.replace("-SRC_FFTW"      , "-I../../fftw-3.2.2/build/include \\")
    text = text.replace("-SRC_DIR_CUDA"  , "")

    file = open("../bin/Makefile", "w")
    file.write(text) 
    file.close() 

if(machine == '3'):

    file = open("Makefile.tmpl", "r")
    text = file.read()
    file.close()

    text = text.replace("COMPILER_C++"   , "CC")
    text = text.replace("COMPILER_C"     , "cc")
    text = text.replace("FLAGOPTIONS"    , "-O3 -Wall -Werror -funroll-loops -ffast-math -mtune=opteron -march=opteron -msse3 -finline-functions")
    text = text.replace("LIBRARIES"      , "../../levmar-2.5/liblevmar.a ../../nfft-3.1.3/build/lib/libnfft3.a ../../fftw-3.2.2/build/lib/libfftw3.a " )
    text = text.replace("MAKE_FFT_THREADS", " -j 4 ")
    text = text.replace("OPTIONS_FFTW"   , "--prefix=`pwd`/build --disable-fortran")
    text = text.replace("OPTIONS_NFFT"   , "--prefix=`pwd`/build --with-fftw3=`pwd`/../fftw-3.2.2/build" )
    text = text.replace("-SRC_DIR_SPRNG" , "")
    text = text.replace("-SRC_FFTW"      , "-I../../fftw-3.2.2/build/include \\")
    text = text.replace("-SRC_DIR_CUDA"  , "")

    file = open("../bin/Makefile", "w")
    file.write(text) 
    file.close() 


if(machine == '4'):

    file = open("Makefile.tmpl", "r")
    text = file.read()
    file.close()

    if(concurrency == '2'):
        text = text.replace("COMPILER_C++"   , "../../../mpich2-1.2.1p1/bin/mpic++")
    else:
        text = text.replace("COMPILER_C++"   , "g++")

    text = text.replace("COMPILER_C"     , "gcc")
    text = text.replace("FLAGOPTIONS"    , "-O3 -Wall -Werror -funroll-loops -ffast-math -msse4 -finline-functions")
    text = text.replace("LIBRARIES"      , "-llapack -lblas -lm ../../levmar-2.5/liblevmar.a ../../nfft-3.1.3/build/lib/libnfft3.a ../../fftw-3.2.2/build/lib/libfftw3.a " )
    text = text.replace("MAKE_FFT_THREADS", " -j 4 ")
    text = text.replace("OPTIONS_FFTW"   , "--prefix=`pwd`/build                                        --without-gcc-arch --with-pic")
    text = text.replace("OPTIONS_NFFT"   , "--prefix=`pwd`/build --with-fftw3=`pwd`/../fftw-3.2.2/build --without-gcc-arch --with-pic" )
    text = text.replace("-SRC_DIR_SPRNG" , "")
    text = text.replace("-SRC_FFTW"      , "-I../../fftw-3.2.2/build/include \\")
    text = text.replace("-SRC_DIR_CUDA"  , "")

    file = open("../bin/Makefile", "w")
    file.write(text) 
    file.close() 

if(machine == '5'):

        file = open("Makefile.tmpl", "r")
        text = file.read()
        file.close()

        text = text.replace("COMPILER_C++"   , "CC")
        text = text.replace("COMPILER_C"     , "cc")
        text = text.replace("FLAGOPTIONS"    , "-O3 -Wall -Werror -funroll-loops -ffast-math -mtune=opteron -march=opteron -msse3 -finline-functions -Wl,-ydgemm_")
        #text = text.replace("LIBRARIES"      , "../../levmar-2.5/liblevmar.a ../../nfft-3.1.3/build/lib/libnfft3.a ../../fftw-3.2.2/build/lib/libfftw3.a ../../../sprng4/lib/libsprng.a" )
        text = text.replace("LIBRARIES"      , "../../levmar-2.5/liblevmar.a ../../nfft-3.1.3/build/lib/libnfft3.a /opt/fftw/3.2.2.1/lib/libfftw3.a ../../../sprng4/lib/libsprng.a" )
        text = text.replace("MAKE_FFT_THREADS", " -j 4 ")
        text = text.replace("OPTIONS_FFTW"   , "--prefix=`pwd`/build --disable-fortran --with-pic")
        text = text.replace("OPTIONS_NFFT"   , "--prefix=`pwd`/build --with-fftw3=/opt/fftw/3.2.2.1" )
        text = text.replace("-SRC_DIR_SPRNG" , "-I../../../sprng4/include \\\n\t-I../../../sprng4/lib ")
        text = text.replace("-SRC_FFTW"      , "-I/opt/fftw/3.2.2.1/include \\")
        text = text.replace("-SRC_DIR_CUDA"  , "")

        file = open("../bin/Makefile", "w")
        file.write(text)
        file.close()

if(machine == '6'):

    file = open("Makefile.tmpl", "r")
    text = file.read()
    file.close()

    if(concurrency == '2' or concurrency == '4'):
        text = text.replace("COMPILER_C++"   , "mpicxx")
    elif(concurrency == '3'):
        text = text.replace("COMPILER_C++"   , "ftmpicxx")
    else:
        text = text.replace("COMPILER_C++"   , "g++")

    text = text.replace("COMPILER_C"     , "gcc")
    text = text.replace("FLAGOPTIONS"    , "-O3 -Wall -Werror -funroll-loops -ffast-math -finline-functions")
    if(gpu == 'y'):
        text = text.replace("LIBRARIES"      , "../../levmar-2.5/liblevmar.a ../../nfft-3.1.3/build/lib/libnfft3.a  ../../fftw-3.2.2/build/lib/libfftw3.a  -L/apps/eiger/Intel-FOR-11.1/mkl/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -L/apps/eiger/Cuda-4.0/cuda/lib64/ -lcublas")
    else:
        text = text.replace("LIBRARIES"      , "../../levmar-2.5/liblevmar.a ../../nfft-3.1.3/build/lib/libnfft3.a  ../../fftw-3.2.2/build/lib/libfftw3.a  -L/apps/eiger/Intel-FOR-11.1/mkl/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread")
    text = text.replace("MAKE_FFT_THREADS", " -j 6 ")
    text = text.replace("OPTIONS_FFTW"   , "--prefix=`pwd`/build                                        --without-gcc-arch --with-pic")
    text = text.replace("OPTIONS_NFFT"   , "--prefix=`pwd`/build --with-fftw3=`pwd`/../fftw-3.2.2/build --without-gcc-arch --with-pic" )
    text = text.replace("-SRC_FFTW"      , "-I../../fftw-3.2.2/build/include \\")
    text = text.replace("-SRC_DIR_SPRNG" , "-I$(SRC_DIR)")
    text = text.replace("-SRC_DIR_CUDA"  , "-I/apps/eiger/Cuda-4.0/cuda/include/ \\")

    file = open("../bin/Makefile", "w")
    file.write(text) 
    file.close() 

if(machine == '7'):

        file = open("Makefile.tmpl", "r")
        text = file.read()
        file.close()

        text = text.replace("COMPILER_C++"   , "CC")
        text = text.replace("COMPILER_C"     , "cc")
        text = text.replace("FLAGOPTIONS"    , "-O3 -Wall -Werror -funroll-loops -ffast-math -mtune=opteron -march=opteron -msse3 -finline-functions -Wl,-ydgemm_")
        text = text.replace("LIBRARIES"      , "../../levmar-2.5/liblevmar.a ../../nfft-3.1.3/build/lib/libnfft3.a  /opt/fftw/3.3.0.0/interlagos/lib/libfftw3.a" )
        text = text.replace("MAKE_FFT_THREADS", " -j 4 ")
        text = text.replace("OPTIONS_FFTW"   , "--prefix=`pwd`/build --disable-fortran --with-pic")
        text = text.replace("OPTIONS_NFFT"   , "--prefix=`pwd`/build --with-fftw3=/opt/fftw/3.3.0.0/interlagos" )
        text = text.replace("-SRC_DIR_SPRNG" , " \ ")
        text = text.replace("-SRC_FFTW"      , "-I/opt/fftw/3.3.0.0/interlagos/include \\")
        text = text.replace("-SRC_DIR_CUDA"  , "")

        file = open("../bin/Makefile", "w")
        file.write(text)
        file.close()


print "\n\t starting to build :: \n\n"

out = raw_input('rebuild FFTW & NFFT ? (n/y) :')

if(out == 'y'):
    if(machine != '5'):
        cmd = "cd ../bin; make FFTW"
        os.system(cmd)

    cmd = "cd ../bin; make NFFT"
    os.system(cmd)

    cmd = "cd ../bin; make LEVMAR"
    os.system(cmd)

if(machine == '3' or machine == '5'):
    cmd = "module load fftw"
    os.system(cmd)

    cmd = "module load xt-papi"
    os.system(cmd)

    cmd = "module load PrgEnv-gnu"
    os.system(cmd)

if(machine == '6'):
    cmd = "module load mkl"
    os.system(cmd)

    if(gpu == 'y'):
        cmd = "module load cuda"
        os.system(cmd)

cmd = "cd ../bin; make"
os.system(cmd)

    


