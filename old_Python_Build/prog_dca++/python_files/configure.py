import commands
import shutil
import os
import sys
import time

class configure:
    """A class that collects the data for the configuration"""

    compiler_dic = {"Mac"      : "$(MPICH__DIR)/build/bin/mpicxx", #"../../../mpich2-1.2.1p1_clang33/bin/mpic++",
                    "Linux"    : " mpic++",
                    "Jaguarpf" : " CC",
                    "Rosa"     : " CC",
                    "Todi"     : " CC",
                    "Titan"    : " CC",
                    "Daint"    : " CC"}

    flags_dic = {"Mac"      : "-g3 -O3    -Wall -Werror -Wno-unused-local-typedefs -Wno-sign-compare #-funroll-loops -finline-functions",
                 "Linux"    : " -lfftw -lblas -llapack",
                 "Jaguarpf" : "-g3 -O3    -Wall -Werror #-funroll-loops -mtune=opteron -march=opteron -finline-functions -Wl,-ydgemm_",
                 "Rosa"     : "-g3 -O3 -Wall -Werror -Wno-unused-local-typedefs -Wno-sign-compare #-funroll-loops -finline-functions -Wl,-ydgemm_ #-mtune=opteron -march=opteron -msse3",
                 "Todi"     : "-g3 -O3 -Wall -Werror -Wno-unused-local-typedefs -Wno-sign-compare #-funroll-loops -finline-functions -Wl,-ydgemm_ #-mtune=opteron -march=opteron -msse3",
                 "Titan"    : "-g3 -O3 -Wall -Werror -Wno-unused-local-typedefs -Wno-sign-compare #-funroll-loops -finline-functions -Wl,-ydgemm_ #-mtune=opteron -march=opteron -msse3",
                 "Daint"    : "-g3 -O3 -Wall -Werror -Wno-unused-local-typedefs -Wno-sign-compare #-funroll-loops -finline-functions -Wl,-ydgemm_ #-mtune=opteron -march=opteron -msse3"}

    include_dic = {"Mac"      : "-I../../fftw-3.3.3/build/include \\\n\t-I$(NFFT___DIR)/build/include \\\n\t-I$(SPGLIB_DIR)/build/include/spglib ",
                   "Linux"    : "-I../../fftw-3.3.3/build/include \\\n\t-I$(NFFT___DIR)/build/include \\\n\t-I$(SPGLIB_DIR)/build/include/spglib ",
                   "Jaguarpf" : "-I/opt/fftw/3.3.0.0/interlagos/include \\\n\t-I$(NFFT___DIR)/build/include \\\n\t-I$(SPGLIB_DIR)/build/include/spglib ",
                   "Rosa"     : "-I/opt/fftw/3.3.0.0/interlagos/include \\\n\t-I$(NFFT___DIR)/build/include \\\n\t-I$(SPGLIB_DIR)/build/include/spglib ",
                   "Todi"     : "-I/opt/fftw/3.3.0.0/interlagos/include \\\n\t-I$(NFFT___DIR)/build/include \\\n\t-I$(SPGLIB_DIR)/build/include/spglib ",
                   "Titan"    : "-I/opt/fftw/3.3.0.0/interlagos/include \\\n\t-I$(NFFT___DIR)/build/include \\\n\t-I$(SPGLIB_DIR)/build/include/spglib ",
                   "Daint"    : "-I/opt/fftw/3.3.0.1/x86_64/include \\\n\t-I$(NFFT___DIR)/build/include \\\n\t-I$(SPGLIB_DIR)/build/include/spglib "}

    mac_libs = "-llapack -lblas -lm \\\n\
	$(HDF5___DIR)/build/lib/libhdf5.8.dylib \\\n\
	$(HDF5___DIR)/build/lib/libhdf5.a \\\n\
	$(HDF5___DIR)/build/lib/libhdf5_cpp.8.dylib \\\n\
	$(HDF5___DIR)/build/lib/libhdf5_cpp.a \\\n\
	$(NFFT___DIR)/build/lib/libnfft3.a \\\n\
	$(FFTW___DIR)/build/lib/libfftw3.a \\\n\
	$(SPGLIB_DIR)/build/lib/libsymspg.a\n"
    
    library_dic = {"Mac"      : mac_libs,#" -llapack -lblas -lm $(SPGLIB_DIR)/build/lib/libsymspg.a $(NFFT___DIR)/build/lib/libnfft3.a ../../fftw-3.2.2/build/lib/libfftw3.a",
                   "Linux"    : mac_libs,#" -llapack -lblas -lm $(SPGLIB_DIR)/build/lib/libsymspg.a $(NFFT___DIR)/build/lib/libnfft3.a ../../fftw-3.2.2/build/lib/libfftw3.a",
                   "Jaguarpf" : "  $(SPGLIB_DIR)/build/lib/libsymspg.a $(NFFT___DIR)/build/lib/libnfft3.a",
                   "Rosa"     : "  $(SPGLIB_DIR)/build/lib/libsymspg.a $(NFFT___DIR)/build/lib/libnfft3.a",
                   "Todi"     : "  $(SPGLIB_DIR)/build/lib/libsymspg.a $(NFFT___DIR)/build/lib/libnfft3.a",
                   "Titan"    : "  $(SPGLIB_DIR)/build/lib/libsymspg.a $(NFFT___DIR)/build/lib/libnfft3.a",
                   "Daint"    : "  $(SPGLIB_DIR)/build/lib/libsymspg.a $(NFFT___DIR)/build/lib/libnfft3.a #-mkl -dynamic",
                   "Eiger"    : "  $(SPGLIB_DIR)/build/lib/libsymspg.a $(NFFT___DIR)/build/lib/libnfft3.a ../../fftw-3.2.2/build/lib/libfftw3.a"}

    GPU_flags_dic = {"Todi"     : "-arch=sm_35",
                     "Titan"    : "-arch=sm_35",
                     "Daint"    : "-arch=sm_35"}

#     GPU_magma_root_dic = {"Todi"     : "/project/s299/rasolca/todi/magmaeigen/",
#                           "Titan"    : "/ccs/home/staarp/magma/magmaeigen/",
#                           "Daint"    : "/project/s299/rasolca/todi/magmaeigen/"}

    GPU_library_dic = {"Todi"     : "  -lmagma -lcuda -lcudart -lcublas -lcblas",
                       "Titan"    : "  -lmagma -lcuda -lcudart -lcublas -lcblas",
                       "Daint"    : "  -lmagma -lcuda -lcudart -lcublas -lcblas"}

    GPU_CC = {"Mac"      : "/Developer/NVIDIA/CUDA-5.0/bin/nvcc",
              "Todi"     : "nvcc",
              "Titan"    : "nvcc",
              "Daint"    : "nvcc"}

    FFTW_THREADS_dic = {"Mac"      : "-j 8",
                        "Linux"    : "-j 1",
                        "Jaguarpf" : "-j 8",
                        "Rosa"     : "-j 8",
                        "Todi"     : "-j 8",
                        "Titan"    : "-j 8",
                        "Daint"    : "-j 8",
                        "Eiger"    : "-j 8"}

    FFTW_dic = {"Mac"      : "--prefix=`pwd`/build --disable-fortran",
                "Linux"    : "--prefix=`pwd`/build --without-gcc-arch --with-pic",
                "Jaguarpf" : "--prefix=`pwd`/build --disable-fortran",
                "Rosa"     : "--prefix=`pwd`/build --disable-fortran --with-pic",
                "Todi"     : "--prefix=`pwd`/build --disable-fortran --with-pic",
                "Titan"    : "--prefix=`pwd`/build --disable-fortran --with-pic",
                "Daint"    : "--prefix=`pwd`/build --disable-fortran --with-pic"}

    NFFT_dic = {"Mac"      : "--prefix=`pwd`/build --with-fftw3=`pwd`/../fftw-3.2.2/build --without-apple-gcc-arch",
                "Linux"    : "--prefix=`pwd`/build --with-fftw3=`pwd`/../fftw-3.2.2/build --without-apple-gcc-arch --with-pic",
                "Jaguarpf" : "--prefix=`pwd`/build --with-fftw3=`pwd`/../fftw-3.2.2/build" ,
                "Rosa"     : "--prefix=`pwd`/build --with-fftw3=/opt/fftw/3.3.0.0/interlagos",
                "Todi"     : "--prefix=`pwd`/build --with-fftw3=/opt/fftw/3.3.0.0/interlagos",
                "Titan"    : "--prefix=`pwd`/build --with-fftw3=/opt/fftw/3.3.0.0/interlagos",
                "Daint"    : "--prefix=`pwd`/build --with-fftw3=/opt/fftw/3.3.0.1/x86_64"}

    SPGLIB_dic = {"Mac"      : "--prefix=`pwd`/build --enable-cxx CC=$(CC) CXX=$(CPPC)",
                  "Rosa"     : "--prefix=`pwd`/build",
                  "Todi"     : "--prefix=`pwd`/build",
                  "Titan"     : "--prefix=`pwd`/build",
                  "Daint"     : "--prefix=`pwd`/build"
                  }

    HDF5_dic = {"Mac"      : "--prefix=`pwd`/build --enable-cxx CC=$(GCC) CXX=$(GPP)",
                "Rosa"     : "--prefix=`pwd`/build",
                "Todi"     : "--prefix=`pwd`/build",
                "Titan"    : "--prefix=`pwd`/build",
                "Daint"    : "--prefix=`pwd`/build"
                }

    MPICH_dic = {"Mac"      : "--prefix=`pwd`/build --enable-cxx CC=$(GCC) CXX=$(GPP) --disable-f77 --disable-fc",
                 "Rosa"     : "--prefix=`pwd`/build",
                 "Todi"     : "--prefix=`pwd`/build",
                 "Titan"    : "--prefix=`pwd`/build",
                 "Daint"    : "--prefix=`pwd`/build"
                 }

    def __init__(self):

        self.version = self.read_version()
        
        self.QMC_type       = "CT-AUX"
        self.concurrency    = "MPI"
        self.use_pthreads   =  False
        self.use_sse        =  False
        self.profiling      = "no-profiling"
        self.auto_tuning    = "False"
        self.rng_type       = "NR"
        
        self.do_BIT                      = False
        self.do_assert                   = False
        self.compute_std_dev             = False
        #self.measure_in_double_precision = False
        self.single_precision_measuring      = True
        self.single_precision_coarsegraining = False
        self.use_full_vertex             = False
        
        self.GPU_support            = False
        self.use_pinned_host_memory = False
        
        self.machine_name         = "Mac"

    def __str__(self):
        str_json = "{\n"
        
        str_json = str_json + "version-stamp : " + self.version + ",\n\n"

        str_json = str_json + "QMC_type     : " + self.QMC_type + ",\n"
        str_json = str_json + "concurrency  : " + self.concurrency + ",\n"
        str_json = str_json + "use_pthreads : " + str(self.use_pthreads) + ",\n"
        str_json = str_json + "profiling    : " + self.profiling + ",\n"
        str_json = str_json + "auto_tuning  : " + self.auto_tuning + ",\n"
        str_json = str_json + "rng_type     : " + self.rng_type + ",\n"

        str_json = str_json + "do_BIT                          : " + str(self.do_BIT)              + ",\n"
        str_json = str_json + "do_assert                       : " + str(self.do_assert)           + ",\n"
        str_json = str_json + "compute_std_dev                 : " + str(self.compute_std_dev)     + ",\n"
        #str_json = str_json + "measure_in_double_precision : " + str(self.measure_in_double_precision) + ",\n"
        str_json = str_json + "single_precision_measuring      : " + str(self.single_precision_measuring     ) + ",\n"
        str_json = str_json + "single_precision_coarsegraining : " + str(self.single_precision_coarsegraining) + ",\n"
        str_json = str_json + "use_full_vertex                 : " + str(self.use_full_vertex)  + ",\n"
        
        str_json = str_json + "machine_name : " + str(self.machine_name) + ",\n"
        str_json = str_json + "GPU_support  : " + str(self.GPU_support) + "\n"
        

        str_json = str_json + "}\n"
        
        return str_json

    def read(self):

        input = raw_input('Please enter Monte-Carlo algorithm : \n\t(1) CT-AUX \n\t(2) SS-Hyb \n\t(3) Cl-Hyb \n\t(4) PCM \n\t(5) High. Temp. Series-expansion\n\n')
        if   input == '1':
            self.QMC_type = "CT-AUX"
        elif input == '2':
            self.QMC_type = "SS-Hyb"
        elif input == '3':
            self.QMC_type = "Cl-Hyb"
        elif input == '4':
            self.QMC_type = "PCM"
        elif input == '5':
            self.QMC_type = "High. Temp. Series-expansion"

        input  = raw_input('Please enter execution-mode : \n\t(1) SERIAL \n\t(2) MPI \n\t(3) MPI_FT \n\t(4) OPENMPI_FT \n\n')
        if   input == '1':
            self.concurrency = "serial"
        elif input == '2':
            self.concurrency = "MPI"
        elif input == '3':
            self.concurrency = "MPI-FT"
        elif input == '4':
            self.concurrency = "OPENMPI-FT"
            
        input  = raw_input('do you wish sse-acceleration? (y/n) \n\n\n')
        if   input == 'y':
            self.use_sse = True
        else:
            self.use_sse = False

        input = raw_input('do you wish to have pthreads? (y/n) \n\n\n')
        if   input == 'y':
            self.use_pthreads = True
        else:
            self.use_pthreads = False
            
        input = raw_input('Please enter profiling-mode : \n\t(1) NONE \n\t(2) COUNTING \n\t(3) PAPI \n\n')
        if   input == '1':
            self.profiling = "no-profiling"
        elif input == '2':
            self.profiling = "counting-profiling"
        elif input == '3':
            self.profiling = "papi-profiling"
            
        if self.profiling == "counting-profiling" or self.profiling == "papi-profiling" :
            input = raw_input('enable profiling for auto-tuning ? (y/n) \n\n\n')

            if input == 'y':
                self.auto_tuning = "True"

        input = raw_input('Please enter random-number generator : \n\t(default) Numerical recipes \n\t(2) SPRNG \n\n')
        if   input == '2':
            self.rng_type = "sprng-type"
        else:
            self.rng_type = "NR"

        input = raw_input('do you wish to do BIT? (y/n) \n\n\n')
        if   input == 'y':
            self.do_BIT = True
        else:
            self.do_BIT = False

        input = raw_input('do you wish to assert? (y/n) \n\n\n')
        if   input == 'y':
            self.do_assert = True
        else:
            self.do_assert = False

        input = raw_input('do you want l1/l2/linf norm of errors? (y/n) \n\n\n')
        if   input == 'y':
            self.compute_std_dev = True
        else:
            self.compute_std_dev = False
         
        input = raw_input('do you wish to measure in single precision? (y/n) \n\n\n')
        if input == 'n':
            self.single_precision_measuring = False
        else:
            self.single_precision_measuring = True

        input = raw_input('do you wish to coarsegrain in single precision? (y/n) \n\n\n')
        if   input == 'n':
            self.single_precision_coarsegraining = False
        else:
            self.single_precision_coarsegraining = True

        input = raw_input('do you want GPU support? (y/n) \n\n\n')
        if   input == 'y':
            self.GPU_support = True

            input = raw_input('do you want pinned-host-memory? (y/n) \n\n\n')
            if   input == 'y':
                self.use_pinned_host_memory = True
            else:
                self.use_pinned_host_memory = False
        else:
            self.GPU_support            = False
            self.use_pinned_host_memory = False

        input = raw_input('Please enter machine : \n\t(1) MAC \n\t(2) LINUX \n\t(3) JAGUAR  \n\t(4) ROSA \n\t(5) TODI \n\t(6) EIGER \n\t(7) DAINT \n\t(8) TITAN  \n\n')
        if   input == '1':
            self.machine_name = "Mac"
        elif input == '2':
            self.machine_name = "Linux"
        elif input == '3':
            self.machine_name = "Jaguarpf"
        elif input == '4':
            self.machine_name = "Rosa"
        elif input == '5':
            self.machine_name = "Todi"
        elif input == '6':
            self.machine_name = "Eiger"
        elif input == '7':
            self.machine_name = "Daint"
        elif input == '8':
            self.machine_name = "Titan"

    def read_version(self):
        print "version \n\n"

        print "\n\t--> print version.txt"
        cmd = "cd ../..; git log -n 1 > version.txt"
        os.system(cmd)

        version_text=""

        for line in open("../../version.txt", "r"):
            line2 = line.strip()
            version_text = version_text + "\n\t\"" + line2 + "\\n\""

        return version_text

    def write_tests(self, dir_name, machine, device):

        if os.path.exists('../'+dir_name):
            shutil.rmtree('../'+dir_name)
        os.makedirs('../'+dir_name)

        self.machine_name = machine

        if(device=="GPU"):
            self.GPU_support            = True
            self.use_pinned_host_memory = True
        else:
            self.GPU_support            = False
            self.use_pinned_host_memory = False

        self.write_Makefile           (dir_name)
        self.write_compiler_directives(dir_name)
        
        self.write_includefile     (dir_name)
        self.write_include_file_GPU(dir_name)

        self.write_typedefs (dir_name)
        self.write_check_cpp(dir_name)

        self.write_tests_input(dir_name, machine, device)

    def write_tests_input(self, dir_name, machine, device):

        if os.path.exists('../'+dir_name):
            
            U_vec  = ["-4", "4"]
            mu_vec = ["-0.1", "0"]

            for l0 in range(0, len(U_vec)):
                for l1 in range(0, len(mu_vec)):

                    file = open("./test_input_files/input_check_repulsive.json", "r")
                    text = file.read()
                    file.close()
                    
                    text = text.replace("HUBBARDU", U_vec[l0])
                    text = text.replace("CHEM_POT", mu_vec[l1])
                    
                    file_name = "../"+dir_name+"/input_"+machine+"_device="+device+"_U="+U_vec[l0]+"_mu="+mu_vec[l1]+".json"

                    print "generating the input-file : " + file_name 
                    file = open(file_name, "w")
                    file.write(text)
                    file.close()


    def write_all(self):
        print "generating the bin-folder \n\n"
        
        if os.path.exists('../bin'):
            shutil.rmtree('../bin')
        os.makedirs('../bin')

        self.write_Makefile()
        
        self.write_compiler_directives()
        
        self.write_includefile()
        self.write_include_file_GPU()

        self.write_typedefs()
        
        self.write_main()

        self.write_analysis()

        self.write_check_cpp()

        self.write_generate_Hamiltonian_cpp()

    def write_Makefile(self, dir_name="bin"):
        print "generating the Makefile in ../"+dir_name 

        file = open("Makefile.tmpl", "r")
        text = file.read()
        file.close()

        if self.concurrency == "MPI-FT":
            text = text.replace("COMPILER_C++"   , "ftmpicxx")
        else:
            text = text.replace("COMPILER_C++"   , self.compiler_dic[self.machine_name])

        if self.GPU_support:
            text = text.replace("GPU_COMPILER"          , self.GPU_CC            [self.machine_name])
            text = text.replace("FLAGOPTIONS_GPU"       , self.GPU_flags_dic     [self.machine_name])
            text = text.replace("LIBRARIES_GPU"         , self.GPU_library_dic   [self.machine_name])
            #text = text.replace("MACHINE_MAGMA_ROOT_DIR", self.GPU_magma_root_dic[self.machine_name])
            text = text.replace("PROG_NAME"    , "dca_gpu")
        else:
            text = text.replace("PROG_NAME"    , "dca_cpu")
            text = text.replace("LIBRARIES_GPU", "")

        text = text.replace("FLAGOPTIONS"       , self.flags_dic   [self.machine_name])
        text = text.replace("LIBRARIES"         , self.library_dic [self.machine_name])
        text = text.replace("OTHER_INCLUDE_DIRS", self.include_dic [self.machine_name])

        text = text.replace("MAKE_FFT_THREADS", self.FFTW_THREADS_dic[self.machine_name])
        text = text.replace("OPTIONS_FFTW"    , self.FFTW_dic[self.machine_name])
        text = text.replace("OPTIONS_NFFT"    , self.NFFT_dic[self.machine_name])
        text = text.replace("OPTIONS_SPGLIB"  , self.SPGLIB_dic[self.machine_name])

        text = text.replace("OPTIONS_HDF5"    , self.HDF5_dic [self.machine_name])
        text = text.replace("OPTIONS_MPICH"   , self.MPICH_dic[self.machine_name])

        file = open("../"+dir_name+"/Makefile", "w")
        file.write(text) 
        file.close() 

    def write_typedefs(self, dir_name="bin"):
         print "generating the type-definitions"

         cmd = "cp ./type_definitions.h ../"+dir_name+"/type_definitions.h"
         os.system(cmd)

    def write_main(self, dir_name="bin"):
        print "generating the main.cpp"

        file = open("main.tmpl", "r")
        text = file.read()
        file.close()

        text = text.replace("\"DEFAULT\"", self.version)

        if  self.GPU_support:
            #text = text.replace("CUDA_FUNCTION"  , "void cuda_organizer(int status);")
            text = text.replace("CUDA_FUNCTION"  , "void print_device_info();")
            text = text.replace("CUDA_INITIALIZE", "cuda_organizer(0);")
            text = text.replace("CUDA_FINALIZE"  , "cuda_organizer(1);")
            text = text.replace("INITIALIZE_MAGMA_0", "void initialize_magma();")
            text = text.replace("INITIALIZE_MAGMA_1", "initialize_magma();")
            text = text.replace("LIN_ALG_DEVICE_TYPE", "LIN_ALG::GPU")
        else:
            text = text.replace("CUDA_FUNCTION"  , "void print_device_info(){}")
            text = text.replace("CUDA_INITIALIZE", "")
            text = text.replace("CUDA_FINALIZE"  , "")
            text = text.replace("INITIALIZE_MAGMA_0", "")
            text = text.replace("INITIALIZE_MAGMA_1", "")
            text = text.replace("void print_device_info();", "")
            text = text.replace("LIN_ALG_DEVICE_TYPE", "LIN_ALG::CPU")

        if   (self.concurrency == "serial"):    
            text = text.replace("PARALLELIZATION_LIBRARY_TYPE"     , "COMP_LIB::SERIAL_LIBRARY")
        elif (self.concurrency == "MPI"):   
            text = text.replace("PARALLELIZATION_LIBRARY_TYPE"     , "COMP_LIB::MPI_LIBRARY")

        if   (self.QMC_type == "CT-AUX"):    
            text = text.replace("CLUSTER_SOLVER_TYPE"    , "DCA::CT_AUX_CLUSTER_SOLVER")
        elif (self.QMC_type == "High. Temp. Series-expansion"):
            text = text.replace("CLUSTER_SOLVER_TYPE"    , "DCA::HIGH_TEMPERATURE_SERIES")

        if self.use_pthreads:
            text = text.replace("//typedef DCA::posix_qmci_integrator<quantum_cluster_solver_type> Monte_Carlo_Integrator_type;", 
                                "typedef DCA::posix_qmci_integrator<quantum_cluster_solver_type> Monte_Carlo_Integrator_type;")
        else:
            text = text.replace("//typedef quantum_cluster_solver_type                             Monte_Carlo_Integrator_type;",
                                "typedef quantum_cluster_solver_type                             Monte_Carlo_Integrator_type;")

        file = open("../"+dir_name+"/main.cpp", "w")
        file.write(text) 
        file.close() 

    def write_analysis(self, dir_name="bin"):
        print "generating the analysis.cpp"

        file = open("analysis.tmpl", "r")
        text = file.read()
        file.close()

        text = text.replace("\"DEFAULT\"", self.version)

        if   (self.concurrency == "serial"):    
            text = text.replace("MPI_LIB_TYPE"     , "SERIAL_LIBRARY")
        elif (self.concurrency == "MPI"):   
            text = text.replace("MPI_LIB_TYPE"     , "MPI_LIBRARY")
        elif (self.concurrency == "MPI-FT"):
            text = text.replace("MPI_LIB_TYPE"     , "MPI_FT_LIBRARY")
        elif (self.concurrency == "OPENMPI-FT"):
            text = text.replace("MPI_LIB_TYPE"     , "OPENMPI_FT_LIBRARY")

        if   (self.QMC_type == "CT-AUX"):    
            text = text.replace("MC_ALG_TYPE"    , "CT_AUX")
        elif (self.QMC_type == "SS-Hyb"):   
            text = text.replace("MC_ALG_TYPE"    , "HYBRIDIZATION")
        elif (self.QMC_type == "Cl-Hyb"):
            text = text.replace("MC_ALG_TYPE"    , "HYBRIDIZATION_FULL")
        elif (self.QMC_type == "PCM"):
            text = text.replace("MC_ALG_TYPE"    , "PCM")

        file = open("../"+dir_name+"/analysis.cpp", "w")
        file.write(text) 
        file.close() 

    def write_includefile(self, dir_name="bin"):
        print "generating the include_file.h"

        file = open("include_files.tmpl", "r")
        text = file.read()
        file.close()

        if (self.do_BIT):
             text = text.replace("QMC_INTEGRATOR_BIT = false;", "QMC_INTEGRATOR_BIT = true;")

        file = open("../"+dir_name+"/include_files.h", "w")
        file.write(text) 
        file.close() 

    def write_include_file_GPU(self, dir_name="bin"):
        print "generating the include_file.h"

        file = open("include_files_GPU.cu", "r")
        text = file.read()
        file.close()

        file = open("../"+dir_name+"/include_files_GPU.cu", "w")
        file.write(text) 
        file.close() 

    def write_check_cpp(self, dir_name="bin"):
        print "generating the ED check"

        file = open("./check.cpp", "r")
        text = file.read()
        file.close()

        file = open("../"+dir_name+"/check.cpp", "w")
        file.write(text) 
        file.close() 

    def write_generate_Hamiltonian_cpp(self, dir_name="bin"):
        print "generating the generate the Hamiltonian check"

        file = open("./generate_hamiltonian.cpp", "r")
        text = file.read()
        file.close()

        file = open("../"+dir_name+"/generate_hamiltonian.cpp", "w")
        file.write(text) 
        file.close() 

    def write_compiler_directives(self, dir_name="bin"):
        print "generating the compiler-directives"

        text = ""

        if not self.do_assert:    
            text = text + "#define NDEBUG" + "\n\n"

        if self.compute_std_dev :
            text = text + "#define MEASURE_ERROR_BARS" + "\n\n"

#         if not self.measure_in_double_precision:
#             text = text + "#define SINGLE_PRECISION_MEASUREMENTS" + "\n\n"

#         if self.single_precision_measuring:
#             text = text + "#define SINGLE_PRECISION_MEASUREMENTS" + "\n\n"

#         if self.single_precision_coarsegraining:
#             text = text + "#define SINGLE_PRECISION_COARSEGRAINING" + "\n\n"

        if not self.use_full_vertex :
            text = text + "#define USE_REDUCED_VERTEX_FUNCTION" + "\n\n"

        if   self.profiling == "counting-profiling":
            text = text + "#define COUNTING_PROFILING" + "\n\n"
        elif self.profiling == "papi-profiling":
            text = text + "#define DCA_PAPI_PROFILING" + "\n\n"
        else:
            text = text + "#define NO_PROFILING" + "\n\n"

        if self.auto_tuning == "True":
            text = text + "#define AUTOTUNING_ENABLED" + "\n\n"

        if self.rng_type == "SPRNG":    
            text = text + "#define RNG_SPRNG" + "\n\n"
        else:
            text = text + "#define RNG_NUMERICAL_RECIPES" + "\n\n"

        if self.GPU_support:
            text = text + "#define USE_GPU" + "\n\n"

        if self.use_sse:
            text = text + "#define USE_SSE_ACCELERATION" + "\n\n"

        if self.use_pinned_host_memory:
            text = text + "#define ENABLE_PINNED_MEMORY_ALLOCATION" + "\n\n"

        file = open("../"+dir_name+"/compiler_directives.h", "w")
        file.write(text) 
        file.close() 
