
SHELL    = /bin/bash
SRC_DIR  = ../src

INCLUDES =  -I. \
	    -I$(SRC_DIR) \
	    -I$(SRC_DIR)/linalg_operations \
	    -I/apps/eiger/Cuda-4.0/cuda/include \
            -I/project/s299/magma/magma-1.1_my-castor-gcc-4.4.5-mkl/ \
	    -I/project/s299/magma/magma-1.1_my-castor-gcc-4.4.5-mkl/include \
	    -L/project/s299/magma/magma-1.1_my-castor-gcc-4.4.5-mkl/lib \
#            -I/project/s299/magma/magma-1.1_my-eiger-gcc-4.3.4-mkl/ \
#	    -I/project/s299/magma/magma-1.1_my-eiger-gcc-4.3.4-mkl/include \
#	    -L/project/s299/magma/magma-1.1_my-eiger-gcc-4.3.4-mkl/lib \
#            -I/project/s299/magma/magma-1.1_my-eiger-gcc-4.3.4-acml \
#            -I/project/s299/magma/magma-1.1_my-eiger-gcc-4.3.4-acml/include


CC           = g++
NVCC         = nvcc

FLAGS        = -O3 -Wall -Werror #-funroll-loops -ffast-math -msse4 -finline-functions
FLAGS_NVCC   = -arch=sm_20

#LIBS         = -lblas -llapack
#GPULIBS         = -L/apps/eiger/Cuda-4.0/cuda/lib64 -lmagma -lmagmablas -lmagma -lcuda -lcudart -lcublas -L/apps/eiger/Intel-FOR-11.1/mkl/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core 
GPULIBS         = -L/apps/castor/CUDA-4.0/cuda/lib64 -lmagma -lmagmablas -lmagma -lcuda -lcudart -lcublas -L/apps/castor/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core 

all: showEnv test_GPU

showEnv:
	@echo 
	@echo ------------------------------------------------ Compiler Version
	$(CC) --version	
	$(NVCC) --version	
	@echo ----------------------------------------------------------------
	@echo

main.o: main.cpp
	$(CC) $(INCLUDES) $(FLAGS) -c main.cpp

test: clean  main.o 
	$(CC) $(INCLUDES) $(FLAGS) main.o  -o test $(LIBS) 

CPU_code.o: main.cpp
	$(CC) $(INCLUDES) $(FLAGS) -c main.cpp

GPU_code.o: kernels.cu
	$(NVCC) $(INCLUDES) $(FLAGS_NVCC) -c kernels.cu

test_GPU: clean_GPU GPU_code.o CPU_code.o 
	$(CC) $(INCLUDES) $(FLAGS) kernels.o main.o -o test $(GPULIBS) 

clean:
	rm -f test *.o 

clean_GPU:
	rm -f test *.o


