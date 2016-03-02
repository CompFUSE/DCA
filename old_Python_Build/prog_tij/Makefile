
SHELL          = /bin/bash
SRC_DIR        = ../../src
MAGMA_ROOT_DIR = MACHINE_MAGMA_ROOT_DIR

LIBS___DIR = ../../libs
FFTW___DIR = ../../libs/fftw-3.3.3
NFFT___DIR = ../../libs/nfft-3.1.3
HDF5___DIR = ../../libs/hdf5-1.8.11
SPGLIB_DIR = ../../libs/spglib-1.4.1
MPICH__DIR = ../../libs/mpich-3.0.4

COMP_DIR = ../../src/comp_libraries
MATH_DIR = ../../src/math_algorithms
PHYS_DIR = ../../src/phys_algorithms/DCA+


#CC   = $(MPICH__DIR)/build/bin/mpicxx
CC = mpicxx-mpich-clang34
#CC = mpicxx-mpich-gcc48
NVCC = GPU_COMPILER

CLANGCC   = /opt/local/bin/clang-mp-3.4
CLANGCPPC = /opt/local/bin/clang++-mp-3.4

GCC = /opt/local/bin/gcc-mp-4.8
GPP = /opt/local/bin/g++-mp-4.8


FLAGS_ANA    = -g3 -O3   -Wall -Werror #-funroll-loops -finline-functions

FLAGS      = -g3 -O3  -Wall -Werror #-Wno-sign-compare -Wno-unused-local-typedefs -fopenmp #-ferror-limit=100  #-funroll-loops -finline-functions
FLAGS_NVCC = -g3 -O3 -Wall -Werror #-funroll-loops -finline-functions_GPU

LIBS     = -llapack -lblas -lm \
	$(HDF5___DIR)/build/lib/libhdf5.8.dylib \
	$(HDF5___DIR)/build/lib/libhdf5.a \
	$(HDF5___DIR)/build/lib/libhdf5_cpp.8.dylib \
	$(HDF5___DIR)/build/lib/libhdf5_cpp.a \
	$(NFFT___DIR)/build/lib/libnfft3.a \
	$(FFTW___DIR)/build/lib/libfftw3.a \
	$(SPGLIB_DIR)/build/lib/libsymspg.a

LIBS_GPU = 

TMP = 	-I$(MPICH__DIR)/build/include \
	-I$(MPICH__DIR)/build/lib \

INCLUDES =  -I. \
	-I$(FFTW___DIR)/build/include \
	-I$(FFTW___DIR)/build/lib \
	-I$(NFFT___DIR)/build/include \
	-I$(NFFT___DIR)/build/lib \
	-I$(HDF5___DIR)/build/include \
	-I$(HDF5___DIR)/build/lib \
	-I$(SPGLIB_DIR)/build/include/spglib \
	-I$(SRC_DIR) \
	-I$(SRC_DIR)/comp_library/JSON \
	-I$(SRC_DIR)/comp_library/JSON/JsonParser \
	-I$(SRC_DIR)/comp_library/IO_operations \
	-I$(SRC_DIR)/comp_library/IO_operations/CSV \
	-I$(SRC_DIR)/comp_library/profiler_library \
	-I$(SRC_DIR)/comp_library/profiler_library/events \
	-I$(SRC_DIR)/comp_library/profiler_library/profilers \
	-I$(SRC_DIR)/comp_library/parallelization_library/ \
	-I$(SRC_DIR)/comp_library/parallelization_library/interfaces\
	-I$(SRC_DIR)/comp_library/type_list \
	-I$(SRC_DIR)/LAPACK_PLANS \
	-I$(SRC_DIR)/LAPACK_PLANS/BLAS \
	-I$(SRC_DIR)/LAPACK_PLANS/CUBLAS \
	-I$(SRC_DIR)/LAPACK_PLANS/LAPACK \
	-I$(SRC_DIR)/LAPACK_PLANS/LAPACK/eigensystem_plans \
	-I$(SRC_DIR)/LAPACK_PLANS/LAPACK/gelss_plans \
	-I$(SRC_DIR)/LAPACK_PLANS/LAPACK/geqr_plans \
	-I$(SRC_DIR)/LAPACK_PLANS/MAGMA \
	-I$(SRC_DIR)/LAPACK_PLANS/BLAS_BLOCK_MATRICES \
	-I$(SRC_DIR)/LAPACK_PLANS/LINEAR_ALGEBRA_PLANS \
	-I$(SRC_DIR)/comp_library/generic_methods_library \
	-I$(SRC_DIR)/comp_library/LIN_ALG \
	-I$(SRC_DIR)/comp_library/LIN_ALG/src \
	-I$(SRC_DIR)/comp_library/LIN_ALG/src/linalg_operations \
	-I$(SRC_DIR)/comp_library/LIN_ALG/src/linalg_structures \
	-I$(SRC_DIR)/comp_library/IO_library/ \
	-I$(SRC_DIR)/comp_library/IO_library/CSV \
	-I$(SRC_DIR)/comp_library/IO_library/HDF5 \
	-I$(SRC_DIR)/comp_library/IO_library/JSON \
	-I$(SRC_DIR)/comp_library/IO_library/JSON/JSON_PARSER \
	-I$(SRC_DIR)/comp_library/function_library\
	-I$(SRC_DIR)/comp_library/function_library/function_operations \
	-I$(SRC_DIR)/comp_library/function_library/domains \
	-I$(SRC_DIR)/comp_library/function_library/domains/special_domains \
	-I$(SRC_DIR)/comp_library/function_library/domains/product_domains \
	-I$(SRC_DIR)/comp_library/function_library/domains/union_domains \
	-I$(SRC_DIR)/comp_library/function_plotting \
	-I$(SRC_DIR)/comp_library/function_plotting/gnuplot \
	-I$(SRC_DIR)/math_library \
	-I$(SRC_DIR)/math_library/statistics \
	-I$(SRC_DIR)/math_library/random_number_library \
	-I$(SRC_DIR)/math_library/geometry_library \
	-I$(SRC_DIR)/math_library/geometry_library/vector_operations \
	-I$(SRC_DIR)/math_library/geometry_library/tetrahedron_mesh/ \
	-I$(SRC_DIR)/math_library/geometry_library/tetrahedron_mesh/tetrahedron_implementation \
	-I$(SRC_DIR)/math_library/geometry_library/tetrahedron_mesh/tetrahedron_mesh_implementation \
	-I$(SRC_DIR)/math_library/geometry_library/gaussian_quadrature \
	-I$(SRC_DIR)/math_library/functional_transforms/ \
	-I$(SRC_DIR)/math_library/functional_transforms/domain_specifications \
	-I$(SRC_DIR)/math_library/functional_transforms/basis_functions \
	-I$(SRC_DIR)/math_library/functional_transforms/basis_functions/hermite_splines \
	-I$(SRC_DIR)/math_library/functional_transforms/basis_transforms \
	-I$(SRC_DIR)/math_library/functional_transforms/domain_transforms \
	-I$(SRC_DIR)/math_library/functional_transforms/function_transforms \
	-I$(SRC_DIR)/math_library/inference_library/ \
	-I$(SRC_DIR)/math_library/interpolation_library/ \
	-I$(SRC_DIR)/math_library/gaussian_processes \
	-I$(SRC_DIR)/math_library/gaussian_processes/covariance_function \
	-I$(SRC_DIR)/math_library/FFT_and_NFFT/ \
	-I$(SRC_DIR)/math_library/FFT_and_NFFT/NFFT/ \
	-I$(SRC_DIR)/math_library/FFT_and_NFFT/NFFT/domains/ \
	-I$(SRC_DIR)/math_library/FFT_and_NFFT/NFFT/basis_functions/ \
	-I$(SRC_DIR)/math_library/FFT_and_NFFT/NFFT/atomic_convolutions/ \
	-I$(SRC_DIR)/algorithms \
	-I$(SRC_DIR)/algorithms/compute_DOS \
	-I$(SRC_DIR)/algorithms/compute_Fermi_surface \
	-I$(SRC_DIR)/algorithms/Bsplines \
	-I$(SRC_DIR)/algorithms/Bsplines/uniform \
	-I$(SRC_DIR)/algorithms/simple_algorithms \
	-I$(SRC_DIR)/algorithms/fourier_transform \
	-I$(SRC_DIR)/algorithms/fourier_transform/FT_type_2_type \
	-I$(SRC_DIR)/algorithms/legendre_transform \
	-I$(SRC_DIR)/algorithms/legendre_transform/legendre_transformation_kernel \
	-I$(SRC_DIR)/algorithms/generic_assert \
	-I$(SRC_DIR)/algorithms/concepts \
	-I$(SRC_DIR)/algorithms/cluster_operations \
	-I$(SRC_DIR)/algorithms/interpolation \
	-I$(SRC_DIR)/algorithms/interpolation/extended_interpolation_domains \
	-I$(SRC_DIR)/algorithms/interpolation/bspline_interpolation_kernel \
	-I$(SRC_DIR)/algorithms/interpolation/hspline_interpolation_kernel \
	-I$(SRC_DIR)/algorithms/interpolation/lspline_interpolation_kernel \
	-I$(SRC_DIR)/algorithms/interpolation/linear_interpolation_kernel \
	-I$(SRC_DIR)/algorithms/interpolation/rbf_interpolation_kernel \
	-I$(SRC_DIR)/algorithms/interpolation/wannier_interpolation_kernel \
	-I$(SRC_DIR)/algorithms/interpolation/projection_kernel \
	-I$(SRC_DIR)/algorithms/smooth_algorithms \
	-I$(SRC_DIR)/algorithms/integration \
	-I$(SRC_DIR)/algorithms/integration/inner_product_kernels \
	-I$(SRC_DIR)/algorithms/coarsegrain \
	-I$(SRC_DIR)/algorithms/coarsegrain/delta_function_integration \
	-I$(SRC_DIR)/algorithms/coarsegrain/delta_function_integration/implementation \
	-I$(SRC_DIR)/algorithms/coarsegrain/isolated_cluster_integration \
	-I$(SRC_DIR)/algorithms/coarsegrain/isolated_cluster_integration/implementation \
	-I$(SRC_DIR)/algorithms/coarsegrain/trapezium_integration \
	-I$(SRC_DIR)/algorithms/coarsegrain/trapezium_integration/implementation \
	-I$(SRC_DIR)/algorithms/coarsegrain/quadrature_integration \
	-I$(SRC_DIR)/algorithms/coarsegrain/quadrature_integration/quadrature_domains \
	-I$(SRC_DIR)/algorithms/coarsegrain/quadrature_integration/quadrature_algorithms \
	-I$(SRC_DIR)/algorithms/coarsegrain/quadrature_integration/quadrature_implementation \
	-I$(SRC_DIR)/algorithms/coarsegrain/tetrahedron_integration \
	-I$(SRC_DIR)/algorithms/coarsegrain/tetrahedron_integration/implementation \
	-I$(SRC_DIR)/algorithms/coarsegrain_inversion \
	-I$(SRC_DIR)/algorithms/coarsegrain_inversion/deconvolute \
	-I$(SRC_DIR)/algorithms/do_cast \
	-I$(SRC_DIR)/algorithms/mesh \
	-I$(SRC_DIR)/algorithms/mesh/delta_function_mesh/ \
	-I$(SRC_DIR)/algorithms/mesh/tetrahedron_mesh/ \
	-I$(SRC_DIR)/algorithms/mesh/tetrahedron_mesh/tetrahedron_mesh_implementation \
	-I$(SRC_DIR)/algorithms/mesh/tetrahedron_mesh/tetrahedron_implementation \
	-I$(SRC_DIR)/algorithms/random_numbers \
	-I$(SRC_DIR)/algorithms/statistics \
	-I$(SRC_DIR)/algorithms/provenance \
	-I$(SRC_DIR)/algorithms/pade \
	-I$(SRC_DIR)/algorithms/pade/basis_functions \
	-I$(SRC_DIR)/algorithms/pade/gradient_method \
	-I$(SRC_DIR)/algorithms/pade/levmar_method \
	-I$(SRC_DIR)/algorithms/pade/conjugate_gradient_method \
	-I$(SRC_DIR)/algorithms/print_data \
	-I$(SRC_DIR)/algorithms/profilers \
	-I$(SRC_DIR)/algorithms/profilers/profilers \
	-I$(SRC_DIR)/algorithms/profilers/events \
	-I$(SRC_DIR)/algorithms/symmetrization \
	-I$(SRC_DIR)/algorithms/compute_moments \
	-I$(SRC_DIR)/algorithms/compute_band_structure \
	-I$(SRC_DIR)/algorithms/compute_chemical_potential \
	-I$(SRC_DIR)/algorithms/Newton_Rapson \
	-I$(SRC_DIR)/phys_library/domains/ \
	-I$(SRC_DIR)/phys_library/domains/cluster \
	-I$(SRC_DIR)/phys_library/domains/cluster/interpolation \
	-I$(SRC_DIR)/phys_library/domains/cluster/interpolation/wannier_interpolation_kernel \
	-I$(SRC_DIR)/phys_library/domains/cluster/interpolation/hspline_interpolation_kernel \
	-I$(SRC_DIR)/phys_library/domains/cluster/symmetries \
	-I$(SRC_DIR)/phys_library/domains/cluster/symmetries/point_groups \
	-I$(SRC_DIR)/phys_library/domains/cluster/symmetries/point_groups/2D \
	-I$(SRC_DIR)/phys_library/domains/cluster/symmetries/point_groups/3D \
	-I$(SRC_DIR)/phys_library/domains/cluster/symmetries/symmetry_operations \
	-I$(SRC_DIR)/phys_library/domains/cluster/symmetries/symmetry_operations/trigoniometric_ops \
	-I$(SRC_DIR)/phys_library/domains/cluster/symmetries/symmetry_operations/2D \
	-I$(SRC_DIR)/phys_library/domains/cluster/symmetries/symmetry_operations/3D \
	-I$(SRC_DIR)/phys_library/domains/cluster/symmetrization_algorithms \
	-I$(SRC_DIR)/phys_library/domains/time_and_frequency \
	-I$(SRC_DIR)/phys_library/domains/Quantum_domain \
	-I$(SRC_DIR)/phys_library/domains/Quantum_domain/Brillouin_zone_cut \
	-I$(SRC_DIR)/phys_library/DFT_connection/VASP/ \
	-I$(SRC_DIR)/phys_library/DFT_connection/VASP/vasp_domains \

all: t_ij #dca_cpu #analysis

showEnv:
	@echo 
	@echo ------------------------------------------------ Compiler Version
	$(CC) --version	
	@echo ----------------------------------------------------------------
	@echo

showEnv_GPU:
	@echo 
	@echo ------------------------------------------------ Compiler Version
	$(NVCC) --version	
	@echo ----------------------------------------------------------------
	@echo


main.o: showEnv main.cpp
	time $(CC)  $(INCLUDES) $(FLAGS) -c main.cpp

t_ij: clean main.o
	time $(CC) $(INCLUDES) $(FLAGS) main.o -o t_ij $(LIBS)

FFTW:
	rm -rf `pwd`/$(FFTW___DIR)
	tar -xvf    $(LIBS___DIR)/fftw-3.3.3.tar
	mv  `pwd`/fftw-3.3.3 `pwd`/$(LIBS___DIR)/fftw-3.3.3
	mkdir `pwd`/$(FFTW___DIR)/build
	cd `pwd`/$(FFTW___DIR); ./configure --prefix=`pwd`/build --disable-fortran #--prefix=`pwd`/build --enable-cxx CC=$(CC) CXX=$(CPPC)
	cd `pwd`/$(FFTW___DIR); make -j 8
	cd `pwd`/$(FFTW___DIR); make install

# FFTW:
# 	rm -rf `pwd`/../../fftw-3.2.2
# 	tar -xvf    ../../fftw-3.2.2.tar
# 	mv  `pwd`/fftw-3.2.2 `pwd`/../../fftw-3.2.2
# 	mkdir `pwd`/../../fftw-3.2.2/build
# 	cd `pwd`/../../fftw-3.2.2; ./configure --prefix=`pwd`/build --disable-fortran
# 	cd `pwd`/../../fftw-3.2.2; make -j 8
# 	cd `pwd`/../../fftw-3.2.2; make install

NFFT:
	rm -rf `pwd`/$(NFFT___DIR)
	tar -xvf  $(LIBS___DIR)/nfft-3.1.3.tar
	mv  `pwd`/nfft-3.1.3 `pwd`/$(NFFT___DIR)
	mkdir `pwd`/$(NFFT___DIR)/build
	cd `pwd`/$(NFFT___DIR); ./configure --prefix=`pwd`/build --with-fftw3=`pwd`/../fftw-3.2.2/build --without-apple-gcc-arch #--prefix=`pwd`/build --with-fftw3=`pwd`/../fftw-3.3.3/build --without-apple-gcc-arch
	cd `pwd`/$(NFFT___DIR); make -j 8
	cd `pwd`/$(NFFT___DIR); make install

# NFFT:
# 	rm -rf `pwd`/../../nfft-3.1.3
# 	tar -xvf ../../nfft-3.1.3.tar
# 	mv  `pwd`/nfft-3.1.3 `pwd`/../../nfft-3.1.3
# 	mkdir `pwd`/../../nfft-3.1.3/build
# 	cd `pwd`/../../nfft-3.1.3; ./configure --prefix=`pwd`/build --with-fftw3=`pwd`/../fftw-3.2.2/build --without-apple-gcc-arch
# 	cd `pwd`/../../nfft-3.1.3; make -j 8
# 	cd `pwd`/../../nfft-3.1.3; make install

SPGLIB:
	rm -rf `pwd`/$(SPGLIB_DIR)
	tar -xvf    $(LIBS___DIR)/spglib-1.4.1.tar
	mv  `pwd`/spglib-1.4.1 `pwd`/$(SPGLIB_DIR)
	mkdir `pwd`/$(SPGLIB_DIR)/build
	cd `pwd`/$(SPGLIB_DIR); ./configure --prefix=`pwd`/build --enable-cxx CC=$(CC) CXX=$(CPPC) #--prefix=`pwd`/build  --enable-cxx CC=$(GCC) CXX=$(GPP)
	cd `pwd`/$(SPGLIB_DIR); make -j 8
	cd `pwd`/$(SPGLIB_DIR); make install

# SPGLIB : 
# 	rm -rf `pwd`/../../spglib-1.4.1
# 	tar -xvf ../../spglib-1.4.1.tar
# 	mv  `pwd`/spglib-1.4.1 `pwd`/../../spglib-1.4.1
# 	mkdir `pwd`/../../spglib-1.4.1/build
# 	cd `pwd`/../../spglib-1.4.1; ./configure --prefix=`pwd`/build
# 	cd `pwd`/../../spglib-1.4.1; make
# 	cd `pwd`/../../spglib-1.4.1; make install

HDF5:
	rm -rf `pwd`/$(HDF5___DIR)
	tar -xvf    $(LIBS___DIR)/hdf5-1.8.11.tar
	mv  `pwd`/hdf5-1.8.11 `pwd`/$(HDF5___DIR)
	mkdir `pwd`/$(HDF5___DIR)/build
	cd `pwd`/$(HDF5___DIR); ./configure --prefix=`pwd`/build --enable-cxx CC=$(CLANGCC) CXX=$(CLANGCPPC)
	cd `pwd`/$(HDF5___DIR); make -j 8
	cd `pwd`/$(HDF5___DIR); make install

# LEVMAR:
# 	rm -rf `pwd`../../levmar-2.5
# 	tar xfvz ../../levmar-2.5-custom.gz
# 	mv  `pwd`/levmar-2.5 `pwd`/../../levmar-2.5
# 	cd `pwd`/../../levmar-2.5/; make

MPICH:
	rm -rf `pwd`/$(MPICH__DIR)
	tar -xvf    $(LIBS___DIR)/mpich-3.0.4.tar
	mv  `pwd`/mpich-3.0.4 `pwd`/$(LIBS___DIR)/
	mkdir `pwd`/$(MPICH__DIR)/build
	cd `pwd`/$(MPICH__DIR); ./configure --prefix=`pwd`/build --enable-cxx CC=$(CLANGCC) CXX=$(CLANGCPPC) --disable-f77 --disable-fc
	cd `pwd`/$(MPICH__DIR); make -j 8
	cd `pwd`/$(MPICH__DIR); make install

clean:
	rm -f dca_cpu dca_gpu main.o

cleananalysis:
	rm -f analysis_cpu analysis_gpu analysis.o

cleancheck:
	rm -f check check.o

