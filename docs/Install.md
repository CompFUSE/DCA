# Installation

### Requirements

#### FFTW3

DCA++ uses the FFTW3 interface to perform FFT tasks. Therefore, in order to build DCA++, FFTW3 or an FFT library with an FFTW3 wrapper (e.g. Intel MKL) has to be provided. When running CMake, the FFTW3(-compatible) library and the path to `fftw3.h` have to be set:

	cmake -DFFTW_INCLUDE_DIR=/path/to/fftw3.h -DFFTW_LIBRARY=/path/to/FFTW3-library ...
	
Note that if you need to link against multiple libraries, they have to be separated by semicolons and enclosed by double quotes:

	-DFFTW_LIBRARY="lib1.a;lib2.a"

If you use a compiler wrapper that already adds the include path for your FFTW3(-compatible) library and automatically links it, use

	cmake -DDCA_HAVE_FFTW=TRUE ...