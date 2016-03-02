//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                       

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <complex>
#include <stdexcept>
#include <assert.h>
#include <limits>
#include <algorithm>
#include <typeinfo>
#include <pthread.h>
#include <omp.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <dirent.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <queue>
#include <bitset>

#define NDEBUG
//#define DEBUG_CUDA
//#define cudaDeviceScheduleBlockingSync 0x04

// using namespace std;

#include "include_linalg.cu.h"


#include "ctaux_G0_matrix_routines_GPU.cu.h"
#include "ctaux_G_matrix_routines_GPU.cu.h"
#include "ctaux_N_matrix_routines_GPU.cu.h"
//#include ""

#include "ctaux_walker_routines_GPU.cu.h"

// #include "SHRINK_ALGORITHMS_GPU.cu.h"
// #include "G_MATRIX_TOOLS_GPU.cu.h"
// #include "N_MATRIX_TOOLS_GPU.cu.h"
// #include "G0_INTERPOLATION_GPU.cu.h"
// #include "GAMMA_MATRIX_TOOLS_GPU.cu.h"

// #include "CT_AUX_WALKER_TOOLS_GPU.cu.h"


