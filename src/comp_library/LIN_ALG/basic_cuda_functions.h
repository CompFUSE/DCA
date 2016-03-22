//-*-C++-*-

#ifndef BASIC_CUDA_FUNCTIONS_H
#define BASIC_CUDA_FUNCTIONS_H

void print_device()
{
  int ndevices;
  cuDeviceGetCount(&ndevices);

  for(int idevice=0; idevice<ndevices; idevice++)
    {
      char   name[200];
      size_t totalMem;

      int      clock;
      CUdevice dev;

      cuDeviceGet( &dev, idevice );
      
      cuDeviceGetName( name, sizeof(name), dev );
      
      cuDeviceTotalMem(&totalMem, dev);
      
      cuDeviceGetAttribute( &clock,
                            CU_DEVICE_ATTRIBUTE_CLOCK_RATE, dev );
      
      printf( "device %d: %s, %.1f MHz clock, %.1f MB memory\n",
              idevice, name, clock/1000.f, totalMem/1024.f/1024.f );
    }
};

void print_device_info()
{  
  // Number of CUDA devices
  int devCount;
  cudaGetDeviceCount(&devCount);
  
  std::cout << "\n\tCUDA Device Query..."<< "\n";
  std::cout << "\tThere are " << devCount << " CUDA devices.\n";
 
  // Iterate through devices
  for(int i=0; i<devCount; ++i)
    {
      std::cout << "\n\tCUDA Device " << i << "\n";
      
      cudaDeviceProp devProp;
      cudaGetDeviceProperties(&devProp, i);
      
      std::cout << "\tMajor revision number:         " <<   devProp.major<< "\n";
      std::cout << "\tMinor revision number:         " <<   devProp.minor<< "\n";
      std::cout << "\tName:                          " <<   devProp.name<< "\n";
      std::cout << "\tTotal global memory:           " <<  devProp.totalGlobalMem << "\n";
      std::cout << "\tTotal shared memory per block: " <<  devProp.sharedMemPerBlock<< "\n";
      std::cout << "\tTotal registers per block:     " <<   devProp.regsPerBlock<< "\n";
      std::cout << "\tWarp size:                     " <<   devProp.warpSize<< "\n";
      std::cout << "\tMaximum memory pitch:          " <<   devProp.memPitch<< "\n";
      std::cout << "\tMaximum threads per block:     " <<   devProp.maxThreadsPerBlock<< "\n";
      
      for(int i = 0; i < 3; ++i)
	std::cout << "\tMaximum dimension " << i << " of block:  " << devProp.maxThreadsDim[i]<< "\n";
      
      for(int i = 0; i < 3; ++i)
	std::cout << "\tMaximum dimension " << i << " of grid:   " << devProp.maxGridSize[i]<< "\n";
      
      std::cout << "\tClock rate:                    " <<   devProp.clockRate<< "\n";
      std::cout << "\tTotal constant memory:         " <<  devProp.totalConstMem<< "\n";
      std::cout << "\tTexture alignment:             " <<  devProp.textureAlignment<< "\n";
      std::cout << "\tConcurrent copy and execution: " <<  (devProp.deviceOverlap ? "Yes" : "No")<< "\n";
      std::cout << "\tNumber of multiprocessors:     " <<   devProp.multiProcessorCount<< "\n";
      std::cout << "\tKernel execution timeout:      " <<  (devProp.kernelExecTimeoutEnabled ? "Yes" : "No")<< "\n";
    }
};

void initialize_magma()
{
  magma_init();
}



#ifdef DEBUG_CUDA

bool cuda_check_for_errors(std::string function_name, std::string file_name, int line)
{
  cudaError_t ret = cudaGetLastError();

  if(ret != cudaSuccess)
    {
      std::stringstream ss;

      ss << "\n\n\t error in function : " << function_name 
	 << "\t file-name : "             << file_name 
	 << "\t at line : "               << line 
	 << "\t error : " << cudaGetErrorString(ret) << "\n\n";

      std::cout << ss.str();

      throw std::logic_error(ss.str());
    }

  return true;
}

bool cuda_check_for_errors_bgn(std::string function_name, std::string file_name, int line)
{
  cudaError_t ret = cudaGetLastError();

  if(ret != cudaSuccess)
    {
      std::stringstream ss;

      ss << "\n\n\t error in function (begin) : " << function_name 
	 << "\t file-name : "             << file_name 
	 << "\t at line : "               << line 
	 << "\t error : " << cudaGetErrorString(ret) << "\n\n";

      std::cout << ss.str();

      throw std::logic_error(ss.str());
    }

  return true;
}

bool cuda_check_for_errors_end(std::string function_name, std::string file_name, int line)
{
  cudaDeviceSynchronize();

  cudaError_t ret = cudaGetLastError();

  if(ret != cudaSuccess)
    {
      std::stringstream ss;

      ss << "\n\n\t error in function (end) : " << function_name 
	 << "\t file-name : "             << file_name 
	 << "\t at line : "               << line 
	 << "\t error : " << cudaGetErrorString(ret) << "\n\n";

      std::cout << ss.str();

      throw std::logic_error(ss.str());
    }

  return true;
}

#endif 

cudaDeviceProp& get_device_properties()
{
    static bool           init=false;
    static cudaDeviceProp prop;

    if(!init){
	int count;

	cudaGetDeviceCount(&count);

	if(count<1)
            throw std::logic_error(__FUNCTION__);

	cudaGetDeviceProperties(&prop, 0);

	init = true;
    }

    return prop;
}

int get_number_of_threads(){
  const static int N_th = get_device_properties().maxThreadsPerBlock;
  return N_th;
}

int get_number_of_blocks(int n){
  const static int N_th = get_device_properties().maxThreadsPerBlock;
  return (n+N_th-1)/N_th;
}

int get_number_of_blocks(int n, int N_th){
  return (n+N_th-1)/N_th;
}

void synchronize_devices()
{
  cudaDeviceSynchronize();
}


#endif
