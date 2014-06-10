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
  
  cout << "\n\tCUDA Device Query..."<< "\n";
  cout << "\tThere are " << devCount << " CUDA devices.\n";
 
  // Iterate through devices
  for(int i=0; i<devCount; ++i)
    {
      cout << "\n\tCUDA Device " << i << "\n";
      
      cudaDeviceProp devProp;
      cudaGetDeviceProperties(&devProp, i);
      
      cout << "\tMajor revision number:         " <<   devProp.major<< "\n";
      cout << "\tMinor revision number:         " <<   devProp.minor<< "\n";
      cout << "\tName:                          " <<   devProp.name<< "\n";
      cout << "\tTotal global memory:           " <<  devProp.totalGlobalMem << "\n";
      cout << "\tTotal shared memory per block: " <<  devProp.sharedMemPerBlock<< "\n";
      cout << "\tTotal registers per block:     " <<   devProp.regsPerBlock<< "\n";
      cout << "\tWarp size:                     " <<   devProp.warpSize<< "\n";
      cout << "\tMaximum memory pitch:          " <<   devProp.memPitch<< "\n";
      cout << "\tMaximum threads per block:     " <<   devProp.maxThreadsPerBlock<< "\n";
      
      for(int i = 0; i < 3; ++i)
	cout << "\tMaximum dimension " << i << " of block:  " << devProp.maxThreadsDim[i]<< "\n";
      
      for(int i = 0; i < 3; ++i)
	cout << "\tMaximum dimension " << i << " of grid:   " << devProp.maxGridSize[i]<< "\n";
      
      cout << "\tClock rate:                    " <<   devProp.clockRate<< "\n";
      cout << "\tTotal constant memory:         " <<  devProp.totalConstMem<< "\n";
      cout << "\tTexture alignment:             " <<  devProp.textureAlignment<< "\n";
      cout << "\tConcurrent copy and execution: " <<  (devProp.deviceOverlap ? "Yes" : "No")<< "\n";
      cout << "\tNumber of multiprocessors:     " <<   devProp.multiProcessorCount<< "\n";
      cout << "\tKernel execution timeout:      " <<  (devProp.kernelExecTimeoutEnabled ? "Yes" : "No")<< "\n";
    }
};

void initialize_magma()
{
  magma_init();
}

/*
void cuda_organizer(int status)
{
  static CUdevice  dev;                                                        
  static CUcontext context;                                                    
  
  if(status == 0)
    {
      if( CUDA_SUCCESS != cuInit(0) ) {                                   
	fprintf(stderr, "CUDA: Not initialized\n" ); exit(-1);                
      }                                                                     
      
      if( CUDA_SUCCESS != cuDeviceGet(&dev, 0) ) {			  
	fprintf(stderr, "CUDA: Cannot get the device\n"); exit(-1);           
      }                                                                     
      
      if( CUDA_SUCCESS != cuCtxCreate( &context, 0, dev ) ) {		  
	fprintf(stderr, "CUDA: Cannot create the context\n"); exit(-1);       
      }                                                                     
      
//       if( CUBLAS_STATUS_SUCCESS != cublasInit() ) {                        
// 	fprintf(stderr, "CUBLAS: Not initialized\n"); exit(-1);               
//       }    
    }
  else
    {
      cuCtxDetach( context ); 
//       cublasShutdown();
    }
};
*/

/*
bool cuda_check_for_errors()
{
  cudaDeviceSynchronize();

  cudaError_t ret = cudaGetLastError();

  if(ret != cudaSuccess)
    {
      std::string error_msg(cudaGetErrorString(ret));

      cout << error_msg << endl;

      throw std::logic_error(__FUNCTION__);
    }

  return true;
}
*/

#ifdef DEBUG_CUDA

bool cuda_check_for_errors(std::string function_name, std::string file_name, int line)
{
  //cout << __FUNCTION__ << endl;

  //cudaDeviceSynchronize();

  cudaError_t ret = cudaGetLastError();

  if(ret != cudaSuccess)
    {
      std::stringstream ss;

      ss << "\n\n\t error in function : " << function_name 
	 << "\t file-name : "             << file_name 
	 << "\t at line : "               << line 
	 << "\t error : " << cudaGetErrorString(ret) << "\n\n";

      cout << ss.str();

      throw std::logic_error(ss.str());
    }

  return true;
}

bool cuda_check_for_errors_bgn(std::string function_name, std::string file_name, int line)
{
  //cout << __FUNCTION__ << endl;

  //cudaDeviceSynchronize();

  cudaError_t ret = cudaGetLastError();

  if(ret != cudaSuccess)
    {
      std::stringstream ss;

      ss << "\n\n\t error in function (begin) : " << function_name 
	 << "\t file-name : "             << file_name 
	 << "\t at line : "               << line 
	 << "\t error : " << cudaGetErrorString(ret) << "\n\n";

      cout << ss.str();

      throw std::logic_error(ss.str());
    }

  return true;
}

bool cuda_check_for_errors_end(std::string function_name, std::string file_name, int line)
{
  //cout << __FUNCTION__ << endl;

  cudaDeviceSynchronize();

  cudaError_t ret = cudaGetLastError();

  if(ret != cudaSuccess)
    {
      std::stringstream ss;

      ss << "\n\n\t error in function (end) : " << function_name 
	 << "\t file-name : "             << file_name 
	 << "\t at line : "               << line 
	 << "\t error : " << cudaGetErrorString(ret) << "\n\n";

      cout << ss.str();

      throw std::logic_error(ss.str());
    }

  return true;
}

#endif 

/*
cublasHandle_t* cublas_organizer(int thread_id, int N)
{
  static cublasHandle_t* handle = ;
  
  switch(status == 0)
    {
    case 0 :
      {
	if(CUDA_SUCCESS != cublasCreate(handle)){                        
	  fprintf(stderr, "CUBLAS: Not initialized\n"); exit(-1);               
	}    
	
	int* version=NULL;
	if(CUDA_SUCCESS != cublasGetVersion(handle, version)){                        
	  fprintf(stderr, "CUBLAS-verssion: Not initialized\n"); exit(-1);               
	}    
	else
	  cout << "\n\n\t\t CUBLAS-VERSION : " << *version << endl;
      }
      break;

    case 1 :
      {
	if(CUDA_SUCCESS != cublasDestroy(handle)){                        
	  fprintf(stderr, "CUBLAS: Not destroyed\n"); exit(-1);               
	}
	
	handle = NULL;
      }
      break;

    default:
    }
  
  return handle;
};
*/

/*
void cuda_thread_synchronize()
{
  cudaThreadSynchronize();	   
}
*/

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
