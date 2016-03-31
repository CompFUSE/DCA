//-*-C++-*-

#ifndef LIN_ALG_DEVICE_MATRICES_H
#define LIN_ALG_DEVICE_MATRICES_H

namespace LIN_ALG {
 
  template<typename scalartype>
  class device_matrix
  {
    typedef typename MATRIX_SCALARTYPE<scalartype, device_name>::new_scalartype matrix_scalartype;
	
  public:

#ifdef __CUDACC__
    __device__ scalartype& operator()(int i, int j)
    {
      assert(i>-1 and i<current_size[0]);
      assert(j>-1 and j<current_size[1]);

      return data[i+j*global_size[0]];
    }
#endif

    int global_size [2];
    int current_size[2];
    
    matrix_scalartype* data;
  };

}

#endif
