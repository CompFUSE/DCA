//-*-C++-*-

#ifndef COPY_SCALARTYPES_IN_FUNCTIONS_H
#define COPY_SCALARTYPES_IN_FUNCTIONS_H

/*! \class  
 *  \author Peter Staar
 *  \date 2011
 *  \version 1.0
 *  \brief This class implements an efficient copy mechanism. If it doesnt know how to copy the type, it will give a compilation error.

 */
template<class whatever_t>
struct copy_from
{};

template<>
struct copy_from<short>
{
public:
  
  static void execute(int size, short* whatever_l, short* whatever_r){
    memcpy(whatever_l, whatever_r, sizeof(short)*size);
  }
};

template<>
struct copy_from<int>
{
public:
  
  static void execute(int size, int* whatever_l, int* whatever_r){
    memcpy(whatever_l, whatever_r, sizeof(int)*size);
  }
};

template<>
struct copy_from<long>
{
public:
  
  static void execute(int size, long* whatever_l, long* whatever_r){
    memcpy(whatever_l, whatever_r, sizeof(long)*size);
  }
};

template<>
struct copy_from<size_t>
{
public:
  
  static void execute(int size, size_t* whatever_l, size_t* whatever_r){
    memcpy(whatever_l, whatever_r, sizeof(size_t)*size);
  }
};

template<>
struct copy_from<float>
{
public:
  
  static void execute(int size, float* whatever_l, float* whatever_r){
    memcpy(whatever_l, whatever_r, sizeof(float)*size);
  }
};

template<>
struct copy_from<double>
{
public:
  
  static void execute(int size, double* whatever_l, double* whatever_r){
    memcpy(whatever_l, whatever_r, sizeof(double)*size);
  }
};

template<>
struct copy_from<std::complex<float> >
{
public:
  
  static void execute(int size, std::complex<float>* whatever_l, std::complex<float>* whatever_r){
    memcpy(whatever_l, whatever_r, sizeof(std::complex<float>)*size);
  }
};

template<>
struct copy_from<std::complex<double> >
{
public:
  
  static void execute(int size, std::complex<double>* whatever_l, std::complex<double>* whatever_r){
    memcpy(whatever_l, whatever_r, sizeof(std::complex<double>)*size);
  }
};


template<class whatever_t>
struct copy_from<std::vector<whatever_t> >
{
public:

  static void execute(int size, std::vector<whatever_t>* whatever_l, std::vector<whatever_t>* whatever_r){
    for(int l=0; l<size; l++){
      whatever_l[l] = whatever_r[l];
    }
  }
};

template<>
struct copy_from<fftw_complex>
{
public:

  static void execute(int size, fftw_complex* whatever_l, fftw_complex* whatever_r){
    memcpy(whatever_l, whatever_r, 2*sizeof(double)*size);
  }
};

#endif
