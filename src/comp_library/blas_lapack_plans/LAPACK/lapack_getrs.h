//-*-C++-*-

#ifndef LAPACK_GETRS_PLAN_H
#define LAPACK_GETRS_PLAN_H

/*!
 *
 *  \author Peter Staar
 */
  
template<typename scalartype>
void lapack_getrs(char& TRANS, int& N, int& NRHS, scalartype* Matrix_A, int& LDA,
                  int* IPIV, scalartype* Matrix_B, int& LDB, int& INFO)
{
  throw std::logic_error(__PRETTY_FUNCTION__);
}

template<>
void lapack_getrs<float>(char& TRANS, int& N, int& NRHS, float* Matrix_A, int& LDA,
                         int* IPIV, float* Matrix_B, int& LDB, int& INFO)
{
  LAPACK::sgetrs_(&TRANS, &N, &NRHS, Matrix_A, &LDA, IPIV, Matrix_B, &LDB, &INFO);

  if(INFO!=0)
    throw std::logic_error(__PRETTY_FUNCTION__); 
}

template<>
void lapack_getrs<double>(char& TRANS, int& N, int& NRHS, double* Matrix_A, int& LDA,
                          int* IPIV, double* Matrix_B, int& LDB, int& INFO)
{
  LAPACK::dgetrs_(&TRANS, &N, &NRHS, Matrix_A, &LDA, IPIV, Matrix_B, &LDB, &INFO);

  if(INFO!=0)
    throw std::logic_error(__PRETTY_FUNCTION__);
}

template<>
void lapack_getrs<std::complex<float> >(char& TRANS, int& N, int& NRHS, std::complex<float>* Matrix_A, int& LDA,
                                        int* IPIV, std::complex<float>* Matrix_B, int& LDB, int& INFO)
{
  LAPACK::cgetrs_(&TRANS, &N, &NRHS, Matrix_A, &LDA, IPIV, Matrix_B, &LDB, &INFO);
  
  if(INFO!=0)
    throw std::logic_error(__PRETTY_FUNCTION__);
}

template<>
void lapack_getrs<std::complex<double> >(char& TRANS, int& N, int& NRHS, std::complex<double>* Matrix_A, int& LDA,
                                         int* IPIV, std::complex<double>* Matrix_B, int& LDB, int& INFO)
{
  LAPACK::zgetrs_(&TRANS, &N, &NRHS, Matrix_A, &LDA, IPIV, Matrix_B, &LDB, &INFO);
  
  if(INFO!=0)
    throw std::logic_error(__PRETTY_FUNCTION__); 
}

#endif
