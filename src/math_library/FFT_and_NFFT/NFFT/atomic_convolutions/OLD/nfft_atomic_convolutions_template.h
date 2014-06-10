//-*-C++-*-

#ifndef MATH_ALGORITHMS_NFFT_ATOMIC_CONVOLUTIONS_1D_H
#define MATH_ALGORITHMS_NFFT_ATOMIC_CONVOLUTIONS_1D_H

namespace MATH_ALGORITHMS
{
  namespace NFFT
  {

    template<int oversampling, int step=1>
    struct nfft_atomic_convolution
    {
      template<typename scalar_type>
      inline static void execute(scalar_type* f, scalar_type tau, scalar_type val, scalar_type* der_phi, scalar_type* phi)
      {
        nfft_atomic_convolution_single_op<2*oversampling+1, 0, step>::execute(f, tau, val, der_phi, phi);
      }

      template<typename scalar_type>
      inline static void execute(scalar_type* f_tau, scalar_type* tau, scalar_type val, scalar_type* alpha, scalar_type* beta, scalar_type* gamma, scalar_type* delta)
      {
        nfft_atomic_convolution_single_op<2*oversampling+1, 0, step>::execute(f_tau, tau, val, alpha, beta, gamma, delta);
      }
    };














    /*
      template<int oversampling, int step=1>
      struct atomic_convolution_square
      {
      template<typename scalar_type>
      inline static void execute(scalar_type* f, scalar_type tau, scalar_type val, scalar_type* der_phi, scalar_type* phi)
      {
      atomic_convolution_single_op_square<2*oversampling+1, 0, step>::execute(f, tau, val, der_phi, phi);
      }

      template<typename scalar_type>
      inline static void execute(scalar_type* f, scalar_type* tau, scalar_type val, scalar_type* alpha, scalar_type* beta, scalar_type* gamma, scalar_type* delta)
      {
      atomic_convolution_single_op_square<2*oversampling+1, 0, step>::execute(f, tau, val, alpha, beta, gamma, delta);
      }
      };
    */
  }

}

#endif
