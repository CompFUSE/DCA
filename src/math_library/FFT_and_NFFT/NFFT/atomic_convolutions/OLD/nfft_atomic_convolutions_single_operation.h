//-*-C++-*-

#ifndef MATH_ALGORITHMS_NFFT_ATOMIC_CONVOLUTIONS_1D_SINGLE_OPERATION_H
#define MATH_ALGORITHMS_NFFT_ATOMIC_CONVOLUTIONS_1D_SINGLE_OPERATION_H

namespace MATH_ALGORITHMS
{
  namespace NFFT
  {

    template<int max_count, int count, int step>
    struct nfft_atomic_convolution_single_op
    {
      template<typename scalar_type>
      inline static void execute(scalar_type* f, scalar_type tau, scalar_type val, scalar_type* der_phi, scalar_type* phi)
      {
        f[count] += val*(der_phi[count*step]*tau+phi[count*step]);
        nfft_atomic_convolution_single_op<max_count, count+1,step>::execute(f, tau, val, der_phi, phi);
      }

      template<typename scalar_type>
      inline static void execute(scalar_type* f, scalar_type* tau, scalar_type val, scalar_type* alpha, scalar_type* beta, scalar_type* gamma, scalar_type* delta)
      {
        f[count] += val*( alpha[count*step]
                          +beta [count*step]*tau[0]
                          +gamma[count*step]*tau[1]
                          +delta[count*step]*tau[2]);

        nfft_atomic_convolution_single_op<max_count, count+1, step>::execute(f, tau, val, alpha, beta, gamma, delta);
      }
    };

    template<int max_count, int step>
    struct nfft_atomic_convolution_single_op<max_count, max_count, step>
    {
      template<typename scalar_type>
      inline static void execute(scalar_type* f, scalar_type tau, scalar_type val, scalar_type* der_phi, scalar_type* phi)
      {}

      template<typename scalar_type>
      inline static void execute(scalar_type* f, scalar_type* tau, scalar_type val, scalar_type* alpha, scalar_type* beta, scalar_type* gamma, scalar_type* delta)
      {}
    };








    /*
      template<int max_count, int count, int step>
      struct atomic_convolution_single_op_square
      {
      template<typename scalar_type>
      inline static void execute(scalar_type* f, scalar_type tau, scalar_type val, scalar_type* der_phi, scalar_type* phi)
      {
      scalar_type new_val =  val*(der_phi[count*step]*tau+phi[count*step]);

      f[count] += (new_val*new_val);

      atomic_convolution_single_op_square<max_count, count+1,step>::execute(f, tau, val, der_phi, phi);
      }

      template<typename scalar_type>
      inline static void execute(scalar_type* f, scalar_type* tau, scalar_type val, scalar_type* alpha, scalar_type* beta, scalar_type* gamma, scalar_type* delta)
      {
      scalar_type new_val = val*( alpha[count*step]
      +beta [count*step]*tau[0]
      +gamma[count*step]*tau[1]
      +delta[count*step]*tau[2]);

      f[count] += (new_val*new_val);

      atomic_convolution_single_op_square<max_count, count+1, step>::execute(f, tau, val, alpha, beta, gamma, delta);
      }
      };

      template<int max_count, int step>
      struct atomic_convolution_single_op_square<max_count, max_count, step>
      {
      template<typename scalar_type>
      inline static void execute(scalar_type* f, scalar_type tau, scalar_type val, scalar_type* der_phi, scalar_type* phi)
      {}

      template<typename scalar_type>
      inline static void execute(scalar_type* f, scalar_type* tau, scalar_type val, scalar_type* alpha, scalar_type* beta, scalar_type* gamma, scalar_type* delta)
      {}
      };
    */
  }

}

#endif
