// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class transforms domains.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_TRANSFORM_DOMAIN_HPP
#define DCA_MATH_FUNCTION_TRANSFORM_TRANSFORM_DOMAIN_HPP

#include <complex>
#include <iostream>

#include <fftw3.h>

#include "dca/function/function.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/math/function_transform/basis_transform/basis_transform.hpp"
#include "dca/math/function_transform/domain_representations.hpp"
#include "dca/math/function_transform/transform_domain_procedure.hpp"

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

//
// Empty template declaration.
//
template <typename type_input, DOMAIN_REPRESENTATIONS DMN_REP_LHS, typename type_output,
          DOMAIN_REPRESENTATIONS DMN_REP_RHS, int DMN_INDEX>
struct TRANSFORM_DOMAIN {};

//
// Specialization for continous --> continous
//
template <typename type_input, typename type_output, int DMN_INDEX>
class TRANSFORM_DOMAIN<type_input, CONTINUOUS, type_output, CONTINUOUS, DMN_INDEX> {
private:
  const static bool VERBOSE = false;

  typedef basis_transformation<type_input, CONTINUOUS, type_output, CONTINUOUS> contraction_transformation_type;
  typedef typename contraction_transformation_type::matrix_type contraction_matrix_type;

  typedef basis_transformation<type_input, CONTINUOUS, type_output, CONTINUOUS> basis_transformation_type;
  typedef typename basis_transformation_type::matrix_type matrix_type;

public:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute(const func::function<scalartype_input, domain_input>& f_input,
                      func::function<scalartype_output, domain_output>& f_output);

  template <typename scalartype, class domain_input, class domain_output>
  static void execute(const func::function<scalartype, domain_input>& f_input,
                      func::function<scalartype, domain_output>& f_output);

  template <typename scalartype, class domain_input, class domain_output>
  static void execute(const func::function<std::complex<scalartype>, domain_input>& f_input,
                      func::function<std::complex<scalartype>, domain_output>& f_output);

private:
  template <typename f_input_t, typename f_output_t>
  static void characterize_transformation(const f_input_t& f_input, const f_output_t& f_output,
                                          int& M, int& K, int& N, int& P);
};

template <typename type_input, typename type_output, int DMN_INDEX>
template <typename scalartype, class domain_input, class domain_output>
void TRANSFORM_DOMAIN<type_input, CONTINUOUS, type_output, CONTINUOUS, DMN_INDEX>::execute(
    const func::function<scalartype, domain_input>& f_input,
    func::function<scalartype, domain_output>& f_output) {
  if (VERBOSE)
    std::cout << "\n\t transform (continuous -> continuous) " << DMN_INDEX << "  "
              << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";

  int M, K, N, P;
  characterize_transformation(f_input, f_output, M, K, N, P);

  matrix_type& T = basis_transformation_type::get_transformation_matrix();

  if (VERBOSE) {
    std::cout << "\n\t M : " << M << ", K : " << K << ", N : " << N << ", P : " << P << "\n\n";

    T.print();

    f_input.print_fingerprint();
    f_output.print_fingerprint();
  }

  for (int l = 0; l < P; l++) {
    int lin_ind_lhs = M * K * l;
    int lin_ind_rhs = M * N * l;

    dca::linalg::blas::gemm("N", "T", M, N, K, scalartype(1), &f_input(lin_ind_lhs), M, &T(0, 0),
                            T.leadingDimension(), scalartype(0), &f_output(lin_ind_rhs), M);
  }
}

template <typename type_input, typename type_output, int DMN_INDEX>
template <typename scalartype, class domain_input, class domain_output>
void TRANSFORM_DOMAIN<type_input, CONTINUOUS, type_output, CONTINUOUS, DMN_INDEX>::execute(
    const func::function<std::complex<scalartype>, domain_input>& f_input,
    func::function<std::complex<scalartype>, domain_output>& f_output) {
  if (VERBOSE)
    std::cout << "\n\t transform (continuous -> continuous) " << DMN_INDEX << "  "
              << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";

  int M, K, N, P;
  characterize_transformation(f_input, f_output, M, K, N, P);

  contraction_matrix_type& T = basis_transformation_type::get_transformation_matrix();

  if (VERBOSE) {
    std::cout << "\n\t M : " << M << ", K : " << K << ", N : " << N << ", P : " << P << "\n\n";

    T.print_fingerprint();

    f_input.print_fingerprint();
    f_output.print_fingerprint();
  }

  scalartype* A = new scalartype[M * K];
  scalartype* C = new scalartype[M * N];

  for (int l = 0; l < P; l++) {
    int lin_ind_lhs = M * K * l;
    int lin_ind_rhs = M * N * l;

    {  // real
      for (int i = 0; i < M * K; i++)
        A[i] = real(f_input(lin_ind_lhs + i));

      dca::linalg::blas::gemm("N", "T", M, N, K, scalartype(1), A, M, &T(0, 0),
                              T.leadingDimension(), scalartype(0), C, M);

      for (int i = 0; i < M * N; i++)
        real(f_output(lin_ind_rhs + i)) = C[i];
    }

    {  // imag
      for (int i = 0; i < M * K; i++)
        A[i] = imag(f_input(lin_ind_lhs + i));

      dca::linalg::blas::gemm('N', 'T', M, N, K, scalartype(1), A, M, &T(0, 0),
                              T.leadingDimension(), scalartype(0), C, M);

      for (int i = 0; i < M * N; i++)
        imag(f_output(lin_ind_rhs + i)) = C[i];
    }
  }

  delete[] A;
  delete[] C;
}

template <typename type_input, typename type_output, int DMN_INDEX>
template <typename f_input_t, typename f_output_t>
void TRANSFORM_DOMAIN<type_input, CONTINUOUS, type_output, CONTINUOUS, DMN_INDEX>::characterize_transformation(
    const f_input_t& f_input, const f_output_t& f_output, int& M, int& K, int& N, int& P) {
  M = 1;
  for (int l = 0; l < DMN_INDEX; l++)
    M *= f_input[l];

  K = f_input[DMN_INDEX];
  N = f_output[DMN_INDEX];

  P = 1;
  for (int l = DMN_INDEX + 1; l < f_input.signature(); l++)
    P *= f_input[l];
}

//
// Specialization for continous --> expansion
//
template <typename type_input, typename type_output, int DMN_INDEX>
class TRANSFORM_DOMAIN<type_input, CONTINUOUS, type_output, EXPANSION, DMN_INDEX>
    : public TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX> {
private:
  const static bool VERBOSE = false;

  typedef basis_transformation<type_input, CONTINUOUS, type_output, EXPANSION> basis_transformation_type;
  typedef typename basis_transformation_type::matrix_type matrix_type;

public:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute(const func::function<scalartype_input, domain_input>& f_input,
                      func::function<scalartype_output, domain_output>& f_output) {
    default_execute(f_input, f_output);
  }

private:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void default_execute(const func::function<scalartype_input, domain_input>& f_input,
                              func::function<scalartype_output, domain_output>& f_output) {
    if (VERBOSE)
      std::cout << "\n\t default-transform (continuous -> expansion) " << DMN_INDEX << "  "
                << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";

    matrix_type& T = basis_transformation_type::get_transformation_matrix();

    TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(f_input, f_output, T);
  }
};

//
// Specialization for discrete --> expansion
//
template <typename type_input, typename type_output, int DMN_INDEX>
class TRANSFORM_DOMAIN<type_input, DISCRETE, type_output, EXPANSION, DMN_INDEX>
    : public TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX> {
private:
  const static bool VERBOSE = false;

  typedef typename type_input::dmn_specifications_type input_specs_type;
  typedef typename type_output::dmn_specifications_type output_specs_type;

  typedef basis_transformation<type_input, DISCRETE, type_output, EXPANSION> basis_transformation_type;
  typedef typename basis_transformation_type::matrix_type matrix_type;

public:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute(const func::function<scalartype_input, domain_input>& f_input,
                      func::function<scalartype_output, domain_output>& f_output) {
    default_execute(f_input, f_output);
  }

private:
  template <typename scalartype>
  static scalartype vector_norm(scalartype x) {
    return abs(x);
  }

  template <typename scalartype>
  static scalartype vector_norm(std::vector<scalartype>& x) {
    scalartype result = 0;

    for (int l = 0; l < x.size(); l++)
      result += x[l] * x[l];

    return result;
  }

  static int find_origin() {
    int index = 0;

    for (int l = 0; l < type_output::get_size(); l++)
      if (vector_norm(type_output::get_elements()[l]) < 1.e-6)
        index = l;

    std::cout << index << "\n";

    return index;
  }

  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void fftw_harmonics_execute(const func::function<scalartype_input, domain_input>& f_input,
                                     func::function<scalartype_output, domain_output>& f_output,
                                     const bool renorm);

  template <typename scalartype, class domain_input, class domain_output>
  static void fftw_harmonics_execute(
      const func::function<std::complex<scalartype>, domain_input>& f_input,
      func::function<std::complex<scalartype>, domain_output>& f_output, const bool renorm) {
    assert(type_input::dmn_specifications_type::DIMENSION ==
           type_output::dmn_specifications_type::DIMENSION);

    if (VERBOSE)
      std::cout << "\n\t ifftw-harmonics-transform (discrete -> expansion) " << DMN_INDEX << "  "
                << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";

    int M, K, N, P;
    characterize_transformation(f_input, f_output, M, K, N, P);

    int rank = type_input::dmn_specifications_type::DIMENSION;
    int* dims = type_input::get_dimensions();

    int how_many = M;

    fftw_complex* in = new fftw_complex[M * K];
    fftw_complex* out = new fftw_complex[M * N];

    int istride = 1;
    int ostride = 1;

    // K=N !
    int idist = K;
    int odist = N;

    const int* inembed = type_input::get_dimensions();
    const int* onembed = type_output::get_dimensions();

    if (false and VERBOSE) {
      std::cout << M << "\t" << K << "\t" << N << "\t" << P << "\n";

      std::cout << rank << "\n";
      for (int i = 0; i < rank; i++)
        std::cout << dims[i] << "\t";
      std::cout << "\n";

      f_input.print_fingerprint();
      f_output.print_fingerprint();
    }

    fftw_plan plan = fftw_plan_many_dft(rank, dims, how_many, in, inembed, istride, idist, out,
                                        onembed, ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE);

    int index = find_origin();
    for (int l = 0; l < P; l++) {
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < K; j++) {
          in[j + K * i][0] = real(f_input(M * K * l + i + j * M));
          in[j + K * i][1] = imag(f_input(M * K * l + i + j * M));
        }
      }

      fftw_execute(plan);

      for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
          int j_0 = (j - index) < 0 ? (j - index + N) : (j - index);

          real(f_output(M * N * l + i + j * M)) = out[j_0 + N * i][0];
          imag(f_output(M * N * l + i + j * M)) = out[j_0 + N * i][1];
        }
      }
    }

    fftw_destroy_plan(plan);

    delete[] in;
    delete[] out;

    if (renorm)
      f_output *= 1. / scalartype(K);
  }

  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void default_execute(const func::function<scalartype_input, domain_input>& f_input,
                              func::function<scalartype_output, domain_output>& f_output) {
    if (VERBOSE)
      std::cout << "\n default-transform (discrete -> expansion) " << DMN_INDEX << "  "
                << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";

    matrix_type& T = basis_transformation_type::get_transformation_matrix();

    // TODO make better
    if constexpr (std::is_same<float, scalartype_input>::value &&
                  std::is_same<float, scalartype_output>::value) {
      func::function<double, domain_input> f_copy;
      f_copy = f_input;
      func::function<double, domain_output> f_out_cpy;
      TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(f_copy, f_out_cpy, T);
      f_output = f_out_cpy;
    }
    else {
      TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(f_input, f_output, T);
    }
  }
};

//
// Specialization for expansion --> continous
//
template <typename type_input, typename type_output, int DMN_INDEX>
class TRANSFORM_DOMAIN<type_input, EXPANSION, type_output, CONTINUOUS, DMN_INDEX>
    : public TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX> {
private:
  const static bool VERBOSE = false;

  typedef basis_transformation<type_input, EXPANSION, type_output, CONTINUOUS> basis_transformation_type;
  typedef typename basis_transformation_type::matrix_type matrix_type;

public:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute(const func::function<scalartype_input, domain_input>& f_input,
                      func::function<scalartype_output, domain_output>& f_output) {
    default_execute(f_input, f_output);
  }

private:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void default_execute(const func::function<scalartype_input, domain_input>& f_input,
                              func::function<scalartype_output, domain_output>& f_output) {
    if (VERBOSE)
      std::cout << "\n\t default-transform (expansion -> continuous) " << DMN_INDEX << "  "
                << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";

    matrix_type& T = basis_transformation_type::get_transformation_matrix();

    TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(f_input, f_output, T);
  }
};

//
// Specialization for expansion --> discrete
//
template <typename type_input, typename type_output, int DMN_INDEX>
class TRANSFORM_DOMAIN<type_input, EXPANSION, type_output, DISCRETE, DMN_INDEX>
    : public TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX> {
private:
  const static bool VERBOSE = false;

  typedef typename type_input::dmn_specifications_type input_specs_type;
  typedef typename type_output::dmn_specifications_type output_specs_type;

  typedef basis_transformation<type_input, EXPANSION, type_output, DISCRETE> basis_transformation_type;
  typedef typename basis_transformation_type::matrix_type matrix_type;

public:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute(const func::function<scalartype_input, domain_input>& f_input,
                      func::function<scalartype_output, domain_output>& f_output) {
    default_execute(f_input, f_output);
  }

private:
  template <typename scalartype>
  static scalartype vector_norm(scalartype x) {
    return abs(x);
  }

  template <typename scalartype>
  static scalartype vector_norm(std::vector<scalartype>& x) {
    scalartype result = 0;

    for (int l = 0; l < x.size(); l++)
      result += x[l] * x[l];

    return result;
  }

  static int find_origin() {
    int index = 0;

    for (int l = 0; l < type_input::get_size(); l++)
      if (vector_norm(type_input::get_elements()[l]) < 1.e-6)
        index = l;

    return index;
  }

  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void fftw_harmonics_execute(const func::function<scalartype_input, domain_input>& f_input,
                                     func::function<scalartype_output, domain_output>& f_output);

  template <typename scalartype, class domain_input, class domain_output>
  static void fftw_harmonics_execute(
      const func::function<std::complex<scalartype>, domain_input>& f_input,
      func::function<std::complex<scalartype>, domain_output>& f_output) {
    assert(type_input::dmn_specifications_type::DIMENSION ==
           type_output::dmn_specifications_type::DIMENSION);

    if (VERBOSE)
      std::cout << "\n\t fftw-harmonics-transform (expansion -> discrete) " << DMN_INDEX << "  "
                << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";

    int M, K, N, P;
    characterize_transformation(f_input, f_output, M, K, N, P);

    int rank = type_input::dmn_specifications_type::DIMENSION;
    int* dims = type_input::get_dimensions();

    int how_many = M;

    fftw_complex* in = new fftw_complex[M * K];
    fftw_complex* out = new fftw_complex[M * N];

    int istride = 1;
    int ostride = 1;

    // K=N !
    int idist = K;
    int odist = N;

    const int* inembed = type_input::get_dimensions();
    const int* onembed = type_output::get_dimensions();

    if (VERBOSE) {
      std::cout << M << "\t" << K << "\t" << N << "\t" << P << "\n";

      f_input.print_fingerprint();
      f_output.print_fingerprint();
    }

    fftw_plan plan = fftw_plan_many_dft(rank, dims, how_many, in, inembed, istride, idist, out,
                                        onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);

    int index = find_origin();
    for (int l = 0; l < P; l++) {
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < K; j++) {
          int j_0 = (j - index) < 0 ? (j - index + K) : (j - index);

          in[j_0 + K * i][0] = real(f_input(M * K * l + i + j * M));
          in[j_0 + K * i][1] = imag(f_input(M * K * l + i + j * M));
        }
      }

      fftw_execute(plan);

      for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
          real(f_output(M * N * l + i + j * M)) = out[j + N * i][0];
          imag(f_output(M * N * l + i + j * M)) = out[j + N * i][1];
        }
      }
    }

    fftw_destroy_plan(plan);

    delete[] in;
    delete[] out;
  }

  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void fftw_cosine_execute(const func::function<scalartype_input, domain_input>& /*f_input*/,
                                  func::function<scalartype_output, domain_output>& /*f_output*/) {
    if (VERBOSE)
      std::cout << "\n\t fftw-cosine-transform (expansion -> discrete) " << DMN_INDEX << "  "
                << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";
  }

  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void fftw_sine_execute(const func::function<scalartype_input, domain_input>& /*f_input*/,
                                func::function<scalartype_output, domain_output>& /*f_output*/) {
    if (VERBOSE)
      std::cout << "\n\t fftw-sine-transform (expansion -> discrete) " << DMN_INDEX << "  "
                << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";
  }

  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void default_execute(const func::function<scalartype_input, domain_input>& f_input,
                              func::function<scalartype_output, domain_output>& f_output) {
    if (VERBOSE)
      std::cout << "\n\t default-transform (expansion -> discrete) " << DMN_INDEX << "  "
                << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";

    matrix_type& T = basis_transformation_type::get_transformation_matrix();

    TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(f_input, f_output, T);
  }
};

//
// Specialization for expansion --> expansion
//
template <typename type_input, typename type_output, int DMN_INDEX>
class TRANSFORM_DOMAIN<type_input, EXPANSION, type_output, EXPANSION, DMN_INDEX>
    : public TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX> {
private:
  const static bool VERBOSE = false;

  typedef basis_transformation<type_input, EXPANSION, type_output, EXPANSION> basis_transformation_type;
  typedef typename basis_transformation_type::matrix_type matrix_type;

public:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute(const func::function<scalartype_input, domain_input>& f_input,
                      func::function<scalartype_output, domain_output>& f_output) {
    default_execute(f_input, f_output);
  }

private:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void default_execute(const func::function<scalartype_input, domain_input>& f_input,
                              func::function<scalartype_output, domain_output>& f_output) {
    if (VERBOSE)
      std::cout << "\n\t default-transform (continuous -> expansion) " << DMN_INDEX << "  "
                << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";

    matrix_type& T = basis_transformation_type::get_transformation_matrix();

    TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(f_input, f_output, T);
  }
};

}  // namespace transform
}  // namespace math
}  // namespace dca

#endif  // DCA_MATH_FUNCTION_TRANSFORM_TRANSFORM_DOMAIN_HPP
