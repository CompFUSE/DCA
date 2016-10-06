// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_DOMAIN_TRANSFORMS_TRANSFORM_DOMAIN_ED_TO_DD_H
#define MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_DOMAIN_TRANSFORMS_TRANSFORM_DOMAIN_ED_TO_DD_H

#include <iostream>
#include <fftw3.h>

#include "dca/function/function.hpp"
#include "comp_library/linalg/linalg.hpp"
#include "math_library/functional_transforms/basis_transforms/basis_transforms.hpp"
#include "math_library/functional_transforms/domain_transforms/transform_domain_template.h"
#include "math_library/functional_transforms/domain_transforms/transformation_characteristics.h"
#include "math_library/functional_transforms/typedefs.hpp"

namespace math_algorithms {
namespace functional_transforms {
// math_algorithms::functional_transforms::

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
  static void execute(func::function<scalartype_input, domain_input>& f_input,
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
  static void fftw_harmonics_execute(func::function<scalartype_input, domain_input>& f_input,
                                     func::function<scalartype_output, domain_output>& f_output);

  template <typename scalartype, class domain_input, class domain_output>
  static void fftw_harmonics_execute(func::function<std::complex<scalartype>, domain_input>& f_input,
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
  static void fftw_cosine_execute(func::function<scalartype_input, domain_input>& /*f_input*/,
                                  func::function<scalartype_output, domain_output>& /*f_output*/) {
    if (VERBOSE)
      std::cout << "\n\t fftw-cosine-transform (expansion -> discrete) " << DMN_INDEX << "  "
                << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";
  }

  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void fftw_sine_execute(func::function<scalartype_input, domain_input>& /*f_input*/,
                                func::function<scalartype_output, domain_output>& /*f_output*/) {
    if (VERBOSE)
      std::cout << "\n\t fftw-sine-transform (expansion -> discrete) " << DMN_INDEX << "  "
                << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";
  }

  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void default_execute(func::function<scalartype_input, domain_input>& f_input,
                              func::function<scalartype_output, domain_output>& f_output) {
    if (VERBOSE)
      std::cout << "\n\t default-transform (expansion -> discrete) " << DMN_INDEX << "  "
                << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";

    matrix_type& T = basis_transformation_type::get_transformation_matrix();

    TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(f_input, f_output, T);
  }
};

}  // functional_transforms
}  // math_algorithms

#endif  // MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_DOMAIN_TRANSFORMS_TRANSFORM_DOMAIN_ED_TO_DD_H
