// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Specialization for expansion --> continous.

template <typename input_type, typename output_type>
class basis_transformation<input_type, EXPANSION, output_type, CONTINUOUS> {
public:
  typedef input_type rh_dmn_type;
  typedef output_type lh_dmn_type;

  typedef typename lh_dmn_type::dmn_specifications_type lh_spec_dmn_type;
  typedef typename rh_dmn_type::dmn_specifications_type rh_spec_dmn_type;

  typedef typename lh_spec_dmn_type::scalar_type lh_scalar_type;
  typedef typename rh_spec_dmn_type::scalar_type rh_scalar_type;

  typedef typename lh_spec_dmn_type::element_type lh_element_type;
  typedef typename rh_spec_dmn_type::element_type rh_element_type;

  typedef basis_function<lh_dmn_type, lh_spec_dmn_type::BASIS_EXPANSION, rh_dmn_type,
                         rh_spec_dmn_type::BASIS_EXPANSION>
      basis_function_type;

  typedef typename basis_function_type::f_scalar_type f_scalar_type;

  typedef dca::linalg::Matrix<f_scalar_type, dca::linalg::CPU> matrix_type;

public:
  static bool& is_initialized() {
    static bool initialized = false;
    return initialized;
  }

  static std::string& get_name() {
    static std::string name = "basis-transformation";
    return name;
  }

  static dca::linalg::Matrix<f_scalar_type, dca::linalg::CPU>& get_transformation_matrix() {
    static dca::linalg::Matrix<f_scalar_type, dca::linalg::CPU> T;

    if (not is_initialized())
      initialize_transformation_matrix();

    return T;
  }

  static void initialize_transformation_matrix() {
    is_initialized() = true;

    int M = lh_dmn_type::get_size();
    int N = rh_dmn_type::get_size();

    assert(M > 0 and N > 0);

    dca::linalg::Matrix<f_scalar_type, dca::linalg::CPU>& T = get_transformation_matrix();

    T.resizeNoCopy(std::pair<int, int>(M, N));

    for (int j = 0; j < N; j++)
      for (int i = 0; i < M; i++)
        T(i, j) = basis_function_type::execute(i, j);
  }
};
