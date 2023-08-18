// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Peter W. Doak (doakpw@ornl.gov)
//
// Specialization for discrete --> expansion.
//
// unlike other transforms this has been altered to not be relay on singleton behavior of static variables
//

template <typename input_type, typename output_type>
class basis_transformation<input_type, DISCRETE, output_type, EXPANSION> {
public:
  typedef input_type rh_dmn_type;
  typedef output_type lh_dmn_type;

  typedef typename lh_dmn_type::dmn_specifications_type lh_spec_dmn_type;
  typedef typename rh_dmn_type::dmn_specifications_type rh_spec_dmn_type;

  typedef typename lh_spec_dmn_type::scalar_type lh_scalar_type;
  typedef typename rh_spec_dmn_type::scalar_type rh_scalar_type;

  typedef typename lh_spec_dmn_type::element_type lh_element_type;
  typedef typename rh_spec_dmn_type::element_type rh_element_type;

  typedef basis_transformation<output_type, EXPANSION, input_type, DISCRETE> inverse_basis_transformation_type;

  typedef typename inverse_basis_transformation_type::matrix_type matrix_type;

  // typedef typename lh_spec_dmn_type::scalar_type f_scalar_type;
  // typedef dca::linalg::Matrix<f_scalar_type, dca::linalg::CPU> matrix_type;

public:

  static std::string get_name() {
    return {"basis-transformation"};
  }

  static matrix_type get_transformation_matrix() {
    matrix_type t_matrix;
    initialize_transformation_matrix(t_matrix);
    return t_matrix;
  }

  static void initialize_transformation_matrix(matrix_type& t_matrix) {
    int M = lh_dmn_type::get_size();
    int N = rh_dmn_type::get_size();

    assert(M > 0 and N > 0);

    t_matrix.resizeNoCopy(std::pair<int, int>(M, N));

    matrix_type inv_t_matrix = inverse_basis_transformation_type::get_transformation_matrix();

    dca::linalg::matrixop::pseudoInverse(inv_t_matrix, t_matrix);
  }
};
