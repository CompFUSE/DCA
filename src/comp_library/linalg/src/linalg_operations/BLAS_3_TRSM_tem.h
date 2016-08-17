//-*-C++-*-

#ifndef LINALG_TRSM_TEM_H
#define LINALG_TRSM_TEM_H

namespace LIN_ALG {

template <device_type device_name>
class TRSM {
public:
  template <typename scalartype>
  inline static void execute(char uplo, char diag, matrix<scalartype, device_name>& A,
                             matrix<scalartype, device_name>& X, int thread_id, int stream_id) {
    assert(uplo == 'U' or uplo == 'L');
    assert(diag == 'U' or diag == 'N');

    dca::linalg::blas::UseDevice<device_name>::trsm(
        "L", &uplo, "N", &diag, X.get_number_of_rows(), X.get_number_of_cols(), scalartype(1),
        A.get_ptr(), A.get_leading_dimension(), X.get_ptr(), X.get_leading_dimension(), thread_id,
        stream_id);
  }
};
}

#endif
