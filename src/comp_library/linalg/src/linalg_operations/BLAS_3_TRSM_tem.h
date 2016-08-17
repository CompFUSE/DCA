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
        "L", &uplo, "N", &diag, X.nrRows(), X.nrCols(), scalartype(1),
        A.ptr(), A.leadingDimension(), X.ptr(), X.leadingDimension(), thread_id,
        stream_id);
  }
};
}

#endif
