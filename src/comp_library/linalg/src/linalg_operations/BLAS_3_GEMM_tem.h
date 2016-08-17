
#ifndef LINALG_GEMM_TEM_H
#define LINALG_GEMM_TEM_H

namespace LIN_ALG {

template <device_type device_name>
class GEMM {
public:
  template <typename scalartype>
  static inline void execute(matrix<scalartype, device_name>& a, matrix<scalartype, device_name>& b,
                             matrix<scalartype, device_name>& c, int thread_id = 0,
                             int stream_id = 0) {
    execute<scalartype>('N', 'N', 1., a, b, 0., c, thread_id, stream_id);
  }

  template <typename scalartype>
  static inline void execute(scalartype alpha, matrix<scalartype, device_name>& a,
                             matrix<scalartype, device_name>& b, scalartype beta,
                             matrix<scalartype, device_name>& c, int thread_id = 0,
                             int stream_id = 0) {
    execute('N', 'N', alpha, a, b, beta, c, thread_id, stream_id);
  }

  template <typename scalartype>
  static void execute(char transa, char transb, matrix<scalartype, device_name>& a,
                      matrix<scalartype, device_name>& b, matrix<scalartype, device_name>& c,
                      int thread_id = 0, int stream_id = 0) {
    execute<scalartype>(transa, transb, 1., a, b, 0., c, thread_id, stream_id);
  }

  template <typename scalartype>
  static void execute(char transa, char transb, scalartype alpha, matrix<scalartype, device_name>& a,
                      matrix<scalartype, device_name>& b, scalartype beta,
                      matrix<scalartype, device_name>& c, int thread_id = 0, int stream_id = 0) {
    int m = c.get_current_size().first;
    int n = c.get_current_size().second;
    int k;

    if (transa == 'N') {
      assert(a.get_current_size().first == m);
      k = a.get_current_size().second;
    }
    else {
      assert(a.get_current_size().second == m);
      k = a.get_current_size().first;
    }

    if (transb == 'N') {
      assert(b.get_current_size().first == k);
      assert(b.get_current_size().second == n);
    }
    else {
      assert(b.get_current_size().second == k);
      assert(b.get_current_size().first == n);
    }

    int lda = a.get_global_size().first;
    int ldb = b.get_global_size().first;
    int ldc = c.get_global_size().first;

    dca::linalg::blas::UseDevice<device_name>::gemm(&transa, &transb, m, n, k, alpha, a.get_ptr(),
                                                    lda, b.get_ptr(), ldb, beta, c.get_ptr(), ldc,
                                                    thread_id, stream_id);
  }

  /*
      template<typename scalartype>
      static void execute(char TRANSA, char TRANSB,
                          matrix<             scalartype , device_name>& A,
                          matrix<std::complex<scalartype>, device_name>& B,
                          matrix<std::complex<scalartype>, device_name>& C)
      {
        assert(test_sizes(TRANSA, TRANSB, A, B, C));

        matrix<scalartype, device_name> B_re("B_re", B.get_current_size(), B.get_global_size());
        matrix<scalartype, device_name> B_im("B_im", B.get_current_size(), B.get_global_size());

        for(int j=0; j<B.get_current_size().second; j++)
          for(int i=0; i<B.get_current_size().first; i++)
            B_re(i,j) = real(B(i,j));

        for(int j=0; j<B.get_current_size().second; j++)
          for(int i=0; i<B.get_current_size().first; i++)
            B_im(i,j) = imag(B(i,j));

        matrix<scalartype, device_name> C_re("C_re", C.get_current_size(), C.get_global_size());
        matrix<scalartype, device_name> C_im("C_im", C.get_current_size(), C.get_global_size());

        execute(TRANSA, TRANSB, A, B_re, C_re, 0, 0);
        execute(TRANSA, TRANSB, A, B_im, C_im, 0, 0);

        std::complex<scalartype> I(0,1);
        for(int j=0; j<C.get_current_size().second; j++)
          for(int i=0; i<C.get_current_size().first; i++)
            C(i,j) = C_re(i,j)+I*C_im(i,j);
      }

      //template<typename scalartype>
      static void execute(char TRANSA, char TRANSB,
                          matrix<std::complex<double>, device_name>& A,
                          matrix<             double , device_name>& B,
                          matrix<std::complex<double>, device_name>& C)
      {
        assert(test_sizes(TRANSA, TRANSB, A, B, C));

        matrix<double, device_name> A_re("A_re", A.get_current_size(), A.get_global_size());
        matrix<double, device_name> A_im("A_im", A.get_current_size(), A.get_global_size());

        for(int j=0; j<A.get_current_size().second; j++)
          for(int i=0; i<A.get_current_size().first; i++)
            A_re(i,j) = real(A(i,j));

        for(int j=0; j<A.get_current_size().second; j++)
          for(int i=0; i<A.get_current_size().first; i++)
            A_im(i,j) = imag(A(i,j));

        matrix<double, device_name> C_re("C_re", C.get_current_size(), C.get_global_size());
        matrix<double, device_name> C_im("C_im", C.get_current_size(), C.get_global_size());

        execute(TRANSA, TRANSB, A_re, B, C_re, 0, 0);
        execute(TRANSA, TRANSB, A_im, B, C_im, 0, 0);

        std::complex<double> I(0,1);
        for(int j=0; j<C.get_current_size().second; j++)
          for(int i=0; i<C.get_current_size().first; i++)
            C(i,j) = C_re(i,j)+I*C_im(i,j);
      }
  */
};
}

#endif
