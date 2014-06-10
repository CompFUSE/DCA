//-*-C++-*-

#ifndef DCA_QMCI_N_MATRIX_ROUTINES_GPU_H
#define DCA_QMCI_N_MATRIX_ROUTINES_GPU_H

namespace DCA
{
  namespace QMCI
  {
    namespace N_MATRIX_TOOLS_GPU_KERNELS {

      void compute_G_cols(int N_i, int N_r, int N_c,
                          int*    p_ptr,
                          double* exp_V_ptr,
                          double* N_ptr,      int N_ld,
                          double* G_ptr,      int G_ld,
                          double* G_cols_ptr, int G_cols_ld,
                          int thread_id, int stream_id);

      void compute_d_vector(int N_i, int* d_ind, double* d_ptr, int* p_ptr, double* N_ptr, int N_ld,
                            int thread_id, int stream_id);
    }

    template<typename parameters_type>
    class N_MATRIX_TOOLS<LIN_ALG::GPU, parameters_type>
    {
      const static int MAX_VERTEX_SINGLETS=4;

      typedef typename parameters_type::concurrency_type concurrency_type;
      typedef typename parameters_type::profiler_type    profiler_t;

    public :

      N_MATRIX_TOOLS(int              id,
                     parameters_type& parameters_ref);

      ~N_MATRIX_TOOLS();

      int* get_permutation();
      void set_permutation(std::vector<int>& p);

      void set_d_vector(std::vector<int>&                      d_index,
                        LIN_ALG::matrix<double, LIN_ALG::GPU>& N,
                        LIN_ALG::vector<double, LIN_ALG::CPU>& d_inv);

      void set_d_vector(LIN_ALG::vector<double, LIN_ALG::CPU>& d_inv);

      void scale_rows(LIN_ALG::matrix<double, LIN_ALG::GPU>& N);


      double* get_device_ptr(LIN_ALG::vector<double, LIN_ALG::CPU>& v);

      void copy_rows(LIN_ALG::matrix<double, LIN_ALG::GPU>& N,
                     LIN_ALG::matrix<double, LIN_ALG::GPU>& N_new_spins);

      void compute_G_cols(std::vector<double>&                    exp_V,
                          LIN_ALG::matrix<double, LIN_ALG::GPU>& N,
                          LIN_ALG::matrix<double, LIN_ALG::GPU>& G,
                          LIN_ALG::matrix<double, LIN_ALG::GPU>& G_cols);

    private:

      int thread_id;
      int stream_id;

      parameters_type&  parameters;
      concurrency_type& concurrency;

      LIN_ALG::vector<int, LIN_ALG::GPU> identity;
      LIN_ALG::vector<int, LIN_ALG::GPU> permutation;

      LIN_ALG::vector<double, LIN_ALG::GPU> tmp;

      LIN_ALG::vector<double, LIN_ALG::GPU> exp_V;

      LIN_ALG::vector<int   , LIN_ALG::GPU> d_ind;
      LIN_ALG::vector<double, LIN_ALG::GPU> d_vec;
    };

    template<typename parameters_type>
    N_MATRIX_TOOLS<LIN_ALG::GPU, parameters_type>::N_MATRIX_TOOLS(int              id,
                                                                  parameters_type& parameters_ref):
      thread_id(id),
      stream_id(0),

      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      identity   ("identity    N_MATRIX_TOOLS<LIN_ALG::GPU>", MAX_VERTEX_SINGLETS*parameters.get_K_PHANI()),
      permutation("permutation N_MATRIX_TOOLS<LIN_ALG::GPU>", MAX_VERTEX_SINGLETS*parameters.get_K_PHANI()),

      tmp  ("tmp   N_MATRIX_TOOLS<LIN_ALG::GPU>", MAX_VERTEX_SINGLETS*parameters.get_K_PHANI()),
      exp_V("exp_V N_MATRIX_TOOLS<LIN_ALG::GPU>", MAX_VERTEX_SINGLETS*parameters.get_K_PHANI()),

      d_ind("d_ind N_MATRIX_TOOLS<LIN_ALG::GPU>", MAX_VERTEX_SINGLETS*parameters.get_K_PHANI()),
      d_vec("d_vec N_MATRIX_TOOLS<LIN_ALG::GPU>", MAX_VERTEX_SINGLETS*parameters.get_K_PHANI())
    {
      {
        identity   .set_thread_and_stream_id(thread_id, stream_id);
        permutation.set_thread_and_stream_id(thread_id, stream_id);

        tmp.set_thread_and_stream_id(thread_id, stream_id);
        exp_V.set_thread_and_stream_id(thread_id, stream_id);

        d_ind.set_thread_and_stream_id(thread_id, stream_id);
        d_vec.set_thread_and_stream_id(thread_id, stream_id);
      }

      {
        std::vector<int> id_tmp(MAX_VERTEX_SINGLETS*parameters.get_K_PHANI());

        for(int l=0; l<MAX_VERTEX_SINGLETS*parameters.get_K_PHANI(); ++l)
          id_tmp[l] = l;

        identity.set(id_tmp);
      }
    }

    template<typename parameters_type>
    N_MATRIX_TOOLS<LIN_ALG::GPU, parameters_type>::~N_MATRIX_TOOLS()
    {}

    template<typename parameters_type>
    int* N_MATRIX_TOOLS<LIN_ALG::GPU, parameters_type>::get_permutation()
    {
      return permutation.get_ptr();
    }

    template<typename parameters_type>
    void N_MATRIX_TOOLS<LIN_ALG::GPU, parameters_type>::set_permutation(std::vector<int>& p)
    {
      permutation.set(p, LIN_ALG::ASYNCHRONOUS);
    }

    template<typename parameters_type>
    void N_MATRIX_TOOLS<LIN_ALG::GPU, parameters_type>::set_d_vector(LIN_ALG::vector<double, LIN_ALG::CPU>& d_inv)
    {
      d_vec.set(d_inv, LIN_ALG::ASYNCHRONOUS);
    }

    template<typename parameters_type>
    void N_MATRIX_TOOLS<LIN_ALG::GPU, parameters_type>::scale_rows(LIN_ALG::matrix<double, LIN_ALG::GPU>& N)
    {
      assert(permutation.size() == d_vec.size());

      int N_i = permutation.size();
      int N_c = N.get_number_of_cols();

      int N_LD = N.get_leading_dimension();

      LIN_ALG::SCALE<LIN_ALG::GPU>::many_rows(N_c, N_i,
                                              permutation.get_ptr(),
                                              d_vec      .get_ptr(),
                                              N          .get_ptr(), N_LD,
                                              thread_id, stream_id);
    }

    template<typename parameters_type>
    double* N_MATRIX_TOOLS<LIN_ALG::GPU, parameters_type>::get_device_ptr(LIN_ALG::vector<double, LIN_ALG::CPU>& v)
    {
      tmp.set(v, LIN_ALG::ASYNCHRONOUS);

      return tmp.get_ptr();
    }

    template<typename parameters_type>
    void N_MATRIX_TOOLS<LIN_ALG::GPU, parameters_type>::copy_rows(LIN_ALG::matrix<double, LIN_ALG::GPU>& N,
                                                                  LIN_ALG::matrix<double, LIN_ALG::GPU>& N_new_spins)
    {
      assert(N_new_spins.get_number_of_cols() == N.get_number_of_cols());
      assert(N_new_spins.get_number_of_rows() == permutation.size());

      int N_i = permutation.size();
      int N_c = N.get_number_of_cols();

      assert(N_i<=identity.size());

      LIN_ALG::COPY<LIN_ALG::GPU>::many_rows(N_c, N_i,
                                             permutation.get_ptr(), N          .get_ptr(), N          .get_leading_dimension(),
                                             identity   .get_ptr(), N_new_spins.get_ptr(), N_new_spins.get_leading_dimension(),
                                             thread_id, stream_id);
    }

    template<typename parameters_type>
    void N_MATRIX_TOOLS<LIN_ALG::GPU, parameters_type>::compute_G_cols(std::vector<double>&                   exp_V_CPU,
                                                                       LIN_ALG::matrix<double, LIN_ALG::GPU>& N,
                                                                       LIN_ALG::matrix<double, LIN_ALG::GPU>& G,
                                                                       LIN_ALG::matrix<double, LIN_ALG::GPU>& G_cols)
    {
      exp_V.set(exp_V_CPU, LIN_ALG::ASYNCHRONOUS);

      assert(N.get_number_of_rows() == G     .get_number_of_rows());
      assert(N.get_number_of_rows() == G_cols.get_number_of_rows());
      assert(exp_V.size()           == permutation.size());

      int N_i = permutation.size();
      int N_r = N.get_number_of_rows();

      int N_c = N.get_number_of_cols()-G.get_number_of_cols();

      N_MATRIX_TOOLS_GPU_KERNELS::compute_G_cols(N_i, N_r, N_c,
                                                 permutation.get_ptr(), exp_V.get_ptr(),
                                                 N     .get_ptr(), N     .get_leading_dimension(),
                                                 G     .get_ptr(), G     .get_leading_dimension(),
                                                 G_cols.get_ptr(), G_cols.get_leading_dimension(),
                                                 thread_id, stream_id);
    }

  }

}

#endif
