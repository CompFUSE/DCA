//-*-C++-*-

#ifndef DCA_QMCI_G_MATRIX_ROUTINES_GPU_H
#define DCA_QMCI_G_MATRIX_ROUTINES_GPU_H

namespace DCA
{
  namespace QMCI
  {
    namespace GPU_KERNELS_G_TOOLS
    {

      void read_G_matrix_elements(int     S,
                                  int     vertex_index,
                                  int*    i_index,
                                  int*    j_index,
                                  bool*   is_Bennett,
                                  double* exp_Vj,
                                  double* N,             int LD_N,
                                  double* G_precomputed, int LD_G,
                                  double* result_ptr,
                                  int     incr);

      void compute_row_on_Gamma_matrix(int     row_index,
                                       int     S,
                                       int     vertex_index,
                                       int*    indices,
                                       double* exp_V,
                                       double* N,             int LD_N,
                                       double* G_precomputed, int LD_G,
                                       double* row_ptr,       int incr);

      void compute_col_on_Gamma_matrix(int     row_index,
                                       int     S,
                                       int     vertex_index,
                                       int*    indices,
                                       double* exp_V,
                                       double* N,             int LD_N,
                                       double* G_precomputed, int LD_G,
                                       double* row_ptr,       int incr);

    }

    template<typename parameters_type>
    class G_MATRIX_TOOLS<LIN_ALG::GPU, parameters_type>
    {
      const static int MAX_VERTEX_SINGLETS=4;

      typedef typename parameters_type::concurrency_type concurrency_type;
      typedef typename parameters_type::profiler_type    profiler_t;

    public:

      G_MATRIX_TOOLS(int              id,
                     parameters_type& parameters_ref):
        thread_id(id),

        parameters(parameters_ref),
        concurrency(parameters.get_concurrency()),

        i_index_gpu   (0, MAX_VERTEX_SINGLETS*parameters_ref.get_K_PHANI()),
        j_index_gpu   (0, MAX_VERTEX_SINGLETS*parameters_ref.get_K_PHANI()),
        is_Bennett_gpu(0, MAX_VERTEX_SINGLETS*parameters_ref.get_K_PHANI()),
        exp_Vj_gpu    (0, MAX_VERTEX_SINGLETS*parameters_ref.get_K_PHANI())
      {}

      ~G_MATRIX_TOOLS()
      {}

      void read_G_matrix_elements(LIN_ALG::vector<int   , LIN_ALG::CPU>& i_index,
                                  LIN_ALG::vector<int   , LIN_ALG::CPU>& j_index,
                                  LIN_ALG::vector<bool  , LIN_ALG::CPU>& is_Bennett,
                                  LIN_ALG::vector<double, LIN_ALG::CPU>& exp_Vj,
                                  LIN_ALG::matrix<double, LIN_ALG::GPU>& N,
                                  LIN_ALG::matrix<double, LIN_ALG::GPU>& G_precomputed,
                                  double*                                result_ptr,
                                  int                                    incr)
      {
        assert(i_index.size() == j_index.size());
        assert(i_index.size() == is_Bennett.size());
        assert(i_index.size() == exp_Vj.size());

        int S = i_index.size();

        i_index_gpu   .resize(S);
        j_index_gpu   .resize(S);
        is_Bennett_gpu.resize(S);
        exp_Vj_gpu    .resize(S);

        LIN_ALG::COPY_FROM<LIN_ALG::CPU, LIN_ALG::GPU>::execute(i_index.get_ptr()   , i_index_gpu.get_ptr()   , i_index.size());
        LIN_ALG::COPY_FROM<LIN_ALG::CPU, LIN_ALG::GPU>::execute(j_index.get_ptr()   , j_index_gpu.get_ptr()   , j_index.size());
        LIN_ALG::COPY_FROM<LIN_ALG::CPU, LIN_ALG::GPU>::execute(is_Bennett.get_ptr(), is_Bennett_gpu.get_ptr(), is_Bennett.size());
        LIN_ALG::COPY_FROM<LIN_ALG::CPU, LIN_ALG::GPU>::execute(exp_Vj.get_ptr()    , exp_Vj_gpu.get_ptr()    , exp_Vj.size());

        int vertex_index = N.get_number_of_cols()-G_precomputed.get_number_of_cols();

        GPU_KERNELS_G_TOOLS::read_G_matrix_elements(S,  vertex_index,
                                                    i_index_gpu.get_ptr(),
                                                    j_index_gpu.get_ptr(),
                                                    is_Bennett_gpu.get_ptr(),
                                                    exp_Vj_gpu.get_ptr(),
                                                    N.get_ptr()            , N.get_leading_dimension(),
                                                    G_precomputed.get_ptr(), G_precomputed.get_leading_dimension(),
                                                    result_ptr, incr);
      }

      void compute_row_on_Gamma_matrix(int                                    row_index,
                                       LIN_ALG::vector<int   , LIN_ALG::GPU>& indices,
                                       LIN_ALG::vector<double, LIN_ALG::GPU>& exp_V,
                                       LIN_ALG::matrix<double, LIN_ALG::GPU>& N,
                                       LIN_ALG::matrix<double, LIN_ALG::GPU>& G_precomputed,
                                       double*                                row_ptr,
                                       int                                    incr)
      {
        int vertex_index = N.get_number_of_cols()-G_precomputed.get_number_of_cols();

        GPU_KERNELS_G_TOOLS::compute_row_on_Gamma_matrix(row_index     ,
                                                         indices.size(),
                                                         vertex_index  ,

                                                         indices.get_ptr(),
                                                         exp_V  .get_ptr(),

                                                         N.get_ptr()            , N.get_leading_dimension(),
                                                         G_precomputed.get_ptr(), G_precomputed.get_leading_dimension(),
                                                         row_ptr                , incr);
      }

      void compute_col_on_Gamma_matrix(int                                    col_index,
                                       LIN_ALG::vector<int   , LIN_ALG::GPU>& indices,
                                       LIN_ALG::vector<double, LIN_ALG::GPU>& exp_V,
                                       LIN_ALG::matrix<double, LIN_ALG::GPU>& N,
                                       LIN_ALG::matrix<double, LIN_ALG::GPU>& G_precomputed,
                                       double*                                col_ptr,
                                       int                                    incr)
      {
        int vertex_index = N.get_number_of_cols()-G_precomputed.get_number_of_cols();

        GPU_KERNELS_G_TOOLS::compute_col_on_Gamma_matrix(col_index     ,
                                                         indices.size(),
                                                         vertex_index  ,

                                                         indices.get_ptr(),
                                                         exp_V  .get_ptr(),

                                                         N.get_ptr()            , N.get_leading_dimension(),
                                                         G_precomputed.get_ptr(), G_precomputed.get_leading_dimension(),
                                                         col_ptr                , incr);
      }

    private:

      int thread_id;

      parameters_type&  parameters;
      concurrency_type& concurrency;

      LIN_ALG::vector<int   , LIN_ALG::GPU> i_index_gpu, j_index_gpu;
      LIN_ALG::vector<bool  , LIN_ALG::GPU> is_Bennett_gpu;
      LIN_ALG::vector<double, LIN_ALG::GPU> exp_Vj_gpu;
    };

  }

}


#endif
