//-*-C++-*-

#ifndef DCA_QMCI_CT_AUX_WALKER_TOOLS_GPU_H
#define DCA_QMCI_CT_AUX_WALKER_TOOLS_GPU_H

namespace DCA
{
  namespace QMCI
  {
    namespace CT_AUX_WALKER_GPU_KERNELS
    {
      void compute_Gamma(double* Gamma, int Gamma_n, int Gamma_ld,
                         double* N, int N_r, int N_c, int N_ld,
                         double* G, int G_r, int G_c, int G_ld,
                         int*    random_vertex_vector,
                         double* exp_V,
                         double* exp_delta_V,
                         int thread_id, int stream_id);
    }

    /*!
     *  \ingroup CT-AUX-WALKER
     *
     *  \author Peter Staar
     *  \brief  ...
     */
    template<>
    class CT_AUX_WALKER_TOOLS<LIN_ALG::GPU>
    {
    public:

      static void compute_Gamma(LIN_ALG::matrix<double, LIN_ALG::GPU>& Gamma,
                                LIN_ALG::matrix<double, LIN_ALG::GPU>& N,
                                LIN_ALG::matrix<double, LIN_ALG::GPU>& G,
                                LIN_ALG::vector<int   , LIN_ALG::GPU>& random_vertex_vector,
                                LIN_ALG::vector<double, LIN_ALG::GPU>& exp_V,
                                LIN_ALG::vector<double, LIN_ALG::GPU>& exp_delta_V,
                                int thread_id, int stream_id)
      {
        Gamma.resize(random_vertex_vector.size());

        assert(Gamma.get_number_of_rows() == Gamma.get_number_of_cols());

        CT_AUX_WALKER_GPU_KERNELS::compute_Gamma(Gamma.get_ptr(), Gamma.get_number_of_rows()                        , Gamma.get_leading_dimension(),
                                                 N.get_ptr(),     N.get_number_of_rows()    , N.get_number_of_cols(), N.get_leading_dimension(),
                                                 G.get_ptr(),     G.get_number_of_rows()    , G.get_number_of_cols(), G.get_leading_dimension(),
                                                 random_vertex_vector.get_ptr(),
                                                 exp_V.get_ptr(),
                                                 exp_delta_V.get_ptr(),
                                                 thread_id, stream_id);
      }

    };

  }

}

#endif
