//-*-C++-*-

#ifndef DCA_QMCI_CT_AUX_WALKER_TOOLS_CPU_H
#define DCA_QMCI_CT_AUX_WALKER_TOOLS_CPU_H

namespace DCA
{
  namespace QMCI
  {
    /*!
     *  \ingroup CT-AUX-WALKER
     *
     *  \author Peter Staar
     *  \brief  ...
     */
    template<>
    class CT_AUX_WALKER_TOOLS<LIN_ALG::CPU>
    {
      const static int BLOCK_SIZE = LIN_ALG::matrix<double, LIN_ALG::CPU>::BLOCK_SIZE;

    public:

      CT_AUX_WALKER_TOOLS(int k_ph);
      ~CT_AUX_WALKER_TOOLS();

      inline static void compute_Gamma(LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma,
                                       LIN_ALG::matrix<double, LIN_ALG::CPU>& N,
                                       LIN_ALG::matrix<double, LIN_ALG::CPU>& G_precomputed,
                                       LIN_ALG::vector<int   , LIN_ALG::CPU>& random_vertex_vector,
                                       LIN_ALG::vector<double, LIN_ALG::CPU>& exp_V,
                                       LIN_ALG::vector<double, LIN_ALG::CPU>& exp_delta_V,
                                       int thread_id, int stream_id);

      inline static void set_to_identity(LIN_ALG::matrix<double, LIN_ALG::CPU>& M, int index);

      inline static void remove_row_and_column(LIN_ALG::matrix<double, LIN_ALG::CPU>& M, int index);

      //inline double solve_Gamma(int n, LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma_LU, double exp_delta_V);
      inline double solve_Gamma        (int n, LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma_LU, double exp_delta_V, double& max, double& min);
      inline double solve_Gamma_blocked(int n, LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma_LU, double exp_delta_V, double& max, double& min);

      inline double apply_bennett_on_Gamma(int k, int n, LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma_LU, double phani_gamma, double& max, double& min);

    private:

      inline void solve_Gamma_slow(int n, LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma_LU);
      inline void solve_Gamma_fast(int n, LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma_LU);
      inline void solve_Gamma_BLAS(int n, LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma_LU);

      inline void solve_Gamma_fast(int n, double* A, int LD);

      inline void solve_Gamma_blocked(int n, LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma_LU);

      bool test_max_min(int n, LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma_LU, double max, double min);

    private:

      LIN_ALG::vector<double, LIN_ALG::CPU> r;
      LIN_ALG::vector<double, LIN_ALG::CPU> c;
      LIN_ALG::vector<double, LIN_ALG::CPU> d;
    };

    CT_AUX_WALKER_TOOLS<LIN_ALG::CPU>::CT_AUX_WALKER_TOOLS(int k_ph):
      r(k_ph),
      c(k_ph),
      d(k_ph)
    {}

    CT_AUX_WALKER_TOOLS<LIN_ALG::CPU>::~CT_AUX_WALKER_TOOLS()
    {}

    void CT_AUX_WALKER_TOOLS<LIN_ALG::CPU>::compute_Gamma(LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma,
                                                          LIN_ALG::matrix<double, LIN_ALG::CPU>& N,
                                                          LIN_ALG::matrix<double, LIN_ALG::CPU>& G_precomputed,
                                                          LIN_ALG::vector<int   , LIN_ALG::CPU>& random_vertex_vector,
                                                          LIN_ALG::vector<double, LIN_ALG::CPU>& exp_V,
                                                          LIN_ALG::vector<double, LIN_ALG::CPU>& exp_delta_V,
                                                          int /*thread_id*/, int /*stream_id*/)
    {
      Gamma.resize(random_vertex_vector.size());

      assert(Gamma.get_number_of_rows() == Gamma.get_number_of_cols());

      int vertex_index = N.get_number_of_cols()-G_precomputed.get_number_of_cols();

      for(int i=0; i<Gamma.get_number_of_rows(); i++){
        for(int j=0; j<Gamma.get_number_of_cols(); j++){

          int configuration_e_spin_index_i = random_vertex_vector[i];
          int configuration_e_spin_index_j = random_vertex_vector[j];

          if(configuration_e_spin_index_j<vertex_index)
            {
              double delta = (configuration_e_spin_index_i == configuration_e_spin_index_j)? 1. : 0.;

              double N_ij  = N(configuration_e_spin_index_i,configuration_e_spin_index_j);

              Gamma(i,j) = (N_ij*exp_V[j] - delta)/(exp_V[j]-1.);
            }
          else
            {
              Gamma(i,j) = G_precomputed(configuration_e_spin_index_i, configuration_e_spin_index_j-vertex_index);
            }

          if(i == j){
            double gamma_k = exp_delta_V[j];
            Gamma(i,j) -= (gamma_k)/(gamma_k-1.);
          }
        }
      }
    }

    void CT_AUX_WALKER_TOOLS<LIN_ALG::CPU>::set_to_identity(LIN_ALG::matrix<double, LIN_ALG::CPU>& M, int index)
    {
      int LD_i = M.get_global_size().first;
      int LD_j = M.get_global_size().second;

      {
        double* M_ptr = &M(0,index);
        for(int i=0; i<LD_i; ++i)
          M_ptr[i] = 0.;
      }

      {
        double* M_ptr = &M(index,0);
        for(int j=0; j<LD_j; ++j)
          M_ptr[j*LD_i] = 0.;
      }

      M(index,index) = 1.;
    }

    void CT_AUX_WALKER_TOOLS<LIN_ALG::CPU>::remove_row_and_column(LIN_ALG::matrix<double, LIN_ALG::CPU>& M, int index)
    {
      LIN_ALG::REMOVE<LIN_ALG::CPU>::row_and_column(M, index);
    }

    bool CT_AUX_WALKER_TOOLS<LIN_ALG::CPU>::test_max_min(int n, LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma_LU, double max_ref, double min_ref)
    {
      double Gamma_val = std::fabs(Gamma_LU(0,0));

      double max = Gamma_val;
      double min = Gamma_val;

      for(int i=1; i<n+1; i++){
        Gamma_val = std::fabs(Gamma_LU(i,i));

        max = Gamma_val>max? Gamma_val : max;
        min = Gamma_val<min? Gamma_val : min;
      }

      if(std::fabs(max_ref-max)<1.e-12 and std::fabs(min_ref-min)<1.e-12)
        return true;
      else
        {
          std::cout << __FUNCTION__ << std::endl;
          std::cout.precision(16);
          std::cout << "\n\t n : " << n << "\n";
          std::cout<<std::scientific;
          std::cout << max << "\t" << max_ref << "\t" <<std::fabs(max_ref-max) << std::endl;
          std::cout << min << "\t" << min_ref << "\t" <<std::fabs(min_ref-min) << std::endl;
          std::cout << std::endl;

	  Gamma_LU.print();

          throw std::logic_error(__FUNCTION__);
        }
    }

    double CT_AUX_WALKER_TOOLS<LIN_ALG::CPU>::solve_Gamma(int     n,
                                                          LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma_LU,
                                                          double  exp_delta_V,
                                                          double& max,
                                                          double& min)
    {
      //solve_Gamma_slow(n, Gamma_LU);
      solve_Gamma_fast(n, Gamma_LU);
      //solve_Gamma_BLAS(n, Gamma_LU);

      double Gamma_LU_n_n = Gamma_LU(n,n);
      double Gamma_val    = std::fabs(Gamma_LU_n_n);

      if(n>0)
        {
          double new_max = (Gamma_val>max)? Gamma_val : max;
          double new_min = (Gamma_val<min)? Gamma_val : min;

          if( (new_max/new_min) > 1.e6)
            return 1.e-16;
          else
            {
              max = new_max;
              min = new_min;
            }
        }
      else
        {
          max = Gamma_val;
          min = Gamma_val;
        }

      assert(test_max_min(n, Gamma_LU, max, min));

      double phani_gamma       = exp_delta_V-1.;
      double determinant_ratio = -phani_gamma*Gamma_LU_n_n;

      return determinant_ratio;
    }

    /*!                 /            |   \          /            |      \ /            |            \
     *                  |            |   |          |            |      | |            |      |
     *  \Gamma_{n+1} := |  \Gamma_n  | s |  ---\    |     L_n    |   0  | |     U_n    |   y  |
     *                  |            |   |  ---/    |            |      | |            |      |
     *                  |------------|---|          |------------|------| |------------|------|
     *                  \     w      | d /          \     x      |   1  / \     0      |\beta /
     *
     *   \Gamma_n = L_n*U_n
     *          s = L_n*y
     *          w = x*U_n
     *          d = -x*y+\beta
     */
    void CT_AUX_WALKER_TOOLS<LIN_ALG::CPU>::solve_Gamma_slow(int n, LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma_LU)
    {
      int LD = Gamma_LU.get_global_size().first;

      {
        double* y = Gamma_LU.get_ptr(0,n);
        double* x = Gamma_LU.get_ptr(n,0);

        {
          if(false)
            {//serial
              for(int i=0; i<n; i++)
                for(int j=0; j<i; j++)
                  y[i] -= Gamma_LU(i,j)*y[j];
            }
          else
            {// parallell
              for(int j=0; j<n; j++)
                for(int i=j+1; i<n; i++)
                  y[i] -= Gamma_LU(i,j)*y[j];
            }
        }

        {
          if(true)
            {// serial
              for(int j=0; j<n; j++){
                for(int i=0; i<j; i++)
                  x[j*LD] -= x[i*LD]*Gamma_LU(i,j);
                x[j*LD] /= Gamma_LU(j,j);
              }
            }
          else
            {// parallell
              for(int i=0; i<n; i++){
                x[i*LD] /= Gamma_LU(i,i);
                for(int j=i+1; j<n; j++)
                  x[j*LD] -= x[i*LD]*Gamma_LU(i,j);
              }
            }
        }

        for(int i=0; i<n; i++)
          Gamma_LU(n,n) -= x[i*LD]*y[i];
      }
    }


    /*!                 /            |   \          /            |      \ /            |            \
     *                  |            |   |          |            |      | |            |      |
     *  \Gamma_{n+1} := |  \Gamma_n  | s |  ---\    |     L_n    |   0  | |     U_n    |   y  |
     *                  |            |   |  ---/    |            |      | |            |      |
     *                  |------------|---|          |------------|------| |------------|------|
     *                  \     w      | d /          \     x      |   1  / \     0      |\beta /
     *
     *   \Gamma_n = L_n*U_n
     *          s = L_n*y
     *          w = x*U_n
     *          d = -x*y+\beta
     */
    void CT_AUX_WALKER_TOOLS<LIN_ALG::CPU>::solve_Gamma_fast(int n, LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma_LU)
    {
      solve_Gamma_fast(n, &Gamma_LU(0,0), Gamma_LU.get_leading_dimension());
    }

    void CT_AUX_WALKER_TOOLS<LIN_ALG::CPU>::solve_Gamma_fast(int n, double* A, int LD)
    {
      {
        double* y = &A[0+n*LD];

        {// parallell
          double  y_val = 0;

          double* y_ptr = NULL;
          double* G_ptr = NULL;

          for(int j=0; j<n; j++)
            {
              y_val = y[j];

              y_ptr = &A[0+n*LD];
              G_ptr = &A[0+j*LD];

              for(int i=j+1; i<n; i++)
                y_ptr[i] -= G_ptr[i]*y_val;
            }
        }

        double* x = &r[0];

        {
          for(int j=0; j<n; j++)
            x[j] = A[n+j*LD];
        }

        {// serial
          double x_val = 0;

          double* x_ptr = NULL;
          double* G_ptr = NULL;

          for(int j=0; j<n; j++){
            x_val = x[j];

            x_ptr = x;
            G_ptr = &A[0+j*LD];

            for(int i=0; i<j; i++)
              x_val -= x_ptr[i]*G_ptr[i];

            x[j] = x_val/G_ptr[j];
          }
        }

        {
          double xy = 0;

          for(int i=0; i<n; i++)
            xy += x[i]*y[i];

          A[n+n*LD] -= xy;
        }

        {
          for(int j=0; j<n; j++)
            A[n+j*LD] = x[j];
        }
      }
    }

    /*!                 /            |   \          /            |      \ /            |            \
     *                  |            |   |          |            |      | |            |      |
     *  \Gamma_{n+1} := |  \Gamma_n  | s |  ---\    |     L_n    |   0  | |     U_n    |   y  |
     *                  |            |   |  ---/    |            |      | |            |      |
     *                  |------------|---|          |------------|------| |------------|------|
     *                  \     w      | d /          \     x      |   1  / \     0      |\beta /
     *
     *   \Gamma_n = L_n*U_n
     *          s = L_n*y
     *          w = x*U_n
     *          d = -x*y+\beta
     */
    void CT_AUX_WALKER_TOOLS<LIN_ALG::CPU>::solve_Gamma_BLAS(int n, LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma_LU/*, double exp_delta_V*/)
    {
      int lda = Gamma_LU.get_global_size() .first;

      {
        LIN_ALG::TRSV<LIN_ALG::CPU>::execute('L', 'N', 'U', n, Gamma_LU.get_ptr(0,0), lda, Gamma_LU.get_ptr(0,n), 1);

        LIN_ALG::TRSV<LIN_ALG::CPU>::execute('U', 'T', 'N', n, Gamma_LU.get_ptr(0,0), lda, Gamma_LU.get_ptr(n,0), lda);

        {
          double xy = 0;
          for(int i=0; i<n; i++)
            xy += Gamma_LU(n,i)*Gamma_LU(i,n);

          Gamma_LU(n,n) -= xy;
        }
      }
    }





    double CT_AUX_WALKER_TOOLS<LIN_ALG::CPU>::solve_Gamma_blocked(int     n,
                                                                  LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma_LU,
                                                                  double  exp_delta_V,
                                                                  double& max,
                                                                  double& min)
    {
      //std::cout << "\t(" << min << ", " << max << " ) ";

      solve_Gamma_blocked(n, Gamma_LU);

      double Gamma_LU_n_n = Gamma_LU(n,n);

      double Gamma_val = std::fabs(Gamma_LU_n_n);

      //std::cout << " --> " << Gamma_val << " --> (";

      if(n>0)
        {
          double new_max = (Gamma_val>max)? Gamma_val : max;
          double new_min = (Gamma_val<min)? Gamma_val : min;

          if( (new_max/new_min) > 1.e6)
            return 1.e-16;
          else
            {
              max = new_max;
              min = new_min;
            }
        }
      else
        {
          max = Gamma_val;
          min = Gamma_val;
        }

      //std::cout << min << ", " << max << ")\t";
      assert(test_max_min(n, Gamma_LU, max, min));

      double phani_gamma       = exp_delta_V-1.;
      double determinant_ratio = -phani_gamma*Gamma_LU_n_n;

      return determinant_ratio;
    }


    /*!                 /            |   \          /            |      \ /            |            \
     *                  |            |   |          |            |      | |            |      |
     *  \Gamma_{n+1} := |  \Gamma_n  | s |  ---\    |     L_n    |   0  | |     U_n    |   y  |
     *                  |            |   |  ---/    |            |      | |            |      |
     *                  |------------|---|          |------------|------| |------------|------|
     *                  \     w      | d /          \     x      |   1  / \     0      |\beta /
     *
     *   \Gamma_n = L_n*U_n
     *          s = L_n*y
     *          w = x*U_n
     *          d = -x*y+\beta
     *
     *
     *
     *   /            |      \ /     \       /     \
     *   |            |      | |     |       |     |
     *   |    L_00    |   0  | | y_0 |  ---  | s_0 |
     *   |            |      | |     |  ---  |     |
     *   |------------|------| |-----|       |-----|
     *   \    L_10    | L_11 / \ y_1 /       \ s_1 /
     *
     *   L_00 y_0            = s_0  -->  y_0 = L_00^{-1} * s_0             // L_00^{-1}*s term is precomputed in dtrsm
     *   L_10 y_0 + L_11 y_1 = s_1  -->  y_1 = L_11^{-1} *(s_1 - L_10 y_0)
     *
     *
     *                         /            |      \
     *                         |            |      |
     *                         |    U_00    | U_01 |  ---
     *                         |            |      |  ---
     *                         |------------|------|
     *  [   x_0    |  x_1  ] * \     0      | U_11 /      [   w_0    | w_1 ]
     *
     *   x_0 U_00            = w_0  --> x_0 =  w_0             * U_00^{-1}  //  w_0 * U_00^{-1} term is precomputed in dtrsm
     *   x_0 U_01 + x_1 U_11 = w_1  --> x_1 = (w_1 - x_0 U_01) * U_11^{-1}
     *
     *
     */
    void CT_AUX_WALKER_TOOLS<LIN_ALG::CPU>::solve_Gamma_blocked(int n,
                                                                LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma_LU)
    {
      assert(n>-1 and n<Gamma_LU.get_current_size().first);

      int Nk = BLOCK_SIZE;

      int N  = Gamma_LU.get_current_size().first;
      int LD = Gamma_LU.get_global_size() .first;

      double* A = &Gamma_LU(0,0);

      {
        int l  = n     % Nk;
        int bl = (n-l) / Nk;

        int Ic = bl*Nk;

        assert(n == l+bl*Nk);

        double* A_00 = &A[(bl+0)*Nk+(bl+0)*Nk*LD];

        // update diagonal block
        if(Ic>0 and l>0)
          {
            LIN_ALG::GEMV<LIN_ALG::CPU>::execute('N', l , Ic, -1., &A[Ic+0 *LD], LD, &A[0+n*LD], 1 , 1., &A[Ic+n *LD], 1 );
            LIN_ALG::GEMV<LIN_ALG::CPU>::execute('T', Ic, l , -1., &A[0 +Ic*LD], LD, &A[n+0*LD], LD, 1., &A[n +Ic*LD], LD);
          }

        solve_Gamma_fast(l, A_00, LD);

        {
          double xy = 0;
          for(int i=0; i<Ic; i++)
            xy += A[n+i*LD]*A[i+n*LD];

          A[n+n*LD] -= xy;
        }

        // update non-diagonal block
        if( l>0 and ((l+1) % Nk)==0 and N-(bl+1)*Nk>0)
          {
            {
              double* A_10 = &A[(bl+1)*Nk+(bl+0)*Nk*LD];

              for(int l=0; l<bl; ++l)
                {
                  double* L_il = &A[(bl+1)*Nk+(l +0)*Nk*LD];
                  double* U_li = &A[(l +0)*Nk+(bl+0)*Nk*LD];

                  LIN_ALG::GEMM<LIN_ALG::CPU>::execute('N','N',N-(bl+1)*Nk,Nk,Nk,-1.,L_il,LD, U_li, LD, 1., A_10,LD);
                }

              LIN_ALG::TRSM<LIN_ALG::CPU>::execute('R','U','N','N',N-(bl+1)*Nk,Nk,1,A_00,LD,A_10,LD);
            }

            {
              double* A_01 = &A[(bl+0)*Nk+(bl+1)*Nk*LD];

              for(int l=0; l<bl; ++l)
                {
                  double* L_il = &A[(bl+0)*Nk+(l +0)*Nk*LD];
                  double* U_li = &A[(l +0)*Nk+(bl+1)*Nk*LD];

                  LIN_ALG::GEMM<LIN_ALG::CPU>::execute('N','N',Nk,N-(bl+1)*Nk,Nk,-1.,L_il,LD, U_li,LD, 1., A_01,LD);
                }

              LIN_ALG::TRSM<LIN_ALG::CPU>::execute('L','L','N','U',Nk,N-(bl+1)*Nk,1.,A_00,LD,A_01,LD);
            }
          }
      }
    }



    double CT_AUX_WALKER_TOOLS<LIN_ALG::CPU>::apply_bennett_on_Gamma(int k, int n, LIN_ALG::matrix<double, LIN_ALG::CPU>& Gamma_LU, double exp_delta_V, double& max, double& min)
    {
      int ld = Gamma_LU.get_global_size() .first;

      double* r_ptr = &r[0];
      double* c_ptr = &c[0];
      double* d_ptr = &d[0];

      {// store previous diagonal
        for(int i=0; i<n; ++i)
          d_ptr[i] = Gamma_LU(i,i);
      }

      {// remove the column
        for(int i=0; i<n; ++i)
          r_ptr[i] = 0.;

        for(int i=0; i<n; ++i)
          c_ptr[i] = 0.;

        r_ptr[k] = -1.;
        c_ptr[k] = -1.;

        for(int i=0; i<n; ++i){

          for(int j=0; i>j and j<k; ++j)
            c_ptr[i] += Gamma_LU(i,j)*Gamma_LU(j,k);

          if(i<=k)
            c_ptr[i] += Gamma_LU(i,k);

          if(i>k)
            c_ptr[i] += Gamma_LU(i,k)*Gamma_LU(k,k);
        }

        LIN_ALG::BENNET<LIN_ALG::CPU>::standard_Bennet(n, ld, &Gamma_LU(0,0), c_ptr, r_ptr);
      }

      {// remove the row
        for(int i=0; i<n; ++i)
          r_ptr[i] = 0.;

        for(int i=0; i<n; ++i)
          c_ptr[i] = 0.;

        r_ptr[k] = -1.;
        c_ptr[k] = -1.;

        for(int i=0; i<n; ++i){

          for(int j=0; k>j and j<i; ++j)
            r_ptr[i] += Gamma_LU(k,j)*Gamma_LU(j,i);

          if(k<=i)
            r_ptr[i] += Gamma_LU(k,i);

          if(k>i)
            r_ptr[i] += Gamma_LU(k,i)*Gamma_LU(i,i);
        }

        LIN_ALG::BENNET<LIN_ALG::CPU>::standard_Bennet(n, ld, &Gamma_LU(0,0), c_ptr, r_ptr);
      }

      double ratio=1.;
      for(int i=0; i<n; ++i)
        ratio *= (Gamma_LU(i,i)/d_ptr[i]);

      {
        double Gamma_val = std::fabs(Gamma_LU(0,0));

        max = Gamma_val;
        min = Gamma_val;

        for(int i=1; i<n; i++){
          Gamma_val = std::fabs(Gamma_LU(i,i));

          max = Gamma_val>max? Gamma_val : max;
          min = Gamma_val<min? Gamma_val : min;
        }

        if( (max/min) > 1.e6)
          return 1.e-16;
      }

      double phani_gamma = exp_delta_V-1.;
      double det_ratio   = -ratio/phani_gamma;

      return det_ratio;
    }

  }

}

#endif
