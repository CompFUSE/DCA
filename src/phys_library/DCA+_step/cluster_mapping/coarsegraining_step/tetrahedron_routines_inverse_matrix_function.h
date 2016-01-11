//-*-C++-*-

#ifndef DCA_TETRAHEDRON_ROUTINES_INVERSE_MATRIX_FUNCTION_H
#define DCA_TETRAHEDRON_ROUTINES_INVERSE_MATRIX_FUNCTION_H

namespace DCA
{

  class tetrahedron_routines_inverse_matrix_function
  {
    inline static double EPSILON() { return 1.e-2; }

    template<typename scalartype>
    struct matrix_element_struct
    {
      int i;

      std::complex<scalartype> e;
      std::complex<scalartype> f;

      std::complex<scalartype> log_min_e;

      std::complex<scalartype> r;
    };

  private:

    template<typename scalartype>
    inline static scalartype Abs(scalartype& val);

    template<typename scalartype>
    inline static std::complex<scalartype> Abs(std::complex<scalartype>& val);

    template<typename scalartype>
    inline static std::complex<scalartype> Power(std::complex<scalartype> val, int n);

    template<typename scalartype>
    inline static std::complex<scalartype> Log(std::complex<scalartype> val);

    template<typename scalartype>
    static bool are_equal(std::complex<scalartype> const& x,
                          std::complex<scalartype> const& y);

    template<typename scalartype>
    static bool not_equal(std::complex<scalartype> const& x,
                          std::complex<scalartype> const& y);

    template<typename scalartype>
    static bool value_comp(matrix_element_struct<scalartype> const& x,
                           matrix_element_struct<scalartype> const& y);

    template<typename scalartype>
    static bool permutation_comp(matrix_element_struct<scalartype> const& x,
                                 matrix_element_struct<scalartype> const& y);

    template<typename scalartype>
    static bool pair_less(std::pair<std::complex<scalartype>, std::complex<scalartype> > const& x,
                          std::pair<std::complex<scalartype>, std::complex<scalartype> > const& y);

    //     template<typename scalartype>
    //     static bool pair_same(std::pair<complex<scalartype>, complex<scalartype> > const& x,
    //                           std::pair<complex<scalartype>, complex<scalartype> > const& y);

  public:

    // 1D
    template<typename scalartype>
    static void execute(int size, scalartype volume,
                        std::complex<scalartype>* G_0,
                        std::complex<scalartype>* G_1,
                        std::complex<scalartype>* f_result,
                        tetrahedron_integration_data<scalartype>& data_obj);

    // 2D
    template<typename scalartype>
    static void execute(int size, scalartype volume,
                        std::complex<scalartype>* G_0,
                        std::complex<scalartype>* G_1,
                        std::complex<scalartype>* G_2,
                        std::complex<scalartype>* f_result,
                        tetrahedron_integration_data<scalartype>& data_obj);

    // 3D
    template<typename scalartype>
    static void execute(int size, scalartype volume,
                        std::complex<scalartype>* G_0,
                        std::complex<scalartype>* G_1,
                        std::complex<scalartype>* G_2,
                        std::complex<scalartype>* G_3,
                        std::complex<scalartype>* f_result,
                        tetrahedron_integration_data<scalartype>& data_obj);

  private:

    /*
      template<typename scalartype>
      static std::complex<scalartype> integrate_matrix_element_1D(std::complex<scalartype>* f,
      std::complex<scalartype>* e);

      template<typename scalartype>
      static std::complex<scalartype> integrate_matrix_element_2D(std::complex<scalartype>* f,
      std::complex<scalartype>* e);

      template<typename scalartype>
      static std::complex<scalartype> integrate_matrix_element_3D(std::complex<scalartype>* f,
      std::complex<scalartype>* e);
    */

    template<typename scalartype>
    static std::complex<scalartype> integrate_matrix_element_1D(std::vector<matrix_element_struct<scalartype> >& vec);

    template<typename scalartype>
    static std::complex<scalartype> integrate_matrix_element_2D(std::vector<matrix_element_struct<scalartype> >& vec);

    template<typename scalartype>
    static void integrate_eigenvalues_2D(eigenvalue_degeneracy_t                          degeneracy,
                                         std::vector<matrix_element_struct<scalartype> >& vec);

    template<typename scalartype>
    static std::complex<scalartype> integrate_matrix_element_3D(std::vector<matrix_element_struct<scalartype> >& vec);

    template<typename scalartype>
    static std::complex<scalartype> integrate_matrix_element_3D(eigenvalue_degeneracy_t                          degeneracy,
                                                                std::vector<matrix_element_struct<scalartype> >& vec);

    template<typename scalartype>
    static void integrate_eigenvalues_3D(eigenvalue_degeneracy_t                          degeneracy,
                                         std::vector<matrix_element_struct<scalartype> >& vec);

    /*
      template<typename scalartype>
      static eigenvalue_degeneracy_t find_degeneracy_1D(std::complex<scalartype>* f,
      std::complex<scalartype>* e);

      template<typename scalartype>
      static eigenvalue_degeneracy_t find_degeneracy_2D(std::complex<scalartype>* f,
      std::complex<scalartype>* e);

      template<typename scalartype>
      static eigenvalue_degeneracy_t find_degeneracy_3D(std::complex<scalartype>* f,
      std::complex<scalartype>* e);
    */

    template<typename scalartype>
    static eigenvalue_degeneracy_t find_degeneracy_1D(std::vector<matrix_element_struct<scalartype> >& vec);

    template<typename scalartype>
    static eigenvalue_degeneracy_t find_degeneracy_2D(std::vector<matrix_element_struct<scalartype> >& vec);

    template<typename scalartype>
    static eigenvalue_degeneracy_t find_degeneracy_3D(std::vector<matrix_element_struct<scalartype> >& vec);

  };

  template<typename scalartype>
  scalartype tetrahedron_routines_inverse_matrix_function::Abs(scalartype& val)
  {
    return std::abs(val);
  }

  template<typename scalartype>
  std::complex<scalartype> tetrahedron_routines_inverse_matrix_function::Abs(std::complex<scalartype>& val)
  {
    return std::abs(val);
  }

  template<typename scalartype>
  std::complex<scalartype> tetrahedron_routines_inverse_matrix_function::Power(std::complex<scalartype> val, int n)
  {
    switch(n)
      {
      case 1:
        return val;
        break;

      case 2:
        return val*val;
        break;

      case 3:
        return val*val*val;
        break;

      case 4:
        return val*val*val*val;
        break;

      case 5:
        return val*val*val*val*val;
        break;

      default:
        return std::pow(val, n);
      }
  }

  template<typename scalartype>
  std::complex<scalartype> tetrahedron_routines_inverse_matrix_function::Log(std::complex<scalartype> val)
  {
    return std::log(val);
  }

  /*
    template<typename scalartype>
    bool tetrahedron_routines_inverse_matrix_function::are_equal(std::complex<scalartype> const& x,
    std::complex<scalartype> const& y)
    {
    scalartype abs_error = std::abs(x-y);
    scalartype rel_error = std::abs(x-y)/(1.e-16+std::max(std::abs(x), std::abs(y)));

    if(abs_error<EPSILON() or rel_error<EPSILON())
    return true;

    return false;
    }
  */

  template<typename scalartype>
  bool tetrahedron_routines_inverse_matrix_function::are_equal(std::complex<scalartype> const& x,
                                                               std::complex<scalartype> const& y)
  {
    scalartype norm_x = std::norm(x);
    scalartype norm_y = std::norm(y);

    scalartype max = norm_x > norm_y? norm_x : norm_y;

    max += 1.e-16;

    scalartype abs_error = std::norm(x-y);
    scalartype rel_error = abs_error/max;//(1.e-16+std::max(std::norm(x), std::norm(y)));

    if(abs_error<EPSILON()*EPSILON() or
       rel_error<EPSILON()*EPSILON())
      return true;

    return false;
  }

  template<typename scalartype>
  bool tetrahedron_routines_inverse_matrix_function::not_equal(std::complex<scalartype> const& x,
                                                               std::complex<scalartype> const& y)
  {
    return (not are_equal(x,y));
  }

  /*
    template<typename scalartype>
    bool tetrahedron_routines_inverse_matrix_function::value_comp(matrix_element_struct<scalartype> const& x,
    matrix_element_struct<scalartype> const& y)
    {
    return std::abs(x.e)>std::abs(y.e);
    }
  */

  template<typename scalartype>
  bool tetrahedron_routines_inverse_matrix_function::value_comp(matrix_element_struct<scalartype> const& x,
                                                                matrix_element_struct<scalartype> const& y)
  {
    return std::norm(x.e)>std::norm(y.e);
  }

  template<typename scalartype>
  bool tetrahedron_routines_inverse_matrix_function::permutation_comp(matrix_element_struct<scalartype> const& x,
                                                                      matrix_element_struct<scalartype> const& y)
  {
    return x.i<y.i;
  }

  //   template<typename scalartype>
  //   bool tetrahedron_routines_inverse_matrix_function::pair_same(std::pair<complex<scalartype>, complex<scalartype> > const& x,
  //                                                                std::pair<complex<scalartype>, complex<scalartype> > const& y)
  //   {
  //     throw std::logic_error(__FUNCTION__);
  //     return false;
  //   }

  template<typename scalartype>
  bool tetrahedron_routines_inverse_matrix_function::pair_less(std::pair<std::complex<scalartype>, std::complex<scalartype> > const& x,
                                                               std::pair<std::complex<scalartype>, std::complex<scalartype> > const& y)
  {
    return std::abs(x.first) > std::abs(y.first);
  }

  //   template<>
  //   bool tetrahedron_routines_inverse_matrix_function::pair_same(std::pair<complex<float>, complex<float> > const& x,
  //                                                                std::pair<complex<float>, complex<float> > const& y)
  //   {
  //     float abs_x = std::abs(x.first);
  //     float abs_y = std::abs(y.first);

  //     if(abs_x < 1. && abs_y < 1.)
  //       {
  //         return std::abs(x.first-y.first)<1.e-2;
  //       }
  //     else
  //       {
  //         float MAX = abs_x>abs_y? abs_x:abs_y;
  //         return std::abs(x.first-y.first)<((1.e-2)*MAX);
  //       }
  //   }

  //   template<>
  //   bool tetrahedron_routines_inverse_matrix_function::pair_same(std::pair<complex<double>, complex<double> > const& x,
  //                                                                std::pair<complex<double>, complex<double> > const& y)
  //   {
  //     //     double abs_x = abs(x.first);
  //     //     double abs_y = abs(y.first);

  //     //     if(abs_x < 1. && abs_y < 1.)
  //     //       {
  //     //         return abs(x.first-y.first)<1.e-6;
  //     //       }
  //     //     else
  //     //       {
  //     //         double MAX = abs_x>abs_y? abs_x:abs_y;
  //     //         return abs(x.first-y.first)<((1.e-6)*MAX);
  //     //       }

  //     double abs_x = norm(x.first);
  //     double abs_y = norm(y.first);

  //     if(abs_x < 1. && abs_y < 1.)
  //       {
  //    return norm(x.first-y.first)<1.e-3;
  //       }
  //     else
  //       {
  //         double MAX = abs_x>abs_y? abs_x:abs_y;
  //         return norm(x.first-y.first)<((1.e-3)*MAX);
  //       }
  //   }

  /************************************
   ***
   ***   1D tetrahedron-integration
   ***
   ************************************/

  template<typename scalartype>
  void tetrahedron_routines_inverse_matrix_function::execute(int N, scalartype volume,
                                                             std::complex<scalartype>* G_0,
                                                             std::complex<scalartype>* G_1,
                                                             std::complex<scalartype>* f_result,
                                                             tetrahedron_integration_data<scalartype>& data_obj)
  {
    //tetrahedron_integration_data<scalartype> data_obj(N);

    {// diagonolize the G-matrices
      {
        memcpy(data_obj.G_inv_0, G_0, sizeof(std::complex<scalartype>)*N*N);
        memcpy(data_obj.G_inv_1, G_1, sizeof(std::complex<scalartype>)*N*N);

        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_0, data_obj.GEINV_IPIV, data_obj.GEINV_WORK, data_obj.GEINV_LWORK, data_obj.INFO); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_0, data_obj.G_inv_0));
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_1, data_obj.GEINV_IPIV, data_obj.GEINV_WORK, data_obj.GEINV_LWORK, data_obj.INFO); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_1, data_obj.G_inv_1));
      }

      {
        //         int INFO  = -1;
        //         int LWORK = 16*std::max(1,2*N-1);

        //         scalartype                RWORK[std::max(1, 3*N-2)];
        //         std::complex<scalartype>*  WORK = new std::complex<scalartype>[LWORK];

        LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'V', N, data_obj.G_inv_0, N, data_obj.W_0, data_obj.VR_inv_0, N, data_obj.VR_0, N, data_obj.GEEV_WORK, data_obj.GEEV_LWORK, data_obj.GEEV_RWORK, data_obj.INFO);
        LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'V', N, data_obj.G_inv_1, N, data_obj.W_1, data_obj.VR_inv_1, N, data_obj.VR_1, N, data_obj.GEEV_WORK, data_obj.GEEV_LWORK, data_obj.GEEV_RWORK, data_obj.INFO);

        //         delete [] WORK;
      }

      {
        memcpy(data_obj.VR_inv_0, data_obj.VR_0, sizeof(std::complex<scalartype>)*N*N);
        memcpy(data_obj.VR_inv_1, data_obj.VR_1, sizeof(std::complex<scalartype>)*N*N);

        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_0, data_obj.GEINV_IPIV, data_obj.GEINV_WORK, data_obj.GEINV_LWORK, data_obj.INFO); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_0, data_obj.VR_inv_0));
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_1, data_obj.GEINV_IPIV, data_obj.GEINV_WORK, data_obj.GEINV_LWORK, data_obj.INFO); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_1, data_obj.VR_inv_1));
      }
    }

    {// integrate G-matrices

      for(int l=0; l<N*N; l++)
        f_result[l] = 0;

      std::vector<matrix_element_struct<scalartype> > vec(2);

      for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){

          for(int l=0; l<N; l++){

            vec[0].i = 0;
            vec[1].i = 1;

            vec[0].e = data_obj.W_0[l];
            vec[1].e = data_obj.W_1[l];

            //             vec[0].f = data_obj.VR_inv_0[i+l*N]*data_obj.VR_0[l+j*N];
            //             vec[1].f = data_obj.VR_inv_1[i+l*N]*data_obj.VR_1[l+j*N];

            vec[0].f = data_obj.VR_0[i+l*N]*data_obj.VR_inv_0[l+j*N];
            vec[1].f = data_obj.VR_1[i+l*N]*data_obj.VR_inv_1[l+j*N];

            f_result[i+j*N] += integrate_matrix_element_1D(vec);
          }
        }
      }

      for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
          f_result[i + j*N] *= volume;
    }
  }

  template<typename scalartype>
  std::complex<scalartype> tetrahedron_routines_inverse_matrix_function::integrate_matrix_element_1D(std::vector<matrix_element_struct<scalartype> >& vec)
  {
    assert(vec.size()==2);

    eigenvalue_degeneracy_t degeneracy = find_degeneracy_1D(vec);

    std::complex<scalartype> r[2];
    std::complex<scalartype> f[3];

    f[0] = vec[0].f;
    f[1] = vec[1].f;

    std::complex<scalartype> e0 = vec[0].e;
    std::complex<scalartype> e1 = vec[1].e;

    switch(degeneracy)
      {
      case NO_DEGENERACY:
        {
          r[0] = (e0 - e1 - e1*Log(-e0) + e1*Log(-e1))/Power(e0 - e1,2);
          r[1] = (-e0 + e1 + e0*Log(-e0) - e0*Log(-e1))/Power(e0 - e1,2);
        }
        break;

      case TWOFOLD_DEGENERACY:
        {
          r[0] = 1/(2.*e0);
          r[1] = 1/(2.*e0);
        }
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }

    std::complex<scalartype> result=0;

    {
      result += f[0]*r[0];
      result += f[1]*r[1];

      assert(result==result); // make sure there is no NAN !
    }

    return result;
  }

  template<typename scalartype>
  eigenvalue_degeneracy_t tetrahedron_routines_inverse_matrix_function::find_degeneracy_1D(std::vector<matrix_element_struct<scalartype> >& vec)
  {
    assert(vec.size()==2);

    if(are_equal(vec[0].e, vec[1].e))
      return TWOFOLD_DEGENERACY;

    return NO_DEGENERACY;
  }

  /************************************
   ***
   ***   2D tetrahedron-integration
   ***
   ************************************/

  template<typename scalartype>
  void tetrahedron_routines_inverse_matrix_function::execute(int N, scalartype volume,
                                                             std::complex<scalartype>* G_0,
                                                             std::complex<scalartype>* G_1,
                                                             std::complex<scalartype>* G_2,
                                                             std::complex<scalartype>* f_result,
                                                             tetrahedron_integration_data<scalartype>& data_obj)
  {
    //tetrahedron_integration_data<scalartype> data_obj(N);

    {// diagonolize the G-matrices
      {
        memcpy(data_obj.G_inv_0, G_0, sizeof(std::complex<scalartype>)*N*N);
        memcpy(data_obj.G_inv_1, G_1, sizeof(std::complex<scalartype>)*N*N);
        memcpy(data_obj.G_inv_2, G_2, sizeof(std::complex<scalartype>)*N*N);

        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_0, data_obj.GEINV_IPIV, data_obj.GEINV_WORK, data_obj.GEINV_LWORK, data_obj.INFO); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_0, data_obj.G_inv_0));
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_1, data_obj.GEINV_IPIV, data_obj.GEINV_WORK, data_obj.GEINV_LWORK, data_obj.INFO); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_1, data_obj.G_inv_1));
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_2, data_obj.GEINV_IPIV, data_obj.GEINV_WORK, data_obj.GEINV_LWORK, data_obj.INFO); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_2, data_obj.G_inv_2));
      }

      {
        //         int INFO  = -1;
        //         int LWORK = 16*std::max(1,2*N-1);

        //         scalartype                RWORK[std::max(1, 3*N-2)];
        //         std::complex<scalartype>*  WORK = new std::complex<scalartype>[LWORK];

        LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'V', N, data_obj.G_inv_0, N, data_obj.W_0, data_obj.VR_inv_0, N, data_obj.VR_0, N, data_obj.GEEV_WORK, data_obj.GEEV_LWORK, data_obj.GEEV_RWORK, data_obj.INFO);
        LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'V', N, data_obj.G_inv_1, N, data_obj.W_1, data_obj.VR_inv_1, N, data_obj.VR_1, N, data_obj.GEEV_WORK, data_obj.GEEV_LWORK, data_obj.GEEV_RWORK, data_obj.INFO);
        LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'V', N, data_obj.G_inv_2, N, data_obj.W_2, data_obj.VR_inv_2, N, data_obj.VR_2, N, data_obj.GEEV_WORK, data_obj.GEEV_LWORK, data_obj.GEEV_RWORK, data_obj.INFO);

        //         delete [] WORK;
      }

      {
        memcpy(data_obj.VR_inv_0, data_obj.VR_0, sizeof(std::complex<scalartype>)*N*N);
        memcpy(data_obj.VR_inv_1, data_obj.VR_1, sizeof(std::complex<scalartype>)*N*N);
        memcpy(data_obj.VR_inv_2, data_obj.VR_2, sizeof(std::complex<scalartype>)*N*N);

        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_0, data_obj.GEINV_IPIV, data_obj.GEINV_WORK, data_obj.GEINV_LWORK, data_obj.INFO); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_0, data_obj.VR_inv_0));
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_1, data_obj.GEINV_IPIV, data_obj.GEINV_WORK, data_obj.GEINV_LWORK, data_obj.INFO); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_1, data_obj.VR_inv_1));
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_2, data_obj.GEINV_IPIV, data_obj.GEINV_WORK, data_obj.GEINV_LWORK, data_obj.INFO); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_2, data_obj.VR_inv_2));
      }
    }

    {// integrate G-matrices

      for(int l=0; l<N*N; l++)
        f_result[l] = 0;

      std::vector<matrix_element_struct<scalartype> > vec(3);

      if(false)
        {
          for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){

              for(int l=0; l<N; l++){

                {
                  vec[0].i = 0;
                  vec[1].i = 1;
                  vec[2].i = 2;

                  vec[0].e = data_obj.W_0[l];
                  vec[1].e = data_obj.W_1[l];
                  vec[2].e = data_obj.W_2[l];

                  vec[0].f = data_obj.VR_0[i+l*N]*data_obj.VR_inv_0[l+j*N];
                  vec[1].f = data_obj.VR_1[i+l*N]*data_obj.VR_inv_1[l+j*N];
                  vec[2].f = data_obj.VR_2[i+l*N]*data_obj.VR_inv_2[l+j*N];

                  f_result[i+j*N] += integrate_matrix_element_2D(vec);
                }
              }
            }
          }
        }
      else
        {
          for(int l=0; l<N; l++){

            vec[0].i = 0;
            vec[1].i = 1;
            vec[2].i = 2;

            vec[0].e = data_obj.W_0[l];
            vec[1].e = data_obj.W_1[l];
            vec[2].e = data_obj.W_2[l];

            vec[0].log_min_e = Log(-vec[0].e);
            vec[1].log_min_e = Log(-vec[1].e);
            vec[2].log_min_e = Log(-vec[2].e);

            eigenvalue_degeneracy_t degeneracy = find_degeneracy_2D(vec);

            integrate_eigenvalues_2D(degeneracy, vec);

            for(int i=0; i<N; i++){
              for(int j=0; j<N; j++){

                vec[0].f = data_obj.VR_0[i+l*N]*data_obj.VR_inv_0[l+j*N];
                vec[1].f = data_obj.VR_1[i+l*N]*data_obj.VR_inv_1[l+j*N];
                vec[2].f = data_obj.VR_2[i+l*N]*data_obj.VR_inv_2[l+j*N];

                f_result[i+j*N] += vec[0].r*vec[0].f;
                f_result[i+j*N] += vec[1].r*vec[1].f;
                f_result[i+j*N] += vec[2].r*vec[2].f;
              }
            }
          }

        }

      for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
          f_result[i + j*N] *= (2.*volume);
    }
  }

  /*
    template<typename scalartype>
    std::complex<scalartype> tetrahedron_routines_inverse_matrix_function::integrate_matrix_element_2D(std::complex<scalartype>* f,
    std::complex<scalartype>* e)
    {
    const static scalartype VALUE[8] = {0,1,2,3,4,5,6,7};

    std::complex<scalartype> r[2+1];

    eigenvalue_degeneracy_t degeneracy = find_degeneracy_2D(f, e);

    switch(degeneracy)
    {
    case NO_DEGENERACY:
    {
    r[0] = (-(e[0]*(e[1] - e[2])*(-VALUE[2]*e[1]*e[2] + e[0]*(e[1] + e[2]))*std::log(-e[0])) + e[0]*(e[1] - e[2])*(-VALUE[2]*e[1]*e[2] + e[0]*(e[1] + e[2]))*std::log(-e[1]) + (e[0] - e[1])*(e[0]*(e[0] - e[2])*(e[1] - e[2]) + (e[0] - e[1])*std::pow(e[2],2)*std::log(e[1]) + (-e[0] + e[1])*std::pow(e[2],2)*std::log(e[2])))/(VALUE[2]*std::pow(e[0] - e[1],2)*std::pow(e[0] - e[2],2)*(e[1] - e[2]));
    r[1] = -(-(std::pow(e[0],2)*std::pow(e[1] - e[2],2)*std::log(-e[0])) + std::pow(e[0],2)*std::pow(e[1] - e[2],2)*std::log(-e[1]) + (e[0] - e[1])*(e[1]*(e[0] - e[2])*(e[1] - e[2]) + (-e[0] + e[1])*std::pow(e[2],2)*std::log(e[1]) + (e[0] - e[1])*std::pow(e[2],2)*std::log(e[2])))/(VALUE[2]*std::pow(e[0] - e[1],2)*(e[0] - e[2])*std::pow(e[1] - e[2],2));
    r[2] = -(-(std::pow(e[0],2)*std::pow(e[1] - e[2],2)*std::log(-e[0])) + std::pow(e[0],2)*std::pow(e[1] - e[2],2)*std::log(-e[1]) + (e[0] - e[1])*e[2]*((e[0] - e[2])*(-e[1] + e[2]) + (VALUE[2]*e[0]*e[1] - e[0]*e[2] - e[1]*e[2])*std::log(e[1]) + (-VALUE[2]*e[0]*e[1] + e[0]*e[2] + e[1]*e[2])*std::log(e[2])))/(VALUE[2]*(e[0] - e[1])*std::pow(e[0] - e[2],2)*std::pow(e[1] - e[2],2));
    }
    break;

    case TWOFOLD_DEGENERACY:
    {
    r[0] = (std::pow(e[0],2) - std::pow(e[1],2) - VALUE[2]*e[0]*e[1]*std::log(-e[0]) + VALUE[2]*e[0]*e[1]*std::log(-e[1]))/(VALUE[2]*std::pow(e[0] - e[1],3));
    r[1] = -(VALUE[3]*std::pow(e[0],2) - VALUE[4]*e[0]*e[1] + std::pow(e[1],2) - VALUE[2]*std::pow(e[0],2)*std::log(-e[0]) + VALUE[2]*std::pow(e[0],2)*std::log(-e[1]))/(VALUE[4]*std::pow(e[0] - e[1],3));
    r[2] = -(VALUE[3]*std::pow(e[0],2) - VALUE[4]*e[0]*e[1] + std::pow(e[1],2) - VALUE[2]*std::pow(e[0],2)*std::log(-e[0]) + VALUE[2]*std::pow(e[0],2)*std::log(-e[1]))/(VALUE[4]*std::pow(e[0] - e[1],3));
    }
    break;

    case THREEFOLD_DEGENERACY:
    r[0] = VALUE[1]/VALUE[6]*VALUE[1]/e[0];
    r[1] = VALUE[1]/VALUE[6]*VALUE[1]/e[0];
    r[2] = VALUE[1]/VALUE[6]*VALUE[1]/e[0];
    break;

    default:
    throw std::logic_error(__FUNCTION__);
    }

    assert(f[0]*r[0]+f[1]*r[1]+f[2]*r[2] == f[0]*r[0]+f[1]*r[1]+f[2]*r[2]);

    return f[0]*r[0]+f[1]*r[1]+f[2]*r[2];
    }

    template<typename scalartype>
    eigenvalue_degeneracy_t tetrahedron_routines_inverse_matrix_function::find_degeneracy_2D(std::complex<scalartype>* f,
    std::complex<scalartype>* e)
    {
    std::vector<std::pair<std::complex<scalartype>, std::complex<scalartype> > > vec(3);

    vec[0].first = e[0]; vec[0].second = f[0];
    vec[1].first = e[1]; vec[1].second = f[1];
    vec[2].first = e[2]; vec[2].second = f[2];

    std::stable_sort            (vec.begin(), vec.end(), tetrahedron_routines_inverse_matrix_function::pair_less<scalartype>);
    int degeneracy = std::unique(vec.begin(), vec.end(), tetrahedron_routines_inverse_matrix_function::pair_same<scalartype>)-vec.begin();

    e[0] = vec[0].first; f[0] = vec[0].second;
    e[1] = vec[1].first; f[1] = vec[1].second;
    e[2] = vec[2].first; f[2] = vec[2].second;

    if(degeneracy == 3)
    return NO_DEGENERACY;

    if(degeneracy == 2)
    return TWOFOLD_DEGENERACY;

    if(degeneracy == 1)
    return THREEFOLD_DEGENERACY;

    throw std::logic_error(__FUNCTION__);
    return NO_DEGENERACY;
    }
  */

  template<typename scalartype>
  std::complex<scalartype> tetrahedron_routines_inverse_matrix_function::integrate_matrix_element_2D(std::vector<matrix_element_struct<scalartype> >& vec)
  {
    assert(vec.size()==3);

    eigenvalue_degeneracy_t degeneracy = find_degeneracy_2D(vec);

    std::complex<scalartype> r[3];
    std::complex<scalartype> f[3];

    f[0] = vec[0].f;
    f[1] = vec[1].f;
    f[2] = vec[2].f;

    std::complex<scalartype> e0 = vec[0].e;
    std::complex<scalartype> e1 = vec[1].e;
    std::complex<scalartype> e2 = vec[2].e;

    switch(degeneracy)
      {
      case NO_DEGENERACY:
        {
          assert(not_equal(vec[0].e, vec[1].e));
          assert(not_equal(vec[1].e, vec[2].e));
          assert(not_equal(vec[2].e, vec[0].e));

          r[0] = (-(e0*(e1 - e2)*(-2.*e1*e2 + e0*(e1 + e2))*Log(-e0)) + Power(e1,2)*Power(e0 - e2,2)*Log(-e1) + (e0 - e1)*(e0*(e0 - e2)*(e1 - e2) + (-e0 + e1)*Power(e2,2)*Log(-e2)))/(2.*Power(e0 - e1,2)*Power(e0 - e2,2)*(e1 - e2));
          r[1] = (Power(e0,2)*Power(e1 - e2,2)*Log(-e0) + e1*(e0 - e2)*(-((e0 - e1)*(e1 - e2)) - (e0*e1 - 2.*e0*e2 + e1*e2)*Log(-e1)) - Power(e0 - e1,2)*Power(e2,2)*Log(-e2))/(2.*Power(e0 - e1,2)*(e0 - e2)*Power(e1 - e2,2));
          r[2] = (Power(e0,2)*Power(e1 - e2,2)*Log(-e0) - Power(e1,2)*Power(e0 - e2,2)*Log(-e1) + (-e0 + e1)*e2*((e0 - e2)*(-e1 + e2) + (-2.*e0*e1 + (e0 + e1)*e2)*Log(-e2)))/(2.*(e0 - e1)*Power(e0 - e2,2)*Power(e1 - e2,2));
        }
        break;

      case TWOFOLD_DEGENERACY:
        {
          assert(not_equal(vec[0].e, vec[1].e));
          assert(are_equal(vec[1].e, vec[2].e));

          r[0] = (Power(e0,2) - Power(e1,2) - 2.*e0*e1*Log(-e0) + 2.*e0*e1*Log(-e1))/(2.*Power(e0 - e1,3));
          r[1] = -(3.*Power(e0,2) - 4.*e0*e1 + Power(e1,2) + 2.*Power(e0,2)*(-Log(-e0) + Log(-e1)))/(4.*Power(e0 - e1,3));
          r[2] = -(3.*Power(e0,2) - 4.*e0*e1 + Power(e1,2) + 2.*Power(e0,2)*(-Log(-e0) + Log(-e1)))/(4.*Power(e0 - e1,3));
        }
        break;

      case THREEFOLD_DEGENERACY:
        {
          assert(are_equal(vec[0].e, vec[1].e));
          assert(are_equal(vec[0].e, vec[2].e));

          r[0] = 1./(6.*e0);
          r[1] = 1./(6.*e0);
          r[2] = 1./(6.*e0);
        }
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }

    std::complex<scalartype> result=0;

    {
      result += f[0]*r[0];
      result += f[1]*r[1];
      result += f[2]*r[2];

      assert(result==result); // make sure there is no NAN !
    }

    return result;
  }

  template<typename scalartype>
  void tetrahedron_routines_inverse_matrix_function::integrate_eigenvalues_2D(eigenvalue_degeneracy_t                          degeneracy,
                                                                              std::vector<matrix_element_struct<scalartype> >& vec)
  {
    assert(vec.size()==3);

    int i0 = vec[0].i;
    int i1 = vec[1].i;
    int i2 = vec[2].i;

    std::complex<scalartype> e0 = vec[0].e;
    std::complex<scalartype> e1 = vec[1].e;
    std::complex<scalartype> e2 = vec[2].e;

    std::complex<scalartype> log_min_e0 = vec[0].log_min_e;
    std::complex<scalartype> log_min_e1 = vec[1].log_min_e;
    std::complex<scalartype> log_min_e2 = vec[2].log_min_e;

    switch(degeneracy)
      {
      case NO_DEGENERACY:
        {
          assert(not_equal(vec[0].e, vec[1].e));
          assert(not_equal(vec[1].e, vec[2].e));
          assert(not_equal(vec[2].e, vec[0].e));

          vec[i0].r = (-(e0*(e1 - e2)*(-2.*e1*e2 + e0*(e1 + e2))*log_min_e0) + Power(e1,2)*Power(e0 - e2,2)*log_min_e1 + (e0 - e1)*(e0*(e0 - e2)*(e1 - e2) + (-e0 + e1)*Power(e2,2)*log_min_e2))/(2.*Power(e0 - e1,2)*Power(e0 - e2,2)*(e1 - e2));
          vec[i1].r = (Power(e0,2)*Power(e1 - e2,2)*log_min_e0 + e1*(e0 - e2)*(-((e0 - e1)*(e1 - e2)) - (e0*e1 - 2.*e0*e2 + e1*e2)*log_min_e1) - Power(e0 - e1,2)*Power(e2,2)*log_min_e2)/(2.*Power(e0 - e1,2)*(e0 - e2)*Power(e1 - e2,2));
          vec[i2].r = (Power(e0,2)*Power(e1 - e2,2)*log_min_e0 - Power(e1,2)*Power(e0 - e2,2)*log_min_e1 + (-e0 + e1)*e2*((e0 - e2)*(-e1 + e2) + (-2.*e0*e1 + (e0 + e1)*e2)*log_min_e2))/(2.*(e0 - e1)*Power(e0 - e2,2)*Power(e1 - e2,2));
        }
        break;

      case TWOFOLD_DEGENERACY:
        {
          assert(not_equal(vec[0].e, vec[1].e));
          assert(are_equal(vec[1].e, vec[2].e));

          vec[i0].r = (Power(e0,2) - Power(e1,2) - 2.*e0*e1*log_min_e0 + 2.*e0*e1*log_min_e1)/(2.*Power(e0 - e1,3));
          vec[i1].r = -(3.*Power(e0,2) - 4.*e0*e1 + Power(e1,2) + 2.*Power(e0,2)*(-log_min_e0 + log_min_e1))/(4.*Power(e0 - e1,3));
          vec[i2].r = -(3.*Power(e0,2) - 4.*e0*e1 + Power(e1,2) + 2.*Power(e0,2)*(-log_min_e0 + log_min_e1))/(4.*Power(e0 - e1,3));
        }
        break;

      case THREEFOLD_DEGENERACY:
        {
          assert(are_equal(vec[0].e, vec[1].e));
          assert(are_equal(vec[0].e, vec[2].e));

          vec[i0].r = 1./(6.*e0);
          vec[i1].r = 1./(6.*e0);
          vec[i2].r = 1./(6.*e0);
        }
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }
  }

  template<typename scalartype>
  eigenvalue_degeneracy_t tetrahedron_routines_inverse_matrix_function::find_degeneracy_2D(std::vector<matrix_element_struct<scalartype> >& vec)
  {
    assert(vec.size()==3);

    if(not_equal(vec[0].e, vec[1].e) and
       not_equal(vec[1].e, vec[2].e) and
       not_equal(vec[2].e, vec[0].e) )
      return NO_DEGENERACY;

    if(are_equal(vec[0].e, vec[1].e) and
       are_equal(vec[0].e, vec[2].e))
      return THREEFOLD_DEGENERACY;

    do
      {
        if(not_equal(vec[0].e, vec[1].e) and
           are_equal(vec[1].e, vec[2].e))
          return TWOFOLD_DEGENERACY;
      }
    while
      ( std::next_permutation(vec.begin(), vec.end(), tetrahedron_routines_inverse_matrix_function::permutation_comp<scalartype>) );

    throw std::logic_error(__FUNCTION__);
    return NO_DEGENERACY;
  }

  /************************************
   ***
   ***   3D tetrahedron-integration
   ***
   ************************************/

  template<typename scalartype>
  void tetrahedron_routines_inverse_matrix_function::execute(int N, scalartype volume,
                                                             std::complex<scalartype>* G_0,
                                                             std::complex<scalartype>* G_1,
                                                             std::complex<scalartype>* G_2,
                                                             std::complex<scalartype>* G_3,
                                                             std::complex<scalartype>* f_result,
                                                             tetrahedron_integration_data<scalartype>& data_obj)
  {
    //     tetrahedron_integration_data<scalartype> data_obj(N);
    //     clock_t t;

    {
      //       t = clock();

      {// obtain G^{-1} from G
        memcpy(data_obj.G_inv_0, G_0, sizeof(std::complex<scalartype>)*N*N);
        memcpy(data_obj.G_inv_1, G_1, sizeof(std::complex<scalartype>)*N*N);
        memcpy(data_obj.G_inv_2, G_2, sizeof(std::complex<scalartype>)*N*N);
        memcpy(data_obj.G_inv_3, G_3, sizeof(std::complex<scalartype>)*N*N);

        //         LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_0); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_0, data_obj.G_inv_0));
        //         LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_1); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_1, data_obj.G_inv_1));
        //         LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_2); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_2, data_obj.G_inv_2));
        //         LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_3); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_3, data_obj.G_inv_3));

        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_0, data_obj.GEINV_IPIV, data_obj.GEINV_WORK, data_obj.GEINV_LWORK, data_obj.INFO); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_0, data_obj.G_inv_0));
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_1, data_obj.GEINV_IPIV, data_obj.GEINV_WORK, data_obj.GEINV_LWORK, data_obj.INFO); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_1, data_obj.G_inv_1));
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_2, data_obj.GEINV_IPIV, data_obj.GEINV_WORK, data_obj.GEINV_LWORK, data_obj.INFO); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_2, data_obj.G_inv_2));
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_3, data_obj.GEINV_IPIV, data_obj.GEINV_WORK, data_obj.GEINV_LWORK, data_obj.INFO); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_3, data_obj.G_inv_3));
      }

      {// diagonolize the G-matrices
        //         int INFO  = -1;
        //         int LWORK = 16*std::max(1,2*N-1);

        //         scalartype                RWORK[std::max(1, 3*N-2)];
        //         std::complex<scalartype>*  WORK = new std::complex<scalartype>[LWORK];

        LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'V', N, data_obj.G_inv_0, N, data_obj.W_0, data_obj.VR_inv_0, N, data_obj.VR_0, N, data_obj.GEEV_WORK, data_obj.GEEV_LWORK, data_obj.GEEV_RWORK, data_obj.INFO);
        LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'V', N, data_obj.G_inv_1, N, data_obj.W_1, data_obj.VR_inv_1, N, data_obj.VR_1, N, data_obj.GEEV_WORK, data_obj.GEEV_LWORK, data_obj.GEEV_RWORK, data_obj.INFO);
        LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'V', N, data_obj.G_inv_2, N, data_obj.W_2, data_obj.VR_inv_2, N, data_obj.VR_2, N, data_obj.GEEV_WORK, data_obj.GEEV_LWORK, data_obj.GEEV_RWORK, data_obj.INFO);
        LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'V', N, data_obj.G_inv_3, N, data_obj.W_3, data_obj.VR_inv_3, N, data_obj.VR_3, N, data_obj.GEEV_WORK, data_obj.GEEV_LWORK, data_obj.GEEV_RWORK, data_obj.INFO);

        //         delete [] WORK;
      }

      {// obtain V^{-1}
        memcpy(data_obj.VR_inv_0, data_obj.VR_0, sizeof(std::complex<scalartype>)*N*N);
        memcpy(data_obj.VR_inv_1, data_obj.VR_1, sizeof(std::complex<scalartype>)*N*N);
        memcpy(data_obj.VR_inv_2, data_obj.VR_2, sizeof(std::complex<scalartype>)*N*N);
        memcpy(data_obj.VR_inv_3, data_obj.VR_3, sizeof(std::complex<scalartype>)*N*N);

        //         LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_0); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_0, data_obj.VR_inv_0));
        //         LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_1); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_1, data_obj.VR_inv_1));
        //         LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_2); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_2, data_obj.VR_inv_2));
        //         LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_3); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_3, data_obj.VR_inv_3));

        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_0, data_obj.GEINV_IPIV, data_obj.GEINV_WORK, data_obj.GEINV_LWORK, data_obj.INFO); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_0, data_obj.VR_inv_0));
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_1, data_obj.GEINV_IPIV, data_obj.GEINV_WORK, data_obj.GEINV_LWORK, data_obj.INFO); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_1, data_obj.VR_inv_1));
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_2, data_obj.GEINV_IPIV, data_obj.GEINV_WORK, data_obj.GEINV_LWORK, data_obj.INFO); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_2, data_obj.VR_inv_2));
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_3, data_obj.GEINV_IPIV, data_obj.GEINV_WORK, data_obj.GEINV_LWORK, data_obj.INFO); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_3, data_obj.VR_inv_3));
      }

      //       cout << " inversion and diag : " << double(clock()-t)/CLOCKS_PER_SEC << "\n";
    }

    {// integrate G-matrices
      //       t = clock();

      for(int j=0; j<N; j++)
        for(int i=0; i<N; i++)
          f_result[i + j*N] = 0.;

      //       std::complex<scalartype> matrix_elements[4];
      //       std::complex<scalartype> eigenvalues    [4];

      std::vector<matrix_element_struct<scalartype> > vec(4);

      /*
        if(true)
        {
        for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){

        for(int l=0; l<N; l++){

        vec[0].i = 0;
        vec[1].i = 1;
        vec[2].i = 2;
        vec[3].i = 3;

        vec[0].e = data_obj.W_0[l];
        vec[1].e = data_obj.W_1[l];
        vec[2].e = data_obj.W_2[l];
        vec[3].e = data_obj.W_3[l];

        vec[0].log_min_e = Log(-vec[0].e);
        vec[1].log_min_e = Log(-vec[1].e);
        vec[2].log_min_e = Log(-vec[2].e);
        vec[3].log_min_e = Log(-vec[3].e);

        vec[0].f = data_obj.VR_0[i+l*N]*data_obj.VR_inv_0[l+j*N];
        vec[1].f = data_obj.VR_1[i+l*N]*data_obj.VR_inv_1[l+j*N];
        vec[2].f = data_obj.VR_2[i+l*N]*data_obj.VR_inv_2[l+j*N];
        vec[3].f = data_obj.VR_3[i+l*N]*data_obj.VR_inv_3[l+j*N];

        f_result[i+j*N] += integrate_matrix_element_3D(vec);
        }
        }
        }
        }
        else
      */
      {
        for(int l=0; l<N; l++){

          vec[0].i = 0;
          vec[1].i = 1;
          vec[2].i = 2;
          vec[3].i = 3;

          vec[0].e = data_obj.W_0[l];
          vec[1].e = data_obj.W_1[l];
          vec[2].e = data_obj.W_2[l];
          vec[3].e = data_obj.W_3[l];

          vec[0].log_min_e = Log(-vec[0].e);
          vec[1].log_min_e = Log(-vec[1].e);
          vec[2].log_min_e = Log(-vec[2].e);
          vec[3].log_min_e = Log(-vec[3].e);

          eigenvalue_degeneracy_t degeneracy = find_degeneracy_3D(vec);

          integrate_eigenvalues_3D(degeneracy, vec);

          for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){

              vec[0].f = data_obj.VR_0[i+l*N]*data_obj.VR_inv_0[l+j*N];
              vec[1].f = data_obj.VR_1[i+l*N]*data_obj.VR_inv_1[l+j*N];
              vec[2].f = data_obj.VR_2[i+l*N]*data_obj.VR_inv_2[l+j*N];
              vec[3].f = data_obj.VR_3[i+l*N]*data_obj.VR_inv_3[l+j*N];

              f_result[i+j*N] += vec[0].r*vec[0].f;
              f_result[i+j*N] += vec[1].r*vec[1].f;
              f_result[i+j*N] += vec[2].r*vec[2].f;
              f_result[i+j*N] += vec[3].r*vec[3].f;
            }
          }
        }
      }

      for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
          f_result[i+j*N] *= (6.*volume);
    }
  }

  /*
    template<typename scalartype>
    std::complex<scalartype> tetrahedron_routines_inverse_matrix_function::integrate_matrix_element_3D(std::complex<scalartype>* f,
    std::complex<scalartype>* e)
    {
    const static scalartype VALUE[40] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
    10,11,12,13,14,15,16,17,18,19,
    20,21,22,23,24,25,26,27,28,29,
    30,31,32,33,34,35,36,37,38,39};

    std::complex<scalartype> r[3+1];

    eigenvalue_degeneracy_t degeneracy = find_degeneracy_3D(f, e);

    switch(degeneracy)
    {
    case NO_DEGENERACY:
    r[0] = (-(std::pow(e[0],2)*(e[1] - e[2])*(e[1] - e[3])*(e[2] - e[3])*(VALUE[3]*e[1]*e[2]*e[3] + std::pow(e[0],2)*(e[1] + e[2] + e[3]) - VALUE[2]*e[0]*(e[2]*e[3] + e[1]*(e[2] + e[3])))*std::log(-e[0])) + std::pow(e[1],3)*std::pow(e[0] - e[2],2)*std::pow(e[0] - e[3],2)*(e[2] - e[3])*std::log(-e[1]) +
    (e[0] - e[1])*(-((e[0] - e[1])*std::pow(e[2],3)*std::pow(e[0] - e[3],2)*(e[1] - e[3])*std::log(-e[2])) + (e[0] - e[2])*(e[1] - e[2])*(std::pow(e[0],2)*(e[0] - e[3])*(-e[1] + e[3])*(-e[2] + e[3]) + (e[0] - e[1])*(e[0] - e[2])*std::pow(e[3],3)*std::log(-e[3]))))/(VALUE[6]*std::pow(e[0] - e[1],2)*std::pow(e[0] - e[2],2)*(e[1] - e[2])*std::pow(e[0] - e[3],2)*(e[1] - e[3])*(e[2] - e[3]));
    r[1] = -(-(std::pow(e[0],3)*std::pow(e[1] - e[2],2)*std::pow(e[1] - e[3],2)*(e[2] - e[3])*std::log(-e[0])) + std::pow(e[1],2)*(e[0] - e[2])*(e[0] - e[3])*(e[2] - e[3])*(e[0]*(std::pow(e[1],2) + VALUE[3]*e[2]*e[3] - VALUE[2]*e[1]*(e[2] + e[3])) + e[1]*(-VALUE[2]*e[2]*e[3] + e[1]*(e[2] + e[3])))*std::log(-e[1]) +
    (e[0] - e[1])*((e[0] - e[1])*std::pow(e[2],3)*(e[0] - e[3])*std::pow(e[1] - e[3],2)*std::log(-e[2]) + (e[0] - e[2])*(e[1] - e[2])*(std::pow(e[1],2)*(e[0] - e[3])*(-e[1] + e[3])*(-e[2] + e[3]) - (e[0] - e[1])*(e[1] - e[2])*std::pow(e[3],3)*std::log(-e[3]))))/(VALUE[6]*std::pow(e[0] - e[1],2)*(e[0] - e[2])*std::pow(e[1] - e[2],2)*(e[0] - e[3])*std::pow(e[1] - e[3],2)*(e[2] - e[3]));
    r[2] = -(std::pow(e[0],3)*std::pow(e[1] - e[2],2)*(e[1] - e[3])*std::pow(e[2] - e[3],2)*std::log(-e[0]) - std::pow(e[1],3)*std::pow(e[0] - e[2],2)*(e[0] - e[3])*std::pow(e[2] - e[3],2)*std::log(-e[1]) +
    (e[0] - e[1])*(std::pow(e[2],2)*(e[0] - e[3])*(-e[1] + e[3])*(e[0]*(-VALUE[2]*e[1]*e[2] + std::pow(e[2],2) + VALUE[3]*e[1]*e[3] - VALUE[2]*e[2]*e[3]) + e[2]*(e[1]*(e[2] - VALUE[2]*e[3]) + e[2]*e[3]))*std::log(-e[2]) + (e[0] - e[2])*(e[1] - e[2])*(std::pow(e[2],2)*(e[0] - e[3])*(-e[1] + e[3])*(-e[2] + e[3]) + (e[0] - e[2])*(e[1] - e[2])*std::pow(e[3],3)*std::log(-e[3]))))/(VALUE[6]*(-e[0] + e[1])*std::pow(e[0] - e[2],2)*std::pow(e[1] - e[2],2)*(e[0] - e[3])*(e[1] - e[3])*std::pow(e[2] - e[3],2));
    r[3] = (std::pow(e[0],3)*(e[1] - e[2])*std::pow(e[1] - e[3],2)*std::pow(e[2] - e[3],2)*std::log(-e[0]) - std::pow(e[1],3)*(e[0] - e[2])*std::pow(e[0] - e[3],2)*std::pow(e[2] - e[3],2)*std::log(-e[1]) +
    (e[0] - e[1])*(std::pow(e[2],3)*std::pow(e[0] - e[3],2)*std::pow(e[1] - e[3],2)*std::log(-e[2]) + (e[0] - e[2])*(-e[1] + e[2])*std::pow(e[3],2)*((e[0] - e[3])*(-e[1] + e[3])*(-e[2] + e[3]) + (e[3]*(-VALUE[2]*e[1]*e[2] + e[1]*e[3] + e[2]*e[3]) + e[0]*(VALUE[3]*e[1]*e[2] - VALUE[2]*e[1]*e[3] - VALUE[2]*e[2]*e[3] + std::pow(e[3],2)))*std::log(-e[3]))))/(VALUE[6]*(e[0] - e[1])*(e[0] - e[2])*(e[1] - e[2])*std::pow(e[0] - e[3],2)*std::pow(e[1] - e[3],2)*std::pow(e[2] - e[3],2));
    break;

    case TWOFOLD_DEGENERACY:
    r[0] = (std::pow(e[0],2)*std::pow(e[1] - e[2],2)*(VALUE[3]*e[1]*e[2] - e[0]*(e[1] + VALUE[2]*e[2]))*std::log(-e[0]) + std::pow(e[1],3)*std::pow(e[0] - e[2],3)*std::log(-e[1]) + (e[0] - e[1])*((e[0] - e[2])*(-e[1] + e[2])*(e[0]*std::pow(e[2],2) - e[1]*std::pow(e[2],2) + std::pow(e[0],2)*(-e[1] + e[2])) - (e[0] - e[1])*std::pow(e[2],2)*(VALUE[3]*e[0]*e[1] - VALUE[2]*e[0]*e[2] - e[1]*e[2])*std::log(-e[2])))/(VALUE[6]*std::pow(e[0] - e[1],2)*std::pow(e[0] - e[2],3)*std::pow(e[1] - e[2],2));
    r[1] = (-(std::pow(e[0],3)*std::pow(e[1] - e[2],3)*std::log(-e[0])) + std::pow(e[1],2)*std::pow(e[0] - e[2],2)*(e[0]*(e[1] - VALUE[3]*e[2]) + VALUE[2]*e[1]*e[2])*std::log(-e[1]) + (e[0] - e[1])*((e[0] - e[2])*(e[1] - e[2])*(-(e[1]*e[2]*(e[1] + e[2])) + e[0]*(std::pow(e[1],2) + std::pow(e[2],2))) + (e[0] - e[1])*std::pow(e[2],2)*(VALUE[3]*e[0]*e[1] - e[0]*e[2] - VALUE[2]*e[1]*e[2])*std::log(-e[2])))/(VALUE[6]*std::pow(e[0] - e[1],2)*std::pow(e[0] - e[2],2)*std::pow(e[1] - e[2],2)*(-e[1] + e[2]));
    r[2] = -(VALUE[2]*std::pow(e[0],3)*std::pow(e[1] - e[2],3)*std::log(-e[0]) - VALUE[2]*std::pow(e[1],3)*std::pow(e[0] - e[2],3)*std::log(-e[1]) + (e[0] - e[1])*e[2]*((e[0] - e[2])*(e[1] - e[2])*(VALUE[5]*e[0]*e[1] - VALUE[3]*e[0]*e[2] - VALUE[3]*e[1]*e[2] + std::pow(e[2],2)) + VALUE[2]*(std::pow(e[1],2)*std::pow(e[2],2) + e[0]*e[1]*e[2]*(-VALUE[3]*e[1] + e[2]) + std::pow(e[0],2)*(VALUE[3]*std::pow(e[1],2) - VALUE[3]*e[1]*e[2] + std::pow(e[2],2)))*std::log(-e[2])))/(VALUE[12]*(-e[0] + e[1])*std::pow(e[0] - e[2],3)*std::pow(e[1] - e[2],3));
    r[3] = -(VALUE[2]*std::pow(e[0],3)*std::pow(e[1] - e[2],3)*std::log(-e[0]) - VALUE[2]*std::pow(e[1],3)*std::pow(e[0] - e[2],3)*std::log(-e[1]) + (e[0] - e[1])*e[2]*((e[0] - e[2])*(e[1] - e[2])*(5.*e[0]*e[1] - VALUE[3]*e[0]*e[2] - VALUE[3]*e[1]*e[2] + std::pow(e[2],2)) + VALUE[2]*(std::pow(e[1],2)*std::pow(e[2],2) + e[0]*e[1]*e[2]*(-VALUE[3]*e[1] + e[2]) + std::pow(e[0],2)*(VALUE[3]*std::pow(e[1],2) - VALUE[3]*e[1]*e[2] + std::pow(e[2],2)))*std::log(-e[2])))/(VALUE[12]*(-e[0] + e[1])*std::pow(e[0] - e[2],3)*std::pow(e[1] - e[2],3));
    break;

    case THREEFOLD_DEGENERACY:
    r[0] = (VALUE[2]*std::pow(e[0],3) + VALUE[3]*std::pow(e[0],2)*e[1] - VALUE[6]*e[0]*std::pow(e[1],2) + std::pow(e[1],3) - VALUE[6]*std::pow(e[0],2)*e[1]*std::log(-e[0]) + VALUE[6]*std::pow(e[0],2)*e[1]*std::log(-e[1]))/(VALUE[12]*std::pow(e[0] - e[1],4));
    r[1] = (-VALUE[11]*std::pow(e[0],3) + VALUE[18]*std::pow(e[0],2)*e[1] - VALUE[9]*e[0]*std::pow(e[1],2) + VALUE[2]*std::pow(e[1],3) + VALUE[6]*std::pow(e[0],3)*std::log(-e[0]) - VALUE[6]*std::pow(e[0],3)*std::log(-e[1]))/(VALUE[36]*std::pow(e[0] - e[1],4));
    r[2] = (-VALUE[11]*std::pow(e[0],3) + VALUE[18]*std::pow(e[0],2)*e[1] - VALUE[9]*e[0]*std::pow(e[1],2) + VALUE[2]*std::pow(e[1],3) + 6.*std::pow(e[0],3)*std::log(-e[0]) - 6.*std::pow(e[0],3)*std::log(-e[1]))/(VALUE[36]*std::pow(e[0] - e[1],4));
    r[3] = (-VALUE[11]*std::pow(e[0],3) + VALUE[18]*std::pow(e[0],2)*e[1] - VALUE[9]*e[0]*std::pow(e[1],2) + VALUE[2]*std::pow(e[1],3) + 6.*std::pow(e[0],3)*std::log(-e[0]) - 6.*std::pow(e[0],3)*std::log(-e[1]))/(VALUE[36]*std::pow(e[0] - e[1],4));
    break;

    case FOURFOLD_DEGENERACY:
    r[0] = VALUE[1]/(VALUE[24]*e[0]);
    r[1] = VALUE[1]/(VALUE[24]*e[0]);
    r[2] = VALUE[1]/(VALUE[24]*e[0]);
    r[3] = VALUE[1]/(VALUE[24]*e[0]);
    break;

    default:
    throw std::logic_error(__FUNCTION__);
    }

    return f[0]*r[0]+f[1]*r[1]+f[2]*r[2]+f[3]*r[3];
    }
  */

  /*
    template<typename scalartype>
    std::complex<scalartype> tetrahedron_routines_inverse_matrix_function::integrate_matrix_element_3D(std::complex<scalartype>* f,
    std::complex<scalartype>* e)
    {
    std::complex<scalartype> r[3+1];

    eigenvalue_degeneracy_t degeneracy = find_degeneracy_3D(f, e);

    std::complex<scalartype> e0 = e[0];
    std::complex<scalartype> e1 = e[1];
    std::complex<scalartype> e2 = e[2];
    std::complex<scalartype> e3 = e[3];

    switch(degeneracy)
    {
    case NO_DEGENERACY:
    {
    assert(abs(e0-e1)>EPSILON());
    assert(abs(e1-e2)>EPSILON());
    assert(abs(e3-e2)>EPSILON());
    assert(abs(e3-e0)>EPSILON());
    assert(abs(e0-e2)>EPSILON());
    assert(abs(e1-e3)>EPSILON());

    r[0] = (-((Power(e0,2)*(3.*e1*e2*e3 + Power(e0,2)*(e1 + e2 + e3) - 2.*e0*(e1*e2 + (e1 + e2)*e3))*Log(-e0))/(Power(e0 - e1,2)*Power(e0 - e2,2)*Power(e0 - e3,2))) + (Power(e1,3)*Log(-e1))/(Power(e0 - e1,2)*(e1 - e2)*(e1 - e3)) + (Power(e2,3)*Log(-e2))/(Power(e0 - e2,2)*(-e1 + e2)*(e2 - e3)) + ((Power(e0,2)*(e0 - e3))/((e0 - e1)*(e0 - e2)) + (Power(e3,3)*Log(-e3))/((-e1 + e3)*(-e2 + e3)))/Power(e0 - e3,2))/6.;
    r[1] = (Power(e1,2)/((-e0 + e1)*(e1 - e2)*(e1 - e3)) + (Power(e0,3)*Log(-e0))/(Power(e0 - e1,2)*(e0 - e2)*(e0 - e3)) - (Power(e1,2)*(e0*(Power(e1,2) + 3.*e2*e3 - 2.*e1*(e2 + e3)) + e1*(-2.*e2*e3 + e1*(e2 + e3)))*Log(-e1))/(Power(e0 - e1,2)*Power(e1 - e2,2)*Power(e1 - e3,2)) + ((Power(e2,3)*Log(-e2))/(Power(e1 - e2,2)*(-e0 + e2)) + (Power(e3,3)*Log(-e3))/((e0 - e3)*Power(e1 - e3,2)))/(e2 - e3))/6.;
    r[2] = (Power(e0,3)*Log(-e0) + (-(Power(e1,3)*Power(e0 - e2,2)*(e0 - e3)*Power(e2 - e3,2)*Log(-e1)) + (e0 - e1)*(Power(e2,2)*(e0 - e3)*(-e1 + e3)*(e0*(-2.*e1*e2 + Power(e2,2) + 3.*e1*e3 - 2.*e2*e3) + e2*(e1*e2 - 2.*e1*e3 + e2*e3))*Log(-e2) + (-e0 + e2)*(-e1 + e2)*(Power(e2,2)*(e0 - e3)*(-e1 + e3)*(-e2 + e3) + (e0 - e2)*(e1 - e2)*Power(e3,3)*Log(-e3))))/(Power(e1 - e2,2)*(e1 - e3)*Power(e2 - e3,2)))/(6.*(e0 - e1)*Power(e0 - e2,2)*(e0 - e3));
    r[3] = (Power(e0,3)*Log(-e0) + (-(Power(e1,3)*(e0 - e2)*Power(e0 - e3,2)*Power(e2 - e3,2)*Log(-e1)) + (e0 - e1)*(Power(e2,3)*Power(e0 - e3,2)*Power(e1 - e3,2)*Log(-e2) + (e0 - e2)*(-e1 + e2)*Power(e3,2)*((e0 - e3)*(-e1 + e3)*(-e2 + e3) + (3.*e0*e1*e2 - 2.*(e0*e1 + (e0 + e1)*e2)*e3 + (e0 + e1 + e2)*Power(e3,2))*Log(-e3))))/((e1 - e2)*Power(e1 - e3,2)*Power(e2 - e3,2)))/(6.*(e0 - e1)*(e0 - e2)*Power(e0 - e3,2));
    }
    break;

    case TWOFOLD_DEGENERACY:
    {
    assert(abs(e0-e1)>EPSILON());
    assert(abs(e1-e2)>EPSILON());
    assert(abs(e2-e0)>EPSILON());
    assert(abs(e2-e3)<EPSILON());

    r[0] = ((Power(e0,2)*(e1 - e2) - e0*Power(e2,2) + e1*Power(e2,2))/((e0 - e1)*Power(e0 - e2,2)*(e1 - e2)) - (Power(e0,2)*(-3.*e1*e2 + e0*(e1 + 2.*e2))*Log(-e0))/(Power(e0 - e1,2)*Power(e0 - e2,3)) + ((Power(e1,3)*Log(-e1))/Power(e0 - e1,2) + (Power(e2,2)*(-3.*e0*e1 + 2.*e0*e2 + e1*e2)*Log(-e2))/Power(e0 - e2,3))/Power(e1 - e2,2))/6.;
    r[1] = (Power(e1,2)/((-e0 + e1)*Power(e1 - e2,2)) + (Power(e0,3)*Log(-e0))/(Power(e0 - e1,2)*Power(e0 - e2,2)) - (Power(e1,2)*(e0*(e1 - 3.*e2) + 2.*e1*e2)*Log(-e1))/(Power(e0 - e1,2)*Power(e1 - e2,3)) + (Power(e2,2)*((e0 - e2)*(e1 - e2) + (3.*e0*e1 - (e0 + 2.*e1)*e2)*Log(-e2)))/(Power(e0 - e2,2)*Power(-e1 + e2,3)))/6.;
    r[2] = -(2.*Power(e0,3)*Power(e1 - e2,3)*Log(-e0) - 2.*Power(e1,3)*Power(e0 - e2,3)*Log(-e1) + (e0 - e1)*e2*((e0 - e2)*(e1 - e2)*(5.*e0*e1 - 3.*(e0 + e1)*e2 + Power(e2,2)) + 2.*(3.*Power(e0,2)*Power(e1,2) - 3.*e0*e1*(e0 + e1)*e2 + (Power(e0,2) + e0*e1 + Power(e1,2))*Power(e2,2))*Log(-e2)))/ (12.*(e0 - e1)*Power(e0 - e2,3)*Power(-e1 + e2,3));
    r[3] = -(2.*Power(e0,3)*Power(e1 - e2,3)*Log(-e0) - 2.*Power(e1,3)*Power(e0 - e2,3)*Log(-e1) + (e0 - e1)*e2*((e0 - e2)*(e1 - e2)*(5.*e0*e1 - 3.*(e0 + e1)*e2 + Power(e2,2)) + 2.*(3.*Power(e0,2)*Power(e1,2) - 3.*e0*e1*(e0 + e1)*e2 + (Power(e0,2) + e0*e1 + Power(e1,2))*Power(e2,2))*Log(-e2)))/ (12.*(e0 - e1)*Power(e0 - e2,3)*Power(-e1 + e2,3));
    }
    break;

    case THREEFOLD_DEGENERACY_A:
    {
    assert(abs(e0-e1)>EPSILON());
    assert(abs(e1-e2)<EPSILON());
    assert(abs(e2-e3)<EPSILON());

    r[0] = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-Log(-e0) + Log(-e1)))/(12.*Power(e0 - e1,4));
    r[1] = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(Log(-e0) - Log(-e1)))/(36.*Power(e0 - e1,4));
    r[2] = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(Log(-e0) - Log(-e1)))/(36.*Power(e0 - e1,4));
    r[3] = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(Log(-e0) - Log(-e1)))/(36.*Power(e0 - e1,4));
    }
    break;

    case THREEFOLD_DEGENERACY_B:
    {
    assert(abs(e0-e1)>EPSILON());
    assert(abs(e0-e2)<EPSILON());
    assert(abs(e1-e3)<EPSILON());

    r[0] = ((e0 - e1)*(Power(e0,2) - 5.*e0*e1 - 2.*Power(e1,2)) + 6.*e0*Power(e1,2)*(Log(-e0) - Log(-e1)))/(12.*Power(e0 - e1,4));
    r[1] = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-Log(-e0) + Log(-e1)))/(12.*Power(e0 - e1,4));
    r[2] = ((e0 - e1)*(Power(e0,2) - 5.*e0*e1 - 2.*Power(e1,2)) + 6.*e0*Power(e1,2)*(Log(-e0) - Log(-e1)))/(12.*Power(e0 - e1,4));
    r[3] = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-Log(-e0) + Log(-e1)))/(12.*Power(e0 - e1,4));
    }
    break;

    case FOURFOLD_DEGENERACY:
    {
    assert(abs(e0-e1)<EPSILON());
    assert(abs(e1-e2)<EPSILON());
    assert(abs(e2-e3)<EPSILON());

    r[0] = 1./(24.*e0);
    r[1] = 1./(24.*e0);
    r[2] = 1./(24.*e0);
    r[3] = 1./(24.*e0);
    }
    break;

    default:
    throw std::logic_error(__FUNCTION__);
    }

    return f[0]*r[0]+f[1]*r[1]+f[2]*r[2]+f[3]*r[3];
    }

    template<typename scalartype>
    eigenvalue_degeneracy_t tetrahedron_routines_inverse_matrix_function::find_degeneracy_3D(std::complex<scalartype>* f,
    std::complex<scalartype>* e)
    {
    std::vector<matrix_element_struct<scalartype> > vec(4);

    vec[0].i = 0; vec[0].e = e[0]; vec[0].f = f[0];
    vec[1].i = 1; vec[1].e = e[1]; vec[1].f = f[1];
    vec[2].i = 2; vec[2].e = e[2]; vec[2].f = f[2];
    vec[3].i = 3; vec[3].e = e[3]; vec[3].f = f[3];

    if(abs(vec[1].e-vec[0].e)<EPSILON() and
    abs(vec[2].e-vec[1].e)<EPSILON() and
    abs(vec[3].e-vec[2].e)<EPSILON())
    {
    //         cout << vec[0].e << "\t";
    //         cout << vec[1].e << "\t";
    //         cout << vec[2].e << "\t";
    //         cout << vec[3].e << "\n";
    //         cout << "\n";
    //         cout << "\tFOURFOLD_DEGENERACY\n";

    return FOURFOLD_DEGENERACY;
    }

    do
    {
    //         cout << vec[0].e << "\t";
    //         cout << vec[1].e << "\t";
    //         cout << vec[2].e << "\t";
    //         cout << vec[3].e << "\n";
    //         //cout << "\n";

    if(abs(vec[1].e-vec[0].e)>EPSILON() and
    abs(vec[2].e-vec[1].e)>EPSILON() and
    abs(vec[0].e-vec[2].e)>EPSILON() and
    abs(vec[3].e-vec[2].e)<EPSILON())
    {
    e[0] = vec[0].e; f[0] = vec[0].f;
    e[1] = vec[1].e; f[1] = vec[1].f;
    e[2] = vec[2].e; f[2] = vec[2].f;
    e[3] = vec[3].e; f[3] = vec[3].f;

    //      cout << "\n\tTWOFOLD_DEGENERACY\n";

    return TWOFOLD_DEGENERACY;
    }

    if(abs(vec[1].e-vec[0].e)>EPSILON() and
    abs(vec[2].e-vec[1].e)<EPSILON() and
    abs(vec[3].e-vec[1].e)<EPSILON())
    {
    e[0] = vec[0].e; f[0] = vec[0].f;
    e[1] = vec[1].e; f[1] = vec[1].f;
    e[2] = vec[2].e; f[2] = vec[2].f;
    e[3] = vec[3].e; f[3] = vec[3].f;

    //      cout << "\n\tTHREEFOLD_DEGENERACY_A\n";

    return THREEFOLD_DEGENERACY_A;
    }

    if(abs(vec[1].e-vec[0].e)>EPSILON() and
    abs(vec[2].e-vec[0].e)<EPSILON() and
    abs(vec[3].e-vec[1].e)<EPSILON())
    {
    e[0] = vec[0].e; f[0] = vec[0].f;
    e[1] = vec[1].e; f[1] = vec[1].f;
    e[2] = vec[2].e; f[2] = vec[2].f;
    e[3] = vec[3].e; f[3] = vec[3].f;

    //      cout << "\n\tTHREEFOLD_DEGENERACY_B\n";

    return THREEFOLD_DEGENERACY_B;
    }
    }
    while
    ( std::next_permutation(vec.begin(), vec.end(), tetrahedron_routines_inverse_matrix_function::permutation_comp<scalartype>) );

    //     cout << "\n\tNO_DEGENERACY\n";

    return NO_DEGENERACY;
    }
  */

  /*
    template<typename scalartype>
    std::complex<scalartype> tetrahedron_routines_inverse_matrix_function::integrate_matrix_element_3D(std::vector<matrix_element_struct<scalartype> >& vec)
    {
    eigenvalue_degeneracy_t degeneracy = find_degeneracy_3D(vec);

    std::complex<scalartype> r[4];
    std::complex<scalartype> f[4];

    f[0] = vec[0].f;
    f[1] = vec[1].f;
    f[2] = vec[2].f;
    f[3] = vec[3].f;

    std::complex<scalartype> e0 = vec[0].e;
    std::complex<scalartype> e1 = vec[1].e;
    std::complex<scalartype> e2 = vec[2].e;
    std::complex<scalartype> e3 = vec[3].e;

    switch(degeneracy)
    {
    case NO_DEGENERACY:
    {
    assert(abs(e0-e1)>EPSILON());
    assert(abs(e1-e2)>EPSILON());
    assert(abs(e2-e3)>EPSILON());
    assert(abs(e3-e0)>EPSILON());
    assert(abs(e0-e2)>EPSILON());
    assert(abs(e1-e3)>EPSILON());

    r[0] = (-((Power(e0,2)*(3.*e1*e2*e3 + Power(e0,2)*(e1 + e2 + e3) - 2.*e0*(e1*e2 + (e1 + e2)*e3))*Log(-e0))/(Power(e0 - e1,2)*Power(e0 - e2,2)*Power(e0 - e3,2))) + (Power(e1,3)*Log(-e1))/(Power(e0 - e1,2)*(e1 - e2)*(e1 - e3)) + (Power(e2,3)*Log(-e2))/(Power(e0 - e2,2)*(-e1 + e2)*(e2 - e3)) + ((Power(e0,2)*(e0 - e3))/((e0 - e1)*(e0 - e2)) + (Power(e3,3)*Log(-e3))/((-e1 + e3)*(-e2 + e3)))/Power(e0 - e3,2))/6.;
    r[1] = (Power(e1,2)/((-e0 + e1)*(e1 - e2)*(e1 - e3)) + (Power(e0,3)*Log(-e0))/(Power(e0 - e1,2)*(e0 - e2)*(e0 - e3)) - (Power(e1,2)*(e0*(Power(e1,2) + 3.*e2*e3 - 2.*e1*(e2 + e3)) + e1*(-2.*e2*e3 + e1*(e2 + e3)))*Log(-e1))/(Power(e0 - e1,2)*Power(e1 - e2,2)*Power(e1 - e3,2)) + ((Power(e2,3)*Log(-e2))/(Power(e1 - e2,2)*(-e0 + e2)) + (Power(e3,3)*Log(-e3))/((e0 - e3)*Power(e1 - e3,2)))/(e2 - e3))/6.;
    r[2] = (Power(e0,3)*Log(-e0) + (-(Power(e1,3)*Power(e0 - e2,2)*(e0 - e3)*Power(e2 - e3,2)*Log(-e1)) + (e0 - e1)*(Power(e2,2)*(e0 - e3)*(-e1 + e3)*(e0*(-2.*e1*e2 + Power(e2,2) + 3.*e1*e3 - 2.*e2*e3) + e2*(e1*e2 - 2.*e1*e3 + e2*e3))*Log(-e2) + (-e0 + e2)*(-e1 + e2)*(Power(e2,2)*(e0 - e3)*(-e1 + e3)*(-e2 + e3) + (e0 - e2)*(e1 - e2)*Power(e3,3)*Log(-e3))))/(Power(e1 - e2,2)*(e1 - e3)*Power(e2 - e3,2)))/(6.*(e0 - e1)*Power(e0 - e2,2)*(e0 - e3));
    r[3] = (Power(e0,3)*Log(-e0) + (-(Power(e1,3)*(e0 - e2)*Power(e0 - e3,2)*Power(e2 - e3,2)*Log(-e1)) + (e0 - e1)*(Power(e2,3)*Power(e0 - e3,2)*Power(e1 - e3,2)*Log(-e2) + (e0 - e2)*(-e1 + e2)*Power(e3,2)*((e0 - e3)*(-e1 + e3)*(-e2 + e3) + (3.*e0*e1*e2 - 2.*(e0*e1 + (e0 + e1)*e2)*e3 + (e0 + e1 + e2)*Power(e3,2))*Log(-e3))))/((e1 - e2)*Power(e1 - e3,2)*Power(e2 - e3,2)))/(6.*(e0 - e1)*(e0 - e2)*Power(e0 - e3,2));
    }
    break;

    case TWOFOLD_DEGENERACY:
    {
    assert(abs(e0-e1)>EPSILON());
    assert(abs(e1-e2)>EPSILON());
    assert(abs(e2-e0)>EPSILON());
    assert(abs(e2-e3)<EPSILON());

    r[0] = ((Power(e0,2)*(e1 - e2) - e0*Power(e2,2) + e1*Power(e2,2))/((e0 - e1)*Power(e0 - e2,2)*(e1 - e2)) - (Power(e0,2)*(-3.*e1*e2 + e0*(e1 + 2.*e2))*Log(-e0))/(Power(e0 - e1,2)*Power(e0 - e2,3)) + ((Power(e1,3)*Log(-e1))/Power(e0 - e1,2) + (Power(e2,2)*(-3.*e0*e1 + 2.*e0*e2 + e1*e2)*Log(-e2))/Power(e0 - e2,3))/Power(e1 - e2,2))/6.;
    r[1] = (Power(e1,2)/((-e0 + e1)*Power(e1 - e2,2)) + (Power(e0,3)*Log(-e0))/(Power(e0 - e1,2)*Power(e0 - e2,2)) - (Power(e1,2)*(e0*(e1 - 3.*e2) + 2.*e1*e2)*Log(-e1))/(Power(e0 - e1,2)*Power(e1 - e2,3)) + (Power(e2,2)*((e0 - e2)*(e1 - e2) + (3.*e0*e1 - (e0 + 2.*e1)*e2)*Log(-e2)))/(Power(e0 - e2,2)*Power(-e1 + e2,3)))/6.;
    r[2] = -(2.*Power(e0,3)*Power(e1 - e2,3)*Log(-e0) - 2.*Power(e1,3)*Power(e0 - e2,3)*Log(-e1) + (e0 - e1)*e2*((e0 - e2)*(e1 - e2)*(5.*e0*e1 - 3.*(e0 + e1)*e2 + Power(e2,2)) + 2.*(3.*Power(e0,2)*Power(e1,2) - 3.*e0*e1*(e0 + e1)*e2 + (Power(e0,2) + e0*e1 + Power(e1,2))*Power(e2,2))*Log(-e2)))/ (12.*(e0 - e1)*Power(e0 - e2,3)*Power(-e1 + e2,3));
    r[3] = -(2.*Power(e0,3)*Power(e1 - e2,3)*Log(-e0) - 2.*Power(e1,3)*Power(e0 - e2,3)*Log(-e1) + (e0 - e1)*e2*((e0 - e2)*(e1 - e2)*(5.*e0*e1 - 3.*(e0 + e1)*e2 + Power(e2,2)) + 2.*(3.*Power(e0,2)*Power(e1,2) - 3.*e0*e1*(e0 + e1)*e2 + (Power(e0,2) + e0*e1 + Power(e1,2))*Power(e2,2))*Log(-e2)))/ (12.*(e0 - e1)*Power(e0 - e2,3)*Power(-e1 + e2,3));
    }
    break;

    case THREEFOLD_DEGENERACY_A:
    {
    assert(abs(e0-e1)>EPSILON());
    assert(abs(e2-e1)<EPSILON());
    assert(abs(e3-e1)<EPSILON());

    r[0] = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-Log(-e0) + Log(-e1)))/(12.*Power(e0 - e1,4));
    r[1] = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(Log(-e0) - Log(-e1)))/(36.*Power(e0 - e1,4));
    r[2] = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(Log(-e0) - Log(-e1)))/(36.*Power(e0 - e1,4));
    r[3] = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(Log(-e0) - Log(-e1)))/(36.*Power(e0 - e1,4));
    }
    break;

    case THREEFOLD_DEGENERACY_B:
    {
    assert(abs(e0-e1)>EPSILON());
    assert(abs(e0-e2)<EPSILON());
    assert(abs(e1-e3)<EPSILON());

    r[0] = ((e0 - e1)*(Power(e0,2) - 5.*e0*e1 - 2.*Power(e1,2)) + 6.*e0*Power(e1,2)*(Log(-e0) - Log(-e1)))/(12.*Power(e0 - e1,4));
    r[1] = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-Log(-e0) + Log(-e1)))/(12.*Power(e0 - e1,4));
    r[2] = ((e0 - e1)*(Power(e0,2) - 5.*e0*e1 - 2.*Power(e1,2)) + 6.*e0*Power(e1,2)*(Log(-e0) - Log(-e1)))/(12.*Power(e0 - e1,4));
    r[3] = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-Log(-e0) + Log(-e1)))/(12.*Power(e0 - e1,4));
    }
    break;

    case FOURFOLD_DEGENERACY:
    {
    assert(abs(e0-e1)<EPSILON());
    assert(abs(e1-e2)<EPSILON());
    assert(abs(e2-e3)<EPSILON());

    r[0] = 1./(24.*e0);
    r[1] = 1./(24.*e0);
    r[2] = 1./(24.*e0);
    r[3] = 1./(24.*e0);
    }
    break;

    default:
    throw std::logic_error(__FUNCTION__);
    }

    return f[0]*r[0]+f[1]*r[1]+f[2]*r[2]+f[3]*r[3];
    }

    template<typename scalartype>
    eigenvalue_degeneracy_t tetrahedron_routines_inverse_matrix_function::find_degeneracy_3D(std::vector<matrix_element_struct<scalartype> >& vec)
    {
    if(abs(vec[1].e-vec[0].e)>EPSILON() and
    abs(vec[2].e-vec[1].e)>EPSILON() and
    abs(vec[3].e-vec[2].e)>EPSILON() and
    abs(vec[0].e-vec[3].e)>EPSILON() and
    abs(vec[3].e-vec[1].e)>EPSILON() and
    abs(vec[2].e-vec[0].e)>EPSILON())
    {
    return NO_DEGENERACY;
    }

    if(abs(vec[1].e-vec[0].e)<EPSILON() and
    abs(vec[2].e-vec[1].e)<EPSILON() and
    abs(vec[3].e-vec[2].e)<EPSILON())
    {
    return FOURFOLD_DEGENERACY;
    }

    vec[0].i = 0; //vec[0].e = e[0]; vec[0].f = f[0];
    vec[1].i = 1; //vec[1].e = e[1]; vec[1].f = f[1];
    vec[2].i = 2; //vec[2].e = e[2]; vec[2].f = f[2];
    vec[3].i = 3; //vec[3].e = e[3]; vec[3].f = f[3];

    do
    {
    if(abs(vec[1].e-vec[0].e)>EPSILON() and
    abs(vec[2].e-vec[1].e)>EPSILON() and
    abs(vec[0].e-vec[2].e)>EPSILON() and
    abs(vec[3].e-vec[2].e)<EPSILON())
    {
    return TWOFOLD_DEGENERACY;
    }

    if(abs(vec[1].e-vec[0].e)>EPSILON() and
    abs(vec[2].e-vec[1].e)<EPSILON() and
    abs(vec[3].e-vec[1].e)<EPSILON())
    {
    return THREEFOLD_DEGENERACY_A;
    }

    if(abs(vec[1].e-vec[0].e)>EPSILON() and
    abs(vec[2].e-vec[0].e)<EPSILON() and
    abs(vec[3].e-vec[1].e)<EPSILON())
    {
    return THREEFOLD_DEGENERACY_B;
    }
    }
    while
    ( std::next_permutation(vec.begin(), vec.end(), tetrahedron_routines_inverse_matrix_function::permutation_comp<scalartype>) );

    throw std::logic_error(__FUNCTION__);
    return NO_DEGENERACY;
    }
  */

  template<typename scalartype>
  std::complex<scalartype> tetrahedron_routines_inverse_matrix_function::integrate_matrix_element_3D(std::vector<matrix_element_struct<scalartype> >& vec)
  {
    assert(vec.size()==4);

    eigenvalue_degeneracy_t degeneracy = find_degeneracy_3D(vec);

    std::complex<scalartype> r[4];
    std::complex<scalartype> f[4];

    f[0] = vec[0].f;
    f[1] = vec[1].f;
    f[2] = vec[2].f;
    f[3] = vec[3].f;

    std::complex<scalartype> e0 = vec[0].e;
    std::complex<scalartype> e1 = vec[1].e;
    std::complex<scalartype> e2 = vec[2].e;
    std::complex<scalartype> e3 = vec[3].e;

    std::complex<scalartype> log_min_e0 = vec[0].log_min_e;
    std::complex<scalartype> log_min_e1 = vec[1].log_min_e;
    std::complex<scalartype> log_min_e2 = vec[2].log_min_e;
    std::complex<scalartype> log_min_e3 = vec[3].log_min_e;

    switch(degeneracy)
      {
      case NO_DEGENERACY:
        {
          assert(not_equal(e0, e1));
          assert(not_equal(e1, e2));
          assert(not_equal(e2, e3));
          assert(not_equal(e3, e0));
          assert(not_equal(e0, e2));
          assert(not_equal(e1, e3));

          /*
            r[0] = (-((Power(e0,2)*(3.*e1*e2*e3 + Power(e0,2)*(e1 + e2 + e3) - 2.*e0*(e1*e2 + (e1 + e2)*e3))*Log(-e0))/(Power(e0 - e1,2)*Power(e0 - e2,2)*Power(e0 - e3,2))) + (Power(e1,3)*Log(-e1))/(Power(e0 - e1,2)*(e1 - e2)*(e1 - e3)) + (Power(e2,3)*Log(-e2))/(Power(e0 - e2,2)*(-e1 + e2)*(e2 - e3)) + ((Power(e0,2)*(e0 - e3))/((e0 - e1)*(e0 - e2)) + (Power(e3,3)*Log(-e3))/((-e1 + e3)*(-e2 + e3)))/Power(e0 - e3,2))/6.;
            r[1] = (Power(e1,2)/((-e0 + e1)*(e1 - e2)*(e1 - e3)) + (Power(e0,3)*Log(-e0))/(Power(e0 - e1,2)*(e0 - e2)*(e0 - e3)) - (Power(e1,2)*(e0*(Power(e1,2) + 3.*e2*e3 - 2.*e1*(e2 + e3)) + e1*(-2.*e2*e3 + e1*(e2 + e3)))*Log(-e1))/(Power(e0 - e1,2)*Power(e1 - e2,2)*Power(e1 - e3,2)) + ((Power(e2,3)*Log(-e2))/(Power(e1 - e2,2)*(-e0 + e2)) + (Power(e3,3)*Log(-e3))/((e0 - e3)*Power(e1 - e3,2)))/(e2 - e3))/6.;
            r[2] = (Power(e0,3)*Log(-e0) + (-(Power(e1,3)*Power(e0 - e2,2)*(e0 - e3)*Power(e2 - e3,2)*Log(-e1)) + (e0 - e1)*(Power(e2,2)*(e0 - e3)*(-e1 + e3)*(e0*(-2.*e1*e2 + Power(e2,2) + 3.*e1*e3 - 2.*e2*e3) + e2*(e1*e2 - 2.*e1*e3 + e2*e3))*Log(-e2) + (-e0 + e2)*(-e1 + e2)*(Power(e2,2)*(e0 - e3)*(-e1 + e3)*(-e2 + e3) + (e0 - e2)*(e1 - e2)*Power(e3,3)*Log(-e3))))/(Power(e1 - e2,2)*(e1 - e3)*Power(e2 - e3,2)))/(6.*(e0 - e1)*Power(e0 - e2,2)*(e0 - e3));
            r[3] = (Power(e0,3)*Log(-e0) + (-(Power(e1,3)*(e0 - e2)*Power(e0 - e3,2)*Power(e2 - e3,2)*Log(-e1)) + (e0 - e1)*(Power(e2,3)*Power(e0 - e3,2)*Power(e1 - e3,2)*Log(-e2) + (e0 - e2)*(-e1 + e2)*Power(e3,2)*((e0 - e3)*(-e1 + e3)*(-e2 + e3) + (3.*e0*e1*e2 - 2.*(e0*e1 + (e0 + e1)*e2)*e3 + (e0 + e1 + e2)*Power(e3,2))*Log(-e3))))/((e1 - e2)*Power(e1 - e3,2)*Power(e2 - e3,2)))/(6.*(e0 - e1)*(e0 - e2)*Power(e0 - e3,2));
          */

          r[0] = (-((Power(e0,2)*(3.*e1*e2*e3 + Power(e0,2)*(e1 + e2 + e3) - 2.*e0*(e1*e2 + (e1 + e2)*e3))*log_min_e0)/(Power(e0 - e1,2)*Power(e0 - e2,2)*Power(e0 - e3,2))) + (Power(e1,3)*log_min_e1)/(Power(e0 - e1,2)*(e1 - e2)*(e1 - e3)) + (Power(e2,3)*log_min_e2)/(Power(e0 - e2,2)*(-e1 + e2)*(e2 - e3)) + ((Power(e0,2)*(e0 - e3))/((e0 - e1)*(e0 - e2)) + (Power(e3,3)*log_min_e3)/((-e1 + e3)*(-e2 + e3)))/Power(e0 - e3,2))/6.;
          r[1] = (Power(e1,2)/((-e0 + e1)*(e1 - e2)*(e1 - e3)) + (Power(e0,3)*log_min_e0)/(Power(e0 - e1,2)*(e0 - e2)*(e0 - e3)) - (Power(e1,2)*(e0*(Power(e1,2) + 3.*e2*e3 - 2.*e1*(e2 + e3)) + e1*(-2.*e2*e3 + e1*(e2 + e3)))*log_min_e1)/(Power(e0 - e1,2)*Power(e1 - e2,2)*Power(e1 - e3,2)) + ((Power(e2,3)*log_min_e2)/(Power(e1 - e2,2)*(-e0 + e2)) + (Power(e3,3)*log_min_e3)/((e0 - e3)*Power(e1 - e3,2)))/(e2 - e3))/6.;
          r[2] = (Power(e0,3)*log_min_e0 + (-(Power(e1,3)*Power(e0 - e2,2)*(e0 - e3)*Power(e2 - e3,2)*log_min_e1) + (e0 - e1)*(Power(e2,2)*(e0 - e3)*(-e1 + e3)*(e0*(-2.*e1*e2 + Power(e2,2) + 3.*e1*e3 - 2.*e2*e3) + e2*(e1*e2 - 2.*e1*e3 + e2*e3))*log_min_e2 + (-e0 + e2)*(-e1 + e2)*(Power(e2,2)*(e0 - e3)*(-e1 + e3)*(-e2 + e3) + (e0 - e2)*(e1 - e2)*Power(e3,3)*log_min_e3)))/(Power(e1 - e2,2)*(e1 - e3)*Power(e2 - e3,2)))/(6.*(e0 - e1)*Power(e0 - e2,2)*(e0 - e3));
          r[3] = (Power(e0,3)*log_min_e0 + (-(Power(e1,3)*(e0 - e2)*Power(e0 - e3,2)*Power(e2 - e3,2)*log_min_e1) + (e0 - e1)*(Power(e2,3)*Power(e0 - e3,2)*Power(e1 - e3,2)*log_min_e2 + (e0 - e2)*(-e1 + e2)*Power(e3,2)*((e0 - e3)*(-e1 + e3)*(-e2 + e3) + (3.*e0*e1*e2 - 2.*(e0*e1 + (e0 + e1)*e2)*e3 + (e0 + e1 + e2)*Power(e3,2))*log_min_e3)))/((e1 - e2)*Power(e1 - e3,2)*Power(e2 - e3,2)))/(6.*(e0 - e1)*(e0 - e2)*Power(e0 - e3,2));
        }
        break;

      case TWOFOLD_DEGENERACY:
        {
          assert(not_equal(e0, e1));
          assert(not_equal(e1, e2));
          assert(not_equal(e2, e0));
          assert(are_equal(e2, e3));

          /*
            r[0] = ((Power(e0,2)*(e1 - e2) - e0*Power(e2,2) + e1*Power(e2,2))/((e0 - e1)*Power(e0 - e2,2)*(e1 - e2)) - (Power(e0,2)*(-3.*e1*e2 + e0*(e1 + 2.*e2))*Log(-e0))/(Power(e0 - e1,2)*Power(e0 - e2,3)) + ((Power(e1,3)*Log(-e1))/Power(e0 - e1,2) + (Power(e2,2)*(-3.*e0*e1 + 2.*e0*e2 + e1*e2)*log_min_e2)/Power(e0 - e2,3))/Power(e1 - e2,2))/6.;
            r[1] = (Power(e1,2)/((-e0 + e1)*Power(e1 - e2,2)) + (Power(e0,3)*Log(-e0))/(Power(e0 - e1,2)*Power(e0 - e2,2)) - (Power(e1,2)*(e0*(e1 - 3.*e2) + 2.*e1*e2)*Log(-e1))/(Power(e0 - e1,2)*Power(e1 - e2,3)) + (Power(e2,2)*((e0 - e2)*(e1 - e2) + (3.*e0*e1 - (e0 + 2.*e1)*e2)*log_min_e2))/(Power(e0 - e2,2)*Power(-e1 + e2,3)))/6.;
            r[2] = -(2.*Power(e0,3)*Power(e1 - e2,3)*Log(-e0) - 2.*Power(e1,3)*Power(e0 - e2,3)*Log(-e1) + (e0 - e1)*e2*((e0 - e2)*(e1 - e2)*(5.*e0*e1 - 3.*(e0 + e1)*e2 + Power(e2,2)) + 2.*(3.*Power(e0,2)*Power(e1,2) - 3.*e0*e1*(e0 + e1)*e2 + (Power(e0,2) + e0*e1 + Power(e1,2))*Power(e2,2))*log_min_e2))/ (12.*(e0 - e1)*Power(e0 - e2,3)*Power(-e1 + e2,3));
            r[3] = -(2.*Power(e0,3)*Power(e1 - e2,3)*Log(-e0) - 2.*Power(e1,3)*Power(e0 - e2,3)*Log(-e1) + (e0 - e1)*e2*((e0 - e2)*(e1 - e2)*(5.*e0*e1 - 3.*(e0 + e1)*e2 + Power(e2,2)) + 2.*(3.*Power(e0,2)*Power(e1,2) - 3.*e0*e1*(e0 + e1)*e2 + (Power(e0,2) + e0*e1 + Power(e1,2))*Power(e2,2))*log_min_e2))/ (12.*(e0 - e1)*Power(e0 - e2,3)*Power(-e1 + e2,3));
          */

          r[0] = ((Power(e0,2)*(e1 - e2) - e0*Power(e2,2) + e1*Power(e2,2))/((e0 - e1)*Power(e0 - e2,2)*(e1 - e2)) - (Power(e0,2)*(-3.*e1*e2 + e0*(e1 + 2.*e2))*log_min_e0)/(Power(e0 - e1,2)*Power(e0 - e2,3)) + ((Power(e1,3)*log_min_e1)/Power(e0 - e1,2) + (Power(e2,2)*(-3.*e0*e1 + 2.*e0*e2 + e1*e2)*log_min_e2)/Power(e0 - e2,3))/Power(e1 - e2,2))/6.;
          r[1] = (Power(e1,2)/((-e0 + e1)*Power(e1 - e2,2)) + (Power(e0,3)*log_min_e0)/(Power(e0 - e1,2)*Power(e0 - e2,2)) - (Power(e1,2)*(e0*(e1 - 3.*e2) + 2.*e1*e2)*log_min_e1)/(Power(e0 - e1,2)*Power(e1 - e2,3)) + (Power(e2,2)*((e0 - e2)*(e1 - e2) + (3.*e0*e1 - (e0 + 2.*e1)*e2)*log_min_e2))/(Power(e0 - e2,2)*Power(-e1 + e2,3)))/6.;
          r[2] = -(2.*Power(e0,3)*Power(e1 - e2,3)*log_min_e0 - 2.*Power(e1,3)*Power(e0 - e2,3)*log_min_e1 + (e0 - e1)*e2*((e0 - e2)*(e1 - e2)*(5.*e0*e1 - 3.*(e0 + e1)*e2 + Power(e2,2)) + 2.*(3.*Power(e0,2)*Power(e1,2) - 3.*e0*e1*(e0 + e1)*e2 + (Power(e0,2) + e0*e1 + Power(e1,2))*Power(e2,2))*log_min_e2))/ (12.*(e0 - e1)*Power(e0 - e2,3)*Power(-e1 + e2,3));
          r[3] = -(2.*Power(e0,3)*Power(e1 - e2,3)*log_min_e0 - 2.*Power(e1,3)*Power(e0 - e2,3)*log_min_e1 + (e0 - e1)*e2*((e0 - e2)*(e1 - e2)*(5.*e0*e1 - 3.*(e0 + e1)*e2 + Power(e2,2)) + 2.*(3.*Power(e0,2)*Power(e1,2) - 3.*e0*e1*(e0 + e1)*e2 + (Power(e0,2) + e0*e1 + Power(e1,2))*Power(e2,2))*log_min_e2))/ (12.*(e0 - e1)*Power(e0 - e2,3)*Power(-e1 + e2,3));
        }
        break;

      case THREEFOLD_DEGENERACY_A:
        {
          assert(not_equal(e0, e1));
          assert(are_equal(e2, e1));
          assert(are_equal(e3, e1));

          /*
            r[0] = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-Log(-e0) + Log(-e1)))/(12.*Power(e0 - e1,4));
            r[1] = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(Log(-e0) - Log(-e1)))/(36.*Power(e0 - e1,4));
            r[2] = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(Log(-e0) - Log(-e1)))/(36.*Power(e0 - e1,4));
            r[3] = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(Log(-e0) - Log(-e1)))/(36.*Power(e0 - e1,4));
          */

          r[0] = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-log_min_e0 + log_min_e1))/(12.*Power(e0 - e1,4));
          r[1] = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(log_min_e0 - log_min_e1))/(36.*Power(e0 - e1,4));
          r[2] = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(log_min_e0 - log_min_e1))/(36.*Power(e0 - e1,4));
          r[3] = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(log_min_e0 - log_min_e1))/(36.*Power(e0 - e1,4));
        }
        break;

      case THREEFOLD_DEGENERACY_B:
        {
          assert(not_equal(e0, e1));
          assert(are_equal(e0, e2));
          assert(are_equal(e1, e3));

          /*
            r[0] = ((e0 - e1)*(Power(e0,2) - 5.*e0*e1 - 2.*Power(e1,2)) + 6.*e0*Power(e1,2)*(Log(-e0) - Log(-e1)))/(12.*Power(e0 - e1,4));
            r[1] = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-Log(-e0) + Log(-e1)))/(12.*Power(e0 - e1,4));
            r[2] = ((e0 - e1)*(Power(e0,2) - 5.*e0*e1 - 2.*Power(e1,2)) + 6.*e0*Power(e1,2)*(Log(-e0) - Log(-e1)))/(12.*Power(e0 - e1,4));
            r[3] = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-Log(-e0) + Log(-e1)))/(12.*Power(e0 - e1,4));
          */

          r[0] = ((e0 - e1)*(Power(e0,2) - 5.*e0*e1 - 2.*Power(e1,2)) + 6.*e0*Power(e1,2)*(log_min_e0 - log_min_e1))/(12.*Power(e0 - e1,4));
          r[1] = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-log_min_e0 + log_min_e1))/(12.*Power(e0 - e1,4));
          r[2] = ((e0 - e1)*(Power(e0,2) - 5.*e0*e1 - 2.*Power(e1,2)) + 6.*e0*Power(e1,2)*(log_min_e0 - log_min_e1))/(12.*Power(e0 - e1,4));
          r[3] = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-log_min_e0 + log_min_e1))/(12.*Power(e0 - e1,4));
        }
        break;

      case FOURFOLD_DEGENERACY:
        {
          assert(are_equal(e0, e1));
          assert(are_equal(e0, e2));
          assert(are_equal(e0, e3));

          r[0] = 1./(24.*e0);
          r[1] = 1./(24.*e0);
          r[2] = 1./(24.*e0);
          r[3] = 1./(24.*e0);
        }
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }

    std::complex<scalartype> result=0;

    {
      result += f[0]*r[0];
      result += f[1]*r[1];
      result += f[2]*r[2];
      result += f[3]*r[3];

      assert(result==result); // make sure there is no NAN
    }

    return result;
  }

  template<typename scalartype>
  std::complex<scalartype> tetrahedron_routines_inverse_matrix_function::integrate_matrix_element_3D(eigenvalue_degeneracy_t                          degeneracy,
                                                                                                     std::vector<matrix_element_struct<scalartype> >& vec)
  {
    assert(vec.size()==4);

    //     eigenvalue_degeneracy_t degeneracy = find_degeneracy_3D(vec);

    std::complex<scalartype> r[4];
    std::complex<scalartype> f[4];

    int i0 = vec[0].i;
    int i1 = vec[1].i;
    int i2 = vec[2].i;
    int i3 = vec[3].i;

    f[0] = vec[i0].f;
    f[1] = vec[i1].f;
    f[2] = vec[i2].f;
    f[3] = vec[i3].f;

    std::complex<scalartype> e0 = vec[0].e;
    std::complex<scalartype> e1 = vec[1].e;
    std::complex<scalartype> e2 = vec[2].e;
    std::complex<scalartype> e3 = vec[3].e;

    std::complex<scalartype> log_min_e0 = vec[0].log_min_e;
    std::complex<scalartype> log_min_e1 = vec[1].log_min_e;
    std::complex<scalartype> log_min_e2 = vec[2].log_min_e;
    std::complex<scalartype> log_min_e3 = vec[3].log_min_e;

    switch(degeneracy)
      {
      case NO_DEGENERACY:
        {
          assert(not_equal(e0, e1));
          assert(not_equal(e1, e2));
          assert(not_equal(e2, e3));
          assert(not_equal(e3, e0));
          assert(not_equal(e0, e2));
          assert(not_equal(e1, e3));

          /*
            r[0] = (-((Power(e0,2)*(3.*e1*e2*e3 + Power(e0,2)*(e1 + e2 + e3) - 2.*e0*(e1*e2 + (e1 + e2)*e3))*Log(-e0))/(Power(e0 - e1,2)*Power(e0 - e2,2)*Power(e0 - e3,2))) + (Power(e1,3)*Log(-e1))/(Power(e0 - e1,2)*(e1 - e2)*(e1 - e3)) + (Power(e2,3)*Log(-e2))/(Power(e0 - e2,2)*(-e1 + e2)*(e2 - e3)) + ((Power(e0,2)*(e0 - e3))/((e0 - e1)*(e0 - e2)) + (Power(e3,3)*Log(-e3))/((-e1 + e3)*(-e2 + e3)))/Power(e0 - e3,2))/6.;
            r[1] = (Power(e1,2)/((-e0 + e1)*(e1 - e2)*(e1 - e3)) + (Power(e0,3)*Log(-e0))/(Power(e0 - e1,2)*(e0 - e2)*(e0 - e3)) - (Power(e1,2)*(e0*(Power(e1,2) + 3.*e2*e3 - 2.*e1*(e2 + e3)) + e1*(-2.*e2*e3 + e1*(e2 + e3)))*Log(-e1))/(Power(e0 - e1,2)*Power(e1 - e2,2)*Power(e1 - e3,2)) + ((Power(e2,3)*Log(-e2))/(Power(e1 - e2,2)*(-e0 + e2)) + (Power(e3,3)*Log(-e3))/((e0 - e3)*Power(e1 - e3,2)))/(e2 - e3))/6.;
            r[2] = (Power(e0,3)*Log(-e0) + (-(Power(e1,3)*Power(e0 - e2,2)*(e0 - e3)*Power(e2 - e3,2)*Log(-e1)) + (e0 - e1)*(Power(e2,2)*(e0 - e3)*(-e1 + e3)*(e0*(-2.*e1*e2 + Power(e2,2) + 3.*e1*e3 - 2.*e2*e3) + e2*(e1*e2 - 2.*e1*e3 + e2*e3))*Log(-e2) + (-e0 + e2)*(-e1 + e2)*(Power(e2,2)*(e0 - e3)*(-e1 + e3)*(-e2 + e3) + (e0 - e2)*(e1 - e2)*Power(e3,3)*Log(-e3))))/(Power(e1 - e2,2)*(e1 - e3)*Power(e2 - e3,2)))/(6.*(e0 - e1)*Power(e0 - e2,2)*(e0 - e3));
            r[3] = (Power(e0,3)*Log(-e0) + (-(Power(e1,3)*(e0 - e2)*Power(e0 - e3,2)*Power(e2 - e3,2)*Log(-e1)) + (e0 - e1)*(Power(e2,3)*Power(e0 - e3,2)*Power(e1 - e3,2)*Log(-e2) + (e0 - e2)*(-e1 + e2)*Power(e3,2)*((e0 - e3)*(-e1 + e3)*(-e2 + e3) + (3.*e0*e1*e2 - 2.*(e0*e1 + (e0 + e1)*e2)*e3 + (e0 + e1 + e2)*Power(e3,2))*Log(-e3))))/((e1 - e2)*Power(e1 - e3,2)*Power(e2 - e3,2)))/(6.*(e0 - e1)*(e0 - e2)*Power(e0 - e3,2));
          */

          r[0] = (-((Power(e0,2)*(3.*e1*e2*e3 + Power(e0,2)*(e1 + e2 + e3) - 2.*e0*(e1*e2 + (e1 + e2)*e3))*log_min_e0)/(Power(e0 - e1,2)*Power(e0 - e2,2)*Power(e0 - e3,2))) + (Power(e1,3)*log_min_e1)/(Power(e0 - e1,2)*(e1 - e2)*(e1 - e3)) + (Power(e2,3)*log_min_e2)/(Power(e0 - e2,2)*(-e1 + e2)*(e2 - e3)) + ((Power(e0,2)*(e0 - e3))/((e0 - e1)*(e0 - e2)) + (Power(e3,3)*log_min_e3)/((-e1 + e3)*(-e2 + e3)))/Power(e0 - e3,2))/6.;
          r[1] = (Power(e1,2)/((-e0 + e1)*(e1 - e2)*(e1 - e3)) + (Power(e0,3)*log_min_e0)/(Power(e0 - e1,2)*(e0 - e2)*(e0 - e3)) - (Power(e1,2)*(e0*(Power(e1,2) + 3.*e2*e3 - 2.*e1*(e2 + e3)) + e1*(-2.*e2*e3 + e1*(e2 + e3)))*log_min_e1)/(Power(e0 - e1,2)*Power(e1 - e2,2)*Power(e1 - e3,2)) + ((Power(e2,3)*log_min_e2)/(Power(e1 - e2,2)*(-e0 + e2)) + (Power(e3,3)*log_min_e3)/((e0 - e3)*Power(e1 - e3,2)))/(e2 - e3))/6.;
          r[2] = (Power(e0,3)*log_min_e0 + (-(Power(e1,3)*Power(e0 - e2,2)*(e0 - e3)*Power(e2 - e3,2)*log_min_e1) + (e0 - e1)*(Power(e2,2)*(e0 - e3)*(-e1 + e3)*(e0*(-2.*e1*e2 + Power(e2,2) + 3.*e1*e3 - 2.*e2*e3) + e2*(e1*e2 - 2.*e1*e3 + e2*e3))*log_min_e2 + (-e0 + e2)*(-e1 + e2)*(Power(e2,2)*(e0 - e3)*(-e1 + e3)*(-e2 + e3) + (e0 - e2)*(e1 - e2)*Power(e3,3)*log_min_e3)))/(Power(e1 - e2,2)*(e1 - e3)*Power(e2 - e3,2)))/(6.*(e0 - e1)*Power(e0 - e2,2)*(e0 - e3));
          r[3] = (Power(e0,3)*log_min_e0 + (-(Power(e1,3)*(e0 - e2)*Power(e0 - e3,2)*Power(e2 - e3,2)*log_min_e1) + (e0 - e1)*(Power(e2,3)*Power(e0 - e3,2)*Power(e1 - e3,2)*log_min_e2 + (e0 - e2)*(-e1 + e2)*Power(e3,2)*((e0 - e3)*(-e1 + e3)*(-e2 + e3) + (3.*e0*e1*e2 - 2.*(e0*e1 + (e0 + e1)*e2)*e3 + (e0 + e1 + e2)*Power(e3,2))*log_min_e3)))/((e1 - e2)*Power(e1 - e3,2)*Power(e2 - e3,2)))/(6.*(e0 - e1)*(e0 - e2)*Power(e0 - e3,2));
        }
        break;

      case TWOFOLD_DEGENERACY:
        {
          assert(not_equal(e0, e1));
          assert(not_equal(e1, e2));
          assert(not_equal(e2, e0));
          assert(are_equal(e2, e3));

          /*
            r[0] = ((Power(e0,2)*(e1 - e2) - e0*Power(e2,2) + e1*Power(e2,2))/((e0 - e1)*Power(e0 - e2,2)*(e1 - e2)) - (Power(e0,2)*(-3.*e1*e2 + e0*(e1 + 2.*e2))*Log(-e0))/(Power(e0 - e1,2)*Power(e0 - e2,3)) + ((Power(e1,3)*Log(-e1))/Power(e0 - e1,2) + (Power(e2,2)*(-3.*e0*e1 + 2.*e0*e2 + e1*e2)*log_min_e2)/Power(e0 - e2,3))/Power(e1 - e2,2))/6.;
            r[1] = (Power(e1,2)/((-e0 + e1)*Power(e1 - e2,2)) + (Power(e0,3)*Log(-e0))/(Power(e0 - e1,2)*Power(e0 - e2,2)) - (Power(e1,2)*(e0*(e1 - 3.*e2) + 2.*e1*e2)*Log(-e1))/(Power(e0 - e1,2)*Power(e1 - e2,3)) + (Power(e2,2)*((e0 - e2)*(e1 - e2) + (3.*e0*e1 - (e0 + 2.*e1)*e2)*log_min_e2))/(Power(e0 - e2,2)*Power(-e1 + e2,3)))/6.;
            r[2] = -(2.*Power(e0,3)*Power(e1 - e2,3)*Log(-e0) - 2.*Power(e1,3)*Power(e0 - e2,3)*Log(-e1) + (e0 - e1)*e2*((e0 - e2)*(e1 - e2)*(5.*e0*e1 - 3.*(e0 + e1)*e2 + Power(e2,2)) + 2.*(3.*Power(e0,2)*Power(e1,2) - 3.*e0*e1*(e0 + e1)*e2 + (Power(e0,2) + e0*e1 + Power(e1,2))*Power(e2,2))*log_min_e2))/ (12.*(e0 - e1)*Power(e0 - e2,3)*Power(-e1 + e2,3));
            r[3] = -(2.*Power(e0,3)*Power(e1 - e2,3)*Log(-e0) - 2.*Power(e1,3)*Power(e0 - e2,3)*Log(-e1) + (e0 - e1)*e2*((e0 - e2)*(e1 - e2)*(5.*e0*e1 - 3.*(e0 + e1)*e2 + Power(e2,2)) + 2.*(3.*Power(e0,2)*Power(e1,2) - 3.*e0*e1*(e0 + e1)*e2 + (Power(e0,2) + e0*e1 + Power(e1,2))*Power(e2,2))*log_min_e2))/ (12.*(e0 - e1)*Power(e0 - e2,3)*Power(-e1 + e2,3));
          */

          r[0] = ((Power(e0,2)*(e1 - e2) - e0*Power(e2,2) + e1*Power(e2,2))/((e0 - e1)*Power(e0 - e2,2)*(e1 - e2)) - (Power(e0,2)*(-3.*e1*e2 + e0*(e1 + 2.*e2))*log_min_e0)/(Power(e0 - e1,2)*Power(e0 - e2,3)) + ((Power(e1,3)*log_min_e1)/Power(e0 - e1,2) + (Power(e2,2)*(-3.*e0*e1 + 2.*e0*e2 + e1*e2)*log_min_e2)/Power(e0 - e2,3))/Power(e1 - e2,2))/6.;
          r[1] = (Power(e1,2)/((-e0 + e1)*Power(e1 - e2,2)) + (Power(e0,3)*log_min_e0)/(Power(e0 - e1,2)*Power(e0 - e2,2)) - (Power(e1,2)*(e0*(e1 - 3.*e2) + 2.*e1*e2)*log_min_e1)/(Power(e0 - e1,2)*Power(e1 - e2,3)) + (Power(e2,2)*((e0 - e2)*(e1 - e2) + (3.*e0*e1 - (e0 + 2.*e1)*e2)*log_min_e2))/(Power(e0 - e2,2)*Power(-e1 + e2,3)))/6.;
          r[2] = -(2.*Power(e0,3)*Power(e1 - e2,3)*log_min_e0 - 2.*Power(e1,3)*Power(e0 - e2,3)*log_min_e1 + (e0 - e1)*e2*((e0 - e2)*(e1 - e2)*(5.*e0*e1 - 3.*(e0 + e1)*e2 + Power(e2,2)) + 2.*(3.*Power(e0,2)*Power(e1,2) - 3.*e0*e1*(e0 + e1)*e2 + (Power(e0,2) + e0*e1 + Power(e1,2))*Power(e2,2))*log_min_e2))/ (12.*(e0 - e1)*Power(e0 - e2,3)*Power(-e1 + e2,3));
          r[3] = -(2.*Power(e0,3)*Power(e1 - e2,3)*log_min_e0 - 2.*Power(e1,3)*Power(e0 - e2,3)*log_min_e1 + (e0 - e1)*e2*((e0 - e2)*(e1 - e2)*(5.*e0*e1 - 3.*(e0 + e1)*e2 + Power(e2,2)) + 2.*(3.*Power(e0,2)*Power(e1,2) - 3.*e0*e1*(e0 + e1)*e2 + (Power(e0,2) + e0*e1 + Power(e1,2))*Power(e2,2))*log_min_e2))/ (12.*(e0 - e1)*Power(e0 - e2,3)*Power(-e1 + e2,3));
        }
        break;

      case THREEFOLD_DEGENERACY_A:
        {
          assert(not_equal(e0, e1));
          assert(are_equal(e2, e1));
          assert(are_equal(e3, e1));

          /*
            r[0] = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-Log(-e0) + Log(-e1)))/(12.*Power(e0 - e1,4));
            r[1] = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(Log(-e0) - Log(-e1)))/(36.*Power(e0 - e1,4));
            r[2] = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(Log(-e0) - Log(-e1)))/(36.*Power(e0 - e1,4));
            r[3] = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(Log(-e0) - Log(-e1)))/(36.*Power(e0 - e1,4));
          */

          r[0] = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-log_min_e0 + log_min_e1))/(12.*Power(e0 - e1,4));
          r[1] = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(log_min_e0 - log_min_e1))/(36.*Power(e0 - e1,4));
          r[2] = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(log_min_e0 - log_min_e1))/(36.*Power(e0 - e1,4));
          r[3] = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(log_min_e0 - log_min_e1))/(36.*Power(e0 - e1,4));
        }
        break;

      case THREEFOLD_DEGENERACY_B:
        {
          assert(not_equal(e0, e1));
          assert(are_equal(e0, e2));
          assert(are_equal(e1, e3));

          /*
            r[0] = ((e0 - e1)*(Power(e0,2) - 5.*e0*e1 - 2.*Power(e1,2)) + 6.*e0*Power(e1,2)*(Log(-e0) - Log(-e1)))/(12.*Power(e0 - e1,4));
            r[1] = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-Log(-e0) + Log(-e1)))/(12.*Power(e0 - e1,4));
            r[2] = ((e0 - e1)*(Power(e0,2) - 5.*e0*e1 - 2.*Power(e1,2)) + 6.*e0*Power(e1,2)*(Log(-e0) - Log(-e1)))/(12.*Power(e0 - e1,4));
            r[3] = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-Log(-e0) + Log(-e1)))/(12.*Power(e0 - e1,4));
          */

          r[0] = ((e0 - e1)*(Power(e0,2) - 5.*e0*e1 - 2.*Power(e1,2)) + 6.*e0*Power(e1,2)*(log_min_e0 - log_min_e1))/(12.*Power(e0 - e1,4));
          r[1] = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-log_min_e0 + log_min_e1))/(12.*Power(e0 - e1,4));
          r[2] = ((e0 - e1)*(Power(e0,2) - 5.*e0*e1 - 2.*Power(e1,2)) + 6.*e0*Power(e1,2)*(log_min_e0 - log_min_e1))/(12.*Power(e0 - e1,4));
          r[3] = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-log_min_e0 + log_min_e1))/(12.*Power(e0 - e1,4));
        }
        break;

      case FOURFOLD_DEGENERACY:
        {
          assert(are_equal(e0, e1));
          assert(are_equal(e0, e2));
          assert(are_equal(e0, e3));

          r[0] = 1./(24.*e0);
          r[1] = 1./(24.*e0);
          r[2] = 1./(24.*e0);
          r[3] = 1./(24.*e0);
        }
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }

    std::complex<scalartype> result=0;

    {
      result += f[0]*r[0];
      result += f[1]*r[1];
      result += f[2]*r[2];
      result += f[3]*r[3];

      assert(result==result); // make sure there is no NAN
    }

    return result;
  }

  template<typename scalartype>
  void tetrahedron_routines_inverse_matrix_function::integrate_eigenvalues_3D(eigenvalue_degeneracy_t                          degeneracy,
                                                                              std::vector<matrix_element_struct<scalartype> >& vec)
  {
    assert(vec.size()==4);

    int i0 = vec[0].i;
    int i1 = vec[1].i;
    int i2 = vec[2].i;
    int i3 = vec[3].i;

    std::complex<scalartype> e0 = vec[0].e;
    std::complex<scalartype> e1 = vec[1].e;
    std::complex<scalartype> e2 = vec[2].e;
    std::complex<scalartype> e3 = vec[3].e;

    std::complex<scalartype> log_min_e0 = vec[0].log_min_e;
    std::complex<scalartype> log_min_e1 = vec[1].log_min_e;
    std::complex<scalartype> log_min_e2 = vec[2].log_min_e;
    std::complex<scalartype> log_min_e3 = vec[3].log_min_e;

    switch(degeneracy)
      {
      case NO_DEGENERACY:
        {
          assert(not_equal(e0, e1));
          assert(not_equal(e1, e2));
          assert(not_equal(e2, e3));
          assert(not_equal(e3, e0));
          assert(not_equal(e0, e2));
          assert(not_equal(e1, e3));

          vec[i0].r = (-((Power(e0,2)*(3.*e1*e2*e3 + Power(e0,2)*(e1 + e2 + e3) - 2.*e0*(e1*e2 + (e1 + e2)*e3))*log_min_e0)/(Power(e0 - e1,2)*Power(e0 - e2,2)*Power(e0 - e3,2))) + (Power(e1,3)*log_min_e1)/(Power(e0 - e1,2)*(e1 - e2)*(e1 - e3)) + (Power(e2,3)*log_min_e2)/(Power(e0 - e2,2)*(-e1 + e2)*(e2 - e3)) + ((Power(e0,2)*(e0 - e3))/((e0 - e1)*(e0 - e2)) + (Power(e3,3)*log_min_e3)/((-e1 + e3)*(-e2 + e3)))/Power(e0 - e3,2))/6.;
          vec[i1].r = (Power(e1,2)/((-e0 + e1)*(e1 - e2)*(e1 - e3)) + (Power(e0,3)*log_min_e0)/(Power(e0 - e1,2)*(e0 - e2)*(e0 - e3)) - (Power(e1,2)*(e0*(Power(e1,2) + 3.*e2*e3 - 2.*e1*(e2 + e3)) + e1*(-2.*e2*e3 + e1*(e2 + e3)))*log_min_e1)/(Power(e0 - e1,2)*Power(e1 - e2,2)*Power(e1 - e3,2)) + ((Power(e2,3)*log_min_e2)/(Power(e1 - e2,2)*(-e0 + e2)) + (Power(e3,3)*log_min_e3)/((e0 - e3)*Power(e1 - e3,2)))/(e2 - e3))/6.;
          vec[i2].r = (Power(e0,3)*log_min_e0 + (-(Power(e1,3)*Power(e0 - e2,2)*(e0 - e3)*Power(e2 - e3,2)*log_min_e1) + (e0 - e1)*(Power(e2,2)*(e0 - e3)*(-e1 + e3)*(e0*(-2.*e1*e2 + Power(e2,2) + 3.*e1*e3 - 2.*e2*e3) + e2*(e1*e2 - 2.*e1*e3 + e2*e3))*log_min_e2 + (-e0 + e2)*(-e1 + e2)*(Power(e2,2)*(e0 - e3)*(-e1 + e3)*(-e2 + e3) + (e0 - e2)*(e1 - e2)*Power(e3,3)*log_min_e3)))/(Power(e1 - e2,2)*(e1 - e3)*Power(e2 - e3,2)))/(6.*(e0 - e1)*Power(e0 - e2,2)*(e0 - e3));
          vec[i3].r = (Power(e0,3)*log_min_e0 + (-(Power(e1,3)*(e0 - e2)*Power(e0 - e3,2)*Power(e2 - e3,2)*log_min_e1) + (e0 - e1)*(Power(e2,3)*Power(e0 - e3,2)*Power(e1 - e3,2)*log_min_e2 + (e0 - e2)*(-e1 + e2)*Power(e3,2)*((e0 - e3)*(-e1 + e3)*(-e2 + e3) + (3.*e0*e1*e2 - 2.*(e0*e1 + (e0 + e1)*e2)*e3 + (e0 + e1 + e2)*Power(e3,2))*log_min_e3)))/((e1 - e2)*Power(e1 - e3,2)*Power(e2 - e3,2)))/(6.*(e0 - e1)*(e0 - e2)*Power(e0 - e3,2));
        }
        break;

      case TWOFOLD_DEGENERACY:
        {
          assert(not_equal(e0, e1));
          assert(not_equal(e1, e2));
          assert(not_equal(e2, e0));
          assert(are_equal(e2, e3));

          vec[i0].r = ((Power(e0,2)*(e1 - e2) - e0*Power(e2,2) + e1*Power(e2,2))/((e0 - e1)*Power(e0 - e2,2)*(e1 - e2)) - (Power(e0,2)*(-3.*e1*e2 + e0*(e1 + 2.*e2))*log_min_e0)/(Power(e0 - e1,2)*Power(e0 - e2,3)) + ((Power(e1,3)*log_min_e1)/Power(e0 - e1,2) + (Power(e2,2)*(-3.*e0*e1 + 2.*e0*e2 + e1*e2)*log_min_e2)/Power(e0 - e2,3))/Power(e1 - e2,2))/6.;
          vec[i1].r = (Power(e1,2)/((-e0 + e1)*Power(e1 - e2,2)) + (Power(e0,3)*log_min_e0)/(Power(e0 - e1,2)*Power(e0 - e2,2)) - (Power(e1,2)*(e0*(e1 - 3.*e2) + 2.*e1*e2)*log_min_e1)/(Power(e0 - e1,2)*Power(e1 - e2,3)) + (Power(e2,2)*((e0 - e2)*(e1 - e2) + (3.*e0*e1 - (e0 + 2.*e1)*e2)*log_min_e2))/(Power(e0 - e2,2)*Power(-e1 + e2,3)))/6.;
          vec[i2].r = -(2.*Power(e0,3)*Power(e1 - e2,3)*log_min_e0 - 2.*Power(e1,3)*Power(e0 - e2,3)*log_min_e1 + (e0 - e1)*e2*((e0 - e2)*(e1 - e2)*(5.*e0*e1 - 3.*(e0 + e1)*e2 + Power(e2,2)) + 2.*(3.*Power(e0,2)*Power(e1,2) - 3.*e0*e1*(e0 + e1)*e2 + (Power(e0,2) + e0*e1 + Power(e1,2))*Power(e2,2))*log_min_e2))/ (12.*(e0 - e1)*Power(e0 - e2,3)*Power(-e1 + e2,3));
          vec[i3].r = -(2.*Power(e0,3)*Power(e1 - e2,3)*log_min_e0 - 2.*Power(e1,3)*Power(e0 - e2,3)*log_min_e1 + (e0 - e1)*e2*((e0 - e2)*(e1 - e2)*(5.*e0*e1 - 3.*(e0 + e1)*e2 + Power(e2,2)) + 2.*(3.*Power(e0,2)*Power(e1,2) - 3.*e0*e1*(e0 + e1)*e2 + (Power(e0,2) + e0*e1 + Power(e1,2))*Power(e2,2))*log_min_e2))/ (12.*(e0 - e1)*Power(e0 - e2,3)*Power(-e1 + e2,3));
        }
        break;

      case THREEFOLD_DEGENERACY_A:
        {
          assert(not_equal(e0, e1));
          assert(are_equal(e2, e1));
          assert(are_equal(e3, e1));

          vec[i0].r = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-log_min_e0 + log_min_e1))/(12.*Power(e0 - e1,4));
          vec[i1].r = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(log_min_e0 - log_min_e1))/(36.*Power(e0 - e1,4));
          vec[i2].r = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(log_min_e0 - log_min_e1))/(36.*Power(e0 - e1,4));
          vec[i3].r = (-((e0 - e1)*(11.*Power(e0,2) - 7.*e0*e1 + 2.*Power(e1,2))) + 6.*Power(e0,3)*(log_min_e0 - log_min_e1))/(36.*Power(e0 - e1,4));
        }
        break;

      case THREEFOLD_DEGENERACY_B:
        {
          assert(not_equal(e0, e1));
          assert(are_equal(e0, e2));
          assert(are_equal(e1, e3));

          vec[i0].r = ((e0 - e1)*(Power(e0,2) - 5.*e0*e1 - 2.*Power(e1,2)) + 6.*e0*Power(e1,2)*(log_min_e0 - log_min_e1))/(12.*Power(e0 - e1,4));
          vec[i1].r = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-log_min_e0 + log_min_e1))/(12.*Power(e0 - e1,4));
          vec[i2].r = ((e0 - e1)*(Power(e0,2) - 5.*e0*e1 - 2.*Power(e1,2)) + 6.*e0*Power(e1,2)*(log_min_e0 - log_min_e1))/(12.*Power(e0 - e1,4));
          vec[i3].r = (2.*Power(e0,3) + 3.*Power(e0,2)*e1 - 6.*e0*Power(e1,2) + Power(e1,3) + 6.*Power(e0,2)*e1*(-log_min_e0 + log_min_e1))/(12.*Power(e0 - e1,4));
        }
        break;

      case FOURFOLD_DEGENERACY:
        {
          assert(are_equal(e0, e1));
          assert(are_equal(e0, e2));
          assert(are_equal(e0, e3));

          vec[i0].r = 1./(24.*e0);
          vec[i1].r = 1./(24.*e0);
          vec[i2].r = 1./(24.*e0);
          vec[i3].r = 1./(24.*e0);
        }
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }
  }

  template<typename scalartype>
  eigenvalue_degeneracy_t tetrahedron_routines_inverse_matrix_function::find_degeneracy_3D(std::vector<matrix_element_struct<scalartype> >& vec)
  {
    assert(vec.size()==4);

    if(not_equal(vec[1].e, vec[0].e) and
       not_equal(vec[2].e, vec[1].e) and
       not_equal(vec[3].e, vec[2].e) and
       not_equal(vec[0].e, vec[3].e) and
       not_equal(vec[3].e, vec[1].e) and
       not_equal(vec[2].e, vec[0].e))
      {
        return NO_DEGENERACY;
      }

    if(are_equal(vec[1].e, vec[0].e) and
       are_equal(vec[2].e, vec[0].e) and
       are_equal(vec[3].e, vec[0].e))
      {
        return FOURFOLD_DEGENERACY;
      }

    vec[0].i = 0;
    vec[1].i = 1;
    vec[2].i = 2;
    vec[3].i = 3;

    do
      {
        if(not_equal(vec[1].e, vec[0].e) and
           not_equal(vec[2].e, vec[1].e) and
           not_equal(vec[0].e, vec[2].e) and
           are_equal(vec[3].e, vec[2].e))
          {
            return TWOFOLD_DEGENERACY;
          }

        if(not_equal(vec[1].e, vec[0].e) and
           are_equal(vec[2].e, vec[1].e) and
           are_equal(vec[3].e, vec[1].e))
          {
            return THREEFOLD_DEGENERACY_A;
          }

        if(not_equal(vec[1].e, vec[0].e) and
           are_equal(vec[2].e, vec[0].e) and
           are_equal(vec[3].e, vec[1].e))
          {
            return THREEFOLD_DEGENERACY_B;
          }
      }
    while
      ( std::next_permutation(vec.begin(), vec.end(), tetrahedron_routines_inverse_matrix_function::permutation_comp<scalartype>) );

    throw std::logic_error(__FUNCTION__);
    return NO_DEGENERACY;
  }

}

#endif
