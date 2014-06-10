//-*-C++-*-

#ifndef DCA_TETRAHEDRON_INTEGRATION_ROUTINES_H
#define DCA_TETRAHEDRON_INTEGRATION_ROUTINES_H

namespace DCA
{

  class tetrahedron_integration_routines
  {
  public:

    // 1D
    template<typename scalartype>
    static void execute(int size, scalartype volume, scalartype* H_0, scalartype* H_1, scalartype* f_result);

    // 2D
    template<typename scalartype>
    static void execute(int size, scalartype volume, std::complex<scalartype>* H_0, std::complex<scalartype>* H_1, std::complex<scalartype>* H_2, std::complex<scalartype>* f_result);

    // 3D
    template<typename scalartype>
    static void execute(int size, scalartype volume, scalartype* H_0, scalartype* H_1, scalartype* H_2, scalartype* H_3, scalartype* f_result);

  private:

    template<typename scalartype>
    static std::complex<scalartype> integrate_matrix_element_2D(std::complex<scalartype>* f, 
								std::complex<scalartype>* e);

    template<typename scalartype>
    static eigenvalue_degeneracy_t find_degeneracy_2D(std::complex<scalartype>* f, 
						      std::complex<scalartype>* e);

    template<typename scalartype>
    static bool pair_same(std::pair<complex<scalartype>, complex<scalartype> > const& x,
                          std::pair<complex<scalartype>, complex<scalartype> > const& y);

  };

  /************************************
   ***
   ***   2D tetrahedron-integration
   ***
   ************************************/

  template<typename scalartype>
  void tetrahedron_integration_routines::execute(int N, scalartype volume,
                                                 std::complex<scalartype>* G_0,
                                                 std::complex<scalartype>* G_1,
                                                 std::complex<scalartype>* G_2,
                                                 std::complex<scalartype>* f_result)
  {
//     std::complex<scalartype> G_inv_0[N*N]; std::complex<scalartype> VR_0[N*N]; std::complex<scalartype> VR_inv_0[N*N]; std::complex<scalartype> W_0[N];
//     std::complex<scalartype> G_inv_1[N*N]; std::complex<scalartype> VR_1[N*N]; std::complex<scalartype> VR_inv_1[N*N]; std::complex<scalartype> W_1[N];
//     std::complex<scalartype> G_inv_2[N*N]; std::complex<scalartype> VR_2[N*N]; std::complex<scalartype> VR_inv_2[N*N]; std::complex<scalartype> W_2[N];

    tetrahedron_integration_data<scalartype> data_obj(N);

    {
      memcpy(data_obj.G_inv_0, G_0, sizeof(std::complex<scalartype>)*N*N);
      memcpy(data_obj.G_inv_1, G_1, sizeof(std::complex<scalartype>)*N*N);
      memcpy(data_obj.G_inv_2, G_2, sizeof(std::complex<scalartype>)*N*N);

      LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_0); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_0, data_obj.G_inv_0));
      LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_1); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_1, data_obj.G_inv_1));
      LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_2); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_2, data_obj.G_inv_2));
    }

    {
      int INFO  = -1;
      int LWORK = 16*std::max(1,2*N-1);

      scalartype                RWORK[std::max(1, 3*N-2)];
      std::complex<scalartype>*  WORK = new std::complex<scalartype>[LWORK];

      LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'V', N, data_obj.G_inv_0, N, data_obj.W_0, data_obj.VR_inv_0, N, data_obj.VR_0, N, WORK, LWORK, RWORK, INFO);
      LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'V', N, data_obj.G_inv_1, N, data_obj.W_1, data_obj.VR_inv_1, N, data_obj.VR_1, N, WORK, LWORK, RWORK, INFO);
      LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'V', N, data_obj.G_inv_2, N, data_obj.W_2, data_obj.VR_inv_2, N, data_obj.VR_2, N, WORK, LWORK, RWORK, INFO);

      delete [] WORK;
    }

    {
      memcpy(data_obj.VR_inv_0, data_obj.VR_0, sizeof(std::complex<scalartype>)*N*N);
      memcpy(data_obj.VR_inv_1, data_obj.VR_1, sizeof(std::complex<scalartype>)*N*N);
      memcpy(data_obj.VR_inv_2, data_obj.VR_2, sizeof(std::complex<scalartype>)*N*N);

      LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_0); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_0, data_obj.VR_inv_0));
      LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_1); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_1, data_obj.VR_inv_1));
      LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_2); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_2, data_obj.VR_inv_2));
    }

    std::complex<scalartype> eigenvalues    [3];
    std::complex<scalartype> matrix_elements[3];

    for(int l=0; l<N*N; l++)
      f_result[l] = 0;

    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){

        for(int l=0; l<N; l++){

          matrix_elements[0] = data_obj.VR_inv_0[i+l*N]*data_obj.VR_0[l+j*N];
          matrix_elements[1] = data_obj.VR_inv_1[i+l*N]*data_obj.VR_1[l+j*N];
          matrix_elements[2] = data_obj.VR_inv_2[i+l*N]*data_obj.VR_2[l+j*N];

          eigenvalues[0] = data_obj.W_0[l];
          eigenvalues[1] = data_obj.W_1[l];
          eigenvalues[2] = data_obj.W_2[l];

          f_result[i + j*N] += integrate_matrix_element_2D(matrix_elements, eigenvalues);
        }
      }
    }

    for(int i=0; i<N; i++)
      for(int j=0; j<N; j++)
        f_result[i + j*N] *= (2.*volume);
  }

  template<typename scalartype>
  std::complex<scalartype> tetrahedron_integration_routines::integrate_matrix_element_2D(std::complex<scalartype>* f, 
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
  eigenvalue_degeneracy_t tetrahedron_integration_routines::find_degeneracy_2D(std::complex<scalartype>* f, 
									       std::complex<scalartype>* e)
  {
    std::vector<std::pair<std::complex<scalartype>, std::complex<scalartype> > > vec(3);

    vec[0].first = e[0]; vec[0].second = f[0];
    vec[1].first = e[1]; vec[1].second = f[1];
    vec[2].first = e[2]; vec[2].second = f[2];

    std::stable_sort(vec.begin(), vec.end(), &pair_less);
    int degeneracy = std::unique(vec.begin(), vec.end(), tetrahedron_integration_routines::pair_same<scalartype>)-vec.begin();

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

  template<typename scalartype>
  bool tetrahedron_integration_routines::pair_same(std::pair<complex<scalartype>, complex<scalartype> > const& x,
                                                   std::pair<complex<scalartype>, complex<scalartype> > const& y)
  {
    scalartype abs_x = abs(x.first);
    scalartype abs_y = abs(y.first);

    if(abs_x < 1. && abs_y < 1.)
      {
        return abs(x.first-y.first)<1.e-6;
      }
    else
      {
        scalartype MAX = abs_x>abs_y? abs_x:abs_y;
        return abs(x.first-y.first)<((1.e-6)*MAX);
      }
  }

}

#endif
