//-*-C++-*-

#ifndef DCA_TETRAHEDRON_INTEGRATION_ROUTINES_H
#define DCA_TETRAHEDRON_INTEGRATION_ROUTINES_H

namespace DCA
{

  class tetrahedron_integration_routines
  {
  private:

    template<typename scalartype>
    static bool pair_same(std::pair<complex<scalartype>, complex<scalartype> > const& x,
                          std::pair<complex<scalartype>, complex<scalartype> > const& y);

  public:

    // 1D
    template<typename scalartype>
    static void execute(int size, scalartype volume,
                        std::complex<scalartype>* G_0,
                        std::complex<scalartype>* G_1,
                        std::complex<scalartype>* f_result);

    // 2D
    template<typename scalartype>
    static void execute(int size, scalartype volume,
                        std::complex<scalartype>* G_0,
                        std::complex<scalartype>* G_1,
                        std::complex<scalartype>* G_2,
                        std::complex<scalartype>* f_result);

    // 3D
    template<typename scalartype>
    static void execute(int size, scalartype volume,
                        std::complex<scalartype>* G_0,
                        std::complex<scalartype>* G_1,
                        std::complex<scalartype>* G_2,
                        std::complex<scalartype>* G_3,
                        std::complex<scalartype>* f_result);

  private:

    template<typename scalartype>
    static std::complex<scalartype> integrate_matrix_element_1D(std::complex<scalartype>* f,
                                                                std::complex<scalartype>* e);

    template<typename scalartype>
    static std::complex<scalartype> integrate_matrix_element_2D(std::complex<scalartype>* f,
                                                                std::complex<scalartype>* e);

    template<typename scalartype>
    static std::complex<scalartype> integrate_matrix_element_3D(std::complex<scalartype>* f,
                                                                std::complex<scalartype>* e);


    template<typename scalartype>
    static eigenvalue_degeneracy_t find_degeneracy_1D(std::complex<scalartype>* f,
                                                      std::complex<scalartype>* e);

    template<typename scalartype>
    static eigenvalue_degeneracy_t find_degeneracy_2D(std::complex<scalartype>* f,
                                                      std::complex<scalartype>* e);

    template<typename scalartype>
    static eigenvalue_degeneracy_t find_degeneracy_3D(std::complex<scalartype>* f,
                                                      std::complex<scalartype>* e);

  };

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

  /************************************
   ***
   ***   1D tetrahedron-integration
   ***
   ************************************/

  template<typename scalartype>
  void tetrahedron_integration_routines::execute(int N, scalartype volume,
                                                 std::complex<scalartype>* G_0,
                                                 std::complex<scalartype>* G_1,
                                                 std::complex<scalartype>* f_result)
  {
    throw std::logic_error(__FUNCTION__);
  }

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
    tetrahedron_integration_data<scalartype> data_obj(N);

    {// diagonolize the G-matrices
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
    }

    {// integrate G-mattices
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

  /************************************
   ***
   ***   3D tetrahedron-integration
   ***
   ************************************/

  template<typename scalartype>
  void tetrahedron_integration_routines::execute(int N, scalartype volume,
                                                 std::complex<scalartype>* G_0,
                                                 std::complex<scalartype>* G_1,
                                                 std::complex<scalartype>* G_2,
                                                 std::complex<scalartype>* G_3,
                                                 std::complex<scalartype>* f_result)
  {
    tetrahedron_integration_data<scalartype> data_obj(N);

    {// diagonolize the G-matrices
      {
        memcpy(data_obj.G_inv_0, G_0, sizeof(std::complex<scalartype>)*N*N);
        memcpy(data_obj.G_inv_1, G_1, sizeof(std::complex<scalartype>)*N*N);
        memcpy(data_obj.G_inv_2, G_2, sizeof(std::complex<scalartype>)*N*N);
        memcpy(data_obj.G_inv_3, G_3, sizeof(std::complex<scalartype>)*N*N);

        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_0); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_0, data_obj.G_inv_0));
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_1); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_1, data_obj.G_inv_1));
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_2); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_2, data_obj.G_inv_2));
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.G_inv_3); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, G_3, data_obj.G_inv_3));
      }

      {
        int INFO  = -1;
        int LWORK = 16*std::max(1,2*N-1);

        scalartype                RWORK[std::max(1, 3*N-2)];
        std::complex<scalartype>*  WORK = new std::complex<scalartype>[LWORK];

        LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'V', N, data_obj.G_inv_0, N, data_obj.W_0, data_obj.VR_inv_0, N, data_obj.VR_0, N, WORK, LWORK, RWORK, INFO);
        LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'V', N, data_obj.G_inv_1, N, data_obj.W_1, data_obj.VR_inv_1, N, data_obj.VR_1, N, WORK, LWORK, RWORK, INFO);
        LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'V', N, data_obj.G_inv_2, N, data_obj.W_2, data_obj.VR_inv_2, N, data_obj.VR_2, N, WORK, LWORK, RWORK, INFO);
        LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'V', N, data_obj.G_inv_3, N, data_obj.W_3, data_obj.VR_inv_3, N, data_obj.VR_3, N, WORK, LWORK, RWORK, INFO);

        delete [] WORK;
      }

      {
        memcpy(data_obj.VR_inv_0, data_obj.VR_0, sizeof(std::complex<scalartype>)*N*N);
        memcpy(data_obj.VR_inv_1, data_obj.VR_1, sizeof(std::complex<scalartype>)*N*N);
        memcpy(data_obj.VR_inv_2, data_obj.VR_2, sizeof(std::complex<scalartype>)*N*N);
        memcpy(data_obj.VR_inv_3, data_obj.VR_3, sizeof(std::complex<scalartype>)*N*N);

        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_0); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_0, data_obj.VR_inv_0));
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_1); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_1, data_obj.VR_inv_1));
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_2); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_2, data_obj.VR_inv_2));
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(N, data_obj.VR_inv_3); assert(LIN_ALG::GEINV<LIN_ALG::CPU>::test(N, data_obj.VR_3, data_obj.VR_inv_3));
      }
    }

    {// integrate G-matrices
      for(int j=0; j<N; j++)
        for(int i=0; i<N; i++)
          f_result[i + j*N] = 0.;

      std::complex<scalartype> matrix_elements[4];
      std::complex<scalartype> eigenvalues[4];

      for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){

          for(int l=0; l<N; l++){

            matrix_elements[0] = data_obj.VR_inv_0[i+l*N]*data_obj.VR_0[l+j*N]; //cout << matrix_elements[0] << endl;
            matrix_elements[1] = data_obj.VR_inv_1[i+l*N]*data_obj.VR_1[l+j*N]; //cout << matrix_elements[1] << endl;
            matrix_elements[2] = data_obj.VR_inv_2[i+l*N]*data_obj.VR_2[l+j*N]; //cout << matrix_elements[2] << endl;
            matrix_elements[3] = data_obj.VR_inv_3[i+l*N]*data_obj.VR_3[l+j*N]; //cout << matrix_elements[3] << endl;

            eigenvalues[0] = data_obj.W_0[l]; //cout << eigenvalues[0] << endl;
            eigenvalues[1] = data_obj.W_1[l]; //cout << eigenvalues[1] << endl;
            eigenvalues[2] = data_obj.W_2[l]; //cout << eigenvalues[2] << endl;
            eigenvalues[3] = data_obj.W_3[l]; //cout << eigenvalues[3] << endl;

            f_result[i+j*N] += integrate_matrix_element_3D(matrix_elements, eigenvalues);
          }
        }
      }

      //   double volume_tetraheder = compute_volume(&mesh[index[0]][0], &mesh[index[1]][0], &mesh[index[2]][0], &mesh[index[3]][0]);
      //   assert(volume_tetraheder > 1.e-6);

      for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
          f_result[i+j*N] *= (6.*volume);//(6.*volume_tetraheder);
    }
  }

  template<typename scalartype>
  std::complex<scalartype> tetrahedron_integration_routines::integrate_matrix_element_3D(std::complex<scalartype>* f,
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
        r[2] = -(VALUE[2]*std::pow(e[0],3)*std::pow(e[1] - e[2],3)*std::log(-e[0]) - VALUE[2]*std::pow(e[1],3)*std::pow(e[0] - e[2],3)*std::log(-e[1]) + (e[0] - e[1])*e[2]*((e[0] - e[2])*(e[1] - e[2])*(5.*e[0]*e[1] - VALUE[3]*e[0]*e[2] - VALUE[3]*e[1]*e[2] + std::pow(e[2],2)) + VALUE[2]*(std::pow(e[1],2)*std::pow(e[2],2) + e[0]*e[1]*e[2]*(-VALUE[3]*e[1] + e[2]) + std::pow(e[0],2)*(VALUE[3]*std::pow(e[1],2) - VALUE[3]*e[1]*e[2] + std::pow(e[2],2)))*std::log(-e[2])))/(VALUE[12]*(-e[0] + e[1])*std::pow(e[0] - e[2],3)*std::pow(e[1] - e[2],3));
        r[3] = r[2];//-(VALUE[2]*std::pow(e[0],3)*std::pow(e[1] - e[2],3)*std::log(-e[0]) - VALUE[2]*std::pow(e[1],3)*std::pow(e[0] - e[2],3)*std::log(-e[1]) + (e[0] - e[1])*e[2]*((e[0] - e[2])*(e[1] - e[2])*(5.*e[0]*e[1] - VALUE[3]*e[0]*e[2] - VALUE[3]*e[1]*e[2] + std::pow(e[2],2)) + VALUE[2]*(std::pow(e[1],2)*std::pow(e[2],2) + e[0]*e[1]*e[2]*(-VALUE[3]*e[1] + e[2]) + std::pow(e[0],2)*(VALUE[3]*std::pow(e[1],2) - VALUE[3]*e[1]*e[2] + std::pow(e[2],2)))*std::log(-e[2])))/(VALUE[12]*(-e[0] + e[1])*std::pow(e[0] - e[2],3)*std::pow(e[1] - e[2],3));
        break;

      case THREEFOLD_DEGENERACY:
        r[0] = (VALUE[2]*std::pow(e[0],3) + VALUE[3]*std::pow(e[0],2)*e[1] - VALUE[6]*e[0]*std::pow(e[1],2) + std::pow(e[1],3) - VALUE[6]*std::pow(e[0],2)*e[1]*std::log(-e[0]) + VALUE[6]*std::pow(e[0],2)*e[1]*std::log(-e[1]))/(VALUE[12]*std::pow(e[0] - e[1],4));
        r[1] = (-VALUE[11]*std::pow(e[0],3) + VALUE[18]*std::pow(e[0],2)*e[1] - VALUE[9]*e[0]*std::pow(e[1],2) + VALUE[2]*std::pow(e[1],3) + VALUE[6]*std::pow(e[0],3)*std::log(-e[0]) - VALUE[6]*std::pow(e[0],3)*std::log(-e[1]))/(VALUE[36]*std::pow(e[0] - e[1],4));
        r[2] = r[1];//(-VALUE[11]*std::pow(e[0],3) + VALUE[18]*std::pow(e[0],2)*e[1] - VALUE[9]*e[0]*std::pow(e[1],2) + VALUE[2]*std::pow(e[1],3) + 6.*std::pow(e[0],3)*std::log(-e[0]) - 6.*std::pow(e[0],3)*std::log(-e[1]))/(VALUE[36]*std::pow(e[0] - e[1],4));
        r[3] = r[1];//(-VALUE[11]*std::pow(e[0],3) + VALUE[18]*std::pow(e[0],2)*e[1] - VALUE[9]*e[0]*std::pow(e[1],2) + VALUE[2]*std::pow(e[1],3) + 6.*std::pow(e[0],3)*std::log(-e[0]) - 6.*std::pow(e[0],3)*std::log(-e[1]))/(VALUE[36]*std::pow(e[0] - e[1],4));
        break;

      case FOURFOLD_DEGENERACY:
        r[0] = VALUE[1]/(VALUE[24]*e[0]);
        r[1] = r[0];//1./24.*1./e[0];
        r[2] = r[0];//1./24.*1./e[0];
        r[3] = r[0];//1./24.*1./e[0];
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }

    return f[0]*r[0]+f[1]*r[1]+f[2]*r[2]+f[3]*r[3];
  }

  template<typename scalartype>
  eigenvalue_degeneracy_t tetrahedron_integration_routines::find_degeneracy_3D(std::complex<scalartype>* f,
                                                                               std::complex<scalartype>* e)
  {
    std::vector<std::pair<std::complex<scalartype>, std::complex<scalartype> > > vec(3+1);

    vec[0].first = e[0]; vec[0].second = f[0];
    vec[1].first = e[1]; vec[1].second = f[1];
    vec[2].first = e[2]; vec[2].second = f[2];
    vec[3].first = e[3]; vec[3].second = f[3];

    stable_sort(vec.begin(), vec.end(), &pair_less);
    int degeneracy = unique(vec.begin(), vec.end(), tetrahedron_integration_routines::pair_same<scalartype>)-vec.begin();

    e[0] = vec[0].first; f[0] = vec[0].second;
    e[1] = vec[1].first; f[1] = vec[1].second;
    e[2] = vec[2].first; f[2] = vec[2].second;
    e[3] = vec[3].first; f[3] = vec[3].second;

    if(degeneracy == 4)
      return NO_DEGENERACY;

    if(degeneracy == 3)
      return TWOFOLD_DEGENERACY;

    if(degeneracy == 2)
      return THREEFOLD_DEGENERACY;

    if(degeneracy == 1)
      return FOURFOLD_DEGENERACY;

    throw std::logic_error(__FUNCTION__);
    return NO_DEGENERACY;
  }

}

#endif
