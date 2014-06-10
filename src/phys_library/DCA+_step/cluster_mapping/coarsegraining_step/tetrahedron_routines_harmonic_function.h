//-*-C++-*-

#ifndef DCA_TETRAHEDRON_ROUTINES_HARMONIC_FUNCTION_H
#define DCA_TETRAHEDRON_ROUTINES_HARMONIC_FUNCTION_H

namespace DCA
{

  class tetrahedron_routines_harmonic_function
  {
  public:

    // 1D
    static std::complex<double> execute(std::vector<double>&             r_vec,
                                        MATH_ALGORITHMS::tetrahedron<1>& tetrahedron);
    // 2D
    static std::complex<double> execute(std::vector<double>&             r_vec,
                                        MATH_ALGORITHMS::tetrahedron<2>& tetrahedron);

    // 3D
    static std::complex<double> execute(std::vector<double>&             r_vec,
                                        MATH_ALGORITHMS::tetrahedron<3>& tetrahedron);

  private:

    // general functionality:
    template<typename scalartype>
    static scalartype Power(scalartype x, int n);

    template<typename scalartype>
    static scalartype Sin(scalartype x);

    template<typename scalartype>
    static scalartype Cos(scalartype x);


    // 2D cases :
    static void permute(MATH_ALGORITHMS::tetrahedron<2>& tetrahedron_new,
                        MATH_ALGORITHMS::tetrahedron<2>& tetrahedron_old);

    static std::complex<double> case_2D   (double dotRD1, double dotRD2, double dotRD2minD1);
    static std::complex<double> case_d1_2D(double dotRD1, double dotRD2, double dotRD2minD1);
    static std::complex<double> case_d2_2D(double dotRD1, double dotRD2, double dotRD2minD1);

    // 3D
    static void permute(MATH_ALGORITHMS::tetrahedron<3>& tetrahedron_new,
                        MATH_ALGORITHMS::tetrahedron<3>& tetrahedron_old);

    static std::complex<double> case_3D(double dotRD1     , double dotRD2     , double dotRD3,
                                        double dotRD2minD1, double dotRD3minD2, double dotRD1minD3);

    static std::complex<double> case_d1_3D(double dotRD1     , double dotRD2     , double dotRD3,
                                           double dotRD2minD1, double dotRD3minD2, double dotRD1minD3);


    static std::complex<double> case_d2_3D(double dotRD1     , double dotRD2     , double dotRD3,
                                           double dotRD2minD1, double dotRD3minD2, double dotRD1minD3);


    static std::complex<double> case_d3_3D(double dotRD1     , double dotRD2     , double dotRD3,
                                           double dotRD2minD1, double dotRD3minD2, double dotRD1minD3);

    static std::complex<double> case_d1_d2_3D(double dotRD1     , double dotRD2     , double dotRD3,
                                              double dotRD2minD1, double dotRD3minD2, double dotRD1minD3);


    static std::complex<double> case_d2_d3_3D(double dotRD1     , double dotRD2     , double dotRD3,
                                              double dotRD2minD1, double dotRD3minD2, double dotRD1minD3);


    static std::complex<double> case_d3_d1_3D(double dotRD1     , double dotRD2     , double dotRD3,
                                              double dotRD2minD1, double dotRD3minD2, double dotRD1minD3);

  };

  template<typename scalartype>
  scalartype tetrahedron_routines_harmonic_function::Power(scalartype x, int n)
  {
    return std::pow(x,n);
  }

  template<typename scalartype>
  scalartype tetrahedron_routines_harmonic_function::Sin(scalartype x)
  {
    return std::sin(x);
  }

  template<typename scalartype>
  scalartype tetrahedron_routines_harmonic_function::Cos(scalartype x)
  {
    return std::cos(x);
  }

  /************************************
   ***
   ***   1D harmonic-integration
   ***
   ************************************/

  std::complex<double> tetrahedron_routines_harmonic_function::execute(std::vector<double>&             r_vec,
                                                                       MATH_ALGORITHMS::tetrahedron<1>& tetrahedron)
  {
    assert(r_vec.size()==1);

    double r = r_vec[0];

    std::vector<double> K0 = tetrahedron.vec_0;
    std::vector<double> K1 = tetrahedron.vec_0;

    double a0 = K0[0];
    double a1 = K1[0];

    std::complex<double> result(0., 0.);

    if(abs(r)<1.e-6)
      {
        real(result) = a1-a0;
        imag(result) = 0;
      }
    else
      {
        real(result) = -Sin(a0*r)/r + Sin(a1*r)/r;
        imag(result) =  Cos(a0*r)/r - Cos(a1*r)/r;
      }

    return result;
  }

  /************************************
   ***
   ***   2D harmonic-integration
   ***
   ************************************/

  /*!
   *   We want to integrate a harmonic function e^{i k r} over an arbitrary triangle, defined by k_0, k_1, k_2.
   *   This comes down to:
   *
   *   \int_0^1 dt_1 \int_0^{1-t_1} dt_2 e^{i r k} * det(d_1, d_2)
   *
   *   with
   *          k   = k_0 + d_1*t_1 + d_2*t_2
   *
   *          d_1 = k_1-k_0,
   *          d_2 = k_2-k_0.
   *          d_3 = d_2-d_1 = k_2-k_1
   *
   *   for convinience, we define
   *
   *          r_dot_di = r \dot d_i
   *
   *   there are 4 cases:
   *
   *    case 0:  |r_dot_d1|>0 and |r_dot_d2|>0  and |r_dot_d3|>0 (general)
   *    case 1:  |r_dot_d1|=0 and |r_dot_d2|>0  and |r_dot_d3|>0 (special 1)
   *    case 2:  |r_dot_d1|>0 and |r_dot_d2|=0  and |r_dot_d3|>0 (special 2)
   *    case 3:  |r_dot_d1|>0 and |r_dot_d2|>0  and |r_dot_d3|=0 (special 3)
   *
   *
   */
  std::complex<double> tetrahedron_routines_harmonic_function::execute(std::vector<double>&             r_vec,
                                                                       MATH_ALGORITHMS::tetrahedron<2>& tetrahedron)
  {
    assert(r_vec.size()==2);

    const static std::complex<double> I(0,1);
    const static double EPSILON = 1.e-6;

    std::complex<double> result(0., 0.);

    if(abs(VECTOR_OPERATIONS::DOT_PRODUCT(r_vec, r_vec))<EPSILON)
      {
        real(result) = tetrahedron.volume;
        imag(result) = 0;
      }
    else
      {
        std::vector<double> K0 = tetrahedron.vec_0;

        std::vector<double> D1        = K0;
        std::vector<double> D2        = K0;
        std::vector<double> D2_min_D1 = K0;

        for(int d=0; d<2; d++){
          D1[d] = tetrahedron.vec_1[d] - tetrahedron.vec_0[d];
          D2[d] = tetrahedron.vec_2[d] - tetrahedron.vec_0[d];

          D2_min_D1[d] = D2[d]-D1[d];
        }

        double dot_R_K0 = VECTOR_OPERATIONS::DOT_PRODUCT(r_vec, K0);
        double dot_R_D1 = VECTOR_OPERATIONS::DOT_PRODUCT(r_vec, D1);
        double dot_R_D2 = VECTOR_OPERATIONS::DOT_PRODUCT(r_vec, D2);

        double dot_R_D2_min_D1 = VECTOR_OPERATIONS::DOT_PRODUCT(r_vec, D2_min_D1);

        double               det   = VECTOR_OPERATIONS::VOLUME(D1, D2);
        std::complex<double> phase = Cos(dot_R_K0) + I*Sin(dot_R_K0);

        if(abs(dot_R_D2_min_D1)>EPSILON)
          {
            if(abs(dot_R_D1)>EPSILON and abs(dot_R_D2)>EPSILON)
              result = case_2D(dot_R_D1, dot_R_D2, dot_R_D2_min_D1)*phase*det;

            if(abs(dot_R_D1)<EPSILON)
              result = case_d1_2D(dot_R_D1, dot_R_D2, dot_R_D2_min_D1)*phase*det;

            if(abs(dot_R_D2)<EPSILON)
              result = case_d2_2D(dot_R_D1, dot_R_D2, dot_R_D2_min_D1)*phase*det;
          }
        else
          //        if(abs(dot_R_D2_min_D1)<EPSILON)
          {
            MATH_ALGORITHMS::tetrahedron<2> tetrahedron_new;

            permute(tetrahedron_new, tetrahedron);

            result = execute(r_vec, tetrahedron_new);
            //result = case_3_2D(dot_R_D1, dot_R_D2, dot_R_D3)*phase*det;
          }
      }

    return result;
  }

  void tetrahedron_routines_harmonic_function::permute(MATH_ALGORITHMS::tetrahedron<2>& tetrahedron_new,
                                                       MATH_ALGORITHMS::tetrahedron<2>& tetrahedron_old)
  {
    tetrahedron_new.vec_1 = tetrahedron_old.vec_0;
    tetrahedron_new.vec_2 = tetrahedron_old.vec_1;
    tetrahedron_new.vec_0 = tetrahedron_old.vec_2;
  }

  std::complex<double> tetrahedron_routines_harmonic_function::case_2D(double dotRD1,
                                                                       double dotRD2,
                                                                       double dotRD2minD1)
  {
    assert(abs(dotRD1) > 1.e-6);
    assert(abs(dotRD2) > 1.e-6);
    assert(abs(dotRD2minD1) > 1.e-6);

    std::complex<double> result(0., 0.);

    real(result) = -(1/(dotRD1*dotRD2minD1)) + 1/(dotRD2*dotRD2minD1) + Cos(dotRD1)/(dotRD1*dotRD2minD1) - Cos(dotRD2)/(dotRD2*dotRD2minD1);
    imag(result) = Sin(dotRD1)/(dotRD1*dotRD2minD1) - Sin(dotRD2)/(dotRD2*dotRD2minD1);

    assert(real(result)==real(result));
    assert(imag(result)==imag(result));

    return result;
  }

  std::complex<double> tetrahedron_routines_harmonic_function::case_d1_2D(double dotRD1,
                                                                          double dotRD2,
                                                                          double dotRD2minD1)
  {
    assert(abs(dotRD1) < 1.e-6);
    assert(abs(dotRD2) > 1.e-6);
    assert(abs(dotRD2minD1) > 1.e-6);

    std::complex<double> result(0., 0.);

    real(result) = 1/(dotRD2*dotRD2minD1) - Cos(dotRD2)/(dotRD2*dotRD2minD1);
    imag(result) = 1/dotRD2minD1 - Sin(dotRD2)/(dotRD2*dotRD2minD1);

    assert(real(result)==real(result));
    assert(imag(result)==imag(result));

    return result;
  }

  std::complex<double> tetrahedron_routines_harmonic_function::case_d2_2D(double dotRD1,
                                                                          double dotRD2,
                                                                          double dotRD2minD1)
  {
    assert(abs(dotRD1) > 1.e-6);
    assert(abs(dotRD2) < 1.e-6);
    assert(abs(dotRD2minD1) > 1.e-6);

    std::complex<double> result(0., 0.);

    real(result) = -(1/(dotRD1*dotRD2minD1)) + Cos(dotRD1)/(dotRD1*dotRD2minD1);
    imag(result) = -(1/dotRD2minD1) + Sin(dotRD1)/(dotRD1*dotRD2minD1);

    assert(real(result)==real(result));
    assert(imag(result)==imag(result));

    return result;
  }

  /*
    std::complex<double> tetrahedron_routines_harmonic_function::case_3_2D(double dotRD1,
    double dotRD2,
    double dotRD3)
    {
    assert(abs(dotRD1) > 1.e-6);
    assert(abs(dotRD2) > 1.e-6);
    assert(abs(dotRD3) < 1.e-6);

    std::complex<double> result(0., 0.);

    real(result) = -Power(dotRD2,-2) + Cos(dotRD2)/Power(dotRD2,2) + Sin(dotRD2)/dotRD2;
    imag(result) = -(Cos(dotRD2)/dotRD2) + Sin(dotRD2)/Power(dotRD2,2);

    return result;
    }
  */

  /************************************
   ***
   ***   3D harmonic-integration
   ***
   ************************************/

  void tetrahedron_routines_harmonic_function::permute(MATH_ALGORITHMS::tetrahedron<3>& tetrahedron_new,
                                                       MATH_ALGORITHMS::tetrahedron<3>& tetrahedron_old)
  {
    tetrahedron_new.vec_1 = tetrahedron_old.vec_0;
    tetrahedron_new.vec_2 = tetrahedron_old.vec_1;
    tetrahedron_new.vec_3 = tetrahedron_old.vec_2;
    tetrahedron_new.vec_0 = tetrahedron_old.vec_3;
  }

  std::complex<double> tetrahedron_routines_harmonic_function::case_3D(double dotRD1     , double dotRD2     , double dotRD3,
                                                                       double dotRD2minD1, double dotRD3minD2, double dotRD1minD3)
  {
    assert(abs(dotRD1) > 1.e-6);
    assert(abs(dotRD2) > 1.e-6);
    assert(abs(dotRD3) > 1.e-6);

    assert(abs(dotRD2minD1) > 1.e-6);
    assert(abs(dotRD3minD2) > 1.e-6);
    assert(abs(dotRD1minD3) > 1.e-6);

    std::complex<double> result(0., 0.);

    real(result) =
      -((dotRD2*Sin(dotRD1))/(dotRD1*dotRD1minD3*dotRD2minD1*dotRD3minD2)) +
      (dotRD3*Sin(dotRD1))/(dotRD1*dotRD1minD3*dotRD2minD1*dotRD3minD2) + Sin(dotRD2)/(dotRD2*dotRD2minD1*dotRD3minD2) +
      Sin(dotRD3)/(dotRD1minD3*dotRD3*dotRD3minD2);

    imag(result) =
      -(1/(dotRD1*dotRD2*dotRD3minD2)) + 1/(dotRD1*dotRD3*dotRD3minD2) +
      (dotRD2*Cos(dotRD1))/(dotRD1*dotRD1minD3*dotRD2minD1*dotRD3minD2) -
      (dotRD3*Cos(dotRD1))/(dotRD1*dotRD1minD3*dotRD2minD1*dotRD3minD2) - Cos(dotRD2)/(dotRD2*dotRD2minD1*dotRD3minD2) -
      Cos(dotRD3)/(dotRD1minD3*dotRD3*dotRD3minD2);

    assert(real(result)==real(result));
    assert(imag(result)==imag(result));

    return result;
  }

  std::complex<double> tetrahedron_routines_harmonic_function::case_d1_3D(double dotRD1     , double dotRD2     , double dotRD3,
                                                                          double dotRD2minD1, double dotRD3minD2, double dotRD1minD3)
  {
    assert(abs(dotRD1) < 1.e-6);
    assert(abs(dotRD2) > 1.e-6);
    assert(abs(dotRD3) > 1.e-6);

    //     assert(abs(dotRD2minD1) > 1.e-6);
    //     assert(abs(dotRD3minD2) > 1.e-6);
    //     assert(abs(dotRD1minD3) > 1.e-6);

    std::complex<double> result(0., 0.);

    if(abs(dotRD3minD2) > 1.e-6)
      {
        real(result) =
          -(1/(dotRD2*dotRD3minD2)) + 1/(dotRD3*dotRD3minD2) + Sin(dotRD2)/(Power(dotRD2,2)*dotRD3minD2) -
          Sin(dotRD3)/(Power(dotRD3,2)*dotRD3minD2);

        imag(result) =
          1/(Power(dotRD2,2)*dotRD3minD2) - 1/(Power(dotRD3,2)*dotRD3minD2) - Cos(dotRD2)/(Power(dotRD2,2)*dotRD3minD2) +
          Cos(dotRD3)/(Power(dotRD3,2)*dotRD3minD2);
      }
    else
      {
	//cout << __FUNCTION__ << " needs implementation\n";

	real(result) =
	  -Power(dotRD2,-2) - Cos(dotRD2)/Power(dotRD2,2) + (2*Sin(dotRD2))/Power(dotRD2,3);

	imag(result) =
	  2/Power(dotRD2,3) - (2*Cos(dotRD2))/Power(dotRD2,3) - Sin(dotRD2)/Power(dotRD2,2);
      }

    assert(real(result)==real(result));
    assert(imag(result)==imag(result));

    return result;
  }

  std::complex<double> tetrahedron_routines_harmonic_function::case_d2_3D(double dotRD1     , double dotRD2     , double dotRD3,
                                                                          double dotRD2minD1, double dotRD3minD2, double dotRD1minD3)
  {
    assert(abs(dotRD1) > 1.e-6);
    assert(abs(dotRD2) < 1.e-6);
    assert(abs(dotRD3) > 1.e-6);

    //     assert(abs(dotRD2minD1) > 1.e-6);
    //     assert(abs(dotRD3minD2) > 1.e-6);
    //     assert(abs(dotRD1minD3) > 1.e-6);

    std::complex<double> result(0., 0.);

    if(abs(dotRD1minD3) > 1.e-6)
      {
        real(result) =
          1/(dotRD1*dotRD1minD3) - 1/(dotRD1minD3*dotRD3) - Sin(dotRD1)/(Power(dotRD1,2)*dotRD1minD3) +
          Sin(dotRD3)/(dotRD1minD3*Power(dotRD3,2));

        imag(result) =
          -(1/(Power(dotRD1,2)*dotRD1minD3)) + 1/(dotRD1minD3*Power(dotRD3,2)) + Cos(dotRD1)/(Power(dotRD1,2)*dotRD1minD3) -
          Cos(dotRD3)/(dotRD1minD3*Power(dotRD3,2));
      }
    else
      {
	//cout << __FUNCTION__ << " needs implementation\n";

	real(result) =
	  -Power(dotRD3,-2) - Cos(dotRD3)/Power(dotRD3,2) + (2*Sin(dotRD3))/Power(dotRD3,3);

	imag(result) =
	  2/Power(dotRD3,3) - (2*Cos(dotRD3))/Power(dotRD3,3) - Sin(dotRD3)/Power(dotRD3,2);
      }

    assert(real(result)==real(result));
    assert(imag(result)==imag(result));

    return result;
  }

  std::complex<double> tetrahedron_routines_harmonic_function::case_d3_3D(double dotRD1     , double dotRD2     , double dotRD3,
                                                                          double dotRD2minD1, double dotRD3minD2, double dotRD1minD3)
  {
    assert(abs(dotRD1) > 1.e-6);
    assert(abs(dotRD2) > 1.e-6);
    assert(abs(dotRD3) < 1.e-6);

    //     assert(abs(dotRD2minD1) > 1.e-6);
    //     assert(abs(dotRD3minD2) > 1.e-6);
    //     assert(abs(dotRD1minD3) > 1.e-6);

    std::complex<double> result(0., 0.);

    if(abs(dotRD2minD1) > 1.e-6)
      {
        real(result) =
          -(1/(dotRD1*dotRD2minD1)) + 1/(dotRD2*dotRD2minD1) + Sin(dotRD1)/(Power(dotRD1,2)*dotRD2minD1) -
          Sin(dotRD2)/(Power(dotRD2,2)*dotRD2minD1);

        imag(result) =
          1/(Power(dotRD1,2)*dotRD2minD1) - 1/(Power(dotRD2,2)*dotRD2minD1) - Cos(dotRD1)/(Power(dotRD1,2)*dotRD2minD1) +
          Cos(dotRD2)/(Power(dotRD2,2)*dotRD2minD1);
      }
    else
      {
	//cout << __FUNCTION__ << " needs implementation\n";

	real(result) = 
	  -Power(dotRD1,-2) - Cos(dotRD1)/Power(dotRD1,2) + (2*Sin(dotRD1))/Power(dotRD1,3);

	imag(result) = 
	  2/Power(dotRD1,3) - (2*Cos(dotRD1))/Power(dotRD1,3) - Sin(dotRD1)/Power(dotRD1,2);
      }

    assert(real(result)==real(result));
    assert(imag(result)==imag(result));

    return result;
  }

  std::complex<double> tetrahedron_routines_harmonic_function::case_d1_d2_3D(double dotRD1     , double dotRD2     , double dotRD3,
                                                                             double dotRD2minD1, double dotRD3minD2, double dotRD1minD3)
  {
    assert(abs(dotRD1) < 1.e-6);
    assert(abs(dotRD2) < 1.e-6);
    assert(abs(dotRD3) > 1.e-6);

    //     assert(abs(dotRD2minD1) > 1.e-6);
    //     assert(abs(dotRD3minD2) > 1.e-6);
    //     assert(abs(dotRD1minD3) > 1.e-6);

    std::complex<double> result(0., 0.);

    real(result) =
      Power(dotRD3,-2) - Sin(dotRD3)/Power(dotRD3,3);

    imag(result) =
      -Power(dotRD3,-3) + 1/(2.*dotRD3) + Cos(dotRD3)/Power(dotRD3,3);

    assert(real(result)==real(result));
    assert(imag(result)==imag(result));

    return result;
  }

  std::complex<double> tetrahedron_routines_harmonic_function::case_d2_d3_3D(double dotRD1     , double dotRD2     , double dotRD3,
                                                                             double dotRD2minD1, double dotRD3minD2, double dotRD1minD3)
  {
    assert(abs(dotRD1) > 1.e-6);
    assert(abs(dotRD2) < 1.e-6);
    assert(abs(dotRD3) < 1.e-6);

    //     assert(abs(dotRD2minD1) > 1.e-6);
    //     assert(abs(dotRD3minD2) > 1.e-6);
    //     assert(abs(dotRD1minD3) > 1.e-6);

    std::complex<double> result(0., 0.);

    real(result) =
      Power(dotRD1,-2) - Sin(dotRD1)/Power(dotRD1,3);

    imag(result) =
      -Power(dotRD1,-3) + 1/(2.*dotRD1) + Cos(dotRD1)/Power(dotRD1,3);

    assert(real(result)==real(result));
    assert(imag(result)==imag(result));

    return result;
  }

  std::complex<double> tetrahedron_routines_harmonic_function::case_d3_d1_3D(double dotRD1     , double dotRD2     , double dotRD3,
                                                                             double dotRD2minD1, double dotRD3minD2, double dotRD1minD3)
  {
    assert(abs(dotRD1) < 1.e-6);
    assert(abs(dotRD2) > 1.e-6);
    assert(abs(dotRD3) < 1.e-6);

    //     assert(abs(dotRD2minD1) > 1.e-6);
    //     assert(abs(dotRD3minD2) > 1.e-6);
    //     assert(abs(dotRD1minD3) > 1.e-6);

    std::complex<double> result(0., 0.);

    real(result) =
      Power(dotRD2,-2) - Sin(dotRD2)/Power(dotRD2,3);

    imag(result) =
      -Power(dotRD2,-3) + 1/(2.*dotRD2) + Cos(dotRD2)/Power(dotRD2,3);

    assert(real(result)==real(result));
    assert(imag(result)==imag(result));

    return result;
  }

  std::complex<double> tetrahedron_routines_harmonic_function::execute(std::vector<double>&             r_vec,
                                                                       MATH_ALGORITHMS::tetrahedron<3>& tetrahedron)
  {
    assert(r_vec.size()==3);

    const static std::complex<double> I(0,1);
    const static double EPSILON = 1.e-6;

    std::complex<double> result(0., 0.);

    if(abs(VECTOR_OPERATIONS::DOT_PRODUCT(r_vec, r_vec))<EPSILON)
      {
        real(result) = tetrahedron.volume;
        imag(result) = 0;
      }
    else
      {
        std::vector<double> K_0 = tetrahedron.vec_0;

        std::vector<double> D1 = K_0;
        std::vector<double> D2 = K_0;
        std::vector<double> D3 = K_0;

        std::vector<double> D2_min_D1 = K_0;
        std::vector<double> D3_min_D2 = K_0;
        std::vector<double> D1_min_D3 = K_0;

        for(int d=0; d<3; d++){
          D1[d] = tetrahedron.vec_1[d] - tetrahedron.vec_0[d];
          D2[d] = tetrahedron.vec_2[d] - tetrahedron.vec_0[d];
          D3[d] = tetrahedron.vec_3[d] - tetrahedron.vec_0[d];

          D2_min_D1[d] = D2[d]-D1[d];
          D3_min_D2[d] = D3[d]-D2[d];
          D1_min_D3[d] = D1[d]-D3[d];
        }

        double dot_R_K0 = VECTOR_OPERATIONS::DOT_PRODUCT(r_vec, K_0);

        double dot_R_D1 = VECTOR_OPERATIONS::DOT_PRODUCT(r_vec, D1);
        double dot_R_D2 = VECTOR_OPERATIONS::DOT_PRODUCT(r_vec, D2);
        double dot_R_D3 = VECTOR_OPERATIONS::DOT_PRODUCT(r_vec, D3);

        double dot_R_D2_min_D1 = VECTOR_OPERATIONS::DOT_PRODUCT(r_vec, D2_min_D1);
        double dot_R_D3_min_D2 = VECTOR_OPERATIONS::DOT_PRODUCT(r_vec, D3_min_D2);
        double dot_R_D1_min_D3 = VECTOR_OPERATIONS::DOT_PRODUCT(r_vec, D1_min_D3);

        double               det   = VECTOR_OPERATIONS::VOLUME(D1, D2, D3);
        std::complex<double> phase = Cos(dot_R_K0) + I*Sin(dot_R_K0);

        if(abs(dot_R_D1)>EPSILON and
           abs(dot_R_D2)>EPSILON and
           abs(dot_R_D3)>EPSILON and
           abs(dot_R_D2_min_D1)>EPSILON and
           abs(dot_R_D3_min_D2)>EPSILON and
           abs(dot_R_D1_min_D3)>EPSILON) // general case
          {
            result = case_3D(dot_R_D1, dot_R_D2, dot_R_D3, dot_R_D2_min_D1, dot_R_D3_min_D2, dot_R_D1_min_D3)*det*phase;
          }
        else
          {
            if(abs(dot_R_D1)<EPSILON or
               abs(dot_R_D2)<EPSILON or
               abs(dot_R_D3)<EPSILON) // special cases where one or two dot-products are zero
              {
                if(abs(dot_R_D1)<EPSILON and
                   abs(dot_R_D2)>EPSILON and
                   abs(dot_R_D3)>EPSILON)
                  result = case_d1_3D(dot_R_D1, dot_R_D2, dot_R_D3, dot_R_D2_min_D1, dot_R_D3_min_D2, dot_R_D1_min_D3)*det*phase;

                if(abs(dot_R_D1)>EPSILON and
                   abs(dot_R_D2)<EPSILON and
                   abs(dot_R_D3)>EPSILON)
                  result = case_d2_3D(dot_R_D1, dot_R_D2, dot_R_D3, dot_R_D2_min_D1, dot_R_D3_min_D2, dot_R_D1_min_D3)*det*phase;

                if(abs(dot_R_D1)>EPSILON and
                   abs(dot_R_D2)>EPSILON and
                   abs(dot_R_D3)<EPSILON)
                  result = case_d3_3D(dot_R_D1, dot_R_D2, dot_R_D3, dot_R_D2_min_D1, dot_R_D3_min_D2, dot_R_D1_min_D3)*det*phase;

                if(abs(dot_R_D1)<EPSILON and
                   abs(dot_R_D2)<EPSILON and
                   abs(dot_R_D3)>EPSILON)
                  result = case_d1_d2_3D(dot_R_D1, dot_R_D2, dot_R_D3, dot_R_D2_min_D1, dot_R_D3_min_D2, dot_R_D1_min_D3)*det*phase;

                if(abs(dot_R_D1)>EPSILON and
                   abs(dot_R_D2)<EPSILON and
                   abs(dot_R_D3)<EPSILON)
                  result = case_d2_d3_3D(dot_R_D1, dot_R_D2, dot_R_D3, dot_R_D2_min_D1, dot_R_D3_min_D2, dot_R_D1_min_D3)*det*phase;

                if(abs(dot_R_D1)<EPSILON and
                   abs(dot_R_D2)>EPSILON and
                   abs(dot_R_D3)<EPSILON)
                  result = case_d3_d1_3D(dot_R_D1, dot_R_D2, dot_R_D3, dot_R_D2_min_D1, dot_R_D3_min_D2, dot_R_D1_min_D3)*det*phase;
              }
            else
              {
//                 cout << "\n\t start permuting\t";
//                 VECTOR_OPERATIONS::PRINT(r_vec);
//                 cout << "\n";
//                 VECTOR_OPERATIONS::PRINT(tetrahedron.vec_0);cout << "\n";
//                 VECTOR_OPERATIONS::PRINT(tetrahedron.vec_1);cout << "\n";
//                 VECTOR_OPERATIONS::PRINT(tetrahedron.vec_2);cout << "\n";
//                 VECTOR_OPERATIONS::PRINT(tetrahedron.vec_3);cout << "\n";

                MATH_ALGORITHMS::tetrahedron<3> tetrahedron_new;

                permute(tetrahedron_new, tetrahedron);

                result = execute(r_vec, tetrahedron_new);
              }
          }
      }

    assert(real(result)==real(result));
    assert(imag(result)==imag(result));

    return result;
  }

}

#endif
