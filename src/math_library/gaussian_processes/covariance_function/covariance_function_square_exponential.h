//-*-C++-*-

#ifndef COVARIANCE_FUNCTION_SQUARED_EXPONENTIAL_H
#define COVARIANCE_FUNCTION_SQUARED_EXPONENTIAL_H

namespace MATH_LIBRARY
{
  /*!
   *
   *    page 19, eqn 2.31, Rasmussen and Williams
   */
  template<typename k_dmn_t>
  class covariance_function<SQUARED_EXPONENTIAL, k_dmn_t>
  {
    const static int DIMENSION;// = k_dmn_t::parameter_type::DIMENSION;

  public:

    covariance_function();
    ~covariance_function();

    double execute(std::vector<double>& x_i);

    void plot();

  public:

    double sigma_f;
    
    LIN_ALG::matrix<double, LIN_ALG::CPU> A;
  };

  template<typename k_dmn_t>
  const int covariance_function<SQUARED_EXPONENTIAL, k_dmn_t>::DIMENSION = k_dmn_t::parameter_type::DIMENSION;

  template<typename k_dmn_t>
  covariance_function<SQUARED_EXPONENTIAL, k_dmn_t>::covariance_function():
    sigma_f(1.),

    A("A", std::pair<int,int>(DIMENSION,DIMENSION))
  {
    for(int li=0; li<DIMENSION; li++)
      A(li,li) = 1.;
  }

  template<typename k_dmn_t>
  covariance_function<SQUARED_EXPONENTIAL, k_dmn_t>::~covariance_function()
  {}

  template<typename k_dmn_t>
  double covariance_function<SQUARED_EXPONENTIAL, k_dmn_t>::execute(std::vector<double>& x_i)
  {
    std::vector<double> y_i(DIMENSION, 0.);

    for(int li=0; li<DIMENSION; li++)
      for(int lj=0; lj<DIMENSION; lj++)
        y_i[li] += k_dmn_t::parameter_type::get_inverse_basis()[li+lj*DIMENSION]*x_i[lj];

    double result=0;

    for(int li=0; li<DIMENSION; li++)
      for(int lj=0; lj<DIMENSION; lj++)
        result += y_i[li]*A(li,lj)*y_i[lj];

    return sigma_f*std::exp(-result);
  }

  template<typename k_dmn_t>
  void covariance_function<SQUARED_EXPONENTIAL, k_dmn_t>::plot()
  {
    std::vector<double> x(0);
    std::vector<double> y(0);
    std::vector<double> z(0);

    std::vector<double> vec(DIMENSION,0);

    switch(DIMENSION)
      {
      case 1 :
        {
	  for(int l=0; l<k_dmn_t::dmn_size(); l++){

	    vec[0] = k_dmn_t::get_elements()[l];

	    x.push_back(vec[0]);	    
	    z.push_back(execute(vec));
	  }

          SHOW::plot_points(x,z);
        }
        break;

//       case 2 :
//         {
// 	  for(int l=0; l<k_dmn_t::dmn_size(); l++){
	    
// 	    vec[0] = k_dmn_t::get_elements()[l][0];	    
// 	    vec[1] = k_dmn_t::get_elements()[l][1];
	    
// 	    x.push_back(vec[0]);
// 	    y.push_back(vec[1]);
	    
// 	    z.push_back(execute(vec));
// 	  }

//           SHOW::heatmap(x,y,z);
//         }
//         break;

      default:
        cout << __FUNCTION__ << endl;
      }

  }

}

#endif
