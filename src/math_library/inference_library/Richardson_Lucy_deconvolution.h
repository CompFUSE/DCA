//-*-C++-*-

#ifndef INFERENCE_RICHARDSON_LUCY_H
#define INFERENCE_RICHARDSON_LUCY_H

namespace INFERENCE
{
  /*! \ingroup INFERENCE
   *
   *  \author Peter Staar
   *  \brief  This class implements the Richardson-Lucy deconvolution.
   *
   */
  template<typename parameters_type, typename k_dmn_t, typename p_dmn_t>
  class Richardson_Lucy_deconvolution
  {

  public:

    Richardson_Lucy_deconvolution(parameters_type& parameters_ref);
    ~Richardson_Lucy_deconvolution();

    void execute(LIN_ALG::matrix<double, LIN_ALG::CPU>& matrix,
                 function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_source,
                 function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_target);

    void execute(LIN_ALG::matrix<double, LIN_ALG::CPU>&      A,
                 function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_source,
                 function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_approx,
                 function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_target);

  private:

    void initialize_matrices(function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_source);

    void initialize_errors(function<bool  , p_dmn_t>& is_finished,
                           function<double, p_dmn_t>& error_function);

    bool update_f_target(function<bool  , p_dmn_t>& is_finished,
                         function<double, p_dmn_t>& error_function,
                         function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_target);

  private:

    parameters_type&  parameters;

    LIN_ALG::matrix<double, LIN_ALG::CPU> c;
    LIN_ALG::matrix<double, LIN_ALG::CPU> d;

    LIN_ALG::matrix<double, LIN_ALG::CPU> d_over_c;

    LIN_ALG::matrix<double, LIN_ALG::CPU> u_t;
    LIN_ALG::matrix<double, LIN_ALG::CPU> u_t_p_1;
  };

  template<typename parameters_type, typename k_dmn_t, typename p_dmn_t>
  Richardson_Lucy_deconvolution<parameters_type, k_dmn_t, p_dmn_t>::Richardson_Lucy_deconvolution(parameters_type& parameters_ref):
    parameters(parameters_ref),

    c("c (Richardson_Lucy_deconvolution)"),
    d("d (Richardson_Lucy_deconvolution)"),

    d_over_c("d/c (Richardson_Lucy_deconvolution)"),

    u_t("u_t (Richardson_Lucy_deconvolution)"),
    u_t_p_1("u_{t+1} (Richardson_Lucy_deconvolution)")
  {}

  template<typename parameters_type, typename k_dmn_t, typename p_dmn_t>
  Richardson_Lucy_deconvolution<parameters_type, k_dmn_t, p_dmn_t>::~Richardson_Lucy_deconvolution()
  {}

  template<typename parameters_type, typename k_dmn_t, typename p_dmn_t>
  void Richardson_Lucy_deconvolution<parameters_type, k_dmn_t, p_dmn_t>::execute(LIN_ALG::matrix<double, LIN_ALG::CPU>&      A,
                                                                                 function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_source,
                                                                                 function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_target)
  {
    //cout << __FUNCTION__ << endl;

    assert(A.get_current_size().first==k_dmn_t::dmn_size());
    assert(A.get_current_size().first==A.get_current_size().second);

    function<bool  , p_dmn_t> is_finished("is_finished");
    function<double, p_dmn_t> error_function("error_function");

    initialize_matrices(f_source);

    // compute c
    LIN_ALG::GEMM<LIN_ALG::CPU>::execute(A, u_t, c);

    initialize_errors(is_finished, error_function);

    for(int l=0; l<parameters.get_deconvolution_iterations(); l++)
      {
        //cout << l << endl;

        for(int j=0; j<p_dmn_t::dmn_size(); j++)
          for(int i=0; i<k_dmn_t::dmn_size(); i++)
            d_over_c(i,j) = d(i,j)/c(i,j);

        // compute u_t_plus_1
        LIN_ALG::GEMM<LIN_ALG::CPU>::execute('T', 'N', A, d_over_c, u_t_p_1);

        for(int j=0; j<p_dmn_t::dmn_size(); j++)
          for(int i=0; i<k_dmn_t::dmn_size(); i++)
            u_t(i,j) = u_t_p_1(i,j)*u_t(i,j);

        // compute c
        LIN_ALG::GEMM<LIN_ALG::CPU>::execute(A, u_t, c);

        bool finished = update_f_target(is_finished, error_function, f_target);

        if(finished)
          break;
      }

    for(int j=0; j<p_dmn_t::dmn_size(); j++)
      if(not is_finished(j))
        for(int i=0; i<k_dmn_t::dmn_size(); i++)
          f_target(i,j) = u_t(i,j);
  }

  template<typename parameters_type, typename k_dmn_t, typename p_dmn_t>
  void Richardson_Lucy_deconvolution<parameters_type, k_dmn_t, p_dmn_t>::execute(LIN_ALG::matrix<double, LIN_ALG::CPU>&      A,
                                                                                 function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_source,
                                                                                 function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_approx,
                                                                                 function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_target)
  {
    //cout << __FUNCTION__ << endl;

    execute(A, f_source, f_target);

    for(int j=0; j<p_dmn_t::dmn_size(); j++)
      for(int i=0; i<k_dmn_t::dmn_size(); i++)
        u_t(i,j) = f_target(i,j);

    LIN_ALG::GEMM<LIN_ALG::CPU>::execute(A, u_t, c);

    for(int j=0; j<p_dmn_t::dmn_size(); j++)
      for(int i=0; i<k_dmn_t::dmn_size(); i++)
        f_approx(i,j) = c(i,j);
  }

  template<typename parameters_type, typename k_dmn_t, typename p_dmn_t>
  void Richardson_Lucy_deconvolution<parameters_type, k_dmn_t, p_dmn_t>::initialize_matrices(function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_source)
  {
    int nr_rows = k_dmn_t::dmn_size();
    int nr_cols = p_dmn_t::dmn_size();

    c.resize_no_copy(std::pair<int,int>(nr_rows, nr_cols));
    d.resize_no_copy(std::pair<int,int>(nr_rows, nr_cols));

    d_over_c.resize_no_copy(std::pair<int,int>(nr_rows, nr_cols));

    u_t    .resize_no_copy(std::pair<int,int>(nr_rows, nr_cols));
    u_t_p_1.resize_no_copy(std::pair<int,int>(nr_rows, nr_cols));

    for(int j=0; j<nr_cols; j++)
      for(int i=0; i<nr_rows; i++)
        d(i,j) = f_source(i,j);

    for(int j=0; j<nr_cols; j++){

      double mean=0;
      for(int i=0; i<nr_rows; i++)
        mean += f_source(i,j);
      mean /= nr_rows;

      for(int i=0; i<nr_rows; i++)
        u_t(i,j) = mean;
    }
  }

  template<typename parameters_type, typename k_dmn_t, typename p_dmn_t>
  void Richardson_Lucy_deconvolution<parameters_type, k_dmn_t, p_dmn_t>::initialize_errors(function<bool  , p_dmn_t>& is_finished,
                                                                                           function<double, p_dmn_t>& error_function)
  {
    for(int j=0; j<p_dmn_t::dmn_size(); j++)
      {
        double error=0;
        for(int i=0; i<k_dmn_t::dmn_size(); i++)
          error += std::pow(c(i,j)-d(i,j), 2);

        is_finished   (j) = false;
        error_function(j) = sqrt(error)/k_dmn_t::dmn_size();
      }
  }

  template<typename parameters_type, typename k_dmn_t, typename p_dmn_t>
  bool Richardson_Lucy_deconvolution<parameters_type, k_dmn_t, p_dmn_t>::update_f_target(function<bool  , p_dmn_t>& is_finished,
                                                                                         function<double, p_dmn_t>& error_function,
                                                                                         function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_target)
  {
    bool all_are_finished = true;

    double epsilon = parameters.get_deconvolution_tolerance();

    for(int j=0; j<p_dmn_t::dmn_size(); j++)
      {
        if(not is_finished(j))
          {
            double error=0;
            for(int i=0; i<k_dmn_t::dmn_size(); i++)
              error += std::pow(c(i,j)-d(i,j), 2);
            error = sqrt(error)/k_dmn_t::dmn_size();

//             if(j==p_dmn_t::dmn_size()/2+1){
//               cout << error << "\t" << error_function(j) << "\t" << (error-error_function(j))/(error+1.e-12) << "\n";
//             }

            if( abs((error-error_function(j))/(error+1.e-12)) < epsilon)
              {
                for(int i=0; i<k_dmn_t::dmn_size(); i++)
                  f_target(i,j) = u_t(i,j);

                is_finished(j) = true;
              }
            else
              {
                all_are_finished = false;

                error_function(j) = error;
              }
          }
      }

    return all_are_finished;
  }

}

#endif
