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
                 FUNC_LIB::function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_source,
                 FUNC_LIB::function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_target);

    void execute(LIN_ALG::matrix<double, LIN_ALG::CPU>&      A,
                 FUNC_LIB::function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_source,
                 FUNC_LIB::function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_approx,
                 FUNC_LIB::function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_target);

  private:

    void initialize_matrices(FUNC_LIB::function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_source);

    void initialize_errors(FUNC_LIB::function<bool  , p_dmn_t>& is_finished,
                           FUNC_LIB::function<double, p_dmn_t>& error_function);

    bool update_f_target(FUNC_LIB::function<bool  , p_dmn_t>& is_finished,
                         FUNC_LIB::function<double, p_dmn_t>& error_function,
                         FUNC_LIB::function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_target);

  private:

    parameters_type&  parameters;
    typename parameters_type::concurrency_type& concurrency;


    LIN_ALG::matrix<double, LIN_ALG::CPU> c;
    LIN_ALG::matrix<double, LIN_ALG::CPU> d;

    LIN_ALG::matrix<double, LIN_ALG::CPU> d_over_c;

    LIN_ALG::matrix<double, LIN_ALG::CPU> u_t;
    LIN_ALG::matrix<double, LIN_ALG::CPU> u_t_p_1;
  };

  template<typename parameters_type, typename k_dmn_t, typename p_dmn_t>
  Richardson_Lucy_deconvolution<parameters_type, k_dmn_t, p_dmn_t>::Richardson_Lucy_deconvolution(parameters_type& parameters_ref):
    parameters(parameters_ref),
    concurrency(parameters.get_concurrency()),

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
                                                                                 FUNC_LIB::function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_source,
                                                                                 FUNC_LIB::function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_target)
  {
    if (concurrency.id() == 0)
      {
        std::cout << "\n\nRichardson_Lucy_deconvolution: " << __FUNCTION__ << std::endl;
        std::cout << "It\tw_ind\tError_Re\tError_Im" << std::endl;
      }

    assert(A.get_current_size().first==k_dmn_t::dmn_size());
    assert(A.get_current_size().first==A.get_current_size().second);

    FUNC_LIB::function<bool  , p_dmn_t> is_finished("is_finished");
    FUNC_LIB::function<double, p_dmn_t> error_function("error_function");

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

        // FIXME: This is hardcoded for 256 fermionic frequencies.
        if (concurrency.id() == 0)
          {
            std::cout << l << "\t" << "253" << "\t" << error_function(0, 0, 0, 0, 253) << "\t" << error_function(1, 0, 0, 0, 253) << std::endl;
            std::cout << l << "\t" << "254" << "\t" << error_function(0, 0, 0, 0, 254) << "\t" << error_function(1, 0, 0, 0, 254) << std::endl;
            std::cout << l << "\t" << "255" << "\t" << error_function(0, 0, 0, 0, 255) << "\t" << error_function(1, 0, 0, 0, 255) << std::endl;
            std::cout << l << "\t" << "256" << "\t" << error_function(0, 0, 0, 0, 256) << "\t" << error_function(1, 0, 0, 0, 256) << std::endl;
            std::cout << l << "\t" << "257" << "\t" << error_function(0, 0, 0, 0, 257) << "\t" << error_function(1, 0, 0, 0, 257) << std::endl;
            std::cout << l << "\t" << "258" << "\t" << error_function(0, 0, 0, 0, 258) << "\t" << error_function(1, 0, 0, 0, 258) << std::endl;
            std::cout << std::endl;
          }

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
                                                                                 FUNC_LIB::function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_source,
                                                                                 FUNC_LIB::function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_approx,
                                                                                 FUNC_LIB::function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_target)
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
  void Richardson_Lucy_deconvolution<parameters_type, k_dmn_t, p_dmn_t>::initialize_matrices(FUNC_LIB::function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_source)
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
  void Richardson_Lucy_deconvolution<parameters_type, k_dmn_t, p_dmn_t>::initialize_errors(FUNC_LIB::function<bool  , p_dmn_t>& is_finished,
                                                                                           FUNC_LIB::function<double, p_dmn_t>& error_function)
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
  bool Richardson_Lucy_deconvolution<parameters_type, k_dmn_t, p_dmn_t>::update_f_target(FUNC_LIB::function<bool  , p_dmn_t>& is_finished,
                                                                                         FUNC_LIB::function<double, p_dmn_t>& error_function,
                                                                                         FUNC_LIB::function<double, dmn_2<k_dmn_t, p_dmn_t> >& f_target)
  {
    bool all_are_finished = true;

    double epsilon = parameters.get_deconvolution_tolerance();

    for(int j=0; j<p_dmn_t::dmn_size(); j++)
      {
        if(not is_finished(j))
          {
            double diff = 0;
            double tot  = 1.e-6;

            for(int i=0; i<k_dmn_t::dmn_size(); i++){
              diff += std::pow(c(i,j)-d(i,j), 2);
              tot  += std::pow(d(i,j), 2);
            }

            error_function(j) = std::sqrt(diff/tot);

            if(error_function(j) < epsilon)
              {
                for(int i=0; i<k_dmn_t::dmn_size(); i++)
                  f_target(i,j) = u_t(i,j);

                is_finished(j) = true;
              }

            else
              {
                all_are_finished = false;
              }
          }
      }

    return all_are_finished;
  }

}

#endif
