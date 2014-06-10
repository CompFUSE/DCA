//-*-C++-*-

#ifndef COMPUTE_SPECTRUM_H
#define COMPUTE_SPECTRUM_H

namespace DCA
{
  /*!
   *  \defgroup SPECTRUM
   *  \ingroup  ANALYSIS
   */

  /*! \class   compute_spectrum
   *  \ingroup SPECTRUM
   *
   *  \author  Peter Staar
   *  \brief   This class organizes the compution of the spectrum, using a CPE analytic continuation.
   *  \version 1.0
   */
  template<class parameters_type, class basis_function_t>
  class compute_spectrum
  {
#include "type_definitions.h"

  public:

    typedef typename parameters_type::profiler_type    profiler_type;
    typedef typename parameters_type::concurrency_type concurrency_type;

  public:

    compute_spectrum(parameters_type& parameters);
    ~compute_spectrum();

    template<typename MOMS_imag_type, typename MOMS_real_type>
    void write(std::string     file_name,
               MOMS_imag_type& MOMS_imag,
               MOMS_real_type& MOMS_real);

    template<IO::FORMAT DATA_FORMAT>
    void write(IO::writer<DATA_FORMAT>& writer);

    template<typename MOMS_imag_type, typename MOMS_real_type>
    void execute(MOMS_imag_type& MOMS_imag,
                 MOMS_real_type& MOMS_real);

  private:

    template<typename MOMS_imag_type, typename MOMS_real_type>
    void execute_without_error_bars(MOMS_imag_type& MOMS_imag,
                                    MOMS_real_type& MOMS_real);

    template<typename MOMS_imag_type, typename MOMS_real_type>
    void execute_with_error_bars(MOMS_imag_type& MOMS_imag,
                                 MOMS_real_type& MOMS_real);


    template<typename MOMS_imag_type, typename MOMS_real_type>
    void test_A_w_versus_G_t(MOMS_imag_type& MOMS_imag,
                             MOMS_real_type& MOMS_real);

  private:

    void generate_f_original(function<std::complex<double>, dmn_4<nu,nu,k_DCA,w> >& S_K_w);

    void generate_new_sample(function<std::complex<double>, dmn_4<nu,nu,k_DCA,w> >& S_K_w_mean,
                             function<std::complex<double>, dmn_4<nu,nu,k_DCA,w> >& S_K_w_stddev,
                             function<std::complex<double>, dmn_4<nu,nu,k_DCA,w> >& S_K_w_sample);

    template<typename k_dmn_t, typename w_imag_dmn_t, typename w_real_dmn_t>
    void perform_analytic_continuation(function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_imag_dmn_t> >& S_k_w_imag,
                                       function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_real_dmn_t> >& S_k_w_real);

    template<typename k_dmn_t, typename w_dmn_t>
    void compute_G_k_w_on_cluster(function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& G0_k_w,
                                  function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& Sigma,
                                  function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& G_k_w);

    template<typename k_host_dmn_t, typename k_cluster_dmn_t, typename w_dmn_t>
    void compute_G_k_w_on_lattice(function<std::complex<double>, dmn_3<nu,nu,k_host_dmn_t           > >& H_k,
                                  function<std::complex<double>, dmn_4<nu,nu,k_cluster_dmn_t,w_dmn_t> >& Sigma,
                                  function<std::complex<double>, dmn_4<nu,nu,k_cluster_dmn_t,w_dmn_t> >& G_k_w);

    template<typename k_dmn_t, typename w_dmn_t>
    void compute_A_w(function<             double ,                     w_dmn_t>  & A_w,
                     function<             double , dmn_3<b,s,          w_dmn_t> >& A_nu_w,
                     function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& G_k_w);

    template<typename k_dmn_t>
    void compute_G_k_beta_over_2(function<double, dmn_3<b, s , k_dmn_t   > >& G_k_beta_over_2,
                                 function<double, dmn_4<nu,nu, k_dmn_t, t> >& G_k_t);

    template<typename k_dmn_t>
    void accumulate_integrated_A_k_div_cosh(function<             double , dmn_3<b , s, k_dmn_t       > >& integrated_A_k_div_cosh,
                                            function<             double , dmn_3<b , s, k_dmn_t       > >& integrated_A_k_div_cosh_stddev,
                                            function<std::complex<double>, dmn_4<nu,nu, k_dmn_t,w_REAL> >& G_k_w);

    template<typename scalar_type, typename dmn_type>
    void accumulate(function<scalar_type, dmn_type>& f,
                    function<scalar_type, dmn_type>& f_average,
                    function<scalar_type, dmn_type>& f_square);

    template<typename scalar_type, typename dmn_type>
    void accumulate(function<std::complex<scalar_type>, dmn_type>& f,
                    function<std::complex<scalar_type>, dmn_type>& f_average,
                    function<std::complex<scalar_type>, dmn_type>& f_square);

    template<typename k_dmn_t, typename w_dmn_t>
    void accumulate_A_w(function<             double,            w_dmn_t>&             A_w,
                        function<             double,            w_dmn_t>&             A_w_stddev,
                        function<             double,  dmn_3<b,s,w_dmn_t> >&           A_nu_w,
                        function<             double,  dmn_3<b,s,w_dmn_t> >&           A_nu_w_stddev,
                        function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& G_k_w);

    template<typename k_dmn_t, typename w_dmn_t>
    void accumulate_f_K_w(function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& f_K_w_sample,
                          function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& f_K_w_mean,
                          function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& f_K_w_square);

    template<typename scalar_type, typename dmn_type>
    void compute_mean_and_stddev(int nb_sample,
                                 function<scalar_type, dmn_type>& f_average,
                                 function<scalar_type, dmn_type>& f_square);

    template<typename scalar_type, typename dmn_type>
    void compute_mean_and_stddev(int nb_sample,
                                 function<std::complex<scalar_type>, dmn_type>& f_average,
                                 function<std::complex<scalar_type>, dmn_type>& f_square);

    template<typename w_dmn_t>
    void compute_mean_and_stddev(int nb_samples,
                                 function<double, w_dmn_t>& A_w,
                                 function<double, w_dmn_t>& A_w_stddev);

    template<typename w_dmn_t>
    void compute_mean_and_stddev(int nb_samples,
                                 function<double, dmn_3<b,s,w_dmn_t> >& A_nu_w,
                                 function<double, dmn_3<b,s,w_dmn_t> >& A_nu_w_stddev);

    template<typename k_dmn_t, typename w_dmn_t>
    void compute_mean_and_stddev(int nb_samples,
                                 function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& f_K_w_mean,
                                 function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& f_K_w_stddev);

    template<typename MOMS_real_type>
    void print_check_sums(MOMS_real_type& MOMS_real);

  private:

    parameters_type&  parameters;
    concurrency_type& concurrency;

    random_number_generator rng;

    //continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD> cpe_obj;
    continuous_pole_expansion<parameters_type, basis_function_t, k_DCA, w_REAL, WEIGHTED_GRADIENT_METHOD> cpe_obj;

    function<double, dmn_3<b, s, k_DCA> > G_k_beta_over_2;
    function<double, dmn_3<b, s, k_DCA> > integrated_A_k_div_cosh;
    function<double, dmn_3<b, s, k_DCA> > integrated_A_k_div_cosh_stddev;

    function<double, dmn_3<b, s, k_DCA> > G0_k_beta_over_2;
    function<double, dmn_3<b, s, k_DCA> > integrated_A0_k_div_cosh;
    function<double, dmn_3<b, s, k_DCA> > integrated_A0_k_div_cosh_stddev;

    function<             double , dmn_3<nu,nu,k_DCA>        > error_function;
    function<             double , dmn_3<nu,nu,k_DCA>        > error_function_stddev;

    function<std::complex<double>, dmn_4<nu,nu,k_DCA,w> > S_approx;
    function<std::complex<double>, dmn_4<nu,nu,k_DCA,w> > S_approx_stddev;

    function<std::complex<double>, dmn_4<nu,nu,k_DCA,w_IMAG> > f_orig;

    function<std::complex<double>, dmn_4<nu,nu,k_DCA,w_IMAG> > f_approx;
    function<std::complex<double>, dmn_4<nu,nu,k_DCA,w_IMAG> > f_approx_stddev;

    function<std::complex<double>, dmn_4<nu,nu,k_DCA,w_IMAG> > f_measured;
    function<std::complex<double>, dmn_4<nu,nu,k_DCA,w_IMAG> > f_measured_stddev;
  };

  template<class parameters_type, class basis_function_t>
  compute_spectrum<parameters_type, basis_function_t>::compute_spectrum(parameters_type&  parameters_ref):
    parameters(parameters_ref),
    concurrency(parameters.get_concurrency()),

    rng(concurrency.id(),
        concurrency.number_of_processors(),
        concurrency.get_seed()),

    cpe_obj(parameters, concurrency),

    G_k_beta_over_2               ("G_k_beta_over_2"),
    integrated_A_k_div_cosh       ("integrated_A_k_div_cosh"),
    integrated_A_k_div_cosh_stddev("integrated_A_k_div_cosh-stddev"),

    G0_k_beta_over_2               ("G0_k_beta_over_2"),
    integrated_A0_k_div_cosh       ("integrated_A0_k_div_cosh"),
    integrated_A0_k_div_cosh_stddev("integrated_A0_k_div_cosh-stddev"),

    error_function       ("L2-CPE-error"),
    error_function_stddev("L2-CPE-error-stddev"),

    S_approx       ("Sigma-approx"),
    S_approx_stddev("Sigma-approx-stddev"),

    f_orig("f-original"),

    f_approx       ("f-approx"),
    f_approx_stddev("f-approx-stddev"),

    f_measured       ("f-measured"),
    f_measured_stddev("f-measured-stddev")
  {}

  template<class parameters_type, class basis_function_t>
  compute_spectrum<parameters_type, basis_function_t>::~compute_spectrum()
  {}

  template<class parameters_type, class basis_function_t>
  template<typename MOMS_imag_type, typename MOMS_real_type>
  void compute_spectrum<parameters_type, basis_function_t>::write(std::string     file_name,
                                                                  MOMS_imag_type& MOMS_imag,
                                                                  MOMS_real_type& MOMS_real)
  {
    IO::FORMAT FORMAT = parameters.get_output_format();

    cout << "\n\n\t\t start writing " << file_name << "\n\n";

    switch(FORMAT)
      {
      case IO::JSON :
        {
          IO::writer<IO::JSON> writer;
          {
            writer.open_file(file_name);

            parameters.write(writer);
            MOMS_real .write(writer);
            cpe_obj   .write(writer);
            this->     write(writer);

            writer.close_file();
          }
        }
        break;

      case IO::HDF5 :
        {
          IO::writer<IO::HDF5> writer;
          {
            writer.open_file(file_name);

            parameters.write(writer);
            MOMS_real .write(writer);
            cpe_obj   .write(writer);
            this->     write(writer);

            writer.close_file();
          }
        }
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }
  }

  template<class parameters_type, class basis_function_t>
  template<IO::FORMAT DATA_FORMAT>
  void compute_spectrum<parameters_type, basis_function_t>::write(IO::writer<DATA_FORMAT>& writer)
  {
    writer.open_group("CPE-functions");

    writer.execute(error_function);
    writer.execute(error_function_stddev);

    writer.execute(S_approx);
    writer.execute(S_approx_stddev);

    writer.execute(f_orig);

    writer.execute(f_approx);
    writer.execute(f_approx_stddev);

    writer.execute(f_measured);
    writer.execute(f_measured_stddev);

    writer.close_group();

    writer.open_group("A-versus-G-functions");

    writer.execute(G_k_beta_over_2);
    writer.execute(integrated_A_k_div_cosh);

    writer.execute(G0_k_beta_over_2);
    writer.execute(integrated_A0_k_div_cosh);

    writer.close_group();
  }

  template<class parameters_type, class basis_function_t>
  template<typename MOMS_imag_type, typename MOMS_real_type>
  void compute_spectrum<parameters_type, basis_function_t>::execute(MOMS_imag_type& MOMS_imag,
                                                                    MOMS_real_type& MOMS_real)
  {
    if(parameters.simulate_gaussian_noise())
      execute_with_error_bars(MOMS_imag, MOMS_real);
    else
      execute_without_error_bars(MOMS_imag, MOMS_real);
  }

  template<class parameters_type, class basis_function_t>
  template<typename MOMS_imag_type, typename MOMS_real_type>
  void compute_spectrum<parameters_type, basis_function_t>::execute_without_error_bars(MOMS_imag_type& MOMS_imag,
                                                                                       MOMS_real_type& MOMS_real)
  {
    if(concurrency.id()==0)
      cout << "\n\n\t start analytic-continuation without error-bars (time = " << print_time() << ")\n";

    if(parameters.compute_free_spectrum())
      {
        MOMS_real.Sigma = 0;
        compute_G_k_w_on_lattice(MOMS_imag.H_HOST, MOMS_real.Sigma, MOMS_real.G0_k_w);

        compute_A_w(MOMS_real.A0_w, MOMS_real.A0_nu_w, MOMS_real.G0_k_w);

        SHOW::execute(MOMS_real.A0_w);
        SHOW::execute_on_bands(MOMS_real.A0_nu_w);

        compute_G_k_beta_over_2(G0_k_beta_over_2, MOMS_imag.G0_k_t);

        accumulate_integrated_A_k_div_cosh(integrated_A0_k_div_cosh,
                                           integrated_A0_k_div_cosh_stddev,
                                           MOMS_real.G0_k_w);
      }

    if(parameters.compute_cluster_spectrum() or
       parameters.compute_lattice_spectrum())
      perform_analytic_continuation(MOMS_imag.Sigma, MOMS_real.Sigma);

    if(parameters.compute_lattice_spectrum())
      compute_G_k_w_on_lattice(MOMS_imag.H_HOST, MOMS_real.Sigma, MOMS_real.G_k_w);

    if(parameters.compute_cluster_spectrum())
      compute_G_k_w_on_cluster(MOMS_real.G0_k_w, MOMS_real.Sigma, MOMS_real.G_k_w);

    if(parameters.compute_cluster_spectrum() or
       parameters.compute_lattice_spectrum())
      {
        {
          compute_G_k_beta_over_2(G_k_beta_over_2, MOMS_imag.G_k_t);

          accumulate_integrated_A_k_div_cosh(integrated_A_k_div_cosh, integrated_A_k_div_cosh_stddev, MOMS_real.G_k_w);

          accumulate_A_w(MOMS_real.A_w,    MOMS_real.A_w_stddev,
                         MOMS_real.A_nu_w, MOMS_real.A_nu_w_stddev,
                         MOMS_real.G_k_w);

          //      accumulate_f_K_w(MOMS_real.Sigma, MOMS_real.Sigma, MOMS_real.Sigma_stddev);
          //      accumulate_f_K_w(MOMS_real.G_k_w, MOMS_real.G_k_w, MOMS_real.G_k_w_stddev);
        }

        {
          SHOW::execute         (MOMS_real.A_w);
          SHOW::execute_on_bands(MOMS_real.A_nu_w);

          print_check_sums(MOMS_real);
        }

        {
          integrated_A_k_div_cosh_stddev = 0.;

          MOMS_real.A_w_stddev    = 0.;
          MOMS_real.A_nu_w_stddev = 0.;

          MOMS_real.Sigma_stddev = 0.;
          MOMS_real.G_k_w_stddev = 0.;

          f_approx_stddev   = 0.;
          f_measured_stddev = 0.;
        }

        if(concurrency.id()==concurrency.last())
          {
            cout.precision(6);
            cout<<scientific;

            cout << "\n\n\t beta/2*G(beta/2) versus 1./(2*T)\\int_{-inf}^{inf} dw A(w)/cosh(w/(2*T)) \n\n";

            for(int b_ind=0; b_ind<b::dmn_size(); b_ind++)
              for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++)
                cout << "\t" << G_k_beta_over_2               (0,0,k_ind)
                     << "\t" << integrated_A_k_div_cosh       (0,0,k_ind)
                     << "\t" << integrated_A_k_div_cosh_stddev(0,0,k_ind)
                     << "\n";
          }

      }
  }

  template<class parameters_type, class basis_function_t>
  template<typename MOMS_imag_type, typename MOMS_real_type>
  void compute_spectrum<parameters_type, basis_function_t>::execute_with_error_bars(MOMS_imag_type& MOMS_imag,
                                                                                    MOMS_real_type& MOMS_real)
  {
    if(parameters.compute_free_spectrum())
      {
        MOMS_real.Sigma = 0;
        compute_G_k_w_on_lattice(MOMS_imag.H_HOST, MOMS_real.Sigma, MOMS_real.G0_k_w);

        compute_A_w(MOMS_real.A0_w, MOMS_real.A0_nu_w, MOMS_real.G0_k_w);

        SHOW::execute         (MOMS_real.A0_w);
        SHOW::execute_on_bands(MOMS_real.A0_nu_w);
        //SHOW::execute_on_bands(MOMS_real.G0_k_w);

        compute_G_k_beta_over_2(G0_k_beta_over_2, MOMS_imag.G0_k_t);
      }

    {
      generate_f_original(MOMS_imag.Sigma);

      compute_G_k_beta_over_2(G_k_beta_over_2, MOMS_imag.G_k_t);
    }

    int nb_samples = parameters.get_nr_of_CPE_samples();
    {
      double magnitude = parameters.get_simulated_CPE_stddev();

      function<std::complex<double>, dmn_4<nu,nu,k_DCA,w     > > S_K_wi;
      function<std::complex<double>, dmn_4<nu,nu,k_DCA,w_REAL> > S_K_wr;

      function<std::complex<double>, dmn_4<nu,nu,k_DCA,w_REAL> > G_K_wr;

      if(concurrency.id()==0)
        cout << "\n\n";

      for(int l=0; l<nb_samples; l++)
        {
          if(concurrency.id()==0)
            cout << "\t start analytic-continuation on sample = " << l << " (time = " << print_time() << ")\n";

          {// generate a new sample that is equal for each MPI-task!

            if(concurrency.id()==0)
              {
                //MOMS_imag.Sigma_stddev = magnitude;

                // uniform error
                for(int l=0; l<MOMS_imag.Sigma_stddev.size(); l++)
                  {
                    std::complex<double> val = MOMS_imag.Sigma(l);

                    real(MOMS_imag.Sigma_stddev(l)) = magnitude*(real(val)*real(val))/(real(val)*real(val)+1.e-6);
                    imag(MOMS_imag.Sigma_stddev(l)) = magnitude*(imag(val)*imag(val))/(imag(val)*imag(val)+1.e-6);
                  }

                generate_new_sample(MOMS_imag.Sigma, MOMS_imag.Sigma_stddev, S_K_wi);
              }
            else
              {
                S_K_wi = 0;
              }

            concurrency.sum(S_K_wi);
          }

          {// do the analytic continuation ...
            if(parameters.compute_cluster_spectrum() or
               parameters.compute_lattice_spectrum())
              perform_analytic_continuation(S_K_wi, S_K_wr);

            if(parameters.compute_lattice_spectrum())
              compute_G_k_w_on_lattice(MOMS_imag.H_HOST, S_K_wr, G_K_wr);

            if(parameters.compute_cluster_spectrum())
              compute_G_k_w_on_cluster(MOMS_real.G0_k_w, S_K_wr, G_K_wr);
          }

          {// accumulate the function and the squared-function ...

            accumulate_integrated_A_k_div_cosh(integrated_A_k_div_cosh, integrated_A_k_div_cosh_stddev, G_K_wr);

            accumulate_A_w(MOMS_real.A_w,    MOMS_real.A_w_stddev,
                           MOMS_real.A_nu_w, MOMS_real.A_nu_w_stddev,
                           G_K_wr);

            accumulate_f_K_w(S_K_wr, MOMS_real.Sigma, MOMS_real.Sigma_stddev);
            accumulate_f_K_w(G_K_wr, MOMS_real.G_k_w, MOMS_real.G_k_w_stddev);
          }

          if(concurrency.id()==0)
            cout << "\n";
        }
    }

    {
      compute_mean_and_stddev(nb_samples, error_function, error_function_stddev);

      compute_mean_and_stddev(nb_samples, integrated_A_k_div_cosh, integrated_A_k_div_cosh_stddev);

      compute_mean_and_stddev(nb_samples, MOMS_real.A_w   , MOMS_real.A_w_stddev   );
      compute_mean_and_stddev(nb_samples, MOMS_real.A_nu_w, MOMS_real.A_nu_w_stddev);

      compute_mean_and_stddev(nb_samples, MOMS_real.Sigma, MOMS_real.Sigma_stddev);
      compute_mean_and_stddev(nb_samples, MOMS_real.G_k_w, MOMS_real.G_k_w_stddev);

      compute_mean_and_stddev(nb_samples, f_approx  , f_approx_stddev);
      compute_mean_and_stddev(nb_samples, f_measured, f_measured_stddev);
    }

    {
      SHOW::plot_error_bars(MOMS_real.A_w, MOMS_real.A_w_stddev);

      print_check_sums(MOMS_real);

      if(concurrency.id()==concurrency.last())
        {
          cout.precision(6);
          cout<<scientific;

          cout << "\n\n\t beta/2*G(beta/2) versus 1./(2*T)\\int_{-inf}^{inf} dw A(w)/cosh(w/(2*T)) \n\n";

          for(int b_ind=0; b_ind<b::dmn_size(); b_ind++)
            for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++)
              cout << "\t" << G_k_beta_over_2               (0,0,k_ind)
                   << "\t" << integrated_A_k_div_cosh       (0,0,k_ind)
                   << "\t" << integrated_A_k_div_cosh_stddev(0,0,k_ind)
                   << "\n";
        }

      //test_A_w_versus_G_t(MOMS_imag, MOMS_real);
    }
  }

  template<class parameters_type, class basis_function_t>
  void compute_spectrum<parameters_type, basis_function_t>::generate_f_original(function<std::complex<double>, dmn_4<nu,nu,k_DCA,w> >& S_K_w)
  {
    for(int w_0=0; w_0<w::dmn_size(); w_0++){
      for(int w_1=0; w_1<w_IMAG::dmn_size(); w_1++){

        if(abs(w::get_elements()[w_0]-w_IMAG::get_elements()[w_1])<1.e-6)
          {
            for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++)
              for(int nu_j=0; nu_j<nu::dmn_size(); nu_j++)
                for(int nu_i=0; nu_i<nu::dmn_size(); nu_i++)
                  f_orig(nu_i, nu_j, k_ind, w_1) = S_K_w(nu_i, nu_j, k_ind, w_0);
          }
      }
    }

    SHOW::execute(f_orig);
  }

  template<class parameters_type, class basis_function_t>
  void compute_spectrum<parameters_type, basis_function_t>::generate_new_sample(function<std::complex<double>, dmn_4<nu,nu,k_DCA,w> >& S_K_w_mean,
                                                                                function<std::complex<double>, dmn_4<nu,nu,k_DCA,w> >& S_K_w_stddev,
                                                                                function<std::complex<double>, dmn_4<nu,nu,k_DCA,w> >& S_K_w_sample)
  {
    for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){
      for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++){
        for(int b_j=0; b_j<b::dmn_size(); b_j++){
          for(int b_i=0; b_i<b::dmn_size(); b_i++){

	    if((b_i == b_j) and parameters.is_an_interacting_band(b_i))
	      {
		double error_re = rng.get_random_number();
		double error_im = rng.get_random_number();
		
		std::complex<double> mean   = S_K_w_mean  (b_i,0,b_j,0,k_ind,w_ind);
		std::complex<double> stddev = S_K_w_stddev(b_i,0,b_j,0,k_ind,w_ind);
		
		std::complex<double> error(gaussian_distribution(error_re, 0, real(stddev)),
					   gaussian_distribution(error_im, 0, imag(stddev)));
		
		assert(error_re==error_re);
		assert(error_im==error_im);
		
		S_K_w_sample(b_i,0,b_j,0,k_ind,w_ind) = mean+error;
		S_K_w_sample(b_i,1,b_j,1,k_ind,w_ind) = mean+error;
	      }
	    else
	      {
		S_K_w_sample(b_i,0,b_j,0,k_ind,w_ind) = 0.;
		S_K_w_sample(b_i,1,b_j,1,k_ind,w_ind) = 0.;
	      }
	  }
	}
      }
    }

    {
      symmetrize::execute(S_K_w_sample);
    }
  }


  template<class parameters_type, class basis_function_t>
  template<typename k_dmn_t, typename w_imag_dmn_t, typename w_real_dmn_t>
  void compute_spectrum<parameters_type, basis_function_t>::perform_analytic_continuation(function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_imag_dmn_t> >& S_k_w_imag,
                                                                                          function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_real_dmn_t> >& S_k_w_real)
  {
    if(concurrency.id()==0)
      cout << "\t\t start CPE (time = " << print_time() << ") --> ";

    //cpe_obj.execute_st(S_k_w_imag, S_k_w_real);
    {
      cpe_obj.execute_mt(S_k_w_imag, S_k_w_real);
      
      accumulate      (cpe_obj.get_error_function(), error_function, error_function_stddev);

      accumulate_f_K_w(cpe_obj.get_S_approx()      , S_approx      , S_approx_stddev);

      accumulate_f_K_w(cpe_obj.get_f_approx()      , f_approx      , f_approx_stddev);
      accumulate_f_K_w(cpe_obj.get_f_measured()    , f_measured    , f_measured_stddev);
    }

    if(concurrency.id()==0)
      cout << " (time = " << print_time() << ")\n";
  }

  template<class parameters_type, class basis_function_t>
  template<typename k_dmn_t, typename w_dmn_t>
  void compute_spectrum<parameters_type, basis_function_t>::compute_G_k_w_on_cluster(function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& G0_k_w,
                                                                                     function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& Sigma,
                                                                                     function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& G_k_w)
  {
    if(concurrency.id()==0)
      cout << "\t\t start AC on G_K_w (time = " << print_time() << ") --> ";

    LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> G_matrix ("G-matrix" , nu::dmn_size());
    LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> S_matrix ("S-matrix" , nu::dmn_size());
    LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> G0_matrix("G0-matrix", nu::dmn_size());

    LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> S_times_G0_matrix("SxG0_matrix", nu::dmn_size());


    LIN_ALG::GEINV<LIN_ALG::CPU>::plan<std::complex<double> > geinv_obj;

    geinv_obj.initialize(S_times_G0_matrix);

    for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); k_ind++){
      for(int w_ind=0; w_ind<w_dmn_t::dmn_size(); w_ind++){

        for(int j=0; j<nu::dmn_size(); j++)
          for(int i=0; i<nu::dmn_size(); i++)
            S_matrix(i,j) = Sigma(i,j,k_ind,w_ind);

        for(int j=0; j<nu::dmn_size(); j++)
          for(int i=0; i<nu::dmn_size(); i++)
            G0_matrix(i,j) = G0_k_w(i,j,k_ind,w_ind);

        LIN_ALG::GEMM<LIN_ALG::CPU>::execute(S_matrix, G0_matrix, S_times_G0_matrix);

        for(int i=0; i<nu::dmn_size(); i++)
          S_times_G0_matrix(i,i) = 1.-S_times_G0_matrix(i,i);

        //LIN_ALG::GEINV<LIN_ALG::CPU>::execute(S_times_G0_matrix);
        geinv_obj.execute(S_times_G0_matrix);

        LIN_ALG::GEMM<LIN_ALG::CPU>::execute(S_times_G0_matrix, G0_matrix, G_matrix);

        for(int j=0; j<nu::dmn_size(); j++)
          for(int i=0; i<nu::dmn_size(); i++)
            G_k_w(i,j,k_ind,w_ind) = G_matrix(i,j);
      }
    }

    if(concurrency.id()==0)
      cout << " (time = " << print_time() << ")\n";
  }

  template<class parameters_type, class basis_function_t>
  template<typename k_host_dmn_t, typename k_cluster_dmn_t, typename w_dmn_t>
  void compute_spectrum<parameters_type, basis_function_t>::compute_G_k_w_on_lattice(function<std::complex<double>, dmn_3<nu,nu,k_host_dmn_t           > >& H_k,
                                                                                     function<std::complex<double>, dmn_4<nu,nu,k_cluster_dmn_t,w_dmn_t> >& Sigma,
                                                                                     function<std::complex<double>, dmn_4<nu,nu,k_cluster_dmn_t,w_dmn_t> >& G_k_w)
  {
    if(concurrency.id()==0)
      cout << "\t\t start TIM (time = " << print_time() << ") --> ";

    DCA::coarsegraining_sp<parameters_type, k_DCA> coarsegraining_sp_obj(parameters);

    coarsegraining_sp_obj.compute_G_K_w_with_TIM(H_k, Sigma, G_k_w);

    if(concurrency.id()==0)
      cout << " (time = " << print_time() << ")\n";
  }

  /*
    template<class parameters_type, class basis_function_t>
    template<typename MOMS_real_type>
    void compute_spectrum<parameters_type, basis_function_t>::compute_G_k_w_on_cluster(MOMS_real_type& MOMS_real)
    {
    int matrix_size = b::dmn_size()*s::dmn_size()*b::dmn_size()*s::dmn_size();
    int matrix_dim  = b::dmn_size()*s::dmn_size();

    std::complex<double>* G_matrix  = new std::complex<double>[matrix_size];
    std::complex<double>* G0_matrix = new std::complex<double>[matrix_size];
    std::complex<double>* S_matrix   = new std::complex<double>[matrix_size];

    std::complex<double>* S_times_G0_matrix     = new std::complex<double>[matrix_size];
    std::complex<double>* S_times_G0_matrix_inv = new std::complex<double>[matrix_size];

    for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++){
    for(int w_ind=0; w_ind<w_REAL::dmn_size(); w_ind++){

    memcpy(G0_matrix, &MOMS_real.G0_k_w(0,0,k_ind,w_ind), sizeof(std::complex<double>)*matrix_size);
    memcpy(S_matrix , &MOMS_real.Sigma (0,0,k_ind,w_ind), sizeof(std::complex<double>)*matrix_size);

    {
    gemm_plan<std::complex<double> > gemm(matrix_dim);
    gemm.A = S_matrix;
    gemm.B = G0_matrix;
    gemm.C = S_times_G0_matrix;

    gemm.execute_plan();
    }

    for(int i=0; i<matrix_dim; i++)
    S_times_G0_matrix[i+i*matrix_dim] = 1.-S_times_G0_matrix[i+i*matrix_dim];

    {
    invert_plan<std::complex<double> > invert_pln(matrix_dim);
    memcpy(invert_pln.Matrix, S_times_G0_matrix, sizeof(std::complex<double>)*matrix_size);
    invert_pln.execute_plan();
    memcpy(S_times_G0_matrix_inv, invert_pln.inverted_matrix, sizeof(std::complex<double>)*matrix_size);
    }

    {
    gemm_plan<std::complex<double> > gemm(matrix_dim);
    gemm.A = S_times_G0_matrix_inv;
    gemm.B = G0_matrix;
    gemm.C = G_matrix;

    gemm.execute_plan();
    }

    memcpy(&MOMS_real.G_k_w(0,0,k_ind,w_ind), G_matrix, sizeof(std::complex<double>)*matrix_size);
    }
    }

    delete [] G_matrix;
    delete [] G0_matrix;
    delete [] S_matrix;

    delete [] S_times_G0_matrix;
    delete [] S_times_G0_matrix_inv;
    }
  */

  template<class parameters_type, class basis_function_t>
  template<typename k_dmn_t, typename w_dmn_t>
  void compute_spectrum<parameters_type, basis_function_t>::compute_A_w(function<             double ,                     w_dmn_t>  & A_w,
                                                                        function<             double , dmn_3<b,s,          w_dmn_t> >& A_nu_w,
                                                                        function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& G_k_w)
  {
    A_w    = 0;
    A_nu_w = 0;

    for(int b=0; b<b::dmn_size(); b++){
      for(int s=0; s<s::dmn_size(); s++){

        for(int w_ind=0; w_ind<w_REAL::dmn_size(); w_ind++){

          for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++){

            double term = -1./M_PI*imag(G_k_w (b,s,b,s,k_ind,w_ind))/double(k_DCA::dmn_size());

            A_w   (    w_ind) += term/s::dmn_size();
            A_nu_w(b,s,w_ind) += term;
          }
        }
      }
    }
  }

  template<class parameters_type, class basis_function_t>
  template<typename k_dmn_t>
  void compute_spectrum<parameters_type, basis_function_t>::compute_G_k_beta_over_2(function<double, dmn_3<b , s, k_dmn_t   > >& G_k_beta_over_2,
                                                                                    function<double, dmn_4<nu,nu, k_dmn_t, t> >& G_k_t)
  {
    int t_ind = 3.*t::dmn_size()/4.;
    assert(abs(t::get_elements()[t_ind]-parameters.get_beta()/2.)<1.e-6);

    for(int b_ind=0; b_ind<b::dmn_size(); b_ind++)
      for(int s_ind=0; s_ind<s::dmn_size(); s_ind++)
        for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); k_ind++)
          G_k_beta_over_2(b_ind,s_ind,k_ind) = -parameters.get_beta()*G_k_t(b_ind,s_ind,b_ind,s_ind,k_ind,t_ind);
  }

  template<class parameters_type, class basis_function_t>
  template<typename scalar_type, typename dmn_type>
  void compute_spectrum<parameters_type, basis_function_t>::accumulate(function<scalar_type, dmn_type>& f,
                                                                       function<scalar_type, dmn_type>& f_average,
                                                                       function<scalar_type, dmn_type>& f_square)
  {
    for(int i=0; i<f.size(); i++)
      {
        f_average(i) += f(i);
        f_square (i) += f(i)*f(i);
      }
  }

  template<class parameters_type, class basis_function_t>
  template<typename scalar_type, typename dmn_type>
  void compute_spectrum<parameters_type, basis_function_t>::accumulate(function<std::complex<scalar_type>, dmn_type>& f,
                                                                       function<std::complex<scalar_type>, dmn_type>& f_average,
                                                                       function<std::complex<scalar_type>, dmn_type>& f_square)
  {
    for(int i=0; i<f.size(); i++)
      {
        f_average(i) += f(i);

        real(f_square(i)) += real(f(i))*real(f(i));
        imag(f_square(i)) += imag(f(i))*imag(f(i));
      }
  }

  template<class parameters_type, class basis_function_t>
  template<typename k_dmn_t>
  void compute_spectrum<parameters_type, basis_function_t>::accumulate_integrated_A_k_div_cosh(function<             double , dmn_3<b , s , k_dmn_t        > >& int_A_k_div_cosh,
                                                                                               function<             double , dmn_3<b , s , k_dmn_t        > >& int_A_k_div_cosh_stddev,
                                                                                               function<std::complex<double>, dmn_4<nu, nu, k_dmn_t, w_REAL> >& G_k_w)
  {
    double twoT = 2./parameters.get_beta();

    double result=0;
    for(int b_ind=0; b_ind<b::dmn_size(); b_ind++){
      for(int s_ind=0; s_ind<s::dmn_size(); s_ind++){
        for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); k_ind++){

          result=0;
          for(int w_ind=0; w_ind<w_REAL::dmn_size(); w_ind++)
            result += -1./(twoT*M_PI)*imag(G_k_w(b_ind,s_ind,b_ind,s_ind,k_ind,w_ind))/(std::cosh(w_REAL::get_elements()[w_ind]/twoT))*basis_function_t::volume(w_ind);

          int_A_k_div_cosh       (b_ind,s_ind,k_ind) += result;
          int_A_k_div_cosh_stddev(b_ind,s_ind,k_ind) += result*result;

          //           result=0;
          //           for(int w_ind=0; w_ind<w_REAL::dmn_size(); w_ind++)
          //             result += -1./M_PI*imag(MOMS_real.G0_k_w(0,0,k_ind,w_ind))/(std::cosh(w_REAL::get_elements()[w_ind]/twoT))*basis_function_t::volume(w_ind)/twoT;

          //           integrated_A0_k_div_cosh(k_ind) = result;
        }
      }
    }
  }

  template<class parameters_type, class basis_function_t>
  template<typename k_dmn_t, typename w_dmn_t>
  void compute_spectrum<parameters_type, basis_function_t>::accumulate_A_w(function<double, w_dmn_t>&                                     A_w,
                                                                           function<double, w_dmn_t>&                                     A_w_stddev,
                                                                           function<double, dmn_3<b,s,w_dmn_t> >&                         A_nu_w,
                                                                           function<double, dmn_3<b,s,w_dmn_t> >&                         A_nu_w_stddev,
                                                                           function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& G_k_w)
  {
    for(int w_ind=0; w_ind<w_REAL::dmn_size(); w_ind++){

      double A_w_result = 0;

      for(int b=0; b<b::dmn_size(); b++){
        for(int s=0; s<s::dmn_size(); s++){

          double A_nu_w_result = 0;

          for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++){

            A_w_result    += -1./M_PI*imag(G_k_w (b,s,b,s,k_ind,w_ind))/double(k_DCA::dmn_size());
            A_nu_w_result += -1./M_PI*imag(G_k_w (b,s,b,s,k_ind,w_ind))/double(k_DCA::dmn_size());
          }

          A_nu_w       (b,s,w_ind) += A_nu_w_result;
          A_nu_w_stddev(b,s,w_ind) += A_nu_w_result*A_nu_w_result;
        }
      }

      A_w_result /= s::dmn_size();

      A_w       (w_ind) += A_w_result;
      A_w_stddev(w_ind) += A_w_result*A_w_result;
    }
  }

  template<class parameters_type, class basis_function_t>
  template<typename k_dmn_t, typename w_dmn_t>
  void compute_spectrum<parameters_type, basis_function_t>::accumulate_f_K_w(function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& f_K_w_sample,
                                                                             function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& f_K_w_mean,
                                                                             function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& f_K_w_square)
  {
    for(int i=0; i<f_K_w_sample.size(); i++)
      {
        f_K_w_mean(i) += f_K_w_sample(i);

        real(f_K_w_square(i)) += real(f_K_w_sample(i))*real(f_K_w_sample(i));
        imag(f_K_w_square(i)) += imag(f_K_w_sample(i))*imag(f_K_w_sample(i));
      }
  }

  template<class parameters_type, class basis_function_t>
  template<typename scalar_type, typename dmn_type>
  void compute_spectrum<parameters_type, basis_function_t>::compute_mean_and_stddev(int nb_samples,
                                                                                    function<scalar_type, dmn_type>& f_average,
                                                                                    function<scalar_type, dmn_type>& f_square)
  {
    f_average /= scalar_type(nb_samples);
    f_square  /= scalar_type(nb_samples);

    for(int i=0; i<f_average.size(); i++)
      f_square(i) = sqrt(abs(f_square(i)-f_average(i)*f_average(i)));
  }

  template<class parameters_type, class basis_function_t>
  template<typename scalar_type, typename dmn_type>
  void compute_spectrum<parameters_type, basis_function_t>::compute_mean_and_stddev(int nb_samples,
                                                                                    function<std::complex<scalar_type>, dmn_type>& f_average,
                                                                                    function<std::complex<scalar_type>, dmn_type>& f_square)
  {
    f_average /= scalar_type(nb_samples);
    f_square  /= scalar_type(nb_samples);

    for(int i=0; i<f_average.size(); i++)
      {
        real(f_square(i)) = sqrt(abs(real(f_square(i))-real(f_average(i))*real(f_average(i))));
        imag(f_square(i)) = sqrt(abs(imag(f_square(i))-imag(f_average(i))*imag(f_average(i))));
      }
  }


  template<class parameters_type, class basis_function_t>
  template<typename w_dmn_t>
  void compute_spectrum<parameters_type, basis_function_t>::compute_mean_and_stddev(int nb_samples,
                                                                                    function<double, w_dmn_t>& A_w,
                                                                                    function<double, w_dmn_t>& A_w_stddev)
  {
    A_w        /= double(nb_samples);
    A_w_stddev /= double(nb_samples);

    for(int i=0; i<A_w.size(); i++)
      A_w_stddev(i) = sqrt(abs(A_w_stddev(i)-A_w(i)*A_w(i)));
  }

  template<class parameters_type, class basis_function_t>
  template<typename w_dmn_t>
  void compute_spectrum<parameters_type, basis_function_t>::compute_mean_and_stddev(int nb_samples,
                                                                                    function<double, dmn_3<b,s,w_dmn_t> >& A_nu_w,
                                                                                    function<double, dmn_3<b,s,w_dmn_t> >& A_nu_w_stddev)
  {
    A_nu_w        /= double(nb_samples);
    A_nu_w_stddev /= double(nb_samples);

    for(int i=0; i<A_nu_w.size(); i++)
      A_nu_w_stddev(i) = sqrt(abs(A_nu_w_stddev(i)-A_nu_w(i)*A_nu_w(i)));
  }

  template<class parameters_type, class basis_function_t>
  template<typename k_dmn_t, typename w_dmn_t>
  void compute_spectrum<parameters_type, basis_function_t>::compute_mean_and_stddev(int nb_samples,
                                                                                    function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& f_K_w,
                                                                                    function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& f_K_w_stddev)
  {
    f_K_w        /= double(nb_samples);
    f_K_w_stddev /= double(nb_samples);

    for(int i=0; i<f_K_w.size(); i++)
      {
        real(f_K_w_stddev(i)) = sqrt(abs(real(f_K_w_stddev(i))-real(f_K_w(i))*real(f_K_w(i))));
        imag(f_K_w_stddev(i)) = sqrt(abs(imag(f_K_w_stddev(i))-imag(f_K_w(i))*imag(f_K_w(i))));
      }
  }

  template<class parameters_type, class basis_function_t>
  template<typename MOMS_real_type>
  void compute_spectrum<parameters_type, basis_function_t>::print_check_sums(MOMS_real_type& MOMS_real)
  {
    if(concurrency.id()==0)
      {
        double result=0;

        {
          cout << "\n\n\t integrated G0 and G : \n\n" ;
          for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++)
            {
              VECTOR_OPERATIONS::PRINT(k_DCA::get_elements()[k_ind]);

              result=0;
              for(int w_ind=0; w_ind<w_REAL::dmn_size()-1; w_ind++){
                double f_w      = imag(MOMS_real.G0_k_w(0,0,0,w_ind+1)+MOMS_real.G0_k_w(0,0,0,w_ind));
                double delta_w  = (w_REAL::get_elements()[w_ind+1]-w_REAL::get_elements()[w_ind]);

                result += -1./M_PI*delta_w*f_w/2;
              }
              cout << result << "\t";

              result=0;
              for(int w_ind=0; w_ind<w_REAL::dmn_size()-1; w_ind++){
                double f_w      = imag(MOMS_real.G_k_w(0,0,0,w_ind+1)+MOMS_real.G_k_w(0,0,0,w_ind));
                double delta_w  = (w_REAL::get_elements()[w_ind+1]-w_REAL::get_elements()[w_ind]);

                result += -1./M_PI*delta_w*f_w/2;
              }
              cout << result << "\n";
            }
        }

        cout << "\n\n";

        {
          cout << "integrated A0 and A: \n\n\t" ;

          result=0;
          for(int w_ind=0; w_ind<w_REAL::dmn_size()-1; w_ind++)
            result += (w_REAL::get_elements()[w_ind+1]-w_REAL::get_elements()[w_ind])*(MOMS_real.A0_w(w_ind)+MOMS_real.A0_w(w_ind+1))/2;

          cout << result << "\t";

          result=0;
          for(int w_ind=0; w_ind<w_REAL::dmn_size()-1; w_ind++)
            result += (w_REAL::get_elements()[w_ind+1]-w_REAL::get_elements()[w_ind])*(MOMS_real.A_w(w_ind)+MOMS_real.A_w(w_ind+1))/2;

          cout << result << "\n";

        }
      }
  }

  template<class parameters_type, class basis_function_t>
  template<typename MOMS_imag_type, typename MOMS_real_type>
  void compute_spectrum<parameters_type, basis_function_t>::test_A_w_versus_G_t(MOMS_imag_type& MOMS_imag,
                                                                                MOMS_real_type& MOMS_real)
  {
    double twoT = 2./parameters.get_beta();

    double result=0;
    for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++){

      result=0;
      for(int w_ind=0; w_ind<w_REAL::dmn_size(); w_ind++)
        result += -1./M_PI*imag(MOMS_real.G_k_w(0,0,k_ind,w_ind))/(std::cosh(w_REAL::get_elements()[w_ind]/twoT))*basis_function_t::volume(w_ind)/twoT;

      integrated_A_k_div_cosh(k_ind) = result;

      result=0;
      for(int w_ind=0; w_ind<w_REAL::dmn_size(); w_ind++)
        result += -1./M_PI*imag(MOMS_real.G0_k_w(0,0,k_ind,w_ind))/(std::cosh(w_REAL::get_elements()[w_ind]/twoT))*basis_function_t::volume(w_ind)/twoT;

      integrated_A0_k_div_cosh(k_ind) = result;
    }

    for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++){
      G_k_beta_over_2 (k_ind) = -parameters.get_beta()*MOMS_imag.G_k_t (0,0,k_ind,3*t::dmn_size()/4);
      G0_k_beta_over_2(k_ind) = -parameters.get_beta()*MOMS_imag.G0_k_t(0,0,k_ind,3*t::dmn_size()/4);
    }

    if(concurrency.id()==concurrency.last())
      {
        cout.precision(6);
        cout<<scientific;

        cout << "\n\n\t G0(beta/2) versus \\int_{-inf}^{inf} dw A0(w)/cosh(w) \n\n";

        for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++)
          cout << "\t" << G0_k_beta_over_2(k_ind)
               << "\t" << integrated_A0_k_div_cosh(k_ind)
               << "\n";

        cout << "\n\n\t G(beta/2) versus \\int_{-inf}^{inf} dw A(w)/cosh(w) \n\n";

        for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++)
          cout << "\t" << G_k_beta_over_2(k_ind)
               << "\t" << integrated_A_k_div_cosh(k_ind)
               << "\n";

        cout << "\n";
      }
  }

}

#endif




























/*
  template<class parameters_type, class basis_function_t>
  void compute_spectrum<parameters_type, basis_function_t>::compute_G_along_cut(function<std::complex<double>, nu_nu_k_LDA>& H_LDA)
  {
  int matrix_dim  = s::dmn_size()*b::dmn_size();
  int matrix_size = square(matrix_dim);

  function<std::complex<double>, dmn_3<nu, nu,  k_dmn_cut_type> > H_k("H_k_interpolated");

  wannier_interpolation<k_LDA, k_dmn_cut_type>::execute(H_LDA, H_k);

  function<int, k_dmn_cut_type> I_k     ("I_k");
  function<int, k_dmn_cut_type> I_k_next("I_k_next");
  {
  for(int k_ind=0; k_ind<k_dmn_cut_type::dmn_size(); k_ind++)
  I_k(k_ind) = k_DCA::parameter_type::get_K_index(k_dmn_cut_type::get_elements()[k_ind]);

  for(int k_ind=0; k_ind<k_dmn_cut_type::dmn_size(); k_ind++){
  for(int l=0; l<k_dmn_cut_type::dmn_size(); l++){

  int i = l%k_dmn_cut_type::dmn_size();

  if(I_k(k_ind) != I_k(i) ){
  I_k_next(k_ind) = i;
  break;
  }

  }
  }
  }

  invert_plan<std::complex<double> > invert_pln(matrix_dim);

  std::complex<double>* H_matrix     = new std::complex<double>[matrix_size];
  std::complex<double>* S_matrix     = new std::complex<double>[matrix_size];

  std::complex<double>* G_matrix     = new std::complex<double>[matrix_size];
  std::complex<double>* G_matrix_inv = new std::complex<double>[matrix_size];

  for(int w_ind=0; w_ind<w_REAL::dmn_size(); w_ind++){

  std::complex<double> z(w_REAL::get_elements()[w_ind]+parameters.get_chemical_potential(), parameters.get_EPSILON());

  for(int k_ind=0; k_ind<k_dmn_cut_type::dmn_size(); k_ind++){

  int K0_ind = I_k     (k_ind);//k_DCA::parameter_type::get_K_index(k_dmn_cut_type::get_elements()[k_ind]);
  //int K1_ind = I_k_next(k_ind);

  memcpy(H_matrix, &H_k  (0,0,0,0,k_ind)      , sizeof(std::complex<double>)*matrix_size);
  memcpy(S_matrix, &Sigma(0,0,0,0,K0_ind,w_ind), sizeof(std::complex<double>)*matrix_size);

  //       {
  //  double LL = std::sqrt(VECTOR_OPERATIONS::L2_NORM(k_DCA::get_elements()[I_k(k_ind)], k_DCA::get_elements()[I_k_next(k_ind)]));
  //  double l0 = std::sqrt(VECTOR_OPERATIONS::L2_NORM(k_DCA::get_elements()[I_k(k_ind)], k_dmn_cut_type::get_elements()[k_ind]));
  //  //double l1 = std::sqrt(VECTOR_OPERATIONS::L2_NORM(k_dmn_cut_type::get_elements()[k_ind], k_DCA::get_elements()[I_k_next(k_ind)]));

  //  for(int j=0; j<matrix_dim; j++)
  //    for(int i=0; i<matrix_dim; i++)
  //      S_matrix[i+matrix_dim*j] = ((LL-l0)/LL*Sigma(i,j,K0_ind,w_ind) + l0/LL*Sigma(i,j,K1_ind,w_ind));
  //       }

  for(int j=0; j<matrix_dim; j++)
  for(int i=0; i<matrix_dim; i++)
  G_matrix_inv[i+matrix_dim*j] = (i==j)? z : 0.;

  for(int j=0; j<matrix_dim; j++)
  for(int i=0; i<matrix_dim; i++)
  G_matrix_inv[i+matrix_dim*j] += -(H_matrix[i+matrix_dim*j] + S_matrix[i+matrix_dim*j]);

  {
  memcpy(invert_pln.Matrix, &G_matrix_inv[0]          , sizeof(std::complex<double>)*matrix_size);
  invert_pln.execute_plan();
  memcpy(&G_matrix[0]     , invert_pln.inverted_matrix, sizeof(std::complex<double>)*matrix_size);
  }

  memcpy(&G_k_dmn_t_w(0,0,0,0,k_ind,w_ind), &G_matrix[0], sizeof(std::complex<double>)*matrix_size);
  }
  }

  //   cout<<scientific;
  //   cout.precision(6);
  //   for(int k_ind=0; k_ind<k_dmn_cut_type::dmn_size(); k_ind++){
  //     for(int w_ind=0; w_ind<w_REAL::dmn_size(); w_ind++){
  //       cout << -imag(G_k_dmn_t_w(0,0,k_ind,w_ind)) << "\t";
  //     }
  //     cout << endl;
  //   }
  }
*/


//   template<class parameters_type, class basis_function_t>
//   void compute_spectrum<parameters_type, basis_function_t>::Fermi_surface(function<std::complex<double>, nu_nu_k_LDA  >& H_LDA,
//                                                                           function<std::complex<double>, nu_nu_k_DCA_w>& Sigma_wn)
//   {
//     typedef typename coarsegrain_type::nu_nu_mesh_k nu_nu_mesh_k;

//     CPE_obj.execute(Sigma_wn, Sigma);

//     int K_index = 0;
//     function<std::complex<double>, nu_nu_mesh_k>& interpolated_H = coarsegrain_obj.interpolate(H_LDA, K_index);

//   }




/*
  template<class parameters_type, class basis_function_t>
  void compute_spectrum<parameters_type, basis_function_t>::execute(function<std::complex<double>, nu_nu_k_LDA  >& H_LDA,
  function<std::complex<double>, nu_nu_k_DCA_w>& Sigma_wn)
  {
  double real_axis_off_set = parameters.get_real_frequencies_off_set();

  if(parameters.compute_free_spectrum())
  {
  concurrency << "\n\n\t start coarsegraining G0 \n\n";

  Sigma = 0;
  coarsegrain_obj.compute_G_from_H_and_Sigma(H_LDA, Sigma, G0_k_w, real_axis_off_set);//parameters.get_EPSILON());

  for(int b=0; b<b::dmn_size(); b++){
  for(int s=0; s<s::dmn_size(); s++){

  for(int w_ind=0; w_ind<w_REAL::dmn_size(); w_ind++){

  for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++){
  A0_nu_w (b,s,w_ind) += -1./M_PI*imag(G0_k_w (b,s,b,s,k_ind,w_ind));

  A0_w (w_ind) += -1./M_PI*imag(G0_k_w (b,s,b,s,k_ind,w_ind));
  }
  }
  }
  }

  A0_nu_w /= (double(k_DCA::dmn_size()));
  A0_w    /= (double(k_DCA::dmn_size())*b::dmn_size()*s::dmn_size());
  }

  {
  concurrency << "\n\n\t start perform analytic continuation \n\n";

  CPE_obj.execute(Sigma_wn, Sigma);

  if(parameters.compute_lattice_spectrum())
  {
  concurrency << "\n\n\t start coarsegraining G \n\n";
  coarsegrain_obj.compute_G_from_H_and_Sigma(H_LDA, Sigma, G_k_w, real_axis_off_set);//parameters.get_EPSILON());
  }

  if(parameters.compute_cluster_spectrum())
  {
  concurrency << "\n\n\t start computing G on the cluster \n\n";

  if(concurrency.id()==0)
  cout << "\n\t start reading : " << parameters.get_output_file_name() << endl;

  std::fstream input_file(&(parameters.get_output_file_name())[0], ios::in);
  G0_k_w.from_JSON(input_file, concurrency);
  input_file.close();

  compute_G_k_w_on_cluster();
  }

  concurrency << "\n\n\t compute spectra \n\n";

  for(int b=0; b<b::dmn_size(); b++){
  for(int s=0; s<s::dmn_size(); s++){

  for(int w_ind=0; w_ind<w_REAL::dmn_size(); w_ind++){

  for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++){
  A_nu_w (b,s,w_ind) += -1./M_PI*imag(G_k_w (b,s,b,s,k_ind,w_ind));

  A_w (w_ind) += -1./M_PI*imag(G_k_w (b,s,b,s,k_ind,w_ind));
  }
  }
  }
  }

  A_nu_w /= (double(k_DCA::dmn_size()));
  A_w    /= (double(k_DCA::dmn_size())*b::dmn_size()*s::dmn_size());

  test_A_w_versus_G_t();

  //compute_G_along_cut(H_LDA);
  }
  }
*/
