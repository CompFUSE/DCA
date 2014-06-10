//-*-C++-*-

#ifndef DCA_UPDATE_CHEMICAL_POTENTIAL_STEP_H
#define DCA_UPDATE_CHEMICAL_POTENTIAL_STEP_H

namespace DCA
{

  template<typename parameters_type, typename MOMS_type, typename coarsegraining_type>
  class update_chemical_potential
  {
#include "type_definitions.h"

    typedef typename parameters_type::concurrency_type concurrency_type;

  public:

    update_chemical_potential(parameters_type&     parameters_ref,
                              MOMS_type&           MOMS_ref,
                              coarsegraining_type& coarsegraining_ref);

    ~update_chemical_potential();

    void execute();

    double compute_density();

  private:

    double get_new_chemical_potential(double d_0,

                                      double mu_lb,
                                      double mu_ub,

                                      double n_lb,
                                      double n_ub);

    double compute_density_correction();

    void compute_density_coefficients(function<             double , dmn_2<nu, k_DCA> >& A,
                                      function<             double , dmn_2<nu, k_DCA> >& B,
                                      function<std::complex<double>, dmn_4<nu, nu, k_DCA, w> >& G);

    void search_bounds(double dens);

    void print_bounds();

  private:

    parameters_type&     parameters;
    concurrency_type&    concurrency;

    MOMS_type&           MOMS;
    coarsegraining_type& coarsegraining;

    std::pair<double, double> lower_bound;
    std::pair<double, double> upper_bound;
  };

  template<typename parameters_type, typename MOMS_type, typename coarsegraining_type>
  update_chemical_potential<parameters_type, MOMS_type, coarsegraining_type>::update_chemical_potential(parameters_type&     parameters_ref,
                                                                                                        MOMS_type&           MOMS_ref,
                                                                                                        coarsegraining_type& coarsegraining_ref):
    parameters(parameters_ref),
    concurrency(parameters.get_concurrency()),

    MOMS(MOMS_ref),
    coarsegraining(coarsegraining_ref)
  {}

  template<typename parameters_type, typename MOMS_type, typename coarsegraining_type>
  update_chemical_potential<parameters_type, MOMS_type, coarsegraining_type>::~update_chemical_potential()
  {}

  template<typename parameters_type, typename MOMS_type, typename coarsegraining_type>
  void update_chemical_potential<parameters_type, MOMS_type, coarsegraining_type>::execute()
  {
    double dens = compute_density();

    if(concurrency.id()==0)
      {
        cout.precision(6);
        cout<<scientific;

        cout << "\n\t\t initial chemical potential : " << parameters.get_chemical_potential() << " (" << dens << ")\n\n";
      }

    if(abs(dens-parameters.get_density())<1.e-3)
      return;

    search_bounds(dens);

    while(true)
      {
        //parameters.get_chemical_potential() = (lower_bound.first+upper_bound.first)/2.;

        double d_0 = parameters.get_density();

        double mu_lb = lower_bound.first;
        double mu_ub = upper_bound.first;

        double n_lb = lower_bound.second;
        double n_ub = upper_bound.second;

        /*
          parameters.get_chemical_potential() = (mu_ub-mu_lb)/(n_ub-n_lb)*(d_0-n_lb) + mu_lb;
        */

        parameters.get_chemical_potential() = get_new_chemical_potential(d_0, mu_lb, mu_ub, n_lb, n_ub);

        dens = compute_density();

        if(abs(dens-parameters.get_density())<1.e-3)
          {
            if(concurrency.id()==0)
              {
                cout.precision(6);
                cout<<scientific;

                cout << "\n\t\t final chemical potential : " << parameters.get_chemical_potential() << " (" << dens << ")\n";
              }

            break;
          }
        else
          {
            if(dens<parameters.get_density())
              {
                lower_bound.first  = parameters.get_chemical_potential();
                lower_bound.second = dens;
              }
            else
              {
                upper_bound.first  = parameters.get_chemical_potential();
                upper_bound.second = dens;
              }
          }

        print_bounds();
      }
  }

  template<typename parameters_type, typename MOMS_type, typename coarsegraining_type>
  double update_chemical_potential<parameters_type, MOMS_type, coarsegraining_type>::get_new_chemical_potential(double d_0,

                                                                                                                double mu_lb,
                                                                                                                double mu_ub,

                                                                                                                double n_lb,
                                                                                                                double n_ub)
  {
    double new_mu = (mu_ub-mu_lb)/(n_ub-n_lb)*(d_0-n_lb) + mu_lb;

    /*
      if(abs(new_mu-mu_lb)/abs(mu_ub-mu_lb) < 0.25 or
      abs(new_mu-mu_ub)/abs(mu_ub-mu_lb) < 0.25)
      {
      return (mu_lb+mu_ub)/2.;
      }
      else
      return new_mu;
    */

    return new_mu;
  }

  template<typename parameters_type, typename MOMS_type, typename coarsegraining_type>
  void update_chemical_potential<parameters_type, MOMS_type, coarsegraining_type>::search_bounds(double dens)
  {
    double factor = 2;
    double delta  = 0.1;

    if(dens<parameters.get_density())
      {
        lower_bound.first  = parameters.get_chemical_potential();
        lower_bound.second = dens;

        while(true)
          {
            parameters.get_chemical_potential() += delta;

            dens = compute_density();

            upper_bound.first  = parameters.get_chemical_potential();
            upper_bound.second = dens;

            print_bounds();

            if(parameters.get_density()<dens)
              break;
            else
              lower_bound = upper_bound;

            delta *= factor;
          }
      }
    else
      {
        upper_bound.first  = parameters.get_chemical_potential();
        upper_bound.second = dens;

        while(true)
          {
            parameters.get_chemical_potential() -= delta;

            dens = compute_density();

            lower_bound.first  = parameters.get_chemical_potential();
            lower_bound.second = dens;

            print_bounds();

            if(dens<parameters.get_density())
              break;
            else
              upper_bound = lower_bound;

            delta *= factor;
          }
      }
  }

  template<typename parameters_type, typename MOMS_type, typename coarsegraining_type>
  double update_chemical_potential<parameters_type, MOMS_type, coarsegraining_type>::compute_density()
  {
    if(parameters.use_interpolated_Self_energy())
      coarsegraining.compute_G_K_w(MOMS.H_HOST, MOMS.Sigma_lattice, MOMS.G_k_w);
    else
      coarsegraining.compute_G_K_w(MOMS.H_HOST, MOMS.Sigma        , MOMS.G_k_w);

    MOMS.G_k_w -= MOMS.G0_k_w;

    {
      MATH_ALGORITHMS::TRANSFORM<w, t>::execute(MOMS.G_k_w, MOMS.G_k_t);

      MOMS.G_k_t += MOMS.G0_k_t;

      MATH_ALGORITHMS::TRANSFORM<k_DCA, r_DCA>::execute(MOMS.G_k_t, MOMS.G_r_t);
    }

    MOMS.G_k_w += MOMS.G0_k_w;

    double result = compute_density_correction();

    for(int i=0; i<nu::dmn_size(); i++)
      result += 1.-MOMS.G_r_t(i,i,r_DCA::parameter_type::origin_index(),0);

    return result;
  }

  /*!
   *  We assume that G_ii(w>>0) ~ 1/(i w_m + A + B/(i w_m))
   */
  template<typename parameters_type, typename MOMS_type, typename coarsegraining_type>
  double update_chemical_potential<parameters_type, MOMS_type, coarsegraining_type>::compute_density_correction()
  {
    std::complex<double> I(0,1);

    double N_k  = k_DCA::dmn_size();
    double beta = parameters.get_beta();

    function<double, dmn_2<nu, k_DCA> > A;
    function<double, dmn_2<nu, k_DCA> > B;

    function<double, dmn_2<nu, k_DCA> > A0;
    function<double, dmn_2<nu, k_DCA> > B0;

    compute_density_coefficients(A , B , MOMS.G_k_w );
    compute_density_coefficients(A0, B0, MOMS.G0_k_w);

    double result = 0;

    for(int k_i=0; k_i<k_DCA::dmn_size(); k_i++){
      for(int nu_i=0; nu_i<nu::dmn_size(); nu_i++){

        double tmp = 0.0;
        double sum = 1.e-16;

        int l=w::dmn_size()/2;

        do
          {
            std::complex<double> I_wn = (M_PI/beta)*(1+2*l)*I;

            std::complex<double> G  = 1./(I_wn + A (nu_i, k_i) + B (nu_i, k_i)/I_wn);
            std::complex<double> G0 = 1./(I_wn + A0(nu_i, k_i) + B0(nu_i, k_i)/I_wn);

            tmp  = real(G - G0);
            sum += tmp;

            l += 1;
          }
        while( abs(tmp/sum) > 1.e-6 and l<1.e6);

        result += sum;
      }
    }

    result *= (2./(beta*N_k));

    return result;
  }

  /*
    template<typename parameters_type, typename MOMS_type, typename coarsegraining_type>
    double update_chemical_potential<parameters_type, MOMS_type, coarsegraining_type>::compute_density_correction()
    {
    double N_k  = k_DCA::dmn_size();
    double beta = parameters.get_beta();

    function<std::complex<double>, dmn_2<nu, k_DCA> > G_diag ("G_diag");
    function<std::complex<double>, dmn_2<nu, k_DCA> > G0_diag("G0_diag");

    for(int k_i=0; k_i<k_DCA::dmn_size(); k_i++)
    for(int nu_i=0; nu_i<nu::dmn_size(); nu_i++)
    G_diag(nu_i, k_i) = MOMS.G_k_w(nu_i, nu_i, k_i, w::dmn_size()-1);

    for(int k_i=0; k_i<k_DCA::dmn_size(); k_i++)
    for(int nu_i=0; nu_i<nu::dmn_size(); nu_i++)
    G0_diag(nu_i, k_i) = MOMS.G0_k_w(nu_i, nu_i, k_i, w::dmn_size()-1);

    std::complex<double> I(0,1);

    double result = 0;

    double wm = w::get_elements()[w::dmn_size()-1];

    for(int k_i=0; k_i<k_DCA::dmn_size(); k_i++){
    for(int nu_i=0; nu_i<nu::dmn_size(); nu_i++){

    double A =   real(1./G_diag(nu_i, k_i));
    double B = -(imag(1./G_diag(nu_i, k_i)) - wm)*wm;

    double A0 =   real(1./G0_diag(nu_i, k_i));
    double B0 = -(imag(1./G0_diag(nu_i, k_i)) - wm)*wm;

    assert(abs((G_diag (nu_i, k_i)) - 1./(I*wm + A  + B /(I*wm))) < 1.e-12);
    assert(abs((G0_diag(nu_i, k_i)) - 1./(I*wm + A0 + B0/(I*wm))) < 1.e-12);

    double tmp = 0.0;
    double sum = 0.0;

    int l=w::dmn_size()/2;

    do
    {
    std::complex<double> I_wn = (M_PI/beta)*(1+2*l)*I;

    tmp  = real(1./(I_wn + A + B/I_wn) - 1./(I_wn + A0 + B0/I_wn));
    sum += tmp;

    l += 1;
    }
    while( abs(tmp/sum) > 1.e-6);

    result += sum;
    }
    }

    result *= (2./(beta*N_k));

    return result;
    }
  */

  template<typename parameters_type, typename MOMS_type, typename coarsegraining_type>
  void update_chemical_potential<parameters_type, MOMS_type, coarsegraining_type>::compute_density_coefficients(function<             double , dmn_2<nu, k_DCA> >& A,
                                                                                                                function<             double , dmn_2<nu, k_DCA> >& B,
                                                                                                                function<std::complex<double>, dmn_4<nu, nu, k_DCA, w> >& G)
  {
    A = 0;
    B = 0;

    int nb_wm = parameters.get_number_of_tail_frequencies();

    if(nb_wm>0)
      {
        for(int k_i=0; k_i<k_DCA::dmn_size(); k_i++){
          for(int nu_i=0; nu_i<nu::dmn_size(); nu_i++){

            for(int w_i=w::dmn_size()-nb_wm; w_i<w::dmn_size(); w_i++){

              double wm = w::get_elements()[w_i];

              A(nu_i, k_i) +=        real(1./G(nu_i, nu_i, k_i, w_i));
              B(nu_i, k_i) += wm*(wm-imag(1./G(nu_i, nu_i, k_i, w_i)));
            }
          }
        }

        A /= nb_wm;
        B /= nb_wm;

        if(nb_wm==1)
          {
            std::complex<double> I(0,1);

            for(int k_i=0; k_i<k_DCA::dmn_size(); k_i++){
              for(int nu_i=0; nu_i<nu::dmn_size(); nu_i++){

                int    w_i = w::dmn_size()-1;
                double wm  = w::get_elements()[w_i];

                if(abs((G(nu_i, nu_i, k_i, w_i)) - 1./(I*wm + A(nu_i, k_i) + B(nu_i, k_i)/(I*wm))) > 1.e-12)
                  throw std::logic_error(__FUNCTION__);
              }
            }
          }
      }
    /*
      if(true)
      {
      function<std::complex<double>, dmn_4<nu, nu, k_DCA, w> > G_diff("G_diff");

      for(int k_i=0; k_i<k_DCA::dmn_size(); k_i++){
      for(int nu_i=0; nu_i<nu::dmn_size(); nu_i++){

      for(int w_i=0; w_i<w::dmn_size(); w_i++){

      std::complex<double> I_wm(0,w::get_elements()[w_i]);

      G_diff(nu_i, nu_i, k_i, w_i) = 1./(I_wm + A (nu_i, k_i) + B (nu_i, k_i)/I_wm) - G(nu_i, nu_i, k_i, w_i);
      }
      }
      }

      SHOW::execute_on_bands(G_diff);
      }
    */
  }

  template<typename parameters_type, typename MOMS_type, typename coarsegraining_type>
  void update_chemical_potential<parameters_type, MOMS_type, coarsegraining_type>::print_bounds()
  {
    if(concurrency.id()==0)
      {
        cout.precision(6);
        cout<<scientific;

        cout << "\t";
        cout << "\t mu : " << lower_bound.first << " (n = " << lower_bound.second << ")";
        cout << "\t mu : " << upper_bound.first << " (n = " << upper_bound.second << ")\t";
        cout << print_time() << "\n";
      }
  }

}

#endif
