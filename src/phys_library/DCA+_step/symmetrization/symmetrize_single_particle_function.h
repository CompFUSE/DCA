//-*-C++-*-

#ifndef SYMMETRIZE_SINGLE_PARTICLE_FUNCTION_H
#define SYMMETRIZE_SINGLE_PARTICLE_FUNCTION_H

/*! \class   symmetrize_single_particle_function
 *  \ingroup SYMMETRIZE-FUNCTIONS
 *
 *  \author Peter Staar
 *  \brief  This class symmetrize_single_particle_function.hs Greens-functions according to cluster-symmetries, matsubara frequencies and band-index symmetries.
 *
 *  \section tau imaginary-time domain
 *
 *   \f{eqnarray*}{
 *     G(\tau) &=& -G(\tau+\beta)
 *   \f}
 *
 *  \section omega matsubara-frequency domain
 *
 *   \f{eqnarray*}{
 *     G(\varpi) &=& \overline{G(-\varpi)}
 *   \f}
 *
 *  \section r_and_k cluster domain
 *
 *   For each symmetry operation \f$\mathcal{S}\f$ of the cluster-domain, we have
 *
 *   \f{eqnarray*}{
 *     G(\vec{r}) &=& G(\mathcal{S}(\vec{r})) \\
 *     G(\vec{k}) &=& G(\mathcal{S}(\vec{k})) \\
 *   \f}
 *
 */
class symmetrize_single_particle_function
{
#include "type_definitions.h"

  typedef dmn_0<electron_band_domain> b_dmn_t;

protected:

  template<typename scalartype, typename nu_dmn_t, typename f_dmn_0, typename f_dmn_1>
  static void execute(FUNC_LIB::function<scalartype, dmn_4<nu_dmn_t, nu_dmn_t, f_dmn_0, f_dmn_1> >& f,
                      FUNC_LIB::function<int       , dmn_2<nu_dmn_t, nu_dmn_t>                   >& H_symmetry,
                      bool do_diff=false);

  /*
    template<typename scalartype, typename cluster_type>
    static void execute(FUNC_LIB::function<scalartype, dmn_0<r_cluster<FULL, cluster_type> > >& f, bool do_diff=false);
  */

  template<typename scalartype, typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S>
  static void execute(FUNC_LIB::function<scalartype, dmn_0<cluster_domain<scalar_type, D, N, REAL_SPACE, S> > >& f, bool do_diff=false);

  /*
    template<typename scalartype, typename cluster_type>
    static void execute(FUNC_LIB::function<scalartype, dmn_0<k_cluster<FULL, cluster_type> > >& f, bool do_diff=false);
  */

  template<typename scalartype, typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S>
  static void execute(FUNC_LIB::function<scalartype, dmn_0<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S> > >& f, bool do_diff=false);

  template<typename scalartype, typename f_dmn_0, typename f_dmn_1>
  static void execute(FUNC_LIB::function<scalartype, dmn_4<b, b, f_dmn_0, f_dmn_1> >& f, bool do_diff=false);

  template<typename scalartype, typename f_dmn_0>
  static void execute(FUNC_LIB::function<scalartype, dmn_3<nu, nu, f_dmn_0> >& f, bool do_diff=false);

  template<typename scalartype, typename f_dmn_0, typename f_dmn_1>
  static void execute(FUNC_LIB::function<scalartype, dmn_4<nu, nu, f_dmn_0, f_dmn_1> >& f, bool do_diff=false);

  /*
    template<typename scalartype, typename cluster_type>
    static void execute(FUNC_LIB::function<scalartype, dmn_4<nu, nu, dmn_0<r_cluster<FULL, cluster_type> >, w> >& f, bool do_diff=false);
  */

private:

  template<typename scalartype>
  static void difference(scalartype val, std::string function_name, std::string dmn_name);

  template<typename scalartype>
  static void difference(scalartype val0, scalartype val1, std::string function_name, std::string dmn_name);

  template<typename scalartype, typename f_dmn_0, typename f_dmn_1>
  static void symmetrize_over_electron_spin(FUNC_LIB::function<scalartype, dmn_4<nu, nu, f_dmn_0, f_dmn_1> >& f, bool do_diff);

  template<typename scalartype>
  static scalartype conjugate(scalartype x);

  template<typename scalartype, typename dmn_t>
  static void execute(FUNC_LIB::function<scalartype, dmn_t>& f, bool do_diff=false);

  template<typename scalartype>
  static void execute(FUNC_LIB::function<scalartype, t>& f, bool do_diff=false);

  template<typename scalartype>
  static void execute(FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, t> >& f, bool do_diff=false);

  template<typename scalartype>
  static void execute(FUNC_LIB::function<scalartype, w>& f, bool do_diff=false);

  template<typename scalartype>
  static void execute(FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, w> >& f, bool do_diff=false);

  template<typename scalartype>
  static void execute(FUNC_LIB::function<scalartype, w_REAL>& f, bool do_diff=false);

  template<typename scalartype>
  static void execute(FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, w_REAL> >& f, bool do_diff=false);

  template<typename scalartype>
  static void execute(FUNC_LIB::function<scalartype, w_VERTEX>& f, bool do_diff=false);

  template<typename scalartype>
  static void execute(FUNC_LIB::function<scalartype, w_VERTEX_EXTENDED>& f, bool do_diff=false);

  /*
    template<typename scalartype, typename cluster_type>
    static void execute(FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, dmn_0<r_cluster<FULL, cluster_type> > > >& f, bool do_diff=false);
  */

  template<typename scalartype, typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S>
  static void execute(FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, dmn_0<cluster_domain<scalar_type, D, N, REAL_SPACE, S> > > >& f, bool do_diff=false);

  /*
    template<typename scalartype, typename cluster_type>
    static void execute(FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, dmn_0<k_cluster<FULL, cluster_type> > > >& f, bool do_diff=false);
  */

  template<typename scalartype, typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S>
  static void execute(FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, dmn_0<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S> > > >& f, bool do_diff=false);

};

template<typename scalartype, typename nu_dmn_t, typename f_dmn_0, typename f_dmn_1>
void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, dmn_4<nu_dmn_t, nu_dmn_t, f_dmn_0, f_dmn_1> >& f,
                                                  FUNC_LIB::function<int       , dmn_2<nu_dmn_t, nu_dmn_t>                   >& /*H_symmetry*/,
                                                  bool do_diff)
{
  execute(f, do_diff);
}

template<typename scalartype, typename f_dmn_0, typename f_dmn_1>
void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, dmn_4<b, b, f_dmn_0, f_dmn_1> >& f, bool do_diff)
{
  {
    FUNC_LIB::function<scalartype, f_dmn_0> f0(f.get_name());

    for(int nu_0=0; nu_0<b::dmn_size(); ++nu_0){
      for(int nu_1=0; nu_1<b::dmn_size(); ++nu_1){
        for(int ind_1=0; ind_1<f_dmn_1::dmn_size(); ++ind_1){

          for(int ind_0=0; ind_0<f_dmn_0::dmn_size(); ++ind_0)
            f0(ind_0) = f(nu_0, nu_1, ind_0, ind_1);

          symmetrize_single_particle_function::execute(f0, do_diff);

          for(int ind_0=0; ind_0<f_dmn_0::dmn_size(); ++ind_0)
            f(nu_0, nu_1, ind_0, ind_1) = f0(ind_0);
        }
      }
    }
  }

  {
    FUNC_LIB::function<scalartype, f_dmn_1> f1(f.get_name());

    for(int nu_0=0; nu_0<b::dmn_size(); ++nu_0){
      for(int nu_1=0; nu_1<b::dmn_size(); ++nu_1){
        for(int ind_0=0; ind_0<f_dmn_0::dmn_size(); ++ind_0){

          for(int ind_1=0; ind_1<f_dmn_1::dmn_size(); ++ind_1)
            f1(ind_1) = f(nu_0, nu_1, ind_0, ind_1);

          symmetrize_single_particle_function::execute(f1, do_diff);

          for(int ind_1=0; ind_1<f_dmn_1::dmn_size(); ++ind_1)
            f(nu_0, nu_1, ind_0, ind_1) = f1(ind_1);
        }
      }
    }
  }
}

template<typename scalartype, typename f_dmn_0>
void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, dmn_3<nu, nu, f_dmn_0> >& f, bool do_diff)
{
  FUNC_LIB::function<scalartype, dmn_3<b, b, f_dmn_0> > f0(f.get_name());

  for(int spin_ind=0; spin_ind<s::dmn_size(); ++spin_ind){

    for(int b_0=0; b_0<b::dmn_size(); ++b_0)
      for(int b_1=0; b_1<b::dmn_size(); ++b_1)
        for(int ind_0=0; ind_0<f_dmn_0::dmn_size(); ++ind_0)
          f0(b_0, b_1, ind_0) = f(b_0, spin_ind, b_1, spin_ind, ind_0);

    symmetrize_single_particle_function::execute(f0, do_diff);

    for(int b_0=0; b_0<b::dmn_size(); ++b_0)
      for(int b_1=0; b_1<b::dmn_size(); ++b_1)
        for(int ind_0=0; ind_0<f_dmn_0::dmn_size(); ++ind_0)
          f(b_0, spin_ind, b_1, spin_ind, ind_0) = f0(b_0, b_1, ind_0);
  }
}

template<typename scalartype, typename f_dmn_0, typename f_dmn_1>
void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, dmn_4<nu, nu, f_dmn_0, f_dmn_1> >& f, bool do_diff)
{
  //   if(do_diff)
  //     cout << "\tsymmetrizing : " << f.get_name() << endl;

  symmetrize_over_electron_spin(f, do_diff);

  {
    FUNC_LIB::function<scalartype, dmn_3<b, b, f_dmn_0> > f0(f.get_name());

    for(int ind_1=0; ind_1<f_dmn_1::dmn_size(); ++ind_1){

      for(int spin_ind=0; spin_ind<s::dmn_size(); ++spin_ind){

        for(int b_0=0; b_0<b::dmn_size(); ++b_0)
          for(int b_1=0; b_1<b::dmn_size(); ++b_1)
            for(int ind_0=0; ind_0<f_dmn_0::dmn_size(); ++ind_0)
              f0(b_0, b_1, ind_0) = f(b_0, spin_ind, b_1, spin_ind, ind_0, ind_1);

        symmetrize_single_particle_function::execute(f0, do_diff);

        for(int b_0=0; b_0<b::dmn_size(); ++b_0)
          for(int b_1=0; b_1<b::dmn_size(); ++b_1)
            for(int ind_0=0; ind_0<f_dmn_0::dmn_size(); ++ind_0)
              f(b_0, spin_ind, b_1, spin_ind, ind_0, ind_1) = f0(b_0, b_1, ind_0);
      }
    }
  }

  {
    FUNC_LIB::function<scalartype, dmn_3<b, b, f_dmn_1> > f1(f.get_name());

    for(int ind_0=0; ind_0<f_dmn_0::dmn_size(); ++ind_0){

      for(int spin_ind=0; spin_ind<s::dmn_size(); ++spin_ind){

        for(int ind_1=0; ind_1<f_dmn_1::dmn_size(); ++ind_1)
          for(int b_1=0; b_1<b::dmn_size(); ++b_1)
            for(int b_0=0; b_0<b::dmn_size(); ++b_0)
              f1(b_0, b_1, ind_1) = f(b_0, spin_ind, b_1, spin_ind, ind_0, ind_1);

        symmetrize_single_particle_function::execute(f1, do_diff);

        for(int ind_1=0; ind_1<f_dmn_1::dmn_size(); ++ind_1)
          for(int b_1=0; b_1<b::dmn_size(); ++b_1)
            for(int b_0=0; b_0<b::dmn_size(); ++b_0)
              f(b_0, spin_ind, b_1, spin_ind, ind_0, ind_1) = f1(b_0, b_1, ind_1);
      }
    }
  }
}

/*
  template<typename scalartype, typename cluster_type>
  void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, dmn_4<nu, nu, dmn_0<r_cluster<FULL, cluster_type> >, w> >& f, bool do_diff)
  {
  typedef dmn_0<r_cluster<FULL, cluster_type> > r_dmn_t;

  //   if(do_diff)
  //     cout << "\tsymmetrizing : " << f.get_name() << endl;

  symmetrize_over_electron_spin(f, do_diff);

  {
  FUNC_LIB::function<scalartype, dmn_3<b, b, r_dmn_t> > f0(f.get_name());

  for(int ind_1=0; ind_1<w::dmn_size(); ++ind_1){

  for(int spin_ind=0; spin_ind<s::dmn_size(); ++spin_ind){

  for(int b_0=0; b_0<b::dmn_size(); ++b_0)
  for(int b_1=0; b_1<b::dmn_size(); ++b_1)
  for(int ind_0=0; ind_0<r_dmn_t::dmn_size(); ++ind_0)
  f0(b_0, b_1, ind_0) = f(b_0, spin_ind, b_1, spin_ind, ind_0, ind_1);

  symmetrize_single_particle_function::execute(f0, do_diff);

  for(int b_0=0; b_0<b::dmn_size(); ++b_0)
  for(int b_1=0; b_1<b::dmn_size(); ++b_1)
  for(int ind_0=0; ind_0<r_dmn_t::dmn_size(); ++ind_0)
  f(b_0, spin_ind, b_1, spin_ind, ind_0, ind_1) = f0(b_0, b_1, ind_0);
  }
  }
  }

  {
  FUNC_LIB::function<scalartype, dmn_4<nu, nu, r_dmn_t, w> > f_new(f.get_name());

  int w_0 = w::dmn_size()-1;
  int r_0 = r_dmn_t::parameter_type::origin_index();

  for(int w_ind_0=0; w_ind_0<w::dmn_size()/2; ++w_ind_0){
  for(int r_ind_0=0; r_ind_0<r_dmn_t::dmn_size(); ++r_ind_0){

  int w_ind_1 = w_0-w_ind_0;
  int r_ind_1 = r_dmn_t::parameter_type::subtract(r_ind_0, r_0);

  for(int b0=0; b0<2*b_dmn_t::dmn_size(); ++b0){
  for(int b1=0; b1<2*b_dmn_t::dmn_size(); ++b1){

  scalartype tmp_0 = f(b0, b1, r_ind_0, w_ind_0);
  scalartype tmp_1 = f(b1, b0, r_ind_1, w_ind_1);

  scalartype tmp = (tmp_0+conjugate(tmp_1))/2.;

  f_new(b0, b1, r_ind_0, w_ind_0) =           tmp ;
  f_new(b1, b0, r_ind_1, w_ind_1) = conjugate(tmp);
  }
  }
  }
  }

  double max=0;
  for(int ind=0; ind<f.size(); ++ind)
  {
  max = std::max(max, abs(f(ind)-f_new(ind)));

  f(ind) = f_new(ind);
  }

  if(do_diff)
  difference(max, f.get_name(), "w-domain of the function : "+f.get_name()+"\n");
  }
  }
*/







template<typename scalartype, typename f_dmn_0, typename f_dmn_1>
void symmetrize_single_particle_function::symmetrize_over_electron_spin(FUNC_LIB::function<scalartype, dmn_4<nu, nu, f_dmn_0, f_dmn_1> >& f, bool /*do_diff*/)
{
  for(int ind_1=0; ind_1<f_dmn_1::dmn_size(); ind_1++){
    for(int ind_0=0; ind_0<f_dmn_0::dmn_size(); ind_0++){

      // spin-symmetry ... --> G_(e_UP, e_DN) == G_(e_DN, e_UP) == 0 !!
      for(int i=0; i<b::dmn_size(); i++){
        for(int j=0; j<b::dmn_size(); j++){
          f(i,0,j,1,ind_0,ind_1) = 0;
          f(i,1,j,0,ind_0,ind_1) = 0;

          scalartype tmp = (f(i,0,j,0,ind_0,ind_1) + f(i,1,j,1,ind_0,ind_1))/2.;

          f(i,0,j,0,ind_0,ind_1) = tmp;
          f(i,1,j,1,ind_0,ind_1) = tmp;
        }
      }
    }
  }
}







template<typename scalartype>
scalartype symmetrize_single_particle_function::conjugate(scalartype x)
{
  return std::conj(x);
}

template<>
double symmetrize_single_particle_function::conjugate(double x)
{
  return x;
}

template<typename scalartype>
void symmetrize_single_particle_function::difference(scalartype val, std::string function_name, std::string dmn_name)
{
  if(std::abs(val)>1.e-6){
    std::cout << "difference detected in : " << dmn_name << "\t" << function_name << "\t" << std::abs(val) << "\n\n";
    //throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

template<typename scalartype>
void symmetrize_single_particle_function::difference(scalartype val0, scalartype val1, std::string function_name, std::string dmn_name)
{
  if(abs(val0-val1)>1.e-6){
    std::cout << "difference detected in : " << dmn_name << "\t" << function_name << "\t" << abs(val0-val1) << "\n\n";
    //throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

template<>
void symmetrize_single_particle_function::difference(float val, std::string function_name, std::string dmn_name)
{
  if(std::abs(val)>1.e-3){
    std::cout << "difference detected in : " << dmn_name << "\t" << function_name << "\t" << std::abs(val) << "\n\n";
    //throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

template<>
void symmetrize_single_particle_function::difference(float val0, float val1, std::string function_name, std::string dmn_name)
{
  if(std::abs(val0-val1)>1.e-3){
    std::cout << "difference detected in : " << dmn_name << "\t" << function_name << "\t" << std::abs(val0-val1) << "\n\n";
    //throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

template<>
void symmetrize_single_particle_function::difference(std::complex<float> val, std::string function_name, std::string dmn_name)
{
  if(abs(val)>1.e-3){
    std::cout << "difference detected in : " << dmn_name << "\t" << function_name << "\t" << abs(val) << "\n\n";
    //throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

template<>
void symmetrize_single_particle_function::difference(std::complex<float> val0, std::complex<float> val1,
                                                     std::string function_name, std::string dmn_name)
{
  if(abs(val0-val1)>1.e-3){
    std::cout << "difference detected in : " << dmn_name << "\t" << function_name << "\t" << abs(val0-val1) << "\n\n";
    //throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

template<typename scalartype>
void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, t>& f, bool do_diff)
{
  int shift = t::dmn_size()/2;

  double max=0;
  for(int i=0; i<t::dmn_size()/2; i++)
  {
    max = std::max(max, abs((f(i) + f(i+shift))/2.));

    scalartype tmp = (f(i) - f(i+shift))/2.;

    f(i)       =  tmp ;
    f(i+shift) = -tmp;
  }

  if(do_diff)
    difference(max, f.get_name(), "tau-domain of the function : "+f.get_name()+"\n");
}

template<typename scalartype>
void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, t> >& f, bool do_diff)
{
  FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, t> > f_new;

  int t_0 = t::dmn_size()/2;

  for(int t_ind=0; t_ind<t::dmn_size()/2; t_ind++){
    for(int b0=0; b0<b_dmn_t::dmn_size(); ++b0){
      for(int b1=0; b1<b_dmn_t::dmn_size(); ++b1){

        scalartype tmp = (f(b0, b1, t_ind) - f(b1, b0, t_ind+t_0))/2.;

        f_new(b0, b1, t_ind)     =  tmp;
        f_new(b1, b0, t_ind+t_0) = -tmp;
      }
    }
  }



  double max=0;
  for(int ind=0; ind<f.size(); ++ind)
  {
    max = std::max(max, std::abs(f(ind)-f_new(ind)));

    f(ind) = f_new(ind);
  }

  if(do_diff)
    difference(max, f.get_name(), "t-domain of the function : "+f.get_name()+"\n");
}

template<typename scalartype>
void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, w>& f, bool do_diff)
{
  double max=0;
  for(int i=0; i<w::dmn_size()/2; i++)
  {
    max = std::max(max, abs( (f(i) - conjugate(f(w::dmn_size()-i-1)))/2. ));

    scalartype tmp = (f(i) + conjugate(f(w::dmn_size()-i-1)))/2.;

    f(i)                 =            tmp ;
    f(w::dmn_size()-1-i) =  conjugate(tmp);
  }

  if(do_diff)
    difference(max, f.get_name(), "w-domain of the function : "+f.get_name()+"\n");
}

template<typename scalartype>
void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, w> >& f, bool do_diff)
{
  FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, w> > f_new;

  int w_0 = w::dmn_size()-1;

  for(int w_ind=0; w_ind<w::dmn_size()/2; ++w_ind){
    for(int b0=0; b0<b_dmn_t::dmn_size(); ++b0){
      for(int b1=0; b1<b_dmn_t::dmn_size(); ++b1){

        scalartype tmp_0 = f(b0, b1,     w_ind);
        scalartype tmp_1 = f(b1, b0, w_0-w_ind);

        scalartype tmp = (tmp_0+conjugate(tmp_1))/2.;

        f_new(b0, b1,     w_ind) =           tmp;
        f_new(b1, b0, w_0-w_ind) = conjugate(tmp);
      }
    }
  }

  double max=0;
  for(int ind=0; ind<f.size(); ++ind)
  {
    max = std::max(max, abs(f(ind)-f_new(ind)));

    f(ind) = f_new(ind);
  }

  if(do_diff)
    difference(max, f.get_name(), "w-domain of the function : "+f.get_name()+"\n");
}

template<typename scalartype>
void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, w_REAL>& /*f*/, bool /*do_diff*/)
{}

template<typename scalartype>
void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, w_REAL> >& /*f*/, bool /*do_diff*/)
{}

template<typename scalartype>
void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, w_VERTEX>& f, bool do_diff)
{
  double max=0;
  for(int i=0; i<w_VERTEX::dmn_size()/2; i++)
  {
    max = std::max(max, abs( (f(i) - conjugate(f(w_VERTEX::dmn_size()-i-1)))/2. ));

    scalartype tmp = (f(i) + conjugate(f(w_VERTEX::dmn_size()-i-1)))/2.;

    f(i)                        =            tmp ;
    f(w_VERTEX::dmn_size()-i-1) =  conjugate(tmp);
  }

  if(do_diff)
    difference(max, "w_VERTEX-domain of the function : "+f.get_name()+"\n");
}

template<typename scalartype>
void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, w_VERTEX_EXTENDED>& f, bool do_diff)
{
  double max=0;
  for(int i=0; i<w_VERTEX_EXTENDED::dmn_size()/2; i++)
  {
    max = std::max(max, abs( (f(i) - conjugate(f(w_VERTEX_EXTENDED::dmn_size()-i-1)))/2. ));

    scalartype tmp = (f(i) + conjugate(f(w_VERTEX_EXTENDED::dmn_size()-i-1)))/2.;

    f(i)                                 =            tmp ;
    f(w_VERTEX_EXTENDED::dmn_size()-i-1) =  conjugate(tmp);
  }

  if(do_diff)
    difference(max, "w_VERTEX_EXTENDED-domain of the function : "+f.get_name()+"\n");
}

/*
  template<typename scalartype, typename cluster_type>
  void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, dmn_0<r_cluster<FULL, cluster_type> > >& f, bool do_diff)
  {
  typedef r_cluster<FULL, cluster_type> r_cluster_type;
  typedef dmn_0<r_cluster_type>         r_dmn_t;

  typedef typename r_cluster_type::sym_super_cell_dmn_t sym_super_cell_dmn_t;

  static FUNC_LIB::function<std::pair<int,int>, dmn_2< dmn_2<r_dmn_t,b_dmn_t>, sym_super_cell_dmn_t > >& r_symmetry_matrix = cluster_symmetry<r_cluster_type>::get_symmetry_matrix();

  static FUNC_LIB::function<scalartype, r_dmn_t> f_new;

  f_new = scalartype(0.);

  for(int S_ind=0; S_ind<sym_super_cell_dmn_t::dmn_size(); ++S_ind){

  for(int r_ind=0; r_ind<r_dmn_t::dmn_size(); ++r_ind){

  int R_new_ind = r_symmetry_matrix(r_ind, 0, S_ind).first;

  f_new(r_ind) += f(R_new_ind);
  }
  }

  if(sym_super_cell_dmn_t::dmn_size()>0)
  f_new /= double(sym_super_cell_dmn_t::dmn_size());
  else
  throw std::logic_error(__FUNCTION__);

  double max=0;
  for(int ind=0; ind<f.size(); ++ind)
  {
  max = std::max(max, abs(f(ind)-f_new(ind)));

  f(ind) = f_new(ind);
  }

  if(do_diff)
  difference(max, f.get_name(), "r-cluster-domain of the function : "+f.get_name()+"\n");
  }
*/

template<typename scalartype, typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S>
void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, dmn_0<cluster_domain<scalar_type, D, N, REAL_SPACE, S> > >& f, bool do_diff)
{
  typedef cluster_domain<scalar_type, D, N, REAL_SPACE, S> r_cluster_type;
  typedef dmn_0<r_cluster_type>                            r_dmn_t;

  typedef typename cluster_symmetry<r_cluster_type>::sym_super_cell_dmn_t sym_super_cell_dmn_t;

  static FUNC_LIB::function<std::pair<int,int>, dmn_2< dmn_2<r_dmn_t,b_dmn_t>, sym_super_cell_dmn_t > >& r_symmetry_matrix = cluster_symmetry<r_cluster_type>::get_symmetry_matrix();

  //typedef typename r_cluster_type::sym_super_cell_dmn_t sym_super_cell_dmn_t;
  //static FUNC_LIB::function<std::pair<int,int>, dmn_2< dmn_2<r_dmn_t,b_dmn_t>, sym_super_cell_dmn_t > >& r_symmetry_matrix = cluster_symmetry<r_cluster_type>::get_symmetry_matrix();

  static FUNC_LIB::function<scalartype, r_dmn_t> f_new;

  f_new = scalartype(0.);

  for(int S_ind=0; S_ind<sym_super_cell_dmn_t::dmn_size(); ++S_ind){

    for(int r_ind=0; r_ind<r_dmn_t::dmn_size(); ++r_ind){

      int R_new_ind = r_symmetry_matrix(r_ind, 0, S_ind).first;

      f_new(r_ind) += f(R_new_ind);
    }
  }

  if(sym_super_cell_dmn_t::dmn_size()>0)
    f_new /= double(sym_super_cell_dmn_t::dmn_size());
  else
    throw std::logic_error(__FUNCTION__);

  double max=0;
  for(int ind=0; ind<f.size(); ++ind)
  {
    max = std::max(max, std::abs(f(ind)-f_new(ind)));

    f(ind) = f_new(ind);
  }

  if(do_diff)
    difference(max, f.get_name(), "r-cluster-domain of the function : "+f.get_name()+"\n");
}

/*
  template<typename scalartype, typename cluster_type>
  void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, dmn_0<r_cluster<FULL, cluster_type> > > >& f, bool do_diff)
  {
  typedef r_cluster<FULL, cluster_type> r_cluster_type;
  typedef dmn_0<r_cluster_type>         r_dmn_t;

  typedef typename cluster_symmetry<r_cluster_type>::sym_super_cell_dmn_t sym_super_cell_dmn_t;

  static FUNC_LIB::function<std::pair<int,int>, dmn_2< dmn_2<r_dmn_t,b_dmn_t>, sym_super_cell_dmn_t > >& r_symmetry_matrix = cluster_symmetry<r_cluster_type>::get_symmetry_matrix();

  static FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, r_dmn_t> > f_new;

  f_new = scalartype(0.);

  for(int S_ind=0; S_ind<sym_super_cell_dmn_t::dmn_size(); ++S_ind){

  for(int b0=0; b0<b_dmn_t::dmn_size(); ++b0){
  for(int b1=0; b1<b_dmn_t::dmn_size(); ++b1){

  for(int r_ind=0; r_ind<r_dmn_t::dmn_size(); ++r_ind){

  int R_new_ind = r_symmetry_matrix(r_ind, 0, S_ind).first;

  int b0_new = r_symmetry_matrix(0    , b0, S_ind).second;
  int b1_new = r_symmetry_matrix(r_ind, b1, S_ind).second;

  f_new(b0, b1, r_ind) += f(b0_new , b1_new, R_new_ind);
  }
  }
  }
  }

  if(sym_super_cell_dmn_t::dmn_size()>0)
  f_new /= double(sym_super_cell_dmn_t::dmn_size());
  else
  throw std::logic_error(__FUNCTION__);

  double max=0;
  for(int ind=0; ind<f.size(); ++ind)
  {
  max = std::max(max, abs(f(ind)-f_new(ind)));

  f(ind) = f_new(ind);
  }

  if(do_diff)
  difference(max, f.get_name(), "r-cluster-domain of the function : "+f.get_name()+"\n");
  }
*/

template<typename scalartype, typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S>
void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, dmn_0<cluster_domain<scalar_type, D, N, REAL_SPACE, S> > > >& f,
                                                  bool do_diff)
{
  typedef cluster_domain<scalar_type, D, N, REAL_SPACE, S> r_cluster_type;
  typedef dmn_0<r_cluster_type>                            r_dmn_t;

  typedef typename cluster_symmetry<r_cluster_type>::sym_super_cell_dmn_t sym_super_cell_dmn_t;

  static FUNC_LIB::function<std::pair<int,int>, dmn_2< dmn_2<r_dmn_t,b_dmn_t>, sym_super_cell_dmn_t > >& r_symmetry_matrix = cluster_symmetry<r_cluster_type>::get_symmetry_matrix();

  static FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, r_dmn_t> > f_new;

  f_new = scalartype(0.);

  for(int S_ind=0; S_ind<sym_super_cell_dmn_t::dmn_size(); ++S_ind){

    for(int b0=0; b0<b_dmn_t::dmn_size(); ++b0){
      for(int b1=0; b1<b_dmn_t::dmn_size(); ++b1){

        for(int r_ind=0; r_ind<r_dmn_t::dmn_size(); ++r_ind){

          int R_new_ind = r_symmetry_matrix(r_ind, 0, S_ind).first;

          int b0_new = r_symmetry_matrix(0    , b0, S_ind).second;
          int b1_new = r_symmetry_matrix(r_ind, b1, S_ind).second;

          f_new(b0, b1, r_ind) += f(b0_new , b1_new, R_new_ind);
        }
      }
    }
  }

  if(sym_super_cell_dmn_t::dmn_size()>0)
    f_new /= double(sym_super_cell_dmn_t::dmn_size());
  else
    throw std::logic_error(__FUNCTION__);

  double max=0;
  for(int ind=0; ind<f.size(); ++ind)
  {
    max = std::max(max, std::abs(f(ind)-f_new(ind)));

    f(ind) = f_new(ind);
  }

  if(do_diff)
    difference(max, f.get_name(), "r-cluster-domain of the function : "+f.get_name()+"\n");
}

/*
  template<typename scalartype, typename cluster_type>
  void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, dmn_0<k_cluster<FULL, cluster_type> > >& f, bool do_diff)
  {
  typedef k_cluster<FULL, cluster_type> k_cluster_type;
  typedef dmn_0<k_cluster_type>         k_dmn_t;

  typedef typename cluster_symmetry<k_cluster_type>::sym_super_cell_dmn_t sym_super_cell_dmn_t;

  static FUNC_LIB::function<std::pair<int,int>, dmn_2< dmn_2<k_dmn_t,b_dmn_t>, sym_super_cell_dmn_t > >& k_symmetry_matrix = cluster_symmetry<k_cluster_type>::get_symmetry_matrix();

  static FUNC_LIB::function<scalartype, k_dmn_t> f_new;

  f_new = scalartype(0.);

  for(int S_ind=0; S_ind<sym_super_cell_dmn_t::dmn_size(); ++S_ind){

  for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); ++k_ind){

  int K_new_ind = k_symmetry_matrix(k_ind, 0, S_ind).first;

  f_new(k_ind) += f(K_new_ind);
  }
  }

  if(sym_super_cell_dmn_t::dmn_size()>0)
  f_new /= double(sym_super_cell_dmn_t::dmn_size());
  else
  throw std::logic_error(__FUNCTION__);

  double max=0;
  for(int ind=0; ind<f.size(); ++ind)
  {
  max = std::max(max, abs(f(ind)-f_new(ind)));

  f(ind) = f_new(ind);
  }

  if(do_diff)
  difference(max, f.get_name(), "k-cluster-domain of the function : "+f.get_name()+"\n");
  }
*/

template<typename scalartype, typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S>
void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, dmn_0<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S> > >& f, bool do_diff)
{
  typedef cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S> k_cluster_type;
  typedef dmn_0<k_cluster_type>                                k_dmn_t;

  typedef typename cluster_symmetry<k_cluster_type>::sym_super_cell_dmn_t sym_super_cell_dmn_t;

  static FUNC_LIB::function<std::pair<int,int>, dmn_2< dmn_2<k_dmn_t,b_dmn_t>, sym_super_cell_dmn_t > >& k_symmetry_matrix = cluster_symmetry<k_cluster_type>::get_symmetry_matrix();

  static FUNC_LIB::function<scalartype, k_dmn_t> f_new;

  f_new = scalartype(0.);

  for(int S_ind=0; S_ind<sym_super_cell_dmn_t::dmn_size(); ++S_ind){

    for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); ++k_ind){

      int K_new_ind = k_symmetry_matrix(k_ind, 0, S_ind).first;

      f_new(k_ind) += f(K_new_ind);
    }
  }

  if(sym_super_cell_dmn_t::dmn_size()>0)
    f_new /= double(sym_super_cell_dmn_t::dmn_size());
  else
    throw std::logic_error(__FUNCTION__);

  double max=0;
  for(int ind=0; ind<f.size(); ++ind)
  {
    max = std::max(max, abs(f(ind)-f_new(ind)));

    f(ind) = f_new(ind);
  }

  if(do_diff)
    difference(max, f.get_name(), "k-cluster-domain of the function : "+f.get_name()+"\n");
}

/*
  template<typename scalartype, typename cluster_type>
  void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, dmn_0<k_cluster<FULL, cluster_type> > > >& f, bool do_diff)
  {
  typedef k_cluster<FULL, cluster_type> k_cluster_type;
  typedef dmn_0<k_cluster_type>         k_dmn_t;

  typedef typename cluster_symmetry<k_cluster_type>::sym_super_cell_dmn_t sym_super_cell_dmn_t;

  static FUNC_LIB::function<std::pair<int,int>, dmn_2< dmn_2<k_dmn_t,b_dmn_t>, sym_super_cell_dmn_t > >& k_symmetry_matrix = cluster_symmetry<k_cluster_type>::get_symmetry_matrix();

  static FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, k_dmn_t> > f_new;

  f_new = scalartype(0.);

  for(int S_ind=0; S_ind<sym_super_cell_dmn_t::dmn_size(); ++S_ind){

  for(int b0=0; b0<b_dmn_t::dmn_size(); ++b0){
  for(int b1=0; b1<b_dmn_t::dmn_size(); ++b1){

  for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); ++k_ind){

  int K_new_ind = k_symmetry_matrix(k_ind, b0, S_ind).first;

  int b0_new = k_symmetry_matrix(0    , b0, S_ind).second;
  int b1_new = k_symmetry_matrix(k_ind, b1, S_ind).second;

  f_new(b0, b1, k_ind) += f(b0_new , b1_new, K_new_ind);
  }
  }
  }
  }

  f_new /= double(sym_super_cell_dmn_t::dmn_size());

  double max=0;
  for(int ind=0; ind<f.size(); ++ind)
  {
  max = std::max(max, abs(f(ind)-f_new(ind)));

  f(ind) = f_new(ind);
  }

  if(do_diff)
  difference(max, f.get_name(), "k-clusterdomain of the function : "+f.get_name()+"\n");
  }
*/

template<typename scalartype, typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S>
void symmetrize_single_particle_function::execute(FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, dmn_0<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S> > > >& f,
                                                  bool do_diff)
{
  typedef cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S> k_cluster_type;
  typedef dmn_0<k_cluster_type>                                k_dmn_t;

  typedef typename cluster_symmetry<k_cluster_type>::sym_super_cell_dmn_t sym_super_cell_dmn_t;

  static FUNC_LIB::function<std::pair<int,int>, dmn_2< dmn_2<k_dmn_t,b_dmn_t>, sym_super_cell_dmn_t > >& k_symmetry_matrix = cluster_symmetry<k_cluster_type>::get_symmetry_matrix();

  static FUNC_LIB::function<scalartype, dmn_3<b_dmn_t, b_dmn_t, k_dmn_t> > f_new;

  f_new = scalartype(0.);

  for(int S_ind=0; S_ind<sym_super_cell_dmn_t::dmn_size(); ++S_ind){

    for(int b0=0; b0<b_dmn_t::dmn_size(); ++b0){
      for(int b1=0; b1<b_dmn_t::dmn_size(); ++b1){

        for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); ++k_ind){

          int K_new_ind = k_symmetry_matrix(k_ind, b0, S_ind).first;  // FIXME: b0 -> b1

          int b0_new = k_symmetry_matrix(0    , b0, S_ind).second;
          int b1_new = k_symmetry_matrix(k_ind, b1, S_ind).second;

          f_new(b0, b1, k_ind) += f(b0_new , b1_new, K_new_ind);
        }
      }
    }
  }

  f_new /= double(sym_super_cell_dmn_t::dmn_size());

  double max=0;
  for(int ind=0; ind<f.size(); ++ind)
  {
    max = std::max(max, std::abs(f(ind)-f_new(ind)));

    f(ind) = f_new(ind);
  }

  if(do_diff)
    difference(max, f.get_name(), "k-clusterdomain of the function : "+f.get_name()+"\n");
}

#endif
