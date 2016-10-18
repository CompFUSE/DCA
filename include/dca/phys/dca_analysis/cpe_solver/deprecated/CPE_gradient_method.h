//-*-C++-*-

#ifndef CONTINUOUS_POLE_EXPANSION_GRADIENT_METHOD_H
#define CONTINUOUS_POLE_EXPANSION_GRADIENT_METHOD_H

namespace DCA
{
  /*!
   *  \ingroup CPE
   *
   *  \author  Peter Staar
   *  \brief   This class implements a CPE analytic continuation, using a gradient method as a minimization-algorithm.
   *  \version 1.0
   */
  template<class parameters_type, class basis_function_t>
  class continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>
  {
#include "type_definitions.h"

  public:

    typedef double                                     scalartype;

    typedef typename parameters_type::profiler_type    profiler_type;
    typedef typename parameters_type::concurrency_type concurrency_type;

    typedef dmn_0<basis_function_t>        alpha_dmn_t;
    typedef dmn_4<nu,nu,k_DCA,alpha_dmn_t> nu_nu_k_DCA_alpha_dmn;


  public:

    continuous_pole_expansion(parameters_type&  parameters,
                              concurrency_type& concurrency,
                              bool              fixed_zero_moment=false,
                              double            zero_moment=0,
                              bool              fixed_first_moment=false,
                              double            first_moment=1);

    ~continuous_pole_expansion();

    template<typename target_dmn_t>
    void execute(FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, w> >&            f_source,
                 FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, target_dmn_t> >& f_target,
                 bool                                                                 use_previous_result=false);

    void write(std::string filename);

    template<IO::FORMAT DATA_FORMAT>
    void write(IO::writer<DATA_FORMAT>& reader);

    template<typename stream_type>
    void domains_to_JSON(stream_type& ss);

    template<typename stream_type>
    void to_JSON(stream_type& ss);

  private:

    void initialize();

    void compute_singular_values();

    bool spline_input_values_to_wn(std::complex<double>* input_values);

    void perform_continuous_pole_expansion();

    void f_wn_distribute(int dmn_number, int* coordinate);

    void compute_gradient_alphas(bool normalize=true);
    void compute_gradient_Sigma_0(bool normalize=true);

    void find_new_Sigma_0(double LAMBDA_MAX=1.);

    int  find_new_alpha  (double LAMBDA_MAX=1.);
    int  find_new_alpha_2(double LAMBDA_MAX=1.);

    double evaluate_Lambda_norm(FUNC_LIB::function<scalartype, alpha_dmn_t>& alpha_real);

    void smooth_alpha();

    void project_from_real_axis_to_imaginary_axis(scalartype* real_values, scalartype* imag_values, scalartype* matrix);
    void project_from_imaginary_axis_to_real_axis(scalartype* imag_values, scalartype* real_values, scalartype* matrix);

    template<typename target_dmn_t>
    void spline_alpha_to_output_values(target_dmn_t& target_dmn, std::complex<double>* output_values);

    void spline_alpha_to_output_values_w_IMAG(std::complex<double>* output_values);
    void spline_alpha_to_output_values_w_REAL(std::complex<double>* output_values);

  private:

    parameters_type&  parameters;
    concurrency_type& concurrency;

    bool static_zero_moment;
    bool project_to_unit_plane;

    scalartype      Sigma_0;
    scalartype grad_Sigma_0;

    FUNC_LIB::function<std::complex<scalartype>, nu_nu_k_DCA>        Sigma_0_moment;

    FUNC_LIB::function<scalartype              , nu_nu_k_DCA_alpha_dmn> alpha;

    FUNC_LIB::function<std::complex<scalartype>, nu_nu_k_DCA_w_IMAG> f_approx;
    FUNC_LIB::function<std::complex<scalartype>, nu_nu_k_DCA_w_IMAG> f_measured;

    FUNC_LIB::function<scalartype, alpha_dmn_t> alpha_ptr;

    FUNC_LIB::function<std::complex<scalartype>, w_IMAG> f_wn;
    FUNC_LIB::function<std::complex<scalartype>, w_IMAG> F_wn;

    FUNC_LIB::function<scalartype, w_IMAG> f_wn_re_ptr;
    FUNC_LIB::function<scalartype, w_IMAG> f_wn_im_ptr;
    FUNC_LIB::function<scalartype, w_IMAG> F_wn_re_ptr;
    FUNC_LIB::function<scalartype, w_IMAG> F_wn_im_ptr;

    FUNC_LIB::function<scalartype, alpha_dmn_t> gradient;
    FUNC_LIB::function<scalartype, alpha_dmn_t> gradient_1;
    FUNC_LIB::function<scalartype, alpha_dmn_t> gradient_2;

    FUNC_LIB::function<scalartype, alpha_dmn_t> a_new    ;
    FUNC_LIB::function<scalartype, alpha_dmn_t> a_new_cpy;

    FUNC_LIB::function<std::complex<scalartype>, dmn_2<w_IMAG, alpha_dmn_t> > A_mat;

    scalartype* A1;
    scalartype* A2;

    int lower_bound, upper_bound, bounds_difference;
  };

  template<class parameters_type, class basis_function_t>
  continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::continuous_pole_expansion(parameters_type&  parameters_ref,
                                                                                                           concurrency_type& concurrency_ref,
                                                                                                           bool              fixed_zeroth_moment,
                                                                                                           double            zeroth_moment_val,
                                                                                                           bool              fixed_first_moment,
                                                                                                           double            first_moment_val):
    parameters(parameters_ref),
    concurrency(concurrency_ref),

    static_zero_moment(fixed_zeroth_moment),
    project_to_unit_plane(fixed_first_moment),

    Sigma_0(zeroth_moment_val),

    alpha("alpha"),
    f_approx("f-approx"),

    f_measured("f-measured")
  {
    basis_function_t::initialize(parameters);

    alpha.reset();

    alpha_ptr.reset();

    gradient.reset();
    gradient_1.reset();
    gradient_2.reset();

    a_new.reset();
    a_new_cpy.reset();

    A_mat.reset();

    initialize();

    //compute_singular_values();
  }

  template<class parameters_type, class basis_function_t>
  continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::~continuous_pole_expansion()
  {
    delete [] A1;
    delete [] A2;
  }

  template<class parameters_type, class basis_function_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::write(std::string file_name)
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
  void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::write(IO::writer<DATA_FORMAT>& writer)
  {
    writer.open_group("CPE-functions");

    writer.execute(alpha);
    writer.execute(f_approx);
    writer.execute(f_measured);

    writer.close_group();
  }

  /*
    template<class parameters_type, class basis_function_t>
    template<typename stream_type>
    void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::to_JSON(stream_type& ss)
    {
    alpha.to_JSON(ss);
    ss << ",";

    f_approx.to_JSON(ss);
    ss << ",";

    f_measured.to_JSON(ss);
    ss << ",";
    }
  */

  template<class parameters_type, class basis_function_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::initialize()
  {
    { // initialize bounds for matrices A_1 and A_2
      //w_REAL              w_REAL_dmn;
      //     alpha_dmn_t         alpha_dmn;
      //     std::pair<int, int> bounds = concurrency.get_bounds(alpha_dmn);

      //     lower_bound = bounds.first;
      //     upper_bound = bounds.second;

      //     bounds_difference = bounds.second-bounds.first;

      lower_bound = 0;
      upper_bound = alpha_dmn_t::dmn_size();

      bounds_difference = upper_bound-lower_bound;

      assert(bounds_difference>=0);
    }

    { // transfer-matrices
      A1 = new scalartype[bounds_difference*w_IMAG::dmn_size()];
      A2 = new scalartype[bounds_difference*w_IMAG::dmn_size()];

      for(int n=0; n<w_IMAG::dmn_size(); n++)
        {
          std::complex<double> z(0., w_IMAG::get_elements()[n]);

          for(int xi=0; xi<bounds_difference; xi++)
            {
              int I = xi+lower_bound;

              A1[xi + bounds_difference*n] = real(basis_function_t::phi(I, z));
              A2[xi + bounds_difference*n] = imag(basis_function_t::phi(I, z));
            }
        }
    }

    {
      for(int wn_ind=0; wn_ind<w_IMAG::dmn_size(); wn_ind++){

        std::complex<double> z(0., w_IMAG::get_elements()[wn_ind]);

        for(int alpha_ind=0; alpha_ind<alpha_dmn_t::dmn_size(); alpha_ind++)
          A_mat(wn_ind, alpha_ind) = basis_function_t::phi(alpha_ind, z);
      }
    }
  }

  /*
    template<class parameters_type, class basis_function_t>
    void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::compute_singular_values()
    {
    cout << __PRETTY_FUNCTION__ << endl;

    for(int i=0; i<alpha_dmn_t::dmn_size(); i++)
    cout << basis_function_t::volume(i) << endl;
    cout << endl;

    std::complex<double> I(0., 1);

    for(int l=0; l<10; l++)
    {
    for(int j=0; j<w_IMAG::dmn_size(); j++){
    std::complex<double> z(0., w_IMAG::get_elements()[j]);
    for(int i=0; i<alpha_dmn_t::dmn_size(); i++)
    A_mat(i,j) = (0.01+l*std::fabs(w_REAL::get_elements()[i]))*basis_function_t::phi(i, z);
    }

    singular_value_decomposition_plan<std::complex<double>, GENERAL> svd_obj(alpha_dmn_t::dmn_size(), w_IMAG::dmn_size(), 'N', 'N');
    memcpy(svd_obj.A, &A_mat(0,0), sizeof(std::complex<double>)*alpha_dmn_t::dmn_size()*w_IMAG::dmn_size() );
    svd_obj.execute_plan();

    for(int l0=0; l0<20; ++l0)
    cout << svd_obj.S[l0] << endl;
    cout << endl;
    }

    throw std::logic_error(__FUNCTION__);
    }
  */

  template<class parameters_type, class basis_function_t>
  template<typename target_dmn_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::execute(FUNC_LIB::function<std::complex<double>,       nu_nu_k_DCA_w             >& f_source,
                                                                                              FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,k_DCA,target_dmn_t> >& f_target,
                                                                                              bool                                                              use_previous_result)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    const static int dmn_number = 5;

    b_k_DCA  parallelized_dmn;
    //nu_k_DCA  parallelized_dmn;
    //nu_nu_k_DCA  parallelized_dmn;
    target_dmn_t target_dmn;

    int Nb_sbdms    = f_source.signature();

    int* coordinate = new int[Nb_sbdms];
    memset(coordinate,0,sizeof(int)*Nb_sbdms);

    std::complex<double>* input_values  = new std::complex<double>[f_source[dmn_number] ];
    std::complex<double>* output_values = new std::complex<double>[f_target[dmn_number] ];

    thread_manager_sum<concurrency_type> sum_manager(concurrency);

    do
      {
        std::pair<int, int> dmn_bounds = sum_manager.get_bounds(parallelized_dmn);

        for(int l=dmn_bounds.first; l<dmn_bounds.second; l++)
          {
            /*
              int linind = l;
              for(int j=Nb_sbdms-1; j>-1; j--)
              {
              if(j != dmn_number)
              {
              coordinate[j] = linind % f_source[j];
              linind = (linind-coordinate[j])/f_source[j];
              }
              }
            */

            int linind = l;
            int coor[2];
            parallelized_dmn.linind_2_subind(linind, coor);

            coordinate[0] = coor[0]; coordinate[1] = 0;
            coordinate[2] = coor[0]; coordinate[3] = 0;
            coordinate[4] = coor[1];
            coordinate[5] = 0;

            f_source.slice(dmn_number, coordinate, input_values);

            bool to_be_done = spline_input_values_to_wn(input_values);

            if(to_be_done && coordinate[0]==coordinate[2] && coordinate[1]==coordinate[3])
              {
                f_measured.distribute(dmn_number, coordinate, &input_values[w::dmn_size()/2]);

                perform_continuous_pole_expansion();

                //            f_wn_distribute(dmn_number, coordinate);
                //            alpha.distribute(dmn_number, coordinate, &alpha_ptr(0));

                spline_alpha_to_output_values(target_dmn, output_values);

                {
                  // spin up
                  coordinate[1] = 0; coordinate[3] = 0;

                  f_wn_distribute(dmn_number, coordinate);
                  alpha.distribute(dmn_number, coordinate, &alpha_ptr(0));
                  f_target.distribute(dmn_number, coordinate, output_values);

                  // spin dn
                  coordinate[1] = 1; coordinate[3] = 1;

                  f_wn_distribute(dmn_number, coordinate);
                  alpha.distribute(dmn_number, coordinate, &alpha_ptr(0));
                  f_target.distribute(dmn_number, coordinate, output_values);
                }
              }
          }
      }
    while(!sum_manager.sum_and_check(f_target));

    concurrency.sum(alpha);
    concurrency.sum(f_approx);
    concurrency.sum(f_measured);

    delete [] coordinate;
    delete [] input_values;
    delete [] output_values;

    symmetrize::execute(f_target);
  }

  template<class parameters_type, class basis_function_t>
  bool continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::spline_input_values_to_wn(std::complex<double>* input_values)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);
    //assert(w_IMAG::dmn_size() == w::dmn_size()/2);

    bool to_be_done = false;
    for(int n=0; n<w_IMAG::dmn_size(); n++){

      F_wn(n) = input_values[w::dmn_size()/2+n];

      F_wn_re_ptr(n) = real(input_values[w::dmn_size()/2+n]);
      F_wn_im_ptr(n) = imag(input_values[w::dmn_size()/2+n]);

      if(abs(input_values[n]) > 1.e-6)
        to_be_done = true;
    }

    return to_be_done;
  }

  template<class parameters_type, class basis_function_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::f_wn_distribute(int dmn_number, int* coordinate)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);
    //assert(w_IMAG::dmn_size() == w::dmn_size()/2);

    //   std::complex<double>* f_wn_ptr = new std::complex<double>[w_IMAG::dmn_size()];

    //   for(int n=0; n<w_IMAG::dmn_size(); n++){
    //     real(f_wn_ptr[n]) = f_wn_re_ptr(n);
    //     imag(f_wn_ptr[n]) = f_wn_im_ptr(n);
    //   }

    //   f_approx.distribute(dmn_number, coordinate, f_wn_ptr);

    //   delete [] f_wn_ptr;

    f_approx.distribute(dmn_number, coordinate, &f_wn(0));
  }

  template<class parameters_type, class basis_function_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::perform_continuous_pole_expansion()
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    Sigma_0      = 0;
    grad_Sigma_0 = 0;

    //   for(int n=0; n<alpha_dmn_t::dmn_size(); n++)
    //     alpha_ptr(n) = 0;

    for(int n=0; n<alpha_dmn_t::dmn_size(); n++)
      alpha_ptr(n) = 1./double(alpha_dmn_t::dmn_size());

    int MAX_ITERATIONS = parameters.get_max_iterations();
    int CPE_iteration  = 0;

    while(MAX_ITERATIONS>CPE_iteration)
      {
        if(!static_zero_moment)
          find_new_Sigma_0();

        int index = find_new_alpha_2();

        if(index == 0)
          break;

        CPE_iteration++;
      }

    cout << CPE_iteration << "\t" << MAX_ITERATIONS << "\t" << evaluate_Lambda_norm(alpha_ptr) << "\n";
  }

  /*
    template<class parameters_type, class basis_function_t>
    void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::compute_gradient_Sigma_0(bool normalize)
    {
    project_from_real_axis_to_imaginary_axis(&alpha_ptr(0), &f_wn_re_ptr(0), A1);

    for(int n=0; n<w_IMAG::dmn_size(); n++)
    f_wn_re_ptr(n) = F_wn_re_ptr(n)-(f_wn_re_ptr(n)+Sigma_0);

    grad_Sigma_0=0.;
    for(int n=0; n<w_IMAG::dmn_size(); n++)
    grad_Sigma_0 += -2.*f_wn_re_ptr(n);
    }

    template<class parameters_type, class basis_function_t>
    void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::find_new_Sigma_0(double LAMBDA_MAX)
    {
    for(int i=0; i<10; i++){
    compute_gradient_Sigma_0();
    Sigma_0 = Sigma_0 - grad_Sigma_0/double(2*w_IMAG::dmn_size());
    }
    }
  */

  template<class parameters_type, class basis_function_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::compute_gradient_Sigma_0(bool normalize)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    project_from_real_axis_to_imaginary_axis(&alpha_ptr(0), &f_wn_re_ptr(0), A1);

    for(int n=w_IMAG::dmn_size()/2; n<w_IMAG::dmn_size(); n++)
      f_wn_re_ptr(n) = F_wn_re_ptr(n)-(f_wn_re_ptr(n)+Sigma_0);

    grad_Sigma_0=0.;
    for(int n=w_IMAG::dmn_size()/2; n<w_IMAG::dmn_size(); n++)
      grad_Sigma_0 += -2.*f_wn_re_ptr(n);
  }

  template<class parameters_type, class basis_function_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::find_new_Sigma_0(double LAMBDA_MAX)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    for(int i=0; i<10; i++){
      compute_gradient_Sigma_0();
      Sigma_0 = Sigma_0 - grad_Sigma_0/scalartype(2.*w_IMAG::dmn_size()/2.);
    }
  }


  template<class parameters_type, class basis_function_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::compute_gradient_alphas(bool normalize)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    //     grad_F = -2*A1*(real(spoints) - transpose(A1)*a) ...
    //               - 2*A2*(imag(spoints) - transpose(A2)*a);

    {// compute the gradient
      {// -2*A1*(real(spoints) - transpose(A1)*a)
        project_from_real_axis_to_imaginary_axis(&alpha_ptr(0), &f_wn_re_ptr(0), A1);

        for(int n=0; n<w_IMAG::dmn_size(); n++)
          f_wn_re_ptr(n) = F_wn_re_ptr(n)-(f_wn_re_ptr(n)+Sigma_0);

        project_from_imaginary_axis_to_real_axis(&f_wn_re_ptr(0), &gradient_1(0), A1);

        for(int n=0; n<alpha_dmn_t::dmn_size(); n++)
          gradient_1(n) = (-2.)*gradient_1(n);
      }

      {// -2*A2*(imag(spoints) - transpose(A2)*a)
        project_from_real_axis_to_imaginary_axis(&alpha_ptr(0), &f_wn_im_ptr(0), A2);

        for(int n=0; n<w_IMAG::dmn_size(); n++)
          f_wn_im_ptr(n) = F_wn_im_ptr(n)-f_wn_im_ptr(n);

        project_from_imaginary_axis_to_real_axis(&f_wn_im_ptr(0), &gradient_2(0), A2);

        for(int n=0; n<alpha_dmn_t::dmn_size(); n++)
          gradient_2(n) = (-2.)*gradient_2(n);
      }
    }

    if(normalize)
      { // normalize
        double norm_gradient = 0;
        for(int n=0; n<alpha_dmn_t::dmn_size(); n++)
          norm_gradient += (gradient_1(n) + gradient_2(n))*(gradient_1(n) + gradient_2(n));

        for(int n=0; n<alpha_dmn_t::dmn_size(); n++)
          gradient(n) = -(gradient_1(n) + gradient_2(n))/sqrt(norm_gradient);
      }
  }

  /*
    template<class parameters_type, class basis_function_t>
    int continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::find_new_alpha(double LAMBDA_MAX)
    {
    cout << scientific;
    cout.precision(12);

    profiler_type profiler(concurrency, __FUNCTION__, __FILE__, __LINE__);

    compute_gradient_alphas();

    double L2_norm   = 1;
    double L2_norm_0 = 0;

    memcpy(&a_new_cpy(0), &alpha_ptr(0), sizeof(scalartype)*alpha_dmn_t::dmn_size());

    double lambda       = 0;
    double delta_lambda = LAMBDA_MAX/100.;

    int index = 0;

    while(true)
    {
    for(int n=0; n<alpha_dmn_t::dmn_size(); n++){
    scalartype value = alpha_ptr(n) + lambda*gradient(n);
    value < 0 ? a_new(n) = 0 : a_new(n) = value;
    }

    L2_norm = 0;

    if(false)
    {
    project_from_real_axis_to_imaginary_axis(&a_new(0), &f_wn_re_ptr(0), A1);
    project_from_real_axis_to_imaginary_axis(&a_new(0), &f_wn_im_ptr(0), A2);

    for(int n=0; n<w_IMAG::dmn_size(); n++)
    f_wn_re_ptr(n) += Sigma_0;

    for(int n=0; n<w_IMAG::dmn_size(); n++)
    L2_norm += square(f_wn_re_ptr(n)-F_wn_re_ptr(n)) + square(f_wn_im_ptr(n)-F_wn_im_ptr(n));
    }
    else
    L2_norm = evaluate_Lambda_norm();

    if(index == 0 || (L2_norm_0 - L2_norm)/L2_norm > 1.e-6){
    L2_norm_0 = L2_norm;
    memcpy(&a_new_cpy(0), &a_new(0), sizeof(scalartype)*alpha_dmn_t::dmn_size());
    }
    else
    break;

    //       if(concurrency.id()==0)
    //  cout << "\t" << index << "\t" << lambda << "\t" << L2_norm << endl;

    lambda += delta_lambda;

    index++;

    if(index>20)
    delta_lambda = LAMBDA_MAX/10.;

    if(index>100)
    delta_lambda = LAMBDA_MAX/1.;
    }

    if(index==0)
    memcpy(&alpha_ptr(0), &a_new_cpy(0), sizeof(scalartype)*alpha_dmn_t::dmn_size());
    else
    smooth_alpha();

    return index;
    }
  */

  template<class parameters_type, class basis_function_t>
  int continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::find_new_alpha_2(double LAMBDA_MAX)
  {
    cout << scientific;
    cout.precision(12);

    profiler_type profiler( __FUNCTION__, __FILE__, __LINE__);

    compute_gradient_alphas();

    memcpy(&a_new_cpy(0), &alpha_ptr(0), sizeof(scalartype)*alpha_dmn_t::dmn_size());

    double lambda       = 0.;
    double delta_lambda = LAMBDA_MAX/10.;

    int index = 0;

    double x[100];
    double y[100];

    for(int tmp=0; tmp<100; tmp++)
      {
        x[tmp] = lambda+tmp*delta_lambda;

        for(int n=0; n<alpha_dmn_t::dmn_size(); n++){
          scalartype value = alpha_ptr(n) + x[tmp]*gradient(n);
          a_new(n)         = value<0.? 0. : value;
        }

        y[tmp] = evaluate_Lambda_norm(a_new);

        if(tmp>1 && y[tmp-1]<y[tmp])
          break;
        else
          index++;
      }

    if(index==100)
      {
        lambda = x[99];
      }
    else
      {
        solve_plan<scalartype> pln(3, 1);

        pln.solved_matrix[0] = y[index-2];
        pln.solved_matrix[1] = y[index-1];
        pln.solved_matrix[2] = y[index-0];

        pln.matrix[0] = square(x[index-2]); pln.matrix[3] = x[index-2]; pln.matrix[6] = 1.;
        pln.matrix[1] = square(x[index-1]); pln.matrix[4] = x[index-1]; pln.matrix[7] = 1.;
        pln.matrix[2] = square(x[index-0]); pln.matrix[5] = x[index-0]; pln.matrix[8] = 1.;

        pln.execute_plan();

        lambda = -pln.solved_matrix[1]/(2.*pln.solved_matrix[0]);

        lambda = lambda>0? lambda : 0;
      }

    for(int n=0; n<alpha_dmn_t::dmn_size(); n++){
      scalartype value = alpha_ptr(n);
      a_new(n)         = value<0? 0 : value;
    }

    scalartype L2_norm_0 = evaluate_Lambda_norm(a_new);

    for(int n=0; n<alpha_dmn_t::dmn_size(); n++){
      scalartype value = alpha_ptr(n) + lambda*gradient(n);
      a_new(n)         = value<0? 0 : value;
    }

    scalartype L2_norm_lambda = evaluate_Lambda_norm(a_new);

    if(L2_norm_lambda<L2_norm_0){

      for(int n=0; n<alpha_dmn_t::dmn_size(); n++){
        scalartype value = alpha_ptr(n) + lambda*gradient(n);
        a_new_cpy(n) = value<0? 0 : value;
      }

      //     if(concurrency.id()==0)
      //       cout << L2_norm_lambda << "\n";
    }
    else
      index=0;

    smooth_alpha();

    return index;
  }


  template<class parameters_type, class basis_function_t>
  double continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::evaluate_Lambda_norm(FUNC_LIB::function<scalartype, alpha_dmn_t>& alpha_real)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    static FUNC_LIB::function<std::complex<scalartype>, alpha_dmn_t> alpha_complex("alpha_complex");
    for(int n=0; n<alpha_dmn_t::dmn_size(); n++)
      alpha_complex(n) = alpha_real(n);

    static gemv_plan<std::complex<scalartype> > gemv_obj(w_IMAG::dmn_size(), alpha_dmn_t::dmn_size());

    gemv_obj.alpha         =  1;
    gemv_obj.beta          =  0.;

    gemv_obj.matrix        = &A_mat(0,0);
    gemv_obj.vector_source = &alpha_complex(0);
    gemv_obj.vector_target = &f_wn(0);

    gemv_obj.execute_plan();

    scalartype L2_norm = 0;
    for(int n=0; n<w_IMAG::dmn_size(); n++){
      f_wn(n) += Sigma_0;
      L2_norm += square(real((f_wn(n)-F_wn(n)))) + square(imag(f_wn(n)-F_wn(n)));
    }

    return L2_norm;
  }

  template<class parameters_type, class basis_function_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::project_from_real_axis_to_imaginary_axis(scalartype* real_values, scalartype* imag_values, scalartype* matrix)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    gemv_plan<scalartype> gemv_obj(bounds_difference, w_IMAG::dmn_size());

    for(int l=0; l<w_IMAG::dmn_size(); l++)
      imag_values[l]=0.;

    gemv_obj.TRANS         = 'T';

    gemv_obj.alpha         =  1;
    gemv_obj.beta          =  0.;

    gemv_obj.matrix        = matrix;
    gemv_obj.vector_source = &real_values[lower_bound];
    gemv_obj.vector_target = imag_values;

    gemv_obj.execute_plan();
  }

  template<class parameters_type, class basis_function_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::project_from_imaginary_axis_to_real_axis(scalartype* imag_values, scalartype* real_values, scalartype* matrix)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    gemv_plan<scalartype> gemv_obj(bounds_difference, w_IMAG::dmn_size());

    for(int l=0; l<alpha_dmn_t::dmn_size(); l++)
      real_values[l]=0.;

    gemv_obj.TRANS         = 'N';

    gemv_obj.alpha         =  1;
    gemv_obj.beta          =  0.;

    gemv_obj.matrix        = matrix;
    gemv_obj.vector_source = imag_values;
    gemv_obj.vector_target = &real_values[lower_bound];

    gemv_obj.execute_plan();
  }

  template<class parameters_type, class basis_function_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::smooth_alpha()
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    int factor = parameters.get_smoothing_factor();

    for(int n=0; n<w_REAL::dmn_size(); n++){

      double N     = 0.;
      alpha_ptr(n) = 0.;

      for(int l=0-factor; l<=0+factor; l++)
        if(n+l>-1 && n+l<w_REAL::dmn_size()){
          alpha_ptr(n) += a_new_cpy(n+l);
          N += 1.;
        }

      alpha_ptr(n) /= N;
    }
  }

  template<class parameters_type, class basis_function_t>
  template<typename target_dmn_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::spline_alpha_to_output_values(target_dmn_t& target_dmn, std::complex<double>* output_values)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    if(IS_EQUAL_TYPE<target_dmn_t, w>::check){
      spline_alpha_to_output_values_w_IMAG(output_values);
      return;
    }

    if(IS_EQUAL_TYPE<target_dmn_t, w_REAL>::check){
      spline_alpha_to_output_values_w_REAL(output_values);
      return;
    }

    throw std::logic_error(__FUNCTION__);
  }

  template<class parameters_type, class basis_function_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::spline_alpha_to_output_values_w_IMAG(std::complex<double>* output_values)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    double sum=0;

    for(int n=0; n<alpha_dmn_t::dmn_size(); n++)
      sum += alpha_ptr(n)*basis_function_t::volume(n);

    int left_index=0;
    int right_index=alpha_dmn_t::dmn_size()-1;
    double num=0;

    while(alpha_ptr(left_index)<1.e-6){
      left_index++;
      num++;
    }

    while(alpha_ptr(right_index)<1.e-6){
      right_index--;
      num++;
    }

    for(int n=0; n<left_index; n++)
      alpha_ptr(n) += (1-sum)/basis_function_t::volume(n)/num;//(left_index+alpha_dmn_t::dmn_size()-right_index);

    for(int n=alpha_dmn_t::dmn_size()-1; n>right_index; n--)
      alpha_ptr(n) += (1-sum)/basis_function_t::volume(n)/num;//(left_index+alpha_dmn_t::dmn_size()-right_index);

    sum=0;

    for(int n=0; n<alpha_dmn_t::dmn_size(); n++)
      sum += alpha_ptr(n)*basis_function_t::volume(n);

    concurrency << "\n\t sum : ";
    concurrency << sum;
    concurrency << "\t  m_0 : ";
    concurrency << Sigma_0;
    concurrency << "\n";

    for(int m=0; m<w::dmn_size(); m++)
      {
        output_values[m] = 0;

        std::complex<double> z(0., w::get_elements()[m]);

        for(int n=0; n<alpha_dmn_t::dmn_size(); n++){
          output_values[m] += double(alpha_ptr(n))*basis_function_t::phi(n, z);//sum;
        }
      }
  }

  template<class parameters_type, class basis_function_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::spline_alpha_to_output_values_w_REAL(std::complex<double>* output_values)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    double real_axis_off_set = parameters.get_real_frequencies_off_set();

    for(int m=0; m<w_REAL::dmn_size(); m++)
      {
        output_values[m] = Sigma_0;

        std::complex<double> z(w_REAL::get_elements()[m], real_axis_off_set);//parameters.get_EPSILON());

        for(int n=0; n<alpha_dmn_t::dmn_size(); n++){
          output_values[m] += double(alpha_ptr(n))*basis_function_t::phi(n, z);
        }
      }
  }

}

#endif


































































// template<class parameters_type, class basis_function_t>
// void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::execute(FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w>& f_source,
//                                                                                                            FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w>& f_target,
//                                                                                                            bool                                           use_previous_result)
// {
//   concurrency << "\n\n\t --> condition_Greens_function \n\n";

//   static_zero_moment    = true;
//   project_to_unit_plane = true;

//   Sigma_0 = 0.;

//   perform_continuous_pole_expansion(f_source, f_target);

//   static_zero_moment    = false;
//   project_to_unit_plane = false;
// }

// template<class parameters_type, class basis_function_t>
// void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::execute(FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w>& f_source,
//                                                                                                            FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w>& f_target,
//                                                                                                            bool                                           use_previous_result)
// {
//   const static int dmn_number = 5;

//   concurrency << "\n\n\t --> condition_Greens_function \n\n";

//   static_zero_moment    = true;
//   project_to_unit_plane = true;

//   Sigma_0 = 0.;

//   w   target_dmn;

//   int Nb_sbdms    = f_source.signature();
//   int Nb_elements = f_source.size();

//   int* coordinate = new int[Nb_sbdms];
//   memset(coordinate,0,sizeof(int)*Nb_sbdms);

//   std::complex<double>* input_values    = new std::complex<double>[f_source[dmn_number] ];
//   std::complex<double>* output_values   = new std::complex<double>[f_target[dmn_number] ];

//   {
//     int Nb_Pade = Nb_elements/f_source[dmn_number];

//     for(int l=0; l<Nb_Pade; l++)
//       {
//      int linind = l;
//      for(int j=Nb_sbdms-1; j>-1; j--)
//        {
//          if(j != dmn_number)
//            {
//              coordinate[j] = linind % f_source[j];
//              linind = (linind-coordinate[j])/f_source[j];
//            }
//        }

//      f_source.slice(dmn_number, coordinate, input_values);

//      bool to_be_done = spline_input_values_to_wn(input_values);

//      if(to_be_done)
//        {
//          S_wn.distribute(dmn_number, coordinate, &input_values[w::dmn_size()/2]);

//          perform_continuous_pole_expansion();

//          a_wn_distribute(dmn_number, coordinate);

//          a_x.distribute(dmn_number, coordinate, &a_x_ptr(0));

//          spline_a_x_to_output_values(target_dmn, output_values);

//          f_target.distribute(dmn_number, coordinate, output_values);
//        }
//       }
//   }

//   delete [] coordinate;
//   delete [] input_values;
//   delete [] output_values;

//   static_zero_moment    = false;
//   project_to_unit_plane = false;
// }

// template<class parameters_type, class basis_function_t>
// void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::condition_Greens_function(FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w>& f_source,
//                                                                                                                              FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w>& f_target)
// {
//   concurrency << "\n\n\t --> condition_Greens_function \n\n";

//   static_zero_moment    = true;
//   project_to_unit_plane = true;

//   Sigma_0 = 0.;

//   perform_continuous_pole_expansion(f_source, f_target);

//   static_zero_moment    = false;
//   project_to_unit_plane = false;
// }

// template<class parameters_type, class basis_function_t>
// void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::perform_continuous_pole_expansion(FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w>& f_source,
//                                                                                                                                      FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w>& f_target)
// {
//   const static int dmn_number = 5;

//   w   target_dmn;

//   int Nb_sbdms    = f_source.signature();
//   int Nb_elements = f_source.size();

//   int* coordinate = new int[Nb_sbdms];
//   memset(coordinate,0,sizeof(int)*Nb_sbdms);

//   std::complex<double>* input_values    = new std::complex<double>[f_source[dmn_number] ];
//   std::complex<double>* output_values   = new std::complex<double>[f_target[dmn_number] ];

//   {
//     int Nb_Pade = Nb_elements/f_source[dmn_number];

//     for(int l=0; l<Nb_Pade; l++)
//       {
//      int linind = l;
//      for(int j=Nb_sbdms-1; j>-1; j--)
//        {
//          if(j != dmn_number)
//            {
//              coordinate[j] = linind % f_source[j];
//              linind = (linind-coordinate[j])/f_source[j];
//            }
//        }

//      f_source.slice(dmn_number, coordinate, input_values);

//      bool to_be_done = spline_input_values_to_wn(input_values);

//      if(to_be_done)
//        {
//          S_wn.distribute(dmn_number, coordinate, &input_values[w::dmn_size()/2]);

//          perform_continuous_pole_expansion();

//          a_wn_distribute(dmn_number, coordinate);

//          a_x.distribute(dmn_number, coordinate, &a_x_ptr(0));

//          spline_a_x_to_output_values(target_dmn, output_values);

//          f_target.distribute(dmn_number, coordinate, output_values);
//        }
//       }
//   }

//   delete [] coordinate;
//   delete [] input_values;
//   delete [] output_values;
// }




// template<class parameters_type, class basis_function_t>
// template<typename MOMS_t>
// void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::compute_spectral_function(MOMS_t& MOMS)
// {
//   concurrency << "\n\n\t start perform continuous-pole-expansion G \n\n";

//   perform_continuous_pole_expansion(MOMS.Sigma, Sigma);

//   concurrency << "\n\n\t start coarsegraining G0 \n\n";

//   FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w_REAL> Sigma_zero;
//   coarsegrain_obj.compute_G_from_H_and_Sigma(MOMS.H_LDA, Sigma_zero, G0_k_w);

//   concurrency << "\n\n\t start coarsegraining G \n\n";

//   coarsegrain_obj.compute_G_from_H_and_Sigma(MOMS.H_LDA, Sigma, G_k_w);

//   concurrency << "\n\n\t compute spectra \n\n";

//   for(int b=0; b<b::dmn_size(); b++){
//     for(int s=0; s<s::dmn_size(); s++){

//       for(int w_ind=0; w_ind<w_REAL::dmn_size(); w_ind++){

//      for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++){
//        A_nu_w (b,s,w_ind) += -1./M_PI*imag(G_k_w (b,s,b,s,k_ind,w_ind));
//        A0_nu_w(b,s,w_ind) += -1./M_PI*imag(G0_k_w(b,s,b,s,k_ind,w_ind));

//        A_w (w_ind) += -1./M_PI*imag(G_k_w (b,s,b,s,k_ind,w_ind));
//        A0_w(w_ind) += -1./M_PI*imag(G0_k_w(b,s,b,s,k_ind,w_ind));
//      }
//       }
//     }
//   }

//   A_nu_w  /= (double(k_DCA::dmn_size()));
//   A0_nu_w /= (double(k_DCA::dmn_size()));

//   A_w  /= (double(k_DCA::dmn_size())*b::dmn_size()*s::dmn_size());
//   A0_w /= (double(k_DCA::dmn_size())*b::dmn_size()*s::dmn_size());
// }

// template<class parameters_type, class basis_function_t>
// void continuous_pole_expansion<parameters_type, basis_function_t, GRADIENT_METHOD>::perform_continuous_pole_expansion(FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w     >& f_source,
//                                                                                                                                      FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w_REAL>& f_target)
// {
//   const static int dmn_number = 5;

//   int Nb_sbdms    = f_source.signature();
//   int Nb_elements = f_source.size();

//   int* coordinate = new int[Nb_sbdms];
//   memset(coordinate,0,sizeof(int)*Nb_sbdms);

//   std::complex<double>* input_values    = new std::complex<double>[f_source[dmn_number] ];
//   std::complex<double>* output_values   = new std::complex<double>[f_target[dmn_number] ];

//   {
//     int Nb_Pade = Nb_elements/f_source[dmn_number];

//     for(int l=0; l<Nb_Pade; l++)
//       {
//      concurrency << "\n\t\t ";
//      concurrency << int(double(l)/double(Nb_Pade)*100);
//      concurrency << " % completed";

//      int linind = l;
//      for(int j=Nb_sbdms-1; j>-1; j--)
//        {
//          if(j != dmn_number)
//            {
//              coordinate[j] = linind % f_source[j];
//              linind = (linind-coordinate[j])/f_source[j];
//            }
//        }

//      f_source.slice(dmn_number, coordinate, input_values);

//      bool to_be_done = spline_input_values_to_wn(input_values);

//      if(to_be_done)
//        {
//          S_wn.distribute(dmn_number, coordinate, &input_values[w::dmn_size()/2]);

//          perform_continuous_pole_expansion();

//          a_wn_distribute(dmn_number, coordinate);

//          a_x.distribute(dmn_number, coordinate, &a_x_ptr(0));

//          spline_a_x_to_output_values(output_values);

//          f_target.distribute(dmn_number, coordinate, output_values);
//        }
//       }
//   }

//   delete [] coordinate;
//   delete [] input_values;
//   delete [] output_values;
// }
