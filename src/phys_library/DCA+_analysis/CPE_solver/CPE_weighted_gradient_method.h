//-*-C++-*-

#ifndef CONTINUOUS_POLE_EXPANSION_WEIGHTED_GRADIENT_METHOD_H
#define CONTINUOUS_POLE_EXPANSION_WEIGHTED_GRADIENT_METHOD_H

namespace DCA
{
  /*!
   *  \ingroup CPE
   *
   *  \author  Peter Staar
   *  \brief   This class implements a CPE analytic continuation, using a gradient method as a minimization-algorithm.
   *  \version 1.0
   */
  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  class continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>
  {
#include "type_definitions.h"

    enum    gradient_minimization_state {BELOW_MAX_ERROR, FOUND_THE_MINIMUM, FOUND_A_MINIMUM, END_OF_PATH};
    typedef gradient_minimization_state gradient_minimization_state_t;

  public:

    typedef double                                     scalartype;

    typedef typename parameters_type::profiler_type    profiler_type;
    typedef typename parameters_type::concurrency_type concurrency_type;

    typedef dmn_0<basis_function_t>          alpha_dmn_t;

    typedef dmn_3<nu,nu,k_dmn_t            > nu_nu_k_dmn;

    typedef dmn_4<nu,nu,k_dmn_t,w     >      nu_nu_k_dmn_w;
    typedef dmn_4<nu,nu,k_dmn_t,w_IMAG>      nu_nu_k_dmn_w_IMAG;

    typedef dmn_4<nu,nu,k_dmn_t,alpha_dmn_t> nu_nu_k_dmn_alpha_dmn;

    typedef CPE_data<scalartype, basis_function_t, k_dmn_t, w_dmn_t> CPE_data_type;

  public:

    continuous_pole_expansion(parameters_type&  parameters,
                              concurrency_type& concurrency,
                              bool              fixed_zero_moment=false,
                              double            zero_moment=0,
                              bool              fixed_first_moment=false,
                              double            first_moment=1);

    ~continuous_pole_expansion();

    FUNC_LIB::function<             scalartype , nu_nu_k_dmn>&        get_error_function() {return error_function;}

    FUNC_LIB::function<std::complex<scalartype>, nu_nu_k_dmn_w_IMAG>& get_f_approx()   { return f_approx;   }
    FUNC_LIB::function<std::complex<scalartype>, nu_nu_k_dmn_w_IMAG>& get_f_measured() { return f_measured; }

    FUNC_LIB::function<std::complex<scalartype>, nu_nu_k_dmn_w>& get_S_approx()   { return S_approx;   }

    void execute_st(FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w      > >& f_source,
                    FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& f_target);

    void execute_mt(FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w      > >& f_source,
                    FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& f_target);

    static void* threaded_analytical_continuation(void* void_ptr);

    void write(std::string filename);

    template<IO::FORMAT DATA_FORMAT>
    void write(IO::writer<DATA_FORMAT>& reader);

  private:

    void initialize();

    static void read_function_values(int nu_ind,
                                     int k_ind,
                                     CPE_data_type& CPE_data_obj);

    static void perform_continuous_pole_expansion_threaded_1(CPE_data_type& CPE_data_obj);
    static void perform_continuous_pole_expansion_threaded_2(CPE_data_type& CPE_data_obj);

    static void write_function_values(int nu_ind,
                                      int k_ind,
                                      CPE_data_type& CPE_data_obj);

    static void compute_gradient_alphas(CPE_data_type& CPE_data_obj);

    static void find_new_Sigma_0        (CPE_data_type& CPE_data_obj);

    static void compute_gradient_Sigma_0(CPE_data_type& CPE_data_obj);

    static gradient_minimization_state find_new_alpha_1(CPE_data_type& CPE_data_obj);
    static int                         find_new_alpha_2(CPE_data_type& CPE_data_obj);

    static scalartype find_minimum_1(int index, CPE_data_type& CPE_data_obj);
    static scalartype find_minimum_2(int index, CPE_data_type& CPE_data_obj);

    static void compute_new_alpha(double lambda, CPE_data_type& CPE_data_obj);

    static double evaluate_Lambda_norm(CPE_data_type& CPE_data_obj);

    static void project_from_real_axis_to_imaginary_axis(LIN_ALG::vector<scalartype, LIN_ALG::CPU>& real_values,
                                                         LIN_ALG::vector<scalartype, LIN_ALG::CPU>& imag_values,
                                                         LIN_ALG::matrix<scalartype, LIN_ALG::CPU>& matrix);

    static void project_from_imaginary_axis_to_real_axis(LIN_ALG::vector<scalartype, LIN_ALG::CPU>& imag_values,
                                                         LIN_ALG::vector<scalartype, LIN_ALG::CPU>& real_values,
                                                         LIN_ALG::matrix<scalartype, LIN_ALG::CPU>& matrix);

  private:

    parameters_type&  parameters;
    concurrency_type& concurrency;

  public:

    FUNC_LIB::function<scalartype, nu_nu_k_dmn>           Sigma_0_moment;
    FUNC_LIB::function<scalartype, nu_nu_k_dmn_alpha_dmn> alpha_function;
    FUNC_LIB::function<scalartype, nu_nu_k_dmn>           error_function;

    FUNC_LIB::function<std::complex<scalartype>, nu_nu_k_dmn_w>      S_approx;

    FUNC_LIB::function<std::complex<scalartype>, nu_nu_k_dmn_w_IMAG> f_approx;
    FUNC_LIB::function<std::complex<scalartype>, nu_nu_k_dmn_w_IMAG> f_measured;

    LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> A_matrix;

    LIN_ALG::matrix<             scalartype , LIN_ALG::CPU> A_matrix_re;
    LIN_ALG::matrix<             scalartype , LIN_ALG::CPU> A_matrix_im;
  };

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::continuous_pole_expansion(parameters_type&  parameters_ref,
                                                                                                                                      concurrency_type& concurrency_ref,
                                                                                                                                      bool              /*fixed_zeroth_moment*/,
                                                                                                                                      double            /*zeroth_moment_val*/,
                                                                                                                                      bool              /*fixed_first_moment*/,
                                                                                                                                      double            /*first_moment_val*/):
    parameters(parameters_ref),
    concurrency(concurrency_ref),

    alpha_function("alpha"),
    error_function("L2-CPE-error"),

    S_approx("Sigma-approx"),

    f_approx("f-approx"),
    f_measured("f-measured"),

    A_matrix("A_matrix"),

    A_matrix_re("A_matrix_re"),
    A_matrix_im("A_matrix_im")
  {
    basis_function_t::initialize(parameters);

    initialize();
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::~continuous_pole_expansion()
  {}

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::write(std::string file_name)
  {
    IO::FORMAT FORMAT = parameters.get_output_format();

    std::cout << "\n\n\t\t start writing " << file_name << "\n\n";

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

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  template<IO::FORMAT DATA_FORMAT>
  void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::write(IO::writer<DATA_FORMAT>& /*writer*/)
  {
    /*
      writer.open_group("CPE-functions");

      writer.execute(error_function);
      writer.execute(alpha_function);

      writer.execute(f_approx);
      writer.execute(f_measured);

      writer.close_group();
    */
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::initialize()
  {
    //     std::cout << "\n\t w_IMAG     ::dmn_size() : " << w_IMAG     ::dmn_size();
    //     std::cout << "\n\t alpha_dmn_t::dmn_size() : " << alpha_dmn_t::dmn_size();

    {
      A_matrix   .resize_no_copy(std::pair<int,int>(w_IMAG::dmn_size(), alpha_dmn_t::dmn_size()));

      A_matrix_re.resize_no_copy(std::pair<int,int>(w_IMAG::dmn_size(), alpha_dmn_t::dmn_size()));
      A_matrix_im.resize_no_copy(std::pair<int,int>(w_IMAG::dmn_size(), alpha_dmn_t::dmn_size()));

      for(int wn_ind=0; wn_ind<w_IMAG::dmn_size(); wn_ind++)
        {
          std::complex<double> z(0., w_IMAG::get_elements()[wn_ind]);

          for(int alpha_ind=0; alpha_ind<alpha_dmn_t::dmn_size(); alpha_ind++)
            {
              A_matrix(wn_ind, alpha_ind) = basis_function_t::phi(alpha_ind, z);

              A_matrix_re(wn_ind, alpha_ind) = real(A_matrix(wn_ind, alpha_ind));
              A_matrix_im(wn_ind, alpha_ind) = imag(A_matrix(wn_ind, alpha_ind));
            }
        }
    }
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::execute_st(FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w      > >& f_source,
                                                                                                                            FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& f_target)
  {
    //std::cout << __FUNCTION__ << "\t" << __LINE__ << "\n";

    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    {// initialize f-measured
      for(int w_ind=0; w_ind<w_IMAG::dmn_size(); w_ind++)
        for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); k_ind++)
          for(int nu_j=0; nu_j<nu::dmn_size(); nu_j++)
            for(int nu_i=0; nu_i<nu::dmn_size(); nu_i++)
              f_measured(nu_i, nu_j, k_ind, w_ind) = f_source(nu_i, nu_j, k_ind, w::dmn_size()/2+w_ind);
    }

    CPE_data_type CPE_data_obj;

    CPE_data_obj.initialize(parameters, f_target, *this);

    //std::cout << "\n\n";
    for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); k_ind++){
      for(int nu_ind=0; nu_ind<nu::dmn_size(); nu_ind++){

        //std::cout << k_ind << "\t" << nu_ind << "\t start reading measured function\n";
        {
          read_function_values(nu_ind, k_ind, CPE_data_obj);
        }

        //std::cout << k_ind << "\t" << nu_ind << "\t start analytic continuation\n";
        {// find the alpha
          perform_continuous_pole_expansion_threaded(CPE_data_obj);
        }

        //std::cout << k_ind << "\t" << nu_ind << "\t start writing approx function\n";
        {
          write_function_values(nu_ind, k_ind, CPE_data_obj);
        }
      }
    }

    symmetrize::execute(f_target);

    //SHOW::execute(f_approx,f_measured);
  }

  /*
    template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
    void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::execute_mt(FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w      > >& f_source,
    FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& f_target)
    {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    {// initialize f-measured
    for(int w_ind=0; w_ind<w_IMAG::dmn_size(); w_ind++)
    for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); k_ind++)
    for(int nu_j=0; nu_j<nu::dmn_size(); nu_j++)
    for(int nu_i=0; nu_i<nu::dmn_size(); nu_i++)
    f_measured(nu_i, nu_j, k_ind, w_ind) = f_source(nu_i, nu_j, k_ind, w::dmn_size()/2+w_ind);
    }

    //SHOW::execute(f_measured);

    int nr_threads = std::min(8, k_dmn_t::dmn_size());

    std::vector<CPE_data_type > CPE_data_vector(nr_threads);

    for(int l=0; l<nr_threads; l++)
    CPE_data_vector[l].initialize(l, parameters, f_target, *this);

    COMP_LIB::parallelization<COMP_LIB::POSIX_LIBRARY> parallelization_obj;

    parallelization_obj.execute(nr_threads, threaded_analytical_continuation, (void*) &CPE_data_vector);

    symmetrize::execute(f_target);

    //SHOW::execute(f_approx,f_measured);
    }

    template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
    void* continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::threaded_analytical_continuation(void* void_ptr)
    {
    COMP_LIB::posix_data*       data_ptr         = static_cast<COMP_LIB::posix_data      *>(void_ptr);
    std::vector<CPE_data_type>* CPE_data_vec_ptr = static_cast<std::vector<CPE_data_type>*>(data_ptr->args);

    int id         = data_ptr->id;
    int nr_threads = data_ptr->nr_threads;

    std::vector<CPE_data_type>& CPE_data_vec = *(CPE_data_vec_ptr);

    k_dmn_t k_dmn;
    std::pair<int, int> k_bounds = COMP_LIB::parallelization<COMP_LIB::POSIX_LIBRARY>::get_bounds(id, nr_threads, k_dmn);

    //     std::stringstream ss;
    //     {
    //       ss << id << "\t" << nr_threads << "\t" << k_bounds.first << "\t" << k_bounds.second << "\t" << CPE_data_vec.size() << "\n";
    //       std::cout << ss.str();
    //     }

    for(int k_ind=k_bounds.first; k_ind<k_bounds.second; k_ind++)
    {
    for(int nu_ind=0; nu_ind<nu::dmn_size(); nu_ind++)
    {
    //             ss << k_ind << "\t" << nu_ind << "\t start reading measured function\n";
    //             std::cout << ss.str();
    {
    read_function_values(nu_ind, k_ind, CPE_data_vec[id]);
    }

    //             ss << k_ind << "\t" << nu_ind << "\t start analytic continuation\n";
    //             std::cout << ss.str();
    {// find the alpha
    perform_continuous_pole_expansion_threaded_1(CPE_data_vec[id]);
    }

    //             ss << k_ind << "\t" << nu_ind << "\t start writing approx function\n";
    //             std::cout << ss.str();
    {
    write_function_values(nu_ind, k_ind, CPE_data_vec[id]);
    }
    }
    }

    return 0;
    }
  */

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::execute_mt(FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w      > >& f_source,
                                                                                                                            FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& f_target)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    {
      S_approx = 0.;

      f_approx = 0;
      f_target = 0;

      f_measured = 0;

      Sigma_0_moment = 0;
      alpha_function = 0;
      error_function = 0;
    }

    symmetrize::execute(f_source);

    {// initialize f-measured
      for(int w_ind=0; w_ind<w_IMAG::dmn_size(); w_ind++)
        for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); k_ind++)
          for(int nu_j=0; nu_j<nu::dmn_size(); nu_j++)
            for(int nu_i=0; nu_i<nu::dmn_size(); nu_i++)
              f_measured(nu_i, nu_j, k_ind, w_ind) = f_source(nu_i, nu_j, k_ind, w::dmn_size()/2+w_ind);
    }

    //SHOW::execute(f_measured);

    {
      //       f_approx = 0;
      //       f_target = 0;

      //       Sigma_0_moment = 0;
      //       alpha_function = 0;
      //       error_function = 0;

      dmn_3<b, s, k_dmn_t> b_s_k_dmn;
      std::pair<int, int>  bounds = concurrency.get_bounds(b_s_k_dmn);

      int nr_tasks = bounds.second-bounds.first;

      if(nr_tasks>0)
        {
          int nr_threads = std::min(8, nr_tasks);

          std::vector<CPE_data_type > CPE_data_vector(nr_threads);

          for(int l=0; l<nr_threads; l++)
            CPE_data_vector[l].initialize(l, bounds, parameters, f_target, *this);

          COMP_LIB::parallelization<COMP_LIB::POSIX_LIBRARY> parallelization_obj;

          parallelization_obj.execute(nr_threads, threaded_analytical_continuation, (void*) &CPE_data_vector);
        }

      // sum
      concurrency.sum(Sigma_0_moment);
      concurrency.sum(alpha_function);
      concurrency.sum(error_function);

      concurrency.sum(S_approx);
      concurrency.sum(f_approx);
      //concurrency.sum(f_measured);

      concurrency.sum(f_target);
    }

    symmetrize::execute(f_target);

    //SHOW::execute(f_approx, f_measured);
    //SHOW::execute(f_target);
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  void* continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::threaded_analytical_continuation(void* void_ptr)
  {
    COMP_LIB::posix_data*       data_ptr         = static_cast<COMP_LIB::posix_data      *>(void_ptr);
    std::vector<CPE_data_type>* CPE_data_vec_ptr = static_cast<std::vector<CPE_data_type>*>(data_ptr->args);

    int id         = data_ptr->id;
    int nr_threads = data_ptr->nr_threads;

    std::vector<CPE_data_type>& CPE_data_vec = *(CPE_data_vec_ptr);

    std::pair<int, int> MPI_bounds = CPE_data_vec[id].bounds;

    dmn_2<b, s>          nu_dmn;
    dmn_3<b, s, k_dmn_t> b_s_k_dmn;

    std::pair<int, int>  bounds = COMP_LIB::parallelization<COMP_LIB::POSIX_LIBRARY>::get_bounds(id, nr_threads, MPI_bounds);

    //     if(id==0)
    //       {
    //         std::stringstream ss;
    //         ss << "\n\n\n\t\t"        << __FUNCTION__ << "\n";
    //  ss << "\t\t MPI-bounds   : " << MPI_bounds.first << "\t" << MPI_bounds.second << "\n";
    //  ss << "\t\t POSIX-bounds : " << id << "\t" << nr_threads << "\t" << bounds.first << "\t" << bounds.second << "\n";
    //  ss << "\n\n\n";
    //         std::cout << ss.str();
    //       }

    int coor[3];
    for(int l=bounds.first; l<bounds.second; l++)
      {
        b_s_k_dmn.linind_2_subind(l, coor);

        int nu_ind = nu_dmn(coor[0], coor[1]);
        int k_ind  = coor[2];

        {
          //         if(id==0)
          //           {
          //             std::stringstream ss;
          //             ss << nu_ind << "\t" << k_ind << "\t start reading measured function\n";
          //             std::cout << ss.str();
          //           }
          {
            read_function_values(nu_ind, k_ind, CPE_data_vec[id]);
          }

          //         if(id==0)
          //           {
          //             std::stringstream ss;
          //             ss << nu_ind << "\t" << k_ind << "\t start analytic continuation\n";
          //             std::cout << ss.str();
          //           }

          if(CPE_data_vec[id].is_interacting_band[nu_ind])
            {// find the alpha
              perform_continuous_pole_expansion_threaded_1(CPE_data_vec[id]);
            }
          else
            {
              for(int n_ind=0; n_ind<alpha_dmn_t::dmn_size(); n_ind++)
                CPE_data_vec[id].alpha_vec_d[n_ind] = 0.;
            }

          //         if(id==0)
          //           {
          //             std::stringstream ss;
          //             ss << nu_ind << "\t" << k_ind << "\t start writing approx function\n";
          //             std::cout << ss.str();
          //           }

          {
            write_function_values(nu_ind, k_ind, CPE_data_vec[id]);
          }

          //         if(id==0)
          //           {
          //             std::stringstream ss;
          //             ss << "\n";
          //             std::cout << ss.str();
          //           }
        }
      }

    return 0;
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::read_function_values(int nu_ind,
                                                                                                                                      int k_ind,
                                                                                                                                      CPE_data_type& CPE_data_obj)
  {
    FUNC_LIB::function<std::complex<scalartype>, dmn_4<nu,nu,k_dmn_t,w_IMAG> >& f_measured_func = *(CPE_data_obj.f_measured_ptr);

    {// read in the values on the imaginary axis.
      for(int w_ind=0; w_ind<w_IMAG::dmn_size(); w_ind++){

        std::complex<double> value = f_measured_func(nu_ind, nu_ind, k_ind,w_ind);

        CPE_data_obj.F_wn_vec[w_ind] = value;

        CPE_data_obj.F_wn_re_vec[w_ind] = real(value);
        CPE_data_obj.F_wn_im_vec[w_ind] = imag(value);
      }
    }
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::write_function_values(int nu_ind,
                                                                                                                                       int k_ind,
                                                                                                                                       CPE_data_type& CPE_data_obj)
  {
    {
      FUNC_LIB::function<scalartype, dmn_3<nu,nu,k_dmn_t>             >& Sigma_0_func = *(CPE_data_obj.Sigma_0_moment_ptr);
      FUNC_LIB::function<scalartype, dmn_4<nu,nu,k_dmn_t,alpha_dmn_t> >& alpha_func   = *(CPE_data_obj.alpha_function_ptr);
      FUNC_LIB::function<scalartype, dmn_3<nu,nu,k_dmn_t>             >& error_func   = *(CPE_data_obj.error_function_ptr);

      Sigma_0_func(nu_ind, nu_ind, k_ind) = CPE_data_obj.Sigma_0;

      for(int n_ind=0; n_ind<alpha_dmn_t::dmn_size(); n_ind++)
        alpha_func(nu_ind, nu_ind, k_ind, n_ind) = CPE_data_obj.alpha_vec_d[n_ind];

      error_func(nu_ind, nu_ind, k_ind) = evaluate_Lambda_norm(CPE_data_obj);
    }

    {// write in the values on the imaginary axis.
      LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& A_matrix = *(CPE_data_obj.A_matrix_ptr);

      FUNC_LIB::function<std::complex<scalartype>, nu_nu_k_DCA_w_IMAG>& f_approx_func = *(CPE_data_obj.f_approx_ptr);

      for(int m_ind=0; m_ind<w_IMAG::dmn_size(); m_ind++){

        std::complex<double> value = CPE_data_obj.Sigma_0;

        for(int n_ind=0; n_ind<alpha_dmn_t::dmn_size(); n_ind++)
          value += A_matrix(m_ind, n_ind)*CPE_data_obj.alpha_vec_d[n_ind];

        f_approx_func(nu_ind, nu_ind, k_ind, m_ind) = value;
      }
    }

    {
      FUNC_LIB::function<std::complex<scalartype>, dmn_4<nu,nu,k_dmn_t,w> >& S_approx_func = *(CPE_data_obj.S_approx_ptr);

      for(int w_ind=w::dmn_size()/2; w_ind<w::dmn_size(); w_ind++)
        {
          std::complex<scalartype> value = CPE_data_obj.Sigma_0;

          double w_re = 0;
          double w_im = w::get_elements()[w_ind];

          std::complex<double> z(w_re, w_im);

          for(int n_ind=0; n_ind<alpha_dmn_t::dmn_size(); n_ind++)
            value += basis_function_t::phi(n_ind, z)*CPE_data_obj.alpha_vec_d[n_ind];

          S_approx_func(nu_ind, nu_ind, k_ind, w_ind) = value;
        }
    }

    {
      FUNC_LIB::function<std::complex<scalartype>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> >& f_target_func = *(CPE_data_obj.f_target_ptr);

      for(int w_ind=0; w_ind<w_dmn_t::dmn_size(); w_ind++)
        {
          std::complex<scalartype> value = CPE_data_obj.Sigma_0;

          double w_re = w_dmn_t::get_elements()[w_ind];
          double w_im = CPE_data_obj.real_axis_off_set;

          std::complex<double> z(w_re, w_im);

          for(int n_ind=0; n_ind<alpha_dmn_t::dmn_size(); n_ind++)
            value += basis_function_t::phi(n_ind, z)*CPE_data_obj.alpha_vec_d[n_ind];

          f_target_func(nu_ind, nu_ind, k_ind, w_ind) = value;
        }
    }
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::perform_continuous_pole_expansion_threaded_1(CPE_data_type& CPE_data_obj)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    CPE_data_obj.Sigma_0      = 0;
    CPE_data_obj.grad_Sigma_0 = 0;

    for(int n=0; n<alpha_dmn_t::dmn_size(); n++)
      CPE_data_obj.alpha_vec_d[n] = 1./double(alpha_dmn_t::dmn_size());

    int MAX_ITERATIONS = CPE_data_obj.max_iterations;
    int CPE_iteration  = 0;

    gradient_minimization_state result = END_OF_PATH;

    while(MAX_ITERATIONS>CPE_iteration)
      {
        find_new_Sigma_0(CPE_data_obj);

        result = find_new_alpha_1(CPE_data_obj);

        if(result == BELOW_MAX_ERROR or
           result == FOUND_THE_MINIMUM)
          break;

        CPE_iteration++;
      }

    //     if(concurrency.id()==0 and CPE_data_obj.id==0)
    //       std::cout << CPE_iteration << "\t" << MAX_ITERATIONS << "\t" << evaluate_Lambda_norm(CPE_data_obj) << "\t" << result << "\n";
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::perform_continuous_pole_expansion_threaded_2(CPE_data_type& CPE_data_obj)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    CPE_data_obj.Sigma_0      = 0;
    CPE_data_obj.grad_Sigma_0 = 0;

    for(int n=0; n<alpha_dmn_t::dmn_size(); n++)
      CPE_data_obj.alpha_vec_d[n] = 1./double(alpha_dmn_t::dmn_size());

    int MAX_ITERATIONS = CPE_data_obj.max_iterations;
    int CPE_iteration  = 0;

    while(MAX_ITERATIONS>CPE_iteration)
      {
        find_new_Sigma_0(CPE_data_obj);

        int index = find_new_alpha_2(CPE_data_obj);

        if(index == 0)
          break;

        CPE_iteration++;
      }

    //     if(concurrency.id()==0 and CPE_data_obj.id==0)
    //       std::cout << CPE_iteration << "\t" << MAX_ITERATIONS << "\t" << evaluate_Lambda_norm(CPE_data_obj) << "\n";
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::compute_gradient_Sigma_0(CPE_data_type& CPE_data_obj)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    project_from_real_axis_to_imaginary_axis(CPE_data_obj.alpha_vec_d,
                                             CPE_data_obj.f_wn_re_vec,
                                             *(CPE_data_obj.A_matrix_re_ptr));

    for(int n=w_IMAG::dmn_size()/2; n<w_IMAG::dmn_size(); n++)
      CPE_data_obj.f_wn_re_vec[n] = CPE_data_obj.F_wn_re_vec[n]-(CPE_data_obj.f_wn_re_vec[n]+CPE_data_obj.Sigma_0);

    CPE_data_obj.grad_Sigma_0=0.;
    for(int n=w_IMAG::dmn_size()/2; n<w_IMAG::dmn_size(); n++)
      CPE_data_obj.grad_Sigma_0 += -2.*CPE_data_obj.f_wn_re_vec[n];
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::find_new_Sigma_0(CPE_data_type& CPE_data_obj)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    for(int i=0; i<10; i++){
      compute_gradient_Sigma_0(CPE_data_obj);
      CPE_data_obj.Sigma_0 = CPE_data_obj.Sigma_0 - CPE_data_obj.grad_Sigma_0/scalartype(2.*w_IMAG::dmn_size()/2.);
    }
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::compute_gradient_alphas(CPE_data_type& CPE_data_obj)
  {
    //std::cout << __FUNCTION__ << "\t" << __LINE__ << "\n";

    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    //     grad_F = -2*A1*(real(spoints) - transpose(A1)*a) ...
    //               - 2*A2*(imag(spoints) - transpose(A2)*a);

    {// compute the gradient
      {// -2*A1*(real(spoints) - transpose(A1)*a)
        project_from_real_axis_to_imaginary_axis(CPE_data_obj.alpha_vec_d,
                                                 CPE_data_obj.f_wn_re_vec,
                                                 *(CPE_data_obj.A_matrix_re_ptr));

        for(int n=0; n<w_IMAG::dmn_size(); n++)
          CPE_data_obj.f_wn_re_vec[n] = CPE_data_obj.F_wn_re_vec[n]-(CPE_data_obj.f_wn_re_vec[n]+CPE_data_obj.Sigma_0);

        project_from_imaginary_axis_to_real_axis(CPE_data_obj.f_wn_re_vec,
                                                 CPE_data_obj.gradient_re,
                                                 *(CPE_data_obj.A_matrix_re_ptr));

        for(int n=0; n<alpha_dmn_t::dmn_size(); n++)
          CPE_data_obj.gradient_re[n] *= -2;
      }

      {// -2*A2*(imag(spoints) - transpose(A2)*a)
        project_from_real_axis_to_imaginary_axis(CPE_data_obj.alpha_vec_d,
                                                 CPE_data_obj.f_wn_im_vec,
                                                 *(CPE_data_obj.A_matrix_im_ptr));

        for(int n=0; n<w_IMAG::dmn_size(); n++)
          CPE_data_obj.f_wn_im_vec[n] = CPE_data_obj.F_wn_im_vec[n]-CPE_data_obj.f_wn_im_vec[n];

        project_from_imaginary_axis_to_real_axis(CPE_data_obj.f_wn_im_vec,
                                                 CPE_data_obj.gradient_im,
                                                 *(CPE_data_obj.A_matrix_im_ptr));

        for(int n=0; n<alpha_dmn_t::dmn_size(); n++)
          CPE_data_obj.gradient_im[n] *= -2;
      }
    }

    //if(normalize)
    { // normalize
      for(int n=0; n<alpha_dmn_t::dmn_size(); n++)
        CPE_data_obj.gradient[n] = -(CPE_data_obj.gradient_re[n]+CPE_data_obj.gradient_im[n]);

      double norm_gradient = 0;
      for(int n=0; n<alpha_dmn_t::dmn_size(); n++)
        norm_gradient += square(CPE_data_obj.gradient[n]);

      for(int n=0; n<alpha_dmn_t::dmn_size(); n++)
        CPE_data_obj.gradient[n] /= sqrt(norm_gradient);
    }
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  typename continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::gradient_minimization_state
  continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::find_new_alpha_1(CPE_data_type& CPE_data_obj)
  {
    profiler_type profiler( __FUNCTION__, __FILE__, __LINE__);

    gradient_minimization_state result = END_OF_PATH;

    compute_gradient_alphas(CPE_data_obj);

    double LAMBDA_MAX = 1;

    // follow the path along the gradient
    int index = 0;
    for( index=0; index<CPE_data_obj.y.size(); index++)
      {
        double lambda       = 0.;
        double delta_lambda = LAMBDA_MAX/10.;

        CPE_data_obj.x[index] = lambda+index*delta_lambda;

        for(int n=0; n<alpha_dmn_t::dmn_size(); n++)
          {
            scalartype value = CPE_data_obj.alpha_vec_d[n] + CPE_data_obj.x[index]*CPE_data_obj.gradient[n];

            CPE_data_obj.alpha_vec_z[n].real(value<0.? 0. : value);
            CPE_data_obj.alpha_vec_z[n].imag(0.);
          }

        CPE_data_obj.y[index] = evaluate_Lambda_norm(CPE_data_obj);

        if(CPE_data_obj.y[index] < CPE_data_obj.max_error)
          return BELOW_MAX_ERROR;

        if(index>1 && CPE_data_obj.y[index-1]<CPE_data_obj.y[index]){
          result = FOUND_A_MINIMUM;
          break;
        }
      }

    double lambda = 0;

    {//find the lambda at the minimum
      switch(result)
        {
        case FOUND_A_MINIMUM:
          {
            assert(index>=2 and index<CPE_data_obj.x.size());

            /*
              {// compute the minimum
              LIN_ALG::matrix<scalartype, LIN_ALG::CPU> V(3);

              V(0,0) = square(CPE_data_obj.x[index-2]); V(0,1) = CPE_data_obj.x[index-2]; V(0,2) = 1.;
              V(1,0) = square(CPE_data_obj.x[index-1]); V(1,1) = CPE_data_obj.x[index-1]; V(1,2) = 1.;
              V(2,0) = square(CPE_data_obj.x[index-0]); V(2,1) = CPE_data_obj.x[index-0]; V(2,2) = 1.;

              LIN_ALG::GEINV<LIN_ALG::CPU>::execute(V);

              scalartype a = V(0,0)*CPE_data_obj.y[index-2] + V(0,1)*CPE_data_obj.y[index-1] + V(0,2)*CPE_data_obj.y[index-0];
              scalartype b = V(1,0)*CPE_data_obj.y[index-2] + V(1,1)*CPE_data_obj.y[index-1] + V(1,2)*CPE_data_obj.y[index-0];
              scalartype c = V(2,0)*CPE_data_obj.y[index-2] + V(2,1)*CPE_data_obj.y[index-1] + V(2,2)*CPE_data_obj.y[index-0];

              assert(abs(CPE_data_obj.y[index-2] - (a*square(CPE_data_obj.x[index-2])+b*CPE_data_obj.x[index-2]+c) )<1.e-6);
              assert(abs(CPE_data_obj.y[index-1] - (a*square(CPE_data_obj.x[index-1])+b*CPE_data_obj.x[index-1]+c) )<1.e-6);
              assert(abs(CPE_data_obj.y[index-0] - (a*square(CPE_data_obj.x[index-0])+b*CPE_data_obj.x[index-0]+c) )<1.e-6);

              lambda = -b/(2*a);
              }
            */

            lambda = find_minimum_1(index, CPE_data_obj);

            if(lambda<0)
              {
                return FOUND_THE_MINIMUM;
              }
            else
              {
                assert(CPE_data_obj.x[index-2]-1.e-6<lambda and lambda<CPE_data_obj.x[index-0]+1.e-6);
              }
          }
          break;

        case END_OF_PATH:
          {
            assert(index==CPE_data_obj.x.size());

            lambda = CPE_data_obj.x[index-1];
          }
          break;

        default:
          throw std::logic_error(__FUNCTION__);
        }
    }

    /*
      {
      // compute alpha at the minimum
      for(int n=0; n<alpha_dmn_t::dmn_size(); n++){
      scalartype value            = CPE_data_obj.alpha_vec_d[n] + lambda*CPE_data_obj.gradient[n];

      real(CPE_data_obj.alpha_vec_z[n]) = value<0.? 0. : value;
      imag(CPE_data_obj.alpha_vec_z[n]) = 0.;
      }

      // smooth alpha out
      {
      int factor = CPE_data_obj.smoothing_factor;

      for(int n=0; n<alpha_dmn_t::dmn_size(); n++){

      double N                    = 0.;
      CPE_data_obj.alpha_vec_d[n] = 0 ;

      for(int l=0-factor; l<=0+factor; l++)
      {
      if(n+l>-1 && n+l<alpha_dmn_t::dmn_size())
      {
      CPE_data_obj.alpha_vec_d[n] += real(CPE_data_obj.alpha_vec_z[n+l]);
      N += 1.;
      }
      }

      CPE_data_obj.alpha_vec_d[n] /= N;
      }
      }
      }
    */

    compute_new_alpha(lambda, CPE_data_obj);

    return FOUND_A_MINIMUM;
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::compute_new_alpha(double lambda, CPE_data_type& CPE_data_obj)
  {
    { // compute alpha at the minimum
      for(int n=0; n<alpha_dmn_t::dmn_size(); n++){
        scalartype value            = CPE_data_obj.alpha_vec_d[n] + lambda*CPE_data_obj.gradient[n];

        CPE_data_obj.alpha_vec_z[n].real(value<0.? 0. : value);
        CPE_data_obj.alpha_vec_z[n].imag(0.);
      }
    }

    { // smooth alpha out
      int factor = CPE_data_obj.smoothing_factor;

      for(int n=0; n<alpha_dmn_t::dmn_size(); n++){

        double N                    = 0.;
        CPE_data_obj.alpha_vec_d[n] = 0 ;

        for(int l=0-factor; l<=0+factor; l++)
          {
            if(n+l>-1 && n+l<alpha_dmn_t::dmn_size())
              {
                CPE_data_obj.alpha_vec_d[n] += real(CPE_data_obj.alpha_vec_z[n+l]);
                N += 1.;
              }
          }

        CPE_data_obj.alpha_vec_d[n] /= N;
      }
    }
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  int continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::find_new_alpha_2(CPE_data_type& CPE_data_obj)
  {
    profiler_type profiler( __FUNCTION__, __FILE__, __LINE__);

    compute_gradient_alphas(CPE_data_obj);

    double LAMBDA_MAX = 1;

    double lambda       = 0.;
    double delta_lambda = LAMBDA_MAX/10.;

    int index = 0;
    for(int tmp=0; tmp<100; tmp++)
      {
        CPE_data_obj.x[tmp] = lambda+tmp*delta_lambda;

        for(int n=0; n<alpha_dmn_t::dmn_size(); n++){
          scalartype value            = CPE_data_obj.alpha_vec_d[n] + CPE_data_obj.x[tmp]*CPE_data_obj.gradient[n];

          real(CPE_data_obj.alpha_vec_z[n]) = value<0.? 0. : value;
          imag(CPE_data_obj.alpha_vec_z[n]) = 0.;
        }

        CPE_data_obj.y[tmp] = evaluate_Lambda_norm(CPE_data_obj);

        if(tmp>1 && CPE_data_obj.y[tmp-1]<CPE_data_obj.y[tmp])
          break;
        else
          index++;
      }

    if(index==100)
      {
        lambda = CPE_data_obj.x[99];
      }
    else
      {
        LIN_ALG::matrix<scalartype, LIN_ALG::CPU> V(3);

        V(0,0) = square(CPE_data_obj.x[index-2]); V(0,1) = CPE_data_obj.x[index-2]; V(0,2) = 1.;
        V(1,0) = square(CPE_data_obj.x[index-1]); V(1,1) = CPE_data_obj.x[index-1]; V(1,2) = 1.;
        V(2,0) = square(CPE_data_obj.x[index-0]); V(2,1) = CPE_data_obj.x[index-0]; V(2,2) = 1.;

        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(V);

        scalartype a = V(0,0)*CPE_data_obj.y[index-2] + V(0,1)*CPE_data_obj.y[index-1] + V(0,2)*CPE_data_obj.y[index-0];
        scalartype b = V(1,0)*CPE_data_obj.y[index-2] + V(1,1)*CPE_data_obj.y[index-1] + V(1,2)*CPE_data_obj.y[index-0];
        //         scalartype c = V(2,0)*CPE_data_obj.y[index-2] + V(2,1)*CPE_data_obj.y[index-1] + V(2,2)*CPE_data_obj.y[index-0];

        //         assert(abs(CPE_data_obj.y[index-2] - (a*square(CPE_data_obj.x[index-2])+b*CPE_data_obj.x[index-2]+c) )<1.e-6);
        //         assert(abs(CPE_data_obj.y[index-1] - (a*square(CPE_data_obj.x[index-1])+b*CPE_data_obj.x[index-1]+c) )<1.e-6);
        //         assert(abs(CPE_data_obj.y[index-0] - (a*square(CPE_data_obj.x[index-0])+b*CPE_data_obj.x[index-0]+c) )<1.e-6);

        lambda = -b/(2*a);

        lambda = lambda>0? lambda : 0;
      }

    {
      for(int n=0; n<alpha_dmn_t::dmn_size(); n++){
        scalartype value            = CPE_data_obj.alpha_vec_d[n];

        real(CPE_data_obj.alpha_vec_z[n]) = value<0.? 0. : value;
        imag(CPE_data_obj.alpha_vec_z[n]) = 0.;
      }

      scalartype L2_norm_0 = evaluate_Lambda_norm(CPE_data_obj);

      for(int n=0; n<alpha_dmn_t::dmn_size(); n++){
        scalartype value            = CPE_data_obj.alpha_vec_d[n] + lambda*CPE_data_obj.gradient[n];

        real(CPE_data_obj.alpha_vec_z[n]) = value<0.? 0. : value;
        imag(CPE_data_obj.alpha_vec_z[n]) = 0.;
      }

      scalartype L2_norm_lambda = evaluate_Lambda_norm(CPE_data_obj);

      if(L2_norm_lambda<L2_norm_0)
        {
          int factor = CPE_data_obj.smoothing_factor;

          for(int n=0; n<alpha_dmn_t::dmn_size(); n++){

            double N                    = 0.;
            CPE_data_obj.alpha_vec_d[n] = 0 ;

            for(int l=0-factor; l<=0+factor; l++)
              {
                if(n+l>-1 && n+l<alpha_dmn_t::dmn_size())
                  {
                    CPE_data_obj.alpha_vec_d[n] += real(CPE_data_obj.alpha_vec_z[n+l]);
                    N += 1.;
                  }
              }

            CPE_data_obj.alpha_vec_d[n] /= N;
          }
        }
      else
        {
          index=0;
        }
    }

    return index;
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  double continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::find_minimum_1(int index, CPE_data_type& CPE_data_obj)
  {
    double lambda = 0;

    scalartype x0 = CPE_data_obj.x[index-2];
    scalartype x1 = CPE_data_obj.x[index-1];
    scalartype x2 = CPE_data_obj.x[index-0];

    scalartype y0 = CPE_data_obj.y[index-2];
    scalartype y1 = CPE_data_obj.y[index-1];
    scalartype y2 = CPE_data_obj.y[index-0];

    scalartype a = (x2*(-y0 + y1) + x1*(y0 - y2) + x0*(-y1 + y2))/((x0 - x1)*(x0 - x2)*(x1 - x2));
    scalartype b = (std::pow(x2,2)*(y0 - y1) + std::pow(x0,2)*(y1 - y2) + std::pow(x1,2)*(-y0 + y2))/((x0 - x1)*(x0 - x2)*(x1 - x2));
    //     scalartype c = (x0*x2*(-x0 + x2)*y1 + std::pow(x1,2)*(x2*y0 - x0*y2) + x1*(-(std::pow(x2,2)*y0) + std::pow(x0,2)*y2))/((x0 - x1)*(x0 - x2)*(x1 - x2));

    //     assert(abs(CPE_data_obj.y[index-2] - (a*square(CPE_data_obj.x[index-2])+b*CPE_data_obj.x[index-2]+c) )<1.e-6);
    //     assert(abs(CPE_data_obj.y[index-1] - (a*square(CPE_data_obj.x[index-1])+b*CPE_data_obj.x[index-1]+c) )<1.e-6);
    //     assert(abs(CPE_data_obj.y[index-0] - (a*square(CPE_data_obj.x[index-0])+b*CPE_data_obj.x[index-0]+c) )<1.e-6);

    lambda = -b/(2*a);

    return lambda;
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  double continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::find_minimum_2(int index, CPE_data_type& CPE_data_obj)
  {
    double lambda = 0;

    LIN_ALG::matrix<scalartype, LIN_ALG::CPU> V(3);

    V(0,0) = square(CPE_data_obj.x[index-2]); V(0,1) = CPE_data_obj.x[index-2]; V(0,2) = 1.;
    V(1,0) = square(CPE_data_obj.x[index-1]); V(1,1) = CPE_data_obj.x[index-1]; V(1,2) = 1.;
    V(2,0) = square(CPE_data_obj.x[index-0]); V(2,1) = CPE_data_obj.x[index-0]; V(2,2) = 1.;

    LIN_ALG::GEINV<LIN_ALG::CPU>::execute(V);

    scalartype a = V(0,0)*CPE_data_obj.y[index-2] + V(0,1)*CPE_data_obj.y[index-1] + V(0,2)*CPE_data_obj.y[index-0];
    scalartype b = V(1,0)*CPE_data_obj.y[index-2] + V(1,1)*CPE_data_obj.y[index-1] + V(1,2)*CPE_data_obj.y[index-0];
    //     scalartype c = V(2,0)*CPE_data_obj.y[index-2] + V(2,1)*CPE_data_obj.y[index-1] + V(2,2)*CPE_data_obj.y[index-0];

    //     assert(abs(CPE_data_obj.y[index-2] - (a*square(CPE_data_obj.x[index-2])+b*CPE_data_obj.x[index-2]+c) )<1.e-6);
    //     assert(abs(CPE_data_obj.y[index-1] - (a*square(CPE_data_obj.x[index-1])+b*CPE_data_obj.x[index-1]+c) )<1.e-6);
    //     assert(abs(CPE_data_obj.y[index-0] - (a*square(CPE_data_obj.x[index-0])+b*CPE_data_obj.x[index-0]+c) )<1.e-6);

    lambda = -b/(2*a);

    return lambda;
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  double continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::evaluate_Lambda_norm(CPE_data_type& CPE_data_obj)
  {
    //std::cout << __FUNCTION__ << "\t" << __LINE__ << "\n";

    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    std::complex<scalartype> MIN_ONE (-1, 0);
    std::complex<scalartype> PLUS_ONE( 1, 0);

    for(int m=0; m<w_IMAG::dmn_size(); m++)
      CPE_data_obj.f_wn_vec[m] = CPE_data_obj.F_wn_vec[m]-CPE_data_obj.Sigma_0;

    LIN_ALG::GEMV<LIN_ALG::CPU>::execute('N', PLUS_ONE, *(CPE_data_obj.A_matrix_ptr), CPE_data_obj.alpha_vec_z, MIN_ONE, CPE_data_obj.f_wn_vec);

    scalartype L2_norm = 0;
    for(int l=0; l<w_IMAG::dmn_size(); l++)
      L2_norm += std::norm(CPE_data_obj.f_wn_vec[l]);

    return sqrt(L2_norm/(w_IMAG::dmn_size()));
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::project_from_real_axis_to_imaginary_axis(LIN_ALG::vector<scalartype, LIN_ALG::CPU>& real_values,
                                                                                                                                                          LIN_ALG::vector<scalartype, LIN_ALG::CPU>& imag_values,
                                                                                                                                                          LIN_ALG::matrix<scalartype, LIN_ALG::CPU>& matrix)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    LIN_ALG::GEMV<LIN_ALG::CPU>::execute('N', matrix, real_values, imag_values);
  }

  template<class parameters_type, class basis_function_t, typename k_dmn_t, typename w_dmn_t>
  void continuous_pole_expansion<parameters_type, basis_function_t, k_dmn_t, w_dmn_t, WEIGHTED_GRADIENT_METHOD>::project_from_imaginary_axis_to_real_axis(LIN_ALG::vector<scalartype, LIN_ALG::CPU>& imag_values,
                                                                                                                                                          LIN_ALG::vector<scalartype, LIN_ALG::CPU>& real_values,
                                                                                                                                                          LIN_ALG::matrix<scalartype, LIN_ALG::CPU>& matrix)
  {
    profiler_type profiler(__FUNCTION__, __FILE__, __LINE__);

    LIN_ALG::GEMV<LIN_ALG::CPU>::execute('T', matrix, imag_values, real_values);
  }

}

#endif

