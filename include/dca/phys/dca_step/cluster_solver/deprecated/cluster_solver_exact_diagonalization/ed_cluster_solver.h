//-*-C++-*-

#ifndef EXACT_DIAGONALIZATION_CLUSTER_SOLVER_H
#define EXACT_DIAGONALIZATION_CLUSTER_SOLVER_H

namespace DCA
{
  /*!
   *  \defgroup EXACT-DIAGONALIZATION
   */

  /*!
   * \class   cluster_solver
   * \ingroup EXACT-DIAGONALIZATION
   * \brief   ED cluster-solver
   * \author  Peter Staar
   * \version 1.0
   */
  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  class cluster_solver<ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>
  {
#include "type_definitions.h"

    typedef DCA_data   <parameters_type> MOMS_w_imag_type;
    typedef MOMS_w_real<parameters_type> MOMS_w_real_type;

  public:

    cluster_solver(parameters_type&   parameters_ref,
                   MOMS_type&         MOMS_ref,
		   MOMS_w_real_type& MOMS_real_ref);

    ~cluster_solver();

    void initialize(int dca_iteration);

    void execute();

    template<typename dca_info_struct_t>
    void finalize(dca_info_struct_t& dca_info_struct);

    void write(std::string file_name);

    template<class stream_type>
    void to_JSON(stream_type& ss);

  private:

    parameters_type&  parameters;
    MOMS_type&        MOMS_imag;
    MOMS_w_real_type& MOMS_real;

    EXACT_DIAGONALIZATION::fermionic_Fock_space <parameters_type, b, s, r_DCA> Fock_obj;
    EXACT_DIAGONALIZATION::fermionic_Hamiltonian<parameters_type, b, s, r_DCA> Ham_obj;

    EXACT_DIAGONALIZATION::fermionic_sp_Greens_function<parameters_type, b, s, r_DCA> Greens_function_obj;
  };

  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  cluster_solver<ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::cluster_solver(parameters_type&  parameters_ref,
                                                                                          MOMS_type&        MOMS_imag_ref,
											  MOMS_w_real_type& MOMS_real_ref):
    parameters(parameters_ref),

    MOMS_imag(MOMS_imag_ref),
    MOMS_real(MOMS_real_ref),

    Fock_obj(parameters),
    Ham_obj (parameters, Fock_obj),

    Greens_function_obj(parameters, Fock_obj, Ham_obj)
  {
    //G_k_w_real.print_fingerprint();
  }

  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  cluster_solver<ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::~cluster_solver()
  {}

  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  void cluster_solver<ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::initialize(int dca_iteration)
  {
    FUNC_LIB::function<std::complex<double>, dmn_3<nu,nu,r_DCA> > H_DCA;
    //FT<k_DCA, r_DCA>::execute(MOMS_imag.H_DCA, H_DCA);
    MATH_ALGORITHMS::TRANSFORM<k_DCA, r_DCA>::execute(MOMS_imag.H_DCA, H_DCA);

    Ham_obj.initialize(H_DCA, MOMS_imag.H_interactions);
  }

  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  void cluster_solver<ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::execute()
  {
    {// non-interacting Greensfunction
      Ham_obj.construct_Hamiltonians(false);

      Ham_obj.diagonolize_Hamiltonians_mt();

      Ham_obj.set_spectrum(MOMS_real.E0_w);

      //Ham_obj.compute_G_k_w(false);
      //Ham_obj.compute_G_k_t(false);

      //Ham_obj.compute_G_k_w_and_G_k_t(false);

      Greens_function_obj.compute_all_sp_functions(MOMS_imag, MOMS_real, false);

      //       Ham_obj            .print_G_k_w();
      //       Greens_function_obj.print_G_k_w(MOMS, false);

      //       throw std::logic_error(__FUNCTION__);
    }

    {// interacting Greensfunction
      Ham_obj.construct_Hamiltonians(true);

      Ham_obj.diagonolize_Hamiltonians_mt();

      Ham_obj.set_spectrum(MOMS_real.E_w);

      //Ham_obj.compute_G_k_w(true);
      //Ham_obj.compute_G_k_t(true);

      //Ham_obj.compute_G_k_w_and_G_k_t(true);
      Greens_function_obj.compute_all_sp_functions(MOMS_imag, MOMS_real, true);

      //       Greens_function_obj.print_G_k_w(MOMS, true);
    }

    {
      MOMS_real.A_w = 0;
      for(int l=0; l<w_REAL::dmn_size(); l++)
        for(int j=0; j<k_DCA::dmn_size(); j++)
          for(int i=0; i<2*b::dmn_size(); i++)
            MOMS_real.A_w(l) -= imag(MOMS_real.G_k_w(i,i,j,l));

      MOMS_real.A0_w = 0;
      for(int l=0; l<w_REAL::dmn_size(); l++)
        for(int j=0; j<k_DCA::dmn_size(); j++)
          for(int i=0; i<2*b::dmn_size(); i++)
            MOMS_real.A0_w(l) -= imag(MOMS_real.G0_k_w(i,i,j,l));

      MOMS_real.A_w  *= 1./double(M_PI*k_DCA::dmn_size()*2*b::dmn_size());
      MOMS_real.A0_w *= 1./double(M_PI*k_DCA::dmn_size()*2*b::dmn_size());
    }

    //Ham_obj.compute_selfenergy();

    //Ham_obj.set_functions(MOMS);
    Greens_function_obj.compute_S_k_w(MOMS_imag.G_k_w, MOMS_imag.G0_k_w, MOMS_imag.Sigma);
    Greens_function_obj.compute_S_k_w(MOMS_real.G_k_w, MOMS_real.G0_k_w, MOMS_real.Sigma);
  }


  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  template<typename dca_info_struct_t>
  void cluster_solver<ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::finalize(dca_info_struct_t& dca_info_struct)
  {
    //     Ham_obj.set_functions(MOMS_imag);

    for(int l=0; l<MOMS_imag.G_r_w.size(); l++)
      MOMS_imag.Sigma_cluster(l) = MOMS_imag.Sigma(l);

    for(int l=0; l<MOMS_imag.G0_k_w.size(); l++)
      MOMS_imag.G0_k_w_cluster_excluded(l) = MOMS_imag.G0_k_w(l);

    for(int l=0; l<MOMS_imag.G0_r_w.size(); l++)
      MOMS_imag.G0_r_w_cluster_excluded(l) = MOMS_imag.G0_r_w(l);

    for(int l=0; l<MOMS_imag.G0_k_t.size(); l++)
      MOMS_imag.G0_k_t_cluster_excluded(l) = MOMS_imag.G0_k_t(l);

    for(int l=0; l<MOMS_imag.G0_r_t.size(); l++)
      MOMS_imag.G0_r_t_cluster_excluded(l) = MOMS_imag.G0_r_t(l);
  }

  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  void cluster_solver<ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::write(std::string file_name)
  {
    IO::FORMAT  FORMAT    = parameters.get_output_format();

    cout << "\n\n\t\t start writing " << file_name << "\n\n";

    switch(FORMAT)
      {
      case IO::JSON :
        {
          IO::writer<IO::JSON> writer;
          {
            writer.open_file(file_name);

            parameters.write(writer);
            MOMS_imag .write(writer);
	    MOMS_real .write(writer);

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
            MOMS_imag .write(writer);
	    MOMS_real .write(writer);

            writer.close_file();
          }
        }
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }
  }


  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  template<class stream_type>
  void cluster_solver<ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::to_JSON(stream_type& ss)
  {
    MOMS_real.A_w .to_JSON(ss);
    ss << ",\n";

    MOMS_real.A0_w.to_JSON(ss);
    ss << ",\n";

    MOMS_real.E_w .to_JSON(ss);
    ss << ",\n";

    MOMS_real.E0_w.to_JSON(ss);
    ss << ",\n";

    MOMS_real.G_r_w.to_JSON(ss);
    ss << ",\n";

    MOMS_real.G_k_w.to_JSON(ss);
    ss << ",\n";

    MOMS_real.G0_r_w.to_JSON(ss);
    ss << ",\n";

    MOMS_real.G0_k_w.to_JSON(ss);
    ss << ",\n";

    MOMS_real.Sigma.to_JSON(ss);
  }

}

#endif
