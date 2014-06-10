//-*-C++-*-

#ifndef ADVANCED_EXACT_DIAGONALIZATION_CLUSTER_SOLVER_H
#define ADVANCED_EXACT_DIAGONALIZATION_CLUSTER_SOLVER_H

namespace DCA
{
  /*!
   *  \defgroup EXACT-DIAGONALIZATION
   */

  /*!
   * \class   cluster_solver
   * \ingroup EXACT-DIAGONALIZATION
   * \brief   ED cluster-solver
   * \author  Urs Haehner
   * \version 1.0
   */
  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  class cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>
  {
#include "type_definitions.h"

    typedef DCA_data   <parameters_type> MOMS_w_imag_type;
    typedef MOMS_w_real<parameters_type> MOMS_w_real_type;

    typedef ADVANCED_EXACT_DIAGONALIZATION::advanced_ed_options<parameters_type> ed_options_type;

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

  private:

    parameters_type&  parameters;
    MOMS_type&        MOMS_imag;
    MOMS_w_real_type& MOMS_real;

    ADVANCED_EXACT_DIAGONALIZATION::Fock_space                  <parameters_type, ed_options_type> Fock_obj;

    ADVANCED_EXACT_DIAGONALIZATION::fermionic_Hamiltonian       <parameters_type, ed_options_type> Ham_obj;

    ADVANCED_EXACT_DIAGONALIZATION::fermionic_overlap_matrices  <parameters_type, ed_options_type> overlap_obj;

    ADVANCED_EXACT_DIAGONALIZATION::fermionic_sp_Greens_function<parameters_type, ed_options_type> sp_Greens_function_obj;
    ADVANCED_EXACT_DIAGONALIZATION::fermionic_tp_Greens_function<parameters_type, ed_options_type> tp_Greens_function_obj;
  };


  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::cluster_solver(parameters_type&  parameters_ref,
                                                                                                   MOMS_type& MOMS_imag_ref,
                                                                                                   MOMS_w_real_type& MOMS_real_ref):
    parameters(parameters_ref),

    MOMS_imag(MOMS_imag_ref),
    MOMS_real(MOMS_real_ref),

    Fock_obj(true, true),
    Ham_obj(parameters),
    overlap_obj(parameters, Ham_obj),

    sp_Greens_function_obj(parameters, Ham_obj, overlap_obj),
    tp_Greens_function_obj(parameters, Ham_obj, overlap_obj)
  {
    cout << "Fock-space without symmetries:" << endl;
    Fock_obj.print_subspaces();

    std::vector<ADVANCED_EXACT_DIAGONALIZATION::Hilbert_space<parameters_type, ed_options_type> >&
      Hilbert_spaces = Fock_obj.get_elements();

    int HS_0 = 0;
    int HS_1 = 0;
    int b_s_r = parameters.get_nu();

    for(int i = 0; i < Hilbert_spaces.size(); ++i){
      if(Hilbert_spaces[i].get_occupation()    == parameters.get_n_0() &&
         Hilbert_spaces[i].get_magnetization() == parameters.get_Sz_0())
        HS_0 = i;

      if(Hilbert_spaces[i].get_occupation()    == parameters.get_n_1() &&
         Hilbert_spaces[i].get_magnetization() == parameters.get_Sz_1())
        HS_1 = i;
    }

    if(false)
      {
        cout << "Print Hilbert-space #" << HS_0 << endl;
        Hilbert_spaces[HS_0].print(true);
        cout << "Print Hilbert-space #" << HS_1 << endl;
        Hilbert_spaces[HS_1].print(true);
        cout << "bsr = " << b_s_r << endl;
      }

    //Fock_obj.apply_rotation_symmetry(parameters.get_symmetries(), parameters.get_ED_method());

    cout << print_time() << endl;
    Fock_obj.apply_translation_symmetry(parameters.get_ED_method());

    cout << print_time() << endl;
    cout << "Create representation" << endl;
    Fock_obj.initialize_rep();
    cout << print_time() << endl;

    cout << "Fock-space with symmetries:" << endl;
    Fock_obj.print_subspaces();

    if(parameters.do_orthogonality_check())
      cout << "subspaces orthogonal: " << Fock_obj.check_orthogonality() << endl;
  }

  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::~cluster_solver()
  {}

  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  void cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::initialize(int dca_iteration)
  {
    //     std::cout << __FUNCTION__ << "\n";

    function<std::complex<double>, dmn_3<nu,nu,r_DCA> > H_DCA;

    //     for(int l=0; l<r_DCA::dmn_size(); l++)
    //       cout << l << "\t" << MOMS_imag.H_DCA(0,0,l) << "\n";
    //     cout << "\n";

    MATH_ALGORITHMS::TRANSFORM<k_DCA, r_DCA>::execute(MOMS_imag.H_DCA, H_DCA);

    //     for(int l=0; l<r_DCA::dmn_size(); l++)
    //       cout << l << "\t" << H_DCA(0,0,l) << "\n";
    //     cout << "\n";

    //assert(false);

    Ham_obj.initialize(H_DCA, MOMS_imag.H_interactions);
  }

  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  void cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::execute()
  {
    std::cout << __FUNCTION__ << "\n";

    {// creation and annihilation matrices
      cout << print_time() << endl;
      overlap_obj.construct_creation_set_all();
      overlap_obj.construct_annihilation_set_all();
      cout << print_time() << endl;
      overlap_obj.construct_creation_set_nonzero_sparse();
      overlap_obj.construct_annihilation_set_nonzero_sparse();
      cout << print_time() << endl;
    }

    {// non-interacting Greensfunction
      cout << print_time() << endl;
      Ham_obj.construct_Hamiltonians(false);

      cout << print_time() << endl;
      Ham_obj.diagonalize_Hamiltonians_st();

      cout << print_time() << endl;
      Ham_obj.set_spectrum(MOMS_real.E0_w);

      // Ham_obj.print_Hamiltonian("data/hamiltonian_nonint.dat");
      // Ham_obj.print_eigen_energies("data/energies_nonint.dat");
      // Ham_obj.print_eigen_states("data/eigenstates_nonint.dat");

      //cout << "Check hermitianess: " << overlap_obj.check_hermitianess() << endl;

      //overlap_obj.print_creation_matrix("data/creation_nonint.dat");

      cout << print_time() << endl;
      //Greens_function_obj.compute_all_sp_functions(MOMS_imag, MOMS_real, false);
      sp_Greens_function_obj.compute_all_sp_functions_slow(MOMS_imag, MOMS_real, false);
      cout << print_time() << endl;

      cout << print_time() << endl;
      tp_Greens_function_obj.compute_two_particle_Greens_function(false);
      cout << print_time() << endl;

      tp_Greens_function_obj.compute_particle_particle_superconducting_A(MOMS_imag.G4_k_k_w_w);
    }

    {// interacting Greensfunction
      cout << print_time() << endl;
      Ham_obj.construct_Hamiltonians(true);

      cout << print_time() << endl;
      Ham_obj.diagonalize_Hamiltonians_st();

      cout << print_time() << endl;
      Ham_obj.set_spectrum(MOMS_real.E_w);

      // Ham_obj.print_Hamiltonian("data/hamiltonian_int.dat");
      // Ham_obj.print_eigen_energies("data/energies_int.dat");
      // Ham_obj.print_eigen_states("data/eigenstates_int.dat");

      //cout << "Check hermitianess: " << overlap_obj.check_hermitianess() << endl;

      //overlap_obj.print_creation_matrix("data/creation_int.dat");

      cout << print_time() << endl;
      //Greens_function_obj.compute_all_sp_functions(MOMS_imag, MOMS_real, true);
      sp_Greens_function_obj.compute_all_sp_functions_slow(MOMS_imag, MOMS_real, true);
      cout << print_time() << endl;

      cout << print_time() << endl;
      tp_Greens_function_obj.compute_two_particle_Greens_function(true);
      cout << print_time() << endl;

      tp_Greens_function_obj.compute_particle_particle_superconducting_A(MOMS_imag.G4_k_k_w_w);

    }

    cout << print_time() << endl;
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

    cout << print_time() << endl;

    sp_Greens_function_obj.compute_S_k_w(MOMS_imag.G_k_w, MOMS_imag.G0_k_w, MOMS_imag.Sigma);
    sp_Greens_function_obj.compute_S_k_w(MOMS_real.G_k_w, MOMS_real.G0_k_w, MOMS_real.Sigma);

//     {

//       cout << print_time() << endl;
//       tp_Greens_function_obj.compute_two_particle_Greens_function(true);
//       cout << print_time() << endl;

//       tp_Greens_function_obj.compute_particle_particle_superconducting(MOMS_imag.G4_k_k_w_w);
//     }

    cout << print_time() << endl;

    //assert(false);
  }

  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  template<typename dca_info_struct_t>
  void cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::finalize(dca_info_struct_t& dca_info_struct)
  {
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
  void cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::write(std::string file_name)
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
            MOMS_imag .write(writer);
            MOMS_real .write(writer);

	    cout << "\n\n\t\t start writing tp-Greens-function\n\n";
	    tp_Greens_function_obj.write(writer);

            writer.close_file();
          }
        }
        break;

      case IO::HDF5 :
        {
          IO::writer<IO::HDF5> writer;
          {
            writer.open_file(file_name);

            cout << "\n\n\t\t start writing parameters\n\n";
            parameters.write(writer);

            cout << "\n\n\t\t start writing MOMS_imag\n\n";
            MOMS_imag .write(writer);

            cout << "\n\n\t\t start writing MOMS_real\n\n";
            MOMS_real .write(writer);

	    cout << "\n\n\t\t start writing tp-Greens-function\n\n";
	    tp_Greens_function_obj.write(writer);

            writer.close_file();
          }
        }
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }
  }

}

#endif
