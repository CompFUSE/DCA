//-*-C++-*-

#ifndef ADVANCED_EXACT_DIAGONALIZATION_CLUSTER_SOLVER_H
#define ADVANCED_EXACT_DIAGONALIZATION_CLUSTER_SOLVER_H

namespace DCA {
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
template <LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
class cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type,
                     MOMS_type> {
#include "type_definitions.h"

  typedef DCA_data<parameters_type> MOMS_w_imag_type;
  typedef MOMS_w_real<parameters_type> MOMS_w_real_type;

  typedef ADVANCED_EXACT_DIAGONALIZATION::advanced_ed_options<parameters_type>
      ed_options_type;

 public:
  cluster_solver(parameters_type& parameters_ref, MOMS_type& MOMS_ref,
                 MOMS_w_real_type& MOMS_real_ref);

  ~cluster_solver();

  void initialize(int dca_iteration);

  void execute();

  template <typename dca_info_struct_t>
  void finalize(dca_info_struct_t& dca_info_struct);

  void write(std::string file_name);

 private:
  parameters_type& parameters;
  typename parameters_type::concurrency_type& concurrency;
  MOMS_type& MOMS_imag;
  MOMS_w_real_type& MOMS_real;

  ADVANCED_EXACT_DIAGONALIZATION::Fock_space<parameters_type, ed_options_type>
      Fock_obj;

  ADVANCED_EXACT_DIAGONALIZATION::fermionic_Hamiltonian<
      parameters_type, ed_options_type> Ham_obj;

  ADVANCED_EXACT_DIAGONALIZATION::fermionic_overlap_matrices<
      parameters_type, ed_options_type> overlap_obj;

  ADVANCED_EXACT_DIAGONALIZATION::fermionic_sp_Greens_function<
      parameters_type, ed_options_type> sp_Greens_function_obj;
  ADVANCED_EXACT_DIAGONALIZATION::fermionic_tp_Greens_function<
      parameters_type, ed_options_type> tp_Greens_function_obj;
};

template <LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type,
               MOMS_type>::cluster_solver(parameters_type& parameters_ref,
                                          MOMS_type& MOMS_imag_ref,
                                          MOMS_w_real_type& MOMS_real_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      MOMS_imag(MOMS_imag_ref),
      MOMS_real(MOMS_real_ref),

      Fock_obj(true, true),
      Ham_obj(parameters),
      overlap_obj(parameters, Ham_obj),

      sp_Greens_function_obj(parameters, Ham_obj, overlap_obj),
      tp_Greens_function_obj(parameters, Ham_obj, overlap_obj) {
#ifdef DEBUG
  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n\nFock-space without symmetries:" << std::endl;
    Fock_obj.print_subspaces();
  }
#endif
  
  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n\n\n"
              << "Apply translation symmetry ..." << std::endl;
  }
  Fock_obj.apply_translation_symmetry(parameters.get_ED_method());

  if (concurrency.id() == concurrency.first()) {
    std::cout << print_time() << std::endl;
    std::cout << "\nCreate representation ..." << std::endl;
  }
  Fock_obj.initialize_rep();

  if (concurrency.id() == concurrency.first()) {
    std::cout << print_time() << std::endl;
  }

#ifdef DEBUG
  if (concurrency.id() == concurrency.first()) {
    std::cout << "\nFock-space with symmetries:" << std::endl;
    Fock_obj.print_subspaces();
  }
#endif

  if (parameters.do_orthogonality_check())
    std::cout << "subspaces orthogonal: " << Fock_obj.check_orthogonality()
              << std::endl;
}

template <LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type,
               MOMS_type>::~cluster_solver() {}

template <LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
void cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type,
                    MOMS_type>::initialize(int /*dca_iteration*/)
{
  
  FUNC_LIB::function<std::complex<double>, dmn_3<nu, nu, r_DCA> > H_DCA;

  MATH_ALGORITHMS::TRANSFORM<k_DCA, r_DCA>::execute(MOMS_imag.H_DCA, H_DCA);

  Ham_obj.initialize(H_DCA, MOMS_imag.H_interactions);
}

template <LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
void cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type,
                    MOMS_type>::execute() {
  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n" << __FUNCTION__ << "\n" << std::endl;
  }

    // creation and annihilation matrices
    overlap_obj.construct_creation_set_all();
    overlap_obj.construct_annihilation_set_all();
    if (concurrency.id() == concurrency.first()) {
      std::cout << "\n" << print_time() << std::endl;
    }
    overlap_obj.construct_creation_set_nonzero_sparse();
    overlap_obj.construct_annihilation_set_nonzero_sparse();
  

  {  // non-interacting Greensfunction
    if (concurrency.id() == concurrency.first()) {
      std::cout << "\n" << print_time() << std::endl;
    }
    Ham_obj.construct_Hamiltonians(false);

    if (concurrency.id() == concurrency.first()) {
      std::cout << "\n" << print_time() << std::endl;
    }
    Ham_obj.diagonalize_Hamiltonians_st();

    if (concurrency.id() == concurrency.first()) {
      std::cout << print_time() << std::endl;
    }
    Ham_obj.set_spectrum(MOMS_real.E0_w);

    sp_Greens_function_obj.compute_all_sp_functions_slow(MOMS_imag, MOMS_real,
                                                         false);

    if (parameters.get_vertex_measurement_type() != NONE) {
      if (concurrency.id() == concurrency.first()) {
        std::cout << "\n" << print_time() << "\n" << std::endl;
      }
      tp_Greens_function_obj.compute_two_particle_Greens_function(false);
      if (concurrency.id() == concurrency.first()) {
        std::cout << "\n" << print_time() << "\n" << std::endl;
      }
      tp_Greens_function_obj.compute_particle_particle_superconducting_A(
          MOMS_imag.G4_k_k_w_w);
    }
  }

    // interacting Greensfunction
    if (concurrency.id() == concurrency.first()) {
      std::cout << "\n" << print_time() << std::endl;
    }
    Ham_obj.construct_Hamiltonians(true);

    if (concurrency.id() == concurrency.first()) {
      std::cout << "\n" << print_time() << std::endl;
    }
    Ham_obj.diagonalize_Hamiltonians_st();

    if (concurrency.id() == concurrency.first()) {
      std::cout << print_time() << std::endl;
    }
    Ham_obj.set_spectrum(MOMS_real.E_w);

    sp_Greens_function_obj.compute_all_sp_functions_slow(MOMS_imag, MOMS_real,
                                                         true);

    if (parameters.get_vertex_measurement_type() != NONE) {
      if (concurrency.id() == concurrency.first()) {
        std::cout << "\n" << print_time() << "\n" << std::endl;
      }
      tp_Greens_function_obj.compute_two_particle_Greens_function(true);
      if (concurrency.id() == concurrency.first()) {
        std::cout << "\n" << print_time() << "\n" << std::endl;
      }
      tp_Greens_function_obj.compute_particle_particle_superconducting_A(
          MOMS_imag.G4_k_k_w_w);
    }

  
    MOMS_real.A_w = 0;
    for (int l = 0; l < w_REAL::dmn_size(); l++)
      for (int j = 0; j < k_DCA::dmn_size(); j++)
        for (int i = 0; i < 2 * b::dmn_size(); i++)
          MOMS_real.A_w(l) -= imag(MOMS_real.G_k_w(i, i, j, l));

    MOMS_real.A0_w = 0;
    for (int l = 0; l < w_REAL::dmn_size(); l++)
      for (int j = 0; j < k_DCA::dmn_size(); j++)
        for (int i = 0; i < 2 * b::dmn_size(); i++)
          MOMS_real.A0_w(l) -= imag(MOMS_real.G0_k_w(i, i, j, l));

    MOMS_real.A_w *= 1. / double(M_PI * k_DCA::dmn_size() * 2 * b::dmn_size());
    MOMS_real.A0_w *= 1. / double(M_PI * k_DCA::dmn_size() * 2 * b::dmn_size());
  

  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n" << print_time() << std::endl;
  }

  sp_Greens_function_obj.compute_S_k_w(MOMS_imag.G_k_w, MOMS_imag.G0_k_w,
                                       MOMS_imag.Sigma);
  sp_Greens_function_obj.compute_S_k_w(MOMS_real.G_k_w, MOMS_real.G0_k_w,
                                       MOMS_real.Sigma);

  if (concurrency.id() == concurrency.first()) {
    std::cout << print_time() << std::endl;
  }
}

template <LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
template <typename dca_info_struct_t>
void cluster_solver<
    ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type,
    MOMS_type>::finalize(dca_info_struct_t& /*dca_info_struct*/) {
  for (int l = 0; l < MOMS_imag.G_r_w.size(); l++)
    MOMS_imag.Sigma_cluster(l) = MOMS_imag.Sigma(l);

  for (int l = 0; l < MOMS_imag.G0_k_w.size(); l++)
    MOMS_imag.G0_k_w_cluster_excluded(l) = MOMS_imag.G0_k_w(l);

  for (int l = 0; l < MOMS_imag.G0_r_w.size(); l++)
    MOMS_imag.G0_r_w_cluster_excluded(l) = MOMS_imag.G0_r_w(l);

  for (int l = 0; l < MOMS_imag.G0_k_t.size(); l++)
    MOMS_imag.G0_k_t_cluster_excluded(l) = MOMS_imag.G0_k_t(l);

  for (int l = 0; l < MOMS_imag.G0_r_t.size(); l++)
    MOMS_imag.G0_r_t_cluster_excluded(l) = MOMS_imag.G0_r_t(l);
}

template <LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
void cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type,
                    MOMS_type>::write(std::string file_name) {
  IO::FORMAT FORMAT = parameters.get_output_format();

  std::cout << "\n\n\t\t start writing " << file_name << "\n\n";

  switch (FORMAT) {
    case IO::JSON: {
      IO::writer<IO::JSON> writer;
      {
        writer.open_file(file_name);

        parameters.write(writer);
        MOMS_imag.write(writer);
        MOMS_real.write(writer);

        if (parameters.get_vertex_measurement_type() != NONE) {
          std::cout << "\n\n\t\t start writing tp-Greens-function\n\n";
          tp_Greens_function_obj.write(writer);
        }

        writer.close_file();
      }
    } break;

    case IO::HDF5: {
      IO::writer<IO::HDF5> writer;
      {
        writer.open_file(file_name);

        std::cout << "\n\n\t\t start writing parameters\n\n";
        parameters.write(writer);

        std::cout << "\n\n\t\t start writing MOMS_imag\n\n";
        MOMS_imag.write(writer);

        std::cout << "\n\n\t\t start writing MOMS_real\n\n";
        MOMS_real.write(writer);

        if (parameters.get_vertex_measurement_type() != NONE) {
          std::cout << "\n\n\t\t start writing tp-Greens-function\n\n";
          tp_Greens_function_obj.write(writer);
        }

        writer.close_file();
      }
    } break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}
}

#endif
