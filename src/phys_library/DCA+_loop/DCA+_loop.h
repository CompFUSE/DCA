//-*-C++-*-

#ifndef DCA_CALCULATION_H
#define DCA_CALCULATION_H

namespace DCA
{
  /*!
   *  \class DCA_Calculation
   *
   *  \author Peter Staar
   *  \brief  This class encapsulates the DCA-loop, using the template programming style.
   */
  template<class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
  class DCA_calculation
  {
#include "type_definitions.h"

    typedef typename parameters_type::profiler_type    profiler_t;
    typedef typename parameters_type::concurrency_type concurrency_type;

    typedef DCA::cluster_exclusion         <parameters_type, MOMS_type> cluster_exclusion_type;
    typedef DCA::double_counting_correction<parameters_type, MOMS_type> double_counting_correction_type;

    typedef DCA::coarsegraining_sp<parameters_type, k_DCA        > coarsegraining_sp_type;
    typedef DCA::lattice_map_sp   <parameters_type, k_DCA, k_HOST> lattice_map_sp_type;

    typedef DCA::update_chemical_potential<parameters_type, MOMS_type, coarsegraining_sp_type> update_chemical_potential_type;

    typedef DCA::cluster_solver<DCA::HIGH_TEMPERATURE_SERIES, LIN_ALG::CPU, parameters_type, MOMS_type> HTS_solver_type;

  public:

    DCA_calculation(parameters_type&   parameters_ref,
                    MOMS_type&         MOMS_ref,
                    concurrency_type&  concurrency_ref);

    ~DCA_calculation();

    void read();

    void write();

    void initialize();

    void execute();

    void finalize();

    //     void iterate_to_consistency();
    //     void iterate_to_consistency_2();

  private:

    void adjust_chemical_potential();

    void perform_cluster_mapping();

    void perform_cluster_mapping_self_energy();

    void perform_cluster_mapping_Greens_function();

    void adjust_coarsegrained_self_energy();

    void perform_cluster_exclusion_step();

    void solve_cluster_problem(int DCA_iteration);

    void adjust_impurity_self_energy();

    void perform_lattice_mapping();

    void perform_lattice_mapping_with_HTS();
    void perform_lattice_mapping_without_HTS();

    void update_DCA_calculation_data_functions(int DCA_iteration);

  private:

    parameters_type&            parameters;
    MOMS_type&                  MOMS;
    concurrency_type&           concurrency;

    DCA_calculation_data   DCA_info_struct;

    cluster_exclusion_type          cluster_exclusion_obj;
    double_counting_correction_type double_counting_correction_obj;

    coarsegraining_sp_type cluster_mapping_obj;
    lattice_map_sp_type    lattice_mapping_obj;

    update_chemical_potential_type update_chemical_potential_obj;

    Monte_Carlo_Integrator_type MonteCarloIntegrator;
  };

  template<class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
  DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::DCA_calculation(parameters_type&   parameters_ref,
                                                                                            MOMS_type&         MOMS_ref,
                                                                                            concurrency_type&  concurrency_ref):
    parameters(parameters_ref),
    MOMS(MOMS_ref),
    concurrency(concurrency_ref),

    DCA_info_struct(),

    cluster_exclusion_obj         (parameters, MOMS),
    double_counting_correction_obj(parameters, MOMS),

    cluster_mapping_obj(parameters),
    lattice_mapping_obj(parameters),

    update_chemical_potential_obj(parameters, MOMS, cluster_mapping_obj),

    MonteCarloIntegrator(parameters_ref, MOMS_ref)
  {
    if(concurrency.id() == concurrency.first())
      cout << "\n\n\t" << __FUNCTION__ << " has started \t" << print_time() << "\n\n";
  }

  template<class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
  DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::~DCA_calculation()
  {
    if(concurrency.id() == concurrency.first())
      cout << "\n\n\t" << __FUNCTION__ << " has finished \t" << print_time() << "\n\n";
  }

  template<class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
  void DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::read()
  {
    if(parameters.get_Sigma_file() != "zero")
      MOMS.read(parameters.get_Sigma_file());
  }

  template<class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
  void DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::write()
  {
    IO::FORMAT  FORMAT    = parameters.get_output_format();
    std::string file_name = parameters.get_directory() + parameters.get_output_file_name();

    cout << "\n\n\t\t start writing " << file_name << "\t" << print_time() << "\n\n";

    switch(FORMAT)
      {
      case IO::JSON :
        {
          IO::writer<IO::JSON> writer;
          {
            writer.open_file(file_name);

            parameters          .write(writer);
            MOMS                .write(writer);
            MonteCarloIntegrator.write(writer);
            DCA_info_struct     .write(writer);

            writer.close_file();
          }
        }
        break;

      case IO::HDF5 :
        {
          IO::writer<IO::HDF5> writer;
          {
            writer.open_file(file_name);

            parameters          .write(writer);
            MOMS                .write(writer);
            MonteCarloIntegrator.write(writer);
            DCA_info_struct     .write(writer);

            writer.close_file();
          }
        }
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }
  }

  /*
    template<class parameters_type,class MOMS_type, class Monte_Carlo_Integrator_type>
    void DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::iterate_to_consistency()
    {
    for(int i=0; i<parameters.get_DCA_iterations(); i++)
    {
    if(concurrency.id() == concurrency.first())
    cout << "\n\n\t" << __FUNCTION__ << " has started with DCA-iteration --> " << i << " \n\n";

    {
    profiler_t profiler("cluster-mapping", "DCA", __LINE__);
    concurrency << "\n\n\t coarsegrain has started";
    MOMS.coarsegrain_functions(i);
    concurrency << "\n\t coarsegrain has ended";
    }

    DCA_info_struct.chemical_potential(i) = parameters.get_chemical_potential();

    {
    profiler_t profiler("Monte Carlo Integration", "DCA", __LINE__);
    MonteCarloIntegrator.initialize(i); // load G0 from MOMS into MCI
    MonteCarloIntegrator.integrate();   // --> does the gather !!
    }

    {
    profiler_t profiler("gather and compute Self-Energy", "DCA", __LINE__);
    MonteCarloIntegrator.finalize(DCA_info_struct); // load G and sigma from MCI into MOMS
    }

    {
    profiler_t profiler("lattice-mapping", "DCA", __LINE__);
    concurrency << "\n\n\t lattice-mapping has started";
    perform_lattice_mapping();
    //perform_lattice_mapping();
    concurrency << "\n\n\t lattice-mapping has ended";
    }

    {
    update_DCA_calculation_data_functions(i);
    }
    }

    if(parameters.use_interpolated_Self_energy())
    {
    MOMS.compute_Sigma_bands();
    MOMS.compute_single_particle_properties();
    }
    }
  */

  template<class parameters_type,class MOMS_type, class Monte_Carlo_Integrator_type>
  void DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::initialize()
  {
    if(parameters.get_Sigma_file() != "zero")
      {
        MOMS.initialize_Sigma();

        perform_lattice_mapping();
      }
  }

  template<class parameters_type,class MOMS_type, class Monte_Carlo_Integrator_type>
  void DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::execute()
  {
    for(int i=0; i<parameters.get_DCA_iterations(); i++)
      {
        adjust_chemical_potential();

        perform_cluster_mapping();

        adjust_coarsegrained_self_energy(); // double-counting-correction

        perform_cluster_exclusion_step();

        solve_cluster_problem(i);

        adjust_impurity_self_energy(); // double-counting-correction

        perform_lattice_mapping();

        update_DCA_calculation_data_functions(i);
      }
  }

  template<class parameters_type,class MOMS_type, class Monte_Carlo_Integrator_type>
  void DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::finalize()
  {
    perform_cluster_mapping_self_energy();

    MOMS.compute_Sigma_bands();

    MOMS.compute_single_particle_properties();
  }

  template<class parameters_type,class MOMS_type, class Monte_Carlo_Integrator_type>
  void DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::adjust_chemical_potential()
  {
    if(parameters.adjust_chemical_potential())
      update_chemical_potential_obj.execute();
  }

  template<class parameters_type,class MOMS_type, class Monte_Carlo_Integrator_type>
  void DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::perform_cluster_mapping()
  {
    perform_cluster_mapping_self_energy();

    perform_cluster_mapping_Greens_function();

    //perform_cluster_exclusion_step();
  }

  template<class parameters_type,class MOMS_type, class Monte_Carlo_Integrator_type>
  void DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::perform_cluster_mapping_self_energy()
  {
    if(concurrency.id()==0)
      cout << "\n\t\t coarsegrain-Selfenergy " << print_time();

    profiler_t profiler("coarsegrain-Selfenergy", "DCA", __LINE__);

    if(parameters.use_interpolated_Self_energy())
      cluster_mapping_obj.compute_S_K_w(MOMS.Sigma_lattice, MOMS.Sigma_cluster);
    else
      MOMS.Sigma_cluster = MOMS.Sigma;

    MOMS.print_Sigma_QMC_versus_Sigma_cg();

    symmetrize::execute(MOMS.Sigma_cluster, MOMS.H_symmetry);
  }

  template<class parameters_type,class MOMS_type, class Monte_Carlo_Integrator_type>
  void DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::perform_cluster_mapping_Greens_function()
  {
    if(concurrency.id()==0)
      cout << "\n\t\t coarsegrain-Greens-function " << print_time();

    profiler_t profiler("coarsegrain-Greens-function", "DCA", __LINE__);

    if(parameters.use_interpolated_Self_energy())
      cluster_mapping_obj.compute_G_K_w(MOMS.H_HOST, MOMS.Sigma_lattice, MOMS.G_k_w);
    else
      cluster_mapping_obj.compute_G_K_w(MOMS.H_HOST, MOMS.Sigma        , MOMS.G_k_w);

    symmetrize::execute(MOMS.G_k_w, MOMS.H_symmetry);
  }

  template<class parameters_type,class MOMS_type, class Monte_Carlo_Integrator_type>
  void DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::adjust_coarsegrained_self_energy()
  {
    double_counting_correction_obj.execute_before_solver();
  }

  template<class parameters_type,class MOMS_type, class Monte_Carlo_Integrator_type>
  void DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::perform_cluster_exclusion_step()
  {
    if(concurrency.id()==0)
      cout << "\n\t\t cluster-exclusion-step " << print_time();

    profiler_t profiler("cluster-exclusion-step", "DCA", __LINE__);

    cluster_exclusion_obj.execute();
  }

  template<class parameters_type,class MOMS_type, class Monte_Carlo_Integrator_type>
  void DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::solve_cluster_problem(int DCA_iteration)
  {
    {
      profiler_t profiler("initialize cluster-solver", "DCA", __LINE__);
      MonteCarloIntegrator.initialize(DCA_iteration);
    }

    {
      profiler_t profiler("Quantum Monte Carlo integration", "DCA", __LINE__);
      MonteCarloIntegrator.integrate();
    }

    {
      profiler_t profiler("finalize cluster-solver", "DCA", __LINE__);
      MonteCarloIntegrator.finalize(DCA_info_struct);
    }
  }

  template<class parameters_type,class MOMS_type, class Monte_Carlo_Integrator_type>
  void DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::adjust_impurity_self_energy()
  {
    double_counting_correction_obj.execute_after_solver();
  }

  template<class parameters_type,class MOMS_type, class Monte_Carlo_Integrator_type>
  void DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::perform_lattice_mapping()
  {
    profiler_t profiler("lattice-mapping", "DCA", __LINE__);

    if(concurrency.id()==0)
      cout << "\n\t\t lattice-mapping " << print_time();

    if(parameters.use_interpolated_Self_energy())
      {
        if(parameters.use_HTS_approximation())
          {
            MOMS_type MOMS_HTS(parameters);

            MOMS_HTS.H_HOST         = MOMS.H_HOST;
            MOMS_HTS.H_interactions = MOMS.H_interactions;

            HTS_solver_type HTS_solver(parameters, MOMS_HTS);

            lattice_mapping_obj.execute_with_HTS_approximation(MOMS_HTS, HTS_solver, cluster_mapping_obj,
                                                               MOMS.Sigma, MOMS.Sigma_lattice_interpolated, MOMS.Sigma_lattice_coarsegrained, MOMS.Sigma_lattice);
          }
        else
          {
            lattice_mapping_obj.execute(MOMS.Sigma, MOMS.Sigma_lattice_interpolated, MOMS.Sigma_lattice_coarsegrained, MOMS.Sigma_lattice);
          }
      }
  }

  template<class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
  void DCA_calculation<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::update_DCA_calculation_data_functions(int i)
  {
    DCA_info_struct.density           (i) = update_chemical_potential_obj.compute_density();
    DCA_info_struct.chemical_potential(i) = parameters.get_chemical_potential();

    if(concurrency.id()==0)
      {
        cout << "\n\n\t\t\t total-density : " << DCA_info_struct.density(i) << "\t (time : " << print_time() << ")\n\n";
      }

    for(int l1=0; l1<b::dmn_size()*s::dmn_size(); l1++)
      DCA_info_struct.orbital_occupancies(l1, i) = MOMS.orbital_occupancy(l1);

    for(int l1=0; l1<b::dmn_size()*s::dmn_size(); l1++)
      for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++)
        DCA_info_struct.n_k(l1, k_ind, i) = 1.-MOMS.G_k_t(l1,l1,k_ind,0);

    for(int l1=0; l1<b::dmn_size()*s::dmn_size(); l1++)
      for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++)
        DCA_info_struct.A_k(l1, k_ind, i) = MOMS.G_k_t(l1,l1,k_ind,parameters.get_sp_time_intervals()/2)*parameters.get_beta()/M_PI;
  }
}

#endif
