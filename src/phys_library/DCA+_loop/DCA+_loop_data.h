//-*-C++-*-

#ifndef DCA_CALCULATION_DATA_H
#define DCA_CALCULATION_DATA_H
#include"phys_library/domain_types.hpp"
using namespace types;

namespace DCA
{
  /*! 
   *  \author Peter Staar
   *  \brief  
   */
  class DCA_calculation_data
  {

    typedef dmn_0<dmn<32, int> > expansion_dmn_t;

  public:

    DCA_calculation_data();
    ~DCA_calculation_data();

    template<IO::FORMAT DATA_FORMAT>
    void read(IO::reader<DATA_FORMAT>& reader);

    template<IO::FORMAT DATA_FORMAT>
    void write(IO::writer<DATA_FORMAT>& reader);


  public:

    FUNC_LIB::function<double, DCA_iteration_domain_type> Gflop_per_mpi_task;
    FUNC_LIB::function<double, DCA_iteration_domain_type> times_per_mpi_task;

    FUNC_LIB::function<double, DCA_iteration_domain_type> thermalization_per_mpi_task;
    FUNC_LIB::function<double, DCA_iteration_domain_type> MC_integration_per_mpi_task;

    FUNC_LIB::function<double, DCA_iteration_domain_type>     Gflops_per_mpi_task;
    FUNC_LIB::function<double, DCA_iteration_domain_type> max_Gflops_per_mpi_task;

    FUNC_LIB::function<double, DCA_iteration_domain_type> sign;

    FUNC_LIB::function<double, DCA_iteration_domain_type> L2_Sigma_difference;

    FUNC_LIB::function<double, dmn_3<nu, k_DCA, DCA_iteration_domain_type> > Sigma_zero_moment;
    FUNC_LIB::function<double, dmn_3<nu, k_DCA, DCA_iteration_domain_type> > standard_deviation;

    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, expansion_dmn_t, DCA_iteration_domain_type> > sigma_lambda;
    
    FUNC_LIB::function<double, dmn_3<nu, k_DCA, DCA_iteration_domain_type> > n_k;
    FUNC_LIB::function<double, dmn_3<nu, k_DCA, DCA_iteration_domain_type> > A_k;
    FUNC_LIB::function<double, dmn_2<nu,        DCA_iteration_domain_type> > orbital_occupancies;

    FUNC_LIB::function<double, DCA_iteration_domain_type> density;
    FUNC_LIB::function<double, DCA_iteration_domain_type> chemical_potential;
    FUNC_LIB::function<double, DCA_iteration_domain_type> average_expansion_order;
  };

  DCA_calculation_data::DCA_calculation_data():
    Gflop_per_mpi_task ("Gflop_per_mpi_task"),
    times_per_mpi_task ("times_per_mpi_task"),

    thermalization_per_mpi_task("thermalization_per_mpi_task"),
    MC_integration_per_mpi_task("MC_integration_per_mpi_task"),

    Gflops_per_mpi_task("Gflops_per_mpi_task"),
    max_Gflops_per_mpi_task("max_Gflops_per_mpi_task"),

    sign                   ("sign"),
    L2_Sigma_difference    ("L2_Sigma_difference"),

    Sigma_zero_moment      ("Sigma_zero_moment"),
    standard_deviation     ("standard_deviation"),

    sigma_lambda           ("sigma_lambda"),

    n_k                    ("n_k"),
    A_k                    ("A_k"),
    orbital_occupancies    ("orbital-occupancies"),

    density                ("density"),
    chemical_potential     ("chemical-potential"),
    average_expansion_order("expansion_order")
  {}

  DCA_calculation_data::~DCA_calculation_data()
  {}

  template<IO::FORMAT DATA_FORMAT>
  void DCA_calculation_data::write(IO::writer<DATA_FORMAT>& writer)
  {
    writer.open_group("DCA-loop-functions");
      
    {
      writer.execute(Gflop_per_mpi_task);    
      writer.execute(times_per_mpi_task);    
      writer.execute(Gflops_per_mpi_task);
      
      writer.execute(sign               );
    
      writer.execute(L2_Sigma_difference);    
      writer.execute(standard_deviation);    
      
      writer.execute(chemical_potential );    
      writer.execute(density            );    
      writer.execute(average_expansion_order);    

      writer.execute(Sigma_zero_moment);    

      writer.execute(n_k);
      writer.execute(A_k);
      writer.execute(orbital_occupancies);
    }

    writer.close_group();
  }


}

#endif

