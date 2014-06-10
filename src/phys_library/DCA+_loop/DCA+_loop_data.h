//-*-C++-*-

#ifndef DCA_CALCULATION_DATA_H
#define DCA_CALCULATION_DATA_H

namespace DCA
{
  /*! 
   *  \author Peter Staar
   *  \brief  
   */
  class DCA_calculation_data
  {
#include "type_definitions.h"

    typedef dmn_0<dmn<32, int> > expansion_dmn_t;

  public:

    DCA_calculation_data();
    ~DCA_calculation_data();

    template<IO::FORMAT DATA_FORMAT>
    void read(IO::reader<DATA_FORMAT>& reader);

    template<IO::FORMAT DATA_FORMAT>
    void write(IO::writer<DATA_FORMAT>& reader);

//     template<class stream_type>
//     void to_JSON(stream_type& ss);

  public:

    function<double, DCA_iteration_domain_type> Gflop_per_mpi_task;
    function<double, DCA_iteration_domain_type> times_per_mpi_task;

    function<double, DCA_iteration_domain_type> thermalization_per_mpi_task;
    function<double, DCA_iteration_domain_type> MC_integration_per_mpi_task;

    function<double, DCA_iteration_domain_type>     Gflops_per_mpi_task;
    function<double, DCA_iteration_domain_type> max_Gflops_per_mpi_task;

    function<double, DCA_iteration_domain_type> sign;

    function<double, DCA_iteration_domain_type> L2_Sigma_difference;

    function<double, dmn_3<nu, k_DCA, DCA_iteration_domain_type> > Sigma_zero_moment;
    function<double, dmn_3<nu, k_DCA, DCA_iteration_domain_type> > standard_deviation;

    function<std::complex<double>, dmn_4<nu, nu, expansion_dmn_t, DCA_iteration_domain_type> > sigma_lambda;
    
    function<double, dmn_3<nu, k_DCA, DCA_iteration_domain_type> > n_k;
    function<double, dmn_3<nu, k_DCA, DCA_iteration_domain_type> > A_k;
    function<double, dmn_2<nu,        DCA_iteration_domain_type> > orbital_occupancies;

    function<double, DCA_iteration_domain_type> density;
    function<double, DCA_iteration_domain_type> chemical_potential;
    function<double, DCA_iteration_domain_type> average_expansion_order;
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

  /*
  template<class stream_type>
  void DCA_calculation_data::to_JSON(stream_type& ss)
  {
    //ss << ",";

    Gflop_per_mpi_task.to_JSON(ss);
    ss << ",";

    times_per_mpi_task.to_JSON(ss);
    ss << ",";

    Gflops_per_mpi_task.to_JSON(ss);
    ss << ",";
    
    max_Gflops_per_mpi_task.to_JSON(ss);
    ss << ",";

    sign               .to_JSON(ss);
    ss << ",";

    L2_Sigma_difference.to_JSON(ss);
    ss << ",";
    
    standard_deviation.to_JSON(ss);
    ss << ",";

    chemical_potential .to_JSON(ss);
    ss << ",";

    density            .to_JSON(ss);
    ss << ",";

    average_expansion_order.to_JSON(ss);
    ss << ",";

    Sigma_zero_moment.to_JSON(ss);
    ss << ",";

    orbital_occupancies.to_JSON(ss);
    ss << ",";

    sigma_lambda.to_JSON(ss);
    ss << ",";

    n_k.to_JSON(ss);
    ss << ",";

    A_k.to_JSON(ss);
  }
  */

}

#endif

