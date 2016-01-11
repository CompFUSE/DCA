//-*-C++-*-

#ifndef DCA_HIGH_TEMPERATURE_SERIES_SOLVER_INTEGRATOR_H
#define DCA_HIGH_TEMPERATURE_SERIES_SOLVER_INTEGRATOR_H

namespace DCA
{
  /*! 
   * \class   cluster_solver<HIGH_TEMPERATURE_SERIES_SOLVER>
   * \ingroup CLUSTER-SOLVER
   * \brief   high temperature series expansion solver
   * \author  Peter Staar
   * \version 1.0
   */
  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  class cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type>
  {
#include "type_definitions.h" 

  public:
    
    cluster_solver(parameters_type&   parameters_ref,
		   MOMS_type&         MOMS_ref);

    ~cluster_solver();

    void initialize(); 
    void initialize(int dca_iteration); 

    void execute();

    void finalize(); 

    template<typename dca_info_struct_t>
    void finalize(dca_info_struct_t& dca_info_struct); 

    void read(std::string filename);
    
    void write(std::string filename);

  private:

    template<IO::FORMAT DATA_FORMAT>
    void write(IO::writer<DATA_FORMAT>& writer);

  private:

    parameters_type& parameters;
    MOMS_type&       MOMS;

    SERIES_EXPANSION::series_expansion<parameters_type, MOMS_type> series_exp_obj;

    int DCA_it;
  };

  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type>::cluster_solver(parameters_type& parameters_ref,
												     MOMS_type&       MOMS_ref):
    parameters(parameters_ref),
    MOMS(MOMS_ref),

    series_exp_obj(parameters, MOMS)
  {}
  
  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type>::~cluster_solver()
  {}

  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  void cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type>::initialize()
  {}
  
  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  void cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type>::initialize(int dca_iteration)
  {
    DCA_it = dca_iteration;
  }

  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  void cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type>::execute()
  {    
    series_exp_obj.execute(false);

    //for(int i=0; i<MOMS.Sigma.size(); ++i)
    //MOMS.Sigma(i) = series_exp_obj.get_Sigma()(i);
    
    for(int i=0; i<MOMS.Sigma_lattice.size(); ++i)
      MOMS.Sigma_lattice(i) = series_exp_obj.get_Sigma()(i);

    MOMS.compute_Sigma_bands();
  }

  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  void cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type>::finalize()
  {
  }

  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  template<typename dca_info_struct_t>
  void cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type>::finalize(dca_info_struct_t& /*dca_info_struct*/)
  {
    finalize();
  }

  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  void cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type>::write(std::string file_name)
  {
    IO::FORMAT  FORMAT    = parameters.get_output_format();
    
    std::cout << "\n\n\t\t start writing " << file_name << "\n\n";
    
    switch(FORMAT)
      {
      case IO::JSON : 
	{
	  IO::writer<IO::JSON> writer;
	  {
	    writer.open_file(file_name);
	    
	    parameters.write(writer);
	    MOMS      .write(writer);
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
	    MOMS      .write(writer);
	    this->     write(writer);
	    
	    writer.close_file();
	  }
	}
	break;
	
      default:
	throw std::logic_error(__FUNCTION__);
      }
  }

  template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  template<IO::FORMAT DATA_FORMAT>
  void cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type>::write(IO::writer<DATA_FORMAT>& /*writer*/)
  {
    //writer.open_group("functions");
    
    //series_exp_obj.write(writer);
  
    //writer.close_group();
  }

}

#endif
