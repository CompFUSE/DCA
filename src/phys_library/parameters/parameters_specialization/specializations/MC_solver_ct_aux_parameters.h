//-*-C++-*-

#ifndef MC_SOLVER_CT_AUX_PARAMETERS_H
#define MC_SOLVER_CT_AUX_PARAMETERS_H

/*!
 *  \author Peter Staar
 */
template<>
class MC_solver_parameters<DCA::CT_AUX_CLUSTER_SOLVER> 
{
public:

  MC_solver_parameters();
  ~MC_solver_parameters();

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

  template<class concurrency_type>
  int  get_buffer_size( concurrency_type& concurrency);

  template<class concurrency_type>
  void pack           ( concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  template<class concurrency_type>
  void unpack         ( concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

  template<class stream_type>
  void to_JSON(stream_type& ss, bool is_end=false);
  
  template<class JSON_reader_type>
  void from_JSON(JSON_reader_type& reader);

  template<class read_write_type>
  void read_write(read_write_type& read_write_obj);

/******************************************
 ***        DATA                        ***
 ******************************************/

  int    get_K_PHANI();
  double get_K_CT_AUX();

  int    get_initial_matrix_size();

private:

  int    submatrix_size;
  int    initial_matrix_size;

  double K_parameter;
};

MC_solver_parameters<DCA::CT_AUX_CLUSTER_SOLVER>::MC_solver_parameters():

  submatrix_size(128),
  initial_matrix_size(128),

  K_parameter(1.)
{}

MC_solver_parameters<DCA::CT_AUX_CLUSTER_SOLVER>::~MC_solver_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<class concurrency_type>
int MC_solver_parameters<DCA::CT_AUX_CLUSTER_SOLVER>::get_buffer_size( concurrency_type& concurrency)
{
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(submatrix_size);
  buffer_size += concurrency.get_buffer_size(initial_matrix_size);
  buffer_size += concurrency.get_buffer_size(K_parameter);

  return buffer_size;
}

template<class concurrency_type>
void MC_solver_parameters<DCA::CT_AUX_CLUSTER_SOLVER> ::pack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.pack(buffer, buffer_size, position, submatrix_size);
  concurrency.pack(buffer, buffer_size, position, initial_matrix_size);
  concurrency.pack(buffer, buffer_size, position, K_parameter);
}

template<class concurrency_type>
void MC_solver_parameters<DCA::CT_AUX_CLUSTER_SOLVER> ::unpack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, submatrix_size);
  concurrency.unpack(buffer, buffer_size, position, initial_matrix_size);
  concurrency.unpack(buffer, buffer_size, position, K_parameter);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<class read_write_type>
void MC_solver_parameters<DCA::CT_AUX_CLUSTER_SOLVER>::read_write(read_write_type& read_write_obj)
{
  try
    {
      read_write_obj.open_group("CT-AUX-solver");
      
      try { read_write_obj.execute("submatrix-size"     , submatrix_size);      } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("initial-matrix-size", initial_matrix_size); } catch(const std::exception& r_e) {}	
      try { read_write_obj.execute("K-parameter"        , K_parameter);         } catch(const std::exception& r_e) {}

      read_write_obj.close_group();
    }
  catch(const std::exception& r_e) 
    { 
      std::cout << "\n\t CT-AUX-solver-parameters defined !!  \n\n";
      throw std::logic_error(__PRETTY_FUNCTION__);
    }
}

/******************************************
 ***        DATA                        ***
 ******************************************/

int MC_solver_parameters<DCA::CT_AUX_CLUSTER_SOLVER>::get_K_PHANI()
{
  return submatrix_size;
}

int MC_solver_parameters<DCA::CT_AUX_CLUSTER_SOLVER>::get_initial_matrix_size()
{
  return initial_matrix_size;
}

double MC_solver_parameters<DCA::CT_AUX_CLUSTER_SOLVER>::get_K_CT_AUX()
{
  return K_parameter;
}

#endif
