//-*-C++-*-

#ifndef MC_SOLVER_PARAMETERS_H
#define MC_SOLVER_PARAMETERS_H

/*!
 *   \ingroup  PARAMETERS
 *
 *   \author   Peter Staar
 *   \brief    ...
 */
template<DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
class MC_solver_parameters 
{
public:

  MC_solver_parameters();
  ~MC_solver_parameters();

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

  template<class concurrency_type>
  int  get_buffer_size(const concurrency_type& concurrency) const;

  template<class concurrency_type>
  void pack           (const concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  template<class concurrency_type>
  void unpack         (const concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

  template<class stream_type>
  void to_JSON(stream_type& ss);
  
  template<class JSON_reader_type>
  void from_JSON(JSON_reader_type& reader);

  template<class read_write_type>
  void read_write(read_write_type& read_write_obj);
};

template<DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
MC_solver_parameters<CLUSTER_SOLVER_NAME>::MC_solver_parameters()
{}

template<DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
MC_solver_parameters<CLUSTER_SOLVER_NAME>::~MC_solver_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
template<class concurrency_type>
int MC_solver_parameters<CLUSTER_SOLVER_NAME>::get_buffer_size(const concurrency_type& concurrency) const
{
  return 0;
}

template<DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
template<class concurrency_type>
void MC_solver_parameters<CLUSTER_SOLVER_NAME>::pack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{}

template<DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
template<class concurrency_type>
void MC_solver_parameters<CLUSTER_SOLVER_NAME>::unpack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
template<class stream_type>
void MC_solver_parameters<CLUSTER_SOLVER_NAME>::to_JSON(stream_type& ss)
{}
  
template<DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
template<class JSON_reader_type>
void MC_solver_parameters<CLUSTER_SOLVER_NAME>::from_JSON(JSON_reader_type& reader)
{}

#endif
