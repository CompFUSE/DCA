//-*-C++-*-

#ifndef MC_PTHREADS_SOLVER_H
#define MC_PTHREADS_SOLVER_H

/*!
 *   \author Peter Staar
 */
class MC_pthreads_solver
{
public:

  MC_pthreads_solver();

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
  void to_JSON(stream_type& ss, bool is_end=false);
  
  template<class JSON_reader_type>
  void from_JSON(JSON_reader_type& reader);

  template<class read_write_type>
  void read_write(read_write_type& read_write_obj);
  
private:

  int nr_walkers;
  int nr_accumulators;
  int additional_steps;

  int HTS_threads;
};

MC_pthreads_solver::MC_pthreads_solver():
  nr_walkers(1),
  nr_accumulators(1),
  additional_steps(1),

  HTS_threads(1)
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<class concurrency_type>
int MC_pthreads_solver::get_buffer_size(const concurrency_type& concurrency) const
{
  int buffer_size = 0;

  buffer_size += concurrency.getBufferSize(nr_walkers);
  buffer_size += concurrency.getBufferSize(nr_accumulators);
  buffer_size += concurrency.getBufferSize(additional_steps);

  buffer_size += concurrency.getBufferSize(HTS_threads);

  return buffer_size;
}

template<class concurrency_type>
void MC_pthreads_solver ::pack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.pack(buffer, buffer_size, position, nr_walkers);
  concurrency.pack(buffer, buffer_size, position, nr_accumulators);
  concurrency.pack(buffer, buffer_size, position, additional_steps);

  concurrency.pack(buffer, buffer_size, position, HTS_threads);
}

template<class concurrency_type>
void MC_pthreads_solver ::unpack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, nr_walkers);
  concurrency.unpack(buffer, buffer_size, position, nr_accumulators);
  concurrency.unpack(buffer, buffer_size, position, additional_steps);

  concurrency.unpack(buffer, buffer_size, position, HTS_threads);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<class stream_type>
void MC_pthreads_solver ::to_JSON(stream_type& ss, bool is_end)
{
  ss << "\"MC-posix-parameters\" :";
  ss << "\n{ \n";

  JSON_writer::write(ss, "nr-walkers"       , nr_walkers);
  JSON_writer::write(ss, "nr-accumulators"  , nr_accumulators);
  JSON_writer::write(ss, "additional-steps" , additional_steps, true);
  
  if(is_end)
    ss << "}\n";
  else
    ss << "},\n";  
}
  
template<class JSON_reader_type>
void MC_pthreads_solver ::from_JSON(JSON_reader_type& reader)
{
  try
    {
      typedef typename JSON_reader_type::JsonAccessor JsonAccessor;
      const JsonAccessor control(reader["MC-posix-parameters"]);
     
      nr_walkers                       <= control["nr-walkers"];
      nr_accumulators                  <= control["nr-accumulators"];
      additional_steps                 <= control["additional-steps"];
      //use_gpu                          <= control["use_gpu"];
    }
  catch(const std::exception& r_e) {}

  cout <<"\n\n"
       << "\t nr_walkers      --> " << nr_walkers      << "\n"
       << "\t nr_accumulators --> " << nr_accumulators << "\n"
       <<"\n\n";

  if(nr_walkers<1 || nr_accumulators<1){
    throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

template<class read_write_type>
void MC_pthreads_solver::read_write(read_write_type& read_write_obj)
{
    try
    {
      read_write_obj.open_group("MC-posix-parameters");

      try { read_write_obj.execute("nr-walkers", nr_walkers);             } catch(const std::exception& r_e) {}

      try { read_write_obj.execute("nr-accumulators", nr_accumulators);   } catch(const std::exception& r_e) {}

      try { read_write_obj.execute("additional-steps", additional_steps); } catch(const std::exception& r_e) {}

      try { read_write_obj.execute("HTS-threads", HTS_threads); } catch(const std::exception& r_e) {}

      read_write_obj.close_group();
    }
  catch(const std::exception& r_e) 
    { 
      cout << "\n\t physics-parameters defined !!  \n\n";
      throw std::logic_error(__PRETTY_FUNCTION__);
    }
}

#endif
