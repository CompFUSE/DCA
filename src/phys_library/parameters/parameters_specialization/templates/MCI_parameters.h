//-*-C++-*-

#ifndef MCI_PARAMETERS_H
#define MCI_PARAMETERS_H

/*!
 *   \ingroup  PARAMETERS
 *
 *   \author   Peter Staar
 *   \brief    ...
 */
class MCI_parameters 
{
public:

  MCI_parameters();
  ~MCI_parameters();

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

  template<class read_write_type>
  void read_write(read_write_type& read_write_obj);

/******************************************
 ***        DATA                        ***
 ******************************************/

  std::string get_Sigma_file();

  int    get_warm_up_sweeps();
  double get_number_of_sweeps_per_measurement();
  int    get_number_of_measurements();

  bool   adjust_self_energy_for_double_counting();


  int    get_nr_walkers();
  int    get_nr_accumulators();
  int    get_additional_steps();

  int    get_nr_HTS_threads();

private:

  std::string  Sigma_file;

  int     warm_up_sweeps;
  double  number_of_sweeps_per_measurement;
  int     measurements;

  std::string do_adaptive_double_counting;


  int    nr_walkers;
  int    nr_accumulators;
  int    additional_steps;

  int    nr_HTS_threads;
};


MCI_parameters::MCI_parameters():

  Sigma_file("zero"),

  warm_up_sweeps                  (20),
  number_of_sweeps_per_measurement(1.),
  measurements                    (100),

  do_adaptive_double_counting("false"),

  nr_walkers(1),
  nr_accumulators(1),
  additional_steps(1),

  nr_HTS_threads(1)
{}

MCI_parameters::~MCI_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<class concurrency_type>
int MCI_parameters::get_buffer_size( concurrency_type& concurrency)
{
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(Sigma_file);

  buffer_size += concurrency.get_buffer_size(warm_up_sweeps);
  buffer_size += concurrency.get_buffer_size(number_of_sweeps_per_measurement);
  buffer_size += concurrency.get_buffer_size(measurements);

  buffer_size += concurrency.get_buffer_size(do_adaptive_double_counting);

  buffer_size += concurrency.get_buffer_size(nr_walkers);
  buffer_size += concurrency.get_buffer_size(nr_accumulators);
  buffer_size += concurrency.get_buffer_size(additional_steps);

  buffer_size += concurrency.get_buffer_size(nr_HTS_threads);

  return buffer_size;
}

template<class concurrency_type>
void MCI_parameters::pack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.pack(buffer, buffer_size, position, Sigma_file);

  
  concurrency.pack(buffer, buffer_size, position, warm_up_sweeps);
  concurrency.pack(buffer, buffer_size, position, number_of_sweeps_per_measurement);
  concurrency.pack(buffer, buffer_size, position, measurements);

  concurrency.pack(buffer, buffer_size, position, do_adaptive_double_counting);

  concurrency.pack(buffer, buffer_size, position, nr_walkers);
  concurrency.pack(buffer, buffer_size, position, nr_accumulators);
  concurrency.pack(buffer, buffer_size, position, additional_steps);

  concurrency.pack(buffer, buffer_size, position, nr_HTS_threads);
}

template<class concurrency_type>
void MCI_parameters::unpack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, Sigma_file);

  
  concurrency.unpack(buffer, buffer_size, position, warm_up_sweeps);
  concurrency.unpack(buffer, buffer_size, position, number_of_sweeps_per_measurement);
  concurrency.unpack(buffer, buffer_size, position, measurements);

  concurrency.unpack(buffer, buffer_size, position, do_adaptive_double_counting);
  
  concurrency.unpack(buffer, buffer_size, position, nr_walkers);
  concurrency.unpack(buffer, buffer_size, position, nr_accumulators);
  concurrency.unpack(buffer, buffer_size, position, additional_steps);

  concurrency.unpack(buffer, buffer_size, position, nr_HTS_threads);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/
template<class read_write_type>
void MCI_parameters::read_write(read_write_type& read_write_obj)
{
  try
    {
      read_write_obj.open_group("Monte-Carlo-Integration");

      try { read_write_obj.execute("Sigma-file"              , Sigma_file);                       } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("warm-up-sweeps"          , warm_up_sweeps);                   } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("sweeps-per-measurement"  , number_of_sweeps_per_measurement); } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("measurements"            , measurements);                     } catch(const std::exception& r_e) {}


      try { read_write_obj.execute("adaptive-double-counting", do_adaptive_double_counting);      } catch(const std::exception& r_e) {}


      {
	read_write_obj.open_group("MC-posix-parameters");

	try { read_write_obj.execute("nr-walkers"      , nr_walkers);       } catch(const std::exception& r_e) {}
	try { read_write_obj.execute("nr-accumulators" , nr_accumulators);  } catch(const std::exception& r_e) {}	
	try { read_write_obj.execute("additional-steps", additional_steps); } catch(const std::exception& r_e) {}

	try { read_write_obj.execute("HTS-threads"     , nr_HTS_threads);   } catch(const std::exception& r_e) {}

	read_write_obj.close_group();
      }

      read_write_obj.close_group();
    }
  catch(const std::exception& r_e) 
    { 
      std::cout << "\n\t MCI-parameters defined !!  \n\n";
      throw std::logic_error(__PRETTY_FUNCTION__);
    }
}

/******************************************
 ***        DATA                        ***
 ******************************************/

std::string MCI_parameters::get_Sigma_file()
{
  return Sigma_file;
}

int MCI_parameters::get_warm_up_sweeps()
{
  return warm_up_sweeps;
}

double MCI_parameters::get_number_of_sweeps_per_measurement()
{
  return number_of_sweeps_per_measurement;
}

int MCI_parameters::get_number_of_measurements()
{
  return measurements;
}

bool MCI_parameters::adjust_self_energy_for_double_counting()
{
  if(do_adaptive_double_counting == "false")
    return false;

  if(do_adaptive_double_counting == "true")
    return true;

  throw std::logic_error("do_adaptive_double_counting needs to be true or false !!! ");
}


int MCI_parameters::get_nr_walkers()
{
  return nr_walkers;
}

int MCI_parameters::get_nr_accumulators()
{
  return nr_accumulators;
}

int MCI_parameters::get_additional_steps()
{
  return additional_steps;
}

int MCI_parameters::get_nr_HTS_threads()
{
  return nr_HTS_threads;
}

#endif
