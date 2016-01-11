//-*-C++-*-

#ifndef CPE_PARAMETERS_H
#define CPE_PARAMETERS_H

/*!
 *   \ingroup  PARAMETERS
 *
 *   \author   Peter Staar
 *   \brief    ...
 */
class CPE_parameters 
{
#include "type_definitions.h"

public:

  CPE_parameters();
  ~CPE_parameters();

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

  template<class concurrency_type>
  int  get_buffer_size( concurrency_type& concurrency) ;

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

  bool    do_CPE();

  int&    get_N_wn();
  double  get_CPE_smoothing_factor();

  int     get_max_CPE_iterations();
  double  get_max_CPE_error();

  bool   simulate_gaussian_noise();
  int    get_nr_of_CPE_samples();
  double get_simulated_CPE_stddev();

  bool compute_free_spectrum();
  bool compute_lattice_spectrum();
  bool compute_cluster_spectrum();

private:

  std::string do_CPE_str;

  int    N_wn;
  int    smoothing_factor;

  int    max_iterations;
  double max_error;

  std::string simulate_gaussian_noise_str;
  int    nr_of_samples;
  double simulated_stddev;

  std::string compute_free_spectrum_str;
  std::string compute_lattice_spectrum_str;
  std::string compute_cluster_spectrum_str;
};

CPE_parameters::CPE_parameters():
  do_CPE_str("false"),

  N_wn(64),

  smoothing_factor(1),

  max_iterations(100),
  max_error(0.),

  simulate_gaussian_noise_str("false"),
  nr_of_samples(1),
  simulated_stddev(0),

  compute_free_spectrum_str   ("false"),
  compute_lattice_spectrum_str("false"),
  compute_cluster_spectrum_str("false")
{}

CPE_parameters::~CPE_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<class concurrency_type>
int CPE_parameters::get_buffer_size( concurrency_type& concurrency) 
{
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(do_CPE_str);

  buffer_size += concurrency.get_buffer_size(N_wn);
  buffer_size += concurrency.get_buffer_size(smoothing_factor);

  buffer_size += concurrency.get_buffer_size(max_iterations);
  buffer_size += concurrency.get_buffer_size(max_error);

  buffer_size += concurrency.get_buffer_size(simulate_gaussian_noise_str);
  buffer_size += concurrency.get_buffer_size(nr_of_samples);
  buffer_size += concurrency.get_buffer_size(simulated_stddev);

  buffer_size += concurrency.get_buffer_size(compute_free_spectrum_str);
  buffer_size += concurrency.get_buffer_size(compute_lattice_spectrum_str);
  buffer_size += concurrency.get_buffer_size(compute_cluster_spectrum_str);

  return buffer_size;
}

template<class concurrency_type>
void CPE_parameters::pack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.pack(buffer, buffer_size, position, do_CPE_str);

  concurrency.pack(buffer, buffer_size, position, N_wn);
  concurrency.pack(buffer, buffer_size, position, smoothing_factor);

  concurrency.pack(buffer, buffer_size, position, max_iterations);
  concurrency.pack(buffer, buffer_size, position, max_error);

  concurrency.pack(buffer, buffer_size, position, simulate_gaussian_noise_str);
  concurrency.pack(buffer, buffer_size, position, nr_of_samples);
  concurrency.pack(buffer, buffer_size, position, simulated_stddev);

  concurrency.pack(buffer, buffer_size, position, compute_free_spectrum_str);
  concurrency.pack(buffer, buffer_size, position, compute_lattice_spectrum_str);
  concurrency.pack(buffer, buffer_size, position, compute_cluster_spectrum_str);
}

template<class concurrency_type>
void CPE_parameters::unpack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, do_CPE_str);

  concurrency.unpack(buffer, buffer_size, position, N_wn);
  concurrency.unpack(buffer, buffer_size, position, smoothing_factor);

  concurrency.unpack(buffer, buffer_size, position, max_iterations);
  concurrency.unpack(buffer, buffer_size, position, max_error);

  concurrency.unpack(buffer, buffer_size, position, simulate_gaussian_noise_str);
  concurrency.unpack(buffer, buffer_size, position, nr_of_samples);
  concurrency.unpack(buffer, buffer_size, position, simulated_stddev);

  concurrency.unpack(buffer, buffer_size, position, compute_free_spectrum_str);
  concurrency.unpack(buffer, buffer_size, position, compute_lattice_spectrum_str);
  concurrency.unpack(buffer, buffer_size, position, compute_cluster_spectrum_str);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<class read_write_type>
void CPE_parameters::read_write(read_write_type& read_write_obj)
{
  try
    {
      read_write_obj.open_group("CPE-parameters");

      try { read_write_obj.execute("do-CPE"                        , do_CPE_str);     } catch(const std::exception& r_e) {}

      try { read_write_obj.execute("number-of-matsubara-freqencies", N_wn);           } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("smoothing-factor"              , smoothing_factor); } catch(const std::exception& r_e) {}

      try { read_write_obj.execute("max-CPE-iterations"            , max_iterations); } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("max-CPE-error"                 , max_error);      } catch(const std::exception& r_e) {}

      try { read_write_obj.execute("simulate-Gaussian-noise"       , simulate_gaussian_noise_str ); } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("nr-of-samples"                 , nr_of_samples               ); } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("simulated-stddev"              , simulated_stddev            ); } catch(const std::exception& r_e) {}

      try { read_write_obj.execute("compute-free-spectrum"   , compute_free_spectrum_str);    } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("compute-lattice-spectrum", compute_lattice_spectrum_str); } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("compute-cluster-spectrum", compute_cluster_spectrum_str); } catch(const std::exception& r_e) {}

      read_write_obj.close_group();
    }
  catch(const std::exception& r_e) 
    { 
      //cout << "\n\t CPE-parameters not well-defined !!  \n\n";
      //throw std::logic_error(__PRETTY_FUNCTION__);
    }
}

/******************************************
 ***        DATA                        ***
 ******************************************/

bool CPE_parameters::do_CPE()
{
  if(do_CPE_str == "true")
    return true;
  
  if(do_CPE_str == "false")
    return false;
  
  throw std::logic_error(__FUNCTION__);
}

int& CPE_parameters::get_N_wn()
{
  return N_wn;
}

double CPE_parameters::get_CPE_smoothing_factor()
{
  return smoothing_factor;
}

int CPE_parameters::get_max_CPE_iterations()
{
  return max_iterations;
}

double CPE_parameters::get_max_CPE_error()
{
  return max_error;
}

bool CPE_parameters::simulate_gaussian_noise()
{
  if(simulate_gaussian_noise_str == "true")
    return true;
  else
    return false;
}

int CPE_parameters::get_nr_of_CPE_samples()
{
  return nr_of_samples;
}

double CPE_parameters::get_simulated_CPE_stddev()
{
  return simulated_stddev;
}

bool CPE_parameters::compute_free_spectrum()
{
  if(compute_free_spectrum_str == "true")
    return true;
  else
    return false;
}

bool CPE_parameters::compute_lattice_spectrum()
{
  if(compute_lattice_spectrum_str == "true")
    return true;
  else
    return false;
}

bool CPE_parameters::compute_cluster_spectrum()
{
  if(compute_cluster_spectrum_str == "true")
    return true;
  else
    return false;
}

#endif

