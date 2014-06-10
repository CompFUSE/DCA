//-*-C++-*-

#ifndef FUNCTION_PARAMETERS_H
#define FUNCTION_PARAMETERS_H

/*!
 *   \ingroup  PARAMETERS
 *
 *   \author   Peter Staar
 *   \brief    ...
 */
class function_parameters 
{
#include "type_definitions.h"

public:

  function_parameters();
  ~function_parameters();

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

//   int    get_number_of_positive_times();
//   int    get_number_of_positive_frequencies();

//   int  get_wn_c();
//   int& get_Delta_wn();
//   int& get_N_Delta();

//   int get_nb_of_legendre_coefficients_single_particle();
//   int get_nb_of_legendre_coefficients_two_particle();
//   int get_max_nb_of_legendre_coefficients();

//   int get_Wn();

  double get_min_real_frequency();
  double get_max_real_frequency();
  int    get_number_of_real_frequencies();
  double get_real_frequencies_off_set();
  
  std::vector<int> get_H_k_grid_size();

  int get_sp_time_intervals();
  
  int get_sp_fermionic_frequencies();
  int get_sp_bosonic_frequencies();
  
  std::vector<std::vector<int> > get_sp_cluster();

  int get_tp_time_intervals();
  
  int get_tp_fermionic_frequencies();
  int get_tp_bosonic_frequencies();
  
  std::vector<std::vector<int> > get_tp_cluster();

private:

//   int time_slices;
//   int frequency_slices;

//   int single_particle_legendre_coeffs;
//   int two_particle_legendre_coeffs;

//   int wn_c;
//   int Delta_wn;
//   int N_Delta;

//   int Wn;

//   double w_max;
//   int    N_w;
//   double idelta;

  std::vector<int> H_k_grid_size;

  int sp_time_intervals;
  
  int sp_fermionic_frequencies;
  int sp_bosonic_frequencies;
  
  std::vector<std::vector<int> > sp_cluster;

  double lower_bound;
  double upper_bound;
  
  int    nr_intervals;
  double real_axis_off_set;

  int tp_time_intervals;
  
  int tp_fermionic_frequencies;
  int tp_bosonic_frequencies;
  
  std::vector<std::vector<int> > tp_cluster;
}; 

function_parameters::function_parameters():
//   time_slices(128),
//   frequency_slices(512),

//   single_particle_legendre_coeffs(0),
//   two_particle_legendre_coeffs(0),

//   wn_c(1),
//   Delta_wn(1),
//   N_Delta(0),

//   Wn(16+1),

//   w_max(8.),
//   N_w  (20),
//   idelta(1.e-2),

  H_k_grid_size(lattice_type::DIMENSION, 8),

  sp_time_intervals(128),
  
  sp_fermionic_frequencies(256),
  sp_bosonic_frequencies(32),
  
  sp_cluster(0),

  lower_bound(-10.),
  upper_bound(10),

  nr_intervals(128),
  real_axis_off_set(0.01),
  
  tp_time_intervals(0),
  
  tp_fermionic_frequencies(0),
  tp_bosonic_frequencies(0),
  
  tp_cluster(0)
{}


function_parameters::~function_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<class concurrency_type>
int function_parameters::get_buffer_size( concurrency_type& concurrency) 
{
  int buffer_size = 0;

//   buffer_size += concurrency.get_buffer_size(time_slices);
//   buffer_size += concurrency.get_buffer_size(frequency_slices);

//   buffer_size += concurrency.get_buffer_size(single_particle_legendre_coeffs);
//   buffer_size += concurrency.get_buffer_size(two_particle_legendre_coeffs   );

//   buffer_size += concurrency.get_buffer_size(wn_c);
//   buffer_size += concurrency.get_buffer_size(Delta_wn);
//   buffer_size += concurrency.get_buffer_size(N_Delta);

//   buffer_size += concurrency.get_buffer_size(Wn);

//   buffer_size += concurrency.get_buffer_size(w_max);
//   buffer_size += concurrency.get_buffer_size(N_w);
//   buffer_size += concurrency.get_buffer_size(idelta);

  {
    buffer_size += concurrency.get_buffer_size(H_k_grid_size);
    
    buffer_size += concurrency.get_buffer_size(sp_time_intervals);
    
    buffer_size += concurrency.get_buffer_size(sp_fermionic_frequencies);
    buffer_size += concurrency.get_buffer_size(sp_bosonic_frequencies);
    
    buffer_size += concurrency.get_buffer_size(sp_cluster);
  }

  {
    buffer_size += concurrency.get_buffer_size(lower_bound);
    buffer_size += concurrency.get_buffer_size(upper_bound);

    buffer_size += concurrency.get_buffer_size(nr_intervals);
    buffer_size += concurrency.get_buffer_size(real_axis_off_set);
  }

  {
    buffer_size += concurrency.get_buffer_size(tp_time_intervals);
    
    buffer_size += concurrency.get_buffer_size(tp_fermionic_frequencies);
    buffer_size += concurrency.get_buffer_size(tp_bosonic_frequencies);
    
    buffer_size += concurrency.get_buffer_size(tp_cluster);
  }

  return buffer_size;
}

template<class concurrency_type>
void function_parameters::pack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
//   concurrency.pack(buffer, buffer_size, position, time_slices);
//   concurrency.pack(buffer, buffer_size, position, frequency_slices);

//   concurrency.pack(buffer, buffer_size, position, single_particle_legendre_coeffs);
//   concurrency.pack(buffer, buffer_size, position, two_particle_legendre_coeffs);

//   concurrency.pack(buffer, buffer_size, position, wn_c);
//   concurrency.pack(buffer, buffer_size, position, Delta_wn);
//   concurrency.pack(buffer, buffer_size, position, N_Delta);

//   concurrency.pack(buffer, buffer_size, position, Wn);

//   concurrency.pack(buffer, buffer_size, position, w_max);
//   concurrency.pack(buffer, buffer_size, position, N_w);
//   concurrency.pack(buffer, buffer_size, position, idelta);

  {
    concurrency.pack(buffer, buffer_size, position, H_k_grid_size);
    
    concurrency.pack(buffer, buffer_size, position, sp_time_intervals);
    
    concurrency.pack(buffer, buffer_size, position, sp_fermionic_frequencies);
    concurrency.pack(buffer, buffer_size, position, sp_bosonic_frequencies);
    
    concurrency.pack(buffer, buffer_size, position, sp_cluster);
  }

  {
    concurrency.pack(buffer, buffer_size, position, lower_bound);
    concurrency.pack(buffer, buffer_size, position, upper_bound);
    concurrency.pack(buffer, buffer_size, position, nr_intervals);
    concurrency.pack(buffer, buffer_size, position, real_axis_off_set);
  }

  {
    concurrency.pack(buffer, buffer_size, position, tp_time_intervals);
  
    concurrency.pack(buffer, buffer_size, position, tp_fermionic_frequencies);
    concurrency.pack(buffer, buffer_size, position, tp_bosonic_frequencies);
  
    concurrency.pack(buffer, buffer_size, position, tp_cluster);
  }
}

template<class concurrency_type>
void function_parameters::unpack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
//   concurrency.unpack(buffer, buffer_size, position, time_slices);
//   concurrency.unpack(buffer, buffer_size, position, frequency_slices);

//   concurrency.unpack(buffer, buffer_size, position, single_particle_legendre_coeffs);
//   concurrency.unpack(buffer, buffer_size, position, two_particle_legendre_coeffs);

//   concurrency.unpack(buffer, buffer_size, position, wn_c);
//   concurrency.unpack(buffer, buffer_size, position, Delta_wn);
//   concurrency.unpack(buffer, buffer_size, position, N_Delta);

//   concurrency.unpack(buffer, buffer_size, position, Wn);

//   concurrency.unpack(buffer, buffer_size, position, w_max);
//   concurrency.unpack(buffer, buffer_size, position, N_w);
//   concurrency.unpack(buffer, buffer_size, position, idelta);

  concurrency.unpack(buffer, buffer_size, position, H_k_grid_size);
  
  concurrency.unpack(buffer, buffer_size, position, sp_time_intervals);
  
  concurrency.unpack(buffer, buffer_size, position, sp_fermionic_frequencies);
  concurrency.unpack(buffer, buffer_size, position, sp_bosonic_frequencies);
  
  concurrency.unpack(buffer, buffer_size, position, sp_cluster);
  
  {
    concurrency.unpack(buffer, buffer_size, position, lower_bound);
    concurrency.unpack(buffer, buffer_size, position, upper_bound);
    concurrency.unpack(buffer, buffer_size, position, nr_intervals);
    concurrency.unpack(buffer, buffer_size, position, real_axis_off_set);
  }

  {
    concurrency.unpack(buffer, buffer_size, position, tp_time_intervals);
  
    concurrency.unpack(buffer, buffer_size, position, tp_fermionic_frequencies);
    concurrency.unpack(buffer, buffer_size, position, tp_bosonic_frequencies);
  
    concurrency.unpack(buffer, buffer_size, position, tp_cluster);
  }
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<class stream_type>
void function_parameters::to_JSON(stream_type& ss, bool is_end)
{
//   ss << "\"function-parameters\" :";
//   ss << "\n{ \n";

//   ss << "\n\t\"single-particle-functions\" : \n\t{\n";
//   JSON_writer::write(ss, "time-intervals"       , time_slices);
//   JSON_writer::write(ss, "matsubara-frequencies", frequency_slices);
//   JSON_writer::write(ss, "bosonic-frequencies"  , Wn);

//   JSON_writer::write(ss, "max-real-frequencies"      , w_max);
//   JSON_writer::write(ss, "number-of-real-frequencies", N_w);
//   JSON_writer::write(ss, "real-frequencies-off-set"  , idelta);
  
//   JSON_writer::write(ss, "nb-of-legendre-coefficients", single_particle_legendre_coeffs, true);
//   ss << "\t},\n";

//   ss << "\n\t\"two-particle-functions\" : \n\t{\n";
//   JSON_writer::write(ss, "wn_c"    , wn_c, true);
// //   JSON_writer::write(ss, "Delta_wn", Delta_wn);
// //   JSON_writer::write(ss, "N_Delta" , N_Delta);

// //   JSON_writer::write(ss, "nb-of-legendre-coefficients", two_particle_legendre_coeffs, true);

//   //JSON_writer::write(ss, "Wn" , Wn, true);
//   ss << "\t}\n";

//   if(is_end)
//     ss << "}\n";
//   else
//     ss << "},\n";
}
  
template<class JSON_reader_type>
void function_parameters::from_JSON(JSON_reader_type& reader)
{
//   typedef typename JSON_reader_type::JsonAccessor JsonAccessor;

//   {
//     const JsonAccessor control(reader["function-parameters"]["single-particle-functions"]);

//     time_slices      <= control["time-intervals"];
//     frequency_slices <= control["matsubara-frequencies"];

//     try
//       { Wn <= control["bosonic-frequencies"]; }
//     catch(const std::exception& r_e) {}

//     try
//       { 
// 	w_max  <= control["max-real-frequencies"]; 
// 	N_w    <= control["number-of-real-frequencies"]; 
// 	idelta <= control["real-frequencies-off-set"]; 
//       }
//     catch(const std::exception& r_e) {}


//     try
//       { single_particle_legendre_coeffs <= control["nb-of-legendre-coefficients"]; }
//     catch(const std::exception& r_e) {}
//   }

//   try
//   {
//     const JsonAccessor control(reader["function-parameters"]["two-particle-functions"]);
//     wn_c     <= control["wn_c"];
// //     Delta_wn <= control["Delta_wn"];
// //     N_Delta  <= control["N_Delta"];
//   }
//   catch(const std::exception& r_e) {}
}

template<class read_write_type>
void function_parameters::read_write(read_write_type& read_write_obj)
{
  try
    {
      read_write_obj.open_group("function-parameters");

      {
	read_write_obj.open_group("single-particle-functions");
	
	try { read_write_obj.execute("H(k) grid-size"       , H_k_grid_size);           } catch(const std::exception& r_e) { cout << "\n not read : H(k) grid-size \n";}

	try { read_write_obj.execute("time-intervals"       , sp_time_intervals);        } catch(const std::exception& r_e) {}

	try { read_write_obj.execute("fermionic-frequencies", sp_fermionic_frequencies); } catch(const std::exception& r_e) {}
	try { read_write_obj.execute("bosonic-frequencies"  , sp_bosonic_frequencies);   } catch(const std::exception& r_e) {}

	try { read_write_obj.execute("sp-cluster"           , sp_cluster);            } catch(const std::exception& r_e) {}

	read_write_obj.close_group();
      }

      {
	read_write_obj.open_group("two-particle-functions");
	
	try { read_write_obj.execute("time-intervals"       , tp_time_intervals);        } catch(const std::exception& r_e) {}

	try { read_write_obj.execute("fermionic-frequencies", tp_fermionic_frequencies); } catch(const std::exception& r_e) {}
	try { read_write_obj.execute("bosonic-frequencies"  , tp_bosonic_frequencies);   } catch(const std::exception& r_e) {}

	try { read_write_obj.execute("tp-cluster"           , tp_cluster);               } catch(const std::exception& r_e) {}

	read_write_obj.close_group();
      }

      {
	read_write_obj.open_group("real-axis-functions");
	
	try { read_write_obj.execute("lower-bound", lower_bound);             } catch(const std::exception& r_e) {}
	try { read_write_obj.execute("upper-bound", upper_bound);             } catch(const std::exception& r_e) {}

	try { read_write_obj.execute("nr-intervals", nr_intervals);           } catch(const std::exception& r_e) {}

	try { read_write_obj.execute("real-axis-off-set", real_axis_off_set); } catch(const std::exception& r_e) {}

	read_write_obj.close_group();
      }

      read_write_obj.close_group();
    }
  catch(const std::exception& r_e) 
    { 
      cout << "\n\t MCI-parameters defined !!  \n\n";
      throw std::logic_error(__PRETTY_FUNCTION__);
    }


}

/******************************************
 ***        DATA                        ***
 ******************************************/

// int function_parameters::get_number_of_positive_times()
// {
//   return time_slices;
// }

// int function_parameters::get_number_of_positive_frequencies()
// {
//   return frequency_slices;
// }

// int function_parameters::get_Wn()
// {
//   return Wn;
// }

double function_parameters::get_min_real_frequency()
{
  return lower_bound;
}

double function_parameters::get_max_real_frequency()
{
  return upper_bound;
}

int function_parameters::get_number_of_real_frequencies()
{
  return nr_intervals;
}

double function_parameters::get_real_frequencies_off_set()
{
  assert(real_axis_off_set>0);
  return real_axis_off_set;
}

// int function_parameters::get_nb_of_legendre_coefficients_single_particle()
// {
//   return single_particle_legendre_coeffs;
// }

// int function_parameters::get_wn_c()
// {
//   return wn_c;
// }

// int& function_parameters::get_Delta_wn()
// {
//   return Delta_wn;
// }

// int& function_parameters::get_N_Delta()
// {
//   return N_Delta;
// }

// // int function_parameters::get_nb_of_legendre_coefficients_two_particle()
// // {
// //   return two_particle_legendre_coeffs;
// // }

// int function_parameters::get_max_nb_of_legendre_coefficients()
// {
//   return max(single_particle_legendre_coeffs, two_particle_legendre_coeffs);
// }

std::vector<int> function_parameters::get_H_k_grid_size()
{
  return H_k_grid_size;
}

int function_parameters::get_sp_time_intervals()
{
  return sp_time_intervals;
}
  
int function_parameters::get_sp_fermionic_frequencies()
{
  return sp_fermionic_frequencies;
}

int function_parameters::get_sp_bosonic_frequencies()
{
  return sp_bosonic_frequencies;
}

std::vector<std::vector<int> > function_parameters::get_sp_cluster()
{
  return sp_cluster;
}

int function_parameters::get_tp_time_intervals()
{
  return tp_time_intervals;
}

int function_parameters::get_tp_fermionic_frequencies()
{
  return tp_fermionic_frequencies;
}

int function_parameters::get_tp_bosonic_frequencies()
{
  return tp_bosonic_frequencies;
}

std::vector<std::vector<int> > function_parameters::get_tp_cluster()
{
  return tp_cluster;
}

#endif
