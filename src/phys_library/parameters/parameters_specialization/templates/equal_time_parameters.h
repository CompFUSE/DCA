//-*-C++-*-

#ifndef EQUAL_TIME_PARAMETERS_H
#define EQUAL_TIME_PARAMETERS_H

/*!
 *  \class  equalt_time_parameters
 *  \ingroup  PARAMETERS
 *
 *  \author Peter Staar
 *  \brief  This class contains all possible parameters that define the two-particle Greens-function and the analysis of it.
 */
class equal_time_parameters
{
#include "type_definitions.h"

public:

  equal_time_parameters();
  ~equal_time_parameters();

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

  template<class read_write_type>
  void read_write(read_write_type& read_write_obj);

  /******************************************
   ***        DATA                        ***
   ******************************************/

  bool do_equal_time_measurements();

private:

  std::string do_equal_time_measurements_str;

};

equal_time_parameters::equal_time_parameters():
  do_equal_time_measurements_str("false")
{}

equal_time_parameters::~equal_time_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<class concurrency_type>
int equal_time_parameters::get_buffer_size(concurrency_type& concurrency) 
{
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(do_equal_time_measurements_str);

  return buffer_size;
}

template<class concurrency_type>
void equal_time_parameters::pack(concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.pack(buffer, buffer_size, position, do_equal_time_measurements_str);
}

template<class concurrency_type>
void equal_time_parameters::unpack(concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, do_equal_time_measurements_str);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<class read_write_type>
void equal_time_parameters::read_write(read_write_type& read_write_obj)
{
  try
    {
      read_write_obj.open_group("equal-time-observables");
      
      try { read_write_obj.execute("do-equal-time-measurements", do_equal_time_measurements_str); } catch(const std::exception& r_e) {}

      read_write_obj.close_group();
    }
  catch(const std::exception& r_e) 
    {}
}

/******************************************
 ***        DATA                        ***
 ******************************************/

bool equal_time_parameters::do_equal_time_measurements()
{
  if(do_equal_time_measurements_str == "true")
    return true;
  
  if(do_equal_time_measurements_str == "false")
    return false;
  
  throw std::logic_error(__FUNCTION__);
}

#endif
