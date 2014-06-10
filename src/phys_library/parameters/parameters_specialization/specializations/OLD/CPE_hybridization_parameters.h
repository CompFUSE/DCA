//-*-C++-*-

/*
 *      Author: peter staar
 */

#ifndef CPE_PARAMETERS_HYBRIDIZATION_H
#define CPE_PARAMETERS_HYBRIDIZATION_H

template<>
class CPE_parameters<HYBRIDIZATION> 
{
#include "type_definitions.h"

public:

  CPE_parameters();
  ~CPE_parameters();

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

/******************************************
 ***        DATA                        ***
 ******************************************/

  bool    do_CPE();
  double  get_CPE_Delta();
  int     get_CPE_N();

private:

  string do_CPE_str;

  double CPE_Delta;
  int    CPE_N;
};

CPE_parameters< HYBRIDIZATION>::CPE_parameters()
{}

CPE_parameters< HYBRIDIZATION>::~CPE_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<class concurrency_type>
int CPE_parameters< HYBRIDIZATION>::get_buffer_size(const concurrency_type& concurrency) const
{
  int buffer_size = 0;

  buffer_size += concurrency.getBufferSize(do_CPE_str);
  buffer_size += concurrency.getBufferSize(CPE_Delta);
  buffer_size += concurrency.getBufferSize(CPE_N);

  return buffer_size;
}

template<class concurrency_type>
void CPE_parameters< HYBRIDIZATION>::pack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.pack(buffer, buffer_size, position, do_CPE_str);
  concurrency.pack(buffer, buffer_size, position, CPE_Delta);
  concurrency.pack(buffer, buffer_size, position, CPE_N);
}

template<class concurrency_type>
void CPE_parameters< HYBRIDIZATION>::unpack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, do_CPE_str);
  concurrency.unpack(buffer, buffer_size, position, CPE_Delta);
  concurrency.unpack(buffer, buffer_size, position, CPE_N);
}


/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<class stream_type>
void CPE_parameters< HYBRIDIZATION>::to_JSON(stream_type& ss, bool is_end)
{
  ss << "\"CPE\" :";
  ss << "\n{ \n";

  JSON_writer::write(ss, "do_CPE"   , do_CPE_str);
  JSON_writer::write(ss, "CPE_Delta", CPE_Delta);
  JSON_writer::write(ss, "CPE_N"    , CPE_N    , true);

  if(is_end)
    ss << "}\n";
  else
    ss << "},\n";  
}
  
template<class JSON_reader_type>
void CPE_parameters< HYBRIDIZATION>::from_JSON(JSON_reader_type& reader)
{
  typedef typename dca::JsonReader::JsonAccessor JsonAccessor;
  const JsonAccessor control(reader["CPE"]);
  
  do_CPE_str <= control["do_CPE"];
  
  CPE_Delta  <= control["CPE_Delta"];
  CPE_N      <= control["CPE_N"];
}

/******************************************
 ***        DATA                        ***
 ******************************************/

bool CPE_parameters< HYBRIDIZATION>::do_CPE()
{
  if(do_CPE_str == "yes")
    return true;

  if(do_CPE_str == "no")
    return false;

  throw std::logic_error(__FUNCTION__);
}

double CPE_parameters< HYBRIDIZATION>::get_CPE_Delta()
{
  return CPE_Delta;
}

int CPE_parameters< HYBRIDIZATION>::get_CPE_N()
{
  return CPE_N;
}

#endif
