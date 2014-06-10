//-*-C++-*-

#ifndef BRILLOUIN_ZONE_PARAMETERS_H
#define BRILLOUIN_ZONE_PARAMETERS_H

/*!
 *   \ingroup  PARAMETERS
 *
 *   \author   Peter Staar
 *   \brief    ...
 */
class brillouin_zone_parameters 
{
#include "type_definitions.h"

public:

  brillouin_zone_parameters();
  ~brillouin_zone_parameters();

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

  template<class concurrency_type>
  int  get_buffer_size(concurrency_type& concurrency) ;

  template<class concurrency_type>
  void pack           (concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  template<class concurrency_type>
  void unpack         (concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

//   template<class stream_type>
//   void to_JSON(stream_type& ss, bool is_end=false);
  
//   template<class JSON_reader_type>
//   void from_JSON(JSON_reader_type& reader);

  template<class read_write_type>
  void read_write(read_write_type& read_write_obj);

/******************************************
 ***        DATA                        ***
 ******************************************/

  std::string&                       get_coordinate_type();

  std::vector<std::string>&          get_coordinate_names();
  std::vector<std::vector<double> >& get_Brillouin_zone_vectors();
    
private:

  std::string                       coordinate_type;

  std::vector<std::string>          coordinate_names;
  std::vector<std::vector<double> > coordinate_vectors;
};

brillouin_zone_parameters::brillouin_zone_parameters():
  coordinate_type("absolute"),

  coordinate_names(0),
  coordinate_vectors(0)
{}

brillouin_zone_parameters::~brillouin_zone_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<class concurrency_type>
int brillouin_zone_parameters::get_buffer_size( concurrency_type& concurrency) 
{
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(coordinate_type);
  buffer_size += concurrency.get_buffer_size(coordinate_names);
  buffer_size += concurrency.get_buffer_size(coordinate_vectors);

  return buffer_size;
}

template<class concurrency_type>
void brillouin_zone_parameters::pack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.pack(buffer, buffer_size, position, coordinate_type);
  concurrency.pack(buffer, buffer_size, position, coordinate_names);
  concurrency.pack(buffer, buffer_size, position, coordinate_vectors);
}

template<class concurrency_type>
void brillouin_zone_parameters::unpack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, coordinate_type);
  concurrency.unpack(buffer, buffer_size, position, coordinate_names);
  concurrency.unpack(buffer, buffer_size, position, coordinate_vectors);
}


/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<class read_write_type>
void brillouin_zone_parameters::read_write(read_write_type& read_write_obj)
{
  try
    {
      read_write_obj.open_group("band-structure-cut");
      
      try { read_write_obj.execute("coordinate-type"       , coordinate_type);    } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("Brillouin-zone-names"  , coordinate_names);   } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("Brillouin-zone-vectors", coordinate_vectors); } catch(const std::exception& r_e) {}
      
      read_write_obj.close_group();
    }
  catch(const std::exception& r_e) 
    {}
}

/******************************************
 ***        DATA                        ***
 ******************************************/

std::string&                       brillouin_zone_parameters::get_coordinate_type()
{
  return coordinate_type;
}

std::vector<std::string>&          brillouin_zone_parameters::get_coordinate_names()
{
  return coordinate_names;
}
 
std::vector<std::vector<double> >& brillouin_zone_parameters::get_Brillouin_zone_vectors()
{
  return coordinate_vectors;
}
  


#endif
