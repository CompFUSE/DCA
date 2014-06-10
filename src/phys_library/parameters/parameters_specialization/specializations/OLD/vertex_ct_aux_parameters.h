//-*-C++-*-

/*
 *      Author: peter staar
 */

#ifndef VERTEX_CT_AUX_PARAMETERS_H
#define VERTEX_CT_AUX_PARAMETERS_H

template<>
class vertex_parameters<CT_AUX> 
{
#include "type_definitions.h"

public:

  vertex_parameters();
  ~vertex_parameters();

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

  vertex_measurement_type get_vertex_measurement_type();

  int                 get_q_channel();
  std::vector<double> get_q_vector();
  int                 get_w_channel();

protected:

  void find_q_channel();

private:

  std::string vertex_measurement_type_str;

  int                 q_channel;
  std::vector<double> q_vector;
  int                 w_channel;
};

vertex_parameters<CT_AUX>::vertex_parameters()
{}

vertex_parameters<CT_AUX>::~vertex_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<class concurrency_type>
int vertex_parameters<CT_AUX>::get_buffer_size(const concurrency_type& concurrency) const
{
  int buffer_size = 0;

  buffer_size += concurrency.getBufferSize(vertex_measurement_type_str);
  buffer_size += concurrency.getBufferSize(q_channel);
  buffer_size += concurrency.getBufferSize(q_vector);
  buffer_size += concurrency.getBufferSize(w_channel);

  return buffer_size;
}

template<class concurrency_type>
void vertex_parameters<CT_AUX>::pack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.pack(buffer, buffer_size, position, vertex_measurement_type_str);
  concurrency.pack(buffer, buffer_size, position, q_channel);
  concurrency.pack(buffer, buffer_size, position, q_vector);
  concurrency.pack(buffer, buffer_size, position, w_channel);
}

template<class concurrency_type>
void vertex_parameters<CT_AUX>::unpack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, vertex_measurement_type_str);
  concurrency.unpack(buffer, buffer_size, position, q_channel);
  concurrency.unpack(buffer, buffer_size, position, q_vector);
  concurrency.unpack(buffer, buffer_size, position, w_channel);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<class stream_type>
void vertex_parameters<CT_AUX>::to_JSON(stream_type& ss, bool is_end)
{
  ss << "\"vertex-measurement\" :";
  ss << "\n{ \n";

  JSON_writer::write(ss, "vertex-measurement-type", vertex_measurement_type_str);
  JSON_writer::write(ss, "q-vector"               , q_vector);
  JSON_writer::write(ss, "q-channel"              , q_channel);
  JSON_writer::write(ss, "w-channel"              , w_channel,true);
  
  if(is_end)
    ss << "}\n";
  else
    ss << "},\n";
}
  
template<class JSON_reader_type>
void vertex_parameters<CT_AUX>::from_JSON(JSON_reader_type& reader)
{
  typedef typename JSON_reader_type::JsonAccessor JsonAccessor;
  const JsonAccessor control(reader["vertex-measurement"]);

  vertex_measurement_type_str   <= control["vertex-measurement-type"];
  q_vector                      <= control["q-channel"];
  w_channel                     <= control["w-channel"];

  if(int(q_vector.size()) != model::DIMENSION)
    throw std::logic_error("int(q_vector.size()) != model::DIMENSION");
}


/******************************************
 ***        DATA                        ***
 ******************************************/

vertex_measurement_type vertex_parameters<CT_AUX>::get_vertex_measurement_type()
{
  if(vertex_measurement_type_str == "PARTICLE_HOLE_MAGNETIC")
    return PARTICLE_HOLE_MAGNETIC;

  if(vertex_measurement_type_str == "PARTICLE_HOLE_CHARGE")
    return PARTICLE_HOLE_CHARGE;

  if(vertex_measurement_type_str == "PARTICLE_PARTICLE_SUPERCONDUCTING")
    return PARTICLE_PARTICLE_SUPERCONDUCTING;

  if(vertex_measurement_type_str == "NONE")
    return NONE;

  throw std::logic_error(__FUNCTION__);
  return NONE;
}

int vertex_parameters<CT_AUX>::get_q_channel()
{
  return q_channel;
}

std::vector<double> vertex_parameters<CT_AUX>::get_q_vector()
{
  return q_vector;
}

int vertex_parameters<CT_AUX>::get_w_channel()
{
  return w_channel;
}

void vertex_parameters<CT_AUX>::find_q_channel()
{
  double MIN = 1.e6;
  for(int l=0; l<DCA_k_cluster_type::get_size(); l++){
    if( L2_norm(q_vector,DCA_k_cluster_type::get_elements()[l]) < MIN ){
      MIN       = L2_norm(q_vector,DCA_k_cluster_type::get_elements()[l]);
      q_channel = l;
    }
  }

//   if( MIN > 1.e-6){
//     concurrency_obj << "\t with \n\t q_channel = \t";
//     concurrency_obj << q_channel; 
//     concurrency_obj << "\n\t Q-vector = ";
//     for(size_t l=0; l<DCA_k_cluster_type::get_elements()[q_channel].size(); l++){
//       concurrency_obj << "\t";
//       concurrency_obj<< DCA_k_cluster_type::get_elements()[q_channel][l];
//     }
//     concurrency_obj << "\n";
//   }
}

#endif
