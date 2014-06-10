//-*-C++-*-

#ifndef CUPRATE_THREE_BAND_MODEL_PARAMETERS_H
#define CUPRATE_THREE_BAND_MODEL_PARAMETERS_H

/*!
 *  \author: peter staar
 */
template<typename dca_point_group_t>
class model_parameters<cuprate_three_band_model<dca_point_group_t> >
{

public:

  model_parameters();
  ~model_parameters();

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

  template<class concurrency_type>
  int  get_buffer_size(const concurrency_type& concurrency) const;

  template<class concurrency_type>
  void pack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  template<class concurrency_type>
  void unpack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

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

  double get_U_dd();
  double get_U_pp();

  std::string       get_compound_file_name();

  std::vector<int>& get_H_k_grid_size();

private:
  
  double U_dd;
  double U_pp;

  std::string      compound_file_name;

  std::vector<int> H_k_grid_size;
};

template<typename dca_point_group_t>
model_parameters<cuprate_three_band_model<dca_point_group_t> >::model_parameters():
  U_dd(0),
  U_pp(0),

  compound_file_name("DEFAULT"),
  H_k_grid_size(2, 16)
{}

template<typename dca_point_group_t>
model_parameters<cuprate_three_band_model<dca_point_group_t> >::~model_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<typename dca_point_group_t>
template<class concurrency_type>
int model_parameters<cuprate_three_band_model<dca_point_group_t> >::get_buffer_size(const concurrency_type& concurrency) const
{
  int buffer_size = 0;

  buffer_size += concurrency.getBufferSize(U_dd);
  buffer_size += concurrency.getBufferSize(U_pp);

  buffer_size += concurrency.getBufferSize(compound_file_name);
  buffer_size += concurrency.getBufferSize(H_k_grid_size);

  return buffer_size;
}

template<typename dca_point_group_t>
template<class concurrency_type>
void model_parameters<cuprate_three_band_model<dca_point_group_t> >::pack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.pack(buffer, buffer_size, position, U_dd);
  concurrency.pack(buffer, buffer_size, position, U_pp);

  concurrency.pack(buffer, buffer_size, position, compound_file_name);
  concurrency.pack(buffer, buffer_size, position, H_k_grid_size);
}

template<typename dca_point_group_t>
template<class concurrency_type>
void model_parameters<cuprate_three_band_model<dca_point_group_t> >::unpack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, U_dd);
  concurrency.unpack(buffer, buffer_size, position, U_pp);

  concurrency.unpack(buffer, buffer_size, position, compound_file_name);
  concurrency.unpack(buffer, buffer_size, position, H_k_grid_size);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<typename dca_point_group_t>
template<class stream_type>
void model_parameters<cuprate_three_band_model<dca_point_group_t> >::to_JSON(stream_type& ss, bool is_end)
{
  ss << "\"model-parameters\" :";
  ss << "\n{ \n";

  JSON_writer::write(ss, "U_dd", U_dd);
  JSON_writer::write(ss, "U_pp", U_pp);

  JSON_writer::write(ss, "compound-name" , compound_file_name);
  JSON_writer::write(ss, "H(k) grid-size", H_k_grid_size, true);

  if(is_end)
    ss << "}\n";
  else
    ss << "},\n";
}
  
template<typename dca_point_group_t>
template<class JSON_reader_type>
void model_parameters<cuprate_three_band_model<dca_point_group_t> >::from_JSON(JSON_reader_type& reader)
{
  typedef typename JSON_reader_type::JsonAccessor JsonAccessor;
  const JsonAccessor control(reader["model-parameters"]);

  U_dd <= control["U_dd"]; 
  U_pp <= control["U_pp"]; 

  compound_file_name <= control["compound-file-name"]; 

  {
    FILE* file_ptr = fopen(&(compound_file_name[0]), "r");
    if(file_ptr == NULL)
      {
	cout << "\t no file with name : " << compound_file_name << endl;
	throw std::logic_error(__FUNCTION__);
      }
  }

  H_k_grid_size <= control["H(k) grid-size"]; 

  if(int(H_k_grid_size.size()) != Bett_cluster_square_2D<dca_point_group_t>::DIMENSION)  
    throw std::logic_error("int(H_k_grid_size.size()) != model::DIMENSION");
}

/******************************************
 ***        DATA                        ***
 ******************************************/

template<typename dca_point_group_t>
double model_parameters<cuprate_three_band_model<dca_point_group_t> >::get_U_dd()
{
  return U_dd;
}

template<typename dca_point_group_t>
double model_parameters<cuprate_three_band_model<dca_point_group_t> >::get_U_pp()
{
  return U_pp;
}

template<typename dca_point_group_t>
std::string model_parameters<cuprate_three_band_model<dca_point_group_t> >::get_compound_file_name()
{
  return compound_file_name;
}

template<typename dca_point_group_t>
std::vector<int>& model_parameters<cuprate_three_band_model<dca_point_group_t> >::get_H_k_grid_size()
{
  return H_k_grid_size;
}

#endif
