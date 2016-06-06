//-*-C++-*-

#ifndef MODEL_FILE_BASED_PARAMETERS_H
#define MODEL_FILE_BASED_PARAMETERS_H

/*!
 *  \author: peter staar
 */
template<>
class model_parameters<Koshevnikov_model>
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

   int     get_LDA_bands();
   string  get_H_LDA_file();

  std::vector<int>& get_H_k_grid_size();

 private:
  
  int    LDA_bands;
  string H_LDA_file;

  std::vector<int> H_k_grid_size;
};

model_parameters<Koshevnikov_model >::model_parameters()
{}

model_parameters<Koshevnikov_model >::~model_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<class concurrency_type>
int model_parameters<Koshevnikov_model >::get_buffer_size(const concurrency_type& concurrency) const
{
  int buffer_size = 0;

   buffer_size += concurrency.getBufferSize(LDA_bands);
   buffer_size += concurrency.getBufferSize(H_LDA_file);
   buffer_size += concurrency.getBufferSize(H_k_grid_size);

  return buffer_size;
}

template<class concurrency_type>
void model_parameters<Koshevnikov_model >::pack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
   concurrency.pack(buffer, buffer_size, position, LDA_bands);
   concurrency.pack(buffer, buffer_size, position, H_LDA_file);
   concurrency.pack(buffer, buffer_size, position, H_k_grid_size);
}

template<class concurrency_type>
void model_parameters<Koshevnikov_model >::unpack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
   concurrency.unpack(buffer, buffer_size, position, LDA_bands);
   concurrency.unpack(buffer, buffer_size, position, H_LDA_file);
   concurrency.unpack(buffer, buffer_size, position, H_k_grid_size);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<class stream_type>
void model_parameters<Koshevnikov_model >::to_JSON(stream_type& ss, bool is_end)
{
  ss << "\"model-parameters\" :";
  ss << "\n{ \n";

   JSON_writer::write(ss, "H_k_bands"    , LDA_bands);
   JSON_writer::write(ss, "H_k_file"   , H_LDA_file);
   JSON_writer::write(ss, "H_k_grid_size", H_k_grid_size, true);

  if(is_end)
    ss << "}\n";
  else
    ss << "},\n";
}
  
template<class JSON_reader_type>
void model_parameters<Koshevnikov_model >::from_JSON(JSON_reader_type& reader)
{
  typedef typename JSON_reader_type::JsonAccessor JsonAccessor;
  const JsonAccessor control(reader["model-parameters"]);
  
  LDA_bands     <= control["H_k_bands"];
  H_LDA_file    <= control["H_k_file"];
  H_k_grid_size <= control["H_k_grid_size"];
  
  {// check that the file exists
    FILE* file_ptr = fopen(&(H_LDA_file[0]), "r");
    if(file_ptr == NULL){
      cout << "\t no file with name : " << H_LDA_file << endl;
      throw std::logic_error(__FUNCTION__);
    }
  }
}


/******************************************
 ***        DATA                        ***
 ******************************************/

int model_parameters<Koshevnikov_model >::get_LDA_bands()
{
  return LDA_bands;
}

string model_parameters<Koshevnikov_model >::get_H_LDA_file()
{
  return H_LDA_file;
}

std::vector<int>& model_parameters<Koshevnikov_model >::get_H_k_grid_size()
{
  return H_k_grid_size;
}

#endif
