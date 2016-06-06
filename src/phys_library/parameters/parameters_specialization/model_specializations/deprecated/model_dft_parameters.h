//-*-C++-*-

#ifndef DFT_MODEL_FILE_BASED_PARAMETERS_H
#define DFT_MODEL_FILE_BASED_PARAMETERS_H

/*!
 *  \author: peter staar
 */
template<int DIMENSION, typename point_group_type>
class model_parameters<dft_model<DIMENSION, point_group_type> >
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

template<int DIMENSION, typename point_group_type>
model_parameters<dft_model<DIMENSION, point_group_type> >::model_parameters()
{}

template<int DIMENSION, typename point_group_type>
model_parameters<dft_model<DIMENSION, point_group_type> >::~model_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<int DIMENSION, typename point_group_type>
template<class concurrency_type>
int model_parameters<dft_model<DIMENSION, point_group_type> >::get_buffer_size(const concurrency_type& concurrency) const
{
  int buffer_size = 0;

   buffer_size += concurrency.getBufferSize(LDA_bands);
   buffer_size += concurrency.getBufferSize(H_LDA_file);
   buffer_size += concurrency.getBufferSize(H_k_grid_size);

  return buffer_size;
}

template<int DIMENSION, typename point_group_type>
template<class concurrency_type>
void model_parameters<dft_model<DIMENSION, point_group_type> >::pack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
   concurrency.pack(buffer, buffer_size, position, LDA_bands);
   concurrency.pack(buffer, buffer_size, position, H_LDA_file);
   concurrency.pack(buffer, buffer_size, position, H_k_grid_size);
}

template<int DIMENSION, typename point_group_type>
template<class concurrency_type>
void model_parameters<dft_model<DIMENSION, point_group_type> >::unpack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
   concurrency.unpack(buffer, buffer_size, position, LDA_bands);
   concurrency.unpack(buffer, buffer_size, position, H_LDA_file);
   concurrency.unpack(buffer, buffer_size, position, H_k_grid_size);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<int DIMENSION, typename point_group_type>
template<class stream_type>
void model_parameters<dft_model<DIMENSION, point_group_type> >::to_JSON(stream_type& ss, bool is_end)
{
  ss << "\"model-parameters\" :";
  ss << "\n{ \n";

  JSON_writer::write(ss, "H_k_bands"    , LDA_bands);
  JSON_writer::write(ss, "H_k_file"     , H_LDA_file);
  JSON_writer::write(ss, "H_k_grid_size", H_k_grid_size, true);
  
  if(is_end)
    ss << "}\n";
  else
    ss << "},\n";
}
  
template<int DIMENSION, typename point_group_type>
template<class JSON_reader_type>
void model_parameters<dft_model<DIMENSION, point_group_type> >::from_JSON(JSON_reader_type& reader)
{
  typedef typename JSON_reader_type::JsonAccessor JsonAccessor;
   const JsonAccessor control(reader["model-parameters"]);
    
   LDA_bands     <= control["H_k_bands"];
   H_LDA_file    <= control["H_k_file"];
   H_k_grid_size <= control["H_k_grid_size"];

   {
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

template<int DIMENSION, typename point_group_type>
int model_parameters<dft_model<DIMENSION, point_group_type> >::get_LDA_bands()
{
  return LDA_bands;
}

template<int DIMENSION, typename point_group_type>
string model_parameters<dft_model<DIMENSION, point_group_type> >::get_H_LDA_file()
{
  return H_LDA_file;
}

template<int DIMENSION, typename point_group_type>
std::vector<int>& model_parameters<dft_model<DIMENSION, point_group_type> >::get_H_k_grid_size()
{
  return H_k_grid_size;
}

#endif
