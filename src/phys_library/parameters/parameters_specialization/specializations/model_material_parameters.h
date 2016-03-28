//-*-C++-*-

#ifndef MATERIAL_MODEL_PARAMETERS_H
#define MATERIAL_MODEL_PARAMETERS_H

/*!
 *  \author: peter staar
 */
template<material_name_type name, typename dca_point_group_t>
class model_parameters<tight_binding_model<material_lattice<name, dca_point_group_t> > >
{
  
public:
  
  model_parameters();
  ~model_parameters();

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

  template<class concurrency_type>
  int  get_buffer_size( concurrency_type& concurrency);

  template<class concurrency_type>
  void pack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  template<class concurrency_type>
  void unpack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

  template<class read_write_type>
  void read_write(read_write_type& read_write_obj);

/******************************************
 ***        DATA                        ***
 ******************************************/

  std::string get_U_ij_file_name();
  std::string get_t_ij_file_name();
  
private:
 
  std::string U_ij_file_name;
  std::string t_ij_file_name;  

};

template<material_name_type name, typename dca_point_group_t>
model_parameters<tight_binding_model<material_lattice<name, dca_point_group_t> > >::model_parameters():
  U_ij_file_name("U_ij_file_name"),
  t_ij_file_name("t_ij_file_name")//,

  //double_counting_correction(0)
{}

template<material_name_type name, typename dca_point_group_t>
model_parameters<tight_binding_model<material_lattice<name, dca_point_group_t> > >::~model_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<material_name_type name, typename dca_point_group_t>
template<class concurrency_type>
int model_parameters<tight_binding_model<material_lattice<name, dca_point_group_t> > >::get_buffer_size( concurrency_type& concurrency)
{
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(U_ij_file_name);
  buffer_size += concurrency.get_buffer_size(t_ij_file_name);

  //buffer_size += concurrency.get_buffer_size(double_counting_correction);

  return buffer_size;
}

template<material_name_type name, typename dca_point_group_t>
template<class concurrency_type>
void model_parameters<tight_binding_model<material_lattice<name, dca_point_group_t> > >::pack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.pack(buffer, buffer_size, position, U_ij_file_name);
  concurrency.pack(buffer, buffer_size, position, t_ij_file_name);

}

template<material_name_type name, typename dca_point_group_t>
template<class concurrency_type>
void model_parameters<tight_binding_model<material_lattice<name, dca_point_group_t> > >::unpack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, U_ij_file_name);
  concurrency.unpack(buffer, buffer_size, position, t_ij_file_name);

}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<material_name_type name, typename dca_point_group_t>
template<class read_write_type>
void  model_parameters<tight_binding_model<material_lattice<name, dca_point_group_t> > >::read_write(read_write_type& read_write_obj)
{
  try
    {
      read_write_obj.open_group("material-model");

      read_write_obj.execute("t_ij-filename", t_ij_file_name);
      read_write_obj.execute("U_ij-filename", U_ij_file_name);

      read_write_obj.close_group();
    }
  catch(const std::exception& r_e) 
    { 
      std::cout << "\n\t material-model parameters defined !!\n\n";
      throw std::logic_error(__PRETTY_FUNCTION__);
    }

  {
    std::stringstream ss;
    
    ss << "\n\n";
    ss << "\tmaterial-model\n";
    ss << "\t--------------\n";
    ss << "\t\tt_ij-filename   : " << t_ij_file_name << "\n";
    ss << "\t\tU_ij-filename   : " << U_ij_file_name << "\n\n";
    
    std::cout << ss.str();
  }
}

/******************************************
 ***        DATA                        ***
 ******************************************/

template<material_name_type name, typename dca_point_group_t>
std::string model_parameters<tight_binding_model<material_lattice<name, dca_point_group_t> > >::get_U_ij_file_name()
{
  return U_ij_file_name;
}

template<material_name_type name, typename dca_point_group_t>
std::string model_parameters<tight_binding_model<material_lattice<name, dca_point_group_t> > >::get_t_ij_file_name()
{
  return t_ij_file_name;
}

#endif
