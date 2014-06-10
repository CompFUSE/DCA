//-*-C++-*-

#ifndef ANALYTIC_MODEL_PARAMETERS_H
#define ANALYTIC_MODEL_PARAMETERS_H

/*!
 *      Author: peter staar
 */
template<typename lattice_t, typename interaction_t>
class model_parameters<tight_binding_model<lattice_t, interaction_t> >
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

  std::vector<int>& get_H_k_grid_size();

  std::vector<std::vector<double> >& get_t_ij();
  std::vector<std::vector<double> >& get_U_ij();

private:

  int U_size;
  int t_size;

  std::vector<std::vector<double> > sparse_U_ij;
  std::vector<std::vector<double> > sparse_t_ij;

  std::vector<int>                  H_k_grid_size;
};

/*
template<typename lattice_t, typename interaction_t>
model_parameters<tight_binding_model<lattice_t, interaction_t> >::model_parameters()
{}

template<typename lattice_t, typename interaction_t>
model_parameters<tight_binding_model<lattice_t, interaction_t> >::~model_parameters()
{}
*/

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

/*
template<typename lattice_t, typename interaction_t>
template<class concurrency_type>
int model_parameters<tight_binding_model<lattice_t, interaction_t> >::get_buffer_size(const concurrency_type& concurrency) const
{
  std::vector<double> tmp_0;
  for(size_t i=0; i<sparse_t_ij.size(); i++)
    for(size_t j=0; j<sparse_t_ij[i].size(); j++)
      tmp_0.push_back(sparse_t_ij[i][j]);

  std::vector<double> tmp_1;
  for(size_t i=0; i<sparse_U_ij.size(); i++)
    for(size_t j=0; j<sparse_U_ij[i].size(); j++)
      tmp_1.push_back(sparse_U_ij[i][j]);

  int buffer_size = 0;

  buffer_size += concurrency.getBufferSize(U_size);
  buffer_size += concurrency.getBufferSize(t_size);
  buffer_size += concurrency.getBufferSize(tmp_0);
  buffer_size += concurrency.getBufferSize(tmp_1);
  buffer_size += concurrency.getBufferSize(H_k_grid_size);

  return buffer_size;
}

template<typename lattice_t, typename interaction_t>
template<class concurrency_type>
void model_parameters<tight_binding_model<lattice_t, interaction_t> >::pack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  std::vector<double> tmp_0;
  for(size_t i=0; i<sparse_t_ij.size(); i++)
    for(size_t j=0; j<sparse_t_ij[i].size(); j++)
      tmp_0.push_back(sparse_t_ij[i][j]);

  std::vector<double> tmp_1;
  for(size_t i=0; i<sparse_U_ij.size(); i++)
    for(size_t j=0; j<sparse_U_ij[i].size(); j++)
      tmp_1.push_back(sparse_U_ij[i][j]);

  concurrency.pack(buffer, buffer_size, position, U_size);
  concurrency.pack(buffer, buffer_size, position, t_size);
  concurrency.pack(buffer, buffer_size, position, tmp_0);
  concurrency.pack(buffer, buffer_size, position, tmp_1);
  concurrency.pack(buffer, buffer_size, position, H_k_grid_size);

}

template<typename lattice_t, typename interaction_t>
template<class concurrency_type>
void model_parameters<tight_binding_model<lattice_t, interaction_t> >::unpack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, U_size);
  concurrency.unpack(buffer, buffer_size, position, t_size);

  std::vector<double> tmp_0(t_size*4);
  std::vector<double> tmp_1(U_size*6);

  concurrency.unpack(buffer, buffer_size, position, tmp_0);
  concurrency.unpack(buffer, buffer_size, position, tmp_1);
  concurrency.unpack(buffer, buffer_size, position, H_k_grid_size);

  {
    assert(tmp_0.size()%4 == 0);
    sparse_t_ij.resize(t_size, std::vector<double>(4,0.));

    int l=0;
    for(size_t i=0; i<sparse_t_ij.size(); i++)
      for(size_t j=0; j<sparse_t_ij[i].size(); j++)
	sparse_t_ij[i][j] = tmp_0[l++];
  }

  {
    assert(tmp_1.size()%6 == 0);
    sparse_U_ij.resize(U_size, std::vector<double>(6,0.));

    int l=0;
    for(size_t i=0; i<sparse_U_ij.size(); i++)
      for(size_t j=0; j<sparse_U_ij[i].size(); j++)
	sparse_U_ij[i][j] = tmp_1[l++];
  }
}
*/

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

/*
template<typename lattice_t, typename interaction_t>
template<class stream_type>
void model_parameters<tight_binding_model<lattice_t, interaction_t> >::to_JSON(stream_type& ss, bool is_end)
{
  ss << "\"model-parameters\" :";
  ss << "\n{ \n";

  JSON_writer::write(ss, "H_k_grid_size", H_k_grid_size);

  JSON_writer::write(ss, "t_ij"         , sparse_t_ij);
  JSON_writer::write(ss, "U_ij"         , sparse_U_ij, true);

  if(is_end)
    ss << "}\n";
  else
    ss << "},\n";
}
  
template<typename lattice_t, typename interaction_t>
template<class JSON_reader_type>
void model_parameters<tight_binding_model<lattice_t, interaction_t> >::from_JSON(JSON_reader_type& reader)
{
  typedef typename JSON_reader_type::JsonAccessor JsonAccessor;
  const JsonAccessor control(reader["model-parameters"]);

  sparse_t_ij   <= control["t_ij"]; 
  sparse_U_ij   <= control["U_ij"]; 

  t_size = sparse_t_ij.size();
  U_size = sparse_U_ij.size();

  H_k_grid_size <= control["H_k_grid_size"]; 

  if(int(H_k_grid_size.size()) != lattice_t::DIMENSION)  
    throw std::logic_error("int(H_k_grid_size.size()) != model::DIMENSION");
}
*/

/******************************************
 ***        DATA                        ***
 ******************************************/

/*
template<typename lattice_t, typename interaction_t>
std::vector<std::vector<double> >& model_parameters<tight_binding_model<lattice_t, interaction_t> >::get_t_ij()
{
  return sparse_t_ij;
}

template<typename lattice_t, typename interaction_t>
std::vector<std::vector<double> >& model_parameters<tight_binding_model<lattice_t, interaction_t> >::get_U_ij()
{
  return sparse_U_ij;
}

template<typename lattice_t, typename interaction_t>
std::vector<int>& model_parameters<tight_binding_model<lattice_t, interaction_t> >::get_H_k_grid_size()
{
  return H_k_grid_size;
}
*/

#endif
