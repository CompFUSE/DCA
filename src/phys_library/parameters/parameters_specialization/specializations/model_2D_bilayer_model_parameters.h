//-*-C++-*-

#ifndef TB_BILAYER_SQUARE_HUBBARD_MODEL_PARAMETERS_H
#define TB_BILAYER_SQUARE_HUBBARD_MODEL_PARAMETERS_H

/*!
 *  \author Peter Staar
 */
template<typename dca_point_group_t, typename interaction_t>
class model_parameters<tight_binding_model<bilayer_lattice<dca_point_group_t>, interaction_t> >
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

  double  get_t();
  double  get_t_prime();

  double  get_t_perp ();

  double  get_U();
  double  get_U_prime();

  double  get_V();
  double  get_V_prime();

private:
  
  double t;
  double t_prime;

  double t_perp;
 
  double U;
  double U_prime;

  double V;
  double V_prime;
};

template<typename dca_point_group_t, typename interaction_t>
model_parameters<tight_binding_model<bilayer_lattice<dca_point_group_t>, interaction_t> >::model_parameters():
  t(1),
  t_prime(0),

  t_perp(0),

  U(4),
  U_prime(0),

  V(0),
  V_prime(0)
{}

template<typename dca_point_group_t, typename interaction_t>
model_parameters<tight_binding_model<bilayer_lattice<dca_point_group_t>, interaction_t> >::~model_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<typename dca_point_group_t, typename interaction_t>
template<class concurrency_type>
int model_parameters<tight_binding_model<bilayer_lattice<dca_point_group_t>, interaction_t> >::get_buffer_size( concurrency_type& concurrency)
{
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(t);
  buffer_size += concurrency.get_buffer_size(t_prime);
  buffer_size += concurrency.get_buffer_size(t_perp);

  buffer_size += concurrency.get_buffer_size(U);
  buffer_size += concurrency.get_buffer_size(U_prime);

  buffer_size += concurrency.get_buffer_size(V);
  buffer_size += concurrency.get_buffer_size(V_prime);

  return buffer_size;
}

template<typename dca_point_group_t, typename interaction_t>
template<class concurrency_type>
void model_parameters<tight_binding_model<bilayer_lattice<dca_point_group_t>, interaction_t> >::pack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.pack(buffer, buffer_size, position, t);
  concurrency.pack(buffer, buffer_size, position, t_prime);
  concurrency.pack(buffer, buffer_size, position, t_perp);

  concurrency.pack(buffer, buffer_size, position, U);
  concurrency.pack(buffer, buffer_size, position, U_prime);

  concurrency.pack(buffer, buffer_size, position, V);
  concurrency.pack(buffer, buffer_size, position, V_prime);
}

template<typename dca_point_group_t, typename interaction_t>
template<class concurrency_type>
void model_parameters<tight_binding_model<bilayer_lattice<dca_point_group_t>, interaction_t> >::unpack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, t);
  concurrency.unpack(buffer, buffer_size, position, t_prime);
  concurrency.unpack(buffer, buffer_size, position, t_perp);

  concurrency.unpack(buffer, buffer_size, position, U);
  concurrency.unpack(buffer, buffer_size, position, U_prime);

  concurrency.unpack(buffer, buffer_size, position, V);
  concurrency.unpack(buffer, buffer_size, position, V_prime);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<typename dca_point_group_t, typename interaction_t>
template<class read_write_type>
void  model_parameters<tight_binding_model<bilayer_lattice<dca_point_group_t>, interaction_t> >::read_write(read_write_type& read_write_obj)
{
  try
    {
      read_write_obj.open_group("bilayer-model");

      try { read_write_obj.execute("t" , t);       } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("t'", t_prime); } catch(const std::exception& r_e) {}	
      try { read_write_obj.execute("tz", t_perp); } catch(const std::exception& r_e) {}	

      try { read_write_obj.execute("U" , U);       } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("U'", U_prime); } catch(const std::exception& r_e) {}

      try { read_write_obj.execute("V" , V);       } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("V'", V_prime); } catch(const std::exception& r_e) {}

      read_write_obj.close_group();
    }
  catch(const std::exception& r_e) 
    { 
      std::cout << "\n\t bilayer-model parameters defined !!\n\n";
      throw std::logic_error(__PRETTY_FUNCTION__);
    }

  {
    std::stringstream ss;

    ss << "\n\n";
    ss << "\tbilayer-model : \n";
    ss << "\t--------------- \n\n";
    ss << "\t\t t  : " << t << "\n";
    ss << "\t\t t' : " << t_prime << "\n";
    ss << "\t\t tz : " << t_perp << "\n";

    ss << "\t\t U  : " << U << "\n";
    ss << "\t\t U' : " << U_prime << "\n";

    ss << "\t\t V  : " << V << "\n";
    ss << "\t\t V' : " << V_prime << "\n";

    ss << "\n\n";

    std::cout << ss.str();
  }
}

/******************************************
 ***        DATA                        ***
 ******************************************/

template<typename dca_point_group_t, typename interaction_t>
double model_parameters<tight_binding_model<bilayer_lattice<dca_point_group_t>, interaction_t> >::get_t()
{
  return t;
}

template<typename dca_point_group_t, typename interaction_t>
double model_parameters<tight_binding_model<bilayer_lattice<dca_point_group_t>, interaction_t> >::get_t_prime()
{
  return t_prime;
}

template<typename dca_point_group_t, typename interaction_t>
double model_parameters<tight_binding_model<bilayer_lattice<dca_point_group_t>, interaction_t> >::get_t_perp()
{
  return t_perp;
}

template<typename dca_point_group_t, typename interaction_t>
double model_parameters<tight_binding_model<bilayer_lattice<dca_point_group_t>, interaction_t> >::get_U()
{
  return U;
}
 
template<typename dca_point_group_t, typename interaction_t>
double model_parameters<tight_binding_model<bilayer_lattice<dca_point_group_t>, interaction_t> >::get_U_prime()
{
  return U_prime;
}

template<typename dca_point_group_t, typename interaction_t>
double model_parameters<tight_binding_model<bilayer_lattice<dca_point_group_t>, interaction_t> >::get_V()
{
  return V;
}
 
template<typename dca_point_group_t, typename interaction_t>
double model_parameters<tight_binding_model<bilayer_lattice<dca_point_group_t>, interaction_t> >::get_V_prime()
{
  return V_prime;
}

#endif
