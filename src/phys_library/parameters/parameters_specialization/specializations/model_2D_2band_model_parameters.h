//-*-C++-*-

#ifndef TB_TWOBAND_SQUARE_HUBBARD_MODEL_PARAMETERS_H
#define TB_TWOBAND_SQUARE_HUBBARD_MODEL_PARAMETERS_H

/*!
 *  \author Peter Staar, Andrei Plamada
 */
template<typename dca_point_group_t>
class model_parameters<tight_binding_model<twoband_lattice<dca_point_group_t> > >
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

  double  get_ei0();
  double  get_eb0();
  double  get_t0();
  double  get_U0();

private:
  
  double ei0;
  double eb0;
  double t0;
  double U0;

};

template<typename dca_point_group_t>
model_parameters<tight_binding_model<twoband_lattice<dca_point_group_t> > >::model_parameters():
  ei0(0),
  eb0(0),
  U0(0),
  t0(0)
{}

template<typename dca_point_group_t>
model_parameters<tight_binding_model<twoband_lattice<dca_point_group_t> > >::~model_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<typename dca_point_group_t>
template<class concurrency_type>
int model_parameters<tight_binding_model<twoband_lattice<dca_point_group_t> > >::get_buffer_size( concurrency_type& concurrency)
{
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(ei0);
  buffer_size += concurrency.get_buffer_size(eb0);
  buffer_size += concurrency.get_buffer_size(t0);
  buffer_size += concurrency.get_buffer_size(U0);

  return buffer_size;
}

template<typename dca_point_group_t>
template<class concurrency_type>
void model_parameters<tight_binding_model<twoband_lattice<dca_point_group_t> > >::pack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.pack(buffer, buffer_size, position, ei0);
  concurrency.pack(buffer, buffer_size, position, eb0);
  concurrency.pack(buffer, buffer_size, position, t0);
  concurrency.pack(buffer, buffer_size, position, U0);
}

template<typename dca_point_group_t>
template<class concurrency_type>
void model_parameters<tight_binding_model<twoband_lattice<dca_point_group_t> > >::unpack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, ei0);
  concurrency.unpack(buffer, buffer_size, position, eb0);
  concurrency.unpack(buffer, buffer_size, position, t0);
  concurrency.unpack(buffer, buffer_size, position, U0);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<typename dca_point_group_t>
template<class read_write_type>
void  model_parameters<tight_binding_model<twoband_lattice<dca_point_group_t> > >::read_write(read_write_type& read_write_obj)
{
  try
    {
      read_write_obj.open_group("twoband-model");

      try { read_write_obj.execute("ei0" , ei0);       } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("eb0" , eb0);       } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("t0" , t0);       } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("U0" , U0);       } catch(const std::exception& r_e) {}

      read_write_obj.close_group();
    }
  catch(const std::exception& r_e) 
    { 
      cout << "\n\t twoband-model parameters defined !!\n\n";
      throw std::logic_error(__PRETTY_FUNCTION__);
    }

  {
    std::stringstream ss;

    ss << "\n\n";
    ss << "\t twoband-model : \n";
    ss << "\t--------------- \n\n";
    ss << "\t\t ei0  : " << ei0 << "\n";
    ss << "\t\t eb0  : " << eb0 << "\n";
    ss << "\t\t  t0  : " << t0 << "\n";
    ss << "\t\t  U0  : " << U0 << "\n";
    ss << "\n\n";

    cout << ss.str();
  }
}

/******************************************
 ***        DATA                        ***
 ******************************************/

template<typename dca_point_group_t>
double model_parameters<tight_binding_model<twoband_lattice<dca_point_group_t> > >::get_ei0()
{
  return ei0;
}

template<typename dca_point_group_t>
double model_parameters<tight_binding_model<twoband_lattice<dca_point_group_t> > >::get_eb0()
{
  return eb0;
}

template<typename dca_point_group_t>
double model_parameters<tight_binding_model<twoband_lattice<dca_point_group_t> > >::get_t0()
{
  return t0;
}

template<typename dca_point_group_t>
double model_parameters<tight_binding_model<twoband_lattice<dca_point_group_t> > >::get_U0()
{
  return U0;
}
 
#endif
