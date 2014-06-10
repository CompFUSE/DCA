//-*-C++-*-

#ifndef PHYSICS_PARAMETERS_H
#define PHYSICS_PARAMETERS_H

/*!
 *   \ingroup  PARAMETERS
 *
 *   \author   Peter Staar
 *   \brief    ...
 */
class physics_parameters 
{
public:

  physics_parameters();
  ~physics_parameters();

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

  template<class concurrency_type>
  int  get_buffer_size( concurrency_type& concurrency);

  template<class concurrency_type>
  void pack           ( concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  template<class concurrency_type>
  void unpack         ( concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

  template<class stream_type>
  void to_JSON(stream_type& ss, bool is_end=false);
  
  template<class JSON_reader_type>
  void from_JSON(JSON_reader_type& reader);

  template<class read_write_type>
  void read_write(read_write_type& read_write_obj);

/******************************************
 ***        DATA                        ***
 ******************************************/

  double  get_beta();

  bool    adjust_chemical_potential();
  double  get_density();
  double& get_chemical_potential();

private:

  double beta;

  int                 beta_index;
  std::vector<double> beta_vector;
  
  std::string adjusting_chemical_potential;

  double density;
  double chemical_potential;
};

physics_parameters::physics_parameters():
  beta(1.),

  beta_index(-1),
  beta_vector(0),

  adjusting_chemical_potential("false"),

  density(1.),
  chemical_potential(0)
{}

physics_parameters::~physics_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<class concurrency_type>
int physics_parameters::get_buffer_size( concurrency_type& concurrency)
{
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(beta);

  buffer_size += concurrency.get_buffer_size(beta_index);
  buffer_size += concurrency.get_buffer_size(beta_vector);

  buffer_size += concurrency.get_buffer_size(adjusting_chemical_potential);

  buffer_size += concurrency.get_buffer_size(density);
  buffer_size += concurrency.get_buffer_size(chemical_potential);

  return buffer_size;
}

template<class concurrency_type>
void physics_parameters::pack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.pack(buffer, buffer_size, position, beta);

  concurrency.pack(buffer, buffer_size, position, beta_index);
  concurrency.pack(buffer, buffer_size, position, beta_vector);

  concurrency.pack(buffer, buffer_size, position, adjusting_chemical_potential);

  concurrency.pack(buffer, buffer_size, position, density);
  concurrency.pack(buffer, buffer_size, position, chemical_potential);
}

template<class concurrency_type>
void physics_parameters::unpack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, beta);

  concurrency.unpack(buffer, buffer_size, position, beta_index);
  concurrency.unpack(buffer, buffer_size, position, beta_vector);

  concurrency.unpack(buffer, buffer_size, position, adjusting_chemical_potential);

  concurrency.unpack(buffer, buffer_size, position, density);
  concurrency.unpack(buffer, buffer_size, position, chemical_potential);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<class stream_type>
void physics_parameters::to_JSON(stream_type& ss, bool is_end)
{
//   ss << "\"physics-parameters\" :";
//   ss << "\n{ \n";

//   JSON_writer::write(ss, "beta"                     , beta);
//   JSON_writer::write(ss, "adjust_chemical_potential", adjusting_chemical_potential);
//   JSON_writer::write(ss, "density"                  , density);
//   JSON_writer::write(ss, "chemical_potential"       , chemical_potential, true);

//   if(is_end)
//     ss << "}\n";
//   else
//     ss << "},\n";
}
  
template<class JSON_reader_type>
void physics_parameters::from_JSON(JSON_reader_type& reader)
{
//   typedef typename JSON_reader_type::JsonAccessor JsonAccessor;
//   const JsonAccessor control(reader["physics-parameters"]);
    
//   beta               <= control["beta"];

//   adjusting_chemical_potential <= control["adjust_chemical_potential"];
//   density                      <= control["density"];
//   chemical_potential           <= control["chemical_potential"];
}

template<class read_write_type>
void physics_parameters::read_write(read_write_type& read_write_obj)
{
  try
    {
      read_write_obj.open_group("physics-parameters");
      
      try { read_write_obj.execute("beta"                     , beta);                         } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("beta-index"               , beta_index);                   } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("beta-vector"              , beta_vector);                  } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("density"                  , density);                      } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("chemical-potential"       , chemical_potential);           } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("adjust-chemical-potential", adjusting_chemical_potential); } catch(const std::exception& r_e) {}

      read_write_obj.close_group();
    }
  catch(const std::exception& r_e) 
    { 
      cout << "\n\t physics-parameters defined !!  \n\n";
      throw std::logic_error(__PRETTY_FUNCTION__);
    }

  if(beta_vector.size()==0){
    beta_index = 0;
    beta_vector.push_back(beta);
  }
}

/******************************************
 ***        DATA                        ***
 ******************************************/

double physics_parameters::get_beta()
{
  return beta;
}

bool physics_parameters::adjust_chemical_potential()
{
  assert(adjusting_chemical_potential == "true" || 
	 adjusting_chemical_potential == "false");

  bool OK;
  adjusting_chemical_potential == "true" ? OK = true : OK = false;
  return OK;
}

double physics_parameters::get_density()
{
  return density;
}

double& physics_parameters::get_chemical_potential()
{
  return chemical_potential;
}

#endif
