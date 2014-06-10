//-*-C++-*-

#ifndef FREQUENCY_DOMAIN_H
#define FREQUENCY_DOMAIN_H

/*!
 *  \author: peter staar
 */
class frequency_domain
{
public:

  const static int DIMENSION = 1;

  typedef double scalar_type;
  typedef double element_type;

  typedef MATH_ALGORITHMS::harmonic_dmn_1D_type dmn_specifications_type;

public:

  static int&         get_size();
  static std::string  get_name();

  static scalar_type* get_basis();
  static scalar_type* get_inverse_basis();

  static std::vector<double>& get_elements();

  static std::vector<int>&    get_integer_wave_vectors();

  static int phase(double frequency);

  template<IO::FORMAT DATA_FORMAT>
  static void read(IO::reader<DATA_FORMAT>& reader);

  template<IO::FORMAT DATA_FORMAT>
  static void write(IO::writer<DATA_FORMAT>& writer);

  template<typename parameters_t>
  static void initialize(parameters_t& parameters);
};

int& frequency_domain::get_size()
{
  static int size;
  return size;
}

std::string frequency_domain::get_name()
{
  static std::string name = "frequency-domain";
  return name;
}

frequency_domain::scalar_type* frequency_domain::get_basis()
{
  static scalar_type basis[DIMENSION];
  return basis;
}

frequency_domain::scalar_type* frequency_domain::get_inverse_basis()
{
  static scalar_type inv_basis[DIMENSION];
  return inv_basis;
}

std::vector<double>& frequency_domain::get_elements()
{
  static std::vector<double> elements;
  return elements;
}

std::vector<int>& frequency_domain::get_integer_wave_vectors()
{
  static std::vector<int> elements;
  return elements;
}

template<IO::FORMAT DATA_FORMAT>
void frequency_domain::write(IO::writer<DATA_FORMAT>& writer)
{
  writer.open_group(get_name());

  writer.execute("elements", get_elements());

  writer.close_group();
}

template<typename parameters_t>
void frequency_domain::initialize(parameters_t& parameters)
{
  get_basis()        [0] = (2.*M_PI)/parameters.get_beta();
  get_inverse_basis()[0] = parameters.get_beta()/(2.*M_PI);

  get_size() = 2*parameters.get_sp_fermionic_frequencies();

  get_elements()            .resize(get_size());
  get_integer_wave_vectors().resize(get_size());

  for(int l=0; l<parameters.get_sp_fermionic_frequencies(); l++)
    {
      get_elements()[get_size()/2+0+l] =  M_PI/parameters.get_beta()*(1+2*l);
      get_elements()[get_size()/2-1-l] = -M_PI/parameters.get_beta()*(1+2*l);

      get_integer_wave_vectors()[get_size()/2+0+l] =  (1+2*l);
      get_integer_wave_vectors()[get_size()/2-1-l] = -(1+2*l);
    }
}

#endif
