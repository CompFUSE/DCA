//-*-C++-*-

#ifndef FREQUENCY_DOMAIN_REAL_AXIS_H
#define FREQUENCY_DOMAIN_REAL_AXIS_H

/*!
 *  \author: peter staar
 */
class frequency_domain_real_axis
{
public:

  typedef double element_type;

public:

  static int         get_size();
  static std::string get_name();

  static std::vector<double>& get_elements();
  
  template<typename stream_type>
  static void to_JSON(stream_type& ss);

  template<IO::FORMAT DATA_FORMAT>
  static void write(IO::writer<DATA_FORMAT>& writer);

  template<typename parameters_t>
  static void initialize(parameters_t& parameters);
};

int frequency_domain_real_axis::get_size()
{
  return get_elements().size();
}

std::string frequency_domain_real_axis::get_name()
{
  static std::string name = "frequency-domain-real-axis";
  return name;
}

std::vector<double >& frequency_domain_real_axis::get_elements()
{
  static std::vector<double > vec_elements(0,0);
  return vec_elements;
} 


template<IO::FORMAT DATA_FORMAT>
void frequency_domain_real_axis::write(IO::writer<DATA_FORMAT>& writer)
{
  writer.open_group(get_name());

  writer.execute("elements", get_elements());

  writer.close_group();
}

template<typename stream_type>
void frequency_domain_real_axis::to_JSON(stream_type& ss)
{
  ss << "\"Pade_real_frequency_domain\" : [\n";
  
  for(int i=0; i<get_size(); i++)
    if(i == get_size()-1)
      ss << get_elements()[i] << "\n";
    else
      ss << get_elements()[i] << ",\n";      
  ss << "]\n";
}

template<typename parameters_t>
void frequency_domain_real_axis::initialize(parameters_t& parameters)
{
  int N = parameters.get_number_of_real_frequencies();

  double min = parameters.get_min_real_frequency();
  double max = parameters.get_max_real_frequency();

  double delta = (max-min)/double(N-1);

  get_elements().resize(N, min);

  for(int l=0; l<get_size(); l++)
    get_elements()[l] += double(l)*delta;
}

#endif
