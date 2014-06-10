//-*-C++-*-

#ifndef DCA_ITERATION_DOMAIN_H
#define DCA_ITERATION_DOMAIN_H

/*!
 *      Author: Peter Staar
 */
class DCA_iteration_domain 
{
public:

  typedef int element_type;

public:

  static int&              get_size();
  static std::string       get_name();

  static std::vector<int>& get_elements();

  template<typename parameters_type>
  static void initialize(parameters_type& parameters);

  template<IO::FORMAT DATA_FORMAT>
  static void read(IO::reader<DATA_FORMAT>& reader);

  template<IO::FORMAT DATA_FORMAT>
  static void write(IO::writer<DATA_FORMAT>& writer);

  template<class stream_type>
  static void to_JSON(stream_type& ss);
  
private:

  static std::vector<int>& initialize_elements();
};

int& DCA_iteration_domain::get_size()  
{
  static int size = 0;
  return size;
}

std::string DCA_iteration_domain::get_name()
{
  static std::string name = "DCA-iteration-domain";
  return name;
}

std::vector<int>& DCA_iteration_domain::get_elements()
{
  static std::vector<int>& v = initialize_elements();
  return v;
}

template<typename parameters_type>
void DCA_iteration_domain::initialize(parameters_type& parameters)
{
  get_size() = parameters.get_DCA_iterations();
}

std::vector<int>& DCA_iteration_domain::initialize_elements()
{
  static std::vector<int> v(get_size());

  for(int i=0; i<get_size(); i++)
    v[i] = i;

  return v;
}

template<IO::FORMAT DATA_FORMAT>
void DCA_iteration_domain::write(IO::writer<DATA_FORMAT>& writer)
{
  writer.open_group(get_name());
  writer.execute("elements", get_elements());
  writer.close_group();
}

template<class stream_type>
void DCA_iteration_domain::to_JSON(stream_type& ss)
{
  ss << "\"DCA_iteration_domain\" : [\n";
    
  for(int i=0; i<get_size(); i++)
    if(i == get_size()-1)
      ss << get_elements()[i] << "\n";
    else
      ss << get_elements()[i] << ",\n";
  
  ss << "]\n";
}

#endif
