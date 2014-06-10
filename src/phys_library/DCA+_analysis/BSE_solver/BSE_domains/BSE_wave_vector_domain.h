//-*-C++-*-

#ifndef BSE_WAVE_VECTOR_DOMAIN_H
#define BSE_WAVE_VECTOR_DOMAIN_H

namespace DCA
{
  /*!
   *  \author: Peter Staar
   */
  template<typename r_dmn_type>
  class wave_vector_domain
  {
  public:

    typedef double                   scalar_type;
    typedef std::vector<scalar_type> element_type;

  public:

    static int&         get_size();
    static std::string  get_name();

    static std::vector<element_type>& get_elements();

    template<typename parameters_type>
    static void initialize(parameters_type& parameters);

    template<IO::FORMAT DATA_FORMAT>
    static void write(IO::writer<DATA_FORMAT>& writer);
  };

  template<typename r_dmn_type>
  int& wave_vector_domain<r_dmn_type>::get_size()
  {
    static int size = 0;
    return size;
  }

  template<typename r_dmn_type>
  std::string wave_vector_domain<r_dmn_type>::get_name()
  {
    static std::string name = "wave_vector_domain";
    return name;
  }

  template<typename r_dmn_type>
  std::vector<typename wave_vector_domain<r_dmn_type>::element_type>& wave_vector_domain<r_dmn_type>::get_elements()
  {
    static std::vector<element_type> elements(0);
    return elements;
  }

  template<typename r_dmn_type>
  template<typename parameters_type>
  void wave_vector_domain<r_dmn_type>::initialize(parameters_type& parameters)
  {
    std::vector<element_type>& elements = get_elements();

    elements origin(r_dmn_type::parameter_type::DIMENSION, 0.);
    for(int l=0; l<r_dmn_type::dmn_size(); l++)
      {
	double dist = cluster_operations::minimal_distance(origin, 
							   r_dmn_type::parameter_type::get_elements()[l], 
							   r_dmn_type::parameter_type::get_basis_vectors());

	if(dist < parameters.get_BSE_radius())
	  elements.push_back(r_dmn_type::parameter_type::get_elements()[l]);
      }

    get_size() = elements.size();
  }

  template<typename r_dmn_type>
  void wave_vector_domain<r_dmn_type>::write(IO::writer<DATA_FORMAT>& writer)
  {
    writer.open_group(get_name());
    writer.execute("elements", get_elements());
    writer.close_group();
  }

}

#endif
