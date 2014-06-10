//-*-C++-*-

#ifndef BRILLOUIN_ZONE_CUT_H
#define BRILLOUIN_ZONE_CUT_H

/*!
 *  \author Peter Staar
 */
template<BRILLOUIN_ZONE_CUT_TYPE brillouin_zone_cut_t>
class brillouin_zone_path_domain
{

public:

  template<class stream_type>
  static void to_JSON(stream_type& ss);

  template<class stream_type>
  static void to_JSON(stream_type&                      ss, 
		      std::string                       name,
		      std::vector<std::vector<double> > elements);
};

#include "BZC_SQUARE_2D_LATTICE.h"
#include "BZC_FERMI_SURFACE_SQUARE_2D_LATTICE.h"

template<BRILLOUIN_ZONE_CUT_TYPE brillouin_zone_cut_t>
template<class stream_type>
void brillouin_zone_path_domain<brillouin_zone_cut_t>::to_JSON(stream_type& ss)
{
  to_JSON(ss, 
	  brillouin_zone_path_domain<SQUARE_2D_LATTICE>::get_name(), 
	  brillouin_zone_path_domain<SQUARE_2D_LATTICE>::get_elements());
  
  ss << ",\n";
  to_JSON(ss, 
	  brillouin_zone_path_domain<FERMI_SURFACE_SQUARE_2D_LATTICE>::get_name(), 
	  brillouin_zone_path_domain<FERMI_SURFACE_SQUARE_2D_LATTICE>::get_elements());
}

template<BRILLOUIN_ZONE_CUT_TYPE brillouin_zone_cut_t>
template<class stream_type>
void brillouin_zone_path_domain<brillouin_zone_cut_t>::to_JSON(stream_type&                      ss, 
							       std::string                       name,
							       std::vector<std::vector<double> > elements)
{
  ss << "\"" << name << "\" : [\n";
  
  for(size_t i=0; i<elements.size(); i++){
    
    ss << "[ ";
    for(size_t z=0; z<elements[i].size(); z++){
      if(z == elements[i].size()-1)
	ss << elements[i][z] << "]";
      else
	ss << elements[i][z] << ", ";
    }
    
    if(i == elements.size()-1)
      ss << "\n";
    else
      ss << ",\n";
  }
  ss << "]\n";
}

#endif
