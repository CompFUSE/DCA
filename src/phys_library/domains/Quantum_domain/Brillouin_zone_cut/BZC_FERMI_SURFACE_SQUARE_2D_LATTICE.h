//-*-C++-*-

#ifndef BRILLOUIN_ZONE_CUT_FERMI_SURFACE_SQUARE_2D_LATTICE_H
#define BRILLOUIN_ZONE_CUT_FERMI_SURFACE_SQUARE_2D_LATTICE_H

/*!
 *  \author Peter Staar
 */
template<>
class brillouin_zone_path_domain<FERMI_SURFACE_SQUARE_2D_LATTICE>
{
public:

  const static int INTERPOLATION_NB = 64;

  typedef std::vector<double>                                         element_type;
  typedef brillouin_zone_path_domain<FERMI_SURFACE_SQUARE_2D_LATTICE> this_type;

public:

  static std::string                get_name();  
  static int                        get_size();
  static std::vector<element_type>& get_elements();

private:

  static std::vector<element_type> initialize_elements();
};

std::string brillouin_zone_path_domain<FERMI_SURFACE_SQUARE_2D_LATTICE>::get_name()  
{
  return std::string("FERMI_SURFACE_SQUARE_2D_LATTICE");
}

int brillouin_zone_path_domain<FERMI_SURFACE_SQUARE_2D_LATTICE>::get_size()  
{
  return get_elements().size();
}

std::vector<std::vector<double> >& brillouin_zone_path_domain<FERMI_SURFACE_SQUARE_2D_LATTICE>::get_elements()
{
  static std::vector<element_type> v = initialize_elements();
  return v;
}

std::vector<std::vector<double> > brillouin_zone_path_domain<FERMI_SURFACE_SQUARE_2D_LATTICE>::initialize_elements()
{
  std::vector<element_type> collection_k_vecs(0);

  std::vector<double> k0(2);
  std::vector<double> k1(2);

  k0[0] = M_PI; k0[1] = 0;
  k1[0] = 0;    k1[1] = M_PI;

  for(int l=0; l<INTERPOLATION_NB; l++){

    std::vector<double> k(2,0);
    
    k[0] = (1.-double(l)/double(INTERPOLATION_NB))*k0[0] + double(l)/double(INTERPOLATION_NB)*k1[0];
    k[1] = (1.-double(l)/double(INTERPOLATION_NB))*k0[1] + double(l)/double(INTERPOLATION_NB)*k1[1];
    
    collection_k_vecs.push_back(k);
  }

  collection_k_vecs.push_back(k1);
  
  return collection_k_vecs;
}

#endif
