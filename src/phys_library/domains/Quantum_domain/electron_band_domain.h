//-*-C++-*-

#ifndef ELECTRON_BAND_DOMAIN_H
#define ELECTRON_BAND_DOMAIN_H
#include <string>
#include <vector>
#include <assert.h>
#include "comp_library/IO_library/template_writer.h"
#include "comp_library/IO_library/template_reader.h"
/*!
 *  \author Peter Staar
 */
struct band_element
{
  band_element();
  ~band_element();

  int                 number;
  int                 flavor;
  std::vector<double> a_vec;
};

band_element::band_element():
  number(-1),
  flavor(-1),
  a_vec(0,0)
{}

band_element::~band_element()
{}

/*!
 *  \author Peter Staar
 */
int n_bands_defined=0;
class electron_band_domain 
{
public:

  typedef band_element element_type;

public:

  static int&                       get_size();
  static std::string                get_name();

  static std::vector<element_type>& get_elements();

  template<IO::FORMAT DATA_FORMAT>
  static void read(IO::reader<DATA_FORMAT>& reader);

  template<IO::FORMAT DATA_FORMAT>
  static void write(IO::writer<DATA_FORMAT>& writer);

  template<typename parameters_type>
  static void initialize(parameters_type& parameters, int Nb_bands, std::vector<int> flavors, std::vector<std::vector<double> > a_vecs);
  
  template<class stream_type>
  static void to_JSON(stream_type& ss);

private:

  static std::vector<element_type>& initialize_elements();
};

int& electron_band_domain::get_size()  
{
  static int size = 0;
  return size;
}

std::string electron_band_domain::get_name()  
{
  static std::string name = "electron-band-domain";
  return name;
}

std::vector<electron_band_domain::element_type>& electron_band_domain::get_elements()
{
  static std::vector<element_type> elements(get_size());
  return elements;
}

template<IO::FORMAT DATA_FORMAT>
void electron_band_domain::write(IO::writer<DATA_FORMAT>& writer)
{
  writer.open_group(get_name());
  writer.execute(get_elements());
  writer.close_group();
}

template<typename parameters_type>
void electron_band_domain::initialize(parameters_type& /*parameters*/, int NB_BANDS, std::vector<int> flavors, std::vector<std::vector<double> > a_vecs)
{
  get_size() = NB_BANDS;

  assert(NB_BANDS == int(flavors.size()));
  assert(NB_BANDS == int(a_vecs.size()));

  for (size_t i = 0; i < get_elements().size(); ++i) {

    get_elements()[i].number = i;
    get_elements()[i].flavor = i;//flavors[i];
    get_elements()[i].a_vec = a_vecs[i];
  }
}

#endif
