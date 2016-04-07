//-*-C++-*-

#ifndef ELECTRON_SPIN_DOMAIN_H
#define ELECTRON_SPIN_DOMAIN_H
#include "comp_library/IO_library/HDF5/HDF5_writer.h"

enum    e_spin_states {e_DN=-1, e_UP=1};
typedef e_spin_states e_spin_states_type;

/*
 *      Author: Peter Staar
 */
class electron_spin_domain
{
public:

  typedef e_spin_states_type element_type;

public:

  static int         get_size();
  static std::string get_name();

  static std::vector<e_spin_states_type>& get_elements();

  template<IO::FORMAT DATA_FORMAT>
  static void write(IO::writer<DATA_FORMAT>& writer);

  static int to_coordinate(element_type spin);

private:

  static std::vector<e_spin_states_type> initialize_elements();
};

int electron_spin_domain::get_size()
{
  return 2;
}

std::string electron_spin_domain::get_name()
{
  static std::string name = "electron-spin-domain";
  return name;
}

std::vector<e_spin_states_type>& electron_spin_domain::get_elements()
{
  static std::vector<e_spin_states_type> v = initialize_elements();
  return v;
}

template<IO::FORMAT DATA_FORMAT>
void electron_spin_domain::write(IO::writer<DATA_FORMAT>& writer)
{
  writer.open_group(get_name());
  writer.execute("elements", get_elements());
  writer.close_group();
}

std::vector<e_spin_states_type> electron_spin_domain::initialize_elements()
{
  static std::vector<e_spin_states_type> v(0);

  v.push_back(e_DN);
  v.push_back(e_UP);

  return v;
}

int electron_spin_domain::to_coordinate(element_type spin)
{
  switch(spin)
    {
    case e_DN:
      return 0;
      break;

    case e_UP:
      return 1;
      break;

    default:
      throw std::logic_error(__FUNCTION__);
    }
}

#endif
