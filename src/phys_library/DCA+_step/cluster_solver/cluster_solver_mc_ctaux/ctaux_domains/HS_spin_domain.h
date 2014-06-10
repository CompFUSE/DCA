//-*-C++-*-

#ifndef DCA_QMCI_HS_SPIN_DOMAIN_H
#define DCA_QMCI_HS_SPIN_DOMAIN_H

// namespace DCA
// {
//   namespace QMCI
//   {

enum    HS_spin_states {HS_DN=-1, HS_ZERO=0, HS_UP=1};
typedef HS_spin_states HS_spin_states_type;

/*!
 *  \author: Peter Staar
 */
class HS_spin_domain
{
public:

  typedef HS_spin_states_type element_type;

public:

  static int                               get_size();
  static std::vector<HS_spin_states_type>& get_elements();

  static int to_coordinate(element_type spin);

private:

  static std::vector<HS_spin_states_type> initialize_elements();
};

int HS_spin_domain::get_size()
{
  return 3;
}

std::vector<HS_spin_states_type>& HS_spin_domain::get_elements()
{
  static std::vector<HS_spin_states_type> v = initialize_elements();
  return v;
}

std::vector<HS_spin_states_type> HS_spin_domain::initialize_elements()
{
  static std::vector<HS_spin_states_type> v(0);

  v.push_back(HS_DN);
  v.push_back(HS_ZERO);
  v.push_back(HS_UP);

  return v;
}

int HS_spin_domain::to_coordinate(element_type spin)
{
  switch(spin)
    {
    case HS_DN:
      return 0;
      break;

    case HS_ZERO:
      return 1;
      break;

    case HS_UP:
      return 2;
      break;

    default:
      throw std::logic_error(__FUNCTION__);
    }
}

//   }

// }

#endif
