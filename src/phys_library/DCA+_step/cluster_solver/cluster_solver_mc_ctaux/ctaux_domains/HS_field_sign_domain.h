//-*-C++-*-

#ifndef DCA_QMCI_HS_FIELD_SIGN_DOMAIN_H
#define DCA_QMCI_HS_FIELD_SIGN_DOMAIN_H

// namespace DCA
// {
//   namespace QMCI
//   {

enum    HS_field_sign {HS_FIELD_DN=-1, HS_FIELD_UP=1};
typedef HS_field_sign HS_field_sign_type;

/*!
 *   \author Peter Staar
 */
class HS_field_sign_domain
{
public:

  typedef HS_field_sign_type element_type;

public:

  static int                              get_size();
  static std::vector<HS_field_sign_type>& get_elements();

  static int to_coordinate(HS_field_sign sign);

private:

  static std::vector<HS_field_sign_type> initialize_elements();
};

int HS_field_sign_domain::get_size()
{
  return 2;
}

std::vector<HS_field_sign_type>& HS_field_sign_domain::get_elements()
{
  static std::vector<HS_field_sign_type> v = initialize_elements();
  return v;
}

std::vector<HS_field_sign_type> HS_field_sign_domain::initialize_elements()
{
  static std::vector<HS_field_sign_type> v(0);

  v.push_back(HS_FIELD_DN);
  v.push_back(HS_FIELD_UP);

  return v;
}


int HS_field_sign_domain::to_coordinate(HS_field_sign sign)
{
  switch(sign)
    {
    case HS_FIELD_DN:
      return 0;
      break;

    case HS_FIELD_UP:
      return 1;
      break;

    default:
      throw std::logic_error(__FUNCTION__);
    }
}


//   }

// }

#endif
