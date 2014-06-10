//-*-C++-*-

#ifndef DCA_QMCI_HS_VERTEX_MOVE_DOMAIN_H
#define DCA_QMCI_HS_VERTEX_MOVE_DOMAIN_H

// namespace DCA
// {
//   namespace QMCI
//   {

enum    HS_vertex_move {ANNIHILATION=-1, STATIC=0, CREATION=1};
typedef HS_vertex_move HS_vertex_move_type;

/*!
 *  \author: Peter Staar
 */
class HS_vertex_move_domain
{
public:

  typedef HS_vertex_move_type element_type;

public:

  static int                               get_size();
  static std::vector<HS_vertex_move_type>& get_elements();

  static int to_coordinate(element_type vertex_move);

private:

  static std::vector<HS_vertex_move_type> initialize_elements();
};

int HS_vertex_move_domain::get_size()
{
  return 3;
}

std::vector<HS_vertex_move_type>& HS_vertex_move_domain::get_elements()
{
  static std::vector<HS_vertex_move_type> v = initialize_elements();
  return v;
}

std::vector<HS_vertex_move_type> HS_vertex_move_domain::initialize_elements()
{
  static std::vector<HS_vertex_move_type> v(0);

  v.push_back(ANNIHILATION);
  v.push_back(STATIC);
  v.push_back(CREATION);

  return v;
}

int HS_vertex_move_domain::to_coordinate(element_type vertex_move)
{
  switch(vertex_move)
    {
    case ANNIHILATION:
      return 0;
      break;

    case STATIC:
      return 1;
      break;

    case CREATION:
      return 2;
      break;

    default:
      throw std::logic_error(__FUNCTION__);
    }
}

//   }

// }

#endif
