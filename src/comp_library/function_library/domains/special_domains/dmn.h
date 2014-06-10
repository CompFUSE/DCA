//-*-C++-*-

#ifndef DMN_H
#define DMN_H

/*! 
 *  \class dmn
 *  \ingroup FUNCTION
 *
 *  \brief concept for each type of domain
 *  \author Peter Staar
 *
 *  \version 1.0
 *  \date    2009-2011
 */
template<int size, class element_t = int>
class dmn
{
public:

  typedef element_t            element_type;
  typedef dmn<size, element_t> this_type;

  static int get_size();
  static std::vector<element_t>& get_elements();
};

template<int size, class element_t>
int dmn<size, element_t>::get_size()
{
  return size;
}

template<int size, class element_t>
std::vector<element_t>& dmn<size, element_t>::get_elements()
{
  static std::vector<element_t> v(size);
  return v;
}

#endif

