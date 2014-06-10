//-*-C++-*-

#ifndef FEYNMAN_EXPANSION_ORDER_DOMAIN_H
#define FEYNMAN_EXPANSION_ORDER_DOMAIN_H

namespace DCA
{
  namespace QMCI
  {
    /*!
     *      Author: Peter Staar
     */
    class Feynman_expansion_order_domain
    {
      const static int MAX_ORDER_SQUARED = 10000;

    public:

      typedef int element_type;

    public:

      static int               get_size();
      static std::vector<int>& get_elements();

    private:

      static std::vector<int>& initialize_elements();
    };

    int Feynman_expansion_order_domain::get_size()
    {
      const static int size = MAX_ORDER_SQUARED;
      return size;
    }

    std::vector<int>& Feynman_expansion_order_domain::get_elements()
    {
      static std::vector<int>& v = initialize_elements();
      return v;
    }

    std::vector<int>& Feynman_expansion_order_domain::initialize_elements()
    {
      static std::vector<int> v(get_size());

      for(int i=0; i<get_size(); i++)
        v[i] = i;

      return v;
    }

  }

}

#endif
