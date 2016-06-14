//-*-C++-*-

#ifndef FERMIONIC_OCCUPATION_NUMBER_DOMAIN_H
#define FERMIONIC_OCCUPATION_NUMBER_DOMAIN_H

namespace DCA
{
  namespace EXACT_DIAGONALIZATION
  {
    template<typename b_dmn, typename s_dmn, typename r_dmn>
    class occupation_number_domain
    {
    public:

      typedef int                                           element_type;
      typedef occupation_number_domain<b_dmn, s_dmn, r_dmn> this_type;

      static int get_size()
      {
        return b_dmn::dmn_size()*s_dmn::dmn_size()*r_dmn::dmn_size()+1;
      }

      static std::vector<element_type>& get_elements()
      {
        static std::vector<element_type>& elements = initialize();
        return elements;
      }

    private:

      static std::vector<element_type>& initialize()
      {
        static std::vector<element_type> elements(0, get_size());

        for(int i=0; i<get_size(); i++)
          elements[i] = i;

        return elements;
      }
    };

  }

}

#endif

