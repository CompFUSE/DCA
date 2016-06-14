//-*-C++-*-

#ifndef FERMIONIC_MAGNETIZATION_DOMAIN_H
#define FERMIONIC_MAGNETIZATION_DOMAIN_H

namespace DCA
{
  namespace EXACT_DIAGONALIZATION
  {
    template<typename b_dmn, typename s_dmn, typename r_dmn>
    class magnetization_domain
    {
    public:

      typedef int                                       element_type;
      typedef magnetization_domain<b_dmn, s_dmn, r_dmn> this_type;

      static int get_size()
      {
        return b_dmn::dmn_size()*s_dmn::dmn_size()*r_dmn::dmn_size()+1;
      }

      static std::vector<element_type>& get_elements()
      {
        static std::vector<element_type>& elements = initialize();
        return elements;
      }

      static int Sz(int i)
      {
        return i-b_dmn::dmn_size()*r_dmn::dmn_size();
      }

      static int*& get_magnetization()
      {
        static int*& ptr = initialize_magnetization();
        return ptr;
      }

    private:

      static std::vector<element_type>& initialize()
      {
        static std::vector<element_type> elements(0, get_size());

        for(int i=0; i<get_size(); i++)
          elements[i] = i-b_dmn::dmn_size()*r_dmn::dmn_size();

        return elements;
      }

      static int*& initialize_magnetization()
      {
        int max_occ = b_dmn::dmn_size()*s_dmn::dmn_size()*r_dmn::dmn_size();
        static int* magnetization = new int[max_occ];

        for(int r=0; r<r_dmn::dmn_size(); ++r)
          for(int i=0; i<s_dmn::dmn_size(); ++i)
            for(int j=0; j<b_dmn::dmn_size(); ++j)
              magnetization[j+b_dmn::dmn_size()*(i+r*s_dmn::dmn_size())] = i==0? 1 : -1;

        return magnetization;
      }
    };

  }

}

#endif
