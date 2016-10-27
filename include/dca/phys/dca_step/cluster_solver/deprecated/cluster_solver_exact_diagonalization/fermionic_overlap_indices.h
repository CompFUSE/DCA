//-*-C++-*-

#ifndef FERMIONIC_OVERLAP_INDICES_DOMAIN_H
#define FERMIONIC_OVERLAP_INDICES_DOMAIN_H

namespace DCA
{
  namespace EXACT_DIAGONALIZATION
  {
    struct overlap_indices
    {
    public:

      int index;

      int lhs;
      int rhs;

      double sign;

    public:

      friend bool operator<(overlap_indices const& lhs,
                            overlap_indices const& rhs)
      {
        if(lhs.index!=rhs.index)
          return lhs.index<rhs.index;

        if(lhs.lhs!=rhs.lhs)
          return lhs.lhs<rhs.lhs;

        return lhs.rhs<rhs.rhs;
      }

    };

  }

}

#endif
