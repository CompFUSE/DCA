// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W. Doak (doakpw@ornl.gov)
//
// Type parameter composition type
//

#ifndef DCA_PHYS_PARAMETERS_NUM_TRAITS_HPP
#define DCA_PHYS_PARAMETERS_NUM_TRAITS_HPP

namespace dca {

template <typename REALPREC, typename SCALAR>
struct NumericalTraits {
  using Real = REALPREC;
  using Scalar = SCALAR;
  using TPAccumPrec = REALPREC;
};

}
#endif
