// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W.  Doak (doakpw@ornl.gov)
//

/** \file
 *  This file provides some "global" types for dca++
 *
 *  These type definitions should be such that they don't force a particular precision,
 *  complex/real, or model regime onto the code using them.
 */
namespace dca {

template<typename Scalar>
using SignType = std::conditional_t<dca::util::IsComplex_t<Scalar>::value, Scalar, std::int8_t>;

}  
