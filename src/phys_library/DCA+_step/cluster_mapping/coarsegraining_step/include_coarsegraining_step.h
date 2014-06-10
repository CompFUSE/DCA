//-*-C++-*-

#ifndef DCA_COARSEGRAIN_DOMAIN_NAMES_H
#define DCA_COARSEGRAIN_DOMAIN_NAMES_H

namespace DCA
{
  enum COARSEGRAIN_DOMAIN_NAMES {ORIGIN, Q_FINE, K, K_PLUS_Q, Q_MINUS_K, TETRAHEDRON_K, TETRAHEDRON_ORIGIN};

  std::string to_str(COARSEGRAIN_DOMAIN_NAMES NAME)
  {
    switch(NAME)
      {
      case ORIGIN:
	return "ORIGIN";

      case Q_FINE:
	return "Q-FINE";

      case K:
	return "K";

      case K_PLUS_Q:
	return "K+Q";

      case Q_MINUS_K:
	return "Q-K";

      case TETRAHEDRON_K:
	return "TETRAHEDRON-K";

      case TETRAHEDRON_ORIGIN:
	return "TETRAHEDRON-ORIGIN";

      default:
	throw std::logic_error(__FUNCTION__);
      }

    return "???";
  }

}

#endif

#include "coarsegraining_domain.h"

#include "coarsegraining_interpolation_matrices.h"

#include "quadrature_integration.h"

#include "tetrahedron_data.h"
#include "tetrahedron_routines_harmonic_function.h"
#include "tetrahedron_routines_inverse_matrix_function.h"
#include "tetrahedron_integration.h"

#include "coarsegraining_routines.h"

#include "coarsegraining_sp.h"
#include "coarsegraining_tp.h"


