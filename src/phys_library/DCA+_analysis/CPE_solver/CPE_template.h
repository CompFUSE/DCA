//-*-C++-*-

#ifndef CONTINUOUS_POLE_EXPANSION_TEMPLATE_H
#define CONTINUOUS_POLE_EXPANSION_TEMPLATE_H

namespace DCA
{
  /*!
   *  \defgroup CPE-BASIS-FUNCTIONS
   *  \ingroup  CPE
   */

  /*! \class   continuous_pole_expansion
   *  \ingroup CPE
   *
   *  \author  Peter Staar
   *  \brief   Empty template class for a CPE analytic continuation.
   *  \version 1.0
   */
  template<class                    parameters_type,
           class                    basis_function_t,
	   class                    k_dmn_t,
	   class                    w_dmn_t,	   
           minimization_method_type minimization_method_t>
  class continuous_pole_expansion
  {};

}

#endif
