//-*-C++-*-

#ifndef BASIS_TRANSFORMATIONS_DMN_2_TO_DMN_0_H
#define BASIS_TRANSFORMATIONS_DMN_2_TO_DMN_0_H

namespace TRAFOR
{
  /*!
   *  \class   TRANSFORM
   *  \ingroup TRANSFORM
   *
   *  \author  Peter Staar
   *  \brief   ...
   */
  template<typename type_input, typename type_output_0, typename type_output_1>
  struct TRANSFORM<dmn_0<type_input>, dmn_2<dmn_0<type_output_0>, dmn_0<type_output_1> > >
  {};

}

#endif
