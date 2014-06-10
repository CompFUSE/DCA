//-*-C++-*-                                                                                                                                                                                                                                                                         

#ifndef LIN_ALG_MATRIX_SCALARTYPE_H
#define LIN_ALG_MATRIX_SCALARTYPE_H

namespace LIN_ALG {

  template<typename scalartype, device_type device_name>
  struct MATRIX_SCALARTYPE
  {
    typedef scalartype new_scalartype;
  };

}

#endif
