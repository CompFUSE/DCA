// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef DCA_IO_JSON_JSON_PARSER_JSON_LEAF_HPP
#define DCA_IO_JSON_JSON_PARSER_JSON_LEAF_HPP

#include <cstdlib>
#include <string>

namespace dca {
namespace io {
// dca::io::

class JSON_leaf {
public:
  JSON_leaf();
  JSON_leaf(int n);
  JSON_leaf(int n, int ld);

  ~JSON_leaf();

  template <typename ss>
  void print(ss& ss_obj);

  template <typename some_type>
  JSON_leaf& operator=(some_type rhs);

  JSON_leaf& operator=(bool rhs);
  JSON_leaf& operator=(char rhs);
  JSON_leaf& operator=(std::string rhs);
  JSON_leaf& operator=(int rhs);
  JSON_leaf& operator=(double rhs);

public:
  std::size_t N;
  std::size_t LD;

  bool* bool_ptr;
  char* char_ptr;
  int* int_ptr;
  double* double_ptr;
};

}  // io
}  // dca

#endif  //  DCA_IO_JSON_JSON_PARSER_JSON_LEAF_HPP
