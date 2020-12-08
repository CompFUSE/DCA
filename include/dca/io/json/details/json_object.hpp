// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// JSON base object.

#ifndef DCA_IO_JSON_DETAILS_JSON_OBJECT_HPP
#define DCA_IO_JSON_DETAILS_JSON_OBJECT_HPP

#include <ostream>

namespace dca::io::details {

class JSONObject {
public:
  JSONObject() = default;
  virtual ~JSONObject() = default;

  virtual void write(std::ostream& stream, int ident) const = 0;
  virtual bool read(std::istream& stream) = 0;
};

}  // namespace dca::io::details

#endif  // DCA_IO_JSON_DETAILS_JSON_OBJECT_HPP
