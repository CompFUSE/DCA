// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Data entry in a JSON file.

#ifndef DCA_IO_JSON_DETAILS_JSON_ENTRY_HPP
#define DCA_IO_JSON_DETAILS_JSON_ENTRY_HPP

#include "dca/io/json/details/json_object.hpp"

#include <cassert>
#include <sstream>
#include <string>
#include <vector>

#include "dca/io/json/details/convert.hpp"

namespace dca::io::details {

class JSONEntry : public JSONObject {
public:
  JSONEntry() = default;

  template <class T>
  JSONEntry(const T& val) {
    std::stringstream ss;
    ss << val;
    data_ = ss.str();
  }

  ~JSONEntry() = default;

  JSONEntry(const std::string& val) {
    data_ = '\"' + val + '\"';
  }
  JSONEntry(const char* val) {
    data_ = '\"' + std::string(val) + '\"';
  }

  void write(std::ostream& stream, int /*ident*/) const override {
    stream << data_;
  }

  // Returns the success of the write operation. A common failure is a type not compatible with the
  // string representation. In case of failure "obj: is unmodified.
  template <class T>
  bool write(T& obj) const {
    return Convert<T>::execute(data_, obj);
  }

  // Return true if this is the last element of a group.
  bool read(std::istream& inp) override;

private:
  std::string data_;
};

}  // namespace dca::io::details

#endif  // DCA_IO_JSON_DETAILS_JSON_GROUP_HPP
