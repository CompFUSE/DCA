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

#ifndef DCA_IO_JSON_JSON_PARSER_JSON_TREE_HPP
#define DCA_IO_JSON_JSON_PARSER_JSON_TREE_HPP

#include <cstdlib>
#include <map>
#include <string>
#include <vector>

#include "dca/io/json/json_parser/json_leaf.hpp"

namespace dca {
namespace io {
// dca::io::

class JSON_tree {
public:
  typedef std::map<std::string, JSON_leaf> map_JSON_leaf_type;
  typedef std::map<std::string, JSON_tree> map_JSON_tree_type;

  typedef std::vector<JSON_tree> vec_arbitrary_JSON_type;

public:
  JSON_tree() : key("no-key"), index(0), parent(NULL), JSON_leaf_obj(), JSON_tree_obj(), vec(0) {}

  std::string& name() {
    return key;
  }

  JSON_tree& operator[](std::string& key) {
    return JSON_tree_obj[key];
  }

  template <typename val_t>
  void set(std::string& key, val_t& val) {
    JSON_leaf_obj[key] = val;
  }

  template <typename val_t>
  void set(std::string& key, std::vector<val_t>& val) {
    JSON_leaf_obj[key] = val;
  }

private:
  std::string key;
  std::size_t index;

  JSON_tree* parent;

  map_JSON_leaf_type JSON_leaf_obj;
  map_JSON_tree_type JSON_tree_obj;

  vec_arbitrary_JSON_type vec;
};

}  // io
}  // dca

#endif  //  DCA_IO_JSON_JSON_PARSER_JSON_TREE_HPP
