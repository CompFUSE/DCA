// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// JSON group.

#ifndef DCA_IO_JSON_DETAILS_JSON_GROUP_HPP
#define DCA_IO_JSON_DETAILS_JSON_GROUP_HPP

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "json_entry.hpp"
#include "json_object.hpp"

namespace dca::io::details {

class JSONGroup : public JSONObject {
public:
  JSONGroup() = default;
  ~JSONGroup() = default;

  JSONGroup* addGroup(const std::string& name);
  JSONGroup* getGroup(const std::string& name);

  template <class T>
  void addEntry(const std::string& name, const T& val) {
    const auto it = objects_.find(name);

    if (it == objects_.end()) {
      auto [new_it, _] = objects_.insert({name, std::make_unique<JSONEntry>(val)});
      order_.push_back(new_it);
    }
    else {
      it->second = std::make_unique<JSONEntry>(val);
    }
  }

  template <class T>
  bool readEntry(const std::string& name, T& val) const noexcept;

  void write(std::ostream& stream, int ident) const override;
  bool read(std::istream& inp) override;

  void clear() {
    objects_.clear();
  }

private:
  using Container = std::unordered_map<std::string, std::unique_ptr<JSONObject>>;
  Container objects_;
  std::vector<Container::const_iterator> order_;
};

template <class T>
bool JSONGroup::readEntry(const std::string& name, T& val) const noexcept {
  JSONEntry* entry = nullptr;
  auto it = objects_.find(name);
  if (it != objects_.end()) {
    entry = dynamic_cast<JSONEntry*>(it->second.get());
  }

  if (!entry)
    return false;

  return entry->write(val);
}

}  // namespace dca::io::details

#endif  // DCA_IO_JSON_DETAILS_JSON_GROUP_HPP
