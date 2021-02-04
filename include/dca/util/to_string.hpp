#include <string>
#include <sstream>

#ifndef DCA_TO_STRING_HPP
#define DCA_TO_STRING_HPP

namespace dca {

template <class T>
inline std::string vectorToString(const std::vector<T>& v) {
  std::ostringstream ss;
  ss << "[";
  for (size_t i = 0; i < v.size(); ++i) {
    if (i != 0)
      ss << ",";
    ss << v[i];
  }
  ss << "]";
  return ss.str();
}

}

#endif
