#include <string>
#include <sstream>

namespace dca {

template <class T>
std::string VectorToString(const std::vector<T>& v) {
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
