#include <string>
#include <sstream>
#include <map>
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

template <class T>
inline std::string vectorToString(const std::vector<std::vector<T>>& v) {
  std::ostringstream ss;
  ss << "[";
  for (size_t i = 0; i < v.size(); ++i) {
    if (i != 0)
      ss << ",";
    ss << vectorToString(v[i]);
  }
  ss << "]";
  return ss.str();
}

/** Useful for dumping possibly nested string to string maps.
 *  safety and more general implementation is an exercise for the reader.
 */
template <class CONTAINER>
std::string mapToString(CONTAINER const& container) {
  std::ostringstream oss;
  auto size = std::size(container);
  oss << "{ ";
  for (auto const& [key, value] : container) {
    if constexpr (std::is_same<std::decay_t<decltype(value)>, std::map<std::string, std::string>>::value)
      oss << '{' << key << ',' << mapToString(value) << (--size ? "}, " : "} ");
    else
      oss << '{' << key << ',' << value << (--size ? "}, " : "} ");
  }
  oss << "}\n";
  return oss.str();
}

}  // namespace dca
#endif
