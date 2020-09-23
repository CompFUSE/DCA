// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Conversion from and to human readable streams and objects.

#ifndef DCA_IO_JSON_DETAILS_CONVERT_HPP
#define DCA_IO_JSON_DETAILS_CONVERT_HPP

#include <array>
#include <string>
#include <sstream>
#include <vector>

namespace dca::io {

template <class T>
std::istream& operator>>(std::istream& stream, std::vector<T>& vec);

template <class T>
struct Convert {
  static T execute(const std::string& val) {
    T ret;
    std::stringstream stream(val);
    stream >> ret;
    return ret;
  }
};

template <>
struct Convert<std::string> {
  static std::string execute(const std::string& val) {
    if (val.size() > 1 && (val[0] == '\"' && val.back() == '\"'))
      return val.substr(1, val.size() - 2);
    else
      throw(std::logic_error("Not a string"));
  }
};

template <>
struct Convert<bool> {
  static bool execute(const std::string& val) {
    if (val.size() >= 4 && val.substr(0, 4) == "true")
      return true;
    else if (val.size() >= 5 && val.substr(0, 5) == "false")
      return false;
    else
      return std::stoi(val);
  }
};

template <class T1, class T2>
struct Convert<std::pair<T1, T2>> {
  static auto execute(const std::string& val) {
    static_assert(std::is_scalar_v<T1> && std::is_scalar_v<T2>,
                  "composite pair members are not supported");

    const auto split = val.find(',');
    if ((val.at(0) != '(' && val.back() != ')') || split == std::string::npos)
      throw(std::logic_error("Not a pair"));

    std::pair<T1, T2> p;
    p.first = Convert<T1>::execute(val.substr(1, split - 1));
    p.second = Convert<T2>::execute(val.substr(split + 1, val.size() - split - 2));

    return p;
  }
};

template <class T1, class T2>
std::ostream& operator<<(std::ostream& stream, const std::pair<T1, T2>& p) {
  return stream << '(' << p.first << ", " << p.second << ')';
}

template <class T>
std::ostream& operator<<(std::ostream& stream, const std::vector<T>& vec) {
  stream << "[";
  for (std::size_t i = 0; i < vec.size(); ++i) {
    if constexpr (std::is_same_v<T, std::string>) {
      stream << "\"" + vec[i] + "\"";
    }
    else {
      stream << vec[i];
    }
    if (i < vec.size() - 1)
      stream << ", ";
  }
  stream << "]";
  return stream;
}

template <class T, std::size_t n>
std::ostream& operator<<(std::ostream& stream, const std::array<T, n>& arr) {
  return stream << std::vector<T>(arr.begin(), arr.end());
}

template <class T, std::size_t n>
std::istream& operator>>(std::istream& stream, std::array<T, n>& arr) {
  std::vector<T> v;
  v.reserve(n);

  stream >> v;
  std::copy(v.begin(), v.end(), arr.begin());
  return stream;
}

template <class T, std::size_t n>
std::istream& operator>>(std::istream& stream, std::vector<std::array<T, n>>& v) {
  std::vector<std::vector<T>> vv;

  stream >> vv;
  v.resize(vv.size());
  for (std::size_t i = 0; i < vv.size(); ++i) {
    std::copy(vv[i].begin(), vv[i].end(), v[i].begin());
  }
  return stream;
}

template <class T>
std::istream& operator>>(std::istream& stream, std::vector<T>& vec) {
  vec.clear();

  char c;
  stream.read(&c, 1);
  if (c != '[') {
    stream.seekg(-1, stream.cur);
    throw(std::logic_error("invalid vector format"));
  }

  std::string value;
  bool quote = false;
  int parentheses = 0;

  while (stream.read(&c, 1)) {
    if (!quote) {
      switch (c) {
        case ']':
          if (value != "")
            vec.push_back(Convert<T>::execute(value));
          return stream;
        case '(':
          ++parentheses;
          value.push_back(c);
          break;
        case ')':
          --parentheses;
          value.push_back(c);
          break;
        case ',':
          if (parentheses == 0) {
            vec.push_back(Convert<T>::execute(value));
            value.clear();
          }
          else {
            value.push_back(c);
          }
          break;
        case '\"':
          quote = true;
          value.push_back(c);
          break;
        case ' ':
          break;
        default:
          value.push_back(c);
      }
    }
    else {
      if (c == '\"')
        quote = false;
      value.push_back(c);
    }
  }

  return stream;
}

template <class T>
std::istream& operator>>(std::istream& stream, std::vector<std::vector<T>>& vec) {
  vec.clear();

  char c;
  stream.read(&c, 1);
  if (c != '[') {
    stream.seekg(-1, stream.cur);
    throw(std::logic_error("invalid vector format"));
  }

  while (true) {
    vec.emplace_back();
    stream >> vec.back();

    // trim whitespaces.
    while (stream.peek() == ' ')
      stream.read(&c, 1);

    stream.read(&c, 1);
    if (c == ']')
      break;

    assert(c == ',');
    while (stream.peek() == ' ')
      stream.read(&c, 1);
  }

  return stream;
}

}  // namespace dca::io

#endif  // DCA_IO_JSON_DETAILS_CONVERT_HPP
