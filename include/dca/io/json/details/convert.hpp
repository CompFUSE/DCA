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

#include "dca/io/json/details/util.hpp"

namespace dca::io {

template <class T>
struct Convert {
  static bool execute(const std::string& inp, T& val) {
    std::stringstream stream(inp);
    try {
      stream >> val;
    }
    catch (...) {
      return false;
    }
    return true;
  }
};

template <>
struct Convert<std::string> {
  static bool execute(std::string_view inp, std::string& val) {
    inp = details::trimSpaces(inp);
    if (inp.size() >= 2 && inp[0] == '\"' && inp.back() == '\"') {
      val = inp.substr(1, inp.size() - 2);
      return true;
    }
    else {
      return false;
    }
  }
};

template <>
struct Convert<bool> {
  static bool execute(std::string_view inp, bool& val) {
    inp = details::trimSpaces(inp);

    if (inp == "true")
      val = true;
    else if (inp == "false")
      val = false;
    else
      val = std::atoi(inp.data());

    return true;
  }
};

template <class T1, class T2>
struct Convert<std::pair<T1, T2>> {
  static bool execute(const std::string& inp, std::pair<T1, T2>& p) {
    static_assert(std::is_scalar_v<T1> && std::is_scalar_v<T2>,
                  "composite pair members are not supported");

    const auto split = inp.find(',');
    if ((inp.at(0) != '(' && inp.back() != ')') || split == std::string::npos) {
      // Not a pair
      return false;
    }

    std::pair<T1, T2> result;
    bool success = true;

    success &= Convert<T1>::execute(inp.substr(1, split - 1), result.first);
    success &= Convert<T2>::execute(inp.substr(split + 1, inp.size() - split - 2), result.second);

    if (success) {
      p = std::move(result);
      return true;
    }
    else {
      return false;
    }
  }
};

template <class T, std::size_t n>
struct Convert<std::array<T, n>> {
  static bool execute(const std::string& inp, std::array<T, n>& arr) {
    std::vector<T> v;
    v.reserve(n);

    // Expansive but easy: delegate to vector conversion.
    const bool success = Convert<std::vector<T>>::execute(inp, v);
    if (!success || v.size() != n)
      return false;

    std::copy(v.begin(), v.end(), arr.begin());
    return true;
  }
};

template <class T>
struct Convert<std::vector<T>> {
  static bool execute(std::string_view inp, std::vector<T>& vec) {
    inp = details::trimSpaces(inp);

    if (inp.size() < 2 || inp[0] != '[' || inp.back() != ']') {
      // Invalid vector format.
      return false;
    }

    std::vector<T> result;  // move result into vec after it has been successfully read.

    bool quote = false;
    int parentheses = 0;
    std::size_t pos = 1;
    std::string entry;

    auto add_element = [&]() {
      bool success = true;
      if (entry != "") {
        result.emplace_back();
        success = Convert<T>::execute(entry, result.back());
        entry.clear();
      }
      return success;
    };

    while (pos < inp.size() - 1) {  // skip first and last square bracket.
      const char c = inp[pos++];

      if (!quote) {
        switch (c) {
          case '[':
          case '(':
            ++parentheses;
            entry.push_back(c);
            break;
          case ')':
          case ']':
            --parentheses;
            if (parentheses < 0) {  // imbalanced parentheses.
              return false;
            }
            entry.push_back(c);
            break;
          case ',':
            if (parentheses == 0) {
              if (!add_element())
                return false;
            }
            else {
              entry.push_back(c);
            }
            break;
          case '\"':
            quote = true;
            entry.push_back(c);
            break;
          default:
            entry.push_back(c);
        }
      }
      else {
        if (c == '\"')
          quote = false;
        entry.push_back(c);
      }
    }

    // Process last element.
    if (!add_element())
      return false;

    vec = std::move(result);
    return true;
  }
};

// Clang requires this forward declaration.
template <class T1, class T2>
std::ostream& operator<<(std::ostream& stream, const std::pair<T1, T2>& p);
template <class T>
std::ostream& operator<<(std::ostream& stream, const std::vector<T>& vec);
template <class T, std::size_t n>
std::ostream& operator<<(std::ostream& stream, const std::array<T, n>& arr);

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

}  // namespace dca::io

#endif  // DCA_IO_JSON_DETAILS_CONVERT_HPP
