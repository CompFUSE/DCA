// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Utilities for JSON reader/writer.

#include "dca/io/json/details/util.hpp"

#include <istream>
#include <stdexcept>

namespace dca::io::details {

std::string readQuotedString(std::istream& inp) {
  trimUntil(inp, '\"');
  const auto pos = inp.tellg();

  std::string result;
  std::getline(inp, result, '\"');

  if (inp.eof()) {
    throw(std::logic_error("Line " + findLine(inp, pos) + ": missing matching \""));
  }

  return result;
}

void skipUntil(std::istream& inp, char target) {
  char c;
  const auto pos = inp.tellg();

  while (inp.read(&c, 1)) {
    if (c == target)
      return;
  }

  throw(std::logic_error("Line " + findLine(inp, pos) + ": parsing ended while looking for \'" +
                         target + "\'"));
}

void trimUntil(std::istream& inp, char target) {
  trimSpaces(inp);
  const auto pos = inp.tellg();
  char c;
  inp.read(&c, 1);

  if (c == target) {
    return;
  }
  else {
    throw(std::logic_error("Line " + findLine(inp, pos) + ": expected next character \'" + target +
                           "\'"));
  }
}

void trimSpaces(std::istream& inp) {
  char c;
  while (true) {
    switch (inp.peek()) {
      case ' ':
      case '\n':
      case '\t':
        break;
      default:
        return;
    }

    inp.read(&c, 1);
  }
}

std::string_view trimSpaces(std::string_view s) {
  auto is_space = [](char c) { return c == ' ' || c == '\n' || c == '\t'; };

  std::size_t start = 0;
  while (start < s.size() && is_space(s[start]))
    ++start;

  std::size_t end = s.size();  // exclusive
  while (end > 0 && is_space(s[end - 1]))
    --end;

  return s.substr(start, end - start);
}

std::string findLine(std::istream& inp) {
  return findLine(inp, inp.tellg());
}

std::string findLine(std::istream& inp, const std::streampos& pos) {
  const auto original_pos = inp.tellg();
  inp.seekg(0, inp.beg);

  int line = 0;
  std::string tmp;

  while (inp.tellg() < pos) {
    ++line;
    std::getline(inp, tmp);
  }

  inp.seekg(original_pos);
  return std::to_string(line);
}

}  // namespace dca::io::details
