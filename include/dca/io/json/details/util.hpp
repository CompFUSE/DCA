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

#ifndef DCA_IO_JSON_DETAILS_UTIL_HPP
#define DCA_IO_JSON_DETAILS_UTIL_HPP

#include <string>

namespace dca::io::details {

std::string readQuotedString(std::istream& inp);

void skipUntil(std::istream& inp, char target);
void trimUntil(std::istream& inp, char target);

void trimSpaces(std::istream& inp);
std::string_view trimSpaces(std::string_view s);

std::string findLine(std::istream& inp);
std::string findLine(std::istream& inp, const std::streampos& pos);

}  // namespace dca::io::details

#endif  // DCA_IO_JSON_DETAILS_UTIL_HPP
