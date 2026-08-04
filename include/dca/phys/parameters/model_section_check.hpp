// Copyright (C) 2026 ETH Zurich
// Copyright (C) 2026 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Tyler Sax (tylersax@gmail.com)
//
// Issue #300: surface model-section input mistakes at parse time.
//
// A DCA++ executable is compiled for exactly one model and reads exactly one "*-model" section.
// This standalone helper inspects the top-level groups of a parsed JSON input file (as 
// ChildGroupStatus records) and reacts to the three situations:
//   * Typo: one or more "*-model" sections are present but none is the one this executable
//     reads -- likely a misspelled section name or the wrong input file. Throws
//     std::invalid_argument, since the run would otherwise silently use default model parameters.
//   * Multi-model file: the executable's section was read, but other "*-model" sections
//     are also present (e.g. a shared input file used by several model builds). Warns, since this
//     is a legitimate pattern currently used in testing.
//   * No model section at all: some non-simulation uses of Parameters (domain/FFT setup) do not
//     need model parameters. Warns rather than throwing, so those uses keep working.
//
// Out: out (warnings are written here; defaults to std::cerr)

#ifndef DCA_PHYS_PARAMETERS_MODEL_SECTION_CHECK_HPP
#define DCA_PHYS_PARAMETERS_MODEL_SECTION_CHECK_HPP

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "dca/io/json/details/json_group.hpp"

namespace dca {
namespace phys {
namespace params {

inline void checkModelSections(
    const std::vector<dca::io::details::ChildGroupStatus>& top_level_groups,
    const std::string& filename, std::ostream& out = std::cerr) {
  constexpr std::string_view suffix = "-model";

  std::vector<std::string> unused;
  bool model_section_read = false;

  for (const auto& [name, accessed] : top_level_groups) {
    const bool is_model_section =
        name.size() >= suffix.size() && std::equal(suffix.rbegin(), suffix.rend(), name.rbegin());
    if (!is_model_section)
      continue;

    if (accessed)
      model_section_read = true;
    else
      unused.push_back(name);
  }

  if (!model_section_read) {
    if (!unused.empty()) {
      // Model sections exist but none matches this executable's model -> typo / wrong file.
      std::string message = "Input '" + filename +
                            "' contains model section(s) [";
      for (std::size_t i = 0; i < unused.size(); ++i)
        message += (i ? ", " : "") + unused[i];
      message +=
          "] but none matches this executable's model. Check for a typo in the model section name.";
      throw std::invalid_argument(message);
    }
    // No model section at all: a non-simulation use of Parameters. Warn, but allow defaults.
    out << "Warning: input '" << filename
        << "' contains no model section; this executable's model will use default parameters.\n";
    return;
  }

  // The executable's model was read; flag any other model sections present but unused.
  for (const auto& name : unused)
    out << "Warning: input model section '" << name
        << "' is not used by this executable and will be ignored.\n";
}

}  // namespace params
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_PARAMETERS_MODEL_SECTION_CHECK_HPP
