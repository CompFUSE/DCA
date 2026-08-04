// Copyright (C) 2026 ETH Zurich
// Copyright (C) 2026 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Tyler Sax (tylersax@gmail.com)
//
// Unit tests for checkModelSections (issue #300 bugs #2 and #3).

#include "dca/phys/parameters/model_section_check.hpp"

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "dca/testing/gtest_h_w_warning_blocking.h"

using dca::phys::params::checkModelSections;
using Groups = std::vector<dca::io::details::ChildGroupStatus>;

// The built model's section was read and nothing else is present: no warning, no throw.
TEST(CheckModelSectionsTest, OnlyBuiltModelRead) {
  const Groups groups{{"domains", true}, {"single-band-Hubbard-model", true}};
  std::ostringstream warnings;
  EXPECT_NO_THROW(checkModelSections(groups, "input.json", warnings));
  EXPECT_TRUE(warnings.str().empty());
}

// The built model was read, but other model sections are present and unused -> warn only.
TEST(CheckModelSectionsTest, UnusedModelSectionWarns) {
  const Groups groups{{"single-band-Hubbard-model", true}, {"bilayer-Hubbard-model", false}};
  std::ostringstream warnings;
  EXPECT_NO_THROW(checkModelSections(groups, "input.json", warnings));
  const std::string out = warnings.str();
  EXPECT_NE(out.find("bilayer-Hubbard-model"), std::string::npos);
  EXPECT_NE(out.find("ignored"), std::string::npos);
  // The section that was actually used must not be warned about.
  EXPECT_EQ(out.find("single-band-Hubbard-model"), std::string::npos);
}

// Model section(s) present but none matches the built model (typo / wrong file)
// -> throw, even though a model section exists.
TEST(CheckModelSectionsTest, WrongModelSectionThrows) {
  const Groups groups{{"bilayer-Hubbard-model", false}, {"domains", true}};
  std::ostringstream warnings;
  try {
    checkModelSections(groups, "input.json", warnings);
    FAIL() << "expected std::invalid_argument";
  }
  catch (const std::invalid_argument& e) {
    const std::string msg = e.what();
    EXPECT_NE(msg.find("input.json"), std::string::npos);
    // The unmatched section is surfaced as a typo hint.
    EXPECT_NE(msg.find("bilayer-Hubbard-model"), std::string::npos);
  }
  // A throw must take precedence over emitting warnings.
  EXPECT_TRUE(warnings.str().empty());
}

// No model section present at all -> warn (do not throw); model-agnostic uses keep working.
TEST(CheckModelSectionsTest, NoModelSectionPresentWarns) {
  const Groups groups{{"domains", true}, {"output", false}};
  std::ostringstream warnings;
  EXPECT_NO_THROW(checkModelSections(groups, "input.json", warnings));
  EXPECT_NE(warnings.str().find("no model section"), std::string::npos);
}

// Sections that merely contain "model" but do not end in "-model" are not model sections.
TEST(CheckModelSectionsTest, SuffixMatchIsExact) {
  const Groups groups{{"model-parameters", false}, {"single-band-Hubbard-model", true}};
  std::ostringstream warnings;
  EXPECT_NO_THROW(checkModelSections(groups, "input.json", warnings));
  EXPECT_TRUE(warnings.str().empty());  // "model-parameters" is not treated as a model section
}
