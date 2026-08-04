// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file provides unit tests or the json reader and writer.

#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "dca/io/json/json_reader.hpp"

#include <complex>

const std::string directory = DCA_SOURCE_DIR "/test/unit/io/";

TEST(ReadTest, All) {
  dca::io::JSONReader reader;

  reader.open_file(directory + "json_reader_input.json");

  double d;
  int i;
  bool b;
  std::string s;
  std::vector<std::string> vs;
  std::vector<std::complex<double>> vc;
  std::vector<std::vector<std::string>> vvs;
  std::vector<std::vector<int>> vvi;

  reader.open_group("group 1");

  // TODO: throw
  bool group_opened = reader.open_group("not a group");
  EXPECT_EQ(group_opened, false);

  group_opened = reader.open_group("group 11");
  EXPECT_EQ(group_opened, true);
  group_opened = reader.execute("scalar", d);
  EXPECT_EQ(group_opened, true);
  EXPECT_DOUBLE_EQ(d, 3.1456);
  reader.execute("int", i);
  EXPECT_EQ(i, 42);
  EXPECT_FALSE(reader.execute("not a field", i));
  EXPECT_FALSE(reader.execute("int", s));
  EXPECT_EQ(i, 42);
  reader.close_group();

  reader.execute("string", s);
  EXPECT_EQ(s, "ciao");

  reader.open_group("group 12");
  reader.execute("vec string", vs);
  const std::vector<std::string> vs_check{"an apple", "a pear"};
  EXPECT_EQ(vs, vs_check);
  reader.close_group();

  reader.close_group();

  reader.open_group("group 2");
  reader.close_group();

  reader.open_group("group 3");
  reader.execute("bool", b);
  EXPECT_EQ(b, false);
  reader.execute("bool 2", b);
  EXPECT_EQ(b, true);
  reader.execute("bool 3", b);
  EXPECT_EQ(b, true);
  reader.execute("bool 4", b);
  EXPECT_EQ(b, false);
  reader.close_group();

  std::vector<int> vi;
  reader.execute("vec int", vi);
  const std::vector<int> vi_check{1, 2, 3};
  EXPECT_EQ(vi, vi_check);
  
  reader.execute("vec vec string", vvs);
  const std::vector<std::vector<std::string>> vvs_check{std::vector<std::string>{"aa", "ab"},
                                                        std::vector<std::string>{"ba", "bc"}};
  EXPECT_EQ(vvs, vvs_check);
  reader.execute("vec vec int", vvi);
  const std::vector<std::vector<int>> vvi_check{std::vector<int>{11, 12}, std::vector<int>{21, 11}};
  EXPECT_EQ(vvi, vvi_check);
  reader.execute("vec cmplx", vc);
  const std::vector<std::complex<double>> vc_check{std::complex<double>(1, -1),
                                                   std::complex<double>(2, 3)};
  EXPECT_EQ(vc, vc_check);
}

// Issue #300: the reader must report which top-level sections were read, so callers can detect
// input that the executable silently ignores.
TEST(ReadTest, TopLevelGroupAccessTracking) {
  dca::io::JSONReader reader;
  reader.open_file(directory + "model_sections_input.json");

  // Before any access, every top-level group is reported as unread.
  auto groups = reader.topLevelGroupAccess();
  // Sorted by name: bilayer-Hubbard-model, domains, single-band-Hubbard-model.
  ASSERT_EQ(groups.size(), 3u);
  EXPECT_EQ(groups[0].name, "bilayer-Hubbard-model");
  EXPECT_EQ(groups[1].name, "domains");
  EXPECT_EQ(groups[2].name, "single-band-Hubbard-model");
  for (const auto& group : groups)
    EXPECT_FALSE(group.accessed) << group.name << " should not be accessed yet";

  // Open one of them, as the executable's ModelParameters would.
  EXPECT_TRUE(reader.open_group("single-band-Hubbard-model"));
  reader.close_group();

  // A failed open (e.g. a typo'd section name) must not mark anything accessed.
  EXPECT_FALSE(reader.open_group("single-band-Hubard-model"));

  groups = reader.topLevelGroupAccess();
  ASSERT_EQ(groups.size(), 3u);
  EXPECT_FALSE(groups[0].accessed);  // bilayer-Hubbard-model: present but never opened
  EXPECT_FALSE(groups[1].accessed);  // domains: never opened
  EXPECT_TRUE(groups[2].accessed);   // single-band-Hubbard-model: opened
}

TEST(ReadTest, InvalidInput) {
  auto test_file = [](const std::string& name, int err_line) {
    dca::io::JSONReader reader;
    try {
      reader.open_file(name);
      EXPECT_TRUE(false);
    }
    catch (const std::logic_error& err) {
      const std::string msg = err.what();
      const std::string check = "Line " + std::to_string(err_line) + ":";
      EXPECT_NE(std::string::npos, msg.find(check));
    }
  };

  test_file(directory + "invalid_input1.json", 8);
  test_file(directory + "invalid_input2.json", 11);
  test_file(directory + "invalid_input3.json", 6);
}
