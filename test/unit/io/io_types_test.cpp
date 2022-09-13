// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file tests io_types

#include "dca/io/io_types.hpp"

#include <string>
#include "gtest/gtest.h"

TEST(IoTypesTest, IoTypeFromExtension) {
  std::string hdf5_file{"dca.hdf5"};
  // This is an unfortunate legacy behavior
  std::string hdf5_autoresume_file{"dca.hdf5.tmp"};
  std::string adios2_file{"dca.bp"};
  std::string json_file{"input.json"};

  EXPECT_EQ(dca::io::IOType::HDF5, dca::io::extensionToIOType(hdf5_file));
}
