// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Conditional inclusion of std::filesystem

#ifndef DCA_IO_FILESYSTEM
#define DCA_IO_FILESYSTEM

#if __has_include(<filesystem>)

#include <filesystem>
namespace filesystem = std::filesystem;

#else

#include <experimental/filesystem>
namespace filesystem = std::experimental::filesystem;

#endif

#endif  // DCA_IO_FILESYSTEM
