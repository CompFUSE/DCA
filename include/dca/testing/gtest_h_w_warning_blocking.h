// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter W. Doak  doakpw@ornl.gov  Oak Ridge National Lab
//
// This file implements a work around for at least older version of google test
// tripping newer compiler warnings. Eventually we will upgrade our bundled version
// of google test or migrate to a different test framework.

#ifndef DCA_GTEST_H_W_WARNING_BLOCKING_H
#define DCA_GTEST_H_W_WARNING_BLOCKING_H

#if defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-warning-option"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#else
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsuggest-override"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

#include "gtest/gtest.h"

#pragma GCC diagnostic_pop

#endif
