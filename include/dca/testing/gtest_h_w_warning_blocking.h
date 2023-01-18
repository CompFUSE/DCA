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
// giving little effort to removing warning trips

#ifndef DCA_GTEST_H_W_WARNING_BLOCKING_
#define DCA_GTEST_H_W_WARNING_BLOCKING_H

#pragma GCC diagnostic push
#include "gtest/gtest.h"
#pragma GCC diagnostic ignored "-Wsuggest-override"
#pragma GCC diagnostic_pop

#endif
