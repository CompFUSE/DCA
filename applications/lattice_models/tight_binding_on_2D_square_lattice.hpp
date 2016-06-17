// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file defines the point group symmetry, the lattice type and the model type for a simulation
// of the tight binding model on a 2D square lattice. This file should be included by main_dca.cpp
// and main_analysis.cpp.

#ifndef DCA_APPLICATIONS_LATTICE_MODELS_TIGHT_BINDING_ON_2D_SQUARE_LATTICE_HPP
#define DCA_APPLICATIONS_LATTICE_MODELS_TIGHT_BINDING_ON_2D_SQUARE_LATTICE_HPP

// Symmetry group
#include "phys_library/domains/cluster/symmetries/point_groups/2D/2D_square.h"
using DcaPointGroupType = D4;

// Lattice type
#include "phys_library/parameters/models/analytic_hamiltonians/lattices/2D_square_lattice.h"
using LatticeType = square_lattice<DcaPointGroupType>;

// Model type
#include "phys_library/parameters/models/tight_binding_model.h"
using ModelType = tight_binding_model<LatticeType>;

#endif  // DCA_APPLICATIONS_LATTICE_MODELS_TIGHT_BINDING_ON_2D_SQUARE_LATTICE_HPP
