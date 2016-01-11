//-*-C++-*-

// template
#include "cluster_solver_template.h"

// HTS
#include "cluster_solver_series_expansion/include_high_temperature_series_expansion_solver.h"

// ED
//#include "include_ed.h"
#include "cluster_solver_exact_diagonalization_advanced/include_advanced_ed_solver.h"

// QMCI-template
#include "cluster_solver_mc_template/include_qmci_templates.h"
#include "cluster_solver_mc_pthread_jacket/include_posix_qmci_integration.h"
#include "cluster_solver_mc_ctaux/include_ctaux_qmci.h"
#include "cluster_solver_ss_hybridization/include_ss_hybridization_qmci.h"
