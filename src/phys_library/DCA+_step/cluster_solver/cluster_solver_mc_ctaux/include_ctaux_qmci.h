//-*-C++-*-

#include "ctaux_typedefinitions.h"

// CT-AUX-domains
#include "ctaux_domains/Feynman_expansion_order_domain.h"
#include "ctaux_domains/HS_spin_domain.h"
#include "ctaux_domains/HS_field_sign_domain.h"
#include "ctaux_domains/HS_vertex_move_domain.h"

// CT-AUX STRUCTURES
#include "ctaux_structs/ctaux_auxilery_field_coefficients.h"

#include "ctaux_structs/ctaux_vertex_singleton.h"
#include "ctaux_structs/ctaux_vertex_pair.h"
#include "ctaux_structs/ctaux_hubbard_stratonovitch_configuration.h"

#include "ctaux_structs/ctaux_walker_data.h"

// CT-AUX WALKER

#include "ctaux_walker/ctaux_walker_tools/ctaux_shrink_routines/ctaux_shrink_routines_TEM.h"
#include "ctaux_walker/ctaux_walker_tools/ctaux_shrink_routines/ctaux_shrink_routines_CPU.h"
#include "ctaux_walker/ctaux_walker_tools/ctaux_shrink_routines/ctaux_shrink_routines_GPU.h"

#include "ctaux_walker/ctaux_walker_tools/ctaux_shrink_routines.h"

#include "ctaux_walker/ctaux_walker_tools/ctaux_N_matrix_routines/ctaux_N_matrix_routines_TEM.h"
#include "ctaux_walker/ctaux_walker_tools/ctaux_N_matrix_routines/ctaux_N_matrix_routines_CPU.h"
#include "ctaux_walker/ctaux_walker_tools/ctaux_N_matrix_routines/ctaux_N_matrix_routines_GPU.h"

#include "ctaux_walker/ctaux_walker_tools/ctaux_N_matrix_routines.h"

#include "ctaux_walker/ctaux_walker_tools/ctaux_G_matrix_routines/ctaux_G_matrix_routines_TEM.h"
#include "ctaux_walker/ctaux_walker_tools/ctaux_G_matrix_routines/ctaux_G_matrix_routines_CPU.h"
#include "ctaux_walker/ctaux_walker_tools/ctaux_G_matrix_routines/ctaux_G_matrix_routines_GPU.h"

#include "ctaux_walker/ctaux_walker_tools/ctaux_G_matrix_routines.h"

#include "ctaux_walker/ctaux_walker_tools/ctaux_G0_matrix_routines/ctaux_G0_matrix_routines.h"

#include "ctaux_walker/ctaux_walker_tools/ctaux_G0_matrix_routines/ctaux_G0_matrix_routines_TEM.h"
#include "ctaux_walker/ctaux_walker_tools/ctaux_G0_matrix_routines/ctaux_G0_matrix_routines_CPU.h"
#include "ctaux_walker/ctaux_walker_tools/ctaux_G0_matrix_routines/ctaux_G0_matrix_routines_GPU.h"

#include "ctaux_walker/ctaux_walker_routines_template.h"
#include "ctaux_walker/ctaux_walker_routines_CPU.h"
#include "ctaux_walker/ctaux_walker_routines_GPU.h"

#include "ctaux_walker/ctaux_walker_build_in_test.h"

#include "ctaux_walker.h"

// CT-AUX ACCUMULATOR

// sp CT-AUX ACCUMULATOR
#include "ctaux_accumulator/sp_accumulator/ctaux_sp_accumulator_nfft.h"

// tp CT-AUX ACCUMULATOR
#include "ctaux_accumulator/tp_accumulator/ctaux_tp_nft.h"

#include "ctaux_accumulator/tp_accumulator/ctaux_accumulator_nonlocal_G.h"
#include "ctaux_accumulator/tp_accumulator/ctaux_accumulator_nonlocal_chi_atomic.h"
#include "ctaux_accumulator/tp_accumulator/ctaux_accumulator_nonlocal_chi.h"

#include "ctaux_accumulator/tp_accumulator/ctaux_accumulator_equal_time_operator.h"

#include "ctaux_accumulator.h"

// CT-AUX SOLVER

#include "ctaux_cluster_solver.h"
