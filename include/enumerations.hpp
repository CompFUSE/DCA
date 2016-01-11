// DCA++ enumerations.

namespace DCA
{
  enum CLUSTER_SOLVER_NAMES {HIGH_TEMPERATURE_SERIES,
                             ED_CLUSTER_SOLVER,
                             ADVANCED_ED_CLUSTER_SOLVER,
                             CT_AUX_CLUSTER_SOLVER,
                             SS_CT_HYB};
}

enum cluster_shape {BETT_CLUSTER,
                    PARALLELEPIPED};
using cluster_shape_type = cluster_shape;

// enum representation {IRREDUCIBLE, FULL, VERTEX};
// typedef representation representation_type;

// enum    legendre_representation {SINGLE_PARTICLE_QUANTITY, TWO_PARTICLE_QUANTITY};
// typedef legendre_representation legendre_representation_type;

// enum vertex_frequency_representation {COMPACT, EXTENDED, COMPACT_POSITIVE, EXTENDED_POSITIVE, COMPACT_SORTED, EXTENDED_SORTED,
//              EXTENDED_BOSONIC, EXTENDED_FERMIONIC, CORE_SORTED, HIGH_FREQUENCY_SORTED};
// typedef vertex_frequency_representation vertex_frequency_representation_type;

enum vertex_measurement {NONE,
                         PARTICLE_HOLE_TRANSVERSE,
                         PARTICLE_HOLE_MAGNETIC,
                         PARTICLE_HOLE_CHARGE,
                         PARTICLE_PARTICLE_SUPERCONDUCTING};
using vertex_measurement_type = vertex_measurement;

enum mesh_shape {PARALLELLOGRAM,
                 FIRST_BRILLOUIN_ZONE};
// PARALLELLOGRAM --> might break the symmetry, be carefull !!!!
const static mesh_shape MESH_SHAPE = FIRST_BRILLOUIN_ZONE;//PARALLELLOGRAM;

// enum    e_spin_states {e_DN=-1, e_UP=1};
// typedef e_spin_states e_spin_states_type;

// enum    HS_spin_states {HS_DN=-1, HS_ZERO=0, HS_UP=1};
// typedef HS_spin_states HS_spin_states_type;

// enum    HS_field_sign {HS_FIELD_DN=-1, HS_FIELD_UP=1};
// typedef HS_field_sign HS_field_sign_type;

// enum    HS_vertex_move {ANNIHILATION=-1, STATIC=0, CREATION=1};
// typedef HS_vertex_move HS_vertex_move_type;

// enum    stop_watch {DEFAULT=-1, START=0, STOP=1, SET=2, CONTINUE=3};
// typedef stop_watch stop_watch_type;

// enum eigenvalue_degeneracy {NO_DEGENERACY, TWOFOLD_DEGENERACY, THREEFOLD_DEGENERACY, FOURFOLD_DEGENERACY};
// typedef eigenvalue_degeneracy eigenvalue_degeneracy_t;

enum integration_method {DELTA_FUNCTION_INTEGRATION,
                         ISOLATED_CLUSTER_INTEGRATION,
                         TRAPEZIUM_INTEGRATION,
                         QUADRATURE_INTEGRATION,
                         TETRAHEDRON_INTEGRATION};
using integration_method_type = integration_method;

enum minimization_method {GRADIENT_METHOD,
                          WEIGHTED_GRADIENT_METHOD,
                          CONJUGATE_GRADIENT_METHOD,
                          LEVMAR_LIBRARY};
using minimization_method_type = minimization_method;

// enum MPI_library {SERIAL_LIBRARY, MPI_LIBRARY, MPI_FT_LIBRARY, OPENMPI_FT_LIBRARY};
// typedef MPI_library MPI_library_type;

enum MC_integration_method {CT_AUX,
                            HYBRIDIZATION,
                            HYBRIDIZATION_FULL,
                            PCM,
                            ANALYSIS,
                            ANALYSIS_INTERPOLATION,
                            ANALYSIS_COMPUTE_REDUCIBLE_VERTEX,
                            HIGH_TEMPERATURE_SERIES_SOLVER,
                            ED_CLUSTER_SOLVER};
using MC_integration_method_type = MC_integration_method;

enum BRILLOUIN_ZONE {BRILLOUIN_ZONE_CUT_TEMPLATE,
                     FERMI_SURFACE_SQUARE_2D_LATTICE,
                     SQUARE_2D_LATTICE,
                     BODY_CENTERED_TETRAGONAL_A,
                     BODY_CENTERED_TETRAGONAL_B,
                     SIMPLE_TETRAGONAL,
                     TRICLINIC,
                     FACE_CENTERED_CUBIC,
                     BODY_CENTERED_CUBIC,
                     SIMPLE_CUBIC,
                     HEXAGONAL,
                     RHOMBOHEDRAL_A,
                     RHOMBOHEDRAL_B,
                     SIMPLE_MONOCLINIC,
                     ONE_FACE_CENTERED_MONOCLINIC_A,
                     ONE_FACE_CENTERED_MONOCLINIC_B,
                     SIMPLE_ORTHOROMBIC,
                     BASE_CENTERED_ORTHORHOMBIC,
                     BODY_CENTERED_ORTHOROMBIC,
                     ALL_FACE_CENTERED_ORTHORHOMBIC_A,
                     ALL_FACE_CENTERED_ORTHORHOMBIC_B};
using BRILLOUIN_ZONE_CUT_TYPE = BRILLOUIN_ZONE;

// enum c_star_algebra_type {ANNIHILATION_OPERATOR, IDENTITY_OPERATOR, CREATION_OPERATOR, NONE_OPERATOR};
// typedef c_star_algebra_type c_star_algebra_t;

// enum MC_accumulator_method {MATSUBARA, LEGENDRE, SERIES};
// typedef MC_accumulator_method MC_accumulator_method_type;

// enum symmetry_operation {PARTICLE_NUMBER, Sz, TOTAL_MOMENTUM};
// typedef symmetry_operation symmetry_operation_type;
