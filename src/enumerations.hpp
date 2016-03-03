// File: enumerations.hpp

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

enum vertex_measurement {NONE,
                         PARTICLE_HOLE_TRANSVERSE,
                         PARTICLE_HOLE_MAGNETIC,
                         PARTICLE_HOLE_CHARGE,
                         PARTICLE_PARTICLE_SUPERCONDUCTING};
using vertex_measurement_type = vertex_measurement;

enum mesh_shape {PARALLELLOGRAM,
                 FIRST_BRILLOUIN_ZONE};
const static mesh_shape MESH_SHAPE = FIRST_BRILLOUIN_ZONE;//PARALLELLOGRAM;


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
