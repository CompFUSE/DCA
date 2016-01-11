//-*-C++-*-

std::string to_str(int DIMENSION)
{
  std::stringstream ss;
  ss << DIMENSION;
  return ss.str();
}

enum CLUSTER_NAMES {CLUSTER, LATTICE_SP, LATTICE_TP, TMP_CLUSTER, VASP_LATTICE};

std::string to_str(CLUSTER_NAMES NAME)
{
  switch(NAME)
    {
    case CLUSTER:
      return "CLUSTER";

    case LATTICE_SP:
      return "LATTICE_SP";

    case LATTICE_TP:
      return "LATTICE_TP";

    case TMP_CLUSTER:
      return "TMP_CLUSTER";

    case VASP_LATTICE:
      return "VASP_LATTICE";

    default:
      throw std::logic_error(__FUNCTION__);
    }
}

enum CLUSTER_REPRESENTATION {REAL_SPACE, MOMENTUM_SPACE};

std::string to_str(CLUSTER_REPRESENTATION NAME)
{
  switch(NAME)
    {
    case REAL_SPACE:
      return "REAL_SPACE";

    case MOMENTUM_SPACE:
      return "MOMENTUM_SPACE";

    default:
      throw std::logic_error(__FUNCTION__);
    }
}

enum CLUSTER_SHAPE {BRILLOUIN_ZONE, PARALLELLEPIPEDUM};

std::string to_str(CLUSTER_SHAPE NAME)
{
  switch(NAME)
    {
    case BRILLOUIN_ZONE:
      return "BRILLOUIN_ZONE";

    case PARALLELLEPIPEDUM:
      return "PARALLELLEPIPEDUM";

    default:
      throw std::logic_error(__FUNCTION__);
    }
}

#include "cluster_operations.h"

#include "cluster_domain.h"
#include "cluster_domain_family.h"

#include "centered_cluster_domain.h"

#include "cluster_domain_initializer.h"

#include "symmetrization_algorithms/set_symmetry_matrices.h"
#include "symmetrization_algorithms/search_maximal_symmetry_group.h"

#include "symmetrization_algorithms/cluster_reduction.h"

#include "cluster_domain_symmetry.h"
#include "cluster_domain_symmetry_initializer.h"

#include "cluster_domain_iterator.h"

#include "interpolation/include_cluster_interpolation_routines.h"
