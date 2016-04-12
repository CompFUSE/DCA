//
// Created by giovanni on 13.04.16.
//
#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_COARSGRAINING_NAMES_HPP
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_COARSGRAINING_NAMES_HPP

namespace DCA {
enum COARSEGRAIN_DOMAIN_NAMES {
  ORIGIN,
  Q_FINE,
  K,
  K_PLUS_Q,
  Q_MINUS_K,
  TETRAHEDRON_K,
  TETRAHEDRON_ORIGIN
};

std::string to_str(COARSEGRAIN_DOMAIN_NAMES NAME) {
  switch (NAME) {
    case ORIGIN:
      return "ORIGIN";

    case Q_FINE:
      return "Q-FINE";

    case K:
      return "K";

    case K_PLUS_Q:
      return "K+Q";

    case Q_MINUS_K:
      return "Q-K";

    case TETRAHEDRON_K:
      return "TETRAHEDRON-K";

    case TETRAHEDRON_ORIGIN:
      return "TETRAHEDRON-ORIGIN";

    default:
      throw std::logic_error(__FUNCTION__);
  }

  return "???";
}
}

#endif