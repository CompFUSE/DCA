#ifndef SUBMATIMPL_HPP
#define SUBMATIMPL_HPP

#include <array>
#include <vector>
#include "dca/linalg/linalg.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {

template<typename SCALAR>
class SubmatImpl {
public:
  struct DelayedMoveType {
    Move move_type;
    std::array<double, 3> removal_rng{1., 1., 1.};
    Real acceptance_rng;
    std::array<int, 2> indices{-1, -1};
  };

private:
  const DelayedMoveType* current_move_;
  std::array<int, 2> Gamma_size_;
  std::array<std::vector<int>, 2> Gamma_indices_;
  std::array<std::vector<int>, 2> sector_indices_;

  std::array<unsigned, 2> nbr_of_indices_;

  std::vector<int> index_;

  std::vector<DelayedMoveType> delayed_moves_;
  std::array<linalg::util::HostVector<int>, 2> move_indices;
  std::array<linalg::util::HostVector<int>, 2> removal_list;
  std::array<linalg::util::HostVector<int>, 2> source_list;
  std::array<std::vector<int>, 2> insertion_list;
  std::array<std::vector<int>, 2> insertion_Gamma_indices;
  std::vector<int> conf_removal_list;
  std::vector<int>::iterator insertion_list_it;

public:
  void mainSubmatrixProcess();
};

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif
