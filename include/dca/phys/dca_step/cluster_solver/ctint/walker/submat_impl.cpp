#include "submat_impl.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {

template <typename SCALAR>
void SubmatImpl<SCALAR>::mainSubmatrixProcess() {
  for (int s = 0; s < 2; ++s) {
    move_indices_[s].clear();
    insertion_list_[s].clear();
    insertion_Gamma_indices_[s].clear();
  }

  std::vector<int> aux_spin_type, new_aux_spin_type, move_band;

  SCALAR acceptance_prob_ = 1.0;
  for (int delay_ind = 0; delay_ind < delayed_moves_.size(); ++delay_ind) {
    current_move_ = &delayed_moves_[delay_ind];
    const auto move_type = current_move_->move_type;
    // there seems to be no need that this be persistent.
    det_ratio_ = 1.;
    std::array<int, 2> indices_array;

    if (move_type == INSERTION) {
      indices_array = current_move_->indices;
    }
    else {  // move_type == REMOVAL
      indices_array = configuration_.randomRemovalCandidate(delayed_moves_[delay_ind].removal_rng);
    }

    index_.clear();
    for (auto idx : indices_array) {
      if (idx >= 0)
        index_.push_back(idx);
    }

    // This leads to a bunch of complex branchy and premature looking optimization
    // The branch predictor must be weeping
    bool at_least_one_recently_added = false;
    bool all_recently_added = false;
    if (move_type == REMOVAL) {
      // Check if the vertex to remove was inserted during the current submatrix update.
      const auto recently_added = [=](int id) { return id >= config_size_init_; };
      all_recently_added = math::util::all(index_, recently_added);
      at_least_one_recently_added = math::util::any(index_, recently_added);
    }

    for (int s = 0; s < 2; ++s) {
      Gamma_indices_[s].clear();
      new_aux_spin_type.clear();
      aux_spin_type.clear();
      move_band.clear();

      findSectorIndices(s);

      if (move_type == INSERTION) {
        for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
          new_aux_spin_type.push_back(
              configuration_.getSector(s).getAuxFieldType(sector_indices_[s][ind]));
          aux_spin_type.push_back(0);
        }
      }
      else if (move_type == REMOVAL) {
        for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
          new_aux_spin_type.push_back(0);
          aux_spin_type.push_back(
              configuration_.getSector(s).getAuxFieldType(sector_indices_[s][ind]));

          if (recentlyAdded(ind, s)) {
            // Find pre-existing Gamma_indices_.
            insertion_list_it_ = std::find(insertion_list_[s].begin(), insertion_list_[s].end(),
                                           sector_indices_[s][ind]);
            Gamma_indices_[s].push_back(insertion_Gamma_indices_[s][std::distance(
                insertion_list_[s].begin(), insertion_list_it_)]);
          }
        }
      }  // endif removal
      for (int index : sector_indices_[s])
        move_band.push_back(configuration_.getSector(s).getLeftB(index));

      for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
        if (move_type == INSERTION || !recentlyAdded(ind, s))
          gamma_[s].push_back(
              gamma_values_[std::make_pair(aux_spin_type[ind], new_aux_spin_type[ind])][move_band[ind]]);
        else {
          // TODO: find instead of adding.
          gamma_[s].push_back(
              gamma_values_[std::make_pair(new_aux_spin_type[ind], aux_spin_type[ind])][move_band[ind]]);
        }
      }

      if (!at_least_one_recently_added)
        computeInsertionMatrices(sector_indices_[s], s);
      else {
        if (all_recently_added)
          computeRemovalMatrix(s);
        else
          computeMixedInsertionAndRemoval(s);
      }

      det_ratio_ *= details::smallDeterminant(s_[s]);
    }  // s loop.

    if (n_ == 0 && move_type == REMOVAL)
      continue;

    // Compute acceptance probability.
    auto [acceptance_prob, mc_weight_ratio] = computeAcceptanceProbability();
    acceptance_prob_ = acceptance_prob;

    // Note: acceptance and rejection can be forced for testing with the appropriate "acceptance_rng".
    const bool accepted =
        delayed_moves_[delay_ind].acceptance_rng < std::min(std::abs(acceptance_prob_), Real(1.));

    // NB: recomputeGammaInv is just a inefficient alternative to updateGammaInv. Only for testing
    // or debbuging.
    // recomputeGammaInv();

    // Update
    if (accepted) {
      ++BaseClass::n_accepted_;

      BaseClass::mc_log_weight_ += std::log(std::abs(mc_weight_ratio));

      // Are we capturing the avg sign properly wrt multiple delayed moves
      phase_.multiply(acceptance_prob_);

      // Update GammaInv if necessary.
      if (!at_least_one_recently_added)
        for (int s = 0; s < 2; ++s)
          updateGammaInv(s);
      else {
        removeRowAndColOfGammaInv();
        for (int s = 0; s < 2; ++s) {
          for (int ind = sector_indices_[s].size() - 1; ind >= 0; --ind) {
            if (!recentlyAdded(ind, s))
              continue;
            insertion_list_[s].erase(std::remove(insertion_list_[s].begin(),
                                                 insertion_list_[s].end(), sector_indices_[s][ind]),
                                     insertion_list_[s].end());
          }
          for (int ind = Gamma_indices_[s].size() - 1; ind >= 0; --ind) {
            insertion_Gamma_indices_[s].erase(std::remove(insertion_Gamma_indices_[s].begin(),
                                                          insertion_Gamma_indices_[s].end(),
                                                          Gamma_indices_[s][ind]),
                                              insertion_Gamma_indices_[s].end());
            gamma_[s].erase(gamma_[s].begin() + Gamma_indices_[s][ind]);

            for (int i = 0; i < insertion_Gamma_indices_[s].size(); ++i) {
              if (insertion_Gamma_indices_[s][i] > Gamma_indices_[s][ind])
                --insertion_Gamma_indices_[s][i];
            }
          }
        }
      }

      for (int s = 0; s < 2; ++s) {
        for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
          if (move_type == INSERTION || !recentlyAdded(ind, s)) {
            move_indices_[s].push_back(sector_indices_[s][ind]);
          }
          else {
            gamma_[s].pop_back();
            move_indices_[s].erase(std::remove(move_indices_[s].begin(), move_indices_[s].end(),
                                               sector_indices_[s][ind]),
                                   move_indices_[s].end());
          }
        }
      }

      if (move_type == INSERTION) {
        n_ += index_.size();
        for (auto idx : index_)
          configuration_.commitInsertion(idx);

        for (int s = 0; s < 2; ++s) {
          for (int ind = 0; ind < nbr_of_indices_[s]; ++ind) {
            insertion_list_[s].push_back(sector_indices_[s][ind]);
            insertion_Gamma_indices_[s].push_back(Gamma_inv_[s].nrRows() - nbr_of_indices_[s] + ind);
          }
        }

        // TODO: cleanup
        for (int idx : index_)
          conf_removal_list_.erase(
              std::remove(conf_removal_list_.begin(), conf_removal_list_.end(), idx),
              conf_removal_list_.end());
      }
      else if (move_type == REMOVAL) {
        n_ -= index_.size();
        for (auto idx : index_)
          configuration_.markForRemoval(idx);

        // TODO: cleanup.
        for (int idx : index_)
          conf_removal_list_.push_back(idx);
      }
    }
    else {  // The move is rejected:
      for (int s = 0; s < 2; ++s) {
        for (int ind = 0; ind < nbr_of_indices_[s]; ++ind)
          gamma_[s].pop_back();
      }
      if (at_least_one_recently_added && !all_recently_added)
        for (int s = 0; s < 2; ++s) {
          Gamma_inv_[s].swap(Gamma_inv_cpy_[s]);
        }
    }
  }

  // This seems correct in that we count all steps
  // just as is done for the non submatrix cpu version
  BaseClass::n_steps_ += delayed_moves_.size();
}

}
}
}
}
