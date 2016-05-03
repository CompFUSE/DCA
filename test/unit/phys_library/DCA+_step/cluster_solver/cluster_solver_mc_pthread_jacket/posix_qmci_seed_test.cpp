// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific
// publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)

#include "mock_solver.hpp"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_pthread_jacket/posix_qmci_cluster_solver.h"
#include "gtest/gtest.h"

void test_posix_initialization(int n_walk, int n_acc) {
  MockRng::reset();
  MockWalker::reset(n_walk);
  MockParameters pars(n_walk, n_acc);
  MOMS_t mom;
  //create the pthread jacke and initialize the rngs
  DCA::posix_qmci_integrator<MockSolver> p_solver(pars, mom);
  p_solver.initialize(0);
  p_solver.integrate();
  std::cout << "With  nr walkers/accumulators: " << n_walk << "  " << n_acc
  << std::endl;
  //check that all the thread used for initializing rngs are different
  MockRng::checkSeeds(n_walk);
  //check if the ids in [0,n_walkers) were used
 // for(int i=0;i<n_walk;i++) std::cout<<pars.available_walker_ids[i]<<std::endl;
  for(int i=0;i<n_walk;i++) ASSERT_EQ(MockWalker::available_walker_ids[i], false);
  long int info;
  p_solver.finalize(info);
}

TEST(posix_solver, seed) {
  //this test make sure that:
  // rng::init_from_id get called with different ids and correct max_id.
  // the correct rng are passed to the walker.
  // the pthread integration runs and accumulator get sent a state before measuring.
  // TODO A test on the final number of measurments could easily be added.
  MockRng::initialize();
  pthread_mutex_init(&MockWalker::walker_lock,NULL);
  test_posix_initialization(3, 3);
  test_posix_initialization(2, 8);
  test_posix_initialization(8, 2);
  MockRng::finalize();
  pthread_mutex_destroy(&MockWalker::walker_lock);
}

int main(int argc, char **argv) {
  // The following line must be executed to initialize Google Mock
  // (and Google Test) before running the tests.
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
