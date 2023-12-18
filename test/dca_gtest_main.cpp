

#include "dca/testing/gtest_h_w_warning_blocking.h"

#ifdef DCA_HAVE_HPX
#include <hpx/hpx_start.hpp>

int hpx_main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  auto ret_gtest = RUN_ALL_TESTS();
  hpx::finalize();
  return ret_gtest;
}

int main(int argc, char** argv) {
  std::cout << "doing HPX gtest main\n";
  hpx::start(argc, argv);
  // returns the value that hpx::main returns
  return hpx::stop();
}

#else
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
#endif
