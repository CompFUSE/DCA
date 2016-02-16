#include "gtest/gtest.h"

TEST(DummyTest1, someTest) {
  int x = 42;
  EXPECT_EQ(42, x);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
