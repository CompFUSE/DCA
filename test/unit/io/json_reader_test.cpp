#include "gtest/gtest.h"

#include "dca/io/json/json_reader.hpp"

#include <complex>

const std::string directory = DCA_SOURCE_DIR "/test/unit/io/";

TEST(ReadTest, All) {
  dca::io::JSONReader reader;

  reader.open_file(directory + "json_reader_input.json");

  double d;
  int i;
  bool b;
  std::string s;
  std::vector<std::string> vs;
  std::vector<std::complex<double>> vc;
  std::vector<std::vector<std::string>> vvs;
  std::vector<std::vector<int>> vvi;

  reader.open_group("group 1");

  // TODO: throw
  EXPECT_NO_THROW(reader.open_group("not a group"));
  reader.close_group();

  reader.open_group("group 11");
  reader.execute("scalar", d);
  EXPECT_DOUBLE_EQ(d, 3.1456);
  reader.execute("int", i);
  EXPECT_EQ(i, 42);
  EXPECT_FALSE(reader.execute("not a field", i));
  EXPECT_FALSE(reader.execute("int", s));
  EXPECT_EQ(i, 42);
  reader.close_group();

  reader.execute("string", s);
  EXPECT_EQ(s, "ciao");

  reader.open_group("group 12");
  reader.execute("vec string", vs);
  const std::vector<std::string> vs_check{"an apple", "a pear"};
  EXPECT_EQ(vs, vs_check);
  reader.close_group();

  reader.close_group();

  reader.open_group("group 2");
  reader.close_group();

  reader.open_group("group 3");
  reader.execute("bool", b);
  EXPECT_EQ(b, false);
  reader.execute("bool 2", b);
  EXPECT_EQ(b, true);
  reader.execute("bool 3", b);
  EXPECT_EQ(b, true);
  reader.execute("bool 4", b);
  EXPECT_EQ(b, false);
  reader.close_group();

  reader.execute("vec vec string", vvs);
  const std::vector<std::vector<std::string>> vvs_check{std::vector<std::string>{"aa", "ab"},
                                                        std::vector<std::string>{"ba", "bc"}};
  EXPECT_EQ(vvs, vvs_check);
  reader.execute("vec vec int", vvi);
  const std::vector<std::vector<int>> vvi_check{std::vector<int>{11, 12}, std::vector<int>{21, 11}};
  EXPECT_EQ(vvi, vvi_check);
  reader.execute("vec cmplx", vc);
  const std::vector<std::complex<double>> vc_check{std::complex<double>(1, -1),
                                                   std::complex<double>(2, 3)};
  EXPECT_EQ(vc, vc_check);
}

TEST(ReadTest, InvalidInput) {
  auto test_file = [](const std::string& name, int err_line) {
    dca::io::JSONReader reader;
    try {
      reader.open_file(name);
      EXPECT_TRUE(false);
    }
    catch (const std::logic_error& err) {
      const std::string msg = err.what();
      const std::string check = "Line " + std::to_string(err_line) + ":";
      EXPECT_NE(std::string::npos, msg.find(check));
    }
  };

  test_file(directory + "invalid_input1.json", 8);
  test_file(directory + "invalid_input2.json", 11);
  test_file(directory + "invalid_input3.json", 6);
}
