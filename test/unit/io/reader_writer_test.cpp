// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides specific tests for the  reader and writer.

#include "dca/io/reader.hpp"
#include "dca/io/writer.hpp"

#include <array>
#include <cctype>
#include <complex>
#include <string>
#include <vector>

#include "gtest/gtest.h"

const std::vector<std::string> types{"JSON", "HDF5"};

std::string toLower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](char c) { return std::tolower(c); });
  return s;
}

TEST(ReaderWriterTest, ReaderDestructorCleanUp) {
  for (auto type : types) {
    std::string test_file_name = type + "_reader_test.hdf5";
    std::string group_name = "magic-numbers";
    std::string object_name = "forty-two";

    // Create test file.
    dca::io::Writer writer(type);
    const int i = 42;

    writer.open_file(test_file_name);
    writer.open_group(group_name);
    writer.execute(object_name, i);
    writer.close_group();
    writer.close_file();

    // Read test file.
    dca::io::Reader reader(type);
    int j;

    reader.open_file(test_file_name);
    reader.open_group(group_name);
    EXPECT_TRUE(reader.execute(object_name, j));

    EXPECT_EQ(i, j);

    // Destructor closes file.
    // reader.close_file();
  }
}

TEST(ReaderWriterTest, WriterDestructorCleanUp) {
  for (auto type : types) {
    std::string test_file_name = "reader_test" + toLower(type);
    std::string group_name_1 = "integers";
    std::string group_name_2 = "magic-numbers";
    std::string object_name = "forty-two";

    dca::io::Writer writer(type);
    const int i = 42;

    writer.open_file(test_file_name);
    writer.open_group(group_name_1);
    writer.open_group(group_name_2);
    writer.execute(object_name, i);

    // Destructor closes groups and file.
    // writer.close_group();
    // writer.close_group();
    // writer.close_file();
  }
}

TEST(ReaderWriterTest, VectorReadWrite) {
  for (auto type : types) {
    const std::string object_name = "a_vector";
    std::string file_name = type + "reader_test" + toLower(type);
    const std::vector<std::complex<double>> a_vector{
        std::complex<double>(1., 0.), std::complex<double>(0., 1.), std::complex<double>(23.4, -1.5)};

    // Create test file.
    dca::io::Writer writer(type);
    writer.open_file(file_name);
    writer.execute(object_name, a_vector);
    writer.close_file();

    // Read test file.
    dca::io::Reader reader(type);
    std::vector<std::complex<double>> vector_read;
    reader.open_file(file_name);
    EXPECT_TRUE(reader.execute(object_name, vector_read));

    ASSERT_EQ(a_vector.size(), vector_read.size());
    for (int i = 0; i < a_vector.size(); ++i) {
      EXPECT_DOUBLE_EQ(std::real(a_vector[i]), std::real(vector_read[i]));
      EXPECT_DOUBLE_EQ(std::imag(a_vector[i]), std::imag(vector_read[i]));
    }

    reader.close_file();
  }
}

TEST(ReaderWriterTest, VectorOfVectorsReadWrite) {
  for (auto type : types) {
    const std::string object_name = "a_vector";
    const std::string file_name = "test_vector_vector." + toLower(type);
    const std::vector<std::vector<double>> data_unequal_size{{0, 0, 2}, {1}, {1, 0}, {}, {0, 0}};

    // Create test file.
    dca::io::Writer writer(type);
    writer.open_file(file_name);
    writer.execute(object_name, data_unequal_size);
    writer.close_file();

    // Read test file.
    dca::io::Reader reader(type);
    std::vector<std::vector<double>> data_read;
    reader.open_file(file_name);
    EXPECT_TRUE(reader.execute(object_name, data_read));

    EXPECT_EQ(data_unequal_size, data_read);
    reader.close_file();
  }
}

TEST(ReaderWriterTest, VectorOfArraysReadWrite) {
  for (auto type : types) {
    const std::string object_name = "obj";
    const std::string file_name = "test_vec_of_arr." + toLower(type);

    std::vector<std::array<int, 3>> data{{-1, 2, 3}, {5, -7, 0}};

    // Create test file.
    dca::io::Writer writer(type);
    writer.open_file(file_name);
    writer.execute(object_name, data);
    writer.close_file();

    // Read test file.
    dca::io::Reader reader(type);
    std::vector<std::array<int, 3>> data_read;
    reader.open_file(file_name);
    EXPECT_TRUE(reader.execute(object_name, data_read));

    EXPECT_EQ(data, data_read);

    reader.close_file();
  }
}

TEST(ReaderWriterTest, StringAndVectorOfStringsReadWrite) {
  for (auto type : types) {
    std::vector<std::string> s_vec1{"foo", "", "baz"};
    std::string s1{"bazinga"};
    const std::string filename = "test_vec_of_strings" + toLower(type);

    // Create test file.
    dca::io::Writer writer(type);
    writer.open_file(filename);
    writer.execute("single-string", s1);
    writer.execute("strings", s_vec1);
    writer.close_file();

    // Read test file.
    dca::io::Reader reader(type);
    reader.open_file(filename);
    //
    std::vector<std::string> s_vec2;
    std::string s2;
    EXPECT_TRUE(reader.execute("strings", s_vec2));
    EXPECT_TRUE(reader.execute("single-string", s2));
    EXPECT_EQ(s_vec1, s_vec2);
    EXPECT_EQ(s1, s2);
  }
}

template <typename Scalar>
class ReaderWriterTest : public ::testing::Test {};
using TestTypes = ::testing::Types<float, std::complex<double>>;
TYPED_TEST_CASE(ReaderWriterTest, TestTypes);

TYPED_TEST(ReaderWriterTest, FunctionReadWrite) {
  for (auto type : types) {
    using Dmn1 = dca::func::dmn_0<dca::func::dmn<5>>;
    using Dmn2 = dca::func::dmn_0<dca::func::dmn<4>>;
    using Dmn3 = dca::func::dmn_0<dca::func::dmn<2>>;
    using Dmn = dca::func::dmn_variadic<Dmn1, Dmn2, Dmn3>;
    using Scalar = TypeParam;

    dca::func::function<Scalar, Dmn> f1("myfunc");
    int val = 0;
    for (auto& x : f1)
      x = ++val;

    dca::io::Writer writer(type);
    writer.open_file("test_func" + toLower(type), true);

    writer.execute(f1);

    writer.close_file();

    // Read test file.
    dca::io::Reader reader(type);
    reader.open_file("test_func" + toLower(type));

    dca::func::function<Scalar, Dmn> f2("myfunc");

    EXPECT_TRUE(reader.execute(f2));

    for (int i = 0; i < f1.size(); ++i)
      EXPECT_EQ(f1(i), f2(i));

    reader.close_file();
  }
}

TYPED_TEST(ReaderWriterTest, MatrixReadWrite) {
  for (auto type : types) {
    using Scalar = TypeParam;
    const std::pair<int, int> m_size(7, 3);
    const std::string filename = "test_mat." + toLower(type);

    dca::linalg::Matrix<Scalar, dca::linalg::CPU> m1(m_size, "mymat");
    int val = 0;
    for (int j = 0; j < m1.nrCols(); ++j)
      for (int i = 0; i < m1.nrRows(); ++i)
        m1(i, j) = ++val;

    dca::io::Writer writer(type);
    writer.open_file(filename, true);
    writer.execute(m1);
    writer.close_file();

    // Read test file.
    dca::io::Reader reader(type);
    reader.open_file(filename);

    dca::linalg::Matrix<Scalar, dca::linalg::CPU> m2;
    EXPECT_TRUE(reader.execute("mymat", m2));
    reader.close_file();

    EXPECT_EQ(m1.get_name(), m2.get_name());
    EXPECT_EQ(m1, m2);
  }
}

TEST(ReaderWriterTest, NonAccessibleFile) {
  for (auto type : types) {
    dca::io::Writer writer(type);
    H5::Exception::dontPrint();
    EXPECT_ANY_THROW(writer.open_file("not_existing_directory/file.txt"));

    dca::io::Reader reader(type);
    EXPECT_THROW(reader.open_file("not_existing_file.txt"), std::runtime_error);
  }
}

TEST(ReaderWriterTest, FunctionNotPresent) {
  for (auto type : types) {
    using Dmn = dca::func::dmn_0<dca::func::dmn<5, int>>;

    dca::func::function<int, Dmn> present("present");
    present = 1;

    dca::io::Writer writer(type);
    writer.open_file("hdf5_missing_func" + toLower(type));
    writer.execute(present);
    writer.close_file();

    dca::func::function<int, Dmn> not_present("not_present");
    present = 0;

    dca::io::Reader reader(type);
    reader.open_file("hdf5_missing_func" + toLower(type));
    EXPECT_FALSE(reader.execute(not_present));
    EXPECT_TRUE(reader.execute(present));

    for (int val : not_present)
      EXPECT_EQ(0, val);
    for (int val : present)
      EXPECT_EQ(1, val);
  }
}

TEST(ReaderWriterTest, GroupOpenclose) {
  for (auto type : types) {
    dca::io::Writer writer(type);
    writer.open_file("group_open_close" + toLower(type));

    writer.open_group("foo");
    writer.execute("a", 0);
    writer.close_group();
    writer.open_group("foo");
    writer.execute("b", 1);
    writer.close_group();
    writer.open_group("bar");
    writer.execute("b2", 1.5);
    writer.close_group();
    writer.open_group("foo");
    writer.execute("c", 2);

    writer.close_file();

    dca::io::Reader reader(type);
    reader.open_file("group_open_close" + toLower(type));

    int i_val;
    double d_val;

    reader.open_group("foo");
    reader.execute("a", i_val);
    EXPECT_EQ(0, i_val);
    reader.execute("b", i_val);
    EXPECT_EQ(1, i_val);
    reader.execute("c", i_val);
    EXPECT_EQ(2, i_val);
    reader.close_group();
    reader.open_group("bar");
    reader.execute("b2", d_val);
    EXPECT_EQ(1.5, d_val);
  }
}

TEST(ReaderWriterTest, Overwrite) {
  for (auto type : types) {
    dca::io::Writer writer(type);
    writer.open_file("test" + toLower(type), true);

    writer.open_group("foo");
    writer.execute("a", 1);
    writer.execute("a", 2);

    // Try to write with different size.
    // EXPECT_THROW(writer.execute("a", std::pair<int, int>(1, 1)), std::length_error);

    writer.close_file();

    dca::io::Reader reader(type);
    reader.open_file("test" + toLower(type));

    int i_val;
    reader.open_group("foo");
    reader.execute("a", i_val);
    EXPECT_EQ(2, i_val);
  }
}
