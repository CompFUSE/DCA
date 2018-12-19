// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides specific tests for the Buffer class.

#include "dca/io/buffer.hpp"

#include <complex>
#include "gtest/gtest.h"

struct PodClass {
  int a;
  double b;
};

struct NonPodClass {
  int a;
  std::unique_ptr<NonPodClass> next;
  float b;
};

using dca::io::Buffer;
Buffer& operator<<(Buffer& buffer, const NonPodClass& obj) {
  buffer << obj.a << obj.next << obj.b;
  return buffer;
}

Buffer& operator>>(Buffer& buffer, NonPodClass& obj) {
  buffer >> obj.a >>  obj.next >> obj.b;
  return buffer;
}

TEST(BufferTest, ReadWritePODTypes) {
  int a = 42;
  double b = 3.14;
  std::complex<float> c(1., 0.5);
  PodClass d{2, 2.7};

  dca::io::Buffer buffer;

  // Write data to buffer.
  EXPECT_EQ(0, buffer.size());
  buffer << a << b << c << d;
  constexpr int expected_size =
      sizeof(int) + sizeof(double) + sizeof(std::complex<float>) + sizeof(PodClass);
  EXPECT_EQ(expected_size, buffer.size());

  // Read data from buffer.
  int a_out;
  double b_out;
  std::complex<float> c_out;
  PodClass d_out;
  buffer >> a_out >> b_out >> c_out >> d_out;
  // The data is still present in the buffer.
  EXPECT_EQ(expected_size, buffer.size());
  EXPECT_EQ(expected_size, buffer.tellg());

  // Check output.
  EXPECT_EQ(a, a_out);
  EXPECT_EQ(b, b_out);
  EXPECT_EQ(c.real(), c_out.real());
  EXPECT_EQ(c.imag(), c_out.imag());
  EXPECT_EQ(d.a, d_out.a);
  EXPECT_EQ(d.b, d_out.b);
}

TEST(BufferTest, ReadWriteNonPODTypes) {
  NonPodClass obj;
  obj.a = 42;
  obj.next.reset(new  NonPodClass{21, nullptr, 3.2});
  obj.b = 3.14;

  dca::io::Buffer buffer;

  buffer << obj;

  NonPodClass obj_out;
  buffer >> obj_out;

  EXPECT_EQ(obj.a, obj_out.a);
  EXPECT_EQ(obj.b, obj_out.b);
  EXPECT_EQ(obj.next->a, obj_out.next->a);
  EXPECT_EQ(obj.next->b, obj_out.next->b);
  EXPECT_FALSE(obj_out.next->next);
}
