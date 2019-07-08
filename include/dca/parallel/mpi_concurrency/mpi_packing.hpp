// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class provides an interface for packing and unpacking with MPI.

#ifndef DCA_PARALLEL_MPI_CONCURRENCY_MPI_PACKING_HPP
#define DCA_PARALLEL_MPI_CONCURRENCY_MPI_PACKING_HPP

#include <string>
#include <vector>
#include <mpi.h>
#include "dca/function/function.hpp"
#include "dca/parallel/mpi_concurrency/mpi_processor_grouping.hpp"
#include "dca/parallel/mpi_concurrency/mpi_type_map.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

class MPIPacking : public virtual MPIProcessorGrouping {
public:
  MPIPacking() {}

  template <typename scalar_type>
  int get_buffer_size(scalar_type item) const;
  template <typename scalar_type>
  int get_buffer_size(const std::basic_string<scalar_type>& str) const;
  template <typename scalar_type>
  int get_buffer_size(const std::vector<scalar_type>& v) const;
  template <typename scalar_type>
  int get_buffer_size(const std::vector<std::basic_string<scalar_type>>& v) const;
  template <typename scalar_type>
  int get_buffer_size(const std::vector<std::vector<scalar_type>>& v) const;
  template <typename scalar_type, class dmn_type>
  int get_buffer_size(const func::function<scalar_type, dmn_type>& f) const;

  template <typename scalar_type>
  void pack(char* buffer, int size, int& off_set, scalar_type item) const;
  template <typename scalar_type>
  void pack(char* buffer, int size, int& off_set, const std::basic_string<scalar_type>& str) const;
  template <typename scalar_type>
  void pack(char* buffer, int size, int& off_set, const std::vector<scalar_type>& v) const;
  template <typename scalar_type>
  void pack(char* buffer, int size, int& off_set,
            const std::vector<std::basic_string<scalar_type>>& v) const;
  template <typename scalar_type>
  void pack(char* buffer, int size, int& off_set,
            const std::vector<std::vector<scalar_type>>& v) const;
  template <typename scalar_type, class dmn_type>
  void pack(char* buffer, int size, int& off_set,
            const func::function<scalar_type, dmn_type>& f) const;

  template <typename scalar_type>
  void unpack(char* buffer, int size, int& off_set, scalar_type& item) const;
  template <typename scalar_type>
  void unpack(char* buffer, int size, int& off_set, std::basic_string<scalar_type>& str) const;
  template <typename scalar_type>
  void unpack(char* buffer, int size, int& off_set, std::vector<scalar_type>& v) const;
  template <typename scalar_type>
  void unpack(char* buffer, int size, int& off_set,
              std::vector<std::basic_string<scalar_type>>& v) const;
  template <typename scalar_type>
  void unpack(char* buffer, int size, int& off_set, std::vector<std::vector<scalar_type>>& v) const;
  template <typename scalar_type, class dmn_type>
  void unpack(char* buffer, int size, int& off_set, func::function<scalar_type, dmn_type>& f) const;
};

template <typename scalar_type>
int MPIPacking::get_buffer_size(scalar_type /*item*/) const {
  int size(0);

  MPI_Pack_size(1, MPITypeMap<scalar_type>::value(), MPIProcessorGrouping::get(), &size);

  return size;
}

template <typename scalar_type>
int MPIPacking::get_buffer_size(const std::basic_string<scalar_type>& str) const {
  /*
    int result = get_buffer_size(str.size());

    int size(0);
    MPI_Pack_size(str.size(), MPI_CHAR, MPIProcessorGrouping::get(), &size);

    result += size;

    return result;
  */

  int result = get_buffer_size(str.size());

  int count = str.size();

  {
    int size(0);
    MPI_Pack_size(count, MPITypeMap<scalar_type>::value(), MPIProcessorGrouping::get(), &size);

    result += size;
  }

  return result;
}

template <typename scalar_type>
int MPIPacking::get_buffer_size(const std::vector<scalar_type>& v) const {
  int result = get_buffer_size(v.size());

  int count = v.size();

  {
    int size(0);
    MPI_Pack_size(count, MPITypeMap<scalar_type>::value(), MPIProcessorGrouping::get(), &size);

    result += size;
  }

  return result;
}

template <typename scalar_type>
int MPIPacking::get_buffer_size(const std::vector<std::basic_string<scalar_type>>& v) const {
  int result = 0;

  std::vector<int> tmp_sizes(0);
  for (int i = 0; i < v.size(); i++)
    tmp_sizes.push_back(v[i].size());

  std::vector<scalar_type> tmp_chars(0);
  for (int i = 0; i < v.size(); i++)
    for (int j = 0; j < v[i].size(); j++)
      tmp_chars.push_back(v[i][j]);

  result += get_buffer_size(v.size());
  result += get_buffer_size(tmp_sizes);
  result += get_buffer_size(tmp_chars);

  return result;
}

template <typename scalar_type>
int MPIPacking::get_buffer_size(const std::vector<std::vector<scalar_type>>& v) const {
  std::vector<scalar_type> tmp;

  tmp.push_back(v.size());
  for (int i = 0; i < v.size(); i++) {
    tmp.push_back(v[i].size());

    for (int j = 0; j < v[i].size(); j++)
      tmp.push_back(v[i][j]);
  }

  return get_buffer_size(tmp);
}

template <typename scalar_type, class dmn_type>
int MPIPacking::get_buffer_size(const func::function<scalar_type, dmn_type>& f) const {
  int result = get_buffer_size(f.size());

  int count = f.size();

  {
    int size = 0;
    MPI_Pack_size(count, MPITypeMap<scalar_type>::value(), MPIProcessorGrouping::get(), &size);

    result += size;
  }

  return result;
}

template <typename scalar_type>
void MPIPacking::pack(char* buffer, int size, int& off_set, scalar_type item) const {
  const scalar_type* tPtr(&item);

  MPI_Pack(tPtr, 1, MPITypeMap<scalar_type>::value(), buffer, size, &off_set,
           MPIProcessorGrouping::get());
}

template <typename scalar_type>
void MPIPacking::pack(char* buffer, int size, int& off_set,
                      const std::basic_string<scalar_type>& str) const {
  /*
  // pack the string's length
  int stringSize(str.size());
  pack(buffer,size,off_set,stringSize);

  MPI_Pack(const_cast<char*>(str.c_str()), stringSize, MPI_CHAR,
  buffer, size, &off_set,
  MPIProcessorGrouping::get());
  */

  // Pack the vector length
  int vectorSize(str.size());
  pack(buffer, size, off_set, vectorSize);

  MPI_Pack(&str[0], vectorSize, MPITypeMap<scalar_type>::value(), buffer, size, &off_set,
           MPIProcessorGrouping::get());
}

template <typename scalar_type>
void MPIPacking::pack(char* buffer, int size, int& off_set, const std::vector<scalar_type>& v) const {
  // Pack the vector length
  int vectorSize(v.size());
  pack(buffer, size, off_set, vectorSize);

  MPI_Pack(&v[0], vectorSize, MPITypeMap<scalar_type>::value(), buffer, size, &off_set,
           MPIProcessorGrouping::get());
}

template <typename scalar_type>
void MPIPacking::pack(char* buffer, int size, int& off_set,
                      const std::vector<std::basic_string<scalar_type>>& v) const {
  std::vector<int> tmp_sizes(0);
  for (int i = 0; i < v.size(); i++)
    tmp_sizes.push_back(v[i].size());

  std::vector<scalar_type> tmp_chars(0);
  for (int i = 0; i < v.size(); i++)
    for (int j = 0; j < v[i].size(); j++)
      tmp_chars.push_back(v[i][j]);

  {
    // Pack the vector length
    int v_size = v.size();
    pack(buffer, size, off_set, v_size);

    pack(buffer, size, off_set, tmp_sizes);
    pack(buffer, size, off_set, tmp_chars);
  }
}

template <typename scalar_type>
void MPIPacking::pack(char* buffer, int size, int& off_set,
                      const std::vector<std::vector<scalar_type>>& v) const {
  std::vector<scalar_type> tmp;

  tmp.push_back(v.size());
  for (int i = 0; i < v.size(); i++) {
    tmp.push_back(v[i].size());

    for (int j = 0; j < v[i].size(); j++)
      tmp.push_back(v[i][j]);
  }

  pack(buffer, size, off_set, tmp);
}

template <typename scalar_type, class dmn_type>
void MPIPacking::pack(char* buffer, int size, int& off_set,
                      const func::function<scalar_type, dmn_type>& f) const {
  // Pack the vector length
  auto function_size = f.size();
  pack(buffer, size, off_set, function_size);

  MPI_Pack(f.values(), function_size, MPITypeMap<scalar_type>::value(), buffer, size, &off_set,
           MPIProcessorGrouping::get());
}

template <typename scalar_type>
void MPIPacking::unpack(char* buffer, int size, int& off_set, scalar_type& item) const {
  scalar_type tmp;

  MPI_Unpack(buffer, size, &off_set, &tmp, 1, MPITypeMap<scalar_type>::value(),
             MPIProcessorGrouping::get());

  item = tmp;
}

template <typename scalar_type>
void MPIPacking::unpack(char* buffer, int size, int& off_set,
                        std::basic_string<scalar_type>& str) const {
  /*
  // Unpack the string length
  int stringSize(0);
  unpack(buffer, size, off_set, stringSize);

  char stringBuffer[stringSize];
  MPI_Unpack(buffer,
  size,
  &off_set,
  stringBuffer,
  stringSize,
  MPI_CHAR,
  MPIProcessorGrouping::get());

  std::basic_string<scalar_type> tmp(stringBuffer,stringSize);

  str = tmp;
  */

  // UnPack the vector length
  int vectorSize(0);
  unpack(buffer, size, off_set, vectorSize);

  str.resize(vectorSize);

  // UnPack the vector
  MPI_Unpack(buffer, size, &off_set, static_cast<scalar_type*>(&str[0]), 1 * vectorSize,
             MPITypeMap<scalar_type>::value(), MPIProcessorGrouping::get());
}

template <typename scalar_type>
void MPIPacking::unpack(char* buffer, int size, int& off_set, std::vector<scalar_type>& v) const {
  // UnPack the vector length
  int vectorSize(0);
  unpack(buffer, size, off_set, vectorSize);

  v.resize(vectorSize);

  // UnPack the vector
  MPI_Unpack(buffer, size, &off_set, static_cast<scalar_type*>(&v[0]), 1 * vectorSize,
             MPITypeMap<scalar_type>::value(), MPIProcessorGrouping::get());
}

template <typename scalar_type>
void MPIPacking::unpack(char* buffer, int size, int& off_set,
                        std::vector<std::basic_string<scalar_type>>& v) const {
  int v_size;

  std::vector<int> tmp_sizes(0);
  std::vector<scalar_type> tmp_chars(0);

  {
    unpack(buffer, size, off_set, v_size);

    unpack(buffer, size, off_set, tmp_sizes);
    unpack(buffer, size, off_set, tmp_chars);
  }

  {
    v.resize(v_size);

    for (int i = 0; i < v.size(); i++)
      v[i].resize(tmp_sizes[i]);

    int index = 0;
    for (int i = 0; i < v.size(); i++) {
      for (int j = 0; j < v[i].size(); j++) {
        v[i][j] = tmp_chars[index];

        index += 1;
      }
    }
  }
}

template <typename scalar_type>
void MPIPacking::unpack(char* buffer, int size, int& off_set,
                        std::vector<std::vector<scalar_type>>& v) const {
  std::vector<scalar_type> tmp;

  unpack(buffer, size, off_set, tmp);

  v.resize(0);

  for (int i = 1; i < tmp.size();) {
    int v_size = tmp[i];

    std::vector<scalar_type> tmp_i(0);
    for (int j = i + 1; j < i + 1 + v_size; j++)
      tmp_i.push_back(tmp[j]);

    v.push_back(tmp_i);

    i += (v_size + 1);
  }
}

template <typename scalar_type, class dmn_type>
void MPIPacking::unpack(char* buffer, int size, int& off_set,
                        func::function<scalar_type, dmn_type>& f) const {
  // UnPack the vector length
  std::size_t function_size = 0;
  unpack(buffer, size, off_set, function_size);

  // UnPack the vector
  MPI_Unpack(buffer, size, &off_set, f.values(), function_size, MPITypeMap<scalar_type>::value(),
             MPIProcessorGrouping::get());
}

}  // namespace parallel
}  // namespace dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_PACKING_HPP
