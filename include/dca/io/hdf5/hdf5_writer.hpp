// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// HDF5 writer.

#ifndef DCA_IO_HDF5_HDF5_WRITER_HPP
#define DCA_IO_HDF5_HDF5_WRITER_HPP

#include <complex>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "H5Cpp.h"

#include "dca/function/domains.hpp"
#include "dca/io/buffer.hpp"
#include "dca/function/function.hpp"
#include "dca/io/hdf5/hdf5_types.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/vector.hpp"

namespace dca {
namespace io {
// dca::io

class HDF5Writer {
public:
  typedef H5::H5File file_type;

public:
  HDF5Writer() : my_file(NULL), my_group(0), my_paths(0) {}
  ~HDF5Writer();

  bool is_reader() {
    return false;
  }
  bool is_writer() {
    return true;
  }

  H5::H5File& open_file(std::string file_name_ref, bool overwrite = true);
  void close_file();

  void open_group(std::string new_path);
  void close_group();

  std::string get_path();

  template <typename arbitrary_struct_t>
  static void to_file(const arbitrary_struct_t& arbitrary_struct, const std::string& file_name);

  template <typename scalar_type>
  void execute(const std::string& name, scalar_type value);

  template <typename scalar_type>
  void execute(const std::string& name, const std::pair<scalar_type, scalar_type>& value);

  template <typename scalar_type>
  void execute(const std::string& name, const std::vector<scalar_type>& value);

  template <typename scalar_type>
  void execute(const std::string& name, const std::vector<std::complex<scalar_type>>& value);

  void execute(const std::string& name, const std::string& value);

  void execute(const std::string& name, const std::vector<std::string>& value);

  template <typename scalar_type>
  void execute(const std::string& name, const std::vector<std::vector<scalar_type>>& value);

  template <typename domain_type>
  void execute(const std::string& name, const func::dmn_0<domain_type>& dmn);

  template <typename scalar_type, typename domain_type>
  void execute(const func::function<scalar_type, domain_type>& f);

  template <typename scalar_type, typename domain_type>
  void execute(const func::function<std::complex<scalar_type>, domain_type>& f);

  template <typename scalar_type, typename domain_type>
  void execute(const std::string& name, const func::function<scalar_type, domain_type>& f);

  template <typename scalar_type, typename domain_type>
  void execute(const std::string& name,
               const func::function<std::complex<scalar_type>, domain_type>& f);

  template <typename scalar_type>
  void execute(const std::string& name, const dca::linalg::Vector<scalar_type, dca::linalg::CPU>& A);

  template <typename scalar_type>
  void execute(const std::string& name,
               const dca::linalg::Vector<std::complex<scalar_type>, dca::linalg::CPU>& A);

  template <typename scalar_type>
  void execute(const std::string& name, const dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& A);

  template <typename scalar_type>
  void execute(const std::string& name,
               const dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU>& A);

  template <class T>
  void execute(const std::string& name, const std::unique_ptr<T>& obj);

  template <class T>
  void execute(const std::unique_ptr<T>& obj);

  void execute(const std::string& name, const io::Buffer& buffer){
      return execute(name, static_cast<io::Buffer::Container>(buffer));
  }

private:
  bool fexists(const char* filename);

  H5::H5File* my_file;

  hid_t file_id;

  std::vector<H5::Group*> my_group;
  std::vector<std::string> my_paths;
};

template <typename arbitrary_struct_t>
void HDF5Writer::to_file(const arbitrary_struct_t& arbitrary_struct, const std::string& file_name) {
  HDF5Writer wr_obj;
  wr_obj.open_file(file_name);
  arbitrary_struct.read_write(wr_obj);
  wr_obj.close_file();
}

template <typename scalar_type>
void HDF5Writer::execute(const std::string& name, scalar_type value) {
  H5::H5File& file = (*my_file);
  std::string path = get_path();

  hsize_t dims[1];

  H5::DataSet* dataset = NULL;
  H5::DataSpace* dataspace = NULL;

  {
    dims[0] = 1;
    dataspace = new H5::DataSpace(1, dims);

    std::string full_name = path + "/" + name;
    dataset = new H5::DataSet(
        file.createDataSet(full_name.c_str(), HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

    H5Dwrite(dataset->getId(), HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL,
             H5P_DEFAULT, &value);
  }

  delete dataset;
  delete dataspace;
}

template <typename scalar_type>
void HDF5Writer::execute(const std::string& name, const std::pair<scalar_type, scalar_type>& value) {
  H5::H5File& file = (*my_file);
  std::string path = get_path();

  hsize_t dims[1];

  H5::DataSet* dataset = NULL;
  H5::DataSpace* dataspace = NULL;

  {
    dims[0] = 2;
    dataspace = new H5::DataSpace(1, dims);

    std::string full_name = path + "/" + name;
    dataset = new H5::DataSet(
        file.createDataSet(full_name.c_str(), HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

    H5Dwrite(dataset->getId(), HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL,
             H5P_DEFAULT, &(value.first));
  }

  delete dataset;
  delete dataspace;
}

template <typename scalar_type>
void HDF5Writer::execute(const std::string& name,
                         const std::vector<scalar_type>& value)  //, H5File& file, std::string path)
{
  if (value.size() > 0) {
    H5::H5File& file = (*my_file);
    std::string path = get_path();

    hsize_t dims[1];

    H5::DataSet* dataset = NULL;
    H5::DataSpace* dataspace = NULL;

    {
      dims[0] = value.size();
      dataspace = new H5::DataSpace(1, dims);

      std::string full_name = path + "/" + name;
      dataset = new H5::DataSet(
          file.createDataSet(full_name.c_str(), HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

      H5Dwrite(dataset->getId(), HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL,
               H5P_DEFAULT, &value[0]);
    }

    delete dataset;
    delete dataspace;
  }
}

template <typename scalar_type>
void HDF5Writer::execute(const std::string& name,
                         const std::vector<std::complex<scalar_type>>& value) {
  H5::H5File& file = (*my_file);
  std::string path = get_path();

  hsize_t dims[2];

  H5::DataSet* dataset = NULL;
  H5::DataSpace* dataspace = NULL;

  {
    dims[0] = 2;
    dims[1] = value.size();
    dataspace = new H5::DataSpace(2, dims);

    std::string full_name = path + "/" + name;
    dataset = new H5::DataSet(
        file.createDataSet(full_name.c_str(), HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

    H5Dwrite(dataset->getId(), HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL,
             H5P_DEFAULT, &value[0]);
  }

  delete dataset;
  delete dataspace;
}

template <typename scalar_type>
void HDF5Writer::execute(const std::string& name, const std::vector<std::vector<scalar_type>>& value) {
  if (value.size() > 0) {
    H5::H5File& file = (*my_file);

    bool all_the_same_size = true;

    std::vector<size_t> dim(2, 0);
    {
      dim[0] = value.size();
      dim[1] = value[0].size();

      for (size_t l = 0; l < dim[0]; l++)
        if (dim[1] != value[l].size())
          all_the_same_size = false;
    }

    open_group(name);

    execute("equal-size", all_the_same_size);

    H5::DataSet* dataset = NULL;
    H5::DataSpace* dataspace = NULL;

    if (all_the_same_size) {
      execute("size", dim);

      hsize_t dims[2];

      {
        dims[0] = value.size();
        dims[1] = value[0].size();

        dataspace = new H5::DataSpace(2, dims);

        std::string full_name = get_path() + "/data";
        dataset = new H5::DataSet(file.createDataSet(
            full_name.c_str(), HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

        scalar_type* tmp = new scalar_type[dims[0] * dims[1]];

        /*
          for(int i=0; i<dims[0]; i++)
          for(int j=0; j<dims[1]; j++)
          tmp[i+j*dims[0]] = value[i][j];
        */

        // hdf5 has row-major ordering!
        for (hsize_t i = 0; i < dims[0]; i++)
          for (hsize_t j = 0; j < dims[1]; j++)
            tmp[i * dims[1] + j] = value[i][j];

        H5Dwrite(dataset->getId(), HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL,
                 H5P_DEFAULT, tmp);

        delete[] tmp;
      }
    }
    else {
      execute("size_i", dim[0]);

      dim.resize(0);
      for (size_t l = 0; l < value.size(); l++)
        dim.push_back(value[l].size());

      execute("size_j", dim);

      hsize_t dims[1];

      open_group("data");
      for (size_t l = 0; l < value.size(); l++) {
        dims[0] = value[l].size();
        dataspace = new H5::DataSpace(1, dims);

        std::stringstream ss;
        ss << get_path() << "/" << l;

        dataset = new H5::DataSet(file.createDataSet(
            ss.str().c_str(), HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

        H5Dwrite(dataset->getId(), HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL,
                 H5P_DEFAULT, &value[l]);
      }
      close_group();
    }

    delete dataset;
    delete dataspace;

    close_group();
  }
}

template <typename domain_type>
void HDF5Writer::execute(const std::string& name, const func::dmn_0<domain_type>& dmn) {
  open_group(name);

  execute("name", dmn.get_name());
  execute("elements", dmn.get_elements());

  close_group();
}

template <typename scalar_type, typename domain_type>
void HDF5Writer::execute(const func::function<scalar_type, domain_type>& f) {
  if (f.size() == 0)
    return;

  std::cout << "\t starts writing function : " << f.get_name() << "\n";

  execute(f.get_name(), f);
}

template <typename scalar_type, typename domain_type>
void HDF5Writer::execute(const func::function<std::complex<scalar_type>, domain_type>& f) {
  if (f.size() == 0)
    return;

  std::cout << "\t starts writing function : " << f.get_name() << "\n";

  execute(f.get_name(), f);
}

template <typename scalar_type, typename domain_type>
void HDF5Writer::execute(const std::string& name, const func::function<scalar_type, domain_type>& f) {
  if (f.size() == 0)
    return;

  H5::H5File& file = (*my_file);

  open_group(name);

  std::string new_path = get_path();

  {
    execute("name", f.get_name());

    {
      std::vector<int> vec(0);

      for (int l = 0; l < f.signature(); l++)
        vec.push_back(f[l]);

      execute("domain-sizes", vec);
    }

    {
      int N_dmns = f.signature();

      std::vector<hsize_t> dims(N_dmns);  // hsize_t dims[N_dmns];

      //      for(int l=0; l<N_dmns; l++)
      //        dims[l] = f[l];

      // be carefull --> HDF5 is by default row-major, while the function-class is column-major !
      for (int l = 0; l < N_dmns; l++)
        dims[N_dmns - 1 - l] = f[l];

      H5::DataSet* dataset = NULL;
      H5::DataSpace* dataspace = NULL;

      {
        dataspace = new H5::DataSpace(N_dmns, &dims[0]);

        std::string full_name = new_path + "/data";
        dataset = new H5::DataSet(file.createDataSet(
            full_name.c_str(), HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

        H5Dwrite(dataset->getId(), HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL,
                 H5P_DEFAULT, &f(0));
      }

      delete dataset;
      delete dataspace;
    }
  }

  close_group();
}

template <typename scalar_type, typename domain_type>
void HDF5Writer::execute(const std::string& name,
                         const func::function<std::complex<scalar_type>, domain_type>& f) {
  if (f.size() == 0)
    return;

  H5::H5File& file = (*my_file);

  open_group(name);

  {
    execute("name", f.get_name());

    {
      std::vector<int> vec(0);

      for (int l = 0; l < f.signature(); l++)
        vec.push_back(f[l]);

      execute("domain-sizes", vec);
    }

    {
      int N_dmns = f.signature();

      std::vector<hsize_t> dims(N_dmns + 1);

      // be carefull --> HDF5 is by default row-major, while the function-class is column-major !
      dims[N_dmns] = 2;
      for (int l = 0; l < N_dmns; l++)
        dims[N_dmns - 1 - l] = f[l];

      H5::DataSet* dataset = NULL;
      H5::DataSpace* dataspace = NULL;

      {
        dataspace = new H5::DataSpace(N_dmns + 1, &dims[0]);

        std::string full_name = get_path() + "/data";
        dataset = new H5::DataSet(file.createDataSet(
            full_name.c_str(), HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

        H5Dwrite(dataset->getId(), HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL,
                 H5P_DEFAULT, &f(0));
      }

      delete dataset;
      delete dataspace;
    }
  }

  close_group();
}

template <typename scalar_type>
void HDF5Writer::execute(const std::string& name,
                         const dca::linalg::Vector<scalar_type, dca::linalg::CPU>& V) {
  H5::H5File& file = (*my_file);

  open_group(name);

  execute("name", V.get_name());
  execute("size", V.get_size());

  hsize_t dims[1];

  dims[0] = V.size();

  H5::DataSet* dataset = NULL;
  H5::DataSpace* dataspace = NULL;

  {
    dataspace = new H5::DataSpace(1, dims);

    std::string full_name = get_path() + "/data";
    dataset = new H5::DataSet(
        file.createDataSet(full_name.c_str(), HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

    H5Dwrite(dataset->getId(), HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL,
             H5P_DEFAULT, &V[0]);
  }

  delete dataset;
  delete dataspace;

  close_group();
}

template <typename scalar_type>
void HDF5Writer::execute(const std::string& name,
                         const dca::linalg::Vector<std::complex<scalar_type>, dca::linalg::CPU>& V) {
  H5::H5File& file = (*my_file);

  open_group(name);

  execute("name", V.get_name());
  execute("size", V.size());

  hsize_t dims[2];

  dims[0] = 2;
  dims[1] = V.size();

  H5::DataSet* dataset = NULL;
  H5::DataSpace* dataspace = NULL;

  {
    dataspace = new H5::DataSpace(2, dims);

    std::string full_name = get_path() + "/data";
    dataset = new H5::DataSet(
        file.createDataSet(full_name.c_str(), HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

    H5Dwrite(dataset->getId(), HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL,
             H5P_DEFAULT, &V[0]);
  }

  delete dataset;
  delete dataspace;

  close_group();
}

template <typename scalar_type>
void HDF5Writer::execute(const std::string& name,
                         const dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& A) {
  H5::H5File& file = (*my_file);

  open_group(name);

  {
    execute("name", A.get_name());

    execute("current-size", A.size());
    execute("global-size", A.capacity());
  }

  hsize_t dims[2];

  dims[0] = A.capacity().first;
  dims[1] = A.capacity().second;

  H5::DataSet* dataset = NULL;
  H5::DataSpace* dataspace = NULL;

  {
    dataspace = new H5::DataSpace(2, dims);

    std::string full_name = get_path() + "/data";
    dataset = new H5::DataSet(
        file.createDataSet(full_name.c_str(), HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

    H5Dwrite(dataset->getId(), HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL,
             H5P_DEFAULT, &A(0, 0));
  }

  delete dataset;
  delete dataspace;

  close_group();
}

template <typename scalar_type>
void HDF5Writer::execute(const std::string& name,
                         const dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU>& A) {
  H5::H5File& file = (*my_file);

  open_group(name);

  {
    execute("name", A.get_name());

    std::vector<int> vec(2, 0);

    vec[0] = A.size().first;
    vec[1] = A.size().second;
    execute("current-size", vec);

    vec[0] = A.capacity().first;
    vec[1] = A.capacity().second;
    execute("global-size", vec);
  }

  hsize_t dims[3];

  dims[0] = 2;
  dims[1] = A.capacity().first;
  dims[2] = A.capacity().second;

  H5::DataSet* dataset = NULL;
  H5::DataSpace* dataspace = NULL;

  {
    dataspace = new H5::DataSpace(3, dims);

    std::string full_name = get_path() + "/data";
    dataset = new H5::DataSet(
        file.createDataSet(full_name.c_str(), HDF5_TYPE<scalar_type>::get_PredType(), *dataspace));

    H5Dwrite(dataset->getId(), HDF5_TYPE<scalar_type>::get(), dataspace->getId(), H5S_ALL,
             H5P_DEFAULT, &A(0, 0));
  }

  delete dataset;
  delete dataspace;

  close_group();
}

template <class T>
void HDF5Writer::execute(const std::string& name, const std::unique_ptr<T>& obj) {
  if (obj)
    execute(name, *obj);
}

template <class T>
void HDF5Writer::execute(const std::unique_ptr<T>& obj) {
  if (obj)
    execute(*obj);
}

}  // io
}  // dca

#endif  // DCA_IO_HDF5_HDF5_WRITER_HPP
