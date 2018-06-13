// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// JSON writer.

#ifndef DCA_IO_JSON_JSON_WRITER_HPP
#define DCA_IO_JSON_JSON_WRITER_HPP

#include <complex>
#include <map>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/json/json_parser/json_context.hpp"
#include "dca/io/json/json_parser/whatever.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/vector.hpp"

namespace dca {
namespace io {
// dca::io::

class JSONWriter {
public:
  typedef std::stringstream file_type;

  typedef Whatever JsonAccessor;
  typedef JSON_context JsonDataType;

public:
  JSONWriter();

  bool is_reader() const {
    return false;
  }
  bool is_writer() const {
    return true;
  }

  file_type& open_file(const std::string& file_name_ref, bool overwrite = true);
  void close_file();

  void open_group(const std::string& new_path);
  void close_group();

  std::string get_path();

  template <typename arbitrary_struct_t>
  static void to_file(arbitrary_struct_t& arbitrary_struct, const std::string& file_name);

  template <typename scalartype>
  void execute(const std::string& name, const scalartype& value);

  template <typename s_t_0, typename s_t_1>
  void execute(const std::string& name, const std::pair<s_t_0, s_t_1>& value);

  template <typename scalartype>
  void execute(const std::string& name, const std::vector<scalartype>& value);

  template <typename scalartype>
  void execute(const std::string& name, const std::vector<std::vector<scalartype>>& value);

  void execute(const std::string& name, const std::string& value);

  void execute(const std::string& name, const std::vector<std::string>& value);

  template <typename domain_type>
  void execute(const std::string& name, const func::dmn_0<domain_type>& dmn);

  template <typename scalar_type, typename domain_type>
  void execute(const func::function<scalar_type, domain_type>& f);

  template <typename scalar_type, typename domain_type>
  void execute(const std::string& name, const func::function<scalar_type, domain_type>& f);

  template <typename scalar_type, typename domain_type>
  void execute(const func::function<std::complex<scalar_type>, domain_type>& f);

  template <typename scalar_type, typename domain_type>
  void execute(const std::string& name,
               const func::function<std::complex<scalar_type>, domain_type>& f);

  template <typename scalar_type>
  void execute(const std::string& name, const linalg::Vector<scalar_type, linalg::CPU>& V);

  template <typename scalar_type>
  void execute(const std::string& name,
               const linalg::Vector<std::complex<scalar_type>, linalg::CPU>& V);

  template <typename scalar_type>
  void execute(const std::string& name, const linalg::Matrix<scalar_type, linalg::CPU>& A);

  template <typename scalar_type>
  void execute(const std::string& name,
               const linalg::Matrix<std::complex<scalar_type>, linalg::CPU>& A);

  template <class T>
  void execute(const std::string& name,
	       const std::unique_ptr<T>& obj);

  template <class T>
  void execute(const std::unique_ptr<T>& obj);

  template <class stream_type>
  static void execute(stream_type& ss, const JsonAccessor& parseResult);

private:
  std::stringstream ss;

  std::string file_name;

  std::string path;

  std::vector<int> elements_in_group;
};

template <typename arbitrary_struct_t>
void JSONWriter::to_file(arbitrary_struct_t& arbitrary_struct, const std::string& file_name) {
  JSONWriter wr_obj;

  wr_obj.open_file(file_name);
  arbitrary_struct.read_write(wr_obj);
  wr_obj.close_file();
}

template <typename scalartype>
void JSONWriter::execute(
    const std::string& name,
    const scalartype& value)  //, file_type& ss), std::string path, bool is_ending)
{
  if (elements_in_group.back() != 0)
    ss << ",\n";

  ss << get_path() << "\"" << name << "\" : " << value;

  elements_in_group.back() += 1;
}

template <typename s_t_0, typename s_t_1>
void JSONWriter::execute(
    const std::string& name,
    const std::pair<s_t_0, s_t_1>& value)  //, file_type& ss), std::string path, bool is_ending)
{
  if (elements_in_group.back() != 0)
    ss << ",\n";

  ss << get_path() << "\"" << name << "\" : [" << value.first << ", " << value.second << "]";

  elements_in_group.back() += 1;
}

template <typename scalartype>
void JSONWriter::execute(
    const std::string& name,
    const std::vector<scalartype>& value)  //, file_type& ss)//, std::string path, bool is_ending)
{
  if (elements_in_group.back() != 0)
    ss << ",\n";

  ss << get_path() << "\"" << name << "\" : [";

  for (size_t i = 0; i < value.size(); i++) {
    ss << value[i];

    if (i < value.size() - 1)
      ss << ", ";
  }

  ss << "]";

  elements_in_group.back() += 1;
}

template <typename scalartype>
void JSONWriter::execute(const std::string& name, const std::vector<std::vector<scalartype>>& value) {
  if (elements_in_group.back() != 0)
    ss << ",\n";

  ss << get_path() << "\"" << name << "\" : [";

  std::string indent = "";
  indent.resize(name.size() + 6, ' ');

  for (size_t i = 0; i < value.size(); i++) {
    ss << "[";

    for (size_t j = 0; j < value[i].size(); j++) {
      ss << value[i][j];

      if (j < value[i].size() - 1)
        ss << ", ";
    }
    ss << "]";

    if (i < value.size() - 1)
      ss << ",\n" << get_path() << indent;
  }
  ss << "]";

  elements_in_group.back() += 1;
}

template <typename domain_type>
void JSONWriter::execute(const std::string& name, const func::dmn_0<domain_type>& dmn) {
  open_group(name);

  execute("name", dmn.get_name());

  execute("elements", dmn.get_elements());

  close_group();

  elements_in_group.back() += 1;
}

template <typename scalar_type, typename domain_type>
void JSONWriter::execute(const func::function<scalar_type, domain_type>& f) {
  std::cout << "\t starts writing function : " << f.get_name() << "\n";

  execute(f.get_name(), f);
}

template <typename scalar_type, typename domain_type>
void JSONWriter::execute(const std::string& name, const func::function<scalar_type, domain_type>& f) {
  open_group(name);

  execute("name", f.get_name());

  {
    std::vector<int> vec(0);

    for (int l = 0; l < f.signature(); l++)
      vec.push_back(f[l]);

    execute("domain-sizes", vec);  //, file, new_path);
  }

  ss << ",\n\n" << get_path() << "\"data\" : [";

  //     ss << std::fixed;
  //     ss.precision(16);

  int* subind = new int[f.signature()];
  for (int i = 0; i < f.size(); i++) {
    ss << "[";

    int dummy = i;
    f.linind_2_subind(dummy, subind);
    for (int j = 0; j < f.signature(); j++) {
      ss << (subind[j]);
      ss << ", ";
    }

    ss << f(i);

    if (i == f.size() - 1)
      ss << "]";
    else
      ss << "],\n" << get_path() << "          ";
  }
  delete[] subind;

  ss << "]";

  close_group();

  elements_in_group.back() += 1;
}

template <typename scalar_type, typename domain_type>
void JSONWriter::execute(const func::function<std::complex<scalar_type>, domain_type>& f) {
  std::cout << "\t starts writing function : " << f.get_name() << "\n";

  execute(f.get_name(), f);
}

template <typename scalar_type, typename domain_type>
void JSONWriter::execute(const std::string& name,
                         const func::function<std::complex<scalar_type>, domain_type>& f) {
  open_group(name);

  execute("name", f.get_name());

  {
    std::vector<int> vec(0);

    for (int l = 0; l < f.signature(); l++)
      vec.push_back(f[l]);

    execute("domain-sizes", vec);
  }

  ss << ",\n\n" << get_path() << "\"data\" : [";

  //     ss << std::fixed;
  //     ss.precision(16);

  int* subind = new int[f.signature()];
  for (int i = 0; i < f.size(); i++) {
    ss << "[";

    int dummy = i;
    f.linind_2_subind(dummy, subind);
    for (int j = 0; j < f.signature(); j++) {
      ss << (subind[j]);
      ss << ", ";
    }

    ss << real(f(i)) << ", " << imag(f(i));

    if (i == f.size() - 1)
      ss << "]";
    else
      ss << "],\n" << get_path() << "          ";
  }
  delete[] subind;

  ss << "]";

  close_group();

  elements_in_group.back() += 1;
}

template <typename scalar_type>
void JSONWriter::execute(const std::string& name, const linalg::Vector<scalar_type, linalg::CPU>& V) {
  open_group(name);

  execute("name", V.get_name());

  execute("current-size", V.size());
  execute("global-size", V.capacity());

  ss << ",\n\n" << get_path() << "\"data\" : [";

  for (int j = 0; j < V.size() - 1; j++)
    ss << V[j] << ", ";

  ss << V[V.size() - 1] << "]";

  close_group();

  elements_in_group.back() += 1;
}

template <typename scalar_type>
void JSONWriter::execute(const std::string& name,
                         const linalg::Vector<std::complex<scalar_type>, linalg::CPU>& V) {
  open_group(name);

  execute("name", V.get_name());

  execute("current-size", V.size());
  execute("global-size", V.capacity());

  ss << ",\n\n" << get_path() << "\"data-real\" : [";
  for (int j = 0; j < V.size() - 1; j++)
    ss << real(V[j]) << ", ";

  ss << real(V[V.size() - 1]) << "]";
  ss << ",\n\n" << get_path() << "\"data-imag\" : [";
  for (int j = 0; j < V.size() - 1; j++)
    ss << imag(V[j]) << ", ";

  ss << imag(V[V.size() - 1]) << "]";

  close_group();

  elements_in_group.back() += 1;
}

template <typename scalar_type>
void JSONWriter::execute(const std::string& name, const linalg::Matrix<scalar_type, linalg::CPU>& A) {
  open_group(name);

  execute("name", A.get_name());

  execute("current-size", A.size());
  execute("global-size", A.capacity());

  ss << ",\n\n" << get_path() << "\"data\" : [";

  for (int i = 0; i < A.size().first; i++) {
    ss << "[";

    for (int j = 0; j < A.size().second - 1; j++)
      ss << A(i, j) << ", ";

    if (i == A.size().first - 1)
      ss << A(i, A.size().second - 1) << "]";
    else
      ss << A(i, A.size().second - 1) << "],\n" << get_path() << "          ";
  }

  ss << "]";

  close_group();

  elements_in_group.back() += 1;
}

template <typename scalar_type>
void JSONWriter::execute(const std::string& name,
                         const linalg::Matrix<std::complex<scalar_type>, linalg::CPU>& A) {
  open_group(name);

  execute("name", A.get_name());

  execute("current-size", A.size());
  execute("current-size", A.capacity());

  ss << ",\n\n" << get_path() << "\"data-real\" : [";

  for (int i = 0; i < A.size().first; i++) {
    ss << "[";

    for (int j = 0; j < A.size().second - 1; j++)
      ss << real(A(i, j)) << ", ";

    if (i == A.size().first - 1)
      ss << real(A(i, A.size().second - 1)) << "]";
    else
      ss << real(A(i, A.size().second - 1)) << "],\n" << get_path() << "               ";
  }

  ss << "]";

  ss << ",\n\n" << get_path() << "\"data-imag\" : [";

  for (int i = 0; i < A.size().first; i++) {
    ss << "[";

    for (int j = 0; j < A.size().second - 1; j++)
      ss << imag(A(i, j)) << ", ";

    if (i == A.size().first - 1)
      ss << imag(A(i, A.size().second - 1)) << "]";
    else
      ss << imag(A(i, A.size().second - 1)) << "],\n" << get_path() << "               ";
  }

  ss << "]";

  close_group();

  elements_in_group.back() += 1;
}

template <class T>
void JSONWriter::execute(const std::string& name,
			 const std::unique_ptr<T>& obj) {
  if (obj)
    execute(name, *obj);
}

template <class T>
void JSONWriter::execute(const std::unique_ptr<T>& obj) {
  if (obj)
    execute(*obj);
}

template <class stream_type>
void JSONWriter::execute(stream_type& os, const JsonAccessor& w) {
  static int level = -1;

  typedef typename std::map<std::wstring, Whatever>::const_iterator WhateverMapItr;
  switch (w.type) {
    case WHATEVER_MAT: {
      std::string wfilename(w.filename.begin(), w.filename.end());
      os << "{ 'fileName': '" << wfilename << "'"
         << ", 'startPos': " << w.startPos << ", 'endPos': " << w.endPos << "}";
      break;
    }

    case WHATEVER_MAP: {
      level += 1;

      os << "\n";
      for (int l = 0; l < level; l++)
        os << "\t";
      os << "{\n";

      int index = 0;
      for (WhateverMapItr itr = w.whateverMap.begin(); itr != w.whateverMap.end(); itr++) {
        const std::wstring& wkey = (*itr).first;
        const std::string key(wkey.begin(), wkey.end());

        for (int l = 0; l < level; l++)
          os << "\t";

        os << "\"" << key << "\" : ";  // << (*itr).second;

        JSONWriter::execute(os, (*itr).second);

        if (int(w.whateverMap.size()) == index + 1)
          os << "";
        else
          os << ",\n";

        index += 1;
      }

      os << "\n";
      for (int l = 0; l < level; l++)
        os << "\t";
      os << "}";

      level -= 1;

      break;
    }

    case WHATEVER_VECTOR: {
      os << "[";
      for (size_t i = 0; i < w.whateverVector.size(); i++) {
        JSONWriter::execute(os, w.whateverVector[i]);
        if (i < w.whateverVector.size() - 1)
          os << ", ";
      }

      os << "]";

      break;
    }

    case WHATEVER_MATRIX:
      os << "WHATEVER_MATRIX";
      break;

    case WHATEVER_STRING: {
      const std::string tmp(w.valueString.begin(), w.valueString.end());
      os << "\"" << tmp << "\"";
    } break;

    case WHATEVER_INTEGER:
      os << w.whateverInteger;
      break;

    case WHATEVER_DOUBLE:
      os << w.whateverDouble;
      break;

    case WHATEVER_BOOL:
      os << w.whateverBool;
      break;

    case WHATEVER_UNKNOWN:
      os << "WHATEVER_UNKNOWN";
      break;

    default:
      throw std::logic_error("typeName given wrong type");
  }
}

}  // io
}  // dca

#endif  // DCA_IO_JSON_JSON_WRITER_HPP
