//-*-C++-*-

#ifndef JSON_READER_HEADER_H
#define JSON_READER_HEADER_H
#include "src/comp_library/IO_library/JSON/JSON_PARSER/what_ever.h"
#include "comp_library/IO_library/JSON/JSON_PARSER/JSON_parser.h"
#include "comp_library/IO_library/JSON/JSON_PARSER/default_context.h"

namespace IO {

/*!
 * \author: based on json-cpp
 */
template <>
class reader<IO::JSON> {
public:
  typedef std::stringstream file_type;

  typedef JSONPARSER::Whatever JsonAccessor;
  typedef JSONPARSER::JSON_context JsonDataType;

public:
  reader();
  ~reader();

  bool is_reader() {
    return true;
  }
  bool is_writer() {
    return false;
  }

  std::string get_current_file_name();

  void open_file(std::string file_name);
  void close_file();

  void open_group(std::string name);
  void close_group();

  std::string get_path();

  template <typename arbitrary_struct_t>
  static void from_file(arbitrary_struct_t& arbitrary_struct, std::string file_name);

  template <typename scalartype>
  void execute(std::string name, scalartype& value);

  template <typename domain_type>
  void execute(std::string name, dmn_0<domain_type>& dmn);

  template <typename scalartype, typename domain_type>
  void execute(FUNC_LIB::function<scalartype, domain_type>& f);

  template <typename scalartype, typename domain_type>
  void execute(std::string name, FUNC_LIB::function<scalartype, domain_type>& f);

  template <typename scalar_type>
  void execute(std::string name, LIN_ALG::vector<scalar_type, LIN_ALG::CPU>& A);

  template <typename scalar_type>
  void execute(std::string name, LIN_ALG::matrix<scalar_type, LIN_ALG::CPU>& A);

private:
  template <typename scalartype>
  void execute(std::string name, scalartype& value, const JsonAccessor& current_result, size_t index);

  template <typename scalartype, typename domain_type>
  void execute(std::string name, FUNC_LIB::function<scalartype, domain_type>& f,
               const JsonAccessor& current_result, size_t index);

  template <typename scalartype, typename domain_type>
  void execute(std::string name, FUNC_LIB::function<std::complex<scalartype>, domain_type>& f,
               const JsonAccessor& current_result, size_t index);

  template <typename scalar_type>
  void execute(std::string name, LIN_ALG::vector<scalar_type, LIN_ALG::CPU>& A,
               const JsonAccessor& current_result, size_t index);

  template <typename scalar_type>
  void execute(std::string name, LIN_ALG::vector<std::complex<scalar_type>, LIN_ALG::CPU>& A,
               const JsonAccessor& current_result, size_t index);

  template <typename scalar_type>
  void execute(std::string name, LIN_ALG::matrix<scalar_type, LIN_ALG::CPU>& A,
               const JsonAccessor& current_result, size_t index);

  template <typename scalar_type>
  void execute(std::string name, LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU>& A,
               const JsonAccessor& current_result, size_t index);

private:
  void parse(std::string& file_name);

  const JsonAccessor& operator[](const std::string& key) const;

  const JsonAccessor& operator[](int index);

  const JsonAccessor& get_JSON_tree();

private:
  std::string current_file_name;

  JSONPARSER::JSON_parser<JsonDataType> parser;
  const JsonAccessor& parse_result;

  std::vector<std::string> my_paths;
};

reader<IO::JSON>::reader()
    : current_file_name("input.json"),

      parser(),
      parse_result(parser.get_JSON_tree().result),

      my_paths(0) {}

reader<IO::JSON>::~reader() {}

std::string reader<IO::JSON>::get_current_file_name() {
  return current_file_name;
}

void reader<IO::JSON>::open_file(std::string file_name) {
  current_file_name = file_name;

  parse(file_name);
}

void reader<IO::JSON>::close_file() {}

void reader<IO::JSON>::open_group(std::string name) {
  my_paths.push_back(name);
}

void reader<IO::JSON>::close_group() {
  my_paths.pop_back();
}

std::string reader<IO::JSON>::get_path() {
  std::string path = "/";

  for (size_t i = 0; i < my_paths.size(); i++) {
    path = path + my_paths[i];

    if (i < my_paths.size() - 1)
      path = path + "/";
  }

  return path;
}

template <typename arbitrary_struct_t>
void reader<IO::JSON>::from_file(arbitrary_struct_t& arbitrary_struct, std::string file_name) {
  reader<IO::JSON> reader_obj;

  reader_obj.open_file(file_name);

  arbitrary_struct.read_write(reader_obj);

  reader_obj.close_file();
}

template <typename scalartype>
void reader<IO::JSON>::execute(std::string name, scalartype& value) {
  size_t index = 0;

  if (index == my_paths.size())
    value <= parse_result[name];
  else {
    const JsonAccessor& new_result(parse_result[my_paths[index]]);

    execute(name, value, new_result, (++index));
  }
}

template <typename scalartype>
void reader<IO::JSON>::execute(std::string name, scalartype& value,
                               const JsonAccessor& current_result, size_t index) {
  if (index == my_paths.size())
    value <= current_result[name];
  else {
    const JsonAccessor& new_result(current_result[my_paths[index]]);

    execute(name, value, new_result, (++index));
  }
}

template <typename domain_type>
void reader<IO::JSON>::execute(std::string /*name*/, dmn_0<domain_type>& /*dmn*/) {}

template <typename scalartype, typename domain_type>
void reader<IO::JSON>::execute(FUNC_LIB::function<scalartype, domain_type>& f) {
  std::cout << "\t starts reading function : " << f.get_name() << "\n";

  execute(f.get_name(), f, parse_result, 0);
}

template <typename scalartype, typename domain_type>
void reader<IO::JSON>::execute(std::string name, FUNC_LIB::function<scalartype, domain_type>& f) {
  execute(name, f, parse_result, 0);
}

template <typename scalartype, typename domain_type>
void reader<IO::JSON>::execute(std::string name, FUNC_LIB::function<scalartype, domain_type>& f,
                               const JsonAccessor& current_result, size_t index) {
  if (index == my_paths.size()) {
    std::vector<std::vector<scalartype>> value;

    value <= current_result[name]["data"];

    int N_dmns = f.signature();
    for (size_t l = 0; l < value.size(); l++)
      f(l) = value[l][N_dmns];
  }
  else {
    const JsonAccessor& new_result(current_result[my_paths[index]]);

    execute(name, f, new_result, (++index));
  }
}

template <typename scalartype, typename domain_type>
void reader<IO::JSON>::execute(std::string name,
                               FUNC_LIB::function<std::complex<scalartype>, domain_type>& f,
                               const JsonAccessor& current_result, size_t index) {
  const std::complex<scalartype> I(0, 1);

  if (index == my_paths.size()) {
    std::vector<std::vector<scalartype>> value;

    value <= current_result[name]["data"];

    int N_dmns = f.signature();
    for (size_t l = 0; l < value.size(); l++)
      f(l) = value[l][N_dmns] + I * value[l][N_dmns + 1];
  }
  else {
    const JsonAccessor& new_result(current_result[my_paths[index]]);

    execute(name, f, new_result, (++index));
  }
}

template <typename scalar_type>
void reader<IO::JSON>::execute(std::string name, LIN_ALG::vector<scalar_type, LIN_ALG::CPU>& V) {
  execute(name, V, parse_result, 0);
}

template <typename scalar_type>
void reader<IO::JSON>::execute(std::string name, LIN_ALG::vector<scalar_type, LIN_ALG::CPU>& V,
                               const JsonAccessor& current_result, size_t index) {
  if (index == my_paths.size()) {
    std::vector<scalar_type> value;

    value <= current_result[name]["data"];

    for (size_t i = 0; i < value.size(); i++)
      V[i] = value[i];
  }
  else {
    const JsonAccessor& new_result(current_result[my_paths[index]]);

    execute(name, V, new_result, (++index));
  }
}

template <typename scalar_type>
void reader<IO::JSON>::execute(std::string name,
                               LIN_ALG::vector<std::complex<scalar_type>, LIN_ALG::CPU>& V,
                               const JsonAccessor& current_result, size_t index) {
  if (index == my_paths.size()) {
    std::vector<scalar_type> value;

    value <= current_result[name]["data-real"];

    for (size_t i = 0; i < value.size(); i++)
      real(V[i]) = value[i];

    value <= current_result[name]["data-imag"];

    for (size_t i = 0; i < value.size(); i++)
      imag(V[i]) = value[i];
  }
  else {
    const JsonAccessor& new_result(current_result[my_paths[index]]);

    execute(name, V, new_result, (++index));
  }
}

template <typename scalar_type>
void reader<IO::JSON>::execute(std::string name, LIN_ALG::matrix<scalar_type, LIN_ALG::CPU>& A) {
  execute(name, A, parse_result, 0);
}

template <typename scalar_type>
void reader<IO::JSON>::execute(std::string name, LIN_ALG::matrix<scalar_type, LIN_ALG::CPU>& A,
                               const JsonAccessor& current_result, size_t index) {
  if (index == my_paths.size()) {
    std::vector<std::vector<scalar_type>> value;

    value <= current_result[name]["data"];

    for (size_t i = 0; i < value.size(); i++)
      for (size_t j = 0; j < value[i].size(); j++)
        A(i, j) = value[i][j];
  }
  else {
    const JsonAccessor& new_result(current_result[my_paths[index]]);

    execute(name, A, new_result, (++index));
  }
}

template <typename scalar_type>
void reader<IO::JSON>::execute(std::string name,
                               LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU>& A,
                               const JsonAccessor& current_result, size_t index) {
  if (index == my_paths.size()) {
    std::vector<std::vector<scalar_type>> value;

    value <= current_result[name]["data-real"];

    for (size_t i = 0; i < value.size(); i++)
      for (size_t j = 0; j < value[i].size(); j++)
        real(A(i, j)) = value[i][j];

    value <= current_result[name]["data-imag"];

    for (size_t i = 0; i < value.size(); i++)
      for (size_t j = 0; j < value[i].size(); j++)
        imag(A(i, j)) = value[i][j];
  }
  else {
    const JsonAccessor& new_result(current_result[my_paths[index]]);

    execute(name, A, new_result, (++index));
  }
}

void reader<IO::JSON>::parse(std::string& file_name_ref) {
  std::string file_name = file_name_ref;

  std::wifstream file(file_name.c_str());

  if (!file or !file.good() or file.bad()) {
    std::cout << "\n\n\tcannot open file : " << file_name << "\n";
    throw std::runtime_error(__FUNCTION__);
  }
  else {
    std::cout << "\n\n\topening file : " << file_name << "\n";
  }

  while (parser.execute(file))
    ;
}

const reader<IO::JSON>::JsonAccessor& reader<IO::JSON>::operator[](const std::string& key) const {
  return parse_result[key];
}

const reader<IO::JSON>::JsonAccessor& reader<IO::JSON>::operator[](int index) {
  return parse_result[index];
}

const reader<IO::JSON>::JsonAccessor& reader<IO::JSON>::get_JSON_tree() {
  return parse_result;
}
}

#endif
