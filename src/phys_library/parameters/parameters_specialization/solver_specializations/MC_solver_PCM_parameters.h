//-*-C++-*-
// INTERNAL UNUSED FILE
#ifndef MC_SOLVER_PCM_PARAMETERS_H
#define MC_SOLVER_PCM_PARAMETERS_H

/*!
 *   \author peter staar
 *   \brief  This class organizes the input-parameters for the PCM-QMC method.
 */
template <>
class MC_solver_parameters<PCM> {
public:
  MC_solver_parameters();
  ~MC_solver_parameters();

  /******************************************
   ***        CONCURRENCY                 ***
   ******************************************/

  template <class concurrency_type>
  int get_buffer_size(const concurrency_type& concurrency) const;

  template <class concurrency_type>
  void pack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  template <class concurrency_type>
  void unpack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  /******************************************
   ***        READ/WRITE                  ***
   ******************************************/

  template <class stream_type>
  void to_JSON(stream_type& ss, bool is_end = false);

  template <class JSON_reader_type>
  void from_JSON(JSON_reader_type& reader);

  /******************************************
   ***        DATA                        ***
   ******************************************/

  int get_K_PHANI();
  double get_K_CT_AUX();

  std::vector<std::vector<double>>& get_phases();
  std::vector<std::vector<int>>& get_PCM_Bett_matrix();

private:
  template <typename scalartype>
  std::vector<scalartype> linearize(const std::vector<std::vector<scalartype>>& vec) const;

  template <typename scalartype>
  void un_linearize(int dimension, std::vector<scalartype>& vec_lin,
                    std::vector<std::vector<scalartype>>& vec);

private:
  int K_PHANI;
  double K_CT_AUX;

  int dimension;

  std::vector<std::vector<int>> Bett_matrix;
  std::vector<std::vector<double>> phases;
};

MC_solver_parameters<PCM>::MC_solver_parameters()
    : K_PHANI(64),
      K_CT_AUX(1),

      dimension(-1),

      Bett_matrix(0),
      phases(0) {}

MC_solver_parameters<PCM>::~MC_solver_parameters() {}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template <class concurrency_type>
int MC_solver_parameters<PCM>::get_buffer_size(const concurrency_type& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.getBufferSize(K_CT_AUX);
  buffer_size += concurrency.getBufferSize(K_PHANI);

  buffer_size += concurrency.getBufferSize(dimension);

  {
    std::vector<int> tmp = linearize(Bett_matrix);
    buffer_size += concurrency.getBufferSize(tmp);
  }

  {
    std::vector<double> tmp = linearize(phases);
    buffer_size += concurrency.getBufferSize(tmp);
  }

  return buffer_size;
}

template <class concurrency_type>
void MC_solver_parameters<PCM>::pack(const concurrency_type& concurrency, int* buffer,
                                     int buffer_size, int& position) {
  concurrency.pack(buffer, buffer_size, position, K_CT_AUX);
  concurrency.pack(buffer, buffer_size, position, K_PHANI);

  concurrency.pack(buffer, buffer_size, position, dimension);

  {
    std::vector<int> tmp = linearize(Bett_matrix);
    concurrency.pack(buffer, buffer_size, position, tmp);
  }

  {
    std::vector<double> tmp = linearize(phases);
    concurrency.pack(buffer, buffer_size, position, tmp);
  }
}

template <class concurrency_type>
void MC_solver_parameters<PCM>::unpack(const concurrency_type& concurrency, int* buffer,
                                       int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, K_CT_AUX);
  concurrency.unpack(buffer, buffer_size, position, K_PHANI);

  concurrency.unpack(buffer, buffer_size, position, dimension);

  {
    std::vector<int> tmp;
    concurrency.unpack(buffer, buffer_size, position, tmp);
    un_linearize(dimension, tmp, Bett_matrix);
  }

  {
    std::vector<double> tmp;
    concurrency.unpack(buffer, buffer_size, position, tmp);
    un_linearize(dimension, tmp, phases);
  }
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template <class stream_type>
void MC_solver_parameters<PCM>::to_JSON(stream_type& ss, bool is_end) {
  //   ss << "\"MC-solver-parameters\" :";
  //   ss << "\n{ \n";

  //   JSON_writer::write(ss, "K_PHANI"    , K_PHANI);
  //   JSON_writer::write(ss, "K_CT_AUX"   , K_CT_AUX);

  //   JSON_writer::write(ss, "PCM-cluster", Bett_matrix);

  //   JSON_writer::write(ss, "phases"     , phases,true);

  //   if(is_end)
  //     ss << "}\n";
  //   else
  //     ss << "},\n";
}

template <class JSON_reader_type>
void MC_solver_parameters<PCM>::from_JSON(JSON_reader_type& reader) {
  typedef typename JSON_reader_type::JsonAccessor JsonAccessor;
  const JsonAccessor control(reader["MC-solver-parameters"]);

  K_PHANI <= control["K_PHANI"];
  K_CT_AUX <= control["K_CT_AUX"];

  Bett_matrix <= control["PCM-cluster"];

  dimension = Bett_matrix.size();

  try {
    const JsonAccessor control2(reader["MC-solver-parameters"]);
    phases <= control2["phases"];
  }
  catch (const std::exception& r_e) {
  }
}

/******************************************
 ***        DATA                        ***
 ******************************************/

int MC_solver_parameters<PCM>::get_K_PHANI() {
  return K_PHANI;
}

double MC_solver_parameters<PCM>::get_K_CT_AUX() {
  return K_CT_AUX;
}

std::vector<std::vector<double>>& MC_solver_parameters<PCM>::get_phases() {
  return phases;
}

std::vector<std::vector<int>>& MC_solver_parameters<PCM>::get_PCM_Bett_matrix() {
  return Bett_matrix;
}

/******************************************
 ***        PRIVATE                     ***
 ******************************************/

template <typename scalartype>
std::vector<scalartype> MC_solver_parameters<PCM>::linearize(
    const std::vector<std::vector<scalartype>>& vec) const {
  std::vector<scalartype> result;

  for (size_t l1 = 0; l1 < vec.size(); l1++)
    for (size_t l2 = 0; l2 < vec[l1].size(); l2++)
      result.push_back(vec[l1][l2]);

  return result;
}

template <typename scalartype>
void MC_solver_parameters<PCM>::un_linearize(int dimension, std::vector<scalartype>& vec_lin,
                                             std::vector<std::vector<scalartype>>& vec) {
  vec.resize(vec_lin.size() / dimension, std::vector<scalartype>(0));

  for (size_t l1 = 0; l1 < vec.size(); l1++)
    for (int l2 = 0; l2 < dimension; l2++)
      vec[l1].push_back(vec_lin[l2 + l1 * dimension]);
}

#endif
