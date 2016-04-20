//-*-C++-*-

#ifndef PARAMETERS_H
#define PARAMETERS_H
#include "phys_library/domain_types.hpp"
using namespace types;

/*!
 *   \ingroup  PARAMETERS
 *
 *   \author   Peter Staar, Raffaele Solca.
 *   \brief    ...
 */
template <class concurrency_t, class model_t, DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
class Parameters : public file_names_parameters,
                   public profiling_parameters<concurrency_t>,
                   public physics_parameters,
                   public model_parameters<model_t>,
                   public DCA_Parameters,

                   public function_parameters,
                   public MCI_parameters,

                   public MC_solver_parameters<CLUSTER_SOLVER_NAME>,
                   public ED_solver_parameters,

                   public vertex_parameters,
                   public equal_time_parameters,
                   public CPE_parameters,

                   public brillouin_zone_parameters,
                   public double_counting_parameters {
public:
  typedef concurrency_t concurrency_type;

#ifdef SINGLE_PRECISION_MEASUREMENTS
  typedef float MC_measurement_scalar_type;
#else
  typedef double MC_measurement_scalar_type;
#endif

#ifdef USE_REDUCED_VERTEX_FUNCTION
  typedef w_VERTEX_EXTENDED_POS G4_w1_dmn_t;
  typedef w_VERTEX_EXTENDED G4_w2_dmn_t;
#else
  typedef w_VERTEX_EXTENDED G4_w1_dmn_t;
  typedef w_VERTEX_EXTENDED G4_w2_dmn_t;
#endif

public:
  Parameters(std::string version_stamp, concurrency_type& concurrency_obj);
  ~Parameters();

  template <IO::FORMAT DATA_FORMAT>
  void read(IO::reader<DATA_FORMAT>& reader);

  template <IO::FORMAT DATA_FORMAT>
  void write(IO::writer<DATA_FORMAT>& writer);

  void read_input_and_broadcast(std::string file_name);

  template <typename read_write_t>
  void read_write(read_write_t& read_write_obj);

  void update_model();
  void update_domains();

  int get_buffer_size(concurrency_type& concurrency);
  void pack(concurrency_type& concurrency, int* buffer, int buffer_size, int& position);
  void unpack(concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  // get_functions
  concurrency_type& get_concurrency();

private:
  static std::string make_python_readable(std::string tmp);

private:
  std::string version_stamp;

  std::string date_str;
  std::string time_str;
  std::string compiler_str;

  concurrency_type& concurrency_obj;
};

template <class concurrency_t, class model_t, DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
Parameters<concurrency_t, model_t, CLUSTER_SOLVER_NAME>::Parameters(std::string version_stamp_str,
                                                                    concurrency_t& concurrency_object)
    : file_names_parameters(),
      physics_parameters(),
      model_parameters<model_t>(),
      DCA_Parameters(),

      function_parameters(),
      MCI_parameters(),

      MC_solver_parameters<CLUSTER_SOLVER_NAME>(),
      vertex_parameters(),
      equal_time_parameters(),
      CPE_parameters(),

      brillouin_zone_parameters(),

      version_stamp(make_python_readable(version_stamp_str)),

      date_str(__DATE__),
      time_str(__TIME__),

      compiler_str("????"),

      concurrency_obj(concurrency_object) {
#ifdef _CRAYC
  compiler_str = _RELEASE_STRING;
#endif

#ifdef __GNUC__
  compiler_str = __VERSION__;
#endif
}

template <class concurrency_t, class model_t, DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
Parameters<concurrency_t, model_t, CLUSTER_SOLVER_NAME>::~Parameters() {}

template <class concurrency_t, class model_t, DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
template <IO::FORMAT DATA_FORMAT>
void Parameters<concurrency_t, model_t, CLUSTER_SOLVER_NAME>::write(IO::writer<DATA_FORMAT>& writer) {
  {
    writer.open_group("parameters");

    this->read_write(writer);

    writer.close_group();
  }

  {
    writer.open_group("domains");

    DCA_cluster_family_type::write(writer);
    HOST_sp_cluster_family_type::write(writer);
    HOST_tp_cluster_family_type::write(writer);

    t::parameter_type::write(writer);
    w::parameter_type::write(writer);

    DCA::vertex_frequency_domain<DCA::EXTENDED_BOSONIC>::write(writer);

    DCA_iteration_domain::write(writer);

    if (get_vertex_measurement_type() != NONE) {
      DCA::vertex_time_domain<DCA::SP_TIME_DOMAIN>::write(writer);
      DCA::vertex_time_domain<DCA::TP_TIME_DOMAIN>::write(writer);
      DCA::vertex_time_domain<DCA::SP_TIME_DOMAIN_POSITIVE>::write(writer);
      DCA::vertex_time_domain<DCA::TP_TIME_DOMAIN_POSITIVE>::write(writer);

      DCA::vertex_frequency_domain<DCA::COMPACT>::write(writer);
      DCA::vertex_frequency_domain<DCA::EXTENDED>::write(writer);

      DCA::vertex_frequency_domain<DCA::COMPACT_POSITIVE>::write(writer);
      DCA::vertex_frequency_domain<DCA::EXTENDED_POSITIVE>::write(writer);
    }

    if (do_CPE()) {
      frequency_domain_real_axis::write(writer);
      frequency_domain_imag_axis::write(writer);
    }

    if (QMC_INTEGRATOR_BIT)
      numerical_error_domain::write(writer);

    writer.close_group();
  }
}

template <class concurrency_t, class model_t, DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
void Parameters<concurrency_t, model_t, CLUSTER_SOLVER_NAME>::read_input_and_broadcast(
    std::string filename) {
  for (bool flag = false; !flag;) {
    if (concurrency_obj.id() == concurrency_obj.first()) {
      file_names_parameters::get_input_file_name() = filename;

      {
        IO::reader<IO::JSON> read_obj;

        read_obj.open_file(filename);

        this->read_write(read_obj);

        read_obj.close_file();
      }
    }

    flag = concurrency_obj.broadcast_object(*this);
  }
}

template <class concurrency_t, class model_t, DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
void Parameters<concurrency_t, model_t, CLUSTER_SOLVER_NAME>::update_model() {
  model::initialize(*this);
}

template <class concurrency_t, class model_t, DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
void Parameters<concurrency_t, model_t, CLUSTER_SOLVER_NAME>::update_domains() {
  DCA_iteration_domain::initialize(*this);
  electron_band_domain::initialize(*this, model::BANDS, model::get_flavors(), model::get_a_vectors());

  // time && frequency-domains
  time_domain::initialize(*this);
  time_domain_left_oriented::initialize(*this);
  frequency_domain::initialize(*this);
  frequency_domain_real_axis::initialize(*this);
  frequency_domain_imag_axis::initialize(*this);

  DCA::vertex_time_domain<DCA::SP_TIME_DOMAIN>::initialize(*this);
  DCA::vertex_time_domain<DCA::TP_TIME_DOMAIN>::initialize(*this);
  DCA::vertex_time_domain<DCA::SP_TIME_DOMAIN_POSITIVE>::initialize(*this);
  DCA::vertex_time_domain<DCA::TP_TIME_DOMAIN_POSITIVE>::initialize(*this);

  DCA::vertex_frequency_domain<DCA::COMPACT>::initialize(*this);
  DCA::vertex_frequency_domain<DCA::EXTENDED>::initialize(*this);

  DCA::vertex_frequency_domain<DCA::COMPACT_POSITIVE>::initialize(*this);
  DCA::vertex_frequency_domain<DCA::EXTENDED_POSITIVE>::initialize(*this);

  DCA::vertex_frequency_domain<DCA::EXTENDED_BOSONIC>::initialize(*this);

  // DCA
  cluster_domain_initializer<r_DCA>::execute(model::get_r_DCA_basis(),
                                             DCA_Parameters::get_DCA_cluster());

  cluster_domain_symmetry_initializer<r_DCA, DCA_point_group_type>::execute();

  if (concurrency_obj.id() == concurrency_obj.last())
    k_DCA::parameter_type::print(std::cout);

  // host
  cluster_domain_initializer<r_HOST>::execute(
      model::get_r_DCA_basis(),  // DCA_lattice_parameters_type::lattice_vectors(),
      function_parameters::get_sp_cluster());

  cluster_domain_symmetry_initializer<r_HOST, DCA_point_group_type>::execute();

  if (concurrency_obj.id() == concurrency_obj.last())
    k_HOST::parameter_type::print(std::cout);

  // host
  cluster_domain_initializer<r_HOST_VERTEX>::execute(
      model::get_r_DCA_basis(),  // DCA_lattice_parameters_type::lattice_vectors(),
      function_parameters::get_tp_cluster());

  cluster_domain_symmetry_initializer<r_HOST_VERTEX, DCA_point_group_type>::execute();

  if (concurrency_obj.id() == concurrency_obj.last())
    k_HOST_VERTEX::parameter_type::print(std::cout);

  // LDA
  cluster_domain_initializer<r_LDA>::execute(
      model::get_r_DCA_basis(),  // DCA_lattice_parameters_type::lattice_vectors(),
      function_parameters::get_H_k_grid_size());

  if (concurrency_obj.id() == concurrency_obj.last())
    k_LDA::parameter_type::print(std::cout);
}

template <class concurrency_t, class model_t, DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
template <typename read_write_t>
void Parameters<concurrency_t, model_t, CLUSTER_SOLVER_NAME>::read_write(read_write_t& read_write_obj) {
  if (read_write_obj.is_writer()) {
    read_write_obj.execute("date", date_str);
    read_write_obj.execute("time", time_str);
    read_write_obj.execute("compiler", compiler_str);
  }

  file_names_parameters::read_write(read_write_obj);

  physics_parameters::read_write(read_write_obj);
  model_parameters<model_t>::read_write(read_write_obj);
  DCA_Parameters::read_write(read_write_obj);

  MCI_parameters::read_write(read_write_obj);
  MC_solver_parameters<CLUSTER_SOLVER_NAME>::read_write(read_write_obj);
  ED_solver_parameters::read_write(read_write_obj);

  function_parameters::read_write(read_write_obj);
  vertex_parameters::read_write(read_write_obj);
  equal_time_parameters::read_write(read_write_obj);

  CPE_parameters::read_write(read_write_obj);
  double_counting_parameters::read_write(read_write_obj);
  brillouin_zone_parameters::read_write(read_write_obj);
}

template <class concurrency_t, class model_t, DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
int Parameters<concurrency_t, model_t, CLUSTER_SOLVER_NAME>::get_buffer_size(concurrency_t& concurrency) {
  int buffer_size = 0;

  buffer_size += file_names_parameters::get_buffer_size(concurrency);

  buffer_size += physics_parameters::get_buffer_size(concurrency);
  buffer_size += model_parameters<model_t>::get_buffer_size(concurrency);
  buffer_size += DCA_Parameters::get_buffer_size(concurrency);

  buffer_size += MCI_parameters::get_buffer_size(concurrency);
  buffer_size += MC_solver_parameters<CLUSTER_SOLVER_NAME>::get_buffer_size(concurrency);
  buffer_size += ED_solver_parameters::get_buffer_size(concurrency);

  buffer_size += function_parameters::get_buffer_size(concurrency);
  buffer_size += vertex_parameters::get_buffer_size(concurrency);
  buffer_size += equal_time_parameters::get_buffer_size(concurrency);
  buffer_size += CPE_parameters::get_buffer_size(concurrency);

  buffer_size += brillouin_zone_parameters::get_buffer_size(concurrency);
  buffer_size += double_counting_parameters::get_buffer_size(concurrency);
  buffer_size += brillouin_zone_parameters::get_buffer_size(concurrency);

  return buffer_size;
}

template <class concurrency_t, class model_t, DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
void Parameters<concurrency_t, model_t, CLUSTER_SOLVER_NAME>::pack(concurrency_t& concurrency,
                                                                   int* buffer, int buffer_size,
                                                                   int& position) {
  file_names_parameters::pack(concurrency, buffer, buffer_size, position);

  physics_parameters::pack(concurrency, buffer, buffer_size, position);
  model_parameters<model_t>::pack(concurrency, buffer, buffer_size, position);
  DCA_Parameters::pack(concurrency, buffer, buffer_size, position);

  MCI_parameters::pack(concurrency, buffer, buffer_size, position);
  MC_solver_parameters<CLUSTER_SOLVER_NAME>::pack(concurrency, buffer, buffer_size, position);
  ED_solver_parameters::pack(concurrency, buffer, buffer_size, position);

  function_parameters::pack(concurrency, buffer, buffer_size, position);
  vertex_parameters::pack(concurrency, buffer, buffer_size, position);
  equal_time_parameters::pack(concurrency, buffer, buffer_size, position);
  CPE_parameters::pack(concurrency, buffer, buffer_size, position);

  brillouin_zone_parameters::pack(concurrency, buffer, buffer_size, position);
  double_counting_parameters::pack(concurrency, buffer, buffer_size, position);
  brillouin_zone_parameters::pack(concurrency, buffer, buffer_size, position);
}

template <class concurrency_t, class model_t, DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
void Parameters<concurrency_t, model_t, CLUSTER_SOLVER_NAME>::unpack(concurrency_t& concurrency,
                                                                     int* buffer, int buffer_size,
                                                                     int& position) {
  file_names_parameters::unpack(concurrency, buffer, buffer_size, position);

  physics_parameters::unpack(concurrency, buffer, buffer_size, position);
  model_parameters<model_t>::unpack(concurrency, buffer, buffer_size, position);
  DCA_Parameters::unpack(concurrency, buffer, buffer_size, position);

  MCI_parameters::unpack(concurrency, buffer, buffer_size, position);
  MC_solver_parameters<CLUSTER_SOLVER_NAME>::unpack(concurrency, buffer, buffer_size, position);
  ED_solver_parameters::unpack(concurrency, buffer, buffer_size, position);

  function_parameters::unpack(concurrency, buffer, buffer_size, position);
  vertex_parameters::unpack(concurrency, buffer, buffer_size, position);
  equal_time_parameters::unpack(concurrency, buffer, buffer_size, position);
  CPE_parameters::unpack(concurrency, buffer, buffer_size, position);

  brillouin_zone_parameters::unpack(concurrency, buffer, buffer_size, position);
  double_counting_parameters::unpack(concurrency, buffer, buffer_size, position);
  brillouin_zone_parameters::unpack(concurrency, buffer, buffer_size, position);
}

template <class concurrency_t, class model_t, DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
concurrency_t& Parameters<concurrency_t, model_t, CLUSTER_SOLVER_NAME>::get_concurrency() {
  return concurrency_obj;
}

template <class concurrency_t, class model_t, DCA::CLUSTER_SOLVER_NAMES CLUSTER_SOLVER_NAME>
std::string Parameters<concurrency_t, model_t, CLUSTER_SOLVER_NAME>::make_python_readable(
    std::string str) {
  {
    std::string tmp("\n\n");
    while (true) {
      int pos = str.find(tmp);
      if (pos == -1)
        break;
      str.replace(pos, tmp.size(), std::string(", "));
    }
  }

  {
    std::string tmp("\n");
    while (true) {
      int pos = str.find(tmp);
      if (pos == -1)
        break;
      str.replace(pos, tmp.size(), std::string(", "));
    }
  }

  {
    std::string tmp("\t");
    while (true) {
      int pos = str.find(tmp);
      if (pos == -1)
        break;
      str.replace(pos, tmp.size(), std::string(" "));
    }
  }

  return str;
}

#endif
