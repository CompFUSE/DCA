//-*-C++-*-

#ifndef DCA_PARAMETERS_H
#define DCA_PARAMETERS_H
#include"phys_library/domain_types.hpp"
using namespace types;

/*!
 *   \ingroup  PARAMETERS
 *
 *   \author   Peter Staar
 *   \brief    ...
 */
class DCA_Parameters
{

public:

  DCA_Parameters();
  ~DCA_Parameters();

  /******************************************
   ***        CONCURRENCY                 ***
   ******************************************/

  template<class concurrency_type>
  int  get_buffer_size( concurrency_type& concurrency);

  template<class concurrency_type>
  void pack           ( concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  template<class concurrency_type>
  void unpack         ( concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  /******************************************
   ***        READ/WRITE                  ***
   ******************************************/

  template<class stream_type>
  void to_JSON(stream_type& ss, bool is_end=false);

  template<class JSON_reader_type>
  void from_JSON(JSON_reader_type& reader);

  template<class read_write_type>
  void read_write(read_write_type& read_write_obj);

  /******************************************
   ***        DATA                        ***
   ******************************************/

  bool is_an_interacting_band(int i);

  std::vector<int>&  get_interacting_bands();

  int    get_DCA_iterations();
  double get_DCA_accuracy();
  double get_DCA_convergence_factor();

  std::vector<std::vector<int> > get_DCA_cluster();

  int get_k_mesh_refinement();
  int get_number_of_periods();

  bool print_phi_k();
  int  get_gaussian_quadrature_rule();

  int get_nr_coarsegraining_threads();
  int get_number_of_tail_frequencies();

  double get_integration_accuracy();
  bool   precompute_H_k();

  bool use_interpolated_Self_energy();
  bool use_wannier_interpolation();

  bool use_HTS_approximation();

  int    get_deconvolution_iterations();
  double get_deconvolution_tolerance();

private:

  template<typename reader_t>
  bool read_boolean(reader_t& control, std::string name);

  template<typename reader_t>
  void read_Bett_vectors(reader_t& control, std::vector<std::vector<int> >& matrix);

private:

  // general
  std::string                     do_DCA_plus;

  std::vector<int>                interacting_bands;
  int                             DCA_iterations;
  double                          DCA_accuracy;
  double                          DCA_mixing_factor;

  std::vector<std::vector<int> >  DCA_cluster;

  // cluster-mapping
  int                             k_mesh_refinement;
  int                             number_of_periods;

  int                             quadrature_rule;
  std::string                     precompute_Hamiltonian;

  int                             nr_coarsegraining_threads;
  int                             nr_tail_frequencies;

  double                          phi_k_integration_accuracy;
  std::string                     print_phi_k_to_python;

  // lattice-mapping
  std::string                     interpolation_method;
  std::string                     HTS_approximation;

  double                          deconvolution_tolerance;
  int                             max_deconvolution_iterations;
};

DCA_Parameters::DCA_Parameters():
  do_DCA_plus("false"),

  interacting_bands(0),
  DCA_iterations(1),
  DCA_accuracy(0.01),
  DCA_mixing_factor(1.),

  DCA_cluster(model::DIMENSION, std::vector<int>(model::DIMENSION,0)),

  // cluster-mapping
  k_mesh_refinement(3),
  number_of_periods(0),

  quadrature_rule(3),
  precompute_Hamiltonian("true"),

  nr_coarsegraining_threads(1),
  nr_tail_frequencies(0),

  phi_k_integration_accuracy(1.e-3),
  print_phi_k_to_python("false"),

  // lattice-mapping
  interpolation_method("wannier-interpolation"),
  HTS_approximation("false"),

  deconvolution_tolerance(0.01),
  max_deconvolution_iterations(16)
{}

DCA_Parameters::~DCA_Parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<class concurrency_type>
int DCA_Parameters::get_buffer_size( concurrency_type& concurrency)
{
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(do_DCA_plus);

  buffer_size += concurrency.get_buffer_size(interacting_bands);
  buffer_size += concurrency.get_buffer_size(DCA_iterations);
  buffer_size += concurrency.get_buffer_size(DCA_accuracy);
  buffer_size += concurrency.get_buffer_size(DCA_mixing_factor);

  buffer_size += concurrency.get_buffer_size(DCA_cluster);

  // cluster-mapping
  buffer_size += concurrency.get_buffer_size(k_mesh_refinement);
  buffer_size += concurrency.get_buffer_size(number_of_periods);

  buffer_size += concurrency.get_buffer_size(quadrature_rule);
  buffer_size += concurrency.get_buffer_size(precompute_Hamiltonian);

  buffer_size += concurrency.get_buffer_size(nr_coarsegraining_threads);
  buffer_size += concurrency.get_buffer_size(nr_tail_frequencies);

  buffer_size += concurrency.get_buffer_size(phi_k_integration_accuracy);
  buffer_size += concurrency.get_buffer_size(print_phi_k_to_python);

  // lattice-mapping
  buffer_size += concurrency.get_buffer_size(interpolation_method);
  buffer_size += concurrency.get_buffer_size(HTS_approximation);
  buffer_size += concurrency.get_buffer_size(deconvolution_tolerance);
  buffer_size += concurrency.get_buffer_size(max_deconvolution_iterations);

  return buffer_size;
}

template<class concurrency_type>
void DCA_Parameters::pack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.pack(buffer, buffer_size, position, do_DCA_plus);

  concurrency.pack(buffer, buffer_size, position, interacting_bands);
  concurrency.pack(buffer, buffer_size, position, DCA_iterations);
  concurrency.pack(buffer, buffer_size, position, DCA_accuracy);
  concurrency.pack(buffer, buffer_size, position, DCA_mixing_factor);

  concurrency.pack(buffer, buffer_size, position, DCA_cluster);

  // cluster-mapping
  concurrency.pack(buffer, buffer_size, position, k_mesh_refinement);
  concurrency.pack(buffer, buffer_size, position, number_of_periods);

  concurrency.pack(buffer, buffer_size, position, quadrature_rule);
  concurrency.pack(buffer, buffer_size, position, precompute_Hamiltonian);

  concurrency.pack(buffer, buffer_size, position, nr_coarsegraining_threads);
  concurrency.pack(buffer, buffer_size, position, nr_tail_frequencies);

  concurrency.pack(buffer, buffer_size, position, phi_k_integration_accuracy);
  concurrency.pack(buffer, buffer_size, position, print_phi_k_to_python);

  // lattice-mapping
  concurrency.pack(buffer, buffer_size, position, interpolation_method);
  concurrency.pack(buffer, buffer_size, position, HTS_approximation);
  concurrency.pack(buffer, buffer_size, position, deconvolution_tolerance);
  concurrency.pack(buffer, buffer_size, position, max_deconvolution_iterations);
}

template<class concurrency_type>
void DCA_Parameters::unpack( concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, do_DCA_plus);

  concurrency.unpack(buffer, buffer_size, position, interacting_bands);
  concurrency.unpack(buffer, buffer_size, position, DCA_iterations);
  concurrency.unpack(buffer, buffer_size, position, DCA_accuracy);
  concurrency.unpack(buffer, buffer_size, position, DCA_mixing_factor);

  concurrency.unpack(buffer, buffer_size, position, DCA_cluster);

  // cluster-mapping
  concurrency.unpack(buffer, buffer_size, position, k_mesh_refinement);
  concurrency.unpack(buffer, buffer_size, position, number_of_periods);

  concurrency.unpack(buffer, buffer_size, position, quadrature_rule);
  concurrency.unpack(buffer, buffer_size, position, precompute_Hamiltonian);

  concurrency.unpack(buffer, buffer_size, position, nr_coarsegraining_threads);
  concurrency.unpack(buffer, buffer_size, position, nr_tail_frequencies);

  concurrency.unpack(buffer, buffer_size, position, phi_k_integration_accuracy);
  concurrency.unpack(buffer, buffer_size, position, print_phi_k_to_python);

  // lattice-mapping
  concurrency.unpack(buffer, buffer_size, position, interpolation_method);
  concurrency.unpack(buffer, buffer_size, position, HTS_approximation);
  concurrency.unpack(buffer, buffer_size, position, deconvolution_tolerance);
  concurrency.unpack(buffer, buffer_size, position, max_deconvolution_iterations);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<class read_write_type>
void DCA_Parameters::read_write(read_write_type& read_write_obj)
{
  try
    {
      read_write_obj.open_group("DCA");

      try { read_write_obj.execute("do-DCA+"          , do_DCA_plus);       } catch(const std::exception& r_e) {}

      try { read_write_obj.execute("interacting-bands", interacting_bands); } catch(const std::exception& r_e) {}

      try { read_write_obj.execute("DCA-iterations"   , DCA_iterations);    } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("DCA-accuracy"     , DCA_accuracy);      } catch(const std::exception& r_e) {}
      try { read_write_obj.execute("DCA-mixing-factor", DCA_mixing_factor); } catch(const std::exception& r_e) {}

      try { read_write_obj.execute("cluster"          , DCA_cluster);       } catch(const std::exception& r_e) {}

      {
        read_write_obj.open_group("cluster-mapping");

        try { read_write_obj.execute("k-mesh-refinement"          , k_mesh_refinement);          } catch(const std::exception& r_e) {}
        try { read_write_obj.execute("quadrature-rule"            , quadrature_rule);            } catch(const std::exception& r_e) {}
        try { read_write_obj.execute("number-of-periods"          , number_of_periods);          } catch(const std::exception& r_e) {}

        try { read_write_obj.execute("number-of-threads"          , nr_coarsegraining_threads);  } catch(const std::exception& r_e) {}
	try { read_write_obj.execute("number-of-tail-frequencies" , nr_tail_frequencies);        } catch(const std::exception& r_e) {}

        try { read_write_obj.execute("precompute-Hamiltonian"     , precompute_Hamiltonian);     } catch(const std::exception& r_e) {}
        try { read_write_obj.execute("phi(k) integration accuracy", phi_k_integration_accuracy); } catch(const std::exception& r_e) {}
        try { read_write_obj.execute("print-phi(k)"               , print_phi_k_to_python);      } catch(const std::exception& r_e) {}

        read_write_obj.close_group();
      }

      {
        read_write_obj.open_group("lattice-mapping");

        try { read_write_obj.execute("interpolation-method"        , interpolation_method);         } catch(const std::exception& r_e) {}
        try { read_write_obj.execute("HTS-approximation"           , HTS_approximation);            } catch(const std::exception& r_e) {}
        try { read_write_obj.execute("deconvolution-tolerance"     , deconvolution_tolerance);      } catch(const std::exception& r_e) {}
        try { read_write_obj.execute("max-deconvolution-iterations", max_deconvolution_iterations); } catch(const std::exception& r_e) {}

        read_write_obj.close_group();
      }

      read_write_obj.close_group();
    }
  catch(const std::exception& r_e)
    {
      std::cout << "\n\t no DCA-parameters defined !!  \n\n";
      throw std::logic_error(__PRETTY_FUNCTION__);
    }
}

/******************************************
 ***        DATA                        ***
 ******************************************/

bool DCA_Parameters::use_interpolated_Self_energy()
{
  if(do_DCA_plus=="true")
    return true;
  else
    return false;
}

bool DCA_Parameters::is_an_interacting_band(int i)
{
  bool result = false;

  for(int l=0; l<interacting_bands.size(); l++)
    if(i==interacting_bands[l])
      result = true;

  return result;
}

std::vector<int>& DCA_Parameters::get_interacting_bands()
{
  return interacting_bands;
}

int DCA_Parameters::get_DCA_iterations()
{
  return DCA_iterations;
}

double DCA_Parameters::get_DCA_accuracy()
{
  return DCA_accuracy;
}

double DCA_Parameters::get_DCA_convergence_factor()
{
  return DCA_mixing_factor;
}

std::vector<std::vector<int> > DCA_Parameters::get_DCA_cluster()
{
  return DCA_cluster;
}

int DCA_Parameters::get_k_mesh_refinement()
{
  return k_mesh_refinement;
}

int DCA_Parameters::get_number_of_periods()
{
  return number_of_periods;
}

int DCA_Parameters::get_gaussian_quadrature_rule()
{
  return quadrature_rule;
}

int DCA_Parameters::get_nr_coarsegraining_threads()
{
  return nr_coarsegraining_threads;
}

int DCA_Parameters::get_number_of_tail_frequencies()
{
  return nr_tail_frequencies;
}

bool DCA_Parameters::print_phi_k()
{
  if(print_phi_k_to_python=="true")
    return true;
  else
    return false;
}

double DCA_Parameters::get_integration_accuracy()
{
  return phi_k_integration_accuracy;
}

bool DCA_Parameters::precompute_H_k()
{
  if(precompute_Hamiltonian=="true")
    return true;
  else
    return false;
}

bool DCA_Parameters::use_wannier_interpolation()
{
  return true;
}

bool DCA_Parameters::use_HTS_approximation()
{
  if(HTS_approximation=="true")
    return true;
  else
    return false;
}

int DCA_Parameters::get_deconvolution_iterations()
{
  return max_deconvolution_iterations;
}

double DCA_Parameters::get_deconvolution_tolerance()
{
  return deconvolution_tolerance;
}


#endif

