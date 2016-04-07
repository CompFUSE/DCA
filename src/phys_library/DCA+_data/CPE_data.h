//-*-C++-*-

#ifndef MOMS_W_REAL_H
#define MOMS_W_REAL_H
#include"phys_library/domain_types.hpp"
using namespace types;

/*!
 * \class   MultiOrbitalMultiSiteStructure
 *
 * \brief   class that contains all functions on the real axis
 *
 * \date    6 Dec, 2013
 * \author  Peter Staar
 * \version 1.0
 */
template<class parameters_type>
class MOMS_w_real
{

  typedef typename parameters_type::profiler_type    profiler_t;
  typedef typename parameters_type::concurrency_type concurrency_type;

public:

  MOMS_w_real(parameters_type& parameters_ref);
  ~MOMS_w_real();

  void read(std::string filename);

  void write(std::string filename);

  template<IO::FORMAT DATA_FORMAT>
  void read(IO::reader<DATA_FORMAT>& reader);

  template<IO::FORMAT DATA_FORMAT>
  void write(IO::writer<DATA_FORMAT>& writer);

private:

  parameters_type&  parameters;
  concurrency_type& concurrency;

public:

  FUNC_LIB::function<double, w_REAL> A_w;
  FUNC_LIB::function<double, w_REAL> A_w_stddev;

  FUNC_LIB::function<double, dmn_3<b, s, w_REAL> > A_nu_w;
  FUNC_LIB::function<double, dmn_3<b, s, w_REAL> > A_nu_w_stddev;

  FUNC_LIB::function<double,             w_REAL  > A0_w;
  FUNC_LIB::function<double, dmn_3<b, s, w_REAL> > A0_nu_w;

  FUNC_LIB::function<double, w_REAL> E_w;
  FUNC_LIB::function<double, w_REAL> E0_w;

  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, w_REAL> > Sigma;
  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, w_REAL> > Sigma_stddev;

  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, w_REAL> > G0_k_w;
  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_DCA, w_REAL> > G0_r_w;

  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, w_REAL> > G_k_w;
  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, w_REAL> > G_k_w_stddev;

  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_DCA, w_REAL> > G_r_w;
};

template<class parameters_type>
MOMS_w_real<parameters_type>::MOMS_w_real(parameters_type& parameters_ref):

  parameters(parameters_ref),
  concurrency(parameters.get_concurrency()),

  A_w       ("spectral-density"),
  A_w_stddev("spectral-density-stddev"),

  A_nu_w       ("spectral-density-per-orbital"),
  A_nu_w_stddev("spectral-density-per-orbital-stddev"),

  A0_w   ("free-spectral-density"),
  A0_nu_w("free-spectral-density-per-orbital"),

  E_w ("E_w"),
  E0_w("E0_w"),

  Sigma       ("self-energy-real-axis"),
  Sigma_stddev("self-energy-real-axis-stddev"),

  G0_k_w("cluster-greens-function-G0_k_w-real-axis"),
  G0_r_w("cluster-greens-function-G0_r_w-real-axis"),

  G_k_w       ("cluster-greens-function-G_k_w-real-axis"),
  G_k_w_stddev("cluster-greens-function-G_k_w-real-axis-stddev"),

  G_r_w ("cluster-greens-function-G_r_w-real-axis")
{}

template<class parameters_type>
MOMS_w_real<parameters_type>::~MOMS_w_real()
{}

template<class parameters_type>
void MOMS_w_real<parameters_type>::read(std::string filename)
{
  if(concurrency.id()==concurrency.first())
    {
      std::cout << "\n\n\t starts reading \n\n";

      IO::FORMAT FORMAT = parameters.get_output_format();

      switch(FORMAT)
        {
        case IO::JSON :
          {
            IO::reader<IO::JSON> reader;

            reader.open_file(filename);

            this->read(reader);

            reader.close_file();
          }
          break;

        case IO::HDF5 :
          {
            IO::reader<IO::HDF5> reader;

            reader.open_file(filename);

            this->read(reader);

            reader.close_file();
          }
          break;

        default:
          throw std::logic_error(__FUNCTION__);
        }
    }

  concurrency.broadcast_object(A_w);
  concurrency.broadcast_object(A0_w);

  concurrency.broadcast_object(E_w);
  concurrency.broadcast_object(E0_w);

  concurrency.broadcast_object(Sigma);

  concurrency.broadcast_object(G_k_w);
  concurrency.broadcast_object(G0_k_w);

  concurrency.broadcast_object(G_r_w);
  concurrency.broadcast_object(G0_r_w);
}

template<class parameters_type>
template<IO::FORMAT DATA_FORMAT>
void MOMS_w_real<parameters_type>::read(IO::reader<DATA_FORMAT>& reader)
{
  reader.open_group("spectral-functions");

  reader.execute(A_w);   
  reader.execute(A0_w);  

  reader.execute(E_w);   
  reader.execute(E0_w);  

  reader.execute(Sigma); 

  reader.execute(G0_k_w);
  reader.execute(G_k_w); 

  reader.execute(G_r_w); 
  reader.execute(G0_r_w);

  reader.close_group();
}

template<class parameters_type>
void MOMS_w_real<parameters_type>::write(std::string filename)
{
  if(concurrency.id()==concurrency.first())
    {
      std::cout << "\n\n\t starts writing \n\n";

      IO::FORMAT FORMAT = parameters.get_output_format();

      switch(FORMAT)
        {
        case IO::JSON :
          {
            IO::writer<IO::JSON> writer;

            writer.open_file(filename);

            this->write(writer);

            writer.close_file();
          }
          break;

        case IO::HDF5 :
          {
            IO::writer<IO::HDF5> writer;

            writer.open_file(filename);

            this->write(writer);

            writer.close_file();
          }
          break;

        default:
          throw std::logic_error(__FUNCTION__);
        }
    }
}

template<class parameters_type>
template<IO::FORMAT DATA_FORMAT>
void MOMS_w_real<parameters_type>::write(IO::writer<DATA_FORMAT>& writer)
{
  writer.open_group("spectral-functions");

  writer.execute(A_w);
  writer.execute(A_w_stddev);

  writer.execute(A_nu_w);
  writer.execute(A_nu_w_stddev);

  writer.execute(A0_w);
  writer.execute(A0_nu_w);

  writer.execute(E_w);
  writer.execute(E0_w);
  
  writer.execute(Sigma);
  writer.execute(Sigma_stddev);

  writer.execute(G_k_w);
  writer.execute(G_k_w_stddev);

  writer.execute(G0_k_w);

  writer.execute(G_r_w);
  writer.execute(G0_r_w);

  writer.close_group();
}

/*
  template<class parameters_type>
  template<class stream_type>
  void MOMS_w_real<parameters_type>::to_JSON(stream_type& ss)
  {
  A_w .to_JSON(ss);
  ss << ",\n";

  A0_w.to_JSON(ss);
  ss << ",\n";

  E_w .to_JSON(ss);
  ss << ",\n";

  E0_w.to_JSON(ss);
  ss << ",\n";

  Sigma.to_JSON(ss);
  ss << ",\n";

  G0_k_w.to_JSON(ss);
  ss << ",\n";

  G_k_w.to_JSON(ss);
  ss << ",\n";

  G_r_w.to_JSON(ss);
  ss << ",\n";

  G0_r_w.to_JSON(ss);
  ss << ",\n";
  }
*/

#endif
