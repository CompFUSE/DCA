//-*-C++-*-

#ifndef MC_SOLVER_ANALYSIS_INTERPOLATION_PARAMETERS_H
#define MC_SOLVER_ANALYSIS_INTERPOLATION_PARAMETERS_H

#include "phys_library/parameters/parameters_specialization/templates/MC_solver_parameters.h"

/*!
 *   \author peter staar
 *   \brief  This class organizes the input-parameters for the ANALYSIS_INTERPOLATION-QMC method.
 */

template<>
class MC_solver_parameters<DCA::ANALYSIS_INTERPOLATION>
{
public:

  MC_solver_parameters();
  ~MC_solver_parameters();

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

  template<class concurrency_type>
  int  get_buffer_size(const concurrency_type& concurrency) const;

  template<class concurrency_type>
  void pack           (const concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  template<class concurrency_type>
  void unpack         (const concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

  template<class stream_type>
  void to_JSON(stream_type& ss, bool is_end=false);
  
  template<class JSON_reader_type>
  void from_JSON(JSON_reader_type& reader);

/******************************************
 ***        DATA                        ***
 ******************************************/

  std::vector<std::vector<int   > >& get_PCM_Bett_matrix();

private:

  template<typename scalartype>
  std::vector<scalartype> linearize(const std::vector<std::vector<scalartype> >& vec) const;

  template<typename scalartype>
  void un_linearize(int                                    dimension,
		    std::vector<scalartype>&		   vec_lin,
		    std::vector<std::vector<scalartype> >& vec);

private:

  int    dimension;

  std::vector<std::vector<int   > > Bett_matrix;
};

MC_solver_parameters<DCA::ANALYSIS_INTERPOLATION>::MC_solver_parameters():
  dimension(-1),
  Bett_matrix(0)
{}

MC_solver_parameters<DCA::ANALYSIS_INTERPOLATION>::~MC_solver_parameters()
{}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template<class concurrency_type>
int MC_solver_parameters<DCA::ANALYSIS_INTERPOLATION> ::get_buffer_size(const concurrency_type& concurrency) const
{
  int buffer_size = 0;

  buffer_size += concurrency.getBufferSize(dimension);

  {
    std::vector<int> tmp = linearize(Bett_matrix);
    buffer_size += concurrency.getBufferSize(tmp);
  }

  return buffer_size;
}

template<class concurrency_type>
void MC_solver_parameters<DCA::ANALYSIS_INTERPOLATION> ::pack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.pack(buffer, buffer_size, position, dimension);

  {
    std::vector<int> tmp = linearize(Bett_matrix);
    concurrency.pack(buffer, buffer_size, position, tmp);
  }
}

template<class concurrency_type>
void MC_solver_parameters<DCA::ANALYSIS_INTERPOLATION> ::unpack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position)
{
  concurrency.unpack(buffer, buffer_size, position, dimension);

  {
    std::vector<int> tmp;
    concurrency.unpack(buffer, buffer_size, position, tmp);
    un_linearize(dimension, tmp, Bett_matrix);
  }
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template<class stream_type>
void MC_solver_parameters<DCA::ANALYSIS_INTERPOLATION> ::to_JSON(stream_type& ss, bool is_end)
{
  ss << "\"MC-solver-parameters\" :";
  ss << "\n{ \n";

  JSON_writer::write(ss, "PCM-cluster", Bett_matrix, true);

  if(is_end)
    ss << "}\n";
  else
    ss << "},\n";  
}
  
template<class JSON_reader_type>
void MC_solver_parameters<DCA::ANALYSIS_INTERPOLATION> ::from_JSON(JSON_reader_type& reader)
{
  typedef typename JSON_reader_type::JsonAccessor JsonAccessor;
  const JsonAccessor control(reader["MC-solver-parameters"]);
     
  Bett_matrix <= control["PCM-cluster"];

  dimension = Bett_matrix.size();
}

/******************************************
 ***        DATA                        ***
 ******************************************/


std::vector<std::vector<int> >& MC_solver_parameters<DCA::ANALYSIS_INTERPOLATION>::get_PCM_Bett_matrix()
{
  return Bett_matrix;
}

/******************************************
 ***        PRIVATE                     ***
 ******************************************/

template<typename scalartype>
std::vector<scalartype> MC_solver_parameters<DCA::ANALYSIS_INTERPOLATION>::linearize(const std::vector<std::vector<scalartype> >& vec) const
{
  std::vector<scalartype> result;

  for(size_t l1=0; l1<vec.size(); l1++)
    for(size_t l2=0; l2<vec[l1].size(); l2++)
      result.push_back(vec[l1][l2]);

  return result;
}

template<typename scalartype>
void MC_solver_parameters<DCA::ANALYSIS_INTERPOLATION>::un_linearize(int                                    dimension,
					     std::vector<scalartype>&		    vec_lin,
					     std::vector<std::vector<scalartype> >& vec)
{
  vec.resize(vec_lin.size()/dimension, std::vector<scalartype>(0));

  for(size_t l1=0; l1<vec.size(); l1++)
    for(int   l2=0; l2<dimension; l2++)
      vec[l1].push_back(vec_lin[l2 + l1*dimension]);
}

#endif
