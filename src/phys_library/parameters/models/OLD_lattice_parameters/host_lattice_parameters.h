//-*-C++-*-

#ifndef HOST_LATTICE_PARAMETERS_H
#define HOST_LATTICE_PARAMETERS_H

/*!
 * \author peter staar
 */
template<class model>
class host_lattice_parameters 
{
 public:

  typedef model                          model_type;
  typedef host_lattice_parameters<model > this_type;

  typedef typename model::DCA_point_group point_group;

  const static int                DIMENSION     = model::DIMENSION;
  const static cluster_shape_type cluster_shape = model::DCA_cluster_shape;

public:

  static double*& lattice_vectors();
  static double*& reciprocal_lattice_vectors();
};

/******************************************
 ***        STATIC METHODS              ***
 ******************************************/

template<class model >
double*& host_lattice_parameters<model >::lattice_vectors()
{
  static double* ptr = model::get_r_DCA_basis();
  return ptr;
}

template<class model >
double*& host_lattice_parameters<model >::reciprocal_lattice_vectors()
{
  static double* ptr = model::get_k_DCA_basis();
  return ptr;
}

#endif

