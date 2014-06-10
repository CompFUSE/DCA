//-*-C++-*-

#ifndef LDA_LATTICE_PARAMETERS_H
#define LDA_LATTICE_PARAMETERS_H

/*!
 *  \author peter staar
 */
template<class model>
class LDA_lattice_parameters {

 public:

  typedef model                          model_type;
  typedef LDA_lattice_parameters<model > this_type;

  typedef typename model::LDA_point_group point_group;
  //typedef typename model::DCA_point_group point_group;

  const static int                DIMENSION     = model::DIMENSION;
  const static cluster_shape_type cluster_shape = model::LDA_cluster_shape;

public:

  static double*&          lattice_vectors();
  static double*&          reciprocal_lattice_vectors();
};


template<class model >
double*& LDA_lattice_parameters<model >::lattice_vectors()
{
  static double* ptr = model::get_r_LDA_basis();
  return ptr;
}

template<class model >
double*& LDA_lattice_parameters<model >::reciprocal_lattice_vectors()
{
  static double* ptr = model::get_k_LDA_basis();
  return ptr;
}

#endif
