//-*-C++-*-

#ifndef CLUSTER_DOMAIN_SYMMETRY_H
#define CLUSTER_DOMAIN_SYMMETRY_H

/*!
 *  \author Peter Staar
 */
template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
class cluster_symmetry<cluster_domain<scalar_type, D, N, R, S> >
{
  const static int DIMENSION = D;

  const static CLUSTER_NAMES          NAME           = N;
  const static CLUSTER_REPRESENTATION REPRESENTATION = R;
  const static CLUSTER_SHAPE          SHAPE          = S;

  const static CLUSTER_REPRESENTATION DUAL_REPRESENTATION = dual_cluster<REPRESENTATION>::REPRESENTATION;

public:

  typedef cluster_domain_family<scalar_type, D, N, S>               cluster_family_type;

  typedef cluster_domain<scalar_type, D, N,      REPRESENTATION, S> this_type;
  typedef cluster_domain<scalar_type, D, N, DUAL_REPRESENTATION, S> dual_type;

  typedef dmn_0<electron_band_domain> b_dmn_t;
  typedef dmn_0<this_type>         c_dmn_t;

  typedef dmn_0<point_group_symmetry_domain<UNIT_CELL , cluster_family_type> > sym_unit_cell_dmn_t;
  typedef dmn_0<point_group_symmetry_domain<SUPER_CELL, cluster_family_type> > sym_super_cell_dmn_t;

  typedef dmn_2< dmn_2<c_dmn_t,b_dmn_t>, sym_super_cell_dmn_t> symmetry_matrix_dmn_t;

public:

  static FUNC_LIB::function<std::pair<int,int>, symmetry_matrix_dmn_t>& get_symmetry_matrix(){
    static FUNC_LIB::function<std::pair<int,int>, symmetry_matrix_dmn_t> symmetry_matrix("symmetry_matrix_super_cell");
    return symmetry_matrix;
  }

};


#endif
