//-*-C++-*-

#ifndef CLUSTER_DOMAIN_FAMILY_H
#define CLUSTER_DOMAIN_FAMILY_H

/*!
 *  \author Peter Staar
 */
template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S>
class cluster_domain_family
{
public:
  
  const static int DIMENSION = D;

  const static CLUSTER_NAMES NAME  = N;
  const static CLUSTER_SHAPE SHAPE = S;

  typedef cluster_domain<scalar_type, D, N, REAL_SPACE    , S> r_cluster_type;
  typedef cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S> k_cluster_type;

public:

  template<IO::FORMAT DATA_FORMAT>
  static void read(IO::reader<DATA_FORMAT>& reader);

  template<IO::FORMAT DATA_FORMAT>
  static void write(IO::writer<DATA_FORMAT>& reader);

};

template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S>
template<IO::FORMAT DATA_FORMAT>
void cluster_domain_family<scalar_type, D, N, S>::write(IO::writer<DATA_FORMAT>& writer)
{
  writer.open_group(to_str(N));

  {
    writer.open_group(to_str(MOMENTUM_SPACE));
    
    writer.execute("basis"      , k_cluster_type::get_basis_vectors());
    writer.execute("super-basis", k_cluster_type::get_super_basis_vectors());
    writer.execute("elements"   , k_cluster_type::get_elements());

    writer.close_group();
  }

  {
    writer.open_group(to_str(REAL_SPACE));

    writer.execute("basis"      , r_cluster_type::get_basis_vectors());
    writer.execute("super-basis", r_cluster_type::get_super_basis_vectors());
    writer.execute("elements"   , r_cluster_type::get_elements());

    writer.close_group();
  }

  writer.close_group();
}



#endif
