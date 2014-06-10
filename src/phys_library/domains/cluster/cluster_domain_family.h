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

//   template<class stream_type>
//   static void to_JSON(stream_type& ss);
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

/*
template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S>
template<class stream_type>
void cluster_domain_family<scalar_type, D, N, S>::to_JSON(stream_type& ss)
{
  ss << "\"cluster\" : ";

  ss << "{\n";

  ss << "\"SUPER-BASIS\" : \n {" ;
  ss << "\"R\" : {\n";
  for(int i=0; i<DIMENSION; i++){
    ss << "\"R_"<<i<<"\" : [";
    for(int j=0; j<DIMENSION; j++){
      ss << r_cluster_type::get_super_basis_vectors()[i][j];
      
      if(j==DIMENSION-1)
	ss << "]";
      else
	ss << ", ";
    }
    
    if(i==DIMENSION-1)
      ss << "\n";
    else
      ss << ",\n";
  }
  ss << "},\n"; 
  ss << "\"K\" : {\n";
  for(int i=0; i<DIMENSION; i++){
    ss << "\"K_"<<i<<"\" : [";
    for(int j=0; j<DIMENSION; j++){
      ss << k_cluster_type::get_super_basis_vectors()[i][j];
      
      if(j==DIMENSION-1)
	ss << "]";
      else
	ss << ", ";
    }
    
    if(i==DIMENSION-1)
      ss << "\n";
    else
      ss << ",\n";
  }
  ss << "}\n"; 
  ss << "},\n";

  ss << "\"BASIS\" : \n {" ;
  ss << "\"R\" : {\n";
  for(int i=0; i<DIMENSION; i++){
    ss << "\"R_"<<i<<"\" : [";
    for(int j=0; j<DIMENSION; j++){
      ss << r_cluster_type::get_basis_vectors()[i][j];
      
      if(j==DIMENSION-1)
	ss << "]";
      else
	ss << ", ";
    }
    
    if(i==DIMENSION-1)
      ss << "\n";
    else
      ss << ",\n";
  }
  ss << "},\n"; 
  ss << "\"K\" : {\n";
  for(int i=0; i<DIMENSION; i++){
    ss << "\"K_"<<i<<"\" : [";
    for(int j=0; j<DIMENSION; j++){
      ss << k_cluster_type::get_basis_vectors()[i][j];
      
      if(j==DIMENSION-1)
	ss << "]";
      else
	ss << ", ";
    }

    if(i==DIMENSION-1)
      ss << "\n";
    else
      ss << ",\n";
  }
  ss << "}\n"; 
  ss << "},\n";

  ss << "\"FULL-CLUSTER\" : \n {";
  ss << "\"R\" : [\n";
  for(int i=0; i<r_cluster_type::get_size(); i++){
    ss << "[";
    for(int j=0; j<DIMENSION; j++)
      if(j == DIMENSION-1)
	if(i == r_cluster_type::get_size()-1)
	  ss << r_cluster_type::get_elements()[i][j] << "]\n";
	else
	  ss << r_cluster_type::get_elements()[i][j] << "],\n";
       else
	ss << r_cluster_type::get_elements()[i][j] << ", ";
  }
  ss << "],\n"; 
  ss << "\"K\" : [\n";
  for(int i=0; i<int(k_cluster_type::get_size()); i++){
    ss << "[";
    for(int j=0; j<DIMENSION; j++)
      if(j == DIMENSION-1)
	if(i == int(k_cluster_type::get_size())-1)
	  ss << k_cluster_type::get_elements()[i][j] << "]\n";
	else
	  ss << k_cluster_type::get_elements()[i][j] << "],\n";
      else
	ss << k_cluster_type::get_elements()[i][j] << ", ";
  }
  ss << "]\n"; 
  ss << "}";

  ss << "}";
}
*/

#endif
