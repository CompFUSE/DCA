//-*-C++-*-

#ifndef CLUSTER_INTERFACE_H_
#define CLUSTER_INTERFACE_H_

/*!
 * \author Peter Staar
 */
template<class parameters >
class cluster
{

public:

  const static  int                        DIMENSION = parameters::DIMENSION;
  typedef parameters                       parameter_type;
  typedef typename parameters::point_group point_group_type;

  const static cluster_shape_type cluster_shape = parameters::cluster_shape;
                              
  typedef cluster<parameters>                                                        this_cluster;
  typedef cluster_initializer<parameters::DIMENSION, this_cluster, point_group_type> initializer;

  typedef point_group_symmetry_domain<UNIT_CELL , this_cluster> sym_unit_cell_t;
  typedef point_group_symmetry_domain<SUPER_CELL, this_cluster> sym_super_cell_t;

  typedef r_cluster<FULL, this_cluster> r_cluster_type;
  typedef k_cluster<FULL, this_cluster> k_cluster_type;

  typedef dmn_0<r_cluster_type> r_dmn_t;
  typedef dmn_0<k_cluster_type> k_dmn_t;

  typedef dmn_0<electron_band_domain> b_dmn_t;

  typedef dmn_0<sym_unit_cell_t>  sym_unit_cell_dmn_t;
  typedef dmn_0<sym_super_cell_t> sym_super_cell_dmn_t;


public:

  static int               get_dimension();
  static point_group_type& get_point_group();

  static std::vector<int>&                   grid_size(); 
  static std::vector<std::vector<int   > >&  get_Bett_matrix();

//   static std::vector<std::vector<double> >&  get_k_basis();
//   static std::vector<std::vector<double> >&  get_k_super_basis();

  static std::vector<std::vector<double> >&  get_k_basis_new();
  static std::vector<std::vector<double> >&  get_k_super_basis_new();

  static std::vector<std::vector<double> >&  get_r_basis();
  static std::vector<std::vector<double> >&  get_r_super_basis();

  static std::vector<std::vector<double> >&  get_r_cluster();
  static std::vector<std::vector<double> >&  get_k_cluster();
  
  /*
  static std::vector<std::vector<double> >&  get_affine_cluster();
  static std::vector<std::vector<int   > >&  get_integer_cluster();

  static std::vector<std::vector<double> >&  get_irreducible_r_cluster();
  static std::vector<std::vector<double> >&  get_irreducible_k_cluster();
  */

  static int get_cluster_size();

  static double& get_r_volume();
  static double& get_k_volume();

  static void print();

  template<class stream_type>
  static void to_JSON(stream_type& ss);

private:

  template<class stream_type>
  static void print_2D(stream_type& ss);

  template<class stream_type>
  static void print_3D(stream_type& ss);

  template<class stream_type>
  static void to_JSON_2D(stream_type& ss);

  template<class stream_type>
  static void to_JSON_3D(stream_type& ss);
};

template<class parameters>
int cluster<parameters>::get_dimension()
{
  const static int dim = parameters::DIMENSION; 
  return dim;
}

template<class parameters>
typename cluster<parameters>::point_group_type& cluster<parameters>::get_point_group()
{
  const static point_group_type symm_list;
  return symm_list;
}

template<class parameters>
std::vector<int>& cluster<parameters>::grid_size()
{
  if(parameters::cluster_shape != PARALLELEPIPED)
    throw std::logic_error(__FUNCTION__);

  static std::vector<int> v(0,DIMENSION);
  return v;
}

template<class parameters>
std::vector<std::vector<int> >& cluster<parameters>::get_Bett_matrix()
{
  if(parameters::cluster_shape != BETT_CLUSTER)
    throw std::logic_error(__FUNCTION__);

  static std::vector<std::vector<int> > Bett(DIMENSION,std::vector<int>(DIMENSION,0.));
  return Bett;
}

template<class parameters>
std::vector<std::vector<double> >& cluster<parameters>::get_k_basis_new()
{
  static std::vector<std::vector<double> > k(DIMENSION,std::vector<double>(DIMENSION,0.));
  return k;
}

template<class parameters>
std::vector<std::vector<double> >& cluster<parameters>::get_k_super_basis_new()
{
  static std::vector<std::vector<double> > k(DIMENSION,std::vector<double>(DIMENSION,0.));
  return k;
}

template<class parameters>
std::vector<std::vector<double> >& cluster<parameters>::get_r_basis()
{
  static std::vector<std::vector<double> > r(DIMENSION,std::vector<double>(DIMENSION,0.));
  return r;
}

template<class parameters>
std::vector<std::vector<double> >& cluster<parameters>::get_r_super_basis()
{
  static std::vector<std::vector<double> > r(DIMENSION,std::vector<double>(DIMENSION,0.));
  return r;
}





template<class parameters>
std::vector<std::vector<double> >& cluster<parameters>::get_r_cluster()
{
  static std::vector<std::vector<double> > cluster(0,std::vector<double>(parameters::DIMENSION,0.));
  return cluster;
}

// template<class parameters>
// std::vector<std::vector<double> >& cluster<parameters>::get_irreducible_r_cluster()
// {
//   static std::vector<std::vector<double> > ircluster(0,std::vector<double>(parameters::DIMENSION,0.));
//   return ircluster;
// }

template<class parameters>
std::vector<std::vector<double> >& cluster<parameters>::get_k_cluster()
{
  static std::vector<std::vector<double> > cluster(0,std::vector<double>(parameters::DIMENSION,0.));
  return cluster;
}

// template<class parameters>
// std::vector<std::vector<double> >& cluster<parameters>::get_irreducible_k_cluster()
// {
//   static std::vector<std::vector<double> > ircluster(0,std::vector<double>(parameters::DIMENSION,0.));
//   return ircluster;
// }

// template<class parameters>
// std::vector<std::vector<double> >& cluster<parameters>::get_affine_cluster()
// {
//   static std::vector<std::vector<double> > affine_cluster(0,std::vector<double>(parameters::DIMENSION,0.));
//   return affine_cluster;
// }

// template<class parameters>
// std::vector<std::vector<int> >& cluster<parameters>::get_integer_cluster()
// {
//   static std::vector<std::vector<int> > affine_cluster(0,std::vector<int>(parameters::DIMENSION,0.));
//   return affine_cluster;
// }

template<class parameters>
double& cluster<parameters>::get_r_volume()
{
  static double r_volume = 0;
  return r_volume;
}

template<class parameters>
double& cluster<parameters>::get_k_volume()
{
  static double k_volume = 0;
  return k_volume;
}

template<class parameters>
int cluster<parameters>::get_cluster_size()
{
  return get_r_cluster().size();
}

template<class parameters>
void cluster<parameters>::print()
{
  std::stringstream ss;

  if(parameters::DIMENSION==2)
    print_2D(ss);

  if(parameters::DIMENSION==3)
    print_3D(ss);

  cout << ss.str() << endl;
}

template<class parameters>
template<class stream_type>
void cluster<parameters>::to_JSON(stream_type& ss)
{
  if(parameters::DIMENSION==2)
    to_JSON_2D(ss);
  
  if(parameters::DIMENSION==3)
    to_JSON_3D(ss);
}

template<class parameters>
template<class stream_type>
void cluster<parameters>::print_2D(stream_type& ss)
{
  ss << scientific;
  ss.precision(6);

  ss << "full-cluster-size        :\t" << get_cluster_size() << "\t" << endl;

//   ss << "SUPER-BASIS" << endl;
//   ss << "A : " << get_r_super_basis()[0][0] << ";" << get_r_super_basis()[0][1] << "\t" << get_k_super_basis()[0][0] << ";" << get_k_super_basis()[0][1] << endl;
//   ss << "B : " << get_r_super_basis()[1][0] << ";" << get_r_super_basis()[1][1] << "\t" << get_k_super_basis()[1][0] << ";" << get_k_super_basis()[1][1] << endl;
//   ss << endl; 
  
  ss << "FULL CLUSTER" << endl;

  for(int i=0; i<int(get_cluster_size()); i++){

    ss << "\t" << i;

    for(int j=0; j<parameters::DIMENSION; j++){
      ss << "\t" << get_r_cluster()[i][j];
    }


    for(int j=0; j<parameters::DIMENSION; j++){
      ss << "\t" << get_k_cluster()[i][j];
    }

    ss << "\t" << endl;
 }
  
  ss << endl;

  ss << "\n\n\t k-space symmetries : \n\n";
  {
    for(int i=0; i<k_dmn_t::dmn_size(); ++i){
      for(int j=0; j<b_dmn_t::dmn_size(); ++j){
	
	ss << "\t" << i << ", " << j << "\t|\t";
	
// 	for(int l=0; l<sym_super_cell_dmn_t::dmn_size(); ++l)
// 	  ss << "\t" << k_cluster_type::get_symmetry_matrix()(i,j, l).first 
// 	       << ", " << k_cluster_type::get_symmetry_matrix()(i,j, l).second;

	for(int l=0; l<sym_super_cell_dmn_t::dmn_size(); ++l)
	  ss << "\t" << cluster_symmetry<k_cluster_type>::get_symmetry_matrix()(i,j, l).first 
	     << ", " << cluster_symmetry<k_cluster_type>::get_symmetry_matrix()(i,j, l).second;
	
	ss << "\n";
      }
    }
    ss << "\n";
  }

  ss << "\n\n\t r-space symmetries : \n\n";
  {
    for(int i=0; i<r_dmn_t::dmn_size(); ++i){
      for(int j=0; j<b_dmn_t::dmn_size(); ++j){
	
	ss << "\t" << i << ", " << j << "\t|\t";
	
// 	for(int l=0; l<sym_super_cell_dmn_t::dmn_size(); ++l)
// 	  ss << "\t" << r_cluster_type::get_symmetry_matrix()(i,j, l).first 
// 	       << ", " << r_cluster_type::get_symmetry_matrix()(i,j, l).second;

	for(int l=0; l<sym_super_cell_dmn_t::dmn_size(); ++l)
	  ss << "\t" << cluster_symmetry<r_cluster_type>::get_symmetry_matrix()(i,j, l).first 
	     << ", " << cluster_symmetry<r_cluster_type>::get_symmetry_matrix()(i,j, l).second;
	
	ss << "\n";
      }
    }
   ss << "\n";
  }
}

template<class parameters>
template<class stream_type>
void cluster<parameters>::print_3D(stream_type& ss)
{
  ss<<scientific;
  ss.precision(6);

  ss << "full-cluster-size        :\t"        << get_cluster_size() << "\t" << endl;
//   ss << "irreducible-cluster-size :\t" << get_irreducible_cluster_size() << "\t" << endl << endl;

//   ss << endl; 
//   ss << "BASIS (rows)" << endl;
//   ss << "  " << get_r_basis()[0][0] << "\t" << get_r_basis()[0][1] << "\t" << get_r_basis()[0][2] << "\t" << get_k_basis()[0][0] << "\t" << get_k_basis()[0][1] << "\t" << get_k_basis()[0][2] << "\t"<< endl;
//   ss << "  " << get_r_basis()[1][0] << "\t" << get_r_basis()[1][1] << "\t" << get_r_basis()[1][2] << "\t" << get_k_basis()[1][0] << "\t" << get_k_basis()[1][1] << "\t" << get_k_basis()[1][2] << "\t"<< endl;
//   ss << "  " << get_r_basis()[2][0] << "\t" << get_r_basis()[2][1] << "\t" << get_r_basis()[2][2] << "\t" << get_k_basis()[2][0] << "\t" << get_k_basis()[2][1] << "\t" << get_k_basis()[2][2] << "\t"<< endl;
//   ss << endl; 

//   ss << endl << "r-volume : " << get_r_volume() << endl; 
//   ss << endl << "k-volume : " << get_k_volume() << endl; 

//   ss << endl; 
//   ss << "SUPER-BASIS (rows)" << endl;
//   ss << "" << get_r_super_basis()[0][0] << "\t" << get_r_super_basis()[0][1] << "\t" << get_r_super_basis()[0][2] << "\t" << get_k_super_basis()[0][0] << "\t" << get_k_super_basis()[0][1] << "\t" << get_k_super_basis()[0][2] << "\t"<< endl;
//   ss << "" << get_r_super_basis()[1][0] << "\t" << get_r_super_basis()[1][1] << "\t" << get_r_super_basis()[1][2] << "\t" << get_k_super_basis()[1][0] << "\t" << get_k_super_basis()[1][1] << "\t" << get_k_super_basis()[1][2] << "\t"<< endl;
//   ss << "" << get_r_super_basis()[2][0] << "\t" << get_r_super_basis()[2][1] << "\t" << get_r_super_basis()[2][2] << "\t" << get_k_super_basis()[2][0] << "\t" << get_k_super_basis()[2][1] << "\t" << get_k_super_basis()[2][2] << "\t"<< endl;
//   ss << endl; 
  
  ss << "FULL CLUSTER" << endl;

  for(int i=0; i<int(get_cluster_size()); i++){

    ss << "\t" << i;

    for(int j=0; j<parameters::DIMENSION; j++){
      ss << "\t" << get_r_cluster()[i][j];
    }

    for(int j=0; j<parameters::DIMENSION; j++){
      ss << "\t" << get_k_cluster()[i][j];
    }
    
    ss << "\t" << endl;
  }
  
  ss << endl;

  ss << "\n\n\t k-space symmetries : \n\n";
  {
    for(int i=0; i<k_dmn_t::dmn_size(); ++i){
      for(int j=0; j<b_dmn_t::dmn_size(); ++j){
	
	ss << "\t" << i << ", " << j << "\t|\t";
	
// 	for(int l=0; l<sym_super_cell_dmn_t::dmn_size(); ++l)
// 	  ss << "\t" << k_cluster_type::get_symmetry_matrix()(i,j, l).first 
// 	       << ", " << k_cluster_type::get_symmetry_matrix()(i,j, l).second;
	
	for(int l=0; l<sym_super_cell_dmn_t::dmn_size(); ++l)
	  ss << "\t" << cluster_symmetry<k_cluster_type>::get_symmetry_matrix()(i,j, l).first 
	     << ", " << cluster_symmetry<k_cluster_type>::get_symmetry_matrix()(i,j, l).second;

	ss << "\n";
      }
    }
    ss << "\n";
  }

  ss << "\n\n\t r-space symmetries : \n\n";
  {
    for(int i=0; i<r_dmn_t::dmn_size(); ++i){
      for(int j=0; j<b_dmn_t::dmn_size(); ++j){
	
	ss << "\t" << i << ", " << j << "\t|\t";
	
// 	for(int l=0; l<sym_super_cell_dmn_t::dmn_size(); ++l)
// 	  ss << "\t" << r_cluster_type::get_symmetry_matrix()(i,j, l).first 
// 	       << ", " << r_cluster_type::get_symmetry_matrix()(i,j, l).second;
	
	for(int l=0; l<sym_super_cell_dmn_t::dmn_size(); ++l)
	  ss << "\t" << cluster_symmetry<r_cluster_type>::get_symmetry_matrix()(i,j, l).first 
	     << ", " << cluster_symmetry<r_cluster_type>::get_symmetry_matrix()(i,j, l).second;

	ss << "\n";
      }
    }
   ss << "\n";
  }
}

template<class parameters>
template<class stream_type>
void cluster<parameters>::to_JSON_2D(stream_type& ss)
{
  ss << "\"cluster\" : ";

  ss << "{\n";
  
//   ss << "\"full-cluster-size\"        : " << get_cluster_size()             << ",\n";

//   ss << "\"SUPER-BASIS\" : \n {" ;
//   ss << "\"R\" : {\n";
//   ss << "\"R_0\" : [[" << get_r_super_basis()[0][0] << ", " << get_r_super_basis()[0][1] << "]],\n";
//   ss << "\"R_1\" : [[" << get_r_super_basis()[1][0] << ", " << get_r_super_basis()[1][1] << "]]"    ;
//   ss << "},\n"; 
//   ss << "\"K\" : {\n";
//   ss << "\"K_0\" : [["  << get_k_super_basis()[0][0] << ", " << get_k_super_basis()[0][1] << "]],\n";
//   ss << "\"K_1\" : [["  << get_k_super_basis()[1][0] << ", " << get_k_super_basis()[1][1] << "]] \n" ;
//   ss << "}\n"; 
//   ss << "},\n";

//   ss << "\"BASIS\" : \n {" ;
//   ss << "\"R\" : {\n";
//   ss << "\"R_0\" : [[" << get_r_basis()[0][0] << ", " << get_r_basis()[0][1] << "]],\n";
//   ss << "\"R_1\" : [[" << get_r_basis()[1][0] << ", " << get_r_basis()[1][1] << "]]"    ;
//   ss << "},\n"; 
//   ss << "\"K\" : {\n";
//   ss << "\"K_0\" : [["  << get_k_basis()[0][0] << ", " << get_k_basis()[0][1] << "]],\n";
//   ss << "\"K_1\" : [["  << get_k_basis()[1][0] << ", " << get_k_basis()[1][1] << "]] \n" ;
//   ss << "}\n"; 
//   ss << "},\n";

  ss << "\"FULL-CLUSTER\" : \n {";
  ss << "\"R\" : [\n";
  for(int i=0; i<int(get_cluster_size()); i++){
    ss << "[";
    for(int j=0; j<parameters::DIMENSION; j++)
      if(j == parameters::DIMENSION-1)
	if(i == int(get_cluster_size())-1)
	  ss << get_r_cluster()[i][j] << "]\n";
	else
	  ss << get_r_cluster()[i][j] << "],\n";
       else
	ss << get_r_cluster()[i][j] << ", ";
  }
  ss << "],\n"; 
  ss << "\"K\" : [\n";
  for(int i=0; i<int(get_cluster_size()); i++){
    ss << "[";
    for(int j=0; j<parameters::DIMENSION; j++)
      if(j == parameters::DIMENSION-1)
	if(i == int(get_cluster_size())-1)
	  ss << get_k_cluster()[i][j] << "]\n";
	else
	  ss << get_k_cluster()[i][j] << "],\n";
      else
	ss << get_k_cluster()[i][j] << ", ";
  }
  ss << "]\n"; 
  ss << "}";

  ss << "}";
}

template<class parameters>
template<class stream_type>
void cluster<parameters>::to_JSON_3D(stream_type& ss)
{
  ss.precision(12);
  ss << fixed;

  ss << "\"cluster\" : ";

  ss << "{\n";
  
//   ss << "\"full-cluster-size\"        : " << get_cluster_size()             << ",\n";

//   ss << "\"BASIS\" : \n {" ;
//   ss << "\"R\" : {\n";
//   ss << "\"R_0\" : [" <<get_r_basis()[0][0] << ", " << get_r_basis()[0][1] << ", " << get_r_basis()[0][2] << "],\n";
//   ss << "\"R_1\" : [" << get_r_basis()[1][0] << ", " << get_r_basis()[1][1] << ", " << get_r_basis()[1][2] << "],\n";
//   ss << "\"R_2\" : [" << get_r_basis()[2][0] << ", " << get_r_basis()[2][1] << ", " << get_r_basis()[2][2] << "] \n"    ;
//   ss << "},\n"; 
//   ss << "\"K\" : {\n";
//   ss << "\"K_0\" : ["  << get_k_basis()[0][0] << ", " << get_k_basis()[0][1] << ", " << get_k_basis()[0][2] << "],\n";
//   ss << "\"K_1\" : ["  << get_k_basis()[1][0] << ", " << get_k_basis()[1][1] << ", " << get_k_basis()[1][2] << "],\n" ;
//   ss << "\"K_2\" : ["  << get_k_basis()[2][0] << ", " << get_k_basis()[2][1] << ", " << get_k_basis()[2][2] << "] \n" ;
//   ss << "}\n"; 
//   ss << "},\n";

//   ss << "\"SUPER-BASIS\" : \n {" ;
//   ss << "\"R\" : {\n";
//   ss << "\"R_0\" : [" << get_r_super_basis()[0][0] << ", " << get_r_super_basis()[0][1] << ", " << get_r_super_basis()[0][2] << "],\n";
//   ss << "\"R_1\" : [" << get_r_super_basis()[1][0] << ", " << get_r_super_basis()[1][1] << ", " << get_r_super_basis()[1][2] << "],\n";
//   ss << "\"R_2\" : [" << get_r_super_basis()[2][0] << ", " << get_r_super_basis()[2][1] << ", " << get_r_super_basis()[2][2] << "] \n";
//   ss << "},\n"; 
//   ss << "\"K\" : {\n";
//   ss << "\"K_0\" : ["  << get_k_super_basis()[0][0] << ", " << get_k_super_basis()[0][1] << ", " << get_k_super_basis()[0][2] << "],\n";
//   ss << "\"K_1\" : ["  << get_k_super_basis()[1][0] << ", " << get_k_super_basis()[1][1] << ", " << get_k_super_basis()[1][2] << "],\n";
//   ss << "\"K_2\" : ["  << get_k_super_basis()[2][0] << ", " << get_k_super_basis()[2][1] << ", " << get_k_super_basis()[2][2] << "] \n";
//   ss << "}\n"; 
//   ss << "},\n";

  ss << "\"FULL-CLUSTER\" : \n {";
  ss << "\"R\" : [\n";
  for(int i=0; i<int(get_cluster_size()); i++){
    ss << "[";
    for(int j=0; j<parameters::DIMENSION; j++)
      if(j == parameters::DIMENSION-1)
	if(i == int(get_cluster_size())-1)
	  ss << std::fixed << std::setprecision(12) << get_r_cluster()[i][j] << "]\n";
	else
	  ss << std::fixed << std::setprecision(12) << get_r_cluster()[i][j] << "],\n";
       else
	ss << std::fixed << std::setprecision(12) << get_r_cluster()[i][j] << ", ";
  }
  ss << "],\n"; 
  ss << "\"K\" : [\n";
  for(int i=0; i<int(get_cluster_size()); i++){
    ss << "[";
    for(int j=0; j<parameters::DIMENSION; j++)
      if(j == parameters::DIMENSION-1)
	if(i == int(get_cluster_size())-1)
	  ss << std::fixed << std::setprecision(12) << get_k_cluster()[i][j] << "]\n";
	else
	  ss << std::fixed << std::setprecision(12) << get_k_cluster()[i][j] << "],\n";
      else
	ss << std::fixed << std::setprecision(12) << get_k_cluster()[i][j] << ", ";
  }
  ss << "]\n"; 
  ss << "}";

  ss << "}";

}








#endif
