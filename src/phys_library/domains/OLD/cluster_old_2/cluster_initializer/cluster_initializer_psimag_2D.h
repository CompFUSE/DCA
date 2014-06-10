//-*-C++-*-

/*
 *      Author: peterstaar
 */

#ifndef CLUSTER_INITIALIZER_PSIMAG_2D_H_
#define CLUSTER_INITIALIZER_PSIMAG_2D_H_

template< class base_cluster>
class cluster_initializer<2 , base_cluster, PsiMag_symmetry_2D>
{
  const static int DIMENSION = 2;
  
public: 
  
  static void initialize(std::vector<std::vector<int> > Bett_matrix,
			 double* r_matrix, 
			 double* k_matrix);

    static void initialize(double* r_matrix, 
			   double* k_matrix, 
			   std::vector<int>& grid_size);

private:

  static void construct_bases(std::vector<int> A, 
			      std::vector<int> B,
			      double* r_matrix, 
			      double* k_matrix);

  template<typename geometry_type>
  static void construct_r_cluster(geometry_type& geometry);

  template<typename geometry_type>
  static void construct_k_cluster(geometry_type& geometry);

  static void test(std::vector<int> A, std::vector<int> B);


};

template< class base_cluster>
void cluster_initializer<2 , base_cluster, PsiMag_symmetry_2D>::initialize(std::vector<std::vector<int> > Bett_matrix,
									   double* r_matrix, 
									   double* k_matrix)
{
  base_cluster::get_Bett_matrix() = Bett_matrix;

  construct_bases(Bett_matrix[0], Bett_matrix[1], r_matrix, k_matrix);

  construct_r_cluster(base_cluster::get_geometry());
  construct_k_cluster(base_cluster::get_geometry());
}

template< class base_cluster>
void cluster_initializer<2 , base_cluster, PsiMag_symmetry_2D>::construct_bases(std::vector<int> A,
										std::vector<int> B,
										double* r_matrix, 
										double* k_matrix)
{
  std::vector<std::vector<double> >& basis = base_cluster::get_r_basis();

  for(int i=0; i<2; i++){
      basis[0][i] = r_matrix[i+0];
      basis[1][i] = r_matrix[i+2];
    }

  std::vector<std::vector<double> >& super_basis = base_cluster::get_r_super_basis();
  // super_basis[0]
  super_basis[0][0] = A[0]*basis[0][0]+A[1]*basis[1][0]; 
  super_basis[0][1] = A[0]*basis[0][1]+A[1]*basis[1][1]; 
  // super_basis[1]
  super_basis[1][0] = B[0]*basis[0][0]+B[1]*basis[1][0]; 
  super_basis[1][1] = B[0]*basis[0][1]+B[1]*basis[1][1];

  std::vector<std::vector<double> >& k_basis     = base_cluster::get_k_basis();

  for(int i=0; i<2; i++){
      k_basis[0][i] = k_matrix[i+0];
      k_basis[1][i] = k_matrix[i+2];
     }

  std::vector<std::vector<double> >& k_super_basis = base_cluster::get_k_super_basis();

  double det = super_basis[0][0]*super_basis[1][1]-super_basis[1][0]*super_basis[0][1];
  // k_super_basis[0]
  k_super_basis[0][0] =  2.*M_PI*super_basis[1][1]/det; 
  k_super_basis[0][1] = -2.*M_PI*super_basis[1][0]/det; 
  // k_super_basis[1]
  k_super_basis[1][0] = -2.*M_PI*super_basis[0][1]/det; 
  k_super_basis[1][1] =  2.*M_PI*super_basis[0][0]/det;

  base_cluster::get_r_volume() = basis[0][0]*basis[1][1]-basis[1][0]*basis[0][1];
  base_cluster::get_k_volume() = k_basis[0][0]*k_basis[1][1]-k_basis[1][0]*k_basis[0][1];
}

template< class base_cluster>
template<typename geometry_type>
void cluster_initializer<2 , base_cluster, PsiMag_symmetry_2D>::construct_r_cluster(geometry_type& geometry)
{
//   cout << __FUNCTION__ << endl;
//   cout << int(geometry.rCrystal.numberOfSites) << endl;

  for(int l=0; l<int(geometry.rCrystal.numberOfSites); l++){
    std::vector<double> r(2,0);

    r[0] = geometry.rCrystal.sites(l,0); 
    r[1] = geometry.rCrystal.sites(l,1); 

    //    cout << r[0] << "\t" << r[1] << endl;

    base_cluster::get_r_cluster().push_back(r);

    if(geometry.rCrystal.wedge.wedge2CrystalContains(l))
      base_cluster::get_irreducible_r_cluster().push_back(r);
  }
}

template< class base_cluster>
template<typename geometry_type>
void cluster_initializer<2 , base_cluster, PsiMag_symmetry_2D>::construct_k_cluster(geometry_type& geometry)
{
//   cout << __FUNCTION__ << endl;
//   cout << int(geometry.kCrystal.numberOfSites) << endl;

  for(int l=0; l<int(geometry.kCrystal.numberOfSites); l++){
    std::vector<double> k(2,0);

    k[0] = geometry.kCrystal.sites(l,0); 
    k[1] = geometry.kCrystal.sites(l,1); 
    
    //    cout << k[0] << "\t" << k[1] << endl;

    base_cluster::get_k_cluster().push_back(k);
    
    if(geometry.kCrystal.wedge.wedge2CrystalContains(l))
      base_cluster::get_irreducible_k_cluster().push_back(k);
  }
}




#endif
