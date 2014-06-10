//-*-C++-*-

#ifndef CLUSTER_INITIALIZER_2D_H
#define CLUSTER_INITIALIZER_2D_H

/*!
 *  \author: Peter Staar
 */
template< class base_cluster, class point_group_type>
class cluster_initializer<2 , base_cluster, point_group_type>
{
  const static int DIMENSION = 2;
 
public: 

  static void initialize(std::vector<std::vector<int> > Bett_matrix,
			 double* r_matrix, 
			 double* k_matrix);

  static void initialize(double* r_matrix, double* k_matrix, std::vector<int>& grid_size);

  static void check();

  static void print();

private:

  static void construct_bases(std::vector<int> A, 
			      std::vector<int> B,
			      double* r_matrix, 
			      double* k_matrix);

  static void construct_r_cluster();
  static void construct_k_cluster();

  static void construct_bases(double* r_matrix, double* k_matrix, std::vector<int>& grid_size);
  static void construct_clusters(std::vector<int>& grid_size);

};

template< class base_cluster, class point_group_type>
void cluster_initializer<2 , base_cluster, point_group_type>::initialize(std::vector<std::vector<int> > Bett_matrix,
									 double* r_matrix, 
									 double* k_matrix)
{
  cluster_initializer<2, base_cluster, point_group_type>::construct_bases(Bett_matrix[0], Bett_matrix[1], r_matrix, k_matrix);
  cluster_initializer<2, base_cluster, point_group_type>::construct_r_cluster();
  cluster_initializer<2, base_cluster, point_group_type>::construct_k_cluster();

  cluster_reduction<base_cluster, point_group_type> cluster_reduction_obj;
  cluster_reduction_obj.execute();
}

template< class base_cluster, class point_group_type>
void cluster_initializer<2 , base_cluster, point_group_type>::initialize(double* r_matrix, double* k_matrix, std::vector<int>& grid_size)
{
  construct_bases(r_matrix, k_matrix, grid_size);
  construct_clusters(grid_size);
}

template< class base_cluster, class point_group_type>
void cluster_initializer<2 , base_cluster, point_group_type>::construct_bases(std::vector<int> A,
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

  super_basis[0][0] = A[0]*basis[0][0]+A[1]*basis[1][0]; 
  super_basis[0][1] = A[0]*basis[0][1]+A[1]*basis[1][1]; 

  super_basis[1][0] = B[0]*basis[0][0]+B[1]*basis[1][0]; 
  super_basis[1][1] = B[0]*basis[0][1]+B[1]*basis[1][1];

  std::vector<std::vector<double> >& k_basis     = base_cluster::get_k_super_basis_new();

  for(int i=0; i<2; i++){
    k_basis[0][i] = k_matrix[i+0];
    k_basis[1][i] = k_matrix[i+2];
  }

  std::vector<std::vector<double> >& k_super_basis = base_cluster::get_k_basis_new();

  double det = super_basis[0][0]*super_basis[1][1]-super_basis[1][0]*super_basis[0][1];

  k_super_basis[0][0] =  2.*M_PI*super_basis[1][1]/det; 
  k_super_basis[0][1] = -2.*M_PI*super_basis[1][0]/det; 

  k_super_basis[1][0] = -2.*M_PI*super_basis[0][1]/det; 
  k_super_basis[1][1] =  2.*M_PI*super_basis[0][0]/det;

  base_cluster::get_r_volume() = fabs(  basis[0][0]*  basis[1][1]-  basis[1][0]*  basis[0][1]);
  base_cluster::get_k_volume() = fabs(k_basis[0][0]*k_basis[1][1]-k_basis[1][0]*k_basis[0][1]);

  double result;

  for(int i=0; i<2; i++){
    for(int j=0; j<2; j++){
      result = 0;
      for(int l=0; l<2; l++)
	result += basis[i][l]*k_basis[j][l]/(2.*M_PI);
      
      if(i==j)
	result -= 1.;
      
      if(fabs(result) > 1.e-6)
	throw std::logic_error(__FUNCTION__);
    }
  }
  
  for(int i=0; i<2; i++){
    for(int j=0; j<2; j++){
      result = 0;
      for(int l=0; l<2; l++)
	result += super_basis[i][l]*k_super_basis[j][l]/(2.*M_PI);
      
      if(i==j)
	 result -= 1.;
      
      if(fabs(result) > 1.e-6)
	throw std::logic_error(__FUNCTION__);
    }
  }
}

template< class base_cluster, class point_group_type>
void cluster_initializer<2 , base_cluster, point_group_type>::construct_r_cluster()
{
  std::vector<std::vector<double> >& basis       = base_cluster::get_r_basis();
  std::vector<std::vector<double> >& super_basis = base_cluster::get_r_super_basis();
  std::vector<std::vector<double> >& r_cluster   = base_cluster::get_r_cluster();

  std::vector<double> rvec(2,0.);

  double x[2]     = {basis[0][0],basis[0][1]};
  double y[2]     = {basis[1][0],basis[1][1]};

  double det = super_basis[0][0]*super_basis[1][1]- super_basis[1][0]*super_basis[0][1];

  // inverse of matrix [A; B], in column-major order !
  double inv_superbasis[2][2] = { {super_basis[1][1]/det, -super_basis[1][0]/det},
				  {-super_basis[0][1]/det,  super_basis[0][0]/det}};

  for(int xind=-100; xind<101; xind++)
    {
      for(int yind=-100; yind<101; yind++)
	{
	  rvec[0] = xind*x[0]+yind*y[0];
	  rvec[1] = xind*x[1]+yind*y[1];

	  double p_X = inv_superbasis[0][0]*rvec[0] + inv_superbasis[0][1]*rvec[1];
	  double p_Y = inv_superbasis[1][0]*rvec[0] + inv_superbasis[1][1]*rvec[1];
	  
	  if( p_X>-1.e-6 && p_X < 1.-1.e-6 && p_Y>-1.e-6 && p_Y < 1.-1.e-6){
	    r_cluster.push_back(rvec);
	  }
	}
    }

  sort(r_cluster.begin(), r_cluster.end(), VECTOR_OPERATIONS::IS_LARGER_VECTOR);
}

template< class base_cluster, class point_group_type>
void cluster_initializer<2 , base_cluster, point_group_type>::construct_k_cluster()
{
  std::vector<std::vector<double> >& k_super_basis      = base_cluster::get_k_basis_new();
  std::vector<std::vector<double> >& k_basis            = base_cluster::get_k_super_basis_new();
  std::vector<std::vector<double> >& k_cluster          = base_cluster::get_k_cluster();

  std::vector<double> k_vec(2,0.);

  double det = k_basis[0][0]*k_basis[1][1] - k_basis[1][0]*k_basis[0][1];

  double inv_basis[2][2] = { {k_basis[1][1]/det, -k_basis[1][0]/det},
			     {-k_basis[0][1]/det,  k_basis[0][0]/det}};
  
  for(int xind=-100; xind<101; xind++)
    {
      for(int yind=-100; yind<101; yind++)
	{
	  k_vec[0] = k_super_basis[0][0]*xind+k_super_basis[1][0]*yind;
	  k_vec[1] = k_super_basis[0][1]*xind+k_super_basis[1][1]*yind;

	  double p_KX = inv_basis[0][0]*k_vec[0] + inv_basis[0][1]*k_vec[1];
	  double p_KY = inv_basis[1][0]*k_vec[0] + inv_basis[1][1]*k_vec[1];
	  
	  if( p_KX>-1.e-6 && p_KX < 1.-1.e-6 && p_KY>-1.e-6 && p_KY < 1.-1.e-6){
	    k_cluster.push_back(k_vec);
	  }
	}
    }

  sort(k_cluster.begin(), k_cluster.end(), VECTOR_OPERATIONS::IS_LARGER_VECTOR);
}





template< class base_cluster, class point_group_type>
void cluster_initializer<2 , base_cluster, point_group_type>::construct_bases(double* r_matrix, double* k_matrix, std::vector<int>& grid_size)
{
  std::vector<std::vector<double> >& basis         = base_cluster::get_r_basis();
  std::vector<std::vector<double> >& super_basis   = base_cluster::get_r_super_basis();
  std::vector<std::vector<double> >& k_basis       = base_cluster::get_k_super_basis_new();
  std::vector<std::vector<double> >& k_super_basis = base_cluster::get_k_basis_new();

  for(int i=0; i<2; i++)
    {
      basis[0][i] = r_matrix[i+0];
      basis[1][i] = r_matrix[i+2];

      super_basis[0][i] = grid_size[0]*r_matrix[i+0];
      super_basis[1][i] = grid_size[1]*r_matrix[i+2];

      k_basis[0][i] = k_matrix[i+0];
      k_basis[1][i] = k_matrix[i+2];

      k_super_basis[0][i] = 1./double(grid_size[0])*k_matrix[i+0];
      k_super_basis[1][i] = 1./double(grid_size[1])*k_matrix[i+2];
     }

  base_cluster::get_r_volume() = fabs(  basis[0][0]*  basis[1][1]-  basis[1][0]*  basis[0][1]);
  base_cluster::get_k_volume() = fabs(k_basis[0][0]*k_basis[1][1]-k_basis[1][0]*k_basis[0][1]);
}

template< class base_cluster, class point_group_type>
void cluster_initializer<2 , base_cluster, point_group_type>::construct_clusters(std::vector<int>& grid_size)
{
  std::vector<std::vector<double> >& basis         = base_cluster::get_r_basis();
  std::vector<std::vector<double> >& k_super_basis = base_cluster::get_k_basis_new();
  
  for(int i=0; i<grid_size[1]; i++)
    {
      for(int j=0; j<grid_size[0]; j++)
	{
	  std::vector<double> r_vec(2,0);
	  r_vec[0] = j*basis[0][0] + i*basis[1][0];
	  r_vec[1] = j*basis[0][1] + i*basis[1][1];
	  
	  base_cluster::get_r_cluster().push_back(r_vec);
	  
	  std::vector<double> k_vec(2,0);
	  k_vec[0] = j*k_super_basis[0][0] + i*k_super_basis[1][0];
	  k_vec[1] = j*k_super_basis[0][1] + i*k_super_basis[1][1];
	  
	  base_cluster::get_k_cluster().push_back(k_vec);
	  
	  std::vector<double> affine_vec(2,0);
	  affine_vec[0] = double(j)/double(grid_size[0]); 
	  affine_vec[1] = double(i)/double(grid_size[1]);
	}
    }
}

template< class base_cluster, class point_group_type>
void cluster_initializer<2 , base_cluster, point_group_type>::check()
{
  for(int i=0; i<int(base_cluster::get_k_cluster().size()); i++)
    {
      std::vector<double > k_vec1 = base_cluster::get_k_cluster()[i];
      std::vector<double > k_vec2 = base_cluster::parameter_type::get_k_vectors()(i);
      
      if(fabs(k_vec1[0]-k_vec2[0]) > 1.e-6)
	cout << "inconsistent LDA-cluster and k_cluster @ " << __PRETTY_FUNCTION__ << endl;
      
      if(fabs(k_vec1[1]-k_vec2[1]) > 1.e-6)
	cout << "inconsistent LDA-cluster and k_cluster @ " << __PRETTY_FUNCTION__ << endl;
    }
}

template< class base_cluster, class point_group_type>
void cluster_initializer<2, base_cluster, point_group_type>::print()
{
  cout << "SUPER-BASIS" << endl;
  cout << "A : " << base_cluster::get_r_super_basis()[0][0] << ";" << base_cluster::get_r_super_basis()[0][1] << "\t";
  cout << base_cluster::get_k_basis_new()[0][0] << ";" << base_cluster::get_k_basis_new()[0][1] << endl;

  cout << "B : " << base_cluster::get_r_super_basis()[1][0] << ";" << base_cluster::get_r_super_basis()[1][1] << "\t";
  cout << base_cluster::get_k_basis_new()[1][0] << ";" << base_cluster::get_k_basis_new()[1][1] << endl;

  cout << endl; 
  
  cout << "FULL CLUSTER" << endl;

  for(int i=0; i<int(base_cluster::get_cluster_size()); i++){
    cout << i  << "\t";

    for(int j=0; j<2; j++){
      cout << "\t" << base_cluster::get_r_cluster()[i][j];
    }

    cout << "\t";
    for(int j=0; j<2; j++){
      cout << "\t" << base_cluster::get_k_cluster()[i][j];
    }

    cout << endl;
  }
}  


#endif
