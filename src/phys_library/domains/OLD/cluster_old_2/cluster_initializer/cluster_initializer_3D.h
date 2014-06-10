//-*-C++-*-

#ifndef CLUSTER_INITIALIZER_3D_H
#define CLUSTER_INITIALIZER_3D_H

/*!
 *  \author: peter staar
 */
template< class base_cluster, class point_group_type>
class cluster_initializer<3, base_cluster, point_group_type>
{
//   typedef r_cluster<FULL, base_cluster > r_cluster_type;
//   typedef typename r_cluster_type::cluster_reduction_type r_crystal_type;

//   typedef k_cluster<FULL, base_cluster > k_cluster_type;
//   typedef typename k_cluster_type::cluster_reduction_type k_crystal_type;

  const static int DIMENSION = 3;

public: 

  static void initialize(std::vector<std::vector<int> > Bett_matrix,
			 double* r_matrix, 
			 double* k_matrix);

  static void initialize(double* r_matrix, double* k_matrix, std::vector<int>& grid_size);

  static void check();

private:

  static void construct_bases(std::vector<int> A, 
			      std::vector<int> B, 
			      std::vector<int> C,
			      double* r_matrix, 
			      double* k_matrix);

  static void construct_r_cluster();
  static void construct_k_cluster();

  static void construct_bases(double* r_matrix, double* k_matrix, std::vector<int>& grid_size);
  static void construct_clusters(std::vector<int>& grid_size);
};


template< class base_cluster, class point_group_type >
void cluster_initializer<3, base_cluster, point_group_type>::initialize(std::vector<std::vector<int> > Bett_matrix,
									double* r_matrix, 
									double* k_matrix)
{
  cluster_initializer<3, base_cluster, point_group_type>::construct_bases(Bett_matrix[0], Bett_matrix[1], Bett_matrix[2], r_matrix, k_matrix);
  cluster_initializer<3, base_cluster, point_group_type>::construct_r_cluster();
  cluster_initializer<3, base_cluster, point_group_type>::construct_k_cluster();

  cluster_reduction<base_cluster, point_group_type> cluster_reduction_obj;
  cluster_reduction_obj.execute();
}

template< class base_cluster, class point_group_type >
void cluster_initializer<3 , base_cluster, point_group_type>::initialize(double* r_matrix, double* k_matrix, std::vector<int>& grid_size)
{
  construct_bases(r_matrix, k_matrix, grid_size);
  construct_clusters(grid_size);
}

template< class base_cluster, class point_group_type >
void cluster_initializer<3 , base_cluster, point_group_type>::construct_bases(std::vector<int> A, 
									      std::vector<int> B, 
									      std::vector<int> C,
									      double* r_matrix, 
									      double* k_matrix)
{
  std::vector<std::vector<double> >& r_basis = base_cluster::get_r_basis();

  for(int i=0; i<3; i++){
    r_basis[0][i] = r_matrix[i+0]; //cout << A[i] << "\t";
    r_basis[1][i] = r_matrix[i+3]; //cout << B[i] << "\t";
    r_basis[2][i] = r_matrix[i+6]; //cout << C[i] << "\n";
  }

  std::vector<std::vector<double> >& super_basis = base_cluster::get_r_super_basis();

  super_basis[0][0] = A[0]*r_basis[0][0]+A[1]*r_basis[1][0]+A[2]*r_basis[2][0]; 
  super_basis[0][1] = A[0]*r_basis[0][1]+A[1]*r_basis[1][1]+A[2]*r_basis[2][1];
  super_basis[0][2] = A[0]*r_basis[0][2]+A[1]*r_basis[1][2]+A[2]*r_basis[2][2]; 
  
  super_basis[1][0] = B[0]*r_basis[0][0]+B[1]*r_basis[1][0]+B[2]*r_basis[2][0]; 
  super_basis[1][1] = B[0]*r_basis[0][1]+B[1]*r_basis[1][1]+B[2]*r_basis[2][1];
  super_basis[1][2] = B[0]*r_basis[0][2]+B[1]*r_basis[1][2]+B[2]*r_basis[2][2]; 

  super_basis[2][0] = C[0]*r_basis[0][0]+C[1]*r_basis[1][0]+C[2]*r_basis[2][0]; 
  super_basis[2][1] = C[0]*r_basis[0][1]+C[1]*r_basis[1][1]+C[2]*r_basis[2][1];
  super_basis[2][2] = C[0]*r_basis[0][2]+C[1]*r_basis[1][2]+C[2]*r_basis[2][2]; 

  std::vector<std::vector<double> >& k_basis     = base_cluster::get_k_basis();

  for(int i=0; i<3; i++){
    k_basis[0][i] = k_matrix[i+0]; 
    k_basis[1][i] = k_matrix[i+3]; 
    k_basis[2][i] = k_matrix[i+6]; 
  }

  std::vector<std::vector<double> >& k_super_basis = base_cluster::get_k_super_basis();

  double det = super_basis[0][0]*super_basis[1][1]*super_basis[2][2]
                + super_basis[1][0]*super_basis[2][1]*super_basis[0][2]
                + super_basis[0][1]*super_basis[1][2]*super_basis[2][0]
              - super_basis[0][2]*super_basis[1][1]*super_basis[2][0]
                - super_basis[1][0]*super_basis[0][1]*super_basis[2][2]
                - super_basis[1][2]*super_basis[2][1]*super_basis[0][0];

  // the k-superbasis is the vector product of the superbasis vectors! 

  k_super_basis[0][0] =  2.*M_PI/det*(super_basis[1][1]*super_basis[2][2]-super_basis[2][1]*super_basis[1][2]);
  k_super_basis[0][1] = -2.*M_PI/det*(super_basis[1][0]*super_basis[2][2]-super_basis[2][0]*super_basis[1][2]);
  k_super_basis[0][2] =  2.*M_PI/det*(super_basis[1][0]*super_basis[2][1]-super_basis[2][0]*super_basis[1][1]);

  k_super_basis[1][0] = -2.*M_PI/det*(super_basis[0][1]*super_basis[2][2]-super_basis[2][1]*super_basis[0][2]);
  k_super_basis[1][1] =  2.*M_PI/det*(super_basis[0][0]*super_basis[2][2]-super_basis[2][0]*super_basis[0][2]);
  k_super_basis[1][2] = -2.*M_PI/det*(super_basis[0][0]*super_basis[2][1]-super_basis[2][0]*super_basis[0][1]);

  k_super_basis[2][0] =  2.*M_PI/det*(super_basis[0][1]*super_basis[1][2]-super_basis[1][1]*super_basis[0][2]);
  k_super_basis[2][1] = -2.*M_PI/det*(super_basis[0][0]*super_basis[1][2]-super_basis[1][0]*super_basis[0][2]);
  k_super_basis[2][2] =  2.*M_PI/det*(super_basis[0][0]*super_basis[1][1]-super_basis[1][0]*super_basis[0][1]);

  base_cluster::get_r_volume() =  fabs(r_basis[0][0]* r_basis[1][1]* r_basis[2][2]
				       +  r_basis[1][0]* r_basis[2][1]* r_basis[0][2]
				       +  r_basis[0][1]* r_basis[1][2]* r_basis[2][0]
				       -  r_basis[0][2]* r_basis[1][1]* r_basis[2][0]
				       -  r_basis[1][0]* r_basis[0][1]* r_basis[2][2]
				       -  r_basis[1][2]* r_basis[2][1]* r_basis[0][0]);

  base_cluster::get_k_volume() = fabs(k_basis[0][0]*k_basis[1][1]*k_basis[2][2]
				      + k_basis[1][0]*k_basis[2][1]*k_basis[0][2]
				      + k_basis[0][1]*k_basis[1][2]*k_basis[2][0]
				      - k_basis[0][2]*k_basis[1][1]*k_basis[2][0]
				      - k_basis[1][0]*k_basis[0][1]*k_basis[2][2]
				      - k_basis[1][2]*k_basis[2][1]*k_basis[0][0]);

   double result = 0;

   for(int i=0; i<3; i++){
     for(int j=0; j<3; j++){
       result = 0;
       for(int l=0; l<3; l++)
	 result += r_basis[i][l]*k_basis[j][l]/(2.*M_PI);
	   
       if(i==j)
	 result -= 1.;

       if(fabs(result) > 1.e-6)
	 throw std::logic_error(__FUNCTION__);
     }
   }

   for(int i=0; i<3; i++){
     for(int j=0; j<3; j++){
       result = 0;
       for(int l=0; l<3; l++)
	 result += super_basis[i][l]*k_super_basis[j][l]/(2.*M_PI);
	   
       if(i==j)
	 result -= 1.;

       if(fabs(result) > 1.e-6)
	 throw std::logic_error(__FUNCTION__);
     }
   }
}



template< class base_cluster, class point_group_type >
void cluster_initializer<3 , base_cluster, point_group_type>::construct_r_cluster()
{
  //cout << __FUNCTION__ << endl;

  std::vector<std::vector<double> >& basis         = base_cluster::get_r_basis();
  std::vector<std::vector<double> >& k_super_basis = base_cluster::get_k_super_basis();
  std::vector<std::vector<double> >& r_cluster     = base_cluster::get_r_cluster();
  

  std::vector<double> rvec(3,0.);
  std::vector<double> affine_vec(3,0.);
  std::vector<int>    int_vec(3,0);

  double x[3]     = {basis[0][0],basis[0][1],basis[0][2]};
  double y[3]     = {basis[1][0],basis[1][1],basis[1][2]};
  double z[3]     = {basis[2][0],basis[2][1],basis[2][2]};

  // inverse of matrix [A; B], in column-major order !
//   double inv_superbasis[3][3] = { {k_super_basis[0][0]/(2.*M_PI), k_super_basis[1][0]/(2.*M_PI), k_super_basis[2][0]/(2.*M_PI)},
// 				  {k_super_basis[0][1]/(2.*M_PI), k_super_basis[1][1]/(2.*M_PI), k_super_basis[2][1]/(2.*M_PI)},
// 				  {k_super_basis[0][2]/(2.*M_PI), k_super_basis[1][2]/(2.*M_PI), k_super_basis[2][2]/(2.*M_PI)}};

  for(int zind=-30; zind<30; zind++)
    {
      for(int yind=-30; yind<30; yind++)
	{
	  for(int xind=-30; xind<30; xind++)	  
	    {
	      rvec[0] = xind*x[0]+yind*y[0]+zind*z[0];
	      rvec[1] = xind*x[1]+yind*y[1]+zind*z[1];
	      rvec[2] = xind*x[2]+yind*y[2]+zind*z[2];

	      double p_X = (k_super_basis[0][0]*rvec[0] + k_super_basis[0][1]*rvec[1] + k_super_basis[0][2]*rvec[2])/(2.*M_PI);
	      double p_Y = (k_super_basis[1][0]*rvec[0] + k_super_basis[1][1]*rvec[1] + k_super_basis[1][2]*rvec[2])/(2.*M_PI);
	      double p_Z = (k_super_basis[2][0]*rvec[0] + k_super_basis[2][1]*rvec[1] + k_super_basis[2][2]*rvec[2])/(2.*M_PI);
	      
	      if( p_X>-1.e-6 && p_X < 1.-1.e-6 && p_Y>-1.e-6 && p_Y < 1.-1.e-6 && p_Z>-1.e-6 && p_Z < 1.-1.e-6){
		r_cluster.push_back(rvec);
	      }
	    }
	}
    }

//   r_crystal_type& r_crystal = r_cluster_type::get_cluster_reduction();
  
//   for(int i=0; i<r_cluster_type::get_size(); i++)
//     if(i == r_crystal.get_irreducible_index(i))
//       base_cluster::get_irreducible_r_cluster().push_back(base_cluster::get_r_cluster()[i]);
}

template< class base_cluster, class point_group_type >
void cluster_initializer<3 , base_cluster, point_group_type>::construct_k_cluster()
{
  //cout << __FUNCTION__ << endl;

  std::vector<std::vector<double> >& r_basis   = base_cluster::get_r_basis();
  std::vector<std::vector<double> >& k_basis   = base_cluster::get_k_super_basis();
  std::vector<std::vector<double> >& k_cluster = base_cluster::get_k_cluster();

  std::vector<double> kvec(3,0.);
  
  double kx[3]     = {k_basis[0][0],k_basis[0][1],k_basis[0][2]};
  double ky[3]     = {k_basis[1][0],k_basis[1][1],k_basis[1][2]};
  double kz[3]     = {k_basis[2][0],k_basis[2][1],k_basis[2][2]};

  for(int kzind=-30; kzind<30; kzind++)
    {
      for(int kyind=-30; kyind<30; kyind++)
	{
	  for(int kxind=-30; kxind<30; kxind++)
	    {
	      kvec[0] = kxind*kx[0]+kyind*ky[0]+kzind*kz[0];
	      kvec[1] = kxind*kx[1]+kyind*ky[1]+kzind*kz[1];
	      kvec[2] = kxind*kx[2]+kyind*ky[2]+kzind*kz[2];
	      
	      double p_KX = (r_basis[0][0]*kvec[0] + r_basis[0][1]*kvec[1] + r_basis[0][2]*kvec[2])/(2.*M_PI);
	      double p_KY = (r_basis[1][0]*kvec[0] + r_basis[1][1]*kvec[1] + r_basis[1][2]*kvec[2])/(2.*M_PI);
	      double p_KZ = (r_basis[2][0]*kvec[0] + r_basis[2][1]*kvec[1] + r_basis[2][2]*kvec[2])/(2.*M_PI);
	      
	      if( p_KX>-1.e-6 && p_KX < 1.-1.e-6 && p_KY>-1.e-6 && p_KY < 1.-1.e-6 && p_KZ>-1.e-6 && p_KZ < 1.-1.e-6){
		k_cluster.push_back(kvec);
	      }
	    }
	}
    }

//   for(size_t l=0; l<k_cluster.size(); l++){
//     for(int d=0; d<DIMENSION; d++)
//       k_cluster[l][d] += DCA_phase[d];
//   }

//   k_crystal_type& k_crystal = k_cluster_type::get_cluster_reduction();

//   for(int i=0; i<k_cluster_type::get_size(); i++)
//     if(i == k_crystal.get_irreducible_index(i))
//       base_cluster::get_irreducible_k_cluster().push_back(base_cluster::get_k_cluster()[i]);
}

template< class base_cluster, class point_group_type >
void cluster_initializer<3 , base_cluster, point_group_type>::construct_bases(double* r_matrix, double* k_matrix, std::vector<int>& grid_size)
{
  std::vector<std::vector<double> >& basis         = base_cluster::get_r_basis();
  std::vector<std::vector<double> >& super_basis   = base_cluster::get_r_super_basis();
  std::vector<std::vector<double> >& k_basis       = base_cluster::get_k_basis();
  std::vector<std::vector<double> >& k_super_basis = base_cluster::get_k_super_basis();

  for(int i=0; i<3; i++)
    {
      basis[0][i] = r_matrix[i+0];
      basis[1][i] = r_matrix[i+3];
      basis[2][i] = r_matrix[i+6];

      super_basis[0][i] = grid_size[0]*r_matrix[i+0];
      super_basis[1][i] = grid_size[1]*r_matrix[i+3];
      super_basis[2][i] = grid_size[2]*r_matrix[i+6];

      k_basis[0][i] = k_matrix[i+0];
      k_basis[1][i] = k_matrix[i+3];
      k_basis[2][i] = k_matrix[i+6];

      k_super_basis[0][i] = 1./double(grid_size[0])*k_matrix[i+0];
      k_super_basis[1][i] = 1./double(grid_size[1])*k_matrix[i+3];
      k_super_basis[2][i] = 1./double(grid_size[2])*k_matrix[i+6];
    }

  base_cluster::get_r_volume() =  fabs(basis[0][0]* basis[1][1]* basis[2][2]
				       +  basis[1][0]* basis[2][1]* basis[0][2]
				       +  basis[0][1]* basis[1][2]* basis[2][0]
				       -  basis[0][2]* basis[1][1]* basis[2][0]
				       -  basis[1][0]* basis[0][1]* basis[2][2]
				       -  basis[1][2]* basis[2][1]* basis[0][0]);

  base_cluster::get_k_volume() = fabs(k_basis[0][0]*k_basis[1][1]*k_basis[2][2]
				      + k_basis[1][0]*k_basis[2][1]*k_basis[0][2]
				      + k_basis[0][1]*k_basis[1][2]*k_basis[2][0]
				      - k_basis[0][2]*k_basis[1][1]*k_basis[2][0]
				      - k_basis[1][0]*k_basis[0][1]*k_basis[2][2]
				      - k_basis[1][2]*k_basis[2][1]*k_basis[0][0]);
}

template< class base_cluster, class point_group_type >
void cluster_initializer<3 , base_cluster, point_group_type>::construct_clusters(std::vector<int>& grid_size)
{
  std::vector<std::vector<double> >& basis         = base_cluster::get_r_basis();
  //std::vector<std::vector<double> >& super_basis   = base_cluster::get_r_super_basis();
  //std::vector<std::vector<double> >& k_basis       = base_cluster::get_k_basis();
  std::vector<std::vector<double> >& k_super_basis = base_cluster::get_k_super_basis();

  for(int i=0; i<grid_size[2]; i++)
    {
      for(int j=0; j<grid_size[1]; j++)
	{
	  for(int l=0; l<grid_size[0]; l++)
	    {
	      std::vector<double> r_vec(3,0);
	      r_vec[0] = l*basis[0][0] + j*basis[1][0] + i*basis[2][0];
	      r_vec[1] = l*basis[0][1] + j*basis[1][1] + i*basis[2][1];
	      r_vec[2] = l*basis[0][2] + j*basis[1][2] + i*basis[2][2];

	      base_cluster::get_r_cluster().push_back(r_vec);

	      std::vector<double> k_vec(3,0);
	      k_vec[0] = l*k_super_basis[0][0] + j*k_super_basis[1][0] + i*k_super_basis[2][0];
	      k_vec[1] = l*k_super_basis[0][1] + j*k_super_basis[1][1] + i*k_super_basis[2][1];
	      k_vec[2] = l*k_super_basis[0][2] + j*k_super_basis[1][2] + i*k_super_basis[2][2];

	      base_cluster::get_k_cluster().push_back(k_vec);

	      std::vector<double> affine_vec(3,0);
	      affine_vec[0] = double(l)/double(grid_size[0]); 
	      affine_vec[1] = double(j)/double(grid_size[1]);
	      affine_vec[2] = double(i)/double(grid_size[2]);
	      
	      //base_cluster::get_affine_cluster().push_back(affine_vec);
	    }
	}
    }
}

template< class base_cluster, class point_group_type >
void cluster_initializer<3 , base_cluster, point_group_type>::check()
{
  for(int i=0; i<int(base_cluster::get_k_cluster().size()); i++)
    {
      std::vector<double > k_vec1 = base_cluster::get_k_cluster()[i];
      std::vector<double > k_vec2 = base_cluster::parameter_type::get_k_vectors()(i);
      
      if(fabs(k_vec1[0]-k_vec2[0]) > 1.e-6)
	cout << "inconsistent LDA-cluster and k_cluster @ k_x --> nb : " << i << "\t" << fabs(k_vec1[0]-k_vec2[0]) << "\t" << __PRETTY_FUNCTION__ << endl;
      
      if(fabs(k_vec1[1]-k_vec2[1]) > 1.e-6)
	cout << "inconsistent LDA-cluster and k_cluster @ k_y --> nb : " << i << "\t" << fabs(k_vec1[1]-k_vec2[1]) << "\t" << __PRETTY_FUNCTION__ << endl;
      
      if(fabs(k_vec1[2]-k_vec2[2]) > 1.e-6)
	cout << "inconsistent LDA-cluster and k_cluster @ k_z --> nb : " << i << "\t" << fabs(k_vec1[2]-k_vec2[2]) << "\t" << __PRETTY_FUNCTION__ << endl;
    }
}
#endif
