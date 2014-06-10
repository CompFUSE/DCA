//-*-C++-*-

/*
 *      Author: peter staar
 */

template<cluster_representation_type representation, typename parameters>
void tetrahedron_mesh<k_cluster<representation , cluster<parameters > > >::initialize_2D()
{
  std::vector<double>& b0 = k_cluster_type::get_super_basis()[0];
  std::vector<double>& b1 = k_cluster_type::get_super_basis()[1];

  simplex<DIMENSION>  s;
  std::vector<double> B(2,0);

  for(int t1=0; t1<2; t1++){
    for(int t0=0; t0<2; t0++){

      int l0 = t0==0? -1:1;
      int l1 = t1==0? -1:1;
      
      B[0] = l0*b0[0] + l1*b1[0];
      B[1] = l0*b0[1] + l1*b1[1];

      s.k_vec = B;

      int i0 = t0;
      int i1 = t1;
      s.index[0] = (i0+1)%2 + 2*i1;
      s.index[1] =  i0      + 2*((i1+1)%2);

      simplices.push_back(s);
    }
  }
}

template<cluster_representation_type representation , typename parameters >
void tetrahedron_mesh<k_cluster<representation , cluster<parameters > > >::make_convex_hull_2D()
{
  std::vector<double> B(2,0);
  std::vector<double> k(2,0);
  std::vector<double> K(2,0);
  
  std::vector<double>& b0 = k_cluster_type::get_super_basis()[0];
  std::vector<double>& b1 = k_cluster_type::get_super_basis()[1];

  bool go_on = false;

  while(!go_on)
    {
      go_on = true;

      for(size_t l=0; l<simplices.size(); l++)
	{
	  k = simplices[l].k_vec;

	  for(int t0=-1; t0<=1; t0++){
	    for(int t1=-1; t1<=1; t1++){
	      if((t0!=0 || t1!=0) && go_on)
		{
		  B[0] = t0*b0[0] + t1*b1[0];
		  B[1] = t0*b0[1] + t1*b1[1];
		  
		  double diff_k_K = L2_norm(k,K);
		  double diff_k_B = L2_norm(k,B);

		  if(diff_k_K > diff_k_B+1.e-6)
		    {
		      go_on = false;
		      update_convex_hull_2D(l, B);
		    }
		}
	    }
	    
	    if(!go_on)
	      break;
	  }
	}
    }

  for(size_t i=0; i<simplices.size(); i++)
    for(size_t j=i+1; j<simplices.size(); j++)
      if(L2_norm(simplices[i].k_vec, simplices[j].k_vec) < 1.e-6)
	remove_simplex(j);
}

template<cluster_representation_type representation , typename parameters >
void tetrahedron_mesh<k_cluster<representation , cluster<parameters > > >::update_convex_hull_2D(int k_ind, 
												std::vector<double> B)
{  
  double lambda_0 = find_lambda_2D(simplices[k_ind].k_vec, simplices[simplices[k_ind].index[0]].k_vec, B);
  double lambda_1 = find_lambda_2D(simplices[k_ind].k_vec, simplices[simplices[k_ind].index[1]].k_vec, B);

  if(lambda_0>1.-1.e-6 && lambda_1>1.-1.e-6){
    remove_simplex(k_ind);
    return;
  }
  
  if((lambda_0<1. && lambda_0>0.) && (lambda_1<1 && lambda_1>0)){
    add_new_simplices_2D(k_ind, B);
    return;
  }

  if(lambda_0<1. && lambda_0>0.){
    shift_simplex(k_ind, lambda_0, 0);
    return;
  }

  if(lambda_1<1 && lambda_1>0){
    shift_simplex(k_ind, lambda_1, 1);
    return;
  }

  plot_simplices();
  throw std::logic_error(__FUNCTION__);
}

template<cluster_representation_type representation , typename parameters >
double tetrahedron_mesh<k_cluster<representation , cluster<parameters > > >::find_lambda_2D(std::vector<double> k,
											    std::vector<double> n,
											    std::vector<double> B)
{
  return -(std::pow(B[0],2) + std::pow(B[1],2) - 2*B[0]*k[0] - 2*B[1]*k[1])/(2.*(B[0]*k[0] + B[1]*k[1] - B[0]*n[0] - B[1]*n[1]));
}


template<cluster_representation_type representation , typename parameters >
void tetrahedron_mesh<k_cluster<representation , cluster<parameters > > >::remove_simplex(int k_ind)
{
  {
    simplex<DIMENSION>& s0 = simplices[simplices[k_ind].index[0]];
    
    if(s0.index[0] == k_ind)
      s0.index[0] = simplices[k_ind].index[1];
    else
      s0.index[1] = simplices[k_ind].index[1];
  }

  {
    simplex<DIMENSION>& s1 = simplices[simplices[k_ind].index[1]];
    
    if(s1.index[0] == k_ind)
      s1.index[0] = simplices[k_ind].index[0];
    else
      s1.index[1] = simplices[k_ind].index[0];
  }

  for(size_t l=0; l<simplices.size(); l++)
    {
      simplices[l].index[0] = simplices[l].index[0]>k_ind? simplices[l].index[0]-1: simplices[l].index[0];
      simplices[l].index[1] = simplices[l].index[1]>k_ind? simplices[l].index[1]-1: simplices[l].index[1];
    }

  simplices.erase(simplices.begin()+k_ind);
}

template<cluster_representation_type representation , typename parameters >
void tetrahedron_mesh<k_cluster<representation , cluster<parameters > > >::shift_simplex(int k_ind, double lambda, int n)
{
  simplices[k_ind].k_vec[0] = simplices[k_ind].k_vec[0] + lambda*(simplices[simplices[k_ind].index[n]].k_vec[0]-simplices[k_ind].k_vec[0]);
  simplices[k_ind].k_vec[1] = simplices[k_ind].k_vec[1] + lambda*(simplices[simplices[k_ind].index[n]].k_vec[1]-simplices[k_ind].k_vec[1]);
}

template<cluster_representation_type representation , typename parameters >
void tetrahedron_mesh<k_cluster<representation , cluster<parameters > > >::add_new_simplices_2D(int k_ind, 
												std::vector<double> B)
{
  std::vector<double> k = simplices[k_ind].k_vec;
  
  std::vector<double> k0 = simplices[simplices[k_ind].index[0]].k_vec;
  std::vector<double> k1 = simplices[simplices[k_ind].index[1]].k_vec;

  simplex<DIMENSION> s_0;
  simplex<DIMENSION> s_1;

  s_0.k_vec.resize(2,0);
  s_1.k_vec.resize(2,0);

  s_0.k_vec[0] = -(-(std::pow(B[0],2)*k[0]) - std::pow(B[1],2)*k[0] + std::pow(B[0],2)*k0[0] + std::pow(B[1],2)*k0[0] - 2*B[1]*k[1]*k0[0] + 2*B[1]*k[0]*k0[1])/(2.*(B[0]*k[0] + B[1]*k[1] - B[0]*k0[0] - B[1]*k0[1]));
  s_0.k_vec[1] = -(-(std::pow(B[0],2)*k[1]) - std::pow(B[1],2)*k[1] + 2*B[0]*k[1]*k0[0] + std::pow(B[0],2)*k0[1] + std::pow(B[1],2)*k0[1] - 2*B[0]*k[0]*k0[1])/(2.*(B[0]*k[0] + B[1]*k[1] - B[0]*k0[0] - B[1]*k0[1]));

  s_0.index[0] = simplices[k_ind].index[0];
  s_0.index[1] = simplices.size()+1;

  if(k_ind == simplices[simplices[k_ind].index[0]].index[0])
    simplices[simplices[k_ind].index[0]].index[0] = simplices.size();
  else
    simplices[simplices[k_ind].index[0]].index[1] = simplices.size();

  s_1.k_vec[0] = -(-(std::pow(B[0],2)*k[0]) - std::pow(B[1],2)*k[0] + std::pow(B[0],2)*k1[0] + std::pow(B[1],2)*k1[0] - 2*B[1]*k[1]*k1[0] + 2*B[1]*k[0]*k1[1])/(2.*(B[0]*k[0] + B[1]*k[1] - B[0]*k1[0] - B[1]*k1[1]));
  s_1.k_vec[1] = -(-(std::pow(B[0],2)*k[1]) - std::pow(B[1],2)*k[1] + 2*B[0]*k[1]*k1[0] + std::pow(B[0],2)*k1[1] + std::pow(B[1],2)*k1[1] - 2*B[0]*k[0]*k1[1])/(2.*(B[0]*k[0] + B[1]*k[1] - B[0]*k1[0] - B[1]*k1[1]));

  s_1.index[0] = simplices[k_ind].index[1];
  s_1.index[1] = simplices.size();

  if(k_ind == simplices[simplices[k_ind].index[1]].index[0])
    simplices[simplices[k_ind].index[1]].index[0] = simplices.size()+1;
  else
    simplices[simplices[k_ind].index[1]].index[1] = simplices.size()+1;

  simplices.push_back(s_0);
  simplices.push_back(s_1);

  for(size_t l=0; l<simplices.size(); l++)
    {
      simplices[l].index[0] = simplices[l].index[0]>k_ind? simplices[l].index[0]-1: simplices[l].index[0];
      simplices[l].index[1] = simplices[l].index[1]>k_ind? simplices[l].index[1]-1: simplices[l].index[1];
    }

  simplices.erase(simplices.begin()+k_ind);
}


template<cluster_representation_type representation , typename parameters >
void tetrahedron_mesh<k_cluster<representation , cluster<parameters > > >::find_facets_2D()
{
  for(size_t l=0; l<simplices.size(); l++)
    {
      facet f0;
      f0.index.push_back(l);
      f0.index.push_back(simplices[l].index[0]);
      sort(f0.index.begin(), f0.index.end());
      
      facet f1;
      f1.index.push_back(l);
      f1.index.push_back(simplices[l].index[1]);
      sort(f1.index.begin(), f1.index.end());

      facets.push_back(f0);
      facets.push_back(f1);
    }

  for(size_t i=0; i<facets.size(); i++)
    {
      for(size_t j=i+1; j<facets.size(); j++)
	if(facets[i].index == facets[j].index)
	  facets.erase(facets.begin()+j);
    }
}

template<cluster_representation_type representation , typename parameters >
void tetrahedron_mesh<k_cluster<representation , cluster<parameters > > >::make_mesh_points_2D()
{
  const static int N=2;

  mesh.resize(1, std::vector<double>(DIMENSION, 0.));

  for(size_t l=0; l<simplices.size(); l++)
    mesh.push_back(simplices[l].k_vec);

  for(size_t l=0; l<facets.size(); l++)
    {
      tetrahedron<DIMENSION> tet;
      tet.index[0] = 0;
      tet.index[1] = facets[l].index[0]+1;
      tet.index[2] = facets[l].index[1]+1;

      tetrahedra.push_back(tet);
    }

  for(int i=0; i<N; i++)
    {
      int n_tet = tetrahedra.size();
      for(int l=0; l<n_tet; l++)
	tetrahedra[l].do_recursion(tetrahedra, mesh);

      for(int l=0; l<n_tet; l++)
	tetrahedra.erase(tetrahedra.begin());
    }

  for(size_t i=0; i<mesh.size(); i++){
    for(size_t j=i+1; j<mesh.size(); j++){
      if(L2_norm(mesh[i], mesh[j]) < 1.e-6){
	mesh.erase(mesh.begin()+j);
	
	for(size_t l=0; l<tetrahedra.size(); l++)
	  for(int z=0; z<DIMENSION+1; z++)
	    if(tetrahedra[l].index[z] == int(j))
	      tetrahedra[l].index[z] = int(i);
	    else
	      tetrahedra[l].index[z] = tetrahedra[l].index[z]>int(j) ? tetrahedra[l].index[z]-1:tetrahedra[l].index[z]; 
      }
    }
  }
}
