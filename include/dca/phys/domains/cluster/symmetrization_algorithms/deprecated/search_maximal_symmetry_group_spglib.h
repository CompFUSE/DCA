//-*-C++-*-

#ifndef SEARCH_MAX_SYMMETRY_GROUP_SPG_LIB_H
#define SEARCH_MAX_SYMMETRY_GROUP_SPG_LIB_H

template<class base_cluster_type, symmetry_group_level_type symmetry_group_level>
class search_maximal_symmetry_group<base_cluster_type, symmetry_package_spg_lib, symmetry_group_level>
{
  typedef point_group_symmetry_domain<symmetry_group_level, base_cluster_type> sym_dmn_t;

  typedef r_cluster<FULL, base_cluster_type> r_cluster_type;
  typedef k_cluster<FULL, base_cluster_type> k_cluster_type;

  typedef dmn_0<r_cluster_type> r_dmn_t;
  typedef dmn_0<k_cluster_type> k_dmn_t;

  const static int MAX_R_a = 1024;
  const static int MAX_OPS = 8192;

public:

  search_maximal_symmetry_group(int DIM, double b[3][3], double B[3][3], int Ns, double r_i[MAX_R_a][3], int Na, int flavors[MAX_R_a], double a_i[MAX_R_a][3]);
  ~search_maximal_symmetry_group();

  void initialize_symmetry_domain();

  void print_on_shell();

private:

  void invert_3D_matrix(double* lattice, double* lattice_inv);

  void generate_cartesian_cluster(int Ns,                         double r_i[MAX_R_a][3], 
				  int Na, int flavor_ai[MAX_R_a], double a_i[MAX_R_a][3]);
  void generate_affine_cluster(double lattice[3][3]);

  void generate_affine_symmetry_operations(double lattice[3][3]);
  void generate_cartesian_symmetry_operations(double lattice[3][3]);

  void generate_symmetry_information(double lattice[3][3]);

  void filter_2D_tranformations();

  void filter_cluster_translations(int Ns, double r_i[MAX_R_a][3]);

private:

  int    DIMENSION;

  int    N_cell;
  int    N_cell_symmetries;

  int    flavors[MAX_R_a];

  double ri_plus_ai_cart[MAX_R_a][3];
  double ri_plus_ai_frac[MAX_R_a][3];

  int    rotations_frac[MAX_OPS][3][3];
  double rotations_cart[MAX_OPS][3][3];

  double translations_frac[MAX_OPS][3];
  double translations_cart[MAX_OPS][3];
};

template<class base_cluster_type, symmetry_group_level_type symmetry_group_level>
search_maximal_symmetry_group<base_cluster_type, symmetry_package_spg_lib, symmetry_group_level>::search_maximal_symmetry_group(int    DIM,
																double b[3][3], 
																double B[3][3], 
																int Ns,                         double r_i[MAX_R_a][3], 
																int Na, int flavor_ai[MAX_R_a], double a_i[MAX_R_a][3]):
  DIMENSION(DIM),

  N_cell           (0),
  N_cell_symmetries(0)
{
  generate_cartesian_cluster(Ns, r_i, Na, flavor_ai, a_i);

  generate_affine_cluster(B);

  generate_symmetry_information(B);

  generate_affine_symmetry_operations(B);

  generate_cartesian_symmetry_operations(B);

  if(DIMENSION==2)
    filter_2D_tranformations();

  filter_cluster_translations(Ns, r_i);

  initialize_symmetry_domain();
}

template<class base_cluster_type, symmetry_group_level_type symmetry_group_level>
search_maximal_symmetry_group<base_cluster_type, symmetry_package_spg_lib, symmetry_group_level>::~search_maximal_symmetry_group()
{}

template<class base_cluster_type, symmetry_group_level_type symmetry_group_level>
void search_maximal_symmetry_group<base_cluster_type, symmetry_package_spg_lib, symmetry_group_level>::invert_3D_matrix(double* lattice, double* lattice_inv)
{
  invert_plan<double> invert_pln(3);
  
  memcpy(invert_pln.Matrix, lattice                   , sizeof(double)*3*3);
  invert_pln.execute_plan();
  memcpy(lattice_inv      , invert_pln.inverted_matrix, sizeof(double)*3*3);
}

template<class base_cluster_type, symmetry_group_level_type symmetry_group_level>
void search_maximal_symmetry_group<base_cluster_type, symmetry_package_spg_lib, symmetry_group_level>::generate_cartesian_cluster(int Ns,                         double r_i[MAX_R_a][3], 
																  int Na, int flavor_ai[MAX_R_a], double a_i[MAX_R_a][3])
{
  N_cell = Ns*Na;

  for(int i=0; i<Ns; ++i){
    for(int j=0; j<Na; ++j){
      
      int I = j+Na*i;
      
      flavors[I] = flavor_ai[j];
	
      for(int d=0; d<3; ++d)
	ri_plus_ai_cart[I][d] = r_i[i][d] + a_i[j][d];
    }
  }
}

template<class base_cluster_type, symmetry_group_level_type symmetry_group_level>
void search_maximal_symmetry_group<base_cluster_type, symmetry_package_spg_lib, symmetry_group_level>::generate_affine_cluster(double lattice[3][3])
{
  double T    [3][3];
  double T_inv[3][3];

  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      T[j][i] = lattice[i][j];

  invert_3D_matrix(&T[0][0], &T_inv[0][0]);

  double t[3];

  for(int l=0; l<N_cell; ++l){
      
    for(size_t i=0; i<3; ++i)
      t[i] = 0.;

    for(size_t i=0; i<3; ++i)
      for(size_t j=0; j<3; ++j)
	t[i] += T_inv[i][j]*ri_plus_ai_cart[l][j];
    
    for(size_t i=0; i<3; ++i){
      
      while(t[i]>(1.-1.e-6))
	t[i] -= 1.;
      
      while(t[i]<(-1.e-6))
	t[i] += 1.;
    }
    
    for(int d=0; d<3; ++d)
      ri_plus_ai_frac[l][d] = t[d];
  }
}

template<class base_cluster_type, symmetry_group_level_type symmetry_group_level>
void search_maximal_symmetry_group<base_cluster_type, symmetry_package_spg_lib, symmetry_group_level>::generate_affine_symmetry_operations(double lattice[3][3])
{
  N_cell_symmetries = SPG_LIB::spg_get_symmetry(rotations_frac,
						translations_frac,
						MAX_OPS,
						lattice,
						ri_plus_ai_frac,
						flavors,
						N_cell,
						1e-3);

//   if(false)
//     {
//       cout << "\n\n";
//       for(int i=0; i<N_cell_symmetries; i++) 
// 	{
// 	  printf("--- %d, ---\n", i+1);
	  
// 	  for (int j=0; j<3; j++)
// 	    cout << rotations_frac[i][j][0] <<"\t"<< rotations_frac[i][j][1] <<"\t"<< rotations_frac[i][j][2] <<"\n";
	  
// 	  cout << "\n\t";
// 	  cout << translations_frac[i][0] <<"\t"<< translations_frac[i][1] <<"\t"<< translations_frac[i][2]<<"\n";
// 	}
//       cout << "\n\n";
//     }
}

template<class base_cluster_type, symmetry_group_level_type symmetry_group_level>
void search_maximal_symmetry_group<base_cluster_type, symmetry_package_spg_lib, symmetry_group_level>::generate_symmetry_information(double lattice[3][3])
{
//   SPG_LIB::SpglibDataset* dataset;

//   dataset = SPG_LIB::spg_get_dataset( lattice,
// 				      ri_plus_ai_frac,
// 				      flavors,
// 				      N_cell,
// 				      1e-5 );
//   cout << "\n\n";
//   cout << "\t\t International : " <<  dataset->international_symbol << "\t" <<   dataset->spacegroup_number << "\n";
//   cout << "\t\t Hall symbol   : " <<  dataset->hall_symbol << "\n\n\n";

//   SPG_LIB::spg_free_dataset(dataset);
}

template<class base_cluster_type, symmetry_group_level_type symmetry_group_level>
void search_maximal_symmetry_group<base_cluster_type, symmetry_package_spg_lib, symmetry_group_level>::generate_cartesian_symmetry_operations(double lattice[3][3])
{
  double T    [3][3];
  double T_inv[3][3];

  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      T[j][i] = lattice[i][j];

  invert_3D_matrix(&T[0][0], &T_inv[0][0]);

  for(int l=0; l<N_cell_symmetries; ++l){
      
    for(size_t i=0; i<3; ++i){
      for(size_t j=0; j<3; ++j){
	rotations_cart[l][i][j] = 0.;
	
	for(size_t l0=0; l0<3; ++l0)
	  for(size_t l1=0; l1<3; ++l1)
	    rotations_cart[l][i][j] += T[i][l0]*rotations_frac[l][l0][l1]*T_inv[l1][j];
      }
    }
    
    for(size_t l0=0; l0<3; ++l0){
      
      std::vector<double> t(DIMENSION, 0);
      for(size_t l1=0; l1<3; ++l1)
	t[l1] += T[l0][l1]*translations_frac[l][l1];
    
      t = r_cluster_type::back_inside_lattice(t);

      for(size_t l1=0; l1<3; ++l1)
	translations_frac[l][l1] = t[l1];
    }
  }
}

template<class base_cluster_type, symmetry_group_level_type symmetry_group_level>
void search_maximal_symmetry_group<base_cluster_type, symmetry_package_spg_lib, symmetry_group_level>::filter_2D_tranformations()
{
  for(int l=0; l<N_cell_symmetries; ){
    
    bool is_a_3D_move = std::fabs(rotations_cart[l][2][2]-1.)>1.e-6 ? true: false;
    
    for(size_t l0=0; l0<2; ++l0)
      if(std::fabs(rotations_cart[l][l0][2]) > 1.e-6 || std::fabs(rotations_cart[l][2][l0]) > 1.e-6)
	is_a_3D_move=true;
    
    if( is_a_3D_move )
      {
	memmove(&rotations_cart   [l][0][0], &rotations_cart   [l+1][0][0], 9*sizeof(double)*(MAX_OPS-l-1));      
	memmove(&translations_cart[l][0]   , &translations_cart[l+1][0]   , 3*sizeof(double)*(MAX_OPS-l-1));
	
	N_cell_symmetries -= 1;
      }
    else
      l++;
  }
}

template<class base_cluster_type, symmetry_group_level_type symmetry_group_level>
void search_maximal_symmetry_group<base_cluster_type, symmetry_package_spg_lib, symmetry_group_level>::filter_cluster_translations(int Ns, double r_i[MAX_R_a][3])
{
  for(int i=0; i<N_cell_symmetries; i++){

    for(int j=i+1; j<N_cell_symmetries; ){

      bool is_the_same = true;

      for(int d0=0; d0<3; d0++)
        if(std::fabs(translations_cart[j][d0]-translations_cart[i][d0]) > 1.e-6)
	  is_the_same = false;

      for(int d0=0; d0<3; d0++)
	for(int d1=0; d1<3; d1++)
	  if(std::fabs(rotations_cart[j][d0][d1]-rotations_cart[i][d0][d1]) > 1.e-6)
	    is_the_same = false;

      if(is_the_same)
	{
	  N_cell_symmetries -= 1;

	  memmove(&rotations_cart   [j][0][0], &rotations_cart   [j+1][0][0], 9*sizeof(double)*(MAX_OPS-j-1));      
	  memmove(&translations_cart[j][0]   , &translations_cart[j+1][0]   , 3*sizeof(double)*(MAX_OPS-j-1));
	}
      else
	j++;
    }
  }
}

template<class base_cluster_type, symmetry_group_level_type symmetry_group_level>
void search_maximal_symmetry_group<base_cluster_type, symmetry_package_spg_lib, symmetry_group_level>::print_on_shell()
{
  for(int i=0; i<N_cell_symmetries; i++) 
    {
      printf("--- %d, ---\n", i+1);

      for (int j=0; j<3; j++)
	printf("%4f %4f %4f\n", rotations_cart[i][j][0], rotations_cart[i][j][1], rotations_cart[i][j][2]);
      
      cout << "\n\t";
      printf("%f %f %f\n", translations_cart[i][0], translations_cart[i][1], translations_cart[i][2]);
    }
}

template<class base_cluster_type, symmetry_group_level_type symmetry_group_level>
void search_maximal_symmetry_group<base_cluster_type, symmetry_package_spg_lib, symmetry_group_level>::initialize_symmetry_domain()
{
  sym_dmn_t::DIMENSION = DIMENSION;
      
  sym_dmn_t::get_size() = N_cell_symmetries;
      
  for(int l=0; l<sym_dmn_t::get_size(); ++l){
    
    for(int i=0; i<DIMENSION; ++i)
      for(int j=0; j<DIMENSION; ++j)
	sym_dmn_t::get_elements()[l].O[i+DIMENSION*j] = rotations_cart[l][i][j];
    
    for(int i=0; i<DIMENSION; ++i)
      sym_dmn_t::get_elements()[l].t[i] = translations_cart[l][i];
  }
}

#endif
