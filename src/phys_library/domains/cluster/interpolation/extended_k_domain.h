//-*-C++-*-

#ifndef EXTENDED_K_DOMAIN_H
#define EXTENDED_K_DOMAIN_H

/*!
 *  \ingroup INTERPOLATION
 *
 *  \author  Peter Staar
 */
template<typename k_cluster_type>
class extended_k_domain
{
  const static int N_delta = 3;

  const static int DIMENSION = k_cluster_type::DIMENSION;  
  //typedef tetrahedron_neighbour_domain<k_cluster_type> tet_dmn_t;

public:

  typedef extended_k_domain<k_cluster_type> this_type;
  typedef std::vector<double>          element_type;

public:

  static int get_size(){
    return get_elements().size();
  }

  static std::vector<element_type>& get_elements()
  {
    static std::vector<element_type> elements = initialize();
    return elements;
  }

  static std::vector<element_type>& get_sorted_elements()
  {
    static std::vector<element_type> elements = initialize_sorted();
    return elements;
  }

  inline static int get_index(std::vector<double>& k){
    static std::vector<std::vector<double> > vecs = get_sorted_elements();
    
    static std::vector<double> k_vec(DIMENSION+1, 0);

    for(int d=0; d<DIMENSION; ++d)
      k_vec[d] = k[d];
    k_vec[DIMENSION] = -1;
    
    int index = int(std::lower_bound(vecs.begin(), vecs.end(), k_vec, VECTOR_OPERATIONS::IS_LARGER_VECTOR<double>) - vecs.begin());

    index = vecs[index][DIMENSION];

    assert(check_consistency(index, k));

    return index;
  }

  inline static bool check_consistency(int index, std::vector<double>& k){
    
    if(index<0 or index>=get_size())
      throw std::logic_error(__FUNCTION__); 
    
    if(VECTOR_OPERATIONS::L2_NORM(k, get_elements()[index])>1.e-6){
      std::cout << "\n\n";
      VECTOR_OPERATIONS::PRINT(k); 
      std::cout << "\t<-->\t"; 
      VECTOR_OPERATIONS::PRINT(get_elements()[index]);
      std::cout << "\n\n";
      throw std::logic_error(__FUNCTION__); 
    }
  
    return true;
  }

  inline static int get_cluster_index(int index){

    static int Nc = k_cluster_type::get_size();

    assert(index>-1 and index<get_size());
    
    return (index % Nc);
  }

private:

  static std::vector<element_type> initialize_sorted(){
    
    std::vector<element_type> sorted_elements = get_elements();
  
    for(size_t l=0; l<sorted_elements.size(); l++)
      sorted_elements[l].push_back(l);
    
    std::sort(sorted_elements.begin(), sorted_elements.end(), VECTOR_OPERATIONS::IS_LARGER_VECTOR<double>);
  
    return sorted_elements;
  }

  static std::vector<element_type> initialize(){

    std::vector<element_type> elements(0);

    switch(DIMENSION)
      {
      case 1:
	elements = initialize_1D();
	break;

      case 2:
	elements = initialize_2D();
	break;

      case 3:
	elements = initialize_3D();
	break;

      default:
	throw std::logic_error(__FUNCTION__);
      }

    return elements;
  }

  static std::vector<element_type> initialize_1D(){
    std::vector<element_type> elements(0);

    element_type K_vec(DIMENSION,0);

    for(int b0=-N_delta; b0<=N_delta; ++b0){
      
      for(int K_ind=0; K_ind<k_cluster_type::get_size(); ++K_ind){
	  
	K_vec = k_cluster_type::get_elements()[K_ind];
	
	for(int d=0; d<DIMENSION; ++d)
	  K_vec[d] += (b0*k_cluster_type::get_super_basis_vectors()[0][d]);
	
	elements.push_back(K_vec);
      }
    }

    return elements;
  }

  static std::vector<element_type> initialize_2D(){
    std::vector<element_type> elements(0);

    element_type K_vec(DIMENSION,0);

    for(int b0=-N_delta; b0<=N_delta; ++b0){
      for(int b1=-N_delta; b1<=N_delta; ++b1){

	for(int K_ind=0; K_ind<k_cluster_type::get_size(); ++K_ind){
	  
	  K_vec = k_cluster_type::get_elements()[K_ind];
	  
	  for(int d=0; d<DIMENSION; ++d)
	    K_vec[d] += (b0*k_cluster_type::get_super_basis_vectors()[0][d] 
			 + b1*k_cluster_type::get_super_basis_vectors()[1][d]);

	  elements.push_back(K_vec);
	}
      }
    }

    return elements;
  }

  static std::vector<element_type> initialize_3D(){
    std::vector<element_type> elements(0);

    element_type K_vec(DIMENSION,0);

    for(int b0=-N_delta; b0<=N_delta; ++b0){
      for(int b1=-N_delta; b1<=N_delta; ++b1){
	for(int b2=-N_delta; b2<=N_delta; ++b2){

	  for(int K_ind=0; K_ind<k_cluster_type::get_size(); ++K_ind){
	  
	    K_vec = k_cluster_type::get_elements()[K_ind];
	    
	    for(int d=0; d<DIMENSION; ++d)
	      K_vec[d] += (b0*k_cluster_type::get_super_basis_vectors()[0][d] 
			   + b1*k_cluster_type::get_super_basis_vectors()[1][d]
			   + b2*k_cluster_type::get_super_basis_vectors()[2][d]);

	    elements.push_back(K_vec);
	  }
	}
      }
    }

    return elements;
  }
};

#endif
