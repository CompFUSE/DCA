//-*-C++-*-

#ifndef TETRAHEDRON_SECOND_NEAREST_NEIGHBOUR_DOMAIN_H
#define TETRAHEDRON_SECOND_NEAREST_NEIGHBOUR_DOMAIN_H


/*! 
 *  \ingroup TETRAHEDRON-MESH
 *
 *  \author Peter Staar
 *  \brief  This class constructs a tetrahedron second nearest neighbour domain in the Brillouin-zone, defined by the k-cluster template.
 */
template<cluster_representation_type representation, class parameters>
class tetrahedron_second_nearest_neighbour_domain<k_cluster<representation, cluster<parameters> > >
{
public:
  
  typedef tetrahedron_second_nearest_neighbour_domain this_type;
  typedef std::vector<double>          element_type;
  
  typedef k_cluster<representation , cluster<parameters > > k_cluster_type;
  typedef cluster<parameters >                              base_cluster_type;

  typedef dmn_0<k_cluster_type> k_dmn_t;

  typedef tetrahedron_neighbour_domain<k_cluster_type> tetrahedron_frist_nn_dmn_t;

  const static int DIMENSION = k_cluster_type::DIMENSION;

  static int get_size(){
    return get_elements().size();
  }

  static std::vector<element_type>& get_elements()
  {
    static std::vector<element_type>& elements = initialize();
    return elements;
  }

  static int get_K_index(int K_ind, int n_ind)
  {
    static function<int, dmn_2<k_dmn_t, dmn_0<this_type> > >& f = initialize_function();
    return f(K_ind, n_ind);
  }
  
private:

  static std::vector<element_type>& initialize()
  {
    static std::vector<element_type> elements(0);

    for(size_t f_ind=0; f_ind<tet_mesh.get_facets().size(); ++f_ind){
      
      std::vector<double> k(DIMENSION,0.);

      std::vector<double> k(DIMENSION,0.);
      
      for(size_t k_ind=0; k_ind<tet_mesh.get_facets()[f_ind].index.size(); ++k_ind)
	k = VECTOR_OPERATIONS::ADD(k, tet_mesh.get_simplices()[tet_mesh.get_facets()[f_ind].index[k_ind]].k_vec);

      elements.push_back(k);
    }

    return elements;
  }

  static function<int, dmn_2<k_dmn_t, dmn_0<this_type> > >& initialize_function()
  {
    std::vector<element_type>& k_vecs = k_dmn_t::get_elements();
    std::vector<element_type>& n_vecs = this_type::get_elements();
    
    static function<int, dmn_2<k_dmn_t, dmn_0<this_type> > > f;

    for(int K_ind=0; K_ind<k_dmn_t::dmn_size(); ++K_ind){
      for(int n_ind=0; n_ind<this_type::get_size(); ++n_ind){
      
	std::vector<double> k(DIMENSION,0.);
	
	k = VECTOR_OPERATIONS::ADD(k_vecs[K_ind], n_vecs[n_ind]);
	
	k = k_cluster_type::back_inside_cluster(k);
	
	int index=-1;
	for(int l=0; l<k_dmn_t::dmn_size(); ++l){
	  if(VECTOR_OPERATIONS::L2_NORM(k, k_vecs[l])<1.e-6){
	    index=l;
	    break;
	  }
	}
	
	assert(index>-1);
	f(K_ind, n_ind) = index;
      }
    }
    
    return f;
  }

};

#endif
