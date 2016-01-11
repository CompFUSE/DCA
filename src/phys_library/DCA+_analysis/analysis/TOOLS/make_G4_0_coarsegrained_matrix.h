//-*-C++-*-

/** \ingroup DCA */

/*@{*/

/*! \file MAKE_G4_0_coarsegrained_MATRIX.h
 *
 * author : peter staar
 */

#ifndef MAKE_G4_0_COARSEGRAINED_MATRIX_H
#define MAKE_G4_0_COARSEGRAINED_MATRIX_H

// namespace dca {

//   template<class parameter_type, class MOMS_type>
//   class make_G4_0_coarsegrained_matrix 
//   {
// #include "type_definitions.h"
    
//   public:
 
//     make_G4_0_coarsegrained_matrix(parameter_type& parameters_ref, MOMS_type& MOMS_ref);
//     ~make_G4_0_coarsegrained_matrix();

//     template<typename analysis_t>
//     void execute(analysis_t& analysis_ref);

//   private:

//     void initialize();

//     void compute_G_k(int K_ind, int w_ind);

//     void compute_G_k_accent(int K_accent_ind, std::pair<int, double>& w_accent);

//     template<typename analysis_t>
//     void compute_bubble(analysis_t& analysis_ref, int K_ind, int w_vertex_index);

//     double                               get_integration_factor();
//     std::pair<int, std::vector<double> > get_k_accent(int K_ind);
//     std::pair<int,             double  > get_w_accent(int w_index, int W);

//   private:    
    
//     parameter_type&  parameters;
//     MOMS_type&       MOMS;

//     b_k_DCA_w_VERTEX_EXTENDED b_k_DCA_w_VERTEX_EXTENDED_domain;
//     b_k_DCA_w_VERTEX          b_k_DCA_w_VERTEX_domain;

//     std::vector<int>  corresponding_extended_index;
//     std::vector<bool> is_compact_frequency;

//     int matrix_size, matrix_dim;

//     std::complex<double>* tmp_matrix;
//     std::complex<double>* Sigma_matrix;
//     std::complex<double>* H_LDA_matrix_k;
//     std::complex<double>* H_LDA_matrix_k_accent;
//     std::complex<double>* G_k;
//     std::complex<double>* G_k_accent;
//   };

//   template<class parameter_type, class MOMS_type>
//   make_G4_0_coarsegrained_matrix<parameter_type, MOMS_type>::make_G4_0_coarsegrained_matrix (parameter_type& parameters_in, 
// 											     MOMS_type& MOMS_in):
//     parameters(parameters_in),
//     MOMS(MOMS_in),

//     corresponding_extended_index(w_VERTEX::dmn_size(),-1),
//     is_compact_frequency(w_VERTEX_EXTENDED::dmn_size(),false)
//   {
//     initialize();

//     tmp_matrix            = new std::complex<double>[matrix_size];
//     Sigma_matrix          = new std::complex<double>[matrix_size];
//     H_LDA_matrix_k        = new std::complex<double>[matrix_size];
//     H_LDA_matrix_k_accent = new std::complex<double>[matrix_size];
//     G_k                   = new std::complex<double>[matrix_size];
//     G_k_accent            = new std::complex<double>[matrix_size];
//   }

//   template<class parameter_type, class MOMS_type>
//   make_G4_0_coarsegrained_matrix<parameter_type, MOMS_type>::~make_G4_0_coarsegrained_matrix()
//   {
//     delete [] tmp_matrix            ;
//     delete [] Sigma_matrix          ;
//     delete [] H_LDA_matrix_k        ;
//     delete [] H_LDA_matrix_k_accent ;
//     delete [] G_k                   ;
//     delete [] G_k_accent            ;
//   }

//   template<class parameter_type, class MOMS_type>
//   void make_G4_0_coarsegrained_matrix<parameter_type, MOMS_type>::initialize()
//   {
//     matrix_dim  = 2*b::dmn_size();
//     matrix_size = square(matrix_dim);

//     for(int i=0; i<w_VERTEX::dmn_size(); i++)
//       for(int j=0; j<w_VERTEX_EXTENDED::dmn_size(); j++)
// 	if(std::fabs(w_VERTEX::parameter_type::get_elements()[i]-w_VERTEX_EXTENDED::parameter_type::get_elements()[j])<1.e-6)
// 	  corresponding_extended_index[i] = j;

//     for(int j=0; j<w_VERTEX_EXTENDED::dmn_size(); j++)
//       for(int i=0; i<w_VERTEX::dmn_size(); i++)
//       	if(std::fabs(w_VERTEX::parameter_type::get_elements()[i]-w_VERTEX_EXTENDED::parameter_type::get_elements()[j])<1.e-6)
// 	  is_compact_frequency[j] = true;	  
//   }

//   template<class parameter_type, class MOMS_type>
//   template<typename analysis_t>
//   void make_G4_0_coarsegrained_matrix<parameter_type, MOMS_type>::execute(analysis_t& analysis_ref)
//   {
//     analysis_ref.full_chi_0_b_k_w__b_k_w = 0.;

//     cout << "Nb_integration_points_per_DCA_cell --> 20" << endl;

//     double                            Nb_interpolation = 20;//parameters.get_Nb_integration_points_per_DCA_cell()[0];
//     std::vector<std::vector<double> > centered_mesh    = Mesh<DCA_k_cluster_type>::execute(int(Nb_interpolation));
//     std::vector<std::vector<double> > mesh             = centered_mesh;
    
//     wannier_interpolation::mesh_k::get_size() = int(mesh.size());
//     wannier_interpolation WIP_k       (MOMS.H_LDA);
//     wannier_interpolation WIP_k_accent(MOMS.H_LDA);

//     invert_plan<std::complex<double> > invert_pln(matrix_dim);

//     for(int K_ind=0; K_ind<DCA_k_cluster_type::get_size(); K_ind++)
//       {
// 	std::vector<double> K = DCA_k_cluster_type::get_elements()[K_ind];
// 	Mesh<DCA_k_cluster_type>::translate_mesh(centered_mesh, mesh, K);
// 	wannier_interpolation::H_k_interpolated_type& H_k        = WIP_k.execute(mesh);

// 	int                 K_accent_ind = get_k_accent(K_ind).first;
// 	std::vector<double> K_accent     = get_k_accent(K_ind).second;
// 	Mesh<DCA_k_cluster_type>::translate_mesh(centered_mesh, mesh, K_accent);
// 	wannier_interpolation::H_k_interpolated_type& H_k_accent = WIP_k_accent.execute(mesh);

// 	// \sum_{k \in K_{k}}
// 	for(int k_ind=0; k_ind<int(mesh.size()); k_ind++)
// 	  {
// 	    memcpy(H_LDA_matrix_k       , &H_k       (0,0,k_ind), sizeof(std::complex<double>)*matrix_size);
// 	    memcpy(H_LDA_matrix_k_accent, &H_k_accent(0,0,k_ind), sizeof(std::complex<double>)*matrix_size);

// 	    for(int w_vertex_index=0; w_vertex_index<w_VERTEX_EXTENDED::dmn_size(); w_vertex_index++)
// 	      {
// // 		double w  = w_VERTEX_EXTENDED::parameter_type::get_elements()[w_vertex_index];
// // 		assert( std::fabs(w -  w::parameter_type::get_elements()[w_ind]) < 1.e-6);

// 		int w_ind = w_vertex_index - w_VERTEX_EXTENDED::dmn_size()/2 + w::dmn_size()/2;

// 		compute_G_k(K_ind, w_vertex_index);

// 		// G(k+Q) --> particle-hole
// 		// G(Q-k) --> particle-particle
// 		{
// 		  std::pair<int, double> w_accent = get_w_accent(w_ind, parameters.get_w_channel());
// 		  compute_G_k_accent(K_accent_ind, w_accent);
// 		}

// 		// G(k)*G(k+Q) --> particle-hole
// 		// G(k)*G(Q-k) --> particle-particle
// 		compute_bubble(analysis_ref, K_ind, w_vertex_index);
// 	      }
// 	  }
//       }

//     double integration_factor = get_integration_factor()/double(centered_mesh.size());
//     analysis_ref.full_chi_0_b_k_w__b_k_w *= integration_factor;
//   }

//   template<class parameter_type, class MOMS_type>
//   void make_G4_0_coarsegrained_matrix<parameter_type, MOMS_type>::compute_G_k(int K_ind, int w_vertex_index)
//   {
//     static invert_plan<std::complex<double> > invert_pln(matrix_dim);

//     double w  = w_VERTEX_EXTENDED::parameter_type::get_elements()[w_vertex_index];
//     int w_ind = w_vertex_index - w_VERTEX_EXTENDED::dmn_size()/2 + w::dmn_size()/2;

//     memcpy(Sigma_matrix, &MOMS.Sigma(0,0,K_ind, w_ind), sizeof(std::complex<double>)*matrix_size);
    
//     {// G^-1 =  -(H_k + Sigma) + i*w + mu
//       for(int index=0; index < matrix_size; index++)
// 	tmp_matrix[index] = -(H_LDA_matrix_k[index] + Sigma_matrix[index]);  
      
//       for(int nu=0; nu<matrix_dim; nu++)
// 	tmp_matrix[nu + matrix_dim*nu] += std::complex<double>(parameters.get_chemical_potential(), w);
//     }
		  
//     {// G(k)^-1 = (-(H(k)+Sigma(k)) + i*w) ==> G(k) = (-(H(k)+Sigma(k)) + i*w)^-1 
//       memcpy(invert_pln.Matrix, tmp_matrix                    , sizeof(std::complex<double>)*matrix_size);
//       invert_pln.execute_plan();
//       memcpy(&G_k[0]          , &invert_pln.inverted_matrix[0], sizeof(std::complex<double>)*matrix_size);
//     }
//   }

//   template<class parameter_type, class MOMS_type>
//   void make_G4_0_coarsegrained_matrix<parameter_type, MOMS_type>::compute_G_k_accent(int K_accent_ind, std::pair<int, double>& w_accent)
//   {
//     static invert_plan<std::complex<double> > invert_pln(matrix_dim);
    
//     memcpy(Sigma_matrix, &MOMS.Sigma(0,0,K_accent_ind, w_accent.first), sizeof(std::complex<double>)*matrix_size);
    
//     {// G^-1 =  -(H_k_accent + Sigma) + i*w + mu
//       for(int index=0; index < matrix_size; index++)
// 	tmp_matrix[index] = -(H_LDA_matrix_k_accent[index] + Sigma_matrix[index]);  
      
//       for(int nu=0; nu<matrix_dim; nu++)
// 	tmp_matrix[nu + matrix_dim*nu] += std::complex<double>(parameters.get_chemical_potential(), w_accent.second);
//     }
    
//     {// G(k+Q)^-1 = (-(H(k+Q)+Sigma(k+Q)) + i*w) ==> G(k+Q) = (-(H(k+Q)+Sigma(k+Q)) + i*w)^-1 
//       memcpy(invert_pln.Matrix, tmp_matrix                     , sizeof(std::complex<double>)*matrix_size);
//       invert_pln.execute_plan();
//       memcpy(&G_k_accent[0]    , &invert_pln.inverted_matrix[0], sizeof(std::complex<double>)*matrix_size);
//     }
//   }

//   template<class parameter_type, class MOMS_type>
//   template<typename analysis_t>
//   void make_G4_0_coarsegrained_matrix<parameter_type, MOMS_type>::compute_bubble(analysis_t& analysis_ref, int K_ind, int w_vertex_index)
//   {
//     for(int b1=0; b1<b::dmn_size(); b1++)
//       for(int b2=0; b2<b::dmn_size(); b2++)
// 	for(int b3=0; b3<b::dmn_size(); b3++)
// 	  analysis_ref.full_chi_0_b_k_w__b_k_w(b1, K_ind, w_vertex_index, b2, K_ind, w_vertex_index) += G_k[b1+b3*matrix_dim] * G_k_accent[b3+b2*matrix_dim];
//   }

//   template<class parameter_type, class MOMS_type>
//   double make_G4_0_coarsegrained_matrix<parameter_type, MOMS_type>::get_integration_factor()
//   {
//     switch(parameters.get_vertex_measurement_type())
//       {
//       case PARTICLE_HOLE_MAGNETIC:
// 	return -1.;
// 	break;
	
//       case PARTICLE_HOLE_CHARGE:
// 	return -2.;
// 	break;
	
//       case PARTICLE_PARTICLE_SUPERCONDUCTING:
// 	return 1;
// 	break;
	
//       default:
// 	throw std::logic_error(__FUNCTION__);
//       }
//   }

//   template<class parameter_type, class MOMS_type>
//   std::pair<int, std::vector<double> > make_G4_0_coarsegrained_matrix<parameter_type, MOMS_type>::get_k_accent(int K_ind)
//   {
//     int Q_ind             = parameters.get_q_channel();
//     std::vector<double> Q = parameters.get_q_vector();

//     std::pair<int, std::vector<double> > k_accent;

//     switch(parameters.get_vertex_measurement_type())
//       {
//       case PARTICLE_HOLE_MAGNETIC:
// 	{
// 	  int index                    = DCA_k_cluster_type::add(K_ind, Q_ind);
// 	  std::vector<double> K_plus_Q = DCA_k_cluster_type::get_elements()[K_ind];

// 	  for(int i=0; i<int(K_plus_Q.size()); i++)
// 	    K_plus_Q[i] += Q[i];

// 	  k_accent.first  = index;
// 	  k_accent.second = K_plus_Q;
// 	}
// 	break;
	
//       case PARTICLE_HOLE_CHARGE:
// 	{
// 	  int index                    = DCA_k_cluster_type::add(K_ind, Q_ind);
// 	  std::vector<double> K_plus_Q = DCA_k_cluster_type::get_elements()[K_ind];

// 	  for(int i=0; i<int(K_plus_Q.size()); i++)
// 	    K_plus_Q[i] += Q[i];

// 	  k_accent.first  = index;
// 	  k_accent.second = K_plus_Q;
// 	}
// 	break;
	
//       case PARTICLE_PARTICLE_SUPERCONDUCTING:
// 	{
// 	  int index                   = DCA_k_cluster_type::subtract(K_ind, Q_ind);
// 	  std::vector<double> Q_min_K = Q;

// 	  for(int i=0; i<int(Q_min_K.size()); i++)
// 	    Q_min_K[i] -= DCA_k_cluster_type::get_elements()[K_ind][i];

// 	  k_accent.first  = index;
// 	  k_accent.second = Q_min_K;
// 	}
// 	break;
	
//       default:
// 	throw std::logic_error(__FUNCTION__);
//       }

//     return k_accent;
//   }

//   template<class parameter_type, class MOMS_type>
//   std::pair<int, double> make_G4_0_coarsegrained_matrix<parameter_type, MOMS_type>::get_w_accent(int w_index, int W)
//   {
//     std::pair<int, double> w_accent;

//     switch(parameters.get_vertex_measurement_type())
//       {
//       case PARTICLE_HOLE_MAGNETIC:
// 	w_accent.first = w_index + W;
// 	assert(w_accent.first >= 0 && w_accent.first < frequency_domain_type::get_size());
// 	w_accent.second = frequency_domain_type::get_elements()[w_accent.first];
// 	break;
	
//       case PARTICLE_HOLE_CHARGE:
// 	w_accent.first = w_index + W;
// 	assert(w_accent.first >= 0 && w_accent.first < frequency_domain_type::get_size());
// 	w_accent.second = frequency_domain_type::get_elements()[w_accent.first];
// 	break;
	
//       case PARTICLE_PARTICLE_SUPERCONDUCTING:
// 	w_accent.first = W + (frequency_domain_type::get_size()-1)-w_index ;
// 	assert(w_accent.first >= 0 && w_accent.first < frequency_domain_type::get_size());
// 	w_accent.second = frequency_domain_type::get_elements()[w_accent.first];
// 	break;
	
//       default:
// 	throw std::logic_error(__FUNCTION__);
//       }

//     return w_accent;
//   }
// }

#endif
