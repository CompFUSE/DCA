//-*-C++-*-

#ifndef DIAGRAMMATIC_SYMMETRIES_H
#define DIAGRAMMATIC_SYMMETRIES_H

/*! 
 *  \ingroup SYMMETRIZE-FUNCTIONS
 *
 *  \author Peter Staar
 *  \brief  This class implements symmetrizations of the 2-particle functions according to the diagrams.
 *
 *  \image html feynman_diagram.pdf
 *
 *
 *  \section all symmetries for all 2-particle quantities.
 *  
 *   All two-particle functions have translation symmetry,
 *
 *   \f{eqnarray*}{
 *     G^{2}(k_1,\varpi_1, k_2,\varpi_1, q,\nu) &=&  G^{2}(k_1+P,\varpi_1+\mu, k_2+P,\varpi_2+\mu, q+P, \nu+\mu)
 *   \f}
 *
 *   Since all two-particle functions are real in real-space and imaginary time, we have an inversion symmetry,
 *
 *   \f{eqnarray*}{
 *     G^{2}(k_1,\varpi_1, k_2,\varpi_1, q,\nu) &=& \overline{G^{2}(-k_1,-\varpi_1, -k_2,-\varpi_2, -q,-\nu)}
 *   \f}
 *
 *   If the crystal has inversion symmetry, then
 *
 *   \f{eqnarray*}{
 *     G^{2}(k_1,\varpi_1, k_2,\varpi_1, q,\nu) &=& \overline{G^{2}(k_1,-\varpi_1, k_2,-\varpi_1, q,-\nu)}
 *   \f}
 *
 *  \section ph particle-hole channel
 *
 *   The particle-hole Greens function \f$G^{ph}\f$ flips the arrow direction under a rotation of \f$\pi\f$ around the horizontal axes (see fig 1).
 *   We therefore get that,
 *
 *   \f{eqnarray*}{
 *   G^{pp}(k_1, k_2, q) &=& G^{ph}(k_1+q, k_1, k_2, k_2+q) \\ 
 *                       &=& \overline{G^{ph}(k_1, k_1+q, k_2+q, k_2)} \\ 
 *			 &=& \overline{G^{ph}(k_1+q, k_2+q, -q)}
 *   \f}
 *
 *   The particle-hole Greens function \f$G^{ph}\f$ flips the arrow direction under a rotation of \f$\pi\f$ around the verical axes (see fig 1).
 *   We therefore get that,
 *
 *   \f{eqnarray*}{
 *   G^{pp}(k_1, k_2, q) &=& G^{ph}(k_1+q, k_1, k_2, k_2+q) \\ 
 *                       &=& \overline{G^{ph}(k_2+q, k_2, k_1, k_1+q)} \\ 
 *                       &=& \overline{G^{ph}(k_2, k_1, q)}
 *   \f}
 *
 *   The particle-hole Greens function \f$G^{ph}\f$ remains invariant under a rotation of \f$\pi\f$ around both the horizontal and vertical axes (see fig 1).
 *   We therefore get that,
 * \f{eqnarray*}{
 *   G^{ph}(k_1, k_2, q) &=& G^{ph}(k_1+q, k_1, k_2, k_2+q) \\ 
 *                       &=& G^{ph}(k_1, k_1+q, k_2+q, k_2) \\ 
 *                       &=& G^{ph}(k_2+q, k_1+q, -q)  
 * \f}
 *
 *
 *
 *  \section pp particle-particle channel
 *
 *   The particle-particle Greens function \f$G^{pp}\f$ remains invariant under a rotation of \f$\pi\f$ around the horizontal axes (see fig 1).
 *   We therefore get that,
 *
 *   \f{eqnarray*}{
 *   G^{pp}(k_1, k_2, q) &=& G^{pp}(q-k_1, k_1, k_2, q-k_2) \\ 
 *                       &=& G^{pp}(k_1, q-k_1, q-k_2, k_2) \\ 
 *                       &=& G^{pp}(q-k_1, q-k_2, q)  
 *   \f}
 *
 *   The particle-particle Greens function \f$G^{pp}\f$ remains invariant under a rotation of \f$\pi\f$ around the vertical axes (see fig 1).
 *   We therefore get that,
 *
 *   \f{eqnarray*}{
 *   G^{pp}(k_1, k_2, q) &=& G^{pp}(q-k_1, k_1, k_2, q-k_2) \\ 
 *                       &=& \overline{G^{pp}(q-k_2, k_2, k_1, q-k_1)} \\ 
 *                       &=& \overline{G^{pp}(k_2, k_1, q)}
 *   \f}
 *
 *   The particle-particle Greens function \f$G^{pp}\f$ remains invariant under a rotation of \f$\pi\f$ around the horizontal and vertical axes (see fig 1).
 *   We therefore get that,
 *
 *   \f{eqnarray*}{
 *   G^{pp}(k_1, k_2, q) &=& G^{pp}(q-k_1, k_1, k_2, q-k_2) \\ 
 *                       &=& \overline{G^{pp}(q-k_2, k_2, k_1, q-k_1)} \\ 
 *                       &=& \overline{G^{pp}(q-k_2, q-k_1, q)}
 *   \f}
 *
 *  Furthermore, \f$G^{pp}\f$ is real for any q-vector, since
 *
 *   \f{eqnarray*}{
 *   G^{pp}(k_1, k_2, q) &=& G^{pp}(q-k_1, k_1, k_2, q-k_2) \\ 
 *                       &=& G^{pp}(k_1, q-k_1, q-k_2, k_2) \\ 
 *                       &=& \overline{G^{pp}(-k_1, k_1-q, k_2-q,-k_2)}  \ \
 *			 &=& \overline{G^{pp}(q-k_1, k_1, k_2,q-k_2)} \	\
 *			 &=& \overline{G^{pp}(k_1, k_2, q)}
 *   \f}
 *
 * 
 */
template<class parameters_type>
class diagrammatic_symmetries
{
#include "type_definitions.h"

public:

  diagrammatic_symmetries(parameters_type& parameters);
  ~diagrammatic_symmetries();

  template<typename scalartype, typename k_dmn, typename w_dmn>
  void execute(function<scalartype, dmn_2<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn> > >& G);

  template<typename scalartype, typename k_dmn, typename w_dmn>
  void execute(function<scalartype, dmn_3<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn>, k_dmn > >& G);

private:

  template<typename scalartype, typename k_dmn, typename w_dmn>
  void symmetrize_over_matsubara_frequencies(function<scalartype, dmn_2<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn> > >& G);

  template<typename scalartype, typename k_dmn, typename w_dmn>
  void symmetrize_over_matsubara_frequencies(function<scalartype, dmn_3<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn>, k_dmn> >& G);



  template<typename scalartype, typename k_dmn, typename w_dmn>
  void symmetrize_over_pi_rotations_ph(function<scalartype, dmn_2<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn> > >& G);

  template<typename scalartype, typename k_dmn, typename w_dmn>
  void symmetrize_over_pi_rotations_ph(function<scalartype, dmn_3<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn>, k_dmn> >& G);



  template<typename scalartype, typename k_dmn, typename w_dmn>
  void symmetrize_over_pi_rotations_pp(function<scalartype, dmn_2<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn> > >& G);
  
  template<typename scalartype, typename k_dmn, typename w_dmn>
  void symmetrize_over_pi_rotations_pp(function<scalartype, dmn_3<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn>, k_dmn> >& G);
  


  template<typename scalartype, typename k_dmn, typename w_dmn>
  void set_real(function<scalartype, dmn_2<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn> > >& G);

  template<typename scalartype, typename k_dmn, typename w_dmn>
  void set_real(function<scalartype, dmn_3<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn>, k_dmn> >& G);

private:

  parameters_type& parameters;

  int                 q_ind;
  std::vector<double> q_vec;

  bool q_vector_is_invertible;
  bool q_vector_is_reciprocal;
};

template<class parameters_type>
diagrammatic_symmetries<parameters_type>::diagrammatic_symmetries(parameters_type& parameters_ref):
  parameters(parameters_ref),

  q_ind(parameters.get_q_channel()),
  q_vec(parameters.get_q_vector()),

  q_vector_is_invertible(false),
  q_vector_is_reciprocal(false)
{
}

template<class parameters_type>
diagrammatic_symmetries<parameters_type>::~diagrammatic_symmetries()
{}

template<class parameters_type>
template<typename scalartype, typename k_dmn, typename w_dmn>
void diagrammatic_symmetries<parameters_type>::execute(function<scalartype, dmn_2<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn> > >& G)
{  
  /*
   *  if you have more bands, make sure that the inversion does not interfere with a permutation of the bands !!!
   */
  if(true)
  {
    std::vector<double> q_rec(q_vec);
    
    for(size_t l=0; l<q_rec.size(); ++l)
      q_rec[l] *= -1.;
    
    //q_rec = k_dmn::parameter_type::back_inside_cluster(q_rec);
    q_rec = cluster_operations::translate_inside_cluster(q_rec, k_dmn::parameter_type::get_super_basis_vectors());

    if(VECTOR_OPERATIONS::L2_NORM(q_rec)<1.e-6){
      q_vector_is_invertible=true;
      //cout << "\n\t q_vec is invertible !!! \n\n";
    }
    else{
      q_vector_is_invertible=false;
      //cout << "\n\t q_vec is NOT invertible !!! \n\n";
    }

    symmetrize_over_matsubara_frequencies(G);
  }

  if(true)
  {
    std::vector<double> q_rec(q_vec);
    
    for(size_t l=0; l<q_rec.size(); ++l)
      q_rec[l] *= 2;
    
    //q_rec = k_dmn::parameter_type::back_inside_cluster(q_rec);
    q_rec = cluster_operations::translate_inside_cluster(q_rec, k_dmn::parameter_type::get_super_basis_vectors());

    if(VECTOR_OPERATIONS::L2_NORM(q_rec)<1.e-6){
      q_vector_is_reciprocal=true;
      //cout << "\n\t q_vec is reciprocal !!! \n\n";
    }
    else{
      q_vector_is_reciprocal=false;
      //cout << "\n\t q_vec is NOT reciprocal !!! \n\n";
    }

    switch(parameters.get_vertex_measurement_type())
      {
      case PARTICLE_HOLE_TRANSVERSE:
	symmetrize_over_pi_rotations_ph(G);
	break;

      case PARTICLE_HOLE_MAGNETIC:
	symmetrize_over_pi_rotations_ph(G);
	break;
      
      case PARTICLE_HOLE_CHARGE:
	symmetrize_over_pi_rotations_ph(G);
	break;
	
      case PARTICLE_PARTICLE_SUPERCONDUCTING:
	set_real(G);
	symmetrize_over_pi_rotations_pp(G);
	break;
	
      default:
	throw std::logic_error(__FUNCTION__);
      }
  }
}

template<class parameters_type>
template<typename scalartype, typename k_dmn, typename w_dmn>
void diagrammatic_symmetries<parameters_type>::execute(function<scalartype, dmn_3<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn>, k_dmn> >& G)
{  
  //cout << __PRETTY_FUNCTION__ << endl;

  symmetrize_over_matsubara_frequencies(G);

  switch(parameters.get_vertex_measurement_type())
    {
    case PARTICLE_HOLE_TRANSVERSE:
      symmetrize_over_pi_rotations_ph(G);
      break;

    case PARTICLE_HOLE_MAGNETIC:
      symmetrize_over_pi_rotations_ph(G);
      break;
      
    case PARTICLE_HOLE_CHARGE:
      symmetrize_over_pi_rotations_ph(G);
      break;
	
    case PARTICLE_PARTICLE_SUPERCONDUCTING:
      //set_real(G);
      symmetrize_over_pi_rotations_pp(G);
      break;
      
    default:
      throw std::logic_error(__FUNCTION__);
    }
}

template<class parameters_type>
template<typename scalartype, typename k_dmn, typename w_dmn>
void diagrammatic_symmetries<parameters_type>::symmetrize_over_matsubara_frequencies(function<scalartype, dmn_2<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn> > >& G)
{
  //cout << __PRETTY_FUNCTION__ << endl;

  if(b::dmn_size()>1)
    throw std::logic_error(__FUNCTION__);

  if(q_vector_is_invertible)
    {
      for(int nu_1=0; nu_1<b::dmn_size(); nu_1++){
	for(int nu_2=0; nu_2<b::dmn_size(); nu_2++){
	  for(int k1=0; k1<k_dmn::dmn_size(); k1++){
	    for(int w1=0; w1<w_dmn::dmn_size(); w1++){
	      
	      for(int mu_1=0; mu_1<b::dmn_size(); mu_1++){
		for(int mu_2=0; mu_2<b::dmn_size(); mu_2++){
		  for(int k2=0; k2<k_dmn::dmn_size(); k2++){
		    for(int w2=0; w2<w_dmn::dmn_size(); w2++){
		      
		      int min_w1 = w_dmn::dmn_size()-1-w1;
		      int min_w2 = w_dmn::dmn_size()-1-w2;
		      
		      std::complex<double> tmp1 = G(nu_1,nu_2,k1,    w1, mu_1,mu_2,k2,     w2);
		      std::complex<double> tmp2 = G(nu_1,nu_2,k1,min_w1, mu_1,mu_2,k2, min_w2);
		      
		      G(nu_1,nu_2,k1,    w1, mu_1,mu_2,k2,     w2) = (     tmp1  + conj(tmp2))/2.;
		      G(nu_1,nu_2,k1,min_w1, mu_1,mu_2,k2, min_w2) = (conj(tmp1) +      tmp2 )/2.;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
}

template<class parameters_type>
template<typename scalartype, typename k_dmn, typename w_dmn>
void diagrammatic_symmetries<parameters_type>::symmetrize_over_matsubara_frequencies(function<scalartype, dmn_3<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn>, k_dmn> >& G)
{
  if(b::dmn_size()>1)
    throw std::logic_error(__FUNCTION__);

  if(parameters.get_w_channel()==0)
    {
      for(int nu_1=0; nu_1<b::dmn_size(); nu_1++){
	for(int nu_2=0; nu_2<b::dmn_size(); nu_2++){
	  for(int k1=0; k1<k_dmn::dmn_size(); k1++){
	    for(int w1=0; w1<w_dmn::dmn_size(); w1++){
	      
	      for(int mu_1=0; mu_1<b::dmn_size(); mu_1++){
		for(int mu_2=0; mu_2<b::dmn_size(); mu_2++){
		  for(int k2=0; k2<k_dmn::dmn_size(); k2++){
		    for(int w2=0; w2<w_dmn::dmn_size(); w2++){
		      
		      for(int q=0; q<k_dmn::dmn_size(); q++){
			int min_w1 = w_dmn::dmn_size()-1-w1;
			int min_w2 = w_dmn::dmn_size()-1-w2;
		      
			std::complex<double> tmp1 = G(nu_1,nu_2,k1,    w1, mu_1,mu_2,k2,     w2, q);
			std::complex<double> tmp2 = G(nu_1,nu_2,k1,min_w1, mu_1,mu_2,k2, min_w2, q);
			
			G(nu_1,nu_2,k1,    w1, mu_1,mu_2,k2,     w2, q) = (     tmp1  + conj(tmp2))/2.;
			G(nu_1,nu_2,k1,min_w1, mu_1,mu_2,k2, min_w2, q) = (conj(tmp1) +      tmp2 )/2.;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
}

template<class parameters_type>
template<typename scalartype, typename k_dmn, typename w_dmn>
void diagrammatic_symmetries<parameters_type>::symmetrize_over_pi_rotations_ph(function<scalartype, dmn_2<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn> > >& G)
{
  //cout << __PRETTY_FUNCTION__ << endl;

  if(b::dmn_size()>1)
    throw std::logic_error(__FUNCTION__);

  for(int nu_1=0; nu_1<b::dmn_size(); nu_1++){
    for(int nu_2=0; nu_2<b::dmn_size(); nu_2++){
      for(int k1=0; k1<k_dmn::dmn_size(); k1++){
	for(int w1=0; w1<w_dmn::dmn_size(); w1++){

	  for(int mu_1=0; mu_1<b::dmn_size(); mu_1++){
	    for(int mu_2=0; mu_2<b::dmn_size(); mu_2++){
	      for(int k2=0; k2<k_dmn::dmn_size(); k2++){
		for(int w2=0; w2<w_dmn::dmn_size(); w2++){

		  int k1_plus_q = k_dmn::parameter_type::add(k1,q_ind);
		  int k2_plus_q = k_dmn::parameter_type::add(k2,q_ind);

// 		  if(q_vector_is_reciprocal)
// 		    {		      
// 		      scalartype tmp1 = G(nu_1,nu_2,k1       ,w1, mu_1,mu_2,k2       ,w2);
// 		      scalartype tmp2 = G(nu_1,nu_2,k2       ,w2, mu_1,mu_2,k1       ,w1);
// 		      scalartype tmp3 = G(nu_1,nu_2,k1_plus_q,w1, mu_1,mu_2,k2_plus_q,w2);
// 		      scalartype tmp4 = G(nu_1,nu_2,k2_plus_q,w2, mu_1,mu_2,k1_plus_q,w1);
		      
// 		      scalartype tmp = (tmp1+conj(tmp2)+conj(tmp3)+tmp4)/4.;
		      
// 		      G(nu_1,nu_2,k1       ,w1, mu_1,mu_2,k2       ,w2) = tmp;
// 		      G(nu_1,nu_2,k2       ,w2, mu_1,mu_2,k1       ,w1) = conj(tmp);
// 		      G(nu_1,nu_2,k1_plus_q,w1, mu_1,mu_2,k2_plus_q,w2) = conj(tmp);
// 		      G(nu_1,nu_2,k2_plus_q,w2, mu_1,mu_2,k1_plus_q,w1) = tmp;
// 		    }
// 		  else

		  //int min_q_ind = k_dmn::parameter_type::subtract(q_ind, k_dmn::parameter_type::get_k_0_index());
		  int min_q_ind = k_dmn::parameter_type::subtract(q_ind, k_dmn::parameter_type::origin_index());

		  if(min_q_ind == q_ind)
		    {
		      scalartype tmp1 =      G(nu_1,nu_2,k1       ,w1, mu_1,mu_2,k2       ,w2);
		      scalartype tmp2 = conj(G(nu_1,nu_2,k2_plus_q,w2, mu_1,mu_2,k1_plus_q,w1));

		      scalartype tmp = (tmp1+tmp2)/2.;

		       G(nu_1,nu_2,k1       ,w1, mu_1,mu_2,k2       ,w2) =      tmp;
		       G(nu_1,nu_2,k2_plus_q,w2, mu_1,mu_2,k1_plus_q,w1) = conj(tmp);
		    }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

template<class parameters_type>
template<typename scalartype, typename k_dmn, typename w_dmn>
void diagrammatic_symmetries<parameters_type>::symmetrize_over_pi_rotations_ph(function<scalartype, dmn_3<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn>, k_dmn> >& G)
{
  if(b::dmn_size()>1)
    throw std::logic_error(__FUNCTION__);

  for(int nu_1=0; nu_1<b::dmn_size(); nu_1++){
    for(int nu_2=0; nu_2<b::dmn_size(); nu_2++){
      for(int k1=0; k1<k_dmn::dmn_size(); k1++){
	for(int w1=0; w1<w_dmn::dmn_size(); w1++){

	  for(int mu_1=0; mu_1<b::dmn_size(); mu_1++){
	    for(int mu_2=0; mu_2<b::dmn_size(); mu_2++){
	      for(int k2=0; k2<k_dmn::dmn_size(); k2++){
		for(int w2=0; w2<w_dmn::dmn_size(); w2++){

		  for(int q=0; q<k_dmn::dmn_size(); q++){

		    int k1_plus_q = k_dmn::parameter_type::add(k1,q);
		    int k2_plus_q = k_dmn::parameter_type::add(k2,q);

		    //int min_q = k_dmn::parameter_type::subtract(q, k_dmn::parameter_type::get_k_0_index());
		    int min_q = k_dmn::parameter_type::subtract(q, k_dmn::parameter_type::origin_index());
		    
		    scalartype tmp1 = G(nu_1,nu_2,k1       ,w1, mu_1,mu_2,k2       ,w2,     q);
		    scalartype tmp2 = G(nu_1,nu_2,k2       ,w2, mu_1,mu_2,k1       ,w1,     q);
		    scalartype tmp3 = G(nu_1,nu_2,k1_plus_q,w1, mu_1,mu_2,k2_plus_q,w2, min_q);
		    scalartype tmp4 = G(nu_1,nu_2,k2_plus_q,w2, mu_1,mu_2,k1_plus_q,w1, min_q);
			
		    scalartype tmp = (tmp1+conj(tmp2)+conj(tmp3)+tmp4)/4.;
			
		    G(nu_1,nu_2,k1       ,w1, mu_1,mu_2,k2       ,w2) = tmp;
		    G(nu_1,nu_2,k2       ,w2, mu_1,mu_2,k1       ,w1) = conj(tmp);
		    G(nu_1,nu_2,k1_plus_q,w1, mu_1,mu_2,k2_plus_q,w2) = conj(tmp);
		    G(nu_1,nu_2,k2_plus_q,w2, mu_1,mu_2,k1_plus_q,w1) = tmp;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

template<class parameters_type>
template<typename scalartype, typename k_dmn, typename w_dmn>
void diagrammatic_symmetries<parameters_type>::symmetrize_over_pi_rotations_pp(function<scalartype, dmn_2<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn> > >& G)
{
  //cout << __PRETTY_FUNCTION__ << endl;

  if(b::dmn_size()>1)
    throw std::logic_error(__FUNCTION__);

  for(int nu_1=0; nu_1<b::dmn_size(); nu_1++){
    for(int nu_2=0; nu_2<b::dmn_size(); nu_2++){
      for(int k1=0; k1<k_dmn::dmn_size(); k1++){
	for(int w1=0; w1<w_dmn::dmn_size(); w1++){

	  for(int mu_1=0; mu_1<b::dmn_size(); mu_1++){
	    for(int mu_2=0; mu_2<b::dmn_size(); mu_2++){
	      for(int k2=0; k2<k_dmn::dmn_size(); k2++){
		for(int w2=0; w2<w_dmn::dmn_size(); w2++){

		  int q_min_k1 = k_dmn::parameter_type::subtract(k1,q_ind);
		  int q_min_k2 = k_dmn::parameter_type::subtract(k2,q_ind);

		  int min_w1 = w_dmn::dmn_size()-1-w1;
		  int min_w2 = w_dmn::dmn_size()-1-w2;

		  scalartype tmp1 = G(nu_1,nu_2,      k1,    w1, mu_1,mu_2,      k2,     w2);
		  scalartype tmp2 = G(nu_1,nu_2,q_min_k1,min_w1, mu_1,mu_2,q_min_k2, min_w2);
		  scalartype tmp3 = G(nu_1,nu_2,      k2,    w2, mu_1,mu_2,      k1,     w1);
		  scalartype tmp4 = G(nu_1,nu_2,q_min_k2,min_w2, mu_1,mu_2,q_min_k1, min_w1);
		      
		  scalartype tmp = (tmp1+tmp2+conj(tmp3)+conj(tmp4))/4.;
		  
		  G(nu_1,nu_2,      k1,    w1, mu_1,mu_2,      k2,     w2) = tmp;
		  G(nu_1,nu_2,q_min_k1,min_w1, mu_1,mu_2,q_min_k2, min_w2) = tmp;
		  G(nu_1,nu_2,      k2,    w2, mu_1,mu_2,      k1,     w1) = conj(tmp);
		  G(nu_1,nu_2,q_min_k2,min_w2, mu_1,mu_2,q_min_k1, min_w1) = conj(tmp);
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

template<class parameters_type>
template<typename scalartype, typename k_dmn, typename w_dmn>
void diagrammatic_symmetries<parameters_type>::symmetrize_over_pi_rotations_pp(function<scalartype, dmn_3<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn>, k_dmn> >& G)
{
  if(b::dmn_size()>1)
    throw std::logic_error(__FUNCTION__);

  for(int nu_1=0; nu_1<b::dmn_size(); nu_1++){
    for(int nu_2=0; nu_2<b::dmn_size(); nu_2++){
      for(int k1=0; k1<k_dmn::dmn_size(); k1++){
	for(int w1=0; w1<w_dmn::dmn_size(); w1++){

	  for(int mu_1=0; mu_1<b::dmn_size(); mu_1++){
	    for(int mu_2=0; mu_2<b::dmn_size(); mu_2++){
	      for(int k2=0; k2<k_dmn::dmn_size(); k2++){
		for(int w2=0; w2<w_dmn::dmn_size(); w2++){

		  for(int q=0; q<k_dmn::dmn_size(); q++){

		    int q_min_k1 = k_dmn::parameter_type::subtract(k1,q);
		    int q_min_k2 = k_dmn::parameter_type::subtract(k2,q);
		    
		    int min_w1 = w_dmn::dmn_size()-1-w1;
		    int min_w2 = w_dmn::dmn_size()-1-w2;
		    
		    scalartype tmp1 = G(nu_1,nu_2,      k1,    w1, mu_1,mu_2,      k2,     w2, q);
		    scalartype tmp2 = G(nu_1,nu_2,q_min_k1,min_w1, mu_1,mu_2,q_min_k2, min_w2, q);
		    scalartype tmp3 = G(nu_1,nu_2,      k2,    w2, mu_1,mu_2,      k1,     w1, q);
		    scalartype tmp4 = G(nu_1,nu_2,q_min_k2,min_w2, mu_1,mu_2,q_min_k1, min_w1, q);
		    
		    scalartype tmp = (tmp1+tmp2+conj(tmp3)+conj(tmp4))/4.;
		    
		    G(nu_1,nu_2,      k1,    w1, mu_1,mu_2,      k2,     w2, q) = tmp;
		    G(nu_1,nu_2,q_min_k1,min_w1, mu_1,mu_2,q_min_k2, min_w2, q) = tmp;
		    G(nu_1,nu_2,      k2,    w2, mu_1,mu_2,      k1,     w1, q) = conj(tmp);
		    G(nu_1,nu_2,q_min_k2,min_w2, mu_1,mu_2,q_min_k1, min_w1, q) = conj(tmp);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

template<class parameters_type>
template<typename scalartype, typename k_dmn, typename w_dmn>
void diagrammatic_symmetries<parameters_type>::set_real(function<scalartype, dmn_2<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn> > >& G)
{
  if(b::dmn_size()>1)
    throw std::logic_error(__FUNCTION__);

  for(int nu_1=0; nu_1<b::dmn_size(); nu_1++){
    for(int nu_2=0; nu_2<b::dmn_size(); nu_2++){
      for(int k1=0; k1<k_dmn::dmn_size(); k1++){
	for(int w1=0; w1<w_dmn::dmn_size(); w1++){

	  for(int mu_1=0; mu_1<b::dmn_size(); mu_1++){
	    for(int mu_2=0; mu_2<b::dmn_size(); mu_2++){
	      for(int k2=0; k2<k_dmn::dmn_size(); k2++){
		for(int w2=0; w2<w_dmn::dmn_size(); w2++){

		  scalartype tmp1 = G(nu_1,nu_2,k1,w1, mu_1,mu_2,k2,w2);

		  G(nu_1,nu_2,k1,w1, mu_1,mu_2,k2,w2) = real(tmp1);
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

template<class parameters_type>
template<typename scalartype, typename k_dmn, typename w_dmn>
void diagrammatic_symmetries<parameters_type>::set_real(function<scalartype, dmn_3<dmn_4<b,b,k_dmn,w_dmn>, dmn_4<b,b,k_dmn,w_dmn>, k_dmn> >& G)
{
  if(b::dmn_size()>1)
    throw std::logic_error(__FUNCTION__);

  for(int nu_1=0; nu_1<b::dmn_size(); nu_1++){
    for(int nu_2=0; nu_2<b::dmn_size(); nu_2++){
      for(int k1=0; k1<k_dmn::dmn_size(); k1++){
	for(int w1=0; w1<w_dmn::dmn_size(); w1++){

	  for(int mu_1=0; mu_1<b::dmn_size(); mu_1++){
	    for(int mu_2=0; mu_2<b::dmn_size(); mu_2++){
	      for(int k2=0; k2<k_dmn::dmn_size(); k2++){
		for(int w2=0; w2<w_dmn::dmn_size(); w2++){

		  for(int q=0; q<k_dmn::dmn_size(); q++){

		    scalartype tmp1 = G(nu_1,nu_2,k1,w1, mu_1,mu_2,k2,w2, q);
		    
		    G(nu_1,nu_2,k1,w1, mu_1,mu_2,k2,w2, q) = real(tmp1);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}


#endif
