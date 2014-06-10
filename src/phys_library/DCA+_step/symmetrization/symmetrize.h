//-*-C++-*-

#ifndef SYMMETRIZE_H
#define SYMMETRIZE_H

/*! \defgroup SYMMETRIZE-FUNCTIONS
*   \ingroup  ALGORITHMS
*/

#include "symmetrize_single_particle_function.h"
#include "symmetrize_two_particle_function.h"

/*! \class   symmetrize
 *  \ingroup SYMMETRIZE-FUNCTIONS
 *
 *  \author Peter Staar
 *  \brief  This class symmetrizes Greens-functions according to cluster-symmetries, matsubara frequencies and band-index symmetries.
 */
class symmetrize : public symmetrize_single_particle_function,
		   public symmetrize_two_particle_function
{
#include "type_definitions.h"

public:

  template<typename scalartype, typename f_dmn_0>
  static void execute(function<scalartype, f_dmn_0>& f, bool do_diff=false); 

  template<typename scalartype, typename nu_dmn_t, typename f_dmn_0>
  static void execute(function<scalartype, dmn_3<nu_dmn_t, nu_dmn_t, f_dmn_0> >& f, bool do_diff=false); 

  template<typename scalartype, typename nu_dmn_t, typename f_dmn_0, typename f_dmn_1>
  static void execute(function<scalartype, dmn_4<nu_dmn_t, nu_dmn_t, f_dmn_0, f_dmn_1> >& f, bool do_diff=false); 

  template<typename scalartype, typename f_dmn_0, typename f_dmn_1>
  static void execute(function<scalartype, dmn_4<nu, nu, f_dmn_0, f_dmn_1> >& f, 
		      function<int       , nu_nu>&                            H_symmetry,
		      bool                                                    do_diff=false);



  
  template<typename scalartype, typename k_dmn_t, typename w_dmn_t>
  static void execute(function<scalartype, dmn_2< dmn_4<b, b, k_dmn_t, w_dmn_t>,  dmn_4<b, b, k_dmn_t, w_dmn_t> > >& f, 
		      std::vector<double>                                                                            Q,
		      bool do_diff=false);

  template<typename scalartype, typename k_dmn_t, typename w_dmn_t>
  static void execute(function<scalartype, dmn_2< dmn_4<b, b, k_dmn_t, w_dmn_t>,  dmn_4<b, b, k_dmn_t, w_dmn_t> > >& f, 
		      function<int       , nu_nu>&  H_symmetry,
		      std::vector<double>                                                                            Q,
		      bool do_diff=false);

  template<typename scalartype, typename k_dmn_t, typename w_dmn_t>
  static void execute(function<scalartype, dmn_8<b, b, b, b, k_dmn_t, k_dmn_t, w_dmn_t, w_dmn_t> >& f, 
		      std::vector<double>                                                           Q,
		      bool do_diff=false);

  template<typename scalartype, typename k_dmn_t, typename w_dmn_t>
  static void execute(function<scalartype, dmn_8<b, b, b, b, k_dmn_t, k_dmn_t, w_dmn_t, w_dmn_t> >& f, 
		      function<int       , nu_nu>&  H_symmetry,
		      std::vector<double>                                                           Q,
		      bool do_diff=false);
  
  template<typename scalartype, typename k_dmn_t, typename w_dmn_t>
  static void execute(function<scalartype, dmn_3< dmn_4<b, b, k_dmn_t, w_dmn_t>,  dmn_4<b, b, k_dmn_t, w_dmn_t>, k_dmn_t> >& f, bool do_diff=false);

};

template<typename scalartype, typename f_dmn_0>
void symmetrize::execute(function<scalartype, f_dmn_0>& f, bool do_diff)
{
  symmetrize_single_particle_function::execute(f, do_diff);
}

template<typename scalartype, typename f_dmn_0, typename f_dmn_1>
void symmetrize::execute(function<scalartype, dmn_4<nu, nu, f_dmn_0, f_dmn_1> >& f, 
			 function<int       , nu_nu>&                            H_symmetry,
			 bool                                                    do_diff)
{
  symmetrize_single_particle_function::execute(f, H_symmetry, do_diff);
}

template<typename scalartype, typename nu_dmn_t, typename f_dmn_0>
void symmetrize::execute(function<scalartype, dmn_3<nu_dmn_t, nu_dmn_t, f_dmn_0> >& f, bool do_diff)
{
  symmetrize_single_particle_function::execute(f, do_diff);
}

template<typename scalartype, typename nu_dmn_t, typename f_dmn_0, typename f_dmn_1>
void symmetrize::execute(function<scalartype, dmn_4<nu_dmn_t, nu_dmn_t, f_dmn_0, f_dmn_1> >& f, bool do_diff)
{
  symmetrize_single_particle_function::execute(f, do_diff);
}

template<typename scalartype, typename k_dmn_t, typename w_dmn_t>
void symmetrize::execute(function<scalartype, dmn_8<b, b, b, b, k_dmn_t, k_dmn_t, w_dmn_t, w_dmn_t> >& f, 
			 function<int       , nu_nu>&                                                  H_symmetry,
			 std::vector<double>                                                           Q,
			 bool                                                                          do_diff)
{
  symmetrize_two_particle_function::execute(f, Q, do_diff);
}

template<typename scalartype, typename k_dmn_t, typename w_dmn_t>
void symmetrize::execute(function<scalartype, dmn_2< dmn_4<b, b, k_dmn_t, w_dmn_t>,  dmn_4<b, b, k_dmn_t, w_dmn_t> > >& f, 
			 function<int       , nu_nu>&                                                                   H_symmetry,
			 std::vector<double>                                                                            Q,
			 bool do_diff)
{
  symmetrize_two_particle_function::execute(f, Q, do_diff);
}

template<typename scalartype, typename k_dmn_t, typename w_dmn_t>
void symmetrize::execute(function<scalartype, dmn_8<b, b, b, b, k_dmn_t, k_dmn_t, w_dmn_t, w_dmn_t> >& f, 
			 std::vector<double>                                                           Q,
			 bool do_diff)
{
  symmetrize_two_particle_function::execute(f, Q, do_diff);
}

template<typename scalartype, typename k_dmn_t, typename w_dmn_t>
void symmetrize::execute(function<scalartype, dmn_2< dmn_4<b, b, k_dmn_t, w_dmn_t>,  dmn_4<b, b, k_dmn_t, w_dmn_t> > >& f, 
			 std::vector<double>                                                                            Q,			                  
			 bool do_diff)
{
  symmetrize_two_particle_function::execute(f, Q, do_diff);
}

template<typename scalartype, typename k_dmn_t, typename w_dmn_t>
void symmetrize::execute(function<scalartype, dmn_3< dmn_4<b, b, k_dmn_t, w_dmn_t>,  dmn_4<b, b, k_dmn_t, w_dmn_t>, k_dmn_t> >& f, bool do_diff)
{
  symmetrize_two_particle_function::execute(f, do_diff);
}


#endif











































































// template<typename scalartype, typename f_dmn_0, typename f_dmn_1>
// void symmetrize::execute(function<scalartype, dmn_4<nu, nu, f_dmn_0, f_dmn_1> >& f, function<int, nu_nu>& H_symmetry,
// 			 bool                                                                             do_diff)
// {
//   //   cout << __PRETTY_FUNCTION__ << endl;

//   symmetrize_over_Hamiltonians(f, H_symmetry, do_diff);

//   if(IS_EQUAL_TYPE<f_dmn_0, t>::check || IS_EQUAL_TYPE<f_dmn_1, t>::check)
//     symmetrize_over_time_domain(f, do_diff);

//   if(IS_EQUAL_TYPE<f_dmn_0, w>::check || IS_EQUAL_TYPE<f_dmn_1, w>::check)
//     symmetrize_over_frequency_domain(f, do_diff);

//   if(IS_EQUAL_TYPE<f_dmn_0, r_DCA>::check || IS_EQUAL_TYPE<f_dmn_1, r_DCA>::check)
//     symmetrize_over_r_DCA(f, do_diff);

//   if(IS_EQUAL_TYPE<f_dmn_0, k_DCA>::check || IS_EQUAL_TYPE<f_dmn_1, k_DCA>::check)
//     symmetrize_over_k_DCA(f, do_diff);
// }


/*
template<typename scalartype>
scalartype symmetrize::conjugate(scalartype x)
{
  return std::conj(x);
}

template<>
int symmetrize::conjugate(int x)
{
  return x;
}


template<>
float symmetrize::conjugate(float x)
{
  return x;
}

template<>
double symmetrize::conjugate(double x)
{
  return x;
}

template<typename scalartype>
bool symmetrize::difference(scalartype x1, scalartype x2, const char* string)
{
  //cout << __FUNCTION__ << "\n";

  bool OK = true;

  if((abs(x1)<1.e-1 && abs(x1-x2)>1.e-6) || (abs(x1)>1.e-1 && abs(x1-x2)/abs(x1)>1.e-6)){
    cout << __FUNCTION__ << "\t" << x1 << "\t" << x2 << "\t" <<abs(x1-x2) << "\t" << string << endl;
    OK = true;
    throw std::logic_error(__FUNCTION__);
  }

  return OK;
}


void symmetrize::permute_band_indices(int l,
				      int    l1, int    l2, int    l3, int    l4,
				      int& p_l1, int& p_l2, int& p_l3, int& p_l4)
{
//   p_l1 = l3;
//   p_l2 = l4;
//   p_l3 = l1;
//   p_l4 = l2;

  std::pair<int,int> initial_state = model::get_orbital_permutations()[l].first;
  std::pair<int,int> final_state   = model::get_orbital_permutations()[l].second;

  p_l1=l1;
  if(l1==initial_state.first)  p_l1 = final_state.first;
  if(l1==initial_state.second) p_l1 = final_state.second;

  p_l2=l2;
  if(l2==initial_state.first)  p_l2 = final_state.first;
  if(l2==initial_state.second) p_l2 = final_state.second;

  p_l3=l3;
  if(l3==initial_state.first)  p_l3 = final_state.first;
  if(l3==initial_state.second) p_l3 = final_state.second;

  p_l4=l4;
  if(l4==initial_state.first)  p_l4 = final_state.first;
  if(l4==initial_state.second) p_l4 = final_state.second;
}

template<typename scalartype>
void symmetrize::symmetrize_over_Hamiltonians(function<scalartype, symmetrize::b_b__b_b__k_DCA_VERTEX_k_DCA_w_VERTEX_w_VERTEX>& f, 
					      function<int       , nu_nu>&                            H_symmetry,
					      bool do_diff)
{
  cout << scientific;
  cout.precision(6);

  scalartype result=0;

  int p_l1, p_l2, p_l3, p_l4;

  for(size_t l=0; l<model::get_orbital_permutations().size(); l++)
    {
      for(int w1=0; w1<w_VERTEX::dmn_size(); w1++){ 
	for(int w2=0; w2<w_VERTEX::dmn_size(); w2++){ 
	  
	  for(int k1=0; k1<k_DCA::dmn_size(); k1++){ 
	    for(int k2=0; k2<k_DCA::dmn_size(); k2++){ 
	      
	      for(int l1=0; l1<b::dmn_size(); l1++){

		for(int l2=0; l2<b::dmn_size(); l2++){

		  for(int l3=0; l3<b::dmn_size(); l3++){

		    for(int l4=0; l4<b::dmn_size(); l4++){

		      permute_band_indices( 0,
					    l1,    l2,    l3,    l4,
					    p_l1,  p_l2,  p_l3,  p_l4);

		      result = (f(l1,l2,l3,l4, k1,k2, w1,w2) + f(p_l1,p_l2, p_l3,p_l4, k1,k2, w1,w2))/2.;
	
		      f(  l1,  l2,   l3,  l4, k1,k2, w1,w2) = result;
		      f(p_l1,p_l2, p_l3,p_l4, k1,k2, w1,w2) = result;
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

template<typename scalartype>
void symmetrize::symmetrize_over_Hamiltonians(function<scalartype, symmetrize::b_b_k_DCA_w_VERTEX_EXTENDED__b_b_k_DCA_w_VERTEX_EXTENDED>& f, 
					      function<int       , nu_nu>&                            H_symmetry,
					      bool do_diff)
{
  cout << scientific;
  cout.precision(6);

  scalartype result=0;

  int p_l1, p_l2, p_l3, p_l4;

  for(size_t l=0; l<model::get_orbital_permutations().size(); l++)
    {
      for(int w1=0; w1<w_VERTEX::dmn_size(); w1++){ 
	for(int w2=0; w2<w_VERTEX::dmn_size(); w2++){ 
	  
	  for(int k1=0; k1<k_DCA::dmn_size(); k1++){ 
	    for(int k2=0; k2<k_DCA::dmn_size(); k2++){ 
	      
	      for(int l1=0; l1<b::dmn_size(); l1++){

		for(int l2=0; l2<b::dmn_size(); l2++){

		  for(int l3=0; l3<b::dmn_size(); l3++){

		    for(int l4=0; l4<b::dmn_size(); l4++){

		      permute_band_indices( 0,
					    l1,    l2,    l3,    l4,
					    p_l1,  p_l2,  p_l3,  p_l4);

		      result = (f(l1,l2,k1,w1, l3,l4,k2,w2) + f(p_l1,p_l2,k1,w1, p_l3,p_l4,k2,w2))/2.;

		      f(  l1,  l2,k1,w1,   l3,  l4,k2,w2) = result;
		      f(p_l1,p_l2,k1,w1, p_l3,p_l4,k2,w2) = result;
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

template<typename scalartype>
void symmetrize::execute(function<scalartype, symmetrize::b_b__b_b__k_DCA_VERTEX_k_DCA_w_VERTEX_w_VERTEX>& f, 
			 function<int       , nu_nu                                                     >& H_symmetry,
			 bool do_diff)
{
  //cout << __PRETTY_FUNCTION__ << endl;

  if(true)
    {// cluster-symmetries G4_k_k_w_w for first k
      int* coor = new int[f.signature()];
      function<std::complex<double>, b_b__b_b__k_DCA_VERTEX_k_DCA_w_VERTEX_w_VERTEX> G4_k_k_w_w;
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	  double weight = DCA_k_cluster_type::get_weight(coor[4]);
	  DCA_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[4], coor[5]);
	    
	  G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]) += f(i)/weight;
	}
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	    
	  DCA_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[4], coor[5]);

	  if(do_diff){
	    assert(difference(f(i), G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]), __FUNCTION__));
	  }

	  f(i) = G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]);
	}

      delete [] coor;
    }

  if(true)
    {// for second k
      int* coor = new int[f.signature()];
      function<std::complex<double>, b_b__b_b__k_DCA_VERTEX_k_DCA_w_VERTEX_w_VERTEX> G4_k_k_w_w;
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	  double weight = DCA_k_cluster_type::get_weight(coor[5]);
	  DCA_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[5], coor[4]);

	  G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]) += f(i)/weight;
	}
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	  DCA_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[5], coor[4]);
	    
	  if(do_diff){
	    assert(difference(f(i), G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]), __FUNCTION__));
	  }

	  f(i) = G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]);
	}
	
      delete [] coor;
    }

  if(true)
    {// G4 is real in r,r,t,t
      int N_w = w_VERTEX::dmn_size()-1;
	
      int* coor = new int[f.signature()];
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);

	  int min_k1 = DCA_k_cluster_type::subtract(coor[4], DCA_k_cluster_type::get_k_0_index());
	  int min_k2 = DCA_k_cluster_type::subtract(coor[5], DCA_k_cluster_type::get_k_0_index());

	  std::complex<double> c = (f(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7])
				    + std::conj(f(coor[0], coor[1], coor[2], coor[3], min_k1, min_k2, N_w - coor[6], N_w - coor[7])))/2.;
	    
// 	  if(true){
// 	    difference(          f(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5],       coor[6],       coor[7]),
// 		       std::conj(f(coor[0], coor[1], coor[2], coor[3], min_k1 , min_k2 , N_w - coor[6], N_w - coor[7])), 
// 			      __FUNCTION__);
// 	  }

	  f(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5],       coor[6],       coor[7]) = c;           
	  f(coor[0], coor[1], coor[2], coor[3], min_k1 , min_k2 , N_w - coor[6], N_w - coor[7]) = std::conj(c);
	}
	
      delete [] coor;
    }

  if(true)
    symmetrize_over_Hamiltonians(f, H_symmetry, do_diff);
}

void symmetrize::execute(function<std::complex<double>, b_b_k_DCA_w_VERTEX_EXTENDED__b_b_k_DCA_w_VERTEX_EXTENDED>& f, 
			 function<int       , nu_nu                                                             >& H_symmetry,
			 bool do_diff)
{
  //cout << __PRETTY_FUNCTION__ << endl;

  if(true)
    {// cluster-symmetries G4_k_k_w_w for first k
      int* coor = new int[f.signature()];
      function<std::complex<double>, b_b_k_DCA_w_VERTEX_EXTENDED__b_b_k_DCA_w_VERTEX_EXTENDED> G4_k_k_w_w;
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	  double weight = DCA_k_cluster_type::get_weight(coor[2]);
	  DCA_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[2], coor[6]);
	    
	  G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]) += f(i)/weight;
	}
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	    
	  DCA_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[2], coor[6]);

// 	  if(true){
// 	    assert(difference(f(i), G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]), __FUNCTION__));
// 	  }

	  f(i) = G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]);
	}

      delete [] coor;
    }

  if(true)
    {// cluster-symmetries G4_k_k_w_w for first k
      int* coor = new int[f.signature()];
      function<std::complex<double>, b_b_k_DCA_w_VERTEX_EXTENDED__b_b_k_DCA_w_VERTEX_EXTENDED> G4_k_k_w_w;
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	  double weight = DCA_k_cluster_type::get_weight(coor[6]);
	  DCA_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[6], coor[2]);
	    
	  G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]) += f(i)/weight;
	}
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	    
	  DCA_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[6], coor[2]);

// 	  if(do_diff){
// 	    assert(difference(f(i), G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]), __FUNCTION__));
// 	  }

	  f(i) = G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]);
	}

      delete [] coor;
    }

  if(true)
    {// G4 is real in r,r,t,t
//       cout << f.get_name() << endl; 

      int N_w = w_VERTEX_EXTENDED::dmn_size()-1;
	
      int* coor = new int[f.signature()];
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);

	  int min_k1 = DCA_k_cluster_type::subtract(coor[2], DCA_k_cluster_type::get_k_0_index());
	  int min_k2 = DCA_k_cluster_type::subtract(coor[6], DCA_k_cluster_type::get_k_0_index());

	  std::complex<double> c = (            f(coor[0], coor[1], coor[2],       coor[3], coor[4], coor[5], coor[6],      coor[7])
				    + std::conj(f(coor[0], coor[1], min_k1,  N_w - coor[3], coor[4], coor[5], min_k2, N_w - coor[7])))/2.;
	    
// 	  if(true){
// 	    difference(          f(coor[0], coor[1], coor[2],       coor[3], coor[4], coor[5], coor[6],       coor[7]),
// 		       std::conj(f(coor[0], coor[1], min_k1 , N_w - coor[3], coor[4], coor[5], min_k2 , N_w - coor[7])), 
// 			      __FUNCTION__);
// 	  }

	  f(coor[0], coor[1], coor[2],       coor[3], coor[4], coor[5], coor[6],       coor[7]) = c;           
	  f(coor[0], coor[1], min_k1 , N_w - coor[3], coor[4], coor[5], min_k2 , N_w - coor[7]) = std::conj(c);
	}
	
      delete [] coor;
    }

  if(true)
    symmetrize_over_Hamiltonians(f, H_symmetry, do_diff);
}








void symmetrize::execute(function<std::complex<double>, b_b_k_PCM_w_VERTEX_EXTENDED__b_b_k_PCM_w_VERTEX_EXTENDED>& f, 
			 function<int       , nu_nu                                                             >& H_symmetry,
			 bool do_diff)
{
  //cout << __PRETTY_FUNCTION__ << endl;

  if(true)
    {// cluster-symmetries G4_k_k_w_w for first k
      int* coor = new int[f.signature()];
      function<std::complex<double>, b_b_k_PCM_w_VERTEX_EXTENDED__b_b_k_PCM_w_VERTEX_EXTENDED> G4_k_k_w_w;
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	  double weight = PCM_k_cluster_type::get_weight(coor[2]);
	  PCM_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[2], coor[6]);
	    
	  G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]) += f(i)/weight;
	}
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	    
	  PCM_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[2], coor[6]);

// 	  if(true){
// 	    assert(difference(f(i), G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]), __FUNCTION__));
// 	  }

	  f(i) = G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]);
	}

      delete [] coor;
    }

  if(true)
    {// cluster-symmetries G4_k_k_w_w for first k
      int* coor = new int[f.signature()];
      function<std::complex<double>, b_b_k_PCM_w_VERTEX_EXTENDED__b_b_k_PCM_w_VERTEX_EXTENDED> G4_k_k_w_w;
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	  double weight = PCM_k_cluster_type::get_weight(coor[6]);
	  PCM_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[6], coor[2]);
	    
	  G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]) += f(i)/weight;
	}
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	    
	  PCM_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[6], coor[2]);

// 	  if(do_diff){
// 	    assert(difference(f(i), G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]), __FUNCTION__));
// 	  }

	  f(i) = G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]);
	}

      delete [] coor;
    }

  if(true)
    {// G4 is real in r,r,t,t
//       cout << f.get_name() << endl; 

      int N_w = w_VERTEX_EXTENDED::dmn_size()-1;
	
      int* coor = new int[f.signature()];
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);

	  int min_k1 = PCM_k_cluster_type::subtract(coor[2], PCM_k_cluster_type::get_k_0_index());
	  int min_k2 = PCM_k_cluster_type::subtract(coor[6], PCM_k_cluster_type::get_k_0_index());

	  std::complex<double> c = (            f(coor[0], coor[1], coor[2],       coor[3], coor[4], coor[5], coor[6],      coor[7])
				    + std::conj(f(coor[0], coor[1], min_k1,  N_w - coor[3], coor[4], coor[5], min_k2, N_w - coor[7])))/2.;
	    
// 	  if(true){
// 	    difference(          f(coor[0], coor[1], coor[2],       coor[3], coor[4], coor[5], coor[6],       coor[7]),
// 		       std::conj(f(coor[0], coor[1], min_k1 , N_w - coor[3], coor[4], coor[5], min_k2 , N_w - coor[7])), 
// 			      __FUNCTION__);
// 	  }

	  f(coor[0], coor[1], coor[2],       coor[3], coor[4], coor[5], coor[6],       coor[7]) = c;           
	  f(coor[0], coor[1], min_k1 , N_w - coor[3], coor[4], coor[5], min_k2 , N_w - coor[7]) = std::conj(c);
	}
	
      delete [] coor;
    }

//   if(true)
//     symmetrize_over_Hamiltonians(f, H_symmetry, do_diff);
}


















void symmetrize::execute(function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX>& f, 
			 function<int       , nu_nu                                           >& H_symmetry,
			 bool do_diff)
{
  //cout << __PRETTY_FUNCTION__ << endl;

  if(true)
    {// cluster-symmetries G4_k_k_w_w for first k
      int* coor = new int[f.signature()];
      function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> G4_k_k_w_w;
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	  double weight = DCA_k_cluster_type::get_weight(coor[2]);
	  DCA_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[2], coor[6]);
	    
	  G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]) += f(i)/weight;
	}
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	    
	  DCA_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[2], coor[6]);

// 	  if(true){
// 	    assert(difference(f(i), G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]), __FUNCTION__));
// 	  }

	  f(i) = G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]);
	}

      delete [] coor;
    }

  if(true)
    {// cluster-symmetries G4_k_k_w_w for first k
      int* coor = new int[f.signature()];
      function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> G4_k_k_w_w;
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	  double weight = DCA_k_cluster_type::get_weight(coor[6]);
	  DCA_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[6], coor[2]);
	    
	  G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]) += f(i)/weight;
	}
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	    
	  DCA_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[6], coor[2]);

// 	  if(do_diff){
// 	    assert(difference(f(i), G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]), __FUNCTION__));
// 	  }

	  f(i) = G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]);
	}

      delete [] coor;
    }

  if(true)
    {// G4 is real in r,r,t,t

      int N_w = w_VERTEX::dmn_size()-1;
	
      int* coor = new int[f.signature()];
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);

	  int min_k1 = DCA_k_cluster_type::subtract(coor[2], DCA_k_cluster_type::get_k_0_index());
	  int min_k2 = DCA_k_cluster_type::subtract(coor[6], DCA_k_cluster_type::get_k_0_index());

	  std::complex<double> c = (            f(coor[0], coor[1], coor[2],       coor[3], coor[4], coor[5], coor[6],      coor[7])
				    + std::conj(f(coor[0], coor[1], min_k1,  N_w - coor[3], coor[4], coor[5], min_k2, N_w - coor[7])))/2.;
	    
// 	  if(true){
// 	    difference(          f(coor[0], coor[1], coor[2],       coor[3], coor[4], coor[5], coor[6],       coor[7]),
// 		       std::conj(f(coor[0], coor[1], min_k1 , N_w - coor[3], coor[4], coor[5], min_k2 , N_w - coor[7])), 
// 			      __FUNCTION__);
// 	  }

	  f(coor[0], coor[1], coor[2],       coor[3], coor[4], coor[5], coor[6],       coor[7]) = c;           
	  f(coor[0], coor[1], min_k1 , N_w - coor[3], coor[4], coor[5], min_k2 , N_w - coor[7]) = std::conj(c);
	}
	
      delete [] coor;
    }

//   if(true)
//     symmetrize_over_Hamiltonians(f, H_symmetry, do_diff);
}

void symmetrize::execute(function<std::complex<double>, b_b_k_PCM_w_VERTEX__b_b_k_PCM_w_VERTEX>& f, 
			 function<int       , nu_nu                                           >& H_symmetry,
			 bool do_diff)
{
  //cout << __PRETTY_FUNCTION__ << endl;

  if(true)
    {// cluster-symmetries G4_k_k_w_w for first k
      int* coor = new int[f.signature()];
      function<std::complex<double>, b_b_k_PCM_w_VERTEX__b_b_k_PCM_w_VERTEX> G4_k_k_w_w;
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	  double weight = PCM_k_cluster_type::get_weight(coor[2]);
	  PCM_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[2], coor[6]);
	    
	  G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]) += f(i)/weight;
	}
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	    
	  PCM_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[2], coor[6]);

// 	  if(true){
// 	    assert(difference(f(i), G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]), __FUNCTION__));
// 	  }

	  f(i) = G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]);
	}

      delete [] coor;
    }

  if(true)
    {// cluster-symmetries G4_k_k_w_w for second k
      int* coor = new int[f.signature()];
      function<std::complex<double>, b_b_k_PCM_w_VERTEX__b_b_k_PCM_w_VERTEX> G4_k_k_w_w;
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	  double weight = PCM_k_cluster_type::get_weight(coor[6]);
	  PCM_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[6], coor[2]);
	    
	  G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]) += f(i)/weight;
	}
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);
	    
	  PCM_k_cluster_type::convert_k1_k2_to_k1_irr_k2(coor[6], coor[2]);

// 	  if(do_diff){
// 	    assert(difference(f(i), G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]), __FUNCTION__));
// 	  }

	  f(i) = G4_k_k_w_w(coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], coor[6], coor[7]);
	}

      delete [] coor;
    }

  if(true)
    {// G4 is real in r,r,t,t

      int N_w = w_VERTEX::dmn_size()-1;
	
      int* coor = new int[f.signature()];
	
      for(int i=0; i<f.size(); i++)
	{
	  memset(coor,0,sizeof(int)*f.signature());
	  int linind = i;
	    
	  f.linind_2_subind(linind, coor);

	  int min_k1 = PCM_k_cluster_type::subtract(coor[2], PCM_k_cluster_type::get_k_0_index());
	  int min_k2 = PCM_k_cluster_type::subtract(coor[6], PCM_k_cluster_type::get_k_0_index());

	  std::complex<double> c = (            f(coor[0], coor[1], coor[2],       coor[3], coor[4], coor[5], coor[6],      coor[7])
				    + std::conj(f(coor[0], coor[1], min_k1,  N_w - coor[3], coor[4], coor[5], min_k2, N_w - coor[7])))/2.;
	    
// 	  if(true){
// 	    difference(          f(coor[0], coor[1], coor[2],       coor[3], coor[4], coor[5], coor[6],       coor[7]),
// 		       std::conj(f(coor[0], coor[1], min_k1 , N_w - coor[3], coor[4], coor[5], min_k2 , N_w - coor[7])), 
// 			      __FUNCTION__);
// 	  }

	  f(coor[0], coor[1], coor[2],       coor[3], coor[4], coor[5], coor[6],       coor[7]) = c;           
	  f(coor[0], coor[1], min_k1 , N_w - coor[3], coor[4], coor[5], min_k2 , N_w - coor[7]) = std::conj(c);
	}
	
      delete [] coor;
    }

//   if(true)
//     symmetrize_over_Hamiltonians(f, H_symmetry, do_diff);
}
*/







// template<typename scalartype, typename f_dmn>
// void symmetrize::symmetrize_over_time_domain(function<scalartype, f_dmn>& f, bool do_diff)
// {
//   //   cout << __FUNCTION__ << endl;

//   typedef typename f_dmn::this_type type_list_dmns;
//   const static int t_index = IndexOf<type_list_dmns, time_domain_type>::value;
//   assert(t_index>=0 && t_index<f.signature());

//   int* coor_1 = new int[f.signature()];
//   int* coor_2 = new int[f.signature()];

//   for(int i=0; i<f.size(); i++)
//     {
//       memset(coor_1,0,sizeof(int)*f.signature());
//       int linind_1 = i;
//       f.linind_2_subind(linind_1, coor_1);

//       memset(coor_2,0,sizeof(int)*f.signature());
//       int linind_2 = i;
//       f.linind_2_subind(linind_2, coor_2);

//       if(coor_1[t_index]<t::dmn_size()/2)
// 	coor_2[t_index] = coor_1[t_index]+t::dmn_size()/2;
//       else
// 	coor_2[t_index] = coor_1[t_index]-t::dmn_size()/2;

//       if(do_diff){
// 	assert(difference(f(coor_1),  (f(coor_1) - f(coor_2))/2., __FUNCTION__));
//       }

//       f(coor_1) =  (f(coor_1) - f(coor_2))/2.;
//       f(coor_2) = -(f(coor_1) - f(coor_2))/2.;
//     }

//   delete [] coor_1;
//   delete [] coor_2;
// }
	
// template<typename scalartype, typename f_dmn_0, typename f_dmn_1>
// void symmetrize::symmetrize_over_frequency_domain(function<scalartype, dmn_4<nu, nu, f_dmn_0, f_dmn_1> >& f, bool do_diff)
// {
//   // H(i w) = conj(transpose(H(-i w)))
  
//   scalartype f_val_1, f_val_2;
//   if(IS_EQUAL_TYPE<f_dmn_1, w>::check)
//     {
//       for(int w_i=0; w_i<w::dmn_size(); w_i++){
// 	for(int f_i=0; f_i<f_dmn_0::dmn_size(); f_i++){
	  
// 	  for(int i=0; i<b::dmn_size()*s::dmn_size(); i++){
// 	    for(int j=0; j<b::dmn_size()*s::dmn_size(); j++){
	      
// 	      f_val_1 = f(i,j,f_i,                w_i);
// 	      f_val_2 = f(j,i,f_i,w::dmn_size()-1-w_i);
	      
// 	      if(do_diff){
// 		assert(difference(f_val_1, conjugate(f_val_2), __FUNCTION__));
// 	      }
	      
// 	      f(i,j,f_i,                w_i) =           (f_val_1 + conjugate(f_val_2)) /2.;
// 	      f(j,i,f_i,w::dmn_size()-1-w_i) = conjugate((f_val_1 + conjugate(f_val_2)))/2.;
// 	    }
// 	  }
	  
// 	}
//       }
//     }
// }
						  
// template<typename scalartype, typename f_dmn>
// void symmetrize::symmetrize_over_r_DCA(function<scalartype, f_dmn>& f, bool do_diff)
// {
//   //   cout << __FUNCTION__ << endl;
  
//   typedef typename f_dmn::this_type type_list_dmns;
//   const static int r_index = IndexOf<type_list_dmns, DCA_r_cluster_type>::value;
//   assert(r_index>=0 && r_index<f.signature());

//   int* coor = new int[f.signature()];
//   function<scalartype, f_dmn> f_new;
      
//   for(int i=0; i<f.size(); i++)
//     f_new(i) = scalartype(0.);

//   for(int i=0; i<f.size(); i++)
//     {
//       memset(coor,0,sizeof(int)*f.signature());
//       int linind = i;

//       f.linind_2_subind(linind, coor);

//       double weight = DCA_r_cluster_type::get_weight(coor[r_index]);
//       coor[r_index] = DCA_r_cluster_type::convert_r_to_r_irr(coor[r_index]);

//       f_new(coor) += f(i)/weight;
//     }
      
//   for(int i=0; i<f.size(); i++)
//     {
//       memset(coor,0,sizeof(int)*f.signature());
//       int linind = i;
	  
//       f.linind_2_subind(linind, coor);

//       coor[r_index] = DCA_r_cluster_type::convert_r_to_r_irr(coor[r_index]);

//       if(do_diff){
// 	difference(f(i), f_new(coor), __FUNCTION__);
//       }

//       f(i) = f_new(coor);
//     }

//   delete [] coor;
// }

// template<typename scalartype, typename f_dmn>
// void symmetrize::symmetrize_over_k_DCA(function<scalartype, f_dmn>& f, bool do_diff)
// {
//   //   cout << __FUNCTION__ << endl;

//   typedef typename f_dmn::this_type type_list_dmns;
//   const static int k_index = IndexOf<type_list_dmns, DCA_k_cluster_type>::value;
//   assert(k_index>=0 && k_index<f.signature());

//   int* coor = new int[f.signature()];
//   function<scalartype, f_dmn> f_new;
      
//   for(int i=0; i<f.size(); i++)
//     f_new(i) = scalartype(0.);

//   for(int i=0; i<f.size(); i++)
//     {
//       memset(coor,0,sizeof(int)*f.signature());
//       int linind = i;

//       f.linind_2_subind(linind, coor);

//       double weight = DCA_k_cluster_type::get_weight(coor[k_index]);
//       coor[k_index] = DCA_k_cluster_type::convert_k_to_k_irr(coor[k_index]);

//       f_new(coor) += f(i)/weight;
//     }
      
//   for(int i=0; i<f.size(); i++)
//     {
//       memset(coor,0,sizeof(int)*f.signature());
//       int linind = i;
	  
//       f.linind_2_subind(linind, coor);

//       coor[k_index] = DCA_k_cluster_type::convert_k_to_k_irr(coor[k_index]);

//       if(do_diff){
// 	difference(f(i), f_new(coor), __FUNCTION__);
//       }

//       f(i) = f_new(coor);
//     }

//   delete [] coor;
// }

// int symmetrize::get_number_of_symmetries(function<int, nu_nu>&  H_symmetry)
// {
//   int Nb_symmetries=-1;
  
//   for(int l=0; l<H_symmetry.size(); l++){
//     //cout << H_symmetry(l) << endl;

//     if(Nb_symmetries<H_symmetry(l))
//       Nb_symmetries = H_symmetry(l);
//   }
//   return Nb_symmetries;
// }

// template<typename scalartype, typename f_dmn_0, typename f_dmn_1 >
// void symmetrize::symmetrize_over_Hamiltonians(function<scalartype, dmn_4<nu, nu, f_dmn_0, f_dmn_1> >& f, 
// 					      function<int       , nu_nu>&                            H_symmetry, 
// 					      bool                                                    do_diff)
// {
//   //   cout << __FUNCTION__ << endl;

//   static int matrix_dim  = b::dmn_size()*s::dmn_size();
//   static int Nb_symmetries = symmetrize::get_number_of_symmetries(H_symmetry);

//   for(int ind_0=0; ind_0<f_dmn_0::dmn_size(); ind_0++){
//     for(int ind_1=0; ind_1<f_dmn_1::dmn_size(); ind_1++){
      
//       // H-symmetries
//       for(int l=0; l<Nb_symmetries+1; l++)
// 	{
// 	  scalartype result = 0;
// 	  double     index  = 0;
	  
// 	  for(int i=0; i<matrix_dim; i++){
// 	    for(int j=0; j<matrix_dim; j++){
// 	      if(l == H_symmetry(i,j)){
// 		result += f(i,j,ind_0,ind_1);
// 		index  += 1;
// 	      }
// 	    }
// 	  }

// 	  assert(index > 1);
// 	  result /= index;
	  
// 	  for(int i=0; i<matrix_dim; i++)
// 	    for(int j=0; j<matrix_dim; j++)
// 	      if(l == H_symmetry(i,j))
// 		f(i,j,ind_0,ind_1) = result;
// 	}
	  
//       // spin-symmetry ... --> G_(e_UP, e_DN) == G_(e_DN, e_UP) == 0 !!
//       for(int i=0; i<b::dmn_size(); i++){
// 	for(int j=0; j<b::dmn_size(); j++){
// 	  f(i,0,j,1,ind_0,ind_1) = 0;
// 	  f(i,1,j,0,ind_0,ind_1) = 0;
// 	}
//       }
//     }
//   }
// }

// template<typename scalartype>
// void symmetrize::symmetrize_over_Hamiltonians(function<scalartype, nu >&   f, 
// 					      function<int       , nu_nu>& H_symmetry, 
// 					      bool                         do_diff)
// {
//   //   cout << __FUNCTION__ << endl;

//   static int matrix_dim  = b::dmn_size()*s::dmn_size();
//   static int Nb_symmetries = symmetrize::get_number_of_symmetries(H_symmetry);
 
//   // H-symmetries
//   for(int l=0; l<Nb_symmetries+1; l++)
//     {
//       double  result = 0;
//       double  index  = 0;
	  
//       for(int i=0; i<matrix_dim; i++){
// 	if(l == H_symmetry(i,i)){
// 	  result += f(i);
// 	  index  += 1;	    
// 	}
//       }

//       assert(index > 1);
//       result /= index;
	  
//       for(int i=0; i<matrix_dim; i++)
// 	if(l == H_symmetry(i,i))
// 	  f(i) = result;
//     }
// }

// void symmetrize::permute_first_pair_band_indices(int l,
// 				      int    l1, int    l3,
// 				      int& p_l1, int& p_l3)
// {
//   std::pair<int,int> initial_state = model::get_orbital_permutations()[l].first;
//   std::pair<int,int> final_state   = model::get_orbital_permutations()[l].second;

//   p_l1=l1;
//   if(l1==initial_state.first)  p_l1 = final_state.first;
//   if(l1==initial_state.second) p_l1 = final_state.second;

//   p_l3=l3;
//   if(l3==initial_state.first)  p_l3 = final_state.first;
//   if(l3==initial_state.second) p_l3 = final_state.second;
// }

// void symmetrize::permute_second_pair_band_indices(int l,
// 						  int    l2, int    l4,
// 						  int& p_l2, int& p_l4)
// {
//   std::pair<int,int> initial_state = model::get_orbital_permutations()[l].first;
//   std::pair<int,int> final_state   = model::get_orbital_permutations()[l].second;

//   p_l2=l2;
//   if(l2==initial_state.first)  p_l2 = final_state.first;
//   if(l2==initial_state.second) p_l2 = final_state.second;

//   p_l4=l4;
//   if(l4==initial_state.first)  p_l4 = final_state.first;
//   if(l4==initial_state.second) p_l4 = final_state.second;
// }
