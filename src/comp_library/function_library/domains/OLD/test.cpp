
#include "Parameter.h"



#include "cluster.h"
#include "r_cluster.h"
#include "k_cluster.h"

#include "D4.h"



#include "time_domain.h"
#include "frequency_domain.h"



#include "function.h"

void test_construction_and_printing()
{
  cout << endl << endl << __FUNCTION__ << endl << endl;


  /*typedef cluster<Parameters> cluster_type;
  cluster_type::initializer::initialize(Parameters::A,Parameters::B);
  
  //reduce< cluster_type >::execute();

  typedef dmn_0<r_cluster<FULL, cluster_type> > r;
  typedef dmn_0<k_cluster<FULL, cluster_type> > k;

  typedef dmn_0<time_domain<FULL, Parameters> > t;
  typedef dmn_0<frequency_domain<FULL, Parameters> > w;

  typedef dmn_2<r,t> r_t;

  function<std::vector<int>, r_t> f;

  f.print_fingerprint();

  std::vector<int> v(3,0);
  f(3,15) = v;*/

}


void test_slize_and_distribute()
{
  cout << endl << endl << __FUNCTION__ << endl << endl;

  /*typedef cluster<Parameters> cluster_type;

  cluster_type::initializer::initialize(Parameters::A,Parameters::B);
  //reduce< cluster_type >::execute();
  cluster_type::print();

  typedef dmn_0<r_cluster<FULL, cluster_type> > r;
  typedef dmn_0<k_cluster<FULL, cluster_type> > k;

  typedef dmn_0<time_domain<FULL, Parameters> > t;
  typedef dmn_0<frequency_domain<FULL, Parameters> > w;

  typedef dmn_3<r,r,t> r_r_t;

  function<double, r_r_t> f;

  f.print_fingerprint();

  
  double data[t::dmn_size()] ;

  for(int i=0; i<t::dmn_size(); i++)
    {
      data[i] = -1.;
      f(2,3,i) = double(i);
    }

  int coor[3] = {2,3,0};
  f.slice(2, &coor[0], &data[0]);


  for(int i=0; i<t::dmn_size(); i++)
    {
      cout << i << " : " << data[i] << "\t" <<  f(2,3,i) << endl;
    }
  


  double data[r::dmn_size()*r::dmn_size()] ;

  for(int i=0; i<r::dmn_size(); i++)
    for(int j=0; j<r::dmn_size(); j++)
      f(i,j,3) = double(i);

  int coor[3] = {0,0,3};
  f.slice(0,1, &coor[0], &data[0]);

  for(int i=0; i<r::dmn_size(); i++){
    for(int j=0; j<r::dmn_size(); j++){
      cout << f(i,j,3) << "\t";
    }
    cout << endl;
    }*/



}

#include "LDA_HUBBARD.h"
#include "DCA_HUBBARD.h"
#include "MultiOrbitalMultiSiteStructure.h"

void do_BIT_test()
{
  typedef tight_binding_models<cluster_8A_2D> model;

  //LDA_HUBBARD<model>::print_fingerprint();
  //DCA_HUBBARD<model>::print_fingerprint();

  typedef LDA_HUBBARD<model> lda_cluster_type;
  typedef DCA_HUBBARD<model> dca_cluster_type;

  MultiOrbitalMultiSiteStructure<dca_cluster_type, lda_cluster_type>::initialize_all_domains();

  MultiOrbitalMultiSiteStructure<dca_cluster_type, lda_cluster_type>::do_BIT_FT_k_2_r();

  MultiOrbitalMultiSiteStructure<dca_cluster_type, lda_cluster_type>::do_BIT_FT_w_2_t();
}

void test_FT_r_2_k()
{
  /*LDA_HUBBARD::print_fingerprint();
  DCA_HUBBARD::print_fingerprint();

  MultiOrbitalMultiSiteStructure<DCA_HUBBARD, LDA_HUBBARD>::initialize_all_domains();
  MultiOrbitalMultiSiteStructure<DCA_HUBBARD, LDA_HUBBARD>::do_BIT_FT_k_2_r();*/


}


void test_FT_t_2_w()
{
  
}

#include "LDA_DFT.h"
#include "DCA_DFT.h"
#include "LDA_parser.h"

#include "MultiOrbitalMultiSiteStructure.h"

void test_DFT()
{
  LDA_parser::read<LDA_DFT>("/Users/peterstaar/Documents/C++TestPrograms/FUNCTION_GENERIC_8/WANN_H_4x4x4.OUT");
  LDA_DFT::print_fingerprint();

  DCA_DFT::initialize();
  DCA_DFT::print_fingerprint();

  MultiOrbitalMultiSiteStructure<DCA_DFT, LDA_DFT>::initialize_all_domains();
  MultiOrbitalMultiSiteStructure<DCA_DFT, LDA_DFT> MOMS;

  MOMS.high_symmetry_line_dispersion_3D();
}

#include "LDA_HUBBARD.h"
#include "DCA_HUBBARD.h"

void test_HUBBARD()
{
  typedef tight_binding_models<cluster_8A_2D> model;
  
  //LDA_HUBBARD<model>::print_fingerprint();
  //DCA_HUBBARD<model>::print_fingerprint();

  typedef LDA_HUBBARD<model> lda_cluster_type;
  typedef DCA_HUBBARD<model> dca_cluster_type;

  MultiOrbitalMultiSiteStructure<dca_cluster_type, lda_cluster_type>::initialize_all_domains();
  MultiOrbitalMultiSiteStructure<dca_cluster_type, lda_cluster_type> MOMS;

  

  //MOMS.high_symmetry_line_dispersion_2D();
}



#include <util.h>
#include <nfft3++.h>

void test_nfft()
{
  { // 1D
    nfft_plan p;
  
    //int N=20;
    int M=20;
    int N = M;

    //nfft_init_1d(&p,N,M);
    nfft_init(&p,1,&N,M);

    for(int i=0;i<M;i++)
      p.x[i] = -0.5 + double(i)/double(M);
	
    for(int i=0;i<M;i++)
      {
	p.f_hat[i][0] = 0.;
	p.f_hat[i][1] = 0.;
      }
  
    //p.f_hat[0][0] = 1.;
	
    for(int i=0;i<M;i++)
      {
	p.f[i][0] = sin(2*3.14159*p.x[i]);
	p.f[i][1] = 0.;
      }
	
    /** precompute psi, the entries of the matrix B */
    if(p.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&p);

    /** direct trafo and show the result */
    //ndft_trafo(&p);
	
    /** approx. adjoint and show the result */
    //ndft_adjoint(&p);

    /** approx. adjoint and show the result */
    nfft_adjoint(&p);
    //nfft_vpr_complex(p.f_hat,p.N_total,"adjoint nfft, vector f_hat");

    /** approx. trafo and show the result */
    nfft_trafo(&p);

    for(int i=0; i<M; i++)
      cout << p.x[i] << "\t" << p.f[i][0] << "\t" << p.f[i][1] << " <---> " << p.f_hat[i][0] << "\t" << p.f_hat[i][1]<< endl;

    nfft_finalize(&p);
  }
  
  { // 2D
    nfft_plan p;
  
    //int N=20;
    int M=25;
    int* N = new int[2];
    N[0] = 5;
    N[1] = 5;

    nfft_init(&p,2,N,M);

    for(int i=0;i<N[0];i++)
      {
	for(int j=0;j<N[1];j++)
	  {
	    p.x[2*(i+N[0]*j)+0] = -0.5 + double(i)/double(N[0]);
	    p.x[2*(i+N[0]*j)+1] = -0.5 + double(j)/double(N[1]);

	    p.f[i+N[0]*j][0]     = sin(2*3.14159*p.x[2*(i+N[0]*j)+0]);;
	    p.f[i+N[0]*j][1]     = 0.;

	    p.f_hat[i+N[0]*j][0] = 0.;
	    p.f_hat[i+N[0]*j][1] = 0.;
	  }
      }

    


    if(p.nfft_flags & PRE_ONE_PSI)
      nfft_precompute_one_psi(&p);


    //ndft_trafo(&p);
	

    //ndft_adjoint(&p);


    //nfft_adjoint(&p);
    //nfft_vpr_complex(p.f_hat,p.N_total,"adjoint nfft, vector f_hat");

    cout << "do 2D nfft" << endl;
    //nfft_trafo(&p);
    ndft_adjoint(&p);

    cout << "px_x :: " << endl;
    for(int i=0;i<N[0];i++)
      {
	for(int j=0;j<N[1];j++)
	  {
	    cout << p.x[2*(i+N[0]*j)+0] << "\t";
	  }
	cout << endl;
      }

    cout << "px_y :: " << endl;
    for(int i=0;i<N[0];i++)
      {
	for(int j=0;j<N[1];j++)
	  {
	    cout << p.x[2*(i+N[0]*j)+1] << "\t";
	  }
	cout << endl;
      }

    cout << "Re[p.f] :: " << endl;
    for(int i=0;i<N[0];i++)
      {
	for(int j=0;j<N[1];j++)
	  {
	    cout << p.f[i+N[0]*j][0] << "\t";
	  }
	cout << endl;
      }

    cout << "Im[p.f] :: " << endl;
    for(int i=0;i<N[0];i++)
      {
	for(int j=0;j<N[1];j++)
	  {
	    cout << p.f[i+N[0]*j][1] << "\t";
	  }
	cout << endl;
      }

    cout << "Re[p.f_hat] :: " << endl;
    for(int i=0;i<N[0];i++)
      {
	for(int j=0;j<N[1];j++)
	  {
	    cout << p.f_hat[i+N[0]*j][0] << "\t";
	  }
	cout << endl;
      }

    cout << "Im[p.f_hat] :: " << endl;
    for(int i=0;i<N[0];i++)
      {
	for(int j=0;j<N[1];j++)
	  {
	    cout << p.f_hat[i+N[0]*j][1] << "\t";
	  }
	cout << endl;
      }


    nfft_finalize(&p);
  }

  assert(false);


}
