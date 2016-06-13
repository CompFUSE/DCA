//-*-C++-*-

#ifndef ANALYSIS_COMPUTE_REDUCIBLE_VERTEX_H
#define ANALYSIS_COMPUTE_REDUCIBLE_VERTEX_H

namespace dca {

  /*! 
   * \author peter staar
   */
  template<class parameter_type, class MOMS_type>
  class analysis<parameter_type, MOMS_type, ANALYSIS_COMPUTE_REDUCIBLE_VERTEX>
  {
#include "type_definitions.h"

    //typedef dmn_0<brillouin_zone_path_domain<SQUARE_2D_LATTICE> > k_dmn_cut_type;
    typedef dmn_0<brillouin_zone_path_domain<FERMI_SURFACE_SQUARE_2D_LATTICE> > k_dmn_cut_type;

  public:
 
    analysis(parameter_type& parameters, MOMS_type& MOMS);
    ~analysis();

    template<class stream_type>
    void to_JSON(stream_type& ss);

    void execute();

    void execute_all();
 
  private:

    void initialize();

    void apply_symmetries();

    void read_this_G4(int q_ind, FUNC_LIB::function<std::complex<double>, dmn_8<b,b,b,b,k_DCA,k_DCA,w_VERTEX,w_VERTEX> >& G4_q);

    void set_into_full_G4(int q_ind, 
			  FUNC_LIB::function<std::complex<double>, dmn_8<b,b,b,b,k_DCA,k_DCA,w_VERTEX,w_VERTEX> >&         G4_q,
			  FUNC_LIB::function<std::complex<double>, dmn_3<b_b_k_DCA_w_VERTEX, b_b_k_DCA_w_VERTEX, k_DCA> >& G4_full);


    void load_G4_b_k_w__b_k_w();
    void load_G4_0_b_k_w__b_k_w();

    void compute_Gamma_b_k_w__b_k_w();

    void compute_full_chi_0_b_k_w__b_k_w();

    void compute_reducible_vertex();

    void test(FUNC_LIB::function<std::complex<double>, dmn_3<b_b_k_DCA_w_VERTEX, b_b_k_DCA_w_VERTEX, k_DCA> >& G4_full);

  private:    

    parameter_type&  parameters;
    MOMS_type&       MOMS;

    std::string      data_filename;

    diagrammatic_symmetries<parameter_type> diagrammatic_symmetries_obj;

  public:

    FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> G4;

    FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> G4_0_b_k_w__b_k_w;

    FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> full_chi_0;
    FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> inverted_full_chi_0;

    FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> reducible_vertex;

    //FUNC_LIB::function<std::complex<double>, dmn_3<b_b_k_DCA, b_b_k_DCA, k_DCA> >                      full_reducible_vertex;
    FUNC_LIB::function<std::complex<double>, dmn_3<b_b_k_DCA_w_VERTEX, b_b_k_DCA_w_VERTEX, k_DCA> >    full_reducible_vertex;
    
//     make_G4_matrix                <parameter_type, MOMS_type> make_G4_obj;
    make_G4_0_matrix              <parameter_type, MOMS_type> make_G4_0_obj;

    b_b_k_DCA_w_VERTEX b_b_k_DCA_w_VERTEX_domain;
//     b_b_k_DCA_w_VERTEX          b_b_k_DCA_w_VERTEX_domain;

//     std::vector<int>  corresponding_extended_index;
//     std::vector<bool> is_compact_frequency;
  };

  template<class parameter_type, class MOMS_type>
  analysis<parameter_type, MOMS_type, ANALYSIS_COMPUTE_REDUCIBLE_VERTEX>::analysis(parameter_type& parameters_in, 
										   MOMS_type& MOMS_in):
    parameters(parameters_in),
    MOMS(MOMS_in),

    data_filename(""),

    diagrammatic_symmetries_obj(parameters),

    G4("G4"),

    full_chi_0("full_chi_0"),
    inverted_full_chi_0("inverted_full_chi_0"),
    
    reducible_vertex("reducible_vertex"),

    full_reducible_vertex("full_reducible_vertex"),
    //full_reducible_vertex_2("full_reducible_vertex_2"),

//     make_G4_obj  (parameters_in, MOMS_in),
    make_G4_0_obj(parameters_in, MOMS_in)

//     corresponding_extended_index(w_VERTEX::dmn_size(),-1),
//     is_compact_frequency(w_VERTEX::dmn_size(),false)
  {
    data_filename                     = parameters.get_output_file_name();
    parameters.get_output_file_name() = parameters.get_output_susceptibilities_file_name();

//     initialize();
  }

  template<class parameter_type, class MOMS_type>
  analysis<parameter_type, MOMS_type, ANALYSIS_COMPUTE_REDUCIBLE_VERTEX>::~analysis()
  {}

//   template<class parameter_type, class MOMS_type>
//   void analysis<parameter_type, MOMS_type, ANALYSIS_COMPUTE_REDUCIBLE_VERTEX>::initialize()
//   {
// //     for(int i=0; i<w_VERTEX::dmn_size(); i++)
// //       for(int j=0; j<w_VERTEX::dmn_size(); j++)
// // 	if(std::fabs(w_VERTEX::parameter_type::get_elements()[i]-w_VERTEX::parameter_type::get_elements()[j])<1.e-6)
// // 	  corresponding_extended_index[i] = j;

// //     for(int j=0; j<w_VERTEX::dmn_size(); j++)
// //       for(int i=0; i<w_VERTEX::dmn_size(); i++)
// //       	if(std::fabs(w_VERTEX::parameter_type::get_elements()[i]-w_VERTEX::parameter_type::get_elements()[j])<1.e-6)
// // 	  is_compact_frequency[j] = true;
//   }

  template<class parameter_type, class MOMS_type>
  template<class stream_type>
  void analysis<parameter_type, MOMS_type, ANALYSIS_COMPUTE_REDUCIBLE_VERTEX>::to_JSON(stream_type& ss)
  {
    //reducible_vertex.to_JSON(ss);

    if(parameters.get_output_file_name() == parameters.get_directory()+"/full_reducible_vertex.json")
      {
	FUNC_LIB::function<std::complex<double>, dmn_7<b,b,k_DCA,  b,b,k_DCA, k_DCA> > reducible_vertex_small("full_reducible_vertex");
	
	for(int b1=0; b1<b::dmn_size(); b1++)
	  for(int b2=0; b2<b::dmn_size(); b2++)
	    for(int k1=0; k1<k_DCA::dmn_size(); k1++)
	      
	      for(int b3=0; b3<b::dmn_size(); b3++)
		for(int b4=0; b4<b::dmn_size(); b4++)
		  for(int k2=0; k2<k_DCA::dmn_size(); k2++)

		    for(int k3=0; k3<k_DCA::dmn_size(); k3++)
		      reducible_vertex_small(b1,b2,k1,  b3,b4,k2, k3) = full_reducible_vertex(b1,b2,k1,w_VERTEX::dmn_size()/2,  b3,b4,k2,w_VERTEX::dmn_size()/2, k3); 

	reducible_vertex_small.to_JSON(ss);
      }
    else
      {
	FUNC_LIB::function<std::complex<double>, dmn_6<b,b,k_DCA,  b,b,k_DCA> > reducible_vertex_small("reducible_vertex");
	
	for(int b1=0; b1<b::dmn_size(); b1++){
	  for(int b2=0; b2<b::dmn_size(); b2++){
	    for(int k1=0; k1<k_DCA::dmn_size(); k1++){
	      
	      for(int b3=0; b3<b::dmn_size(); b3++){
		for(int b4=0; b4<b::dmn_size(); b4++){
		  for(int k2=0; k2<k_DCA::dmn_size(); k2++){
		    reducible_vertex_small(b1,b2,k1,  b3,b4,k2) = reducible_vertex(b1,b2,k1,w_VERTEX::dmn_size()/2,  b3,b4,k2,w_VERTEX::dmn_size()/2); 
		  }
		}
	      }
	    }
	  }
	}

	reducible_vertex_small.to_JSON(ss);
      }
  }

  template<class parameter_type, class MOMS_type>
  void analysis<parameter_type, MOMS_type, ANALYSIS_COMPUTE_REDUCIBLE_VERTEX>::execute()
  {
    cout << __FUNCTION__ << endl;

    apply_symmetries();

    double renorm = 1./(parameters.get_beta()*k_DCA::dmn_size());
    
    {
      int* coor_1 = new int[     G4        .signature()];
      int* coor_2 = new int[MOMS.G4_k_k_w_w.signature()];
      
      for(int i=0; i<MOMS.G4_k_k_w_w.size(); i++)
	{
	  MOMS.G4_k_k_w_w.linind_2_subind(i, coor_2);
	  
	  coor_1[0] = coor_2[0];
	  coor_1[1] = coor_2[1];
	  coor_1[2] = coor_2[4];//k_1
	  coor_1[3] = coor_2[6];//w_1
	  coor_1[4] = coor_2[2];
	  coor_1[5] = coor_2[3];
	  coor_1[6] = coor_2[5];//k_2
	  coor_1[7] = coor_2[7];//w_2
	  
	  G4(coor_1[0], coor_1[1], coor_1[2], coor_1[3], coor_1[4], coor_1[5], coor_1[6], coor_1[7])
	    = MOMS.G4_k_k_w_w(coor_2[0], coor_2[1], coor_2[2], coor_2[3], coor_2[4], coor_2[5], coor_2[6], coor_2[7]);
	}
      
      delete [] coor_1;
      delete [] coor_2;

//       diagrammatic_symmetries_obj.execute(G4);
//       symmetrize::execute(G4, parameters.get_q_vector());
    }

    {
      make_G4_0_obj.execute(*this);
      for(int l=0; l<full_chi_0.size(); l++)
	full_chi_0(l) = G4_0_b_k_w__b_k_w(l);

//       diagrammatic_symmetries_obj.execute(full_chi_0);
//       symmetrize::execute(full_chi_0, parameters.get_q_vector());
    }

    G4         *= renorm;
    full_chi_0 *= renorm;

    compute_reducible_vertex();
  }

  template<class parameter_type, class MOMS_type>
  void analysis<parameter_type, MOMS_type, ANALYSIS_COMPUTE_REDUCIBLE_VERTEX>::execute_all()
  {
    cout << __FUNCTION__ << endl;

    parameters.get_output_file_name() = parameters.get_directory()+"/full_reducible_vertex.json";

    apply_symmetries();

    FUNC_LIB::function<std::complex<double>, dmn_8<b,b,b,b,k_DCA,k_DCA,w_VERTEX,w_VERTEX>         > G4_k_k_w_w("G4_k_k_w_w");
    FUNC_LIB::function<std::complex<double>, dmn_3<b_b_k_DCA_w_VERTEX, b_b_k_DCA_w_VERTEX, k_DCA> > G4_full   ("G4_k_w_k_w");
    //FUNC_LIB::function<std::complex<double>, dmn_3<b_b_k_DCA_w_VERTEX, b_b_k_DCA_w_VERTEX, k_DCA> > G4_0_full ("G4_0_k_w_k_w");

    for(int q_ind=0; q_ind<k_DCA::dmn_size(); q_ind++)
      {
	parameters.get_q_channel() = q_ind;
	parameters.get_q_vector()  = k_DCA::get_elements()[q_ind];

	cout << parameters.get_q_channel() << "\t";
	VECTOR_OPERATIONS::PRINT(parameters.get_q_vector());
	cout << "\t";

	read_this_G4(q_ind, G4_k_k_w_w);

	set_into_full_G4(q_ind, G4_k_k_w_w, G4_full);

// 	make_G4_0_obj.execute(*this);
// 	for(int j=0; j<sqrt(full_chi_0.size()); ++j)
// 	  for(int i=0; i<sqrt(full_chi_0.size()); ++i)
// 	    G4_0_full(i,j,q_ind) = G4_0_b_k_w__b_k_w(i,j);
      }

    symmetrize::                execute(G4_full);
    //diagrammatic_symmetries_obj.execute(G4_full);

    double renorm = 1./(parameters.get_beta()*k_DCA::dmn_size());

    int N = square(b::dmn_size())*k_DCA::dmn_size()*w_VERTEX::dmn_size();    
    for(int q_ind=0; q_ind<k_DCA::dmn_size(); q_ind++)
      {
	parameters.get_q_channel() = q_ind;
	parameters.get_q_vector()  = k_DCA::get_elements()[q_ind];

	cout << parameters.get_q_channel() << "\t";
	VECTOR_OPERATIONS::PRINT(parameters.get_q_vector());
	cout << "\n";

	{
	  make_G4_0_obj.execute(*this);

	  for(int l=0; l<full_chi_0.size(); l++)
	    full_chi_0(l) = G4_0_b_k_w__b_k_w(l);
	}

	memcpy(&G4(0), &G4_full(0,0,q_ind), sizeof(std::complex<double>)*square(N));

	G4         *= renorm;
	full_chi_0 *= renorm;

	compute_reducible_vertex();

	memcpy(&full_reducible_vertex(0,0,q_ind), &reducible_vertex(0), sizeof(std::complex<double>)*square(N));
      }

    symmetrize::                execute(full_reducible_vertex);
    //diagrammatic_symmetries_obj.execute(full_reducible_vertex);
  }


  template<class parameter_type, class MOMS_type>
  void analysis<parameter_type, MOMS_type, ANALYSIS_COMPUTE_REDUCIBLE_VERTEX>::apply_symmetries()
  {
    //cout << __FUNCTION__ << endl;

    symmetrize::execute(MOMS.Sigma, MOMS.H_symmetry);

    symmetrize::execute(MOMS.G_k_w, MOMS.H_symmetry);
  }

  template<class parameter_type, class MOMS_type>
  void analysis<parameter_type, MOMS_type, ANALYSIS_COMPUTE_REDUCIBLE_VERTEX>::read_this_G4(int q_ind, 
											    FUNC_LIB::function<std::complex<double>, dmn_8<b,b,b,b,k_DCA,k_DCA,w_VERTEX,w_VERTEX> >& G4_q)
  {
    std::string filename = data_filename;

    size_t loc = filename.find( string("q=0") );
    filename.erase(loc, 3);

    std::ostringstream os;
    os << "q=";
    os << q_ind;
    std::string q_str(os.str());

    filename.insert(loc, q_str);

    cout << filename << endl;

    {
//       if(! ifstream(&(filename)[0])){
// 	cout << "\t no file with name : " << &(filename)[0] << endl;
// 	throw std::logic_error(__FUNCTION__);
//       }

//       FILE* file_ptr = fopen(&(filename)[0], "r");
//       if(file_ptr == NULL){
// 	cout << "\t no file with name : " << &(filename)[0] << endl;
// 	throw std::logic_error(__FUNCTION__);
//       }
//       fclose(file_ptr);
    }

    {
      std::fstream input_file(&(filename)[0], ios::in);
      G4_q.from_JSON(input_file);
      input_file.close();
    }
  }
    
  template<class parameter_type, class MOMS_type>
  void analysis<parameter_type, MOMS_type, ANALYSIS_COMPUTE_REDUCIBLE_VERTEX>::set_into_full_G4(int q_ind, 
												FUNC_LIB::function<std::complex<double>, dmn_8<b,b,b,b,k_DCA,k_DCA,w_VERTEX,w_VERTEX> >&                          G4_q,
												FUNC_LIB::function<std::complex<double>, dmn_3<b_b_k_DCA_w_VERTEX,b_b_k_DCA_w_VERTEX, k_DCA> >& G4_full)
  {
    int* coor_1 = new int[G4_full.signature()];
    int* coor_2 = new int[G4_q   .signature()];
    
    for(int i=0; i<G4_q.size(); i++)
      {
	G4_q.linind_2_subind(i, coor_2);
	
	coor_1[0] = coor_2[0];
	coor_1[1] = coor_2[1];
	coor_1[2] = coor_2[4];//k_1
	coor_1[3] = coor_2[6];//w_1
	coor_1[4] = coor_2[2];
	coor_1[5] = coor_2[3];
	coor_1[6] = coor_2[5];//k_2
	coor_1[7] = coor_2[7];//w_2
	
	G4_full(coor_1[0], coor_1[1], coor_1[2], coor_1[3], coor_1[4], coor_1[5], coor_1[6], coor_1[7], q_ind)
	  = G4_q(coor_2[0], coor_2[1], coor_2[2], coor_2[3], coor_2[4], coor_2[5], coor_2[6], coor_2[7]);
      }
    
    delete [] coor_1;
    delete [] coor_2;
  }

  template<class parameter_type, class MOMS_type>
  void analysis<parameter_type, MOMS_type, ANALYSIS_COMPUTE_REDUCIBLE_VERTEX>::compute_reducible_vertex()
  {
    cout << __FUNCTION__ << endl;
    cout << scientific;
    cout.precision(6);

    int N = square(b::dmn_size())*k_DCA::dmn_size()*w_VERTEX::dmn_size();    

    {
      invert_plan<std::complex<double> > invert_pln(N);
      
      memcpy(&invert_pln.Matrix[0], &full_chi_0(0), sizeof(std::complex<double>)*N*N);
      invert_pln.execute_plan();
      memcpy(&inverted_full_chi_0(0), invert_pln.inverted_matrix, sizeof(std::complex<double>)*N*N);

      //cout << "symmetrize chi_0_b_k_w__b_k_w" << endl;
      //symmetrize::execute(inverted_full_chi_0, MOMS.H_symmetry, true);
    }

    for(int l=0; l<G4.size(); l++)
      G4(l) -= full_chi_0(l);

    FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> G4_min_full_chi_0__times_inverted_full_chi_0("tmp");

    {
      gemm_plan<std::complex<double> > gemm_pln(N);
      
      gemm_pln.A = &G4(0);
      gemm_pln.B = &inverted_full_chi_0(0);
      gemm_pln.C = &G4_min_full_chi_0__times_inverted_full_chi_0(0);
      
      gemm_pln.execute_plan();
    }

    {
      gemm_plan<std::complex<double> > gemm_pln(N);
      
      gemm_pln.A = &inverted_full_chi_0(0);
      gemm_pln.B = &G4_min_full_chi_0__times_inverted_full_chi_0(0);
      gemm_pln.C = &reducible_vertex(0);
      
      gemm_pln.execute_plan();
    }

  }

  template<class parameter_type, class MOMS_type>
  void analysis<parameter_type, MOMS_type, ANALYSIS_COMPUTE_REDUCIBLE_VERTEX>::test(FUNC_LIB::function<std::complex<double>, dmn_3<b_b_k_DCA_w_VERTEX, b_b_k_DCA_w_VERTEX, k_DCA> >& G4_full)
  {
    cout << scientific;
    cout.precision(6);

    if(parameters.get_w_channel()==-1)
      for(int k3=0; k3<k_DCA::dmn_size(); k3++)

	for(int b1=0; b1<b::dmn_size(); b1++)
	  for(int b2=0; b2<b::dmn_size(); b2++)
	    for(int k1=0; k1<k_DCA::dmn_size(); k1++)
	      
	      for(int b3=0; b3<b::dmn_size(); b3++){
		for(int b4=0; b4<b::dmn_size(); b4++)
		  //for(int k2=0; k2<k_DCA::dmn_size(); k2++)
		    {
		      
		      int k2=k1;

		      int k1_plus_k2      = k_DCA::parameter_type::add(k1,k2);
		      int q_min_k1_min_k2 = k_DCA::parameter_type::subtract(k1_plus_k2, k3);
		    
		      cout << real(G4_full(b1,b2,k1,w_VERTEX::dmn_size()/2,  b3,b4,k2,w_VERTEX::dmn_size()/2, q_min_k1_min_k2)) << "\t"; 
		    }
		cout << "\n";
	      }
    
    if(parameters.get_w_channel()==0)
      for(int k3=0; k3<k_DCA::dmn_size(); k3++)

	for(int b1=0; b1<b::dmn_size(); b1++)
	  for(int b2=0; b2<b::dmn_size(); b2++)
	  for(int k1=0; k1<k_DCA::dmn_size(); k1++)
	    
	    for(int b3=0; b3<b::dmn_size(); b3++){
	      for(int b4=0; b4<b::dmn_size(); b4++)
		//for(int k2=0; k2<k_DCA::dmn_size(); k2++)
		{
		  int k2=k1;

		  cout << real(G4_full(b1,b2,k1,w_VERTEX::dmn_size()/2,  b3,b4,k2,w_VERTEX::dmn_size()/2, k3)) << "\t"; 
		}
	      cout << "\n";
	    }


  }
}

#endif
