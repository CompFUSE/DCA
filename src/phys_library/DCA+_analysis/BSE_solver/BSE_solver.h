//-*-C++-*-

#ifndef BSE_SOLVER_H
#define BSE_SOLVER_H
#include"phys_library/domain_types.hpp"
using namespace types;

namespace DCA
{
  /*!
   *  \author Peter Staar
   */
  using types::k_HOST_VERTEX;
  using types::w_VERTEX;
  template<class parameters_type, class MOMS_type>
  class BSE_solver
  {

    typedef double scalartype;

    const static int N_LAMBDAS = 10;
    typedef dmn_0<dmn<N_LAMBDAS, int> > lambda_dmn_type;

    const static int N_HARMONICS = 3;
    typedef dmn_0<dmn<N_HARMONICS, int> > harmonics_dmn_type;

    typedef dmn_0<brillouin_zone_path_domain<SQUARE_2D_LATTICE> > k_dmn_cut_type;

    typedef typename parameters_type::profiler_type    profiler_t;
    typedef typename parameters_type::concurrency_type concurrency_t;

    typedef dmn_4<b,b,k_DCA        , w_VERTEX> cluster_eigenvector_dmn_t;
    typedef dmn_4<b,b,k_HOST_VERTEX, w_VERTEX> lattice_eigenvector_dmn_t;

    typedef dmn_2<dmn_4<b,b,k_DCA        , w_VERTEX>, dmn_4<b,b,k_DCA        , w_VERTEX> >   DCA_DCA_matrix_dmn_t;
    typedef dmn_2<dmn_4<b,b,k_HOST_VERTEX, w_VERTEX>, dmn_4<b,b,k_DCA        , w_VERTEX> >  HOST_DCA_matrix_dmn_t;
    typedef dmn_2<dmn_4<b,b,k_DCA        , w_VERTEX>, dmn_4<b,b,k_HOST_VERTEX, w_VERTEX> >  DCA_HOST_matrix_dmn_t;
    typedef dmn_2<dmn_4<b,b,k_HOST_VERTEX, w_VERTEX>, dmn_4<b,b,k_HOST_VERTEX, w_VERTEX> > HOST_HOST_matrix_dmn_t;

    typedef dmn_2<dmn_4<b,b,k_DCA        , w_VERTEX>, dmn_4<b,b,k_DCA        , w_VERTEX> > DCA_matrix_dmn_t;
    typedef dmn_2<dmn_4<b,b,k_HOST_VERTEX, w_VERTEX>, dmn_4<b,b,k_HOST_VERTEX, w_VERTEX> > HOST_matrix_dmn_t;

  public:

    BSE_solver(parameters_type& parameters, MOMS_type& MOMS);
    ~BSE_solver();

    void write();

    template<IO::FORMAT DATA_FORMAT>
    void write(IO::writer<DATA_FORMAT>& reader);

    template<class stream_type>
    void to_JSON(stream_type& ss);

    void calculate_cuts();

    void calculate_susceptibilities_1();
    void calculate_susceptibilities_2();

    FUNC_LIB::function<std::complex<scalartype>, lambda_dmn_type>& get_leading_eigenvalues() { return BSE_lattice_solver_obj.get_leading_eigenvalues(); };

  private:

    void initialize_wave_functions();

    void apply_symmetries();

    void load_the_matrices();

    void compute_Gamma_cluster();

    void compute_chi_0();
    void compute_chi_0_2();

    void compute_Gamma_lattice();

    void diagonolize_Gamma_times_chi_0();

    void diagonolize_Gamma_times_chi_0_general();

    void diagonolize_Gamma_times_chi_0_symmetric_0();
    void diagonolize_Gamma_times_chi_0_symmetric_1();

    void compute_P_q_cluster();
    void compute_P_q_lattice();

    void find_harmonic_expansion();

    void write_on_shell();

  private:

    parameters_type& parameters;
    concurrency_t&   concurrency;

    MOMS_type&       MOMS;

    BSE_cluster_solver<parameters_type, MOMS_type> BSE_cluster_solver_obj;
    BSE_lattice_solver<parameters_type, MOMS_type> BSE_lattice_solver_obj;

    cluster_eigenvector_dmn_t cluster_eigenvector_dmn;
    lattice_eigenvector_dmn_t lattice_eigenvector_dmn;

    diagrammatic_symmetries<parameters_type> diagrammatic_symmetries_obj;

    FUNC_LIB::function<std::string             ,                      harmonics_dmn_type>   wave_functions_names;
    FUNC_LIB::function<std::complex<scalartype>, dmn_2<k_HOST_VERTEX, harmonics_dmn_type> > harmonics;

    FUNC_LIB::function<std::complex<double>, DCA_matrix_dmn_t> G4;
    FUNC_LIB::function<std::complex<double>, DCA_matrix_dmn_t> G4_0;

    FUNC_LIB::function<std::complex<scalartype>, DCA_matrix_dmn_t > Gamma_cluster;
    FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t> Gamma_lattice;

    FUNC_LIB::function<std::complex<double>, dmn_4<b_b, b_b, k_HOST_VERTEX, w_VERTEX> > chi_0_function;//("phi");

    FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t> chi;
    FUNC_LIB::function<std::complex<scalartype>, HOST_matrix_dmn_t> chi_0;

    FUNC_LIB::function<std::complex<scalartype>, dmn_0<dmn<1, int> > > chi_q;

    FUNC_LIB::function<std::complex<scalartype>, harmonics_dmn_type  > P_q_cluster;
    FUNC_LIB::function<std::complex<scalartype>, harmonics_dmn_type  > P_q_lattice;

    FUNC_LIB::function<std::complex<scalartype>,       lambda_dmn_type>                              leading_eigenvalues;
    FUNC_LIB::function<std::complex<scalartype>,       lambda_dmn_type>                              leading_phi_t_chi_0_phi;
    FUNC_LIB::function<std::complex<scalartype>,       lambda_dmn_type>                              leading_phi_t_Gamma_phi;

    FUNC_LIB::function<std::complex<scalartype>, dmn_2<lambda_dmn_type, harmonics_dmn_type> >        leading_symmetries;
    FUNC_LIB::function<std::complex<scalartype>, dmn_2<lambda_dmn_type, lattice_eigenvector_dmn_t> > leading_eigenvectors;

    FUNC_LIB::function<std::complex<scalartype>, k_dmn_cut_type> S_k_cut;
    FUNC_LIB::function<std::complex<scalartype>, k_dmn_cut_type> a_k_cut;

    FUNC_LIB::function<std::complex<scalartype>, dmn_2<lambda_dmn_type, cluster_eigenvector_dmn_t> > leading_U_K;
    FUNC_LIB::function<std::complex<scalartype>, dmn_2<lambda_dmn_type, cluster_eigenvector_dmn_t> > leading_Vt_K;

    FUNC_LIB::function<std::complex<scalartype>, dmn_2<lambda_dmn_type, lattice_eigenvector_dmn_t> > leading_U_k;
    FUNC_LIB::function<std::complex<scalartype>, dmn_2<lambda_dmn_type, lattice_eigenvector_dmn_t> > leading_Vt_k;
  };

  template<class parameters_type, class MOMS_type>
  BSE_solver<parameters_type, MOMS_type>::BSE_solver(parameters_type& parameters_ref, MOMS_type& MOMS_ref):
    parameters(parameters_ref),
    concurrency(parameters.get_concurrency()),

    MOMS(MOMS_ref),

    BSE_cluster_solver_obj(parameters, MOMS),
    BSE_lattice_solver_obj(parameters, MOMS),

    cluster_eigenvector_dmn(),
    lattice_eigenvector_dmn(),

    //eigensystem_pln(lattice_eigenvector_dmn.get_size()),
    diagrammatic_symmetries_obj(parameters),

    harmonics("harmonics"),

    G4("G4"),
    G4_0("G4_0"),

    Gamma_cluster("Gamma_cluster"),
    Gamma_lattice("Gamma_lattice"),

    chi_0_function("chi_0_function"),

    chi("chi"),
    chi_0("chi_0"),

    chi_q("chi_q"),

    P_q_cluster("P_q_cluster"),
    P_q_lattice("P_q_lattice"),

    leading_eigenvalues("leading_eigenvalues"),
    leading_phi_t_chi_0_phi("leading_phi_t_chi_0_phi"),
    leading_phi_t_Gamma_phi("leading_phi_t_Gamma_phi"),

    leading_symmetries("leading_symmetries"),
    leading_eigenvectors("leading_eigenvectors"),

    S_k_cut("S_k_cut"),
    a_k_cut("a_k_cut"),

    leading_U_K("leading_U_K"),
    leading_Vt_K("leading_Vt_K"),

    leading_U_k("leading_U_k_interpolated"),
    leading_Vt_k("leading_Vt_k_interpolated")
  {
    initialize_wave_functions();

    {
      profiler_t prof("compute-H(k)", "input", __LINE__);

      wannier_interpolation<k_LDA, k_DCA >::execute(MOMS.H_LDA, MOMS.H_DCA);
      wannier_interpolation<k_LDA, k_HOST>::execute(MOMS.H_LDA, MOMS.H_HOST);
    }

    {
      profiler_t prof("compute-band-structure", "input", __LINE__);
      compute_band_structure::execute(parameters, MOMS.H_LDA, MOMS.band_structure);
    }

  }

  template<class parameters_type, class MOMS_type>
  BSE_solver<parameters_type, MOMS_type>::~BSE_solver()
  {}

  template<class parameters_type, class MOMS_type>
  void BSE_solver<parameters_type, MOMS_type>::write()
  {
    IO::FORMAT  FORMAT    = parameters.get_output_format();
    std::string file_name = parameters.get_directory() + parameters.get_susceptibilities_file_name();

    std::cout << "\n\n\t\t start writing " << file_name << "\n\n";

    switch(FORMAT)
      {
      case IO::JSON :
        {
          IO::writer<IO::JSON> writer;
          {
            writer.open_file(file_name);

            parameters.write(writer);
	    this->     write(writer);

            writer.close_file();
          }
        }
        break;

      case IO::HDF5 :
        {
          IO::writer<IO::HDF5> writer;
          {
            writer.open_file(file_name);

            parameters.write(writer);
            //MOMS      .write(writer);
            this->     write(writer);

            writer.close_file();
          }
        }
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }
  }

  template<class parameters_type, class MOMS_type>
  template<IO::FORMAT DATA_FORMAT>
  void BSE_solver<parameters_type, MOMS_type>::write(IO::writer<DATA_FORMAT>& writer)
  {
    writer.open_group("analysis-functions");

    {
      BSE_lattice_solver_obj.write(writer);
      BSE_cluster_solver_obj.write(writer);
    }

    writer.close_group();
  }

  
  template<class parameters_type, class MOMS_type>
  void BSE_solver<parameters_type, MOMS_type>::initialize_wave_functions()
  {
    wave_functions_names = "no-name";

    wave_functions_names(0) = "s-wave";
    wave_functions_names(1) = "p-wave";
    wave_functions_names(2) = "d-wave";
    
    {// s-wave
      std::complex<scalartype> norm_psi=0;

      for(int k_ind=0; k_ind<k_HOST_VERTEX::dmn_size(); k_ind++)
        harmonics(k_ind, 0) = 1.;

      for(int k_ind=0; k_ind<k_HOST_VERTEX::dmn_size(); k_ind++)
        norm_psi += harmonics(k_ind, 0)*conj(harmonics(k_ind, 0));

      for(int k_ind=0; k_ind<k_HOST_VERTEX::dmn_size(); k_ind++)
        harmonics(k_ind, 0) /= std::sqrt(norm_psi);
    }

    {// p-wave
      std::complex<scalartype> norm_psi=0;

      scalartype alpha_x = 1;//host_vertex_cluster_type::get_r_basis()[0][0];
      scalartype alpha_y = 1;//host_vertex_cluster_type::get_r_basis()[1][1];

      for(int k_ind=0; k_ind<k_HOST_VERTEX::dmn_size(); k_ind++)
        harmonics(k_ind, 1) = (cos(alpha_x*k_HOST_VERTEX::get_elements()[k_ind][0])
                               + cos(alpha_y*k_HOST_VERTEX::get_elements()[k_ind][1]));

      for(int k_ind=0; k_ind<k_HOST_VERTEX::dmn_size(); k_ind++)
        norm_psi += harmonics(k_ind, 1)*conj(harmonics(k_ind, 1));

      for(int k_ind=0; k_ind<k_HOST_VERTEX::dmn_size(); k_ind++)
        harmonics(k_ind, 1) /= std::sqrt(norm_psi);
    }

    {// d-wave
      std::complex<scalartype> norm_psi=0;

      scalartype alpha_x = 1;//host_vertex_cluster_type::get_r_basis()[0][0];
      scalartype alpha_y = 1;//host_vertex_cluster_type::get_r_basis()[1][1];

      for(int k_ind=0; k_ind<k_HOST_VERTEX::dmn_size(); k_ind++)
        harmonics(k_ind, 2) = (cos(alpha_x*k_HOST_VERTEX::get_elements()[k_ind][0])
                               - cos(alpha_y*k_HOST_VERTEX::get_elements()[k_ind][1]));

      for(int k_ind=0; k_ind<k_HOST_VERTEX::dmn_size(); k_ind++)
        norm_psi += harmonics(k_ind, 2)*conj(harmonics(k_ind, 2));

      for(int k_ind=0; k_ind<k_HOST_VERTEX::dmn_size(); k_ind++)
        harmonics(k_ind, 2) /= std::sqrt(norm_psi);
    }
  }


  template<class parameters_type, class MOMS_type>
  void BSE_solver<parameters_type, MOMS_type>::calculate_susceptibilities_2()
  {
    if(concurrency.id()==concurrency.last())
      std::cout << "\t" << __FUNCTION__ << std::endl;

    if(true)
      {
        BSE_cluster_solver_obj.compute_Gamma_cluster();

        Gamma_cluster = BSE_cluster_solver_obj.get_Gamma_matrix();
      }
    else
      {
        apply_symmetries();

        load_the_matrices();

        compute_Gamma_cluster();
      }

    {
      BSE_lattice_solver_obj.compute_Gamma_lattice_3(Gamma_cluster);

      BSE_lattice_solver_obj.compute_chi_0_lattice(chi_0);

      Gamma_lattice = BSE_lattice_solver_obj.get_Gamma_lattice();

      BSE_lattice_solver_obj.diagonolize_Gamma_chi_0(Gamma_lattice, chi_0);
    }
  }


  template<class parameters_type, class MOMS_type>
  void BSE_solver<parameters_type, MOMS_type>::apply_symmetries()
  {
    profiler_t prof(__FUNCTION__, __FILE__, __LINE__);

    if(concurrency.id()==concurrency.last())
      std::cout << "\t" << __FUNCTION__ << std::endl << std::endl;

    symmetrize::execute(MOMS.Sigma, MOMS.H_symmetry);

    symmetrize::execute(MOMS.G_k_w, MOMS.H_symmetry);
  }

}
#endif
