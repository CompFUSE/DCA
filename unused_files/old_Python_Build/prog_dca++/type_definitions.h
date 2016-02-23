//-*-C++-*-

//typedef no_symmetry<2> DCA_point_group_type;
//typedef no_symmetry<3> DCA_point_group_type;

//typedef symmetry_package_spg_lib DCA_point_group_type;

//typedef S4               DCA_point_group_type;
//typedef S4_plus               DCA_point_group_type;
typedef D4               DCA_point_group_type;
//typedef C4               DCA_point_group_type;
//typedef D4_trial         DCA_point_group_type;

//typedef C3               DCA_point_group_type;
//typedef S3               DCA_point_group_type;
//typedef D3               DCA_point_group_type;

//typedef C6               DCA_point_group_type;
//typedef S6               DCA_point_group_type;
//typedef D6               DCA_point_group_type;

//typedef Oh               DCA_point_group_type;

//typedef spg_symmetry_package DCA_point_group_type;

typedef square_lattice<DCA_point_group_type> lattice_type;
//typedef triangular_lattice<DCA_point_group_type> lattice_type;
//typedef bilayer_lattice<DCA_point_group_type> lattice_type;
//typedef cubic_lattice<DCA_point_group_type> lattice_type;

//typedef material_lattice<NiO, DCA_point_group_type> lattice_type;
//typedef material_lattice<CuO2, DCA_point_group_type> lattice_type;

//typedef Bett_cluster_square_2D  <DCA_point_group_type> lattice_type;
//typedef cuprate_three_band_model<DCA_point_group_type> lattice_type;

//typedef Bett_cluster_square_2D_singlets<DCA_point_group_type> lattice_type;
//typedef Bett_cluster_square_3D<DCA_point_group_type> lattice_type;
//typedef Bett_cluster_triangular_2D<DCA_point_group_type> lattice_type;
//typedef Andersen_Hamiltonians<DCA_point_group_type> lattice_type;
//typedef cuprate_model_3_bands<DCA_point_group_type> lattice_type;
//typedef Bett_cluster_square_2D_layered<DCA_point_group_type> lattice_type;
//typedef Bett_cluster_square_2D_layered_non_interacting<DCA_point_group_type> lattice_type;
//typedef bipartite_square_model<DCA_point_group_type> lattice_type;
//typedef cuprate_2_band_model<DCA_point_group_type> lattice_type;
//typedef cuprate_3_band_model<DCA_point_group_type> lattice_type;
//typedef cuprate_3_band_model_kent<DCA_point_group_type> lattice_type;
//typedef non_interacting_cluster_square_2D<DCA_point_group_type> lattice_type;

typedef on_site_u   interaction_type;
typedef tight_binding_model<lattice_type, interaction_type> model;

//typedef dft_model<3, DCA_point_group_type>                    model;
//typedef Koshevnikov_model                                     model;

//typedef cuprate_single_band_model<DCA_point_group_type> model;
//typedef cuprate_three_band_model<DCA_point_group_type> model;

//typedef bilayer_model<DCA_point_group_type> model;

// typedef LDA_lattice_parameters<model>         LDA_lattice_parameters_type;//LDA_parameters_type;
// typedef DCA_lattice_parameters<model>         DCA_lattice_parameters_type;//DCA_parameters_type;
// typedef PCM_lattice_parameters<model>         PCM_lattice_parameters_type;//DCA_parameters_type;
// typedef host_lattice_parameters<model>        host_lattice_parameters_type;
// typedef host_vertex_lattice_parameters<model> host_vertex_lattice_parameters_type;

//typedef cluster_domain_family<double, lattice_type::DIMENSION, LATTICE_SP, MOMENTUM_SPACE, PARALLELLEPIPEDUM> LDA_cluster_family_type;
typedef cluster_domain_family<double, lattice_type::DIMENSION, CLUSTER, BRILLOUIN_ZONE>    DCA_cluster_family_type;
typedef cluster_domain_family<double, lattice_type::DIMENSION, LATTICE_SP, BRILLOUIN_ZONE> HOST_sp_cluster_family_type;
typedef cluster_domain_family<double, lattice_type::DIMENSION, LATTICE_TP, BRILLOUIN_ZONE> HOST_tp_cluster_family_type;

// domain representations
// const static cluster_representation_type r_DCA_cluster_representation = FULL;//IRREDUCIBLE; 
// const static cluster_representation_type k_DCA_cluster_representation = FULL;//IRREDUCIBLE; 
// const static cluster_representation_type r_LDA_cluster_representation = FULL;
// const static cluster_representation_type k_LDA_cluster_representation = FULL;

//typedef cluster<DCA_lattice_parameters_type>         DCA_cluster_type;
//typedef cluster<LDA_lattice_parameters_type>         LDA_cluster_type;
//typedef cluster<PCM_lattice_parameters_type>         PCM_cluster_type;
//typedef cluster<host_lattice_parameters_type>        host_cluster_type;
//typedef cluster<host_vertex_lattice_parameters_type> host_vertex_cluster_type;

//typedef k_cluster<k_LDA_cluster_representation, LDA_cluster_type>  LDA_k_cluster_type;
typedef cluster_domain<double, lattice_type::DIMENSION, LATTICE_SP, MOMENTUM_SPACE, PARALLELLEPIPEDUM> LDA_k_cluster_type;
//typedef k_cluster<k_DCA_cluster_representation, DCA_cluster_type>  DCA_k_cluster_type;
typedef cluster_domain<double, lattice_type::DIMENSION, CLUSTER, MOMENTUM_SPACE, BRILLOUIN_ZONE> DCA_k_cluster_type;
//typedef k_cluster<k_DCA_cluster_representation, PCM_cluster_type>  PCM_k_cluster_type;
typedef cluster_domain<double, lattice_type::DIMENSION, CLUSTER, MOMENTUM_SPACE, BRILLOUIN_ZONE> PCM_k_cluster_type;
//typedef k_cluster<k_DCA_cluster_representation, host_cluster_type> host_k_cluster_type;
typedef cluster_domain<double, lattice_type::DIMENSION, LATTICE_SP, MOMENTUM_SPACE, BRILLOUIN_ZONE> host_k_cluster_type;

//typedef k_cluster<k_DCA_cluster_representation, host_vertex_cluster_type> host_vertex_k_cluster_type;
typedef cluster_domain<double, lattice_type::DIMENSION, LATTICE_TP       , MOMENTUM_SPACE, BRILLOUIN_ZONE> host_vertex_k_cluster_type;
//typedef cluster_domain<double, lattice_type::DIMENSION, LATTICE_TP_INTERP, MOMENTUM_SPACE, BRILLOUIN_ZONE> host_vertex_k_cluster_type;

//typedef k_cluster<VERTEX_cluster_representation, DCA_cluster_type>        DCA_k_VERTEX_cluster_type;

//typedef r_cluster<r_LDA_cluster_representation, LDA_cluster_type>  LDA_r_cluster_type;
typedef cluster_domain<double, lattice_type::DIMENSION, LATTICE_SP, REAL_SPACE, PARALLELLEPIPEDUM> LDA_r_cluster_type;
//typedef r_cluster<r_DCA_cluster_representation, DCA_cluster_type>  DCA_r_cluster_type;
typedef cluster_domain<double, lattice_type::DIMENSION, CLUSTER, REAL_SPACE, BRILLOUIN_ZONE> DCA_r_cluster_type;
//typedef r_cluster<r_DCA_cluster_representation, PCM_cluster_type>  PCM_r_cluster_type;
typedef cluster_domain<double, lattice_type::DIMENSION, CLUSTER, REAL_SPACE, BRILLOUIN_ZONE> PCM_r_cluster_type;
//typedef r_cluster<r_DCA_cluster_representation, host_cluster_type> host_r_cluster_type;
typedef cluster_domain<double, lattice_type::DIMENSION, LATTICE_SP, REAL_SPACE, BRILLOUIN_ZONE> host_r_cluster_type;

//typedef r_cluster<r_DCA_cluster_representation, host_vertex_cluster_type> host_vertex_r_cluster_type;
typedef cluster_domain<double, lattice_type::DIMENSION, LATTICE_TP, REAL_SPACE, BRILLOUIN_ZONE> host_vertex_r_cluster_type;

//typedef r_cluster<VERTEX_cluster_representation, DCA_cluster_type> DCA_r_VERTEX_cluster_type;

typedef dmn_0<LDA_k_cluster_type>         k_LDA;
typedef dmn_0<DCA_k_cluster_type>         k_DCA;
typedef dmn_0<PCM_k_cluster_type>         k_PCM;
typedef dmn_0<host_k_cluster_type>        k_HOST;
typedef dmn_0<host_vertex_k_cluster_type> k_HOST_VERTEX;
//typedef dmn_0<DCA_k_VERTEX_cluster_type>  k_DCA_VERTEX;

typedef dmn_0<LDA_r_cluster_type>         r_LDA;
typedef dmn_0<DCA_r_cluster_type>         r_DCA;
typedef dmn_0<PCM_r_cluster_type>         r_PCM;
typedef dmn_0<host_r_cluster_type>        r_HOST;
typedef dmn_0<host_vertex_r_cluster_type> r_HOST_VERTEX;
//typedef dmn_0<DCA_r_VERTEX_cluster_type>  r_DCA_VERTEX;

typedef centered_cluster_domain<r_HOST_VERTEX::parameter_type> crystal_harmonics_expansion;
typedef dmn_0<crystal_harmonics_expansion>                     crystal_harmonics_expansion_dmn_t;

typedef electron_spin_domain        electron_spin_domain_type;
typedef electron_band_domain        electron_band_domain_type;

// typedef HS_spin_domain              HS_spin_domain_type;
// typedef HS_field_sign_domain        HS_field_sign_domain_type;

typedef dmn_0<electron_spin_domain_type>  s;
typedef dmn_0<electron_band_domain_type>  b;
// typedef dmn_0<HS_spin_domain_type>        HS_s;
// typedef dmn_0<HS_field_sign_domain_type>  HS_f;

typedef particle_number_domain<b, r_DCA, s> pn__b_r_DCA_e_spin_domain_type;
typedef particle_number_domain<b, k_DCA, s> pn__b_k_DCA_e_spin_domain_type;

typedef dmn_0<pn__b_r_DCA_e_spin_domain_type> pn__b_r_DCA_s;
typedef dmn_0<pn__b_r_DCA_e_spin_domain_type> pn__b_k_DCA_s;

typedef time_domain           time_domain_type;
typedef frequency_domain frequency_domain_type;

// typedef legendre_domain<time_domain_type, SINGLE_PARTICLE_QUANTITY> legendre_sp_expansion_type;
// typedef dmn_0<legendre_sp_expansion_type>                           legendre_sp_dmn_t;

typedef DCA::vertex_time_domain<DCA::SP_TIME_DOMAIN         > sp_vertex_time_domain_type;
typedef DCA::vertex_time_domain<DCA::TP_TIME_DOMAIN         > tp_vertex_time_domain_type;
typedef DCA::vertex_time_domain<DCA::SP_TIME_DOMAIN_POSITIVE> sp_vertex_time_domain_pos_type;
typedef DCA::vertex_time_domain<DCA::TP_TIME_DOMAIN_POSITIVE> tp_vertex_time_domain_pos_type;

typedef DCA::vertex_frequency_domain<DCA::COMPACT>           compact_vertex_frequency_domain_type;
typedef DCA::vertex_frequency_domain<DCA::EXTENDED>          extended_vertex_frequency_domain_type;
typedef DCA::vertex_frequency_domain<DCA::COMPACT_POSITIVE>  compact_positive_vertex_frequency_domain_type;
typedef DCA::vertex_frequency_domain<DCA::EXTENDED_POSITIVE> extended_positive_vertex_frequency_domain_type;
typedef DCA::vertex_frequency_domain<DCA::EXTENDED_BOSONIC>  bosonic_vertex_frequency_domain_type;

typedef dmn_0<time_domain>      t;
typedef dmn_0<frequency_domain> w;

typedef dmn_0<frequency_domain_real_axis> w_REAL;
typedef dmn_0<frequency_domain_imag_axis> w_IMAG;

typedef dmn_0<sp_vertex_time_domain_type>     sp_time_dmn_t;
typedef dmn_0<tp_vertex_time_domain_type>     tp_time_dmn_t;
typedef dmn_0<sp_vertex_time_domain_pos_type> sp_time_pos_dmn_t;
typedef dmn_0<tp_vertex_time_domain_pos_type> tp_time_pos_dmn_t;

typedef dmn_0<compact_vertex_frequency_domain_type>           w_VERTEX;
//typedef dmn_0<compact_sorted_vertex_frequency_domain_type>    w_VERTEX_SORTED;

typedef dmn_0<extended_vertex_frequency_domain_type>          w_VERTEX_EXTENDED;
typedef dmn_0<compact_positive_vertex_frequency_domain_type>  w_VERTEX_POS;
typedef dmn_0<extended_positive_vertex_frequency_domain_type> w_VERTEX_EXTENDED_POS;
//typedef dmn_0<extended_sorted_vertex_frequency_domain_type>   w_VERTEX_EXTENDED_SORTED;
typedef dmn_0<bosonic_vertex_frequency_domain_type>           w_VERTEX_BOSONIC;
// typedef dmn_0<fermionic_vertex_frequency_domain_type>         w_VERTEX_FERMIONIC;
// typedef dmn_0<vertex_frequency_domain_core_sorted_type>           w_VERTEX_CORE_SORTED;
// typedef dmn_0<vertex_frequency_domain_high_frequency_sorted_type> w_VERTEX_HF_SORTED;
// typedef dmn_2<w_VERTEX_EXTENDED_SORTED,
// 	      w_VERTEX_EXTENDED_SORTED>  w_VERTEX_EXTENDED_SORTED_w_VERTEX_EXTENDED_SORTED;


typedef dmn_2<b,s>                     nu;
typedef dmn_2<b,r_DCA>                 b_r_DCA;
typedef dmn_2<b,k_DCA>                 b_k_DCA;
typedef dmn_2<nu,nu>                   nu_nu;
typedef dmn_2<nu_nu,nu_nu>             nu_nu__nu_nu;
typedef dmn_2<nu,r_DCA>                nu_r_DCA;
typedef dmn_2<nu,k_DCA>                nu_k_DCA;
typedef dmn_2<nu,w_REAL>               nu_w_REAL;
typedef dmn_2<b,b>                     b_b;
typedef dmn_2<b_b,b_b>                 b_b__b_b;
typedef dmn_2<s,s>                     s_s;
typedef dmn_2<w,k_DCA>                 w_k_DCA;
typedef dmn_2<k_DCA,w>                 k_DCA_w;
typedef dmn_2<k_DCA,t>                 k_DCA_t;
typedef dmn_2<t,k_DCA>                 t_k_DCA;
typedef dmn_2<k_DCA,w_REAL>            k_DCA_w_REAL;
typedef dmn_2<w_REAL,k_DCA>            w_REAL_k_DCA;
typedef dmn_2<k_DCA,w_VERTEX>          k_DCA_w_VERTEX;
typedef dmn_2<k_DCA,w_VERTEX_EXTENDED> k_DCA_w_VERTEX_EXTENDED;
typedef dmn_2<nu,k_LDA>                nu_k_LDA;
typedef dmn_2<r_DCA,r_DCA>             r_DCA_r_DCA;
typedef dmn_2<r_DCA,k_DCA>             r_DCA_k_DCA;
typedef dmn_2<k_DCA,k_DCA>             k_DCA_k_DCA;
// typedef dmn_2<k_DCA,k_PCM>             k_DCA_k_PCM;
// typedef dmn_2<k_PCM,k_DCA>             k_PCM_k_DCA;
// typedef dmn_2<r_PCM,r_PCM>             r_PCM_r_PCM;
// typedef dmn_2<k_PCM,k_PCM>             k_PCM_k_PCM;
// typedef dmn_2<r_PCM,k_PCM>             r_PCM_k_PCM;
// typedef dmn_2<r_DCA,r_PCM>             r_DCA_r_PCM;
// typedef dmn_2<k_PCM,r_PCM>             k_PCM_r_PCM;
typedef dmn_2<k_DCA,r_DCA>             k_DCA_r_DCA;
// typedef dmn_2<r_DCA_VERTEX,r_DCA>      r_DCA_VERTEX_r_DCA;
// typedef dmn_2<k_DCA_VERTEX,k_DCA>      k_DCA_VERTEX_k_DCA;

typedef dmn_3<b,b,r_DCA>            b_b_r_DCA;
typedef dmn_3<b,b,k_DCA>            b_b_k_DCA;
// typedef dmn_3<b,b,r_PCM>            b_b_r_PCM;
// typedef dmn_3<b,b,k_PCM>            b_b_k_PCM;
typedef dmn_3<b,r_DCA,s>            b_r_DCA_s;
typedef dmn_3<b,k_DCA,s>            b_k_DCA_s;
typedef dmn_3<b,s,k_LDA>            b_s__k_LDA;
typedef dmn_3<b_b,s,k_LDA>          b_b__s__k_LDA;
typedef dmn_3<nu,nu,k_LDA>          nu_nu_k_LDA;
typedef dmn_3<nu,nu,r_LDA>          nu_nu_r_LDA;
typedef dmn_3<nu,nu,r_DCA>          nu_nu_r_DCA;
typedef dmn_3<nu,nu,k_DCA>          nu_nu_k_DCA;
// typedef dmn_3<nu,nu,r_PCM>          nu_nu_r_PCM;
// typedef dmn_3<nu,nu,k_PCM>          nu_nu_k_PCM;
typedef dmn_3<nu,nu,w>              nu_nu_w;
typedef dmn_3<nu,nu,w_REAL>         nu_nu_w_REAL;

typedef dmn_4<b ,b,r_DCA,r_DCA>     b_b_r_DCA_r_DCA;

typedef dmn_4<nu,nu,k_DCA,w>        nu_nu_k_DCA_w;
typedef dmn_4<nu,nu,k_DCA,t>        nu_nu_k_DCA_t;
typedef dmn_4<nu,nu,r_DCA,t>        nu_nu_r_DCA_t;
typedef dmn_4<nu,nu,r_DCA,w>        nu_nu_r_DCA_w;

typedef dmn_4<nu,nu,k_PCM,w>        nu_nu_k_PCM_w;
typedef dmn_4<nu,nu,k_PCM,t>        nu_nu_k_PCM_t;
typedef dmn_4<nu,nu,r_PCM,t>        nu_nu_r_PCM_t;
typedef dmn_4<nu,nu,r_PCM,w>        nu_nu_r_PCM_w;

typedef dmn_4<nu,nu,k_HOST,w>        nu_nu_k_HOST_w;
typedef dmn_4<nu,nu,k_HOST,t>        nu_nu_k_HOST_t;
typedef dmn_4<nu,nu,r_HOST,t>        nu_nu_r_HOST_t;
typedef dmn_4<nu,nu,r_HOST,w>        nu_nu_r_HOST_w;

typedef dmn_4<t,r_DCA,nu,nu>        t_r_DCA_nu_nu;


typedef dmn_4<b,b,b,b>              b_b_b_b;
typedef dmn_4<nu,nu,k_DCA,w_REAL>   nu_nu_k_DCA_w_REAL;
typedef dmn_4<nu,nu,k_DCA,w_IMAG>   nu_nu_k_DCA_w_IMAG;

// typedef dmn_5<nu,nu,HS_s,HS_f,r_DCA>       nu_nu_HS_s_HS_f_r_DCA;
// typedef dmn_5<nu,nu,HS_s,HS_f,r_PCM>       nu_nu_HS_s_HS_f_r_PCM;
typedef dmn_5<nu,nu, nu,nu, k_DCA>         nu_nu__nu_nu__k_DCA;
typedef dmn_5<b,b,r_DCA,w_VERTEX,w_VERTEX> b_b_r_DCA_w_VERTEX_w_VERTEX;

// typedef dmn_6<nu,nu,HS_s,HS_s,HS_f,r_DCA>                   nu_nu_HS_s_HS_s_HS_f_r_DCA;
//typedef dmn_6<nu,nu,HS_s,HS_s,HS_f,r_PCM>                   nu_nu_HS_s_HS_s_HS_f_r_PCM;
typedef dmn_6<nu,nu,nu,nu,k_DCA,w_VERTEX_BOSONIC>           nu_nu__nu_nu__k_DCA_w_VERTEX_BOSONIC;
typedef dmn_6<b ,b ,b ,b, k_DCA,w_VERTEX_BOSONIC>           b_b__b_b__k_DCA_w_VERTEX_BOSONIC;

// typedef dmn_6<b,b,r_DCA_VERTEX,r_DCA,w_VERTEX,w_VERTEX>     b_b_r_DCA_VERTEX_r_DCA_w_VERTEX_w_VERTEX;
// typedef dmn_6<b,b,k_DCA_VERTEX,k_DCA,w_VERTEX,w_VERTEX>     b_b_k_DCA_VERTEX_k_DCA_w_VERTEX_w_VERTEX;

// typedef dmn_6<b,b,r_PCM,r_PCM,w_VERTEX,w_VERTEX>            b_b_r_PCM_r_PCM_w_VERTEX_w_VERTEX;
// typedef dmn_6<b,b,k_PCM,k_PCM,w_VERTEX,w_VERTEX>            b_b_k_PCM_k_PCM_w_VERTEX_w_VERTEX;

// typedef dmn_6<b,b,r_PCM,r_PCM,w_VERTEX_POS,w_VERTEX>        b_b_r_PCM_r_PCM_w_VERTEX_POS_w_VERTEX;
// typedef dmn_6<b,b,k_PCM,k_PCM,w_VERTEX_POS,w_VERTEX>        b_b_k_PCM_k_PCM_w_VERTEX_POS_w_VERTEX;

// typedef dmn_6<b,b,r_DCA_VERTEX,r_DCA,w_VERTEX_POS,w_VERTEX> b_b_r_DCA_VERTEX_r_DCA_w_VERTEX_POS_w_VERTEX;
// typedef dmn_6<b,b,k_DCA_VERTEX,k_DCA,w_VERTEX_POS,w_VERTEX> b_b_k_DCA_VERTEX_k_DCA_w_VERTEX_POS_w_VERTEX;

// typedef dmn_6<b,b,w_VERTEX_EXTENDED_SORTED,
// 	      b,b,w_VERTEX_EXTENDED_SORTED>                 b_b_w_VERTEX_EXTENDED_SORTED__b_b_w_VERTEX_EXTENDED_SORTED;
typedef dmn_6<b,b,k_DCA,b,b,k_DCA>                          b_b_k_DCA_b_b_k_DCA;

//typedef dmn_8<b,b,b,b,k_DCA_VERTEX,k_DCA,w_VERTEX,w_VERTEX> b_b__b_b__k_DCA_VERTEX_k_DCA_w_VERTEX_w_VERTEX;
//typedef dmn_8<b,b,b,b,k_PCM,k_PCM,w_VERTEX,w_VERTEX>        b_b__b_b__k_PCM_k_PCM_w_VERTEX_w_VERTEX;


// band-structure plot

typedef brillouin_zone_cut_domain<101>        brillouin_zone_cut_domain_type;
typedef dmn_0<brillouin_zone_cut_domain_type> k_domain_cut_dmn_type;  
typedef dmn_2<nu, k_domain_cut_dmn_type>      nu_k_cut;

// average Feynman expansion order histogram 

// static const int MAX_MEASURED_DIAGRAM_ORDER = 10000;
// typedef Feynman_expansion_order_domain<MAX_MEASURED_DIAGRAM_ORDER> diagram_order_type;
// typedef dmn_0<diagram_order_type>                                  diagram_order_domain_type;

// functions over the DMFT-iteration steps

typedef dmn_0<DCA_iteration_domain> DCA_iteration_domain_type;

// information structures

//typedef QMC::QMC_information_structure QMC_information_structure_type;
//typedef dca::DCA_information_structure DCA_information_structure_type;

// random number generator

typedef random_number_generator random_number_generator; 

// analysis

typedef dmn_4<b,b,r_DCA,w_VERTEX>                     b_b_r_DCA_w_VERTEX;
typedef dmn_4<b,b,k_DCA,w_VERTEX>                     b_b_k_DCA_w_VERTEX;

// typedef dmn_4<b,b,r_PCM,w_VERTEX>                     b_b_r_PCM_w_VERTEX;
// typedef dmn_4<b,b,k_PCM,w_VERTEX>                     b_b_k_PCM_w_VERTEX;


//typedef dmn_4<b,b,k_DCA,w_VERTEX_SORTED>              b_b_k_DCA_w_VERTEX_SORTED;
typedef dmn_2<b_b_k_DCA_w_VERTEX, b_b_k_DCA_w_VERTEX> b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX;

// typedef dmn_2<b_b_k_PCM_w_VERTEX, b_b_k_PCM_w_VERTEX> b_b_k_PCM_w_VERTEX__b_b_k_PCM_w_VERTEX;

typedef dmn_4<b,b,k_DCA,w_VERTEX_EXTENDED>                              b_b_k_DCA_w_VERTEX_EXTENDED;
typedef dmn_2<b_b_k_DCA_w_VERTEX_EXTENDED, b_b_k_DCA_w_VERTEX_EXTENDED> b_b_k_DCA_w_VERTEX_EXTENDED__b_b_k_DCA_w_VERTEX_EXTENDED;

// typedef dmn_4<b,b,k_PCM,w_VERTEX_EXTENDED>                              b_b_k_PCM_w_VERTEX_EXTENDED;
// typedef dmn_4<b,b,r_PCM,w_VERTEX_EXTENDED>                              b_b_r_PCM_w_VERTEX_EXTENDED;
// typedef dmn_2<b_b_k_PCM_w_VERTEX_EXTENDED, b_b_k_PCM_w_VERTEX_EXTENDED> b_b_k_PCM_w_VERTEX_EXTENDED__b_b_k_PCM_w_VERTEX_EXTENDED;

// typedef dmn_4<b, b, k_DCA, w_VERTEX_CORE_SORTED>                              b_b_k_DCA_w_VERTEX_CORE_SORTED;
// typedef dmn_2<b_b_k_DCA_w_VERTEX_CORE_SORTED, b_b_k_DCA_w_VERTEX_CORE_SORTED> b_b_k_DCA_w_VERTEX_CORE_SORTED__b_b_k_DCA_w_VERTEX_CORE_SORTED;

// typedef dmn_4<b, b, k_DCA, w_VERTEX_EXTENDED_SORTED>                                  b_b_k_DCA_w_VERTEX_EXTENDED_SORTED;
// typedef dmn_2<b_b_k_DCA_w_VERTEX_EXTENDED_SORTED, b_b_k_DCA_w_VERTEX_EXTENDED_SORTED> b_b_k_DCA_w_VERTEX_EXTENDED_SORTED__b_b_k_DCA_w_VERTEX_EXTENDED_SORTED;


