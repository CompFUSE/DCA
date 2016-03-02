//-*-C++-*-

typedef cluster_domain<double, 3, VASP_LATTICE, MOMENTUM_SPACE, PARALLELLEPIPEDUM> vasp_k_cluster_type;
typedef cluster_domain<double, 3, VASP_LATTICE, REAL_SPACE    , PARALLELLEPIPEDUM> vasp_r_cluster_type;

typedef dmn_0<vasp_k_cluster_type> k_vasp;
typedef dmn_0<vasp_r_cluster_type> r_vasp;

typedef dmn_0<DFT::VASP::dmft_band_domain   > b_dmft;
typedef dmn_0<DFT::VASP::dmft_orbital_domain> o_dmft;

typedef dmn_0<DFT::VASP::vasp_band_domain>    b_vasp;
typedef dmn_0<DFT::VASP::vasp_orbital_domain> o_vasp;

typedef dmn_0<DFT::VASP::vasp_brillouin_zone_cut_domain> bz_cut;
