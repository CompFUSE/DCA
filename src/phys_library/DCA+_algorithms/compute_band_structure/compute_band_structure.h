//-*-C++-*-

#ifndef COMPUTE_BAND_STRUCTURE_H
#define COMPUTE_BAND_STRUCTURE_H

/*!
 *  \author peter staar
 */
class compute_band_structure
{
#include "type_definitions.h"

  // we give it a prime number such that it is positioned on an edge or corner of a patch!
  static const int INTERPOLATION_POINTS_BAND_STRUCTURE = brillouin_zone_cut_domain_type::INTERPOLATION_POINTS;//101;

public:

  template<typename K_dmn_t, typename parameter_type>
  static void execute(parameter_type&                                        parameters,
                      function<std::complex<double>, dmn_3<nu,nu,K_dmn_t> >& H_LDA,
                      function<double, nu_k_cut>&                            bands);

private:

  static void construct_path(std::string                        coordinate_type,
                             std::vector<std::vector<double> >  path_vecs,
                             std::vector<std::vector<double> >& collection_k_vecs);

  static void high_symmetry_line_1D(std::vector<std::vector<double> >& collection_k_vecs);
  static void high_symmetry_line_2D(std::vector<std::vector<double> >& collection_k_vecs);
  static void high_symmetry_line_3D(std::vector<std::vector<double> >& collection_k_vecs);

  static void high_symmetry_line_1D(std::string& name, std::vector<std::vector<double> >& k_vecs);

  static void high_symmetry_line_2D(std::string& name, std::vector<std::vector<double> >& k_vecs);

  static void high_symmetry_line_3D(std::string& name, std::vector<std::vector<double> >& k_vecs);
};

template<typename K_dmn_t, typename parameter_type>
void compute_band_structure::execute(parameter_type&                                        parameters,
                                     function<std::complex<double>, dmn_3<nu,nu,K_dmn_t> >& H_LDA,
                                     function<double, nu_k_cut>&                            band_structure)
{
  //int matrix_dim = s::dmn_size()*b::dmn_size();

  std::vector<std::vector<double> > collection_k_vecs;

  { // construct the path in the Brilluoin zone ...
    if(parameters.get_Brillouin_zone_vectors().size() == 0)
      {
        if(model::DIMENSION == 1)
          high_symmetry_line_1D(parameters.get_coordinate_type(),
                                parameters.get_Brillouin_zone_vectors());

        if(model::DIMENSION == 2)
          high_symmetry_line_2D(parameters.get_coordinate_type(),
                                parameters.get_Brillouin_zone_vectors());

        if(model::DIMENSION == 3)
          high_symmetry_line_3D(parameters.get_coordinate_type(),
                                parameters.get_Brillouin_zone_vectors());
      }

    construct_path(parameters.get_coordinate_type(),
                   parameters.get_Brillouin_zone_vectors(),
                   collection_k_vecs);

    brillouin_zone_cut_domain_type::get_size()     = collection_k_vecs.size();
    brillouin_zone_cut_domain_type::get_elements() = collection_k_vecs;

    band_structure.reset();
  }

  function<std::complex<double>, dmn_3<nu, nu, k_domain_cut_dmn_type> > H_k("H_k_interpolated");

  { // get H(k)
    wannier_interpolation<K_dmn_t, k_domain_cut_dmn_type>::execute(H_LDA, H_k);
  }

  { // compute the bands ...

    LIN_ALG::vector<             double , LIN_ALG::CPU> L_vec(nu::dmn_size());
    LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> H_mat(nu::dmn_size());
    LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> V_mat(nu::dmn_size());

    for(int l=0; l<int(collection_k_vecs.size()); l++)
      {
        for(int i=0; i<nu::dmn_size(); i++)
          for(int j=0; j<nu::dmn_size(); j++)
            H_mat(i,j) = H_k(i,j,l);

        LIN_ALG::GEEV<LIN_ALG::CPU>::execute('N', 'U', H_mat, L_vec, V_mat);

        for(int i=0; i<b::dmn_size(); i++)
          for(int j=0; j<s::dmn_size(); j++)
            band_structure(i, j, l) = L_vec[2*i+j];
      }
  }
}

void compute_band_structure::construct_path(std::string                        coordinate_type,
                                            std::vector<std::vector<double> >  path_vecs,
                                            std::vector<std::vector<double> >& collection_k_vecs)
{
  //cout << __FUNCTION__ << endl;

  int Nb_interpolation = INTERPOLATION_POINTS_BAND_STRUCTURE;

  collection_k_vecs.resize(0);

  std::vector<double> k(model::DIMENSION);

  for(int i=0; i<int(path_vecs.size())-1; i++)
    {
      std::vector<double>& p0 = path_vecs[i];
      std::vector<double>& p1 = path_vecs[i+1];

      //       VECTOR_OPERATIONS::PRINT(p0); cout << endl;
      //       VECTOR_OPERATIONS::PRINT(p1); cout << endl;

      assert(model::DIMENSION == int(p0.size()));
      assert(model::DIMENSION == int(p1.size()));

      std::vector<double> p0_tmp(model::DIMENSION, 0.);
      std::vector<double> p1_tmp(model::DIMENSION, 0.);

      if(coordinate_type == "relative")
        {
          for(int di=0; di<model::DIMENSION; di++)
            for(int dj=0; dj<model::DIMENSION; dj++)
              p0_tmp[dj] += p0[di]*k_DCA::parameter_type::get_super_basis_vectors()[di][dj];

          for(int di=0; di<model::DIMENSION; di++)
            for(int dj=0; dj<model::DIMENSION; dj++)
              p1_tmp[dj] += p1[di]*k_DCA::parameter_type::get_super_basis_vectors()[di][dj];
        }
      else
        {
          p0_tmp = p0;
          p1_tmp = p1;
        }

      //       VECTOR_OPERATIONS::PRINT(p0_tmp); cout << endl;
      //       VECTOR_OPERATIONS::PRINT(p1_tmp); cout << endl << endl;

      for(int l=0; l<Nb_interpolation; l++)
        {

          for(int z=0; z<model::DIMENSION; z++)
            k[z] = (1.-double(l)/double(Nb_interpolation))*p0_tmp[z] + double(l)/double(Nb_interpolation)*p1_tmp[z];

          collection_k_vecs.push_back(k);
        }
    }
}

void compute_band_structure::high_symmetry_line_1D(std::string& name,
                                                   std::vector<std::vector<double> >& k_vecs)
{
  name = "absolute";

  std::vector<double>& b0 = k_DCA::parameter_type::get_super_basis_vectors()[0];

  std::vector<double> k0(1);
  std::vector<double> k1(2);

  for(int i=0; i<1; i++)
    {
      k0[i]  = 0    *b0[i];//+0.   *b1[i];
      k1[i]  = 1./2.*b0[i];//+1./2.*b1[i];
    }

  k_vecs.resize(0);

  k_vecs.push_back(k0);
  k_vecs.push_back(k1);
}

void compute_band_structure::high_symmetry_line_2D(std::string& name,
                                                   std::vector<std::vector<double> >& k_vecs)
{
  name = "absolute";

  std::vector<double>& b0 = k_DCA::parameter_type::get_super_basis_vectors()[0];
  std::vector<double>& b1 = k_DCA::parameter_type::get_super_basis_vectors()[1];

  std::vector<double> k0(2);
  std::vector<double> k1(2);
  std::vector<double> k2(2);
  std::vector<double> k3(2);

  for(int i=0; i<2; i++)
    {
      k0[i]  = 0    *b0[i]+0.   *b1[i];
      k1[i]  = 1./2.*b0[i]+1./2.*b1[i];
      k2[i]  = 1./2.*b0[i]+0.   *b1[i];
      k3[i]  = 0.   *b0[i]+1./2.*b1[i];
    }

  k_vecs.resize(0);

  k_vecs.push_back(k0);
  k_vecs.push_back(k1);
  k_vecs.push_back(k2);
  k_vecs.push_back(k3);
  k_vecs.push_back(k0);
}

void compute_band_structure::high_symmetry_line_3D(std::string& name,
                                                   std::vector<std::vector<double> >& k_vecs)
{
  name = "absolute";

  std::vector<double>& b0 = k_DCA::parameter_type::get_super_basis_vectors()[0];
  std::vector<double>& b1 = k_DCA::parameter_type::get_super_basis_vectors()[1];
  std::vector<double>& b2 = k_DCA::parameter_type::get_super_basis_vectors()[2];

  std::vector<double> k0(3);
  std::vector<double> k1(3);
  std::vector<double> k2(3);
  std::vector<double> k3(3);

  for(int i=0; i<3; i++)
    {
      k0[i]  = 0    *b0[i]+0.   *b1[i]+0.   *b2[i];
      k1[i]  = 1./2.*b0[i]+1./2.*b1[i]+1./2.*b2[i];
      k2[i]  = 1./2.*b0[i]+1./2.*b1[i]+0.   *b2[i];
      k3[i]  = 1./2.*b0[i]+0.   *b1[i]+0.   *b2[i];
    }

  k_vecs.resize(0);

  k_vecs.push_back(k0);
  k_vecs.push_back(k1);
  k_vecs.push_back(k2);
  k_vecs.push_back(k3);
  k_vecs.push_back(k0);
}

void compute_band_structure::high_symmetry_line_1D(std::vector<std::vector<double> >& collection_k_vecs)
{
  //cout << __FUNCTION__ << endl;

  int Nb_interpolation = INTERPOLATION_POINTS_BAND_STRUCTURE;

  std::vector<double>& b0 = k_DCA::parameter_type::get_super_basis_vectors()[0];

  for(int l=0; l<Nb_interpolation; l++)
    {
      std::vector<double> k(1,0);

      k[0] = double(l)/double(Nb_interpolation)*b0[0];

      collection_k_vecs.push_back(k);
    }
}

void compute_band_structure::high_symmetry_line_2D(std::vector<std::vector<double> >& collection_k_vecs)
{
  //cout << __FUNCTION__ << endl;

  int Nb_interpolation = INTERPOLATION_POINTS_BAND_STRUCTURE;

  std::vector<double>& b0 = k_DCA::parameter_type::get_super_basis_vectors()[0];
  std::vector<double>& b1 = k_DCA::parameter_type::get_super_basis_vectors()[1];

  std::vector<double> k0(2);
  std::vector<double> k1(2);
  std::vector<double> k2(2);
  std::vector<double> k3(2);
  std::vector<double> k4(2);

  for(int i=0; i<2; i++)
    {
      k0[i]  = 0    *b0[i]+0.   *b1[i];
      k1[i]  = 1./2.*b0[i]+1./2.*b1[i];
      k2[i]  = 1./2.*b0[i]+0.   *b1[i];
      k3[i]  = 0.   *b0[i]+1./2.*b1[i];
      k4[i]  = 0.   *b0[i]+0.   *b1[i];
    }

  for(int l=0; l<Nb_interpolation; l++)
    {
      std::vector<double> k(2,0);

      k[0] = (1.-double(l)/double(Nb_interpolation))*k0[0] + double(l)/double(Nb_interpolation)*k1[0];
      k[1] = (1.-double(l)/double(Nb_interpolation))*k0[1] + double(l)/double(Nb_interpolation)*k1[1];

      collection_k_vecs.push_back(k);
    }


  for(int l=0; l<Nb_interpolation; l++)
    {
      std::vector<double> k(2,0);

      k[0] = (1.-double(l)/double(Nb_interpolation))*k1[0] + double(l)/double(Nb_interpolation)*k2[0];
      k[1] = (1.-double(l)/double(Nb_interpolation))*k1[1] + double(l)/double(Nb_interpolation)*k2[1];

      collection_k_vecs.push_back(k);
    }

  for(int l=0; l<Nb_interpolation; l++)
    {
      std::vector<double> k(2,0);

      k[0] = (1.-double(l)/double(Nb_interpolation))*k2[0] + double(l)/double(Nb_interpolation)*k3[0];
      k[1] = (1.-double(l)/double(Nb_interpolation))*k2[1] + double(l)/double(Nb_interpolation)*k3[1];

      collection_k_vecs.push_back(k);
    }

  for(int l=0; l<Nb_interpolation; l++)
    {
      std::vector<double> k(2,0);

      k[0] = (1.-double(l)/double(Nb_interpolation))*k3[0] + double(l)/double(Nb_interpolation)*k4[0];
      k[1] = (1.-double(l)/double(Nb_interpolation))*k3[1] + double(l)/double(Nb_interpolation)*k4[1];

      collection_k_vecs.push_back(k);
    }
}

void compute_band_structure::high_symmetry_line_3D(std::vector<std::vector<double> >& collection_k_vecs)
{
  //cout << __FUNCTION__ << endl;

  int Nb_interpolation = INTERPOLATION_POINTS_BAND_STRUCTURE;

  //   std::vector<double>& b0 = DCA_cluster_type::get_k_super_basis_new()[0];
  //   std::vector<double>& b1 = DCA_cluster_type::get_k_super_basis_new()[1];
  //   std::vector<double>& b2 = DCA_cluster_type::get_k_super_basis_new()[2];

  std::vector<double>& b0 = k_DCA::parameter_type::get_super_basis_vectors()[0];
  std::vector<double>& b1 = k_DCA::parameter_type::get_super_basis_vectors()[1];
  std::vector<double>& b2 = k_DCA::parameter_type::get_super_basis_vectors()[2];

  std::vector<double> k0(3);
  std::vector<double> k1(3);
  std::vector<double> k2(3);
  std::vector<double> k3(3);

  std::vector<double> kA(3);
  std::vector<double> kB(3);

  for(int i=0; i<3; i++)
    {
      k0[i]  = 0    *b0[i]+0.   *b1[i]+0.   *b2[i];
      k1[i]  = 1./2.*b0[i]+1./2.*b1[i]+1./2.*b2[i];
      k2[i]  = 1./2.*b0[i]+1./2.*b1[i]+0.   *b2[i];
      k3[i]  = 1./2.*b0[i]+0.   *b1[i]+0.   *b2[i];

      kA[i]  = -3.*b0[i]  -3.   *b1[i]-3.   *b2[i];
      kB[i]  =  3.*b0[i]  +3.   *b1[i]+3.   *b2[i];
    }

  for(int l=0; l<Nb_interpolation; l++)
    {
      std::vector<double> k(3,0);

      k[0] = (1.-double(l)/double(Nb_interpolation))*k0[0] + double(l)/double(Nb_interpolation)*k1[0];
      k[1] = (1.-double(l)/double(Nb_interpolation))*k0[1] + double(l)/double(Nb_interpolation)*k1[1];
      k[2] = (1.-double(l)/double(Nb_interpolation))*k0[2] + double(l)/double(Nb_interpolation)*k1[2];

      collection_k_vecs.push_back(k);
    }

  for(int l=0; l<Nb_interpolation; l++)
    {
      std::vector<double> k(3,0);

      k[0] = (1.-double(l)/double(Nb_interpolation))*k1[0] + double(l)/double(Nb_interpolation)*k2[0];
      k[1] = (1.-double(l)/double(Nb_interpolation))*k1[1] + double(l)/double(Nb_interpolation)*k2[1];
      k[2] = (1.-double(l)/double(Nb_interpolation))*k1[2] + double(l)/double(Nb_interpolation)*k2[2];

      collection_k_vecs.push_back(k);
    }

  for(int l=0; l<Nb_interpolation; l++)
    {
      std::vector<double> k(3,0);

      k[0] = (1.-double(l)/double(Nb_interpolation))*k2[0] + double(l)/double(Nb_interpolation)*k3[0];
      k[1] = (1.-double(l)/double(Nb_interpolation))*k2[1] + double(l)/double(Nb_interpolation)*k3[1];
      k[2] = (1.-double(l)/double(Nb_interpolation))*k2[2] + double(l)/double(Nb_interpolation)*k3[2];

      collection_k_vecs.push_back(k);
    }

  for(int l=0; l<Nb_interpolation; l++)
    {
      std::vector<double> k(3,0);

      k[0] = (1.-double(l)/double(Nb_interpolation))*k3[0] + double(l)/double(Nb_interpolation)*k0[0];
      k[1] = (1.-double(l)/double(Nb_interpolation))*k3[1] + double(l)/double(Nb_interpolation)*k0[1];
      k[2] = (1.-double(l)/double(Nb_interpolation))*k3[2] + double(l)/double(Nb_interpolation)*k0[2];

      collection_k_vecs.push_back(k);
    }
}

#endif
