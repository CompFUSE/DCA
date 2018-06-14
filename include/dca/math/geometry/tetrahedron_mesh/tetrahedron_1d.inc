// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Full template specialization for a 1D tetrahedron.

template <>
struct tetrahedron<1> {
public:
  tetrahedron();
  ~tetrahedron();

  std::vector<double> compute_cm();

  void translate(const std::vector<double>& q);

  void plot(::dca::util::Plot& plot);
  void plot_q_vecs(::dca::util::Plot& plot);

  template <typename scalartype>
  void update_gaussian_domain(int& size, std::vector<scalartype>& weights,
                              std::vector<std::vector<scalartype>>& k_vecs);

  template <typename scalartype>
  void update_tetrahedron_domain(int& size, std::vector<scalartype>& weights,
                                 std::vector<std::vector<scalartype>>& k_vecs);

  template <typename mesh_t>
  void do_recursion(std::vector<tetrahedron<1>>& tetrahedra, mesh_t& mesh);

public:
  int index[2];

  double volume;

  std::vector<double> vec_0;
  std::vector<double> vec_1;

  std::vector<double> normal;

  int N_q;

  double* q_w;
  double* q_vecs;
};

template <typename mesh_t>
void tetrahedron<1>::do_recursion(std::vector<tetrahedron<1>>& tetrahedra, mesh_t& mesh) {
  std::vector<double> k0 = mesh[index[0]];
  std::vector<double> k1 = mesh[index[1]];

  static std::vector<double> km(1, 0);
  km[0] = (k0[0] + k1[0]) / 2.;

  int ind = mesh.size();
  mesh.push_back(km);

  static tetrahedron<1> tet_l;
  {
    tet_l.index[0] = index[0];
    tet_l.index[1] = ind;

    tet_l.vec_0 = k0;
    tet_l.vec_1 = km;

    tet_l.normal = this->normal;
    tet_l.volume = std::abs(km[0] - k0[0]);
  }

  static tetrahedron<1> tet_r;
  {
    tet_r.index[0] = ind;
    tet_r.index[1] = index[1];

    tet_r.vec_0 = km;
    tet_r.vec_1 = k1;

    tet_r.volume = std::abs(k1[0] - km[0]);
    tet_r.normal = normal;
  }

  tetrahedra.push_back(tet_l);
  tetrahedra.push_back(tet_r);
}
