//-*-C++-*-

#ifndef MATERIAL_LATTICE_TEMPLATE_H
#define MATERIAL_LATTICE_TEMPLATE_H

/*!
 *  \author Peter Staar
 */
template<material_name_type name, typename point_group_type>
class material_lattice
{
public:

  const static int DIMENSION = -1;
  const static int BANDS     = -1;

  typedef no_symmetry<DIMENSION> LDA_point_group;
  typedef point_group_type       DCA_point_group;

  static double* initialize_r_DCA_basis();
  static double* initialize_r_LDA_basis();

  static std::vector<int>                  get_flavors();
  static std::vector<std::vector<double> > get_a_vectors();

  static std::vector<std::pair<std::pair<int,int>, std::pair<int,int> > > get_orbital_permutations();

  template<class domain, class parameters_type>
  static void initialize_H_interaction(FUNC_LIB::function<double , domain >& H_interaction,
                                       parameters_type&            parameters);

  template<class domain>
  static void initialize_H_symmetry(FUNC_LIB::function<int, domain>& H_symmetry);

  template<class parameters_type>
  static std::complex<double> get_LDA_Hamiltonians(parameters_type& parameters,
                                                   std::vector<double> k, int b1, int s1, int b2, int s2);
};

#endif
