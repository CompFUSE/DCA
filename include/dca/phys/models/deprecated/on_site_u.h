// INTERNAL: Keep this file for later reference.

//-*-C++-*-

#ifndef ON_SITE_U_H
#define ON_SITE_U_H

/*!
 *   \ingroup MODELS
 *
 *   \author Peter Staar
 *   \brief  Simple template that implements on-site interaction, for different MC-solvers.
 */
class on_site_u {
public:
  template <class vertex_pair_type, class parameters_type, class concurrency_type,
            class H_interactions_function_type>
  static void set_vertex(vertex_pair_type& vertex, parameters_type& parameters,
                         concurrency_type& concurrency, H_interactions_function_type& H_interactions) {
    typedef typename H_interactions_function_type::this_domain_type interaction_domain_type;

    typedef typename interaction_domain_type::template domain_typelist<0> dmn_typelist_0;
    typedef typename interaction_domain_type::template domain_typelist<2> dmn_typelist_2;

    typedef typename TypeAt<1, dmn_typelist_0>::type spin_domain_type;

    typedef typename TypeAt<0, dmn_typelist_2>::type r_dmn_type;

    typedef typename parameters_type::nu nu;

    int BANDS = parameters.get_interacting_orbitals().size();
    int origin = r_dmn_type::origin_index();

    do {
      int band_ind_1 = concurrency.get_random_number() * BANDS;
      // int band_ind_2 = concurrency.get_random_number()*BANDS;  // old version: gives the wrong
      // probabilities

      // In order to not double-count the intra-band interaction we need different probabilities
      // for chosing the second band to be the same/different to the first band. Each band
      // different to the first band should have twice the probablility to be chosen than the same
      // band.
      int band_ind_2 = band_ind_1;

      // Hardcoded (efficient) version for 2 bands
      // assert(BANDS == 2);
      // double p = 2./3;
      // if (concurrency.get_random_number() < p) {
      //   band_ind_2 = 1 - band_ind_1;
      // }

      // General version
      double p = 1. - 1. / (2. * BANDS - 1.);  // probablity to choose any different band
      if (concurrency.get_random_number() < p) {
        while (band_ind_2 == band_ind_1) {
          band_ind_2 = concurrency.get_random_number() * BANDS;
        }
      }

      vertex.get_bands().first = parameters.get_interacting_orbitals()[band_ind_1];
      vertex.get_bands().second = parameters.get_interacting_orbitals()[band_ind_2];

      vertex.get_e_spins().first =
          spin_domain_type::get_elements()[int(concurrency.get_random_number() * 2.)];

      if (vertex.get_bands().first == vertex.get_bands().second) {
        if (vertex.get_e_spins().first == e_UP)
          vertex.get_e_spins().second = e_DN;
        else
          vertex.get_e_spins().second = e_UP;
      }
      else
        vertex.get_e_spins().second =
            spin_domain_type::get_elements()[int(concurrency.get_random_number() * 2.)];

      vertex.get_spin_orbitals().first = domains::convert<int, nu>::spin_orbital(
          vertex.get_bands().first, vertex.get_e_spins().first);
      vertex.get_spin_orbitals().second = domains::convert<int, nu>::spin_orbital(
          vertex.get_bands().second, vertex.get_e_spins().second);
    } while (std::fabs(H_interactions(vertex.get_spin_orbitals().first,
                                      vertex.get_spin_orbitals().second, origin)) < 1.e-3);

    int r_site = int(r_dmn_type::get_size() * concurrency.get_random_number());
    vertex.get_r_sites().first = r_site;
    vertex.get_r_sites().second = r_site;
  }
};

#endif
