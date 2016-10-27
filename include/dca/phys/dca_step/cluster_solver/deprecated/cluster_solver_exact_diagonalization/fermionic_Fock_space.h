//-*-C++-*-

#ifndef FERMIONIC_FOCK_SPACE_H
#define FERMIONIC_FOCK_SPACE_H

namespace DCA
{
  namespace EXACT_DIAGONALIZATION
  {
    /*!
     *   \author Peter Staar
     */
    template<typename parameters_type, typename b_dmn, typename s_dmn, typename r_dmn>
    class fermionic_Fock_space
    {
    public:

      typedef size_t int_t;

      typedef dmn_3<b_dmn, s_dmn, r_dmn> b_s_r;

      typedef dmn_0<occupation_number_domain<b_dmn, s_dmn, r_dmn> > occ_dmn;
      typedef dmn_0<magnetization_domain    <b_dmn, s_dmn, r_dmn> > mag_dmn;

      typedef dmn_2<occ_dmn, mag_dmn> occ_mag_dmn;

      typedef LIN_ALG::matrix<int, LIN_ALG::CPU> matrix_type;

      typedef typename parameters_type::profiler_type    profiler_t;
      typedef typename parameters_type::concurrency_type concurrency_type;

    public:

      fermionic_Fock_space(parameters_type& parameters_ref);

      ~fermionic_Fock_space();

    private:

      void initialize();

      void compute_n_occupation_states();

      void print_occupation_matrix();

      void compute_occupation_states();

      void symmetrize_occupation_states();

      void compute_creation_overlap();

      void compute_annihilation_overlap();

      void print_creation_overlap();

      void print_annihilation_overlap();

      int_t number_of_states(int_t n, int_t l0, int_t l1);

      int_t factorial(int_t l0, int_t l1);

      int_t factorial(int_t n);

    public:

      parameters_type&  parameters;
      concurrency_type& concurrency;

      int max_occ;

      //     int* magnetization;
      //     int* fermionic_state;

      LIN_ALG::vector<int, LIN_ALG::CPU> magnetization;
      LIN_ALG::vector<int, LIN_ALG::CPU> fermionic_state;

      FUNC_LIB::function<int        , occ_mag_dmn > n_occupation_states;

      //FUNC_LIB::function<int*     , occ_mag_dmn > occupation_states;
      FUNC_LIB::function<matrix_type, occ_mag_dmn > occupation_states;

      FUNC_LIB::function<std::vector<overlap_indices>, dmn_2<occ_mag_dmn, occ_mag_dmn> > creation_overlap_of_states;
      FUNC_LIB::function<std::vector<overlap_indices>, dmn_2<occ_mag_dmn, occ_mag_dmn> > annihilation_overlap_of_states;

      FUNC_LIB::function<std::vector<int>, dmn_2<occ_mag_dmn, occ_mag_dmn> > creation_overlap_break_points;
      FUNC_LIB::function<std::vector<int>, dmn_2<occ_mag_dmn, occ_mag_dmn> > annihilation_overlap_break_points;
    };

    template<typename parameters_type, typename b_dmn, typename s_dmn, typename r_dmn>
    fermionic_Fock_space<parameters_type, b_dmn, s_dmn, r_dmn>::fermionic_Fock_space(parameters_type& parameters_ref):
      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      max_occ(b_dmn::dmn_size()*s_dmn::dmn_size()*r_dmn::dmn_size()),

      magnetization  ("magnetization"  , max_occ),
      fermionic_state("fermionic_state", max_occ),

      n_occupation_states("n_occupation_states"),
      occupation_states("occupation_states"),

      creation_overlap_of_states("creation_overlap_of_states"),
      annihilation_overlap_of_states("annihilation_overlap_of_states")
    {
      int index = 0;
      for(int r=0; r<r_dmn::dmn_size(); ++r){
        for(int i=0; i<s_dmn::dmn_size(); ++i){
          for(int j=0; j<b_dmn::dmn_size(); ++j){
            magnetization[index] = i==0? 1 : -1;

            index += 1;
          }
        }
      }

      initialize();
    }

    template<typename parameters_type, typename b_dmn, typename s_dmn, typename r_dmn>
    fermionic_Fock_space<parameters_type, b_dmn, s_dmn, r_dmn>::~fermionic_Fock_space()
    {}

    template<typename parameters_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Fock_space<parameters_type, b_dmn, s_dmn, r_dmn>::initialize()
    {
      compute_n_occupation_states();

      if(concurrency.id()==0)
        print_occupation_matrix();

      compute_occupation_states();

      //symmetrize_occupation_states();

      compute_creation_overlap();

      compute_annihilation_overlap();
    }

    template<typename parameters_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Fock_space<parameters_type, b_dmn, s_dmn, r_dmn>::compute_n_occupation_states()
    {
      n_occupation_states = 0;

      for(int N_occ=0; N_occ<=max_occ; ++N_occ)
        {
          for(int i=0; i<max_occ-N_occ; ++i)
            fermionic_state[i] = 0;

          for(int i=max_occ-N_occ; i<max_occ; ++i)
            fermionic_state[i] = 1;

          do
            {
              int Sz = b_dmn::dmn_size()*r_dmn::dmn_size();
              for(int i=0; i<max_occ; ++i)
                Sz += magnetization[i]*fermionic_state[i];

              n_occupation_states(N_occ, Sz) += 1;
            }
          while(std::next_permutation(&fermionic_state[0], &fermionic_state[0]+max_occ) );
        }
    }

    template<typename parameters_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Fock_space<parameters_type, b_dmn, s_dmn, r_dmn>::print_occupation_matrix()
    {
      cout << "\n\n\tN \\ Sz\t";
      for(int j=0; j<mag_dmn::dmn_size(); j++)
        cout << mag_dmn::parameter_type::Sz(j) << "\t";

      cout << "\n";
      for(int i=0; i<occ_dmn::dmn_size(); i++){

        cout << "\tN = "<< i << "\t";

        for(int j=0; j<mag_dmn::dmn_size(); j++)
          cout << n_occupation_states(i,j) << "\t";

        cout << "\n";
      }
      cout << "\n";
    }

    template<typename parameters_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Fock_space<parameters_type, b_dmn, s_dmn, r_dmn>::compute_occupation_states()
    {
      for(int i=0; i<occ_dmn::dmn_size(); i++){
        for(int j=0; j<mag_dmn::dmn_size(); j++){
          if(n_occupation_states(i,j)==0)
            {
              occupation_states(i,j).resize_no_copy(0);
            }
          else
            {
              occupation_states(i,j).resize_no_copy(std::pair<int, int>(max_occ, n_occupation_states(i,j)));
            }
        }
      }

      FUNC_LIB::function<int , dmn_2<occ_dmn, mag_dmn> > tmp;

      tmp = 0;

      for(int N_occ=0; N_occ<=max_occ; ++N_occ)
        {
          for(int i=0; i<max_occ-N_occ; ++i)
            fermionic_state[i] = 0;

          for(int i=max_occ-N_occ; i<max_occ; ++i)
            fermionic_state[i] = 1;

          do
            {
              int Sz = b_dmn::dmn_size()*r_dmn::dmn_size();
              for(int i=0; i<max_occ; ++i)
                Sz += magnetization[i]*fermionic_state[i];

              int j = tmp(N_occ, Sz);
              for(int i=0; i<max_occ; ++i)
                occupation_states(N_occ,Sz)(i,j) = fermionic_state[i];

              tmp(N_occ, Sz) += 1;
            }
          while(std::next_permutation(&fermionic_state[0], &fermionic_state[0]+max_occ) );
        }
    }

    template<typename parameters_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Fock_space<parameters_type, b_dmn, s_dmn, r_dmn>::symmetrize_occupation_states()
    {
      //     typedef typename r_dmn::parameter_type                r_cluster_type;
      //     typedef typename r_cluster_type::sym_super_cell_dmn_t sym_super_cell_dmn_t;

      //     static FUNC_LIB::function<std::pair<int,int>, dmn_2< dmn_2<r_dmn_t,b_dmn_t>, sym_super_cell_dmn_t > >& r_symmetry_matrix = r_cluster_type::get_symmetry_matrix();

      //     for(int S_ind=0; S_ind<sym_super_cell_dmn_t::dmn_size(); S_ind++)
      //       {

      //       }
    }

    template<typename parameters_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Fock_space<parameters_type, b_dmn, s_dmn, r_dmn>::compute_creation_overlap()
    {
      if(concurrency.id()==0)
        cout << "\n\n\t" << __FUNCTION__ << endl;

      overlap_indices overlap_obj;

      //     b_s_r   b_s_r_obj;
      //     int     coor[3];

      for(int n_0=0; n_0<occ_dmn::dmn_size(); n_0++){
        for(int Sz_0=0; Sz_0<mag_dmn::dmn_size(); Sz_0++){
          for(int n_1=0; n_1<occ_dmn::dmn_size(); n_1++){
            for(int Sz_1=0; Sz_1<mag_dmn::dmn_size(); Sz_1++){

              int N_0 = n_occupation_states(n_0, Sz_0);
              int N_1 = n_occupation_states(n_1, Sz_1);

              std::vector<overlap_indices>& overlap = creation_overlap_of_states(n_1, Sz_1, n_0, Sz_0);

              overlap.resize(0);

              if(N_0>0 and N_1>0 and n_1==n_0+1)
                {
                  matrix_type& ptr_0 = occupation_states(n_0, Sz_0);
                  matrix_type& ptr_1 = occupation_states(n_1, Sz_1);

                  for(int l_0=0; l_0<N_0; l_0++){
                    for(int l_1=0; l_1<N_1; l_1++){

                      int diff=0;
                      for(int d=0; d<max_occ; d++)
                        diff += square(ptr_0(d,l_0)-ptr_1(d,l_1));

                      if(diff==1)
                        {
                          int index=-1;
                          for(int d=0; d<max_occ; d++)
                            if(ptr_1(d,l_1) == 1 and ptr_0(d,l_0) == 0)
                              index=d;

                          assert(index>-1);
                          assert(ptr_0(index,l_0) == 0);
                          assert(ptr_1(index,l_1) == 1);

                          double sign = 1.;
                          for(int d=0; d<index; d++)
                            if(ptr_0(d,l_0) == 1)
                              sign *= -1.;

                          overlap_obj.index = index;

                          //                    b_s_r_obj.linind_2_subind(index, coor);
                          //                    overlap_obj.b_i = coor[0];
                          //                    overlap_obj.s_i = coor[1];
                          //                    overlap_obj.r_i = coor[2];

                          overlap_obj.lhs = l_1;
                          overlap_obj.rhs = l_0;

                          overlap_obj.sign = sign;

                          //creation_overlap_of_states(n_1, Sz_1, n_0, Sz_0).push_back(overlap_obj);

                          overlap.push_back(overlap_obj);
                        }
                    }
                  }
                }

              if(overlap.size() > 0)
                {
                  std::sort(overlap.begin(), overlap.end());

                  std::vector<int>& break_points = creation_overlap_break_points(n_1, Sz_1, n_0, Sz_0);

                  break_points.resize(0);
                  break_points.push_back(0);

                  for(size_t l=1; l<overlap.size(); l++)
                    if(overlap[l-1].index != overlap[l].index)
                      break_points.push_back(l);

                  break_points.push_back(overlap.size());

                  //            cout << "\n\n";
                  //            for(size_t l=0; l<overlap.size(); l++)
                  //              cout << "\t" << overlap[l].index << "\t" << overlap[l].lhs << "\t"<< overlap[l].rhs << "\n";
                  //            cout << "\n\n";
                }

              //          if(creation_overlap_of_states(n_0, Sz_0, n_1, Sz_1).size()>0)
              //            cout << n_0 << "\t" << mag_dmn::parameter_type::Sz(Sz_0) << "\t"
              //                 << n_1 << "\t" << mag_dmn::parameter_type::Sz(Sz_1) << "\t"
              //                 << creation_overlap_of_states(n_0, Sz_0, n_1, Sz_1).size()
              //                 << "\n";
            }
          }
        }
      }

      //     if(concurrency.id()==0)
      //       print_creation_overlap();
    }

    template<typename parameters_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Fock_space<parameters_type, b_dmn, s_dmn, r_dmn>::compute_annihilation_overlap()
    {
      if(concurrency.id()==0)
        cout << "\n\n\t" << __FUNCTION__ << endl;

      overlap_indices overlap_obj;

      //     b_s_r   b_s_r_obj;
      //     int     coor[3];

      for(int n_0=0; n_0<occ_dmn::dmn_size(); n_0++){
        for(int Sz_0=0; Sz_0<mag_dmn::dmn_size(); Sz_0++){
          for(int n_1=0; n_1<occ_dmn::dmn_size(); n_1++){
            for(int Sz_1=0; Sz_1<mag_dmn::dmn_size(); Sz_1++){

              int N_0 = n_occupation_states(n_0, Sz_0);
              int N_1 = n_occupation_states(n_1, Sz_1);

              std::vector<overlap_indices>& overlap = annihilation_overlap_of_states(n_1, Sz_1, n_0, Sz_0);

              overlap.resize(0);

              if(N_0>0 and N_1>0 and n_1==n_0-1)
                {
                  matrix_type& ptr_0 = occupation_states(n_0, Sz_0);
                  matrix_type& ptr_1 = occupation_states(n_1, Sz_1);

                  for(int l_0=0; l_0<N_0; l_0++){
                    for(int l_1=0; l_1<N_1; l_1++){

                      int diff=0;
                      for(int d=0; d<max_occ; d++)
                        diff += square(ptr_0(d,l_0)-ptr_1(d,l_1));

                      if(diff==1)
                        {
                          int index=-1;
                          for(int d=0; d<max_occ; d++)
                            if(ptr_1(d,l_1)==0 and ptr_0(d,l_0)==1)
                              index=d;

                          assert(index>-1);
                          assert(ptr_0(index,l_0)==1);
                          assert(ptr_1(index,l_1)==0);

                          double sign = 1.;
                          for(int d=0; d<index; d++)
                            if(ptr_1(d,l_1)==1)
                              sign *= -1.;

                          overlap_obj.index = index;

                          //                    b_s_r_obj.linind_2_subind(index, coor);
                          //                    overlap_obj.b_i = coor[0];
                          //                    overlap_obj.s_i = coor[1];
                          //                    overlap_obj.r_i = coor[2];

                          overlap_obj.lhs = l_1;
                          overlap_obj.rhs = l_0;

                          overlap_obj.sign = sign;

                          //annihilation_overlap_of_states(n_1, Sz_1, n_0, Sz_0).push_back(overlap_obj);

                          overlap.push_back(overlap_obj);
                        }
                    }
                  }
                }

              if(overlap.size() > 0)
                {
                  std::sort(overlap.begin(), overlap.end());

                  std::vector<int>& break_points = annihilation_overlap_break_points(n_1, Sz_1, n_0, Sz_0);

                  break_points.resize(0);
                  break_points.push_back(0);

                  for(size_t l=1; l<overlap.size(); l++)
                    if(overlap[l-1].index != overlap[l].index)
                      break_points.push_back(l);

                  break_points.push_back(overlap.size());
                }

              //          if(annihilation_overlap_of_states(n_0, Sz_0, n_1, Sz_1).size()>0)
              //            cout << n_0 << "\t" << mag_dmn::parameter_type::Sz(Sz_0) << "\t"
              //                 << n_1 << "\t" << mag_dmn::parameter_type::Sz(Sz_1) << "\t"
              //                 << annihilation_overlap_of_states(n_0, Sz_0, n_1, Sz_1).size()
              //                 << "\n";
            }
          }
        }
      }


    }

    template<typename parameters_type, typename b_dmn, typename s_dmn, typename r_dmn>
    size_t fermionic_Fock_space<parameters_type, b_dmn, s_dmn, r_dmn>::number_of_states(int_t l0, int_t l1, int_t l2)
    {
      return factorial(l0, l1)/factorial(l2);
    }

    template<typename parameters_type, typename b_dmn, typename s_dmn, typename r_dmn>
    size_t fermionic_Fock_space<parameters_type, b_dmn, s_dmn, r_dmn>::factorial(int_t n, int_t m)
    {
      switch(n)
        {
        case 0:
          return 1;

        case 1:
          return 1;

        default:
          return (n == m) ? 1 : factorial(n - 1, m) * n;
        }
    }

    template<typename parameters_type, typename b_dmn, typename s_dmn, typename r_dmn>
    size_t fermionic_Fock_space<parameters_type, b_dmn, s_dmn, r_dmn>::factorial(int_t n)
    {
      return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
    }


    template<typename parameters_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Fock_space<parameters_type, b_dmn, s_dmn, r_dmn>::print_creation_overlap()
    {
      cout << "\n\n\t creation_overlap : \n\n";

      for(int n_0=0; n_0<occ_dmn::dmn_size(); n_0++)
        for(int Sz_0=0; Sz_0<mag_dmn::dmn_size(); Sz_0++)
          for(int n_1=0; n_1<occ_dmn::dmn_size(); n_1++)
            for(int Sz_1=0; Sz_1<mag_dmn::dmn_size(); Sz_1++)
              if(creation_overlap_of_states(n_0, Sz_0, n_1, Sz_1).size()>0)
                {
                  cout << "\t" << n_0 << " , " << mag_dmn::parameter_type::Sz(Sz_0) << "\t -->  ";
                  cout << "\t" << n_1 << " , " << mag_dmn::parameter_type::Sz(Sz_1) << "\t :  ";
                  cout << creation_overlap_of_states(n_0, Sz_0, n_1, Sz_1).size()   << "\n";
                }

      cout << "\n\n";
    }

  }

}

#endif
