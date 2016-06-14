//-*-C++-*-

#ifndef FERMIONIC_HAMILTONIAN_H
#define FERMIONIC_HAMILTONIAN_H

namespace DCA
{
  namespace EXACT_DIAGONALIZATION
  {
    template<typename scalar_type>
    struct V_struct
    {
    public:

      long long index;

      scalar_type value;
    };

    template<typename scalar_type>
    struct t_struct
    {
    public:

      long long lhs;
      long long rhs;

      scalar_type value;
      //std::complex<scalar_type> value;
    };

    template<typename scalar_type>
    struct U_struct
    {
    public:

      long long lhs;
      long long rhs;

      scalar_type value;
      //std::complex<scalar_type> value;
    };

    /*!
     *   \author Peter Staar
     */
    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    class fermionic_Hamiltonian
    {
#include "type_definitions.h"

    public:

      typedef ED_type_definitions<parameter_type, b_dmn, s_dmn, r_dmn> ED_type_def;

      typedef typename ED_type_def::profiler_t       profiler_t;
      typedef typename ED_type_def::concurrency_type concurrency_type;

      typedef typename ED_type_def::scalar_type  scalar_type;
      typedef typename ED_type_def::complex_type complex_type;

      typedef typename ED_type_def::vector_type         vector_type;
      typedef typename ED_type_def::matrix_type         matrix_type;
      typedef typename ED_type_def::int_matrix_type int_matrix_type;

      typedef typename ED_type_def::occ_dmn occ_dmn;
      typedef typename ED_type_def::mag_dmn mag_dmn;

      typedef typename ED_type_def::occ_mag_dmn occ_mag_dmn;

      typedef fermionic_Fock_space<parameter_type, b_dmn, s_dmn, r_dmn> fermionic_Fock_space_type;

    public:

      fermionic_Hamiltonian(parameter_type&            parameters_ref,
                            fermionic_Fock_space_type& Fock_space_ref);

      ~fermionic_Hamiltonian();

      void initialize(FUNC_LIB::function<std::complex<double>, dmn_3<dmn_2<b_dmn, s_dmn>, dmn_2<b_dmn, s_dmn>, r_dmn> >& H_0,
                      FUNC_LIB::function<double              , dmn_3<dmn_2<b_dmn, s_dmn>, dmn_2<b_dmn, s_dmn>, r_dmn> >& H_i);

      void construct_Hamiltonians(bool interacting);

      void diagonolize_Hamiltonians_st();
      void diagonolize_Hamiltonians_mt();

      void set_spectrum(FUNC_LIB::function<double, w_REAL>& A_w);
      void print_spectrum();

      FUNC_LIB::function<vector_type, dmn_2<occ_dmn, mag_dmn> >& get_eigen_energies() { return eigen_energies; }

      FUNC_LIB::function<matrix_type, dmn_2<occ_dmn, mag_dmn> >& get_eigen_states  () { return eigen_states; }

    private:

      void initialize_t_ij_and_U_ij(FUNC_LIB::function<std::complex<double>, dmn_3<dmn_2<b_dmn, s_dmn>, dmn_2<b_dmn, s_dmn>, r_dmn> >& H_0,
                                    FUNC_LIB::function<double              , dmn_3<dmn_2<b_dmn, s_dmn>, dmn_2<b_dmn, s_dmn>, r_dmn> >& H_i);

      void shift_the_energies();

      void add_V_to_Hamiltonian(int N, matrix_type& H, int_matrix_type& Psi);
      void add_T_to_Hamiltonian(int N, matrix_type& H, int_matrix_type& Psi);
      void add_U_to_Hamiltonian(int N, matrix_type& H, int_matrix_type& Psi);

      void print_processor_lay_out();

    private:

      parameter_type&   parameters;
      concurrency_type& concurrency;

      fermionic_Fock_space_type& Fock_space;

      FUNC_LIB::function<int            , dmn_2<occ_dmn, mag_dmn> >& n_occupation_states;
      FUNC_LIB::function<int_matrix_type, dmn_2<occ_dmn, mag_dmn> >&   occupation_states;

      FUNC_LIB::function<std::vector<overlap_indices>, dmn_2<occ_mag_dmn, occ_mag_dmn> >& creation_overlap_of_states;
      FUNC_LIB::function<std::vector<overlap_indices>, dmn_2<occ_mag_dmn, occ_mag_dmn> >& annihilation_overlap_of_states;

      FUNC_LIB::function<std::vector<int>, dmn_2<occ_mag_dmn, occ_mag_dmn> >& creation_overlap_break_points;
      FUNC_LIB::function<std::vector<int>, dmn_2<occ_mag_dmn, occ_mag_dmn> >& annihilation_overlap_break_points;

      double CUT_OFF;

      std::vector<V_struct<scalar_type> >  V_i;
      std::vector<t_struct<complex_type> > t_ij;
      std::vector<U_struct<complex_type> > U_ij;

      FUNC_LIB::function<matrix_type, dmn_2<occ_dmn, mag_dmn> > Hamiltonians;

      FUNC_LIB::function<vector_type, dmn_2<occ_dmn, mag_dmn> > eigen_energies;
      FUNC_LIB::function<matrix_type, dmn_2<occ_dmn, mag_dmn> > eigen_states;
    };

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    fermionic_Hamiltonian<parameter_type, b_dmn, s_dmn, r_dmn>::fermionic_Hamiltonian(parameter_type&            parameters_ref,
                                                                                      fermionic_Fock_space_type& Fock_space_ref):
      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      Fock_space(Fock_space_ref),

      n_occupation_states(Fock_space.n_occupation_states),
      occupation_states  (Fock_space.  occupation_states),

      creation_overlap_of_states    (Fock_space.creation_overlap_of_states),
      annihilation_overlap_of_states(Fock_space.annihilation_overlap_of_states),

      creation_overlap_break_points    (Fock_space.creation_overlap_break_points),
      annihilation_overlap_break_points(Fock_space.annihilation_overlap_break_points),

      CUT_OFF(parameters.get_eigenvalue_cut_off()),

      V_i (0),
      t_ij(0),
      U_ij(0),

      Hamiltonians  ("Hamiltonians"),

      eigen_energies("eigen_energies"),
      eigen_states  ("eigen_states")
    {}

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    fermionic_Hamiltonian<parameter_type, b_dmn, s_dmn, r_dmn>::~fermionic_Hamiltonian()
    {}

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Hamiltonian<parameter_type, b_dmn, s_dmn, r_dmn>::initialize(FUNC_LIB::function<std::complex<double>, dmn_3<dmn_2<b_dmn, s_dmn>, dmn_2<b_dmn, s_dmn>, r_dmn> >& H_0,
                                                                                FUNC_LIB::function<double              , dmn_3<dmn_2<b_dmn, s_dmn>, dmn_2<b_dmn, s_dmn>, r_dmn> >& H_i)
    {
      initialize_t_ij_and_U_ij(H_0, H_i);
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Hamiltonian<parameter_type, b_dmn, s_dmn, r_dmn>::initialize_t_ij_and_U_ij(FUNC_LIB::function<std::complex<double>, dmn_3<dmn_2<b_dmn, s_dmn>, dmn_2<b_dmn, s_dmn>, r_dmn> >& H_0,
                                                                                              FUNC_LIB::function<double              , dmn_3<dmn_2<b_dmn, s_dmn>, dmn_2<b_dmn, s_dmn>, r_dmn> >& H_i)
    {
      {
        for(int r_j=0; r_j<r_dmn::dmn_size(); r_j++){
          for(int s_j=0; s_j<s_dmn::dmn_size(); s_j++){
            for(int b_j=0; b_j<b_dmn::dmn_size(); b_j++){

              V_struct<scalar_type> V_obj;

              V_obj.index = b_j+b_dmn::dmn_size()*(s_j+s_dmn::dmn_size()*r_j);

              V_obj.value = -parameters.get_chemical_potential();

              V_i.push_back(V_obj);
            }
          }
        }
      }

      for(int r_j=0; r_j<r_dmn::dmn_size(); r_j++){
        for(int r_i=0; r_i<r_dmn::dmn_size(); r_i++){

          int delta_r = r_dmn::parameter_type::subtract(r_j, r_i);

          for(int s_j=0; s_j<s_dmn::dmn_size(); s_j++){
            for(int b_j=0; b_j<b_dmn::dmn_size(); b_j++){

              for(int s_i=0; s_i<s_dmn::dmn_size(); s_i++){
                for(int b_i=0; b_i<b_dmn::dmn_size(); b_i++){

                  if(abs(H_0(b_i,s_i, b_j,s_j, delta_r))>1.e-3)
                    {
                      t_struct<complex_type> t_obj;

                      t_obj.lhs = b_i+b_dmn::dmn_size()*(s_i+s_dmn::dmn_size()*r_i);
                      t_obj.rhs = b_j+b_dmn::dmn_size()*(s_j+s_dmn::dmn_size()*r_j);

                      t_obj.value = real(H_0(b_i,s_i, b_j,s_j, delta_r)/2.);

                      t_ij.push_back(t_obj);
                    }

                  if(abs(H_i(b_i,s_i, b_j,s_j, delta_r))>1.e-3)
                    {
                      U_struct<complex_type> U_obj;

                      U_obj.lhs = b_i+b_dmn::dmn_size()*(s_i+s_dmn::dmn_size()*r_i);
                      U_obj.rhs = b_j+b_dmn::dmn_size()*(s_j+s_dmn::dmn_size()*r_j);

                      U_obj.value = H_i(b_i,s_i, b_j,s_j, delta_r)/2.;

                      U_ij.push_back(U_obj);
                    }
                }
              }
            }
          }
        }
      }
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Hamiltonian<parameter_type, b_dmn, s_dmn, r_dmn>::construct_Hamiltonians(bool interacting)
    {
      if(concurrency.id()==0)
        cout << "\n\t" << __FUNCTION__ << endl;

      {// resize the matrices
        for(int i=0; i<occ_dmn::dmn_size(); i++){
          for(int j=0; j<mag_dmn::dmn_size(); j++){

            int N = n_occupation_states(i,j);

            if(N==0)
              Hamiltonians(i,j).resize_no_copy(0);
            else
              Hamiltonians(i,j).resize_no_copy(N);
          }
        }
      }

      {
        for(int i=0; i<occ_dmn::dmn_size(); i++){
          for(int j=0; j<mag_dmn::dmn_size(); j++){

            int N = n_occupation_states(i,j);

            if(N>0){

              matrix_type&     H   = Hamiltonians     (i,j);
              int_matrix_type& Psi = occupation_states(i,j);

              for(int l1=0; l1<N; l1++)
                for(int l0=0; l0<N; l0++)
                  H(l0,l1) = 0;

              add_V_to_Hamiltonian(N, H, Psi);
              add_T_to_Hamiltonian(N, H, Psi);

              if(interacting)
                add_U_to_Hamiltonian(N, H, Psi);
            }
          }
        }
      }
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Hamiltonian<parameter_type, b_dmn, s_dmn, r_dmn>::add_V_to_Hamiltonian(int N, matrix_type& H, int_matrix_type& Psi)
    {
      int Ns = b_dmn::dmn_size()*s_dmn::dmn_size()*r_dmn::dmn_size();

      int* Psi_lhs = new int[Ns];
      int* Psi_rhs = new int[Ns];

      for(int l0=0; l0<N; l0++)
        {
          for(int l=0; l<Ns; l++)
            Psi_lhs[l] = Psi(l,l0);

          for(int l1=0; l1<N; l1++)
            {
              for(int l=0; l<Ns; l++)
                Psi_rhs[l] = Psi(l,l1);

              int diff=0;
              for(int d=0; d<Ns; d++)
                diff += square(Psi_lhs[d]-Psi_rhs[d]);

              if(diff==0)
                {
                  for(size_t l=0; l<V_i.size(); l++)
                    H(l0,l1) += V_i[l].value*Psi_rhs[V_i[l].index];
                }
            }
        }

      delete [] Psi_lhs;
      delete [] Psi_rhs;
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Hamiltonian<parameter_type, b_dmn, s_dmn, r_dmn>::add_T_to_Hamiltonian(int N, matrix_type& H, int_matrix_type& Psi)
    {
      int Ns = b_dmn::dmn_size()*s_dmn::dmn_size()*r_dmn::dmn_size();

      int* Psi_lhs = new int[Ns];
      int* Psi_rhs = new int[Ns];
      int* Psi_new = new int[Ns];

      for(int l0=0; l0<N; l0++)
        {
          for(int l=0; l<Ns; l++)
            Psi_lhs[l] = Psi(l,l0);

          for(int l1=0; l1<N; l1++)
            {
              for(int l=0; l<Ns; l++)
                Psi_rhs[l] = Psi(l,l1);

              int diff=0;
              for(int d=0; d<Ns; d++)
                diff += square(Psi_lhs[d]-Psi_rhs[d]);

              if(diff==2)
                {
                  for(size_t l=0; l<t_ij.size(); l++)
                    {
                      if(Psi_lhs[t_ij[l].lhs] == 1 and Psi_lhs[t_ij[l].rhs] == 0 and
                         Psi_rhs[t_ij[l].lhs] == 0 and Psi_rhs[t_ij[l].rhs] == 1)
                        {
                          scalar_type phase = 1.;

                          for(int d=0; d<Ns; d++)
                            Psi_new[d] = Psi_rhs[d];

                          {// apply annihilation operator
                            for(int d=0; d<t_ij[l].rhs; d++)
                              if(Psi_new[d] == 1)
                                phase *= scalar_type(-1.);

                            Psi_new[t_ij[l].rhs] = 0;
                          }

                          {// apply creation operator
                            for(int d=0; d<t_ij[l].lhs; d++)
                              if(Psi_new[d] == 1)
                                phase *= scalar_type(-1.);

                            Psi_new[t_ij[l].lhs] = 1;
                          }

                          if(false)
                            {
                              int res=0;
                              for(int d=0; d<Ns; d++)
                                res += square(Psi_lhs[d]-Psi_new[d]);

                              if(res != 0)
                                throw std::logic_error(__FUNCTION__);
                            }

                          H(l0,l1) += phase*t_ij[l].value;
                        }
                    }
                }
            }
        }

      delete [] Psi_lhs;
      delete [] Psi_rhs;
      delete [] Psi_new;
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Hamiltonian<parameter_type, b_dmn, s_dmn, r_dmn>::add_U_to_Hamiltonian(int N, matrix_type& H, int_matrix_type& Psi)
    {
      int Ns = b_dmn::dmn_size()*s_dmn::dmn_size()*r_dmn::dmn_size();

      int* Psi_lhs = new int[Ns];
      int* Psi_rhs = new int[Ns];

      for(int l0=0; l0<N; l0++)
        {
          for(int l=0; l<Ns; l++)
            Psi_lhs[l] = Psi(l,l0);

          for(int l1=0; l1<N; l1++)
            {
              for(int l=0; l<Ns; l++)
                Psi_rhs[l] = Psi(l,l1);

              int diff=0;
              for(int d=0; d<Ns; d++)
                diff += square(Psi_lhs[d]-Psi_rhs[d]);

              if(diff==0)
                {
                  for(size_t l=0; l<U_ij.size(); l++)
                    {
                      H(l0,l1) += U_ij[l].value/scalar_type(4.);

                      if(abs(Psi_rhs[U_ij[l].lhs]-scalar_type(1.))<1.e-3)
                        H(l0,l1) -= U_ij[l].value/scalar_type(2.);

                      if(abs(Psi_rhs[U_ij[l].rhs]-scalar_type(1.))<1.e-3)
                        H(l0,l1) -= U_ij[l].value/scalar_type(2.);

                      if(abs(Psi_rhs[U_ij[l].lhs]-scalar_type(1.))<1.e-3 and
                         abs(Psi_rhs[U_ij[l].rhs]-scalar_type(1.))<1.e-3)
                        H(l0,l1) += U_ij[l].value;
                    }
                }

            }
        }

      delete [] Psi_lhs;
      delete [] Psi_rhs;
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Hamiltonian<parameter_type, b_dmn, s_dmn, r_dmn>::diagonolize_Hamiltonians_st()
    {
      if(concurrency.id()==0)
        cout << "\n\t" << __FUNCTION__ << "\n\n";

      int start = clock();

      for(int i=0; i<occ_dmn::dmn_size(); i++){
        for(int j=0; j<mag_dmn::dmn_size(); j++){

          int N = n_occupation_states(i,j);

          if(N==0)
            {
              eigen_energies(i,j).resize        (0);
              eigen_states  (i,j).resize_no_copy(0);
            }
          else
            {
              eigen_energies(i,j).resize        (N);
              eigen_states  (i,j).resize_no_copy(N);

              {
                if(concurrency.id()==0)
                  cout << "\t N : " << i << ", Sz : " << mag_dmn::parameter_type::Sz(j)
                       << ", \t size : " << N << ", \t time : ";

                int start = clock();
                LIN_ALG::GEEV<LIN_ALG::CPU>::execute('V', 'U', Hamiltonians(i,j), eigen_energies(i,j), eigen_states(i,j));
                int end = clock();

                if(concurrency.id()==0)
                  cout << double(end-start)/double(CLOCKS_PER_SEC) << "\n";
              }
            }
        }
      }

      int end = clock();

      if(concurrency.id()==0)
        {
          cout << "\n\t" << __FUNCTION__ << "\t total time : " << double(end-start)/double(CLOCKS_PER_SEC) << "\n\n";

          print_spectrum();
        }

      shift_the_energies();
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Hamiltonian<parameter_type, b_dmn, s_dmn, r_dmn>::diagonolize_Hamiltonians_mt()
    {
      if(concurrency.id()==0)
        cout << "\n\t" << __FUNCTION__ << "\n\n";

      int start = clock();

      for(int i=0; i<occ_dmn::dmn_size(); i++){
        for(int j=0; j<mag_dmn::dmn_size(); j++){

          int N = n_occupation_states(i,j);

          if(N==0)
            {
              eigen_energies(i,j).resize        (0);
              eigen_states  (i,j).resize_no_copy(0);
            }
          else
            {
              eigen_energies(i,j).resize        (N);
              eigen_states  (i,j).resize_no_copy(N);

              for(int l=0; l<N; l++)
                eigen_energies(i,j)[l] = 0;

              for(int l_1=0; l_1<N; l_1++)
                for(int l_0=0; l_0<N; l_0++)
                  eigen_states(i,j)(l_0,l_1) = 0;
            }
        }
      }

      print_processor_lay_out();

      int index  = 0;
      int N_proc = concurrency.number_of_processors();

      for(int i=0; i<occ_dmn::dmn_size(); i++){
        for(int j=0; j<mag_dmn::dmn_size(); j++){

          int N = n_occupation_states(i,j);

          if(N!=0)
            {
              if((index % N_proc) == concurrency.id())
                LIN_ALG::GEEV<LIN_ALG::CPU>::execute('V', 'U', Hamiltonians(i,j), eigen_energies(i,j), eigen_states(i,j));

              index += 1;
            }
        }
      }

      for(int i=0; i<occ_dmn::dmn_size(); i++){
        for(int j=0; j<mag_dmn::dmn_size(); j++){

          int N = n_occupation_states(i,j);

          if(N!=0)
            {
              concurrency.sum(eigen_energies(i,j));
              concurrency.sum(eigen_states  (i,j));
            }
        }
      }

      int end = clock();

      if(concurrency.id()==0)
        {
          cout << "\n\t" << __FUNCTION__ << "\t total time : " << double(end-start)/double(CLOCKS_PER_SEC) << "\n\n";

          print_spectrum();
        }

      shift_the_energies();
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Hamiltonian<parameter_type, b_dmn, s_dmn, r_dmn>::print_processor_lay_out()
    {
      if(concurrency.id()==0)
        {
          int index  = 0;
          int N_proc = concurrency.number_of_processors();

          cout << "\n\n\tN \\ Sz\t";
          for(int j=0; j<mag_dmn::dmn_size(); j++)
            cout << mag_dmn::parameter_type::Sz(j) << "\t";

          cout << "\n";
          for(int i=0; i<occ_dmn::dmn_size(); i++){

            cout << "\tN = "<< i << "\t";

            for(int j=0; j<mag_dmn::dmn_size(); j++)
              {
                //cout << n_occupation_states(i,j) << "\t";

                int N = n_occupation_states(i,j);
                if(N!=0)
                  {
                    cout << (index % N_proc) << "\t";
                    index += 1;
                  }
                else
                  cout << " "               << "\t";

              }
            cout << "\n";
          }
          cout << "\n";
        }
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Hamiltonian<parameter_type, b_dmn, s_dmn, r_dmn>::print_spectrum()
    {
      double N_BINS = 20;

      double E_MIN = 0.;
      double E_MAX = 0.;

      for(int i=0; i<occ_dmn::dmn_size(); i++)
        for(int j=0; j<mag_dmn::dmn_size(); j++)
          for(int n=0; n<n_occupation_states(i,j); n++)
            E_MIN = E_MIN > eigen_energies(i,j)[n]? eigen_energies(i,j)[n] : E_MIN;

      for(int i=0; i<occ_dmn::dmn_size(); i++)
        for(int j=0; j<mag_dmn::dmn_size(); j++)
          for(int n=0; n<n_occupation_states(i,j); n++)
            E_MAX = E_MAX < eigen_energies(i,j)[n]? eigen_energies(i,j)[n] : E_MAX;

      if(concurrency.id()==0)
        cout << "\n\t E_min : " << E_MIN << "\t E_max : " << E_MAX << "\n\n";

      if(E_MAX-E_MIN>1.e-2)
	{
	  size_t NUMBER_OF_LAMBDAS = 0;

	  for(int i=0; i<occ_dmn::dmn_size(); i++)
	    for(int j=0; j<mag_dmn::dmn_size(); j++)
	      NUMBER_OF_LAMBDAS += eigen_energies(i,j).size();
	  
	  double delta = (E_MAX-E_MIN)/(N_BINS-1.);
	  std::vector<size_t> y(N_BINS, 0);

	  for(int i=0; i<occ_dmn::dmn_size(); i++)
	    for(int j=0; j<mag_dmn::dmn_size(); j++)
	      for(int n=0; n<n_occupation_states(i,j); n++)
		y[int((eigen_energies(i,j)[n]-E_MIN)/delta)] += 1;

	  if(concurrency.id()==0)
	    {
	      cout << "\n\t distribution of the energies : \n\n";
	      for(int l=0; l<N_BINS; l++)
		cout << "\t" << E_MIN+delta/2. + l*delta << "\t" << std::string(int(double(y[l])/double(NUMBER_OF_LAMBDAS)*400.),'*') << endl;
	      cout << "\n\n";
	    }
	}
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Hamiltonian<parameter_type, b_dmn, s_dmn, r_dmn>::set_spectrum(FUNC_LIB::function<double, w_REAL>& A_w)
    {
      A_w = 0;

      std::vector<double>& w_elem = w_REAL::get_elements();

      for(int i=0; i<occ_dmn::dmn_size(); i++)
        for(int j=0; j<mag_dmn::dmn_size(); j++)
          for(int n=0; n<n_occupation_states(i,j); n++)
            for(int w_ind=0; w_ind<w_REAL::dmn_size()-1; w_ind++)
              if(w_elem[w_ind] <= eigen_energies(i,j)[n] and eigen_energies(i,j)[n] < w_elem[w_ind+1])
                A_w(w_ind) += 1.;

      double total = 0;
      for(int w_ind=0; w_ind<w_REAL::dmn_size()-1; w_ind++)
        total += A_w(w_ind);

      A_w/= total;
    }

    template<typename parameter_type, typename b_dmn, typename s_dmn, typename r_dmn>
    void fermionic_Hamiltonian<parameter_type, b_dmn, s_dmn, r_dmn>::shift_the_energies()
    {
      double E_0 = 0.;
      for(int i=0; i<occ_dmn::dmn_size(); i++)
        for(int j=0; j<mag_dmn::dmn_size(); j++)
          for(int n=0; n<n_occupation_states(i,j); n++)
            E_0 = E_0 > eigen_energies(i,j)[n]? eigen_energies(i,j)[n] : E_0;

      for(int i=0; i<occ_dmn::dmn_size(); i++)
        for(int j=0; j<mag_dmn::dmn_size(); j++)
          for(int n=0; n<n_occupation_states(i,j); n++)
            eigen_energies(i,j)[n] -= E_0;

      {
        double beta = parameters.get_beta();

        int number_of_eigenvalues = 0;
        for(int i=0; i<occ_dmn::dmn_size(); i++)
          for(int j=0; j<mag_dmn::dmn_size(); j++)
            for(int n=0; n<n_occupation_states(i,j); n++)
              if(std::exp(-beta*eigen_energies(i,j)[n]) > CUT_OFF)
                number_of_eigenvalues += 1;

        int total_of_eigenvalues = 0;
        for(int i=0; i<occ_dmn::dmn_size(); i++)
          for(int j=0; j<mag_dmn::dmn_size(); j++)
            for(int n=0; n<n_occupation_states(i,j); n++)
              total_of_eigenvalues += 1;

        if(concurrency.id()==0)
          cout << "\n\n\t number of eigenvalues exp(-beta*lambda) > CUT_OFF: " << number_of_eigenvalues << ", " << total_of_eigenvalues << "\n\n";
      }
    }

  }

}

#endif
