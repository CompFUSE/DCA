//-*-C++-*-

#ifndef ADVANCED_FERMIONIC_HILBERT_SPACE_PHI_REPRESENTATION_H
#define ADVANCED_FERMIONIC_HILBERT_SPACE_PHI_REPRESENTATION_H

namespace DCA
{
  namespace ADVANCED_EXACT_DIAGONALIZATION
  {
    /*
    template<typename parameter_type, typename ed_options>       // N: size of bitset sequence
    struct rep_struct
    {
    public:

      typedef typename ed_options::scalar_type  scalar_type;
      typedef typename ed_options::complex_type complex_type;
 
      typedef typename ed_options::phi_type phi_type;
 
    public:

      void sort();

    public:

      phi_type                  phi;

      std::vector<int>          index;
      std::vector<complex_type> alpha;
    };

    template<typename parameter_type, typename ed_options>
    void rep_struct<parameter_type, ed_options>::sort()
    {
      std::vector<int>          sorted_index;
      std::vector<complex_type> sorted_alpha;

      for(int i = 0; i < index.size(); ++i)
	{
	  int idx = 0;
	
	  while(idx < sorted_index.size() && index[i] >= sorted_index[idx])
	    {
	      ++idx;
	    }
	  
	  sorted_index.insert(sorted_index.begin()+idx, index[i]);
	  sorted_alpha.insert(sorted_alpha.begin()+idx, alpha[i]);
	}
      
      index.swap(sorted_index);
      alpha.swap(sorted_alpha);
    }


    template<typename parameter_type, typename ed_options>
    bool operator< (const rep_struct<parameter_type, ed_options>& phi_obj1, 
		    const rep_struct<parameter_type, ed_options>& phi_obj2)
    {
      return phi_obj1.phi.to_ulong() < phi_obj2.phi.to_ulong();
    }
    */

//     template<typename parameter_type, typename ed_options>       // N: size of bitset sequence
//     class Hilbert_space;

    template<typename parameter_type, typename ed_options>       // N: size of bitset sequence
    class Hilbert_space_rep
    {
      typedef typename ed_options::scalar_type  scalar_type;
      typedef typename ed_options::complex_type complex_type;
 
      typedef typename ed_options::phi_type phi_type;
 
    public:

      Hilbert_space_rep() {}

      template<typename Hilbert_space_type>
      void initialize(Hilbert_space_type& HS);

      int size() { return rep.size(); }

      phi_type&                  get_phi(int i)     { assert(i>=0 && i<size()); return rep[i].phi; }
      std::vector<int>&          get_indices(int i) { assert(i>=0 && i<size()); return rep[i].index; }
      std::vector<complex_type>& get_alphas(int i)  { assert(i>=0 && i<size()); return rep[i].alpha; }

      int find(const phi_type& phi_);

    private:

      void sort();

      void merge();

      void sort_phis();

      //std::vector< typename psi_state<parameter_type, ed_options>::phi_type > phis;
      //std::vector< std::vector< std::pair<int,complex_type> > > idx_coeff;
      std::vector< phi_state<parameter_type, ed_options, PHI_MULTIPLET> > rep;
    };

    template<typename parameter_type, typename ed_options>
    template<typename Hilbert_space_type>
    void Hilbert_space_rep<parameter_type, ed_options>::initialize(Hilbert_space_type& HS)//Hilbert_space<parameter_type, ed_options>& HS)
    {
      for(int i = 0; i < HS.size(); ++i)
	{	
	  psi_state<parameter_type, ed_options>& Psi = HS.get_element(i);

	  for (int j = 0; j < Psi.size(); ++j)
	    {	  
	      phi_state<parameter_type, ed_options, PHI_MULTIPLET> tmp;

	      tmp.phi   = Psi.get_phi(j);
	      tmp.index = std::vector<int>(1,i);
	      tmp.alpha = std::vector<complex_type>(1,Psi.get_alpha(j));

	      rep.push_back(tmp);
	    }
	}

      sort();
      
      merge();
      
      sort_phis();
    }

    template<typename parameter_type, typename ed_options>
    void Hilbert_space_rep<parameter_type, ed_options>::sort()
    {
      std::sort(rep.begin(), rep.end());
    }

    template<typename parameter_type, typename ed_options>
    void Hilbert_space_rep<parameter_type, ed_options>::merge()
    {
      std::vector< phi_state<parameter_type, ed_options, PHI_MULTIPLET> > merged_rep;

      merged_rep.push_back(rep[0]);

      for(int i = 1; i < size(); ++i){
	
	if(rep[i].phi == merged_rep.back().phi){
	  merged_rep.back().index.insert(merged_rep.back().index.end(), rep[i].index.begin(), rep[i].index.end());
	  merged_rep.back().alpha.insert(merged_rep.back().alpha.end(), rep[i].alpha.begin(), rep[i].alpha.end());
	}

	else{
	  merged_rep.push_back(rep[i]);
	}
      }

      rep.swap(merged_rep);
    }

    template<typename parameter_type, typename ed_options>
    void Hilbert_space_rep<parameter_type, ed_options>::sort_phis()
    {
      for(int i = 0; i < size(); ++i)
	rep[i].sort();
    }

    template<typename parameter_type, typename ed_options>
    int Hilbert_space_rep<parameter_type, ed_options>::find(const phi_type& phi_)
    {
      phi_state<parameter_type, ed_options, PHI_MULTIPLET> tmp;
      tmp.phi = phi_;

      typename std::vector< phi_state<parameter_type, ed_options, PHI_MULTIPLET> >::iterator it = std::lower_bound(rep.begin(), rep.end(), tmp);

      if(it != rep.end() && !(tmp < *it))
	return it - rep.begin();

      else
	return rep.size();
    }


    /*
    template<typename parameter_type, typename ed_options>       // N: size of bitset sequence
    class Hilbert_space
    {
      typedef typename ed_options::int_type     int_type;

      typedef typename ed_options::scalar_type  scalar_type;
      typedef typename ed_options::complex_type complex_type;
 
      typedef typename ed_options::phi_type phi_type;
 
    public:

      Hilbert_space(const std::vector< std::string >& names_, const std::vector< int >& eigenvalues_);
      
      Hilbert_space(const Hilbert_space<parameter_type, ed_options>& HS);


      void  insert(const psi_state<parameter_type, ed_options>& psi);
      
      void  remove(int i){ psi_states.erase(psi_states.begin()+i); }
      bool  remove(const psi_state<parameter_type, ed_options>& psi);

      bool  mark_state(const psi_state<parameter_type, ed_options>& psi);

      typename std::vector< psi_state<parameter_type, ed_options> >::iterator
      binary_find(const psi_state<parameter_type, ed_options>& psi);

      void  print(bool full=false);

      void  sort_wrt_magnetization();
      void  sort();

      int   size() const     { return psi_states.size(); }
      int   rep_size() const { return rep.size(); }

      bool  empty()          { return psi_states.empty(); }

      psi_state<parameter_type, ed_options>&   get_element(int state) { return psi_states[state]; }
      int                                        get_N_occ()            { return psi_states[0].occupation_number; }
      std::vector< std::string >                 get_name()             { return names; }
      std::vector< int >                         get_eigenvalues()      { return eigenvalues; }


      // assert or boolean to check wether occupation and magnetization is defined
      int                                        get_occupation()       { return eigenvalues[0]; }
      int                                        get_magnetization()    { return eigenvalues[1]; }

      Hilbert_space_rep<parameter_type, ed_options>& get_rep()            { return rep; }

      typename std::vector< psi_state<parameter_type, ed_options> >::iterator psi_states_begin()
      { return psi_states.begin(); }
      typename std::vector< psi_state<parameter_type, ed_options> >::iterator psi_states_end()
      { return psi_states.end(); }

      //Hilbert_space_rep<parameter_type, ed_options> get_compact_rep();

      void initialize_rep() { rep.initialize(*this); }

    private:

      std::vector< std::string >                             names;
      std::vector< int >                                     eigenvalues;
      std::vector< psi_state<parameter_type, ed_options> > psi_states;

      Hilbert_space_rep<parameter_type, ed_options>              rep;

    };


    template<typename parameter_type, typename ed_options>
    Hilbert_space<parameter_type, ed_options>::Hilbert_space(const std::vector< std::string >& names_,
							 const std::vector< int  >& eigenvalues_)
      : names(names_), eigenvalues(eigenvalues_)
    {}

    template<typename parameter_type, typename ed_options>
    Hilbert_space<parameter_type, ed_options>::Hilbert_space(const Hilbert_space<parameter_type, ed_options>& HS)
      : names(HS.names), eigenvalues(HS.eigenvalues), psi_states(HS.psi_states)
    {}

    template<typename parameter_type, typename ed_options>
    void Hilbert_space<parameter_type, ed_options>::insert(const psi_state<parameter_type, ed_options>& psi)
    {
      psi_states.push_back(psi);
    }

    template<typename parameter_type, typename ed_options>
    bool Hilbert_space<parameter_type, ed_options>::remove(const psi_state<parameter_type, ed_options>& psi)
    {
      typename std::vector< psi_state<parameter_type, ed_options> >::iterator
	it = std::find(psi_states.begin(), psi_states.end(), psi);

      if (it != psi_states.end()){
        psi_states.erase(it);
        return true;
      }

      return false;
    }

    template<typename parameter_type, typename ed_options>
    void Hilbert_space<parameter_type, ed_options>::print(bool full)
    {
      std::cout << size() << "\t" << rep.size();
      for (int_type i = 0; i < eigenvalues.size(); ++i){
        std::cout << "\t" << eigenvalues[i];
      }
      std::cout << "\n";
      if (full){
        for (int_type i = 0; i < psi_states.size(); ++i){
          psi_states[i].print();
          std::cout << "\n";
        }
      }
    }

    template<typename parameter_type, typename ed_options>
    void Hilbert_space<parameter_type, ed_options>::sort_wrt_magnetization()
    {
      std::sort(psi_states.begin(), psi_states.end(), less_magnetization<parameter_type, ed_options>);
    }


    template<typename parameter_type, typename ed_options>
    void Hilbert_space<parameter_type, ed_options>::sort()
    {
      std::sort(psi_states.begin(), psi_states.end());
    }

    template<typename parameter_type, typename ed_options>
    bool Hilbert_space<parameter_type, ed_options>::mark_state(const psi_state<parameter_type, ed_options>& psi)
    {
      // for (int i = 0; i < size(); ++i){
      //        if (psi_states[i] == psi){
      //          psi_states[i].marker = true;
      //          return true;
      //        }
      // }
      //        return false;

      typename std::vector< psi_state<parameter_type, ed_options> >::iterator it = binary_find(psi);

      if (it != psi_states.end()){
        it->mark();
        return true;
      }
      else
        return false;
    }

    // Make sure that the subspace is sorted
    template<typename parameter_type, typename ed_options>
    typename std::vector< psi_state<parameter_type, ed_options> >::iterator
    Hilbert_space<parameter_type, ed_options>::binary_find(const psi_state<parameter_type, ed_options>& psi)
    {
      typename std::vector< psi_state<parameter_type, ed_options> >::iterator
	it = std::lower_bound(psi_states.begin(), psi_states.end(), psi);

      while (it != psi_states.end() && !(operator< <parameter_type, ed_options>(psi,*it))){
        if (abs(abs(scalar_product<parameter_type, ed_options>(*it, psi))-scalar_type(1.)) < 1.e-10)
          return it;
        ++it;
      }

      return psi_states.end();
    }
    */

  }
}

#endif
