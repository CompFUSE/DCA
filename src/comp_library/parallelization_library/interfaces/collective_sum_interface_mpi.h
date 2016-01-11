//-*-C++-*-

#ifndef COLLECTIVE_SUM_INTERFACE_MPI_H
#define COLLECTIVE_SUM_INTERFACE_MPI_H

namespace COMP_LIB
{
  /*!
   *  \author Peter Staar
   */
  template<>
  class collective_sum_interface<MPI_LIBRARY>
  {
  public:

    collective_sum_interface(processor_grouping<MPI_LIBRARY>& grouping_ref);
    ~collective_sum_interface();

    template<typename scalar_type>
    void sum(scalar_type& value);

    template<typename scalar_type>
    void sum(std::vector<scalar_type>& m);

    template<typename scalartype>
    void sum(std::map<std::string, std::vector<scalartype> >& m);

    template<typename scalar_type, class domain>
    void sum(FUNC_LIB::function<scalar_type, domain>& f);

    template<typename scalar_type, class domain>
    void sum(FUNC_LIB::function<scalar_type, domain>& f, FUNC_LIB::function<scalar_type, domain>& f_target);

    template<typename scalar_type, class domain>
    void sum(FUNC_LIB::function<std::vector<scalar_type>, domain>& f);

    template<typename scalar_type>
    void sum(LIN_ALG::vector<scalar_type, LIN_ALG::CPU>& f);

    template<typename scalar_type>
    void sum(LIN_ALG::matrix<scalar_type, LIN_ALG::CPU>& f);

    template<typename some_type>
    void sum_and_average(some_type& obj, int size);

    template<typename scalar_type, class domain>
    void average_and_compute_stddev(FUNC_LIB::function<scalar_type, domain>& f_mean, 
				    FUNC_LIB::function<scalar_type, domain>& f_stddev, 
				    size_t size);

    template<typename scalar_type, class domain>
    void average_and_compute_stddev(FUNC_LIB::function<std::complex<scalar_type>, domain>& f_mean, 
				    FUNC_LIB::function<std::complex<scalar_type>, domain>& f_stddev, 
				    size_t size);

  private:

    processor_grouping<MPI_LIBRARY>& grouping;
  };

  collective_sum_interface<MPI_LIBRARY>::collective_sum_interface(processor_grouping<MPI_LIBRARY>& grouping_ref):
    grouping(grouping_ref)
  {}

  collective_sum_interface<MPI_LIBRARY>::~collective_sum_interface()
  {}

  template<typename scalar_type>
  void collective_sum_interface<MPI_LIBRARY>::sum(scalar_type& value)
  {
    //cout << __PRETTY_FUNCTION__ << endl;

    scalar_type result;

    MPI_Allreduce(&value,
                  &result,
                  type_map_interface<MPI_LIBRARY, scalar_type>::factor(),
                  type_map_interface<MPI_LIBRARY, scalar_type>::value(),
                  MPI_SUM,
                  grouping.get());

    value = result;
  }

  template<typename scalar_type>
  void collective_sum_interface<MPI_LIBRARY>::sum(std::vector<scalar_type>& m)
  {
    //cout << __PRETTY_FUNCTION__ << endl;

    std::vector<scalar_type> result(m.size(),scalar_type(0));

    MPI_Allreduce(&(m[0]),
                  &(result[0]),
                  type_map_interface<MPI_LIBRARY, scalar_type>::factor()*m.size(),
                  type_map_interface<MPI_LIBRARY, scalar_type>::value(),
                  MPI_SUM,
                  grouping.get());

    for (size_t i=0;i<m.size();i++)
      m[i]=result[i];
  }

  template<typename scalar_type>
  void collective_sum_interface<MPI_LIBRARY>::sum(std::map<std::string, std::vector<scalar_type> >& m)
  {
    //cout << __PRETTY_FUNCTION__ << endl;

    typedef typename std::map <std::string, std::vector<scalar_type> >::iterator iterator_type;

    ////cout << id() << "\t" << m.size() << "\n";

    iterator_type it = m.begin();

    for( ; it != m.end(); ++it )
      {
        ////cout << it->first  << "\n";

        std::vector<scalar_type> values((it->second).size());

        for(size_t l=0; l<(it->second).size(); l++)
          values[l] = (it->second)[l];

        sum(values);

        for(size_t l=0; l<(it->second).size(); l++)
          (it->second)[l] = values[l];
      }
  }

  template<typename scalar_type, class domain>
  void collective_sum_interface<MPI_LIBRARY>::sum(FUNC_LIB::function<scalar_type, domain>& f)
  {
    //cout << __PRETTY_FUNCTION__ << "\t" << f.get_name() << endl;

    FUNC_LIB::function<scalar_type, domain> F;

    MPI_Allreduce(&f(0),
                  &F(0),
                  type_map_interface<MPI_LIBRARY, scalar_type>::factor()*f.size(),
                  type_map_interface<MPI_LIBRARY, scalar_type>::value(),
                  MPI_SUM,
                  grouping.get());

    for(int i=0;i<F.size();i++)
      f(i)=F(i);

    for(int i=0;i<F.size();i++)
      {
	if(f(i)!=f(i))
	  {
	    std::stringstream ss;
	    ss << i << "\t" << f.get_name() << "\n";
      std::cout << ss.str();

	    throw std::logic_error(__FUNCTION__);
	  }
      }
  }

  template<typename scalar_type, class domain>
  void collective_sum_interface<MPI_LIBRARY>::sum(FUNC_LIB::function<std::vector<scalar_type>, domain>& f)
  {
    //cout << __PRETTY_FUNCTION__ << endl;

    int Nr = f(0).size();
    int Nc = f   .size();

    LIN_ALG::matrix<scalar_type, LIN_ALG::CPU> M("M", std::pair<int,int>(Nr,Nc));

    for(int j=0; j<Nc; j++)
      for(int i=0; i<Nr; i++)
        M(i,j) = f(j)[i];

    sum(M);

    for(int j=0; j<Nc; j++)
      for(int i=0; i<Nr; i++)
        f(j)[i] = M(i,j);
  }

  template<typename scalar_type, class domain>
  void collective_sum_interface<MPI_LIBRARY>::sum(FUNC_LIB::function<scalar_type, domain>& f,
                                                  FUNC_LIB::function<scalar_type, domain>& f_target)
  {
    //cout << __PRETTY_FUNCTION__ << endl;

    MPI_Allreduce(&f(0),
                  &f_target(0),
                  type_map_interface<MPI_LIBRARY, scalar_type>::factor()*f.size(),
                  type_map_interface<MPI_LIBRARY, scalar_type>::value(),
                  MPI_SUM,
                  grouping.get());
  }

  template<typename scalar_type>
  void collective_sum_interface<MPI_LIBRARY>::sum(LIN_ALG::vector<scalar_type, LIN_ALG::CPU>& f)
  {
    ////cout << __PRETTY_FUNCTION__ << endl;

    LIN_ALG::vector<scalar_type, LIN_ALG::CPU> F("F", f.get_current_size());

    MPI_Allreduce(&f[0],
                  &F[0],
                  type_map_interface<MPI_LIBRARY, scalar_type>::factor()*f.get_current_size(),
                  type_map_interface<MPI_LIBRARY, scalar_type>::value(),
                  MPI_SUM,
                  grouping.get());

    for(int i=0;i<F.get_current_size();i++)
      f[i] = F[i];

    for(int i=0;i<F.get_current_size();i++)
      if(f[i]!=f[i])
	throw std::logic_error(__FUNCTION__);
  }

  template<typename scalar_type>
  void collective_sum_interface<MPI_LIBRARY>::sum(LIN_ALG::matrix<scalar_type, LIN_ALG::CPU>& f)
  {
    LIN_ALG::matrix<scalar_type, LIN_ALG::CPU> F("F", f.get_current_size(), f.get_global_size());

    assert(f.get_global_size().first  == F.get_global_size().first);
    assert(f.get_global_size().second == F.get_global_size().second);

    int Nr = f.get_global_size().first;
    int Nc = f.get_global_size().second;

    MPI_Allreduce(&f(0,0),
                  &F(0,0),
                  type_map_interface<MPI_LIBRARY, scalar_type>::factor()*Nr*Nc,
                  type_map_interface<MPI_LIBRARY, scalar_type>::value(),
                  MPI_SUM,
                  grouping.get());

    for(int j=0; j<F.get_current_size().second; j++)
      for(int i=0; i<F.get_current_size().first; i++)
        f(i,j) = F(i,j);
  }

  template<typename some_type>
  void collective_sum_interface<MPI_LIBRARY>::sum_and_average(some_type& obj, int size)
  {
    //cout << __PRETTY_FUNCTION__ << endl;

    sum(obj);

    double one_over_N = 1./(size*grouping.get_Nr_threads());

    obj *= one_over_N;
  }

  template<typename scalar_type, class domain>
  void collective_sum_interface<MPI_LIBRARY>::average_and_compute_stddev(FUNC_LIB::function<scalar_type, domain>& f_mean, 
									 FUNC_LIB::function<scalar_type, domain>& f_stddev, 
									 size_t size)
  {
//     if(grouping.get_id()==0)
//       cout << "\n\t\t average_and_compute_stddev " << f_mean.get_name() << "\t" /*<< print_time()*/ << "\n\n";

    scalar_type factor = 1./(size*grouping.get_Nr_threads());

    FUNC_LIB::function<scalar_type, domain> f_sum ("f-sum");
    FUNC_LIB::function<scalar_type, domain> f_diff("f-diff");

    {
      sum(f_mean, f_sum);
    }

    {
      f_sum *= factor;
      
      for(int i=0; i<f_sum.size(); i++)
	f_diff(i) = (f_mean(i) - f_sum(i))*(f_mean(i) - f_sum(i));
      
      for(int i=0; i<f_sum.size(); i++)
	f_mean(i) = f_sum(i);
    }

    {
      sum(f_diff, f_stddev);
    }

    {
      f_stddev *= factor;

      for(int i=0; i<f_sum.size(); i++)
	f_stddev(i) = std::sqrt(f_stddev(i));

      f_stddev /= std::sqrt(grouping.get_Nr_threads());
    }
  }

  template<typename scalar_type, class domain>
  void collective_sum_interface<MPI_LIBRARY>::average_and_compute_stddev(FUNC_LIB::function<std::complex<scalar_type>, domain>& f_mean, 
									 FUNC_LIB::function<std::complex<scalar_type>, domain>& f_stddev, 
									 size_t size)
  {
//     if(grouping.get_id()==0)
//       cout << "\n\t\t average_and_compute_stddev " << f_mean.get_name() << "\t" /*<< print_time()*/ << "\n\n";

    scalar_type factor = 1./(size*grouping.get_Nr_threads());

    FUNC_LIB::function<std::complex<scalar_type>, domain> f_sum ("f-sum");
    FUNC_LIB::function<std::complex<scalar_type>, domain> f_diff("f-diff");

    {
      sum(f_mean, f_sum);
    }

    {
      f_sum *= factor;
      
      for(int i=0; i<f_sum.size(); i++){
	f_diff(i).real(real(f_mean(i) - f_sum(i))*real(f_mean(i) - f_sum(i)));
	f_diff(i).imag(imag(f_mean(i) - f_sum(i))*imag(f_mean(i) - f_sum(i)));
      }

      for(int i=0; i<f_sum.size(); i++)
	f_mean(i) = f_sum(i);
    }

    {
      sum(f_diff, f_stddev);
    }

    {
      f_stddev *= factor;

      for(int i=0; i<f_sum.size(); i++){
	f_stddev(i).real(std::sqrt(real(f_stddev(i))));
	f_stddev(i).imag(std::sqrt(imag(f_stddev(i))));
      }

      f_stddev /= std::sqrt(grouping.get_Nr_threads());
    }
  }

}

#endif
