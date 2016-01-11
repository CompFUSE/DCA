//-*-C++-*-

#ifndef FREQUENCY_DOMAIN_COMPACT_H
#define FREQUENCY_DOMAIN_COMPACT_H

namespace DCA
{
  enum VERTEX_FREQUENCY_NAME {COMPACT, EXTENDED, COMPACT_POSITIVE, EXTENDED_POSITIVE, COMPACT_SORTED, EXTENDED_SORTED,
                              EXTENDED_BOSONIC, EXTENDED_FERMIONIC, CORE_SORTED, HIGH_FREQUENCY_SORTED};

  std::string to_str(VERTEX_FREQUENCY_NAME NAME)
  {
    switch(NAME)
      {
      case COMPACT:
        return "COMPACT";
        break;

      case EXTENDED:
        return "EXTENDED";
        break;

      case COMPACT_POSITIVE:
        return "COMPACT_POSITIVE";
        break;

      case EXTENDED_POSITIVE:
        return "EXTENDED_POSITIVE";
        break;

      case EXTENDED_BOSONIC:
        return "EXTENDED_BOSONIC";
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }

    return "NONE";
  }

  /*!
   * \class  vertex_frequency_domain
   *
   * \author Peter Staar
   *
   * \brief  This class implements various types and orderings of the matsubara-frequency domain, via a templates.
   */
  template<VERTEX_FREQUENCY_NAME NAME>
  class vertex_frequency_domain
  {
  public:

    const static int DIMENSION = 1;

    typedef double scalar_type;
    typedef double element_type;

  public:

    static int&                 get_size();
    static std::string          get_name();

    static scalar_type* get_basis();
    static scalar_type* get_inverse_basis();

    static std::vector<double>& get_elements();

    static std::vector<int>&    get_corresponding_frequency_domain_index();

    template<IO::FORMAT DATA_FORMAT>
    static void write(IO::writer<DATA_FORMAT>& writer);

    template<typename parameters_t>
    static void initialize(parameters_t& parameters);
  };

  template<VERTEX_FREQUENCY_NAME NAME>
  int& vertex_frequency_domain<NAME>::get_size()
  {
    static int size;
    return size;
  }

  template<VERTEX_FREQUENCY_NAME NAME>
  std::string vertex_frequency_domain<NAME>::get_name()
  {
    static std::string name = "vertex-frequency-domain (" + to_str(NAME) + ")";
    return name;
  }

  template<VERTEX_FREQUENCY_NAME NAME>
  typename vertex_frequency_domain<NAME>::scalar_type* vertex_frequency_domain<NAME>::get_basis()
  {
    static scalar_type basis[DIMENSION];
    return basis;
  }

  template<VERTEX_FREQUENCY_NAME NAME>
  typename vertex_frequency_domain<NAME>::scalar_type* vertex_frequency_domain<NAME>::get_inverse_basis()
  {
    static scalar_type inv_basis[DIMENSION];
    return inv_basis;
  }

  template<VERTEX_FREQUENCY_NAME NAME>
  std::vector<typename vertex_frequency_domain<NAME>::element_type>& vertex_frequency_domain<NAME>::get_elements()
  {
    static std::vector<element_type> elements;
    return elements;
  }

  template< VERTEX_FREQUENCY_NAME NAME>
  std::vector<int>& vertex_frequency_domain<NAME>::get_corresponding_frequency_domain_index()
  {
    static std::vector<int> elements;
    return elements;
  }

  template<VERTEX_FREQUENCY_NAME NAME>
  template<IO::FORMAT DATA_FORMAT>
  void vertex_frequency_domain<NAME>::write(IO::writer<DATA_FORMAT>& writer)
  {
    writer.open_group(get_name());

    writer.execute("elements", get_elements());

    writer.close_group();
  }

  template<>
  template<typename parameters_t>
  void vertex_frequency_domain<COMPACT>::initialize(parameters_t& parameters)
  {
    get_basis()        [0] = (2.*M_PI)/parameters.get_beta();
    get_inverse_basis()[0] = parameters.get_beta()/(2.*M_PI);

    get_size() = 2*parameters.get_tp_fermionic_frequencies();//wn_c();

    get_elements().resize(get_size());

    for(int l=0; l<parameters.get_tp_fermionic_frequencies(); l++)
      {
        get_elements()[get_size()/2  +l] =  M_PI/parameters.get_beta()*(1+2*l);
        get_elements()[get_size()/2-1-l] = -M_PI/parameters.get_beta()*(1+2*l);
      }

    get_corresponding_frequency_domain_index().resize(get_size(),-1);

    std::vector<double>& wn = frequency_domain::get_elements();

    for(int i=0; i<get_size(); i++)
      for(size_t j=0; j<wn.size(); j++)
        if(std::fabs(wn[j]-get_elements()[i])<1.e-6)
          get_corresponding_frequency_domain_index()[i] = j;

    for(int i=0; i<get_size(); i++)
      if(get_corresponding_frequency_domain_index()[i] == -1 ||
         std::fabs(wn[get_corresponding_frequency_domain_index()[i]]-get_elements()[i])>1.e-6)
        throw std::logic_error(__FUNCTION__);
  }

  template<>
  template<typename parameters_t>
  void vertex_frequency_domain<COMPACT_POSITIVE>::initialize(parameters_t& parameters)
  {
    get_basis()        [0] = (2.*M_PI)/parameters.get_beta();
    get_inverse_basis()[0] = parameters.get_beta()/(2.*M_PI);

    get_size() = parameters.get_tp_fermionic_frequencies();//wn_c();

    get_elements().resize(get_size());

    for(int l=0; l<parameters.get_tp_fermionic_frequencies(); l++)
      get_elements()[l] =  M_PI/parameters.get_beta()*(1+2*l);

    get_corresponding_frequency_domain_index().resize(get_size(),-1);

    std::vector<double>& wn = frequency_domain::get_elements();

    for(int i=0; i<get_size(); i++)
      for(size_t j=0; j<wn.size(); j++)
        if(std::fabs(wn[j]-get_elements()[i])<1.e-6)
          get_corresponding_frequency_domain_index()[i] = j;

    for(int i=0; i<get_size(); i++)
      if(get_corresponding_frequency_domain_index()[i] == -1||
         std::fabs(wn[get_corresponding_frequency_domain_index()[i]]-get_elements()[i])>1.e-6)
        throw std::logic_error(__FUNCTION__);

    assert(get_elements().back() == vertex_frequency_domain<COMPACT>::get_elements().back());
  }

  template<>
  template<typename parameters_t>
  void vertex_frequency_domain<EXTENDED>::initialize(parameters_t& parameters)
  {
    get_basis()        [0] = (2.*M_PI)/parameters.get_beta();
    get_inverse_basis()[0] = parameters.get_beta()/(2.*M_PI);

    get_size() = 2*(parameters.get_tp_fermionic_frequencies() + abs(parameters.get_w_channel()));

    get_elements().resize(get_size());

    for(int l=0; l<get_size()/2; l++)
      {
        get_elements()[get_size()/2  +l] =  M_PI/parameters.get_beta()*(1+2*l);
        get_elements()[get_size()/2-1-l] = -M_PI/parameters.get_beta()*(1+2*l);
      }

    get_corresponding_frequency_domain_index().resize(get_size(),-1);

    std::vector<double>& wn = frequency_domain::get_elements();

    for(int i=0; i<get_size(); i++)
      for(size_t j=0; j<wn.size(); j++)
        if(std::fabs(wn[j]-get_elements()[i])<1.e-6)
          get_corresponding_frequency_domain_index()[i] = j;

    for(int i=0; i<get_size(); i++)
      if(get_corresponding_frequency_domain_index()[i] == -1 ||
         std::fabs(wn[get_corresponding_frequency_domain_index()[i]]-get_elements()[i])>1.e-6)
        throw std::logic_error(__FUNCTION__);
  }

  template<>
  template<typename parameters_t>
  void vertex_frequency_domain<EXTENDED_POSITIVE>::initialize(parameters_t& parameters)
  {
    get_basis()        [0] = (2.*M_PI)/parameters.get_beta();
    get_inverse_basis()[0] = parameters.get_beta()/(2.*M_PI);

    get_size() = parameters.get_tp_fermionic_frequencies() + abs(parameters.get_w_channel());

    get_elements().resize(get_size());

    for(int l=0; l<get_size(); l++)
      get_elements()[l] =  M_PI/parameters.get_beta()*(1+2*l);

    get_corresponding_frequency_domain_index().resize(get_size(),-1);

    std::vector<double>& wn = frequency_domain::get_elements();

    for(int i=0; i<get_size(); i++)
      for(size_t j=0; j<wn.size(); j++)
        if(std::fabs(wn[j]-get_elements()[i])<1.e-6)
          get_corresponding_frequency_domain_index()[i] = j;

    for(int i=0; i<get_size(); i++)
      if(get_corresponding_frequency_domain_index()[i] == -1 ||
         std::fabs(wn[get_corresponding_frequency_domain_index()[i]]-get_elements()[i])>1.e-6)
        throw std::logic_error(__FUNCTION__);

    assert(get_elements().back() == vertex_frequency_domain<EXTENDED>::get_elements().back());
  }

  template<>
  template<typename parameters_t>
  void vertex_frequency_domain<EXTENDED_BOSONIC>::initialize(parameters_t& parameters)
  {
    get_basis()        [0] = (2.*M_PI)/parameters.get_beta();
    get_inverse_basis()[0] = parameters.get_beta()/(2.*M_PI);

    get_size() = 2*parameters.get_sp_bosonic_frequencies()+1;

    get_elements().resize(get_size(), -2.*M_PI/parameters.get_beta()*int(get_size()/2));

    for(int l=0; l<get_size(); l++){
      get_elements()[l] += l*2.*M_PI/parameters.get_beta();
      //cout << get_elements()[l] << endl;
    }
  }

}

#endif
