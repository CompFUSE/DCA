// Copyright (C) 2010 Philipp Werner
//
// Integrated into DCA++ by Peter Staar (taa@zurich.ibm.com) and Bart Ydens.

#ifndef HYBRIDIZATION_HILBERT_DOMAIN_H
#define HYBRIDIZATION_HILBERT_DOMAIN_H
namespace QMC 
{
  class Hybridization_Hilbert_domain 
  {
  public:

    typedef int element_type;

  public:

    static void initialize();

    static int&              get_size();
    static std::vector<int>& get_elements();

  private:
    
    static bool initialized;

  };

  bool Hybridization_Hilbert_domain::initialized = false;

  void Hybridization_Hilbert_domain::initialize()
  {
    Hybridization_Hilbert_space::initialize_Hilbert_space(3);

    Hybridization_symmetries<PARTICLE_NUMBER>::execute_symmetrie();
    Hybridization_symmetries<Sz>             ::execute_symmetrie();
    Hybridization_symmetries<TOTAL_MOMENTUM> ::execute_symmetrie();

    Hybridization_Hilbert_space::sort();

    initialized = true;
  }

  int& Hybridization_Hilbert_domain::get_size()  
  {
    if(!initialized)
      initialize();
    static int size = Hybridization_Hilbert_space::size();
    return size;
  }

  std::vector<int>& Hybridization_Hilbert_domain::get_elements()
  {
    if(!initialized)
      initialize();
    static std::vector<int>& v = Hybridization_Hilbert_space::sizes();
    return v;
  }
}
#endif
