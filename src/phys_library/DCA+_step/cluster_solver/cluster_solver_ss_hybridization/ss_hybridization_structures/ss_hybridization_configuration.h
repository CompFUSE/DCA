//-*-C++-*-

#ifndef SS_HYBRIDIZATION_CONFIGURATION_H
#define SS_HYBRIDIZATION_CONFIGURATION_H

namespace DCA 
{
  /*! 
   *  \brief   This class organizes the configuration space in the single-site hybridization QMC.
   *  \author  Peter Staar
   *  \author  Bart Ydens
   *  \version 1.0
   *
   */ 
  class SS_CT_HYB_configuration
  {
#include "type_definitions.h" 

  public:

    typedef SS_CT_HYB_configuration this_type;

    typedef std::vector<Hybridization_vertex> orbital_configuration_type;

  public:

    SS_CT_HYB_configuration();

    ~SS_CT_HYB_configuration();

    int size();

    void initialize();

    orbital_configuration_type& get_vertices(int i);
    bool&                       get_full_line(int i);

    void copy_from(this_type& other_configuration);

    void print();    

  private:

    function<orbital_configuration_type, nu> vertices;
    function<bool                      , nu> has_full_line;

    int N_spin_orbitals;
  };

  SS_CT_HYB_configuration::SS_CT_HYB_configuration():
    vertices("SS_CT_HYB_vertices"),
    has_full_line("SS_CT_HYB_lines")
  {
    N_spin_orbitals = s::dmn_size()*b::dmn_size();

    for(int i=0; i<N_spin_orbitals; i++)
      has_full_line(i) = false;	
  }
  
  SS_CT_HYB_configuration::~SS_CT_HYB_configuration()
  {}

  int SS_CT_HYB_configuration::size()
  {
    int size=0;

    for(int l=0; l<N_spin_orbitals; l++)
      size += get_vertices(l).size();

    return size;
  }

  void SS_CT_HYB_configuration::initialize()
  {
    for(int i=0; i<has_full_line.size(); i++)
      has_full_line(i) = false;

    for(int i=0; i<vertices.size(); i++)
      vertices(i).resize(0);
  }

  std::vector<Hybridization_vertex>& SS_CT_HYB_configuration::get_vertices(int i)
  {
    return vertices(i);
  }

  bool& SS_CT_HYB_configuration::get_full_line(int i)
  {
    return has_full_line(i);
  }

  void SS_CT_HYB_configuration::copy_from(this_type& other_configuration)
  {
    for(int l=0; l<nu::dmn_size(); l++){

      orbital_configuration_type& other_vertices = other_configuration.get_vertices(l);

      vertices(l).resize(other_vertices.size());

      for(int i=0; i<other_vertices.size(); i++){
	vertices(l)[i].set_t_end  (other_vertices[i].t_end  ());
	vertices(l)[i].set_t_start(other_vertices[i].t_start());
      }
    }

    for(int l=0; l<nu::dmn_size(); l++)
      has_full_line(l) = other_configuration.get_full_line(l);
  }

  void SS_CT_HYB_configuration::print()
  {
    for(int i=0; i<N_spin_orbitals; i++)
      {
	cout << i;
	if(vertices(i).size() >0)
	  for(size_t j=0; j<vertices(i).size(); j++)
	    cout << "\t {" << vertices(i)[j].t_start() << " , " << vertices(i)[j].t_end() << " }"; 
	else if(has_full_line(i))
	  cout << "\t full line";
	cout << "\n";
      }
    cout << "\n";
  } 


}

#endif

