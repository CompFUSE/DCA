//-*-C++-*-

#ifndef TIME_DOMAIN_LEFT_ORIENTED_H
#define TIME_DOMAIN_LEFT_ORIENTED_H

/*! 
 *  \ingroup DOMAINS
 *
 *  \author Peter Staar
 *  \brief  This class parametrizes the time-domain with N intervals, but it leaves out the right edge.
 */
struct time_domain_left_oriented
{
  typedef double                    element_type;

  typedef time_domain_left_oriented this_type;
  
  static int& get_size(){
    static int SIZE = time_domain::get_size()-2;
    return SIZE;
  }
  
  static std::vector<element_type>& get_elements(){
    static std::vector<element_type> elements(get_size(), 0);
    return elements;
  }

  template<class parameters_t>
  static void initialize(parameters_t& parameters){

    for(int t_ind=0; t_ind<time_domain::get_size()/2-1; t_ind++)      
      time_domain_left_oriented::get_elements()[t_ind] = time_domain::get_elements()[t_ind];

    for(int t_ind=time_domain::get_size()/2; t_ind<time_domain::get_size()-1; t_ind++)
      time_domain_left_oriented::get_elements()[t_ind-1] = t_ind==time_domain::get_size()/2? 0 : time_domain::get_elements()[t_ind];     
  }

};

#endif
