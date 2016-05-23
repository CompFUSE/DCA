//-*-C++-*-

#ifndef COLLECTIVE_SUM_INTERFACE_PTHREADS_H
#define COLLECTIVE_SUM_INTERFACE_PTHREADS_H

namespace COMP_LIB
{
  /*!
   *  \author Peter Staar
   */
  template<>
  class collective_sum_interface<POSIX_LIBRARY>
  {
  public:

    template<typename scalar_type>
    static void sum(scalar_type&     value,
		    scalar_type&     result,
		    pthread_mutex_t& mutex);

    template<typename scalar_type, class domain>
    static void sum(FUNC_LIB::function<scalar_type, domain>& f,
		    FUNC_LIB::function<scalar_type, domain>& f_result,
		    pthread_mutex_t&               mutex);

  private:

  };

  template<typename scalar_type>
  void collective_sum_interface<POSIX_LIBRARY>::sum(scalar_type&     value, 
						    scalar_type&     result,
						    pthread_mutex_t& mutex)
  {
    pthread_mutex_lock( &mutex);
    
    value += result;
    
    pthread_mutex_unlock(&mutex);
  }

  template<typename scalar_type, class domain>
  void collective_sum_interface<POSIX_LIBRARY>::sum(FUNC_LIB::function<scalar_type, domain>& f,
						    FUNC_LIB::function<scalar_type, domain>& f_result,
						    pthread_mutex_t&               mutex)
  {
    pthread_mutex_lock( &mutex);
    
    for(int l=0; l<f.size(); l++)
      f_result(l) += f(l);
    
    pthread_mutex_unlock(&mutex);
  }

}

#endif
