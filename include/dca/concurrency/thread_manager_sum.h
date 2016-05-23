//-*-C++-*-

#ifndef THREAD_MANAGER_SUM_H
#define THREAD_MANAGER_SUM_H

/*!
 *  \author Raffaele Solca' <rasolca@phys.ethz.ch>
 *  \author  Peter Staar <staarp@phys.ethz.ch>   
 */
template<typename concurrency_t>
class thread_manager_sum
{
public:

  thread_manager_sum(concurrency_t& concurrency_ref);
  ~thread_manager_sum();

  template<typename domain_t>
  std::pair<int, int> get_bounds(domain_t& dmn);
    
  template<typename T>
  bool sum_and_check(T& result);

private:
    
  concurrency_t& concurrency_obj;
};

template<typename concurrency_t>
thread_manager_sum<concurrency_t>::thread_manager_sum(concurrency_t& concurrency_ref):
  concurrency_obj(concurrency_ref)
{}

template<typename concurrency_t>
thread_manager_sum<concurrency_t>::~thread_manager_sum()
{}

template<typename concurrency_t>
template<typename domain_t>
std::pair<int, int> thread_manager_sum<concurrency_t>::get_bounds(domain_t& dmn)
{
  return concurrency_obj.get_bounds(dmn);
}

template<typename concurrency_t>
template<typename T>
bool thread_manager_sum<concurrency_t>::sum_and_check(T& result)
{
  concurrency_obj.sum(result);
  return true;    
}

#endif
