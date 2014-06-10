//-*-C++-*-

#ifndef PROFILING_PARAMETERS_H
#define PROFILING_PARAMETERS_H

/*!
 *   \ingroup  PARAMETERS
 *
 *   \author   Peter Staar
 *   \brief    ...
 */
template<class concurrency_type>
class profiling_parameters 
{
public:

#ifdef NO_PROFILING
  typedef PROFILER::no_profiling_mode profiler_type;
#endif
  
#ifdef COUNTING_PROFILING
  typedef PROFILER::time_event      <size_t>     event_type;
  typedef PROFILER::CountingProfiler<event_type> profiler_type;
#endif

#ifdef DCA_CRAYPAT_PROFILING
  typedef size_t counter_type;
  typedef dca::CraypatProfiler<counter_type> profiler_type;
#endif

#ifdef DCA_PAPI_PROFILING
  typedef PROFILER::papi_and_time_event<long long>  event_type;
  typedef PROFILER::CountingProfiler   <event_type> profiler_type;
#endif
};

#endif
