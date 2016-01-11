//-*-C++-*-

/*!
 *  \defgroup PROFILING  
 *  \ingroup  ALGORITHMS
 */

std::string print_time()
{
  time_t rawtime;
  struct tm * timeinfo;
  
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  
  //return asctime(timeinfo);

  char buffer[80];
  strftime(buffer,80,"%d-%m-%Y %I:%M:%S",timeinfo);
  std::string str(buffer);

  return str;
}

#include "events/time_file_name_changed.h"
#include "events/time_events.h"

#ifdef NO_PROFILING

#include "profilers/null_profiler.h"	

#endif

#ifdef COUNTING_PROFILING

#include "profilers/counting_profiler.h"

#endif

#ifdef DCA_CRAYPAT_PROFILING

#include "CraypatProfiler.h"

#endif

#ifdef DCA_PAPI_PROFILING

#include "events/papi_events.h"
#include "profilers/counting_profiler.h"

#endif
