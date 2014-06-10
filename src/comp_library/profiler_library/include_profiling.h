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

#include "time_file_name_changed.h"
#include "time_events.h"

#ifdef NO_PROFILING

#include "null_profiler.h"	

#endif

#ifdef COUNTING_PROFILING

#include "counting_profiler.h"

#endif

#ifdef DCA_CRAYPAT_PROFILING

#include "CraypatProfiler.h"

#endif

#ifdef DCA_PAPI_PROFILING

#include "papi_events.h"
#include "counting_profiler.h"

#endif
