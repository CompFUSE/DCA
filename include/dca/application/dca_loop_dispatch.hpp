
// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// Dispatches based on runtime selectable template parameters.
// Obviously this must be limited as it results in combinatorial amounts of
// code being compiled.

#ifndef DCA_APPLICATION_DCA_LOOP_DISPATCH_HPP
#define DCA_APPLICATION_DCA_LOOP_DISPATCH_HPP
#include "dca/config/dca.hpp"

template <dca::DistType DT>
class DCALoopDispatch {
public:
  void operator()(ParametersType& parameters, Concurrency& concurrency) {
    DcaDataType<DT> dca_data(parameters);
    dca_data.initialize();
    DcaLoopType<DT> dca_loop(parameters, dca_data, concurrency);
    {
      Profiler profiler(__FUNCTION__, __FILE__, __LINE__);

      // Create and initialize the DCA data object.
    
      try {
	dca_loop.initialize();
      } catch (const std::exception& exc) {
	std::cout << "unhandled exception in dca_loop.initialize(): " << exc.what() << std::endl;
	throw exc;
      }
      try {
	dca_loop.execute();
      } catch (const std::exception& exc) {
	std::cout << "unhandled exception in dca_loop.execute(): " << exc.what() << std::endl;
	throw exc;
      }
      try {
      dca_loop.finalize();
      Profiler::stop(concurrency, parameters.get_filename_profiling());
      } catch (const std::exception& exc) {
	std::cout << "unhandled exception in dca_loop.execute(): " << exc.what() << std::endl;
	throw exc;
      }
 
      // Whether writing is on all ranks or not is now controlled at the writer level.
      //      if (concurrency.id() == concurrency.first()) {
      //std::cout << "\nProcessor " << concurrency.id() << " is writing data." << std::endl;
      dca_loop.write();

      std::cout << "\nFinish time: " << dca::util::print_time() << "\n" << std::endl;
      //      }
    }
  }
};

#endif
