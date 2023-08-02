# Profiling

add_library(profiling STATIC events/time.cpp)

if(DCA_PROFILER STREQUAL "PAPI")
	find_library(PAPI_LIB papi)
    add_library(papi_profiling STATIC events/papi_and_time_event.cpp)
	if(PAPI_LIB)
		target_link_libraries(papi_profiling PUBLIC papi)
	else()
		message(WARNING "Unable to find -lpapi. Be sure to include manually")
	endif()
    target_link_libraries(profiling PUBLIC papi_profiling)
endif()
