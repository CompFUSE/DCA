# Profiling

add_library(profiling STATIC events/time.cpp)

find_library(PAPI_LIB papi)
if(PAPI_LIB)
    add_library(papi_profiling STATIC events/papi_and_time_event.cpp)
    target_link_libraries(papi_profiling PUBLIC papi)
    target_link_libraries(profiling PUBLIC papi_profiling)
endif()
