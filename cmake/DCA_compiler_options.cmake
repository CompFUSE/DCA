################################################################################
# Compiler options and tweaks
#
# TODO: - Use target_compile_options().
#       - Use set_property(TARGET target PROPERTY CXX_STANDARD 11) instead of
#         -std=c++11.
################################################################################

set(WARNINGS "")
set(WARNINGS_RELEASE "-Wall -Wextra -Wpedantic -Werror -Wno-sign-compare -Wno-unused-parameter")
set(FLAGS "-std=c++14")  # -funroll-loops -finline-functions
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNINGS} ${FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${WARNINGS_RELEASE}")
