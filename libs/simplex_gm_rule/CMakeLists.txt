# Simplex GM rule

add_library(simplex_gm_rule STATIC simplex_gm_rule.cpp)

# Silence compiler warnings.
target_compile_options(simplex_gm_rule PRIVATE "-Wno-unused-parameter")
if (CMAKE_COMPILER_IS_GNUCXX)
  target_compile_options(simplex_gm_rule PRIVATE "-Wno-unused-but-set-variable")
endif()

set(SIMPLEX_GM_RULE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
set(SIMPLEX_GM_RULE_LIBRARY simplex_gm_rule CACHE INTERNAL "")
