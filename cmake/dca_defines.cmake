################################################################################
# Author: John Biddiscombe (john.biddiscombe@cscs.ch)
#         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# Contains global variables with lists of preprocessor definitions.
# DCA_CONFIG_DEFINITIONS holds preprocessor definitions to configure the applications, e.g.
# DCA_WITH_CT_AUX.
# DCA_HAVES_DEFINITIONS holds preprocessor defintions about availability of certain libraries, e.g.
# MPI or Pthreads.

################################################################################
# Initialize global variables to empty strings.
set_property(GLOBAL PROPERTY DCA_CONFIG_DEFINITIONS "")
set_property(GLOBAL PROPERTY DCA_HAVES_DEFINITIONS "")

################################################################################
# Adds a definition to the global 'config definitions' string.
function(dca_add_config_define definition)
  if(ARGN)
    set_property(GLOBAL APPEND PROPERTY DCA_CONFIG_DEFINITIONS "${definition} ${ARGN}")
  else()
    set_property(GLOBAL APPEND PROPERTY DCA_CONFIG_DEFINITIONS "${definition}")
  endif()
   # debugging only
   # get_property(DCA_CONFIG_DEFINITIONS_VAR GLOBAL PROPERTY DCA_CONFIG_DEFINITIONS)
   # message("Config vars are\n" ${DCA_CONFIG_DEFINITIONS_VAR})
endfunction()

################################################################################
# Adds a definition to the global 'have definitions' string.
function(dca_add_haves_define definition)
  if(ARGN)
    set_property(GLOBAL APPEND PROPERTY DCA_HAVES_DEFINITIONS "${definition} ${ARGN}")
  else()
    set_property(GLOBAL APPEND PROPERTY DCA_HAVES_DEFINITIONS "${definition}")
  endif()
endfunction()

################################################################################
# Generates in the build directory the config_defines.hpp that contains all 'config preprocessor
# definitions'.
function(dca_write_config_definitions_file)
  get_property(DCA_CONFIG_DEFINITIONS_VAR GLOBAL PROPERTY DCA_CONFIG_DEFINITIONS)

  list(SORT DCA_CONFIG_DEFINITIONS_VAR)
  list(REMOVE_DUPLICATES DCA_CONFIG_DEFINITIONS_VAR)
  list(REMOVE_ITEM DCA_CONFIG_DEFINITIONS_VAR "")

  set(dca_config_defines "")
  foreach(def ${DCA_CONFIG_DEFINITIONS_VAR})
    set(dca_config_defines "${dca_config_defines}#define ${def} ${${def}_define}\n")
  endforeach()

  configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/config_defines.hpp.in"
    "${CMAKE_BINARY_DIR}/include/dca/config/config_defines.hpp"
    @ONLY)
endfunction()

################################################################################
# Generates in the build directory the haves_defines.hpp that contains all 'haves preprocessor
# definitions'.
function(dca_write_haves_definitions_file)
  get_property(DCA_HAVES_DEFINITIONS_VAR GLOBAL PROPERTY DCA_HAVES_DEFINITIONS)

  list(SORT DCA_HAVES_DEFINITIONS_VAR)
  list(REMOVE_DUPLICATES DCA_HAVES_DEFINITIONS_VAR)
  list(REMOVE_ITEM DCA_HAVES_DEFINITIONS_VAR "")

  set(dca_haves_defines "")
  foreach(def ${DCA_HAVES_DEFINITIONS_VAR})
    string(CONCAT dca_haves_defines
        "${dca_haves_defines}"
        "#ifndef ${def}\n"
        " #define ${def} ${${def}_define}\n"
        "#endif\n\n")
  endforeach()

  configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/haves_defines.hpp.in"
    "${CMAKE_BINARY_DIR}/include/dca/config/haves_defines.hpp"
    @ONLY)
endfunction()
