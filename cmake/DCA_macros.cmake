################################################################################
# initializa a global variable to empty 'definitions' string
################################################################################
set_property(GLOBAL PROPERTY DCA_CONFIG_DEFINITIONS "")

################################################################################
# function to add a definition to the global 'definitions' string
################################################################################
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
# Configure the header to include all compile definitions
################################################################################
function(write_global_definitions_file filename)
  get_property(DCA_CONFIG_DEFINITIONS_VAR GLOBAL PROPERTY DCA_CONFIG_DEFINITIONS)
  
  list(SORT DCA_CONFIG_DEFINITIONS_VAR)
  list(REMOVE_DUPLICATES DCA_CONFIG_DEFINITIONS_VAR)
  list(REMOVE_ITEM DCA_CONFIG_DEFINITIONS_VAR "")

  set(dca_config_defines "// Generated from CMake definitons\n\n")
  foreach(def ${DCA_CONFIG_DEFINITIONS_VAR})
    set(dca_config_defines "${dca_config_defines}#define ${def} ${${def}_define}\n")
  endforeach()

  # Generate a defines.hpp to be used in the build directory ...
  set(DCA_DEFINES_PREFIX ${DCA_BUILD_PREFIX})
  configure_file("${PROJECT_SOURCE_DIR}/cmake/templates/config_defines.hpp.in"
                 "${PROJECT_BINARY_DIR}/${filename}"
                 @ONLY)

  # Generate a defines.hpp to be used in the install directory ...
  set(DCA_DEFINES_PREFIX ${DCA_PREFIX})
  configure_file("${PROJECT_SOURCE_DIR}/cmake/templates/config_defines.hpp.in"
                 "${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/${filename}"
                 @ONLY)
endfunction()
