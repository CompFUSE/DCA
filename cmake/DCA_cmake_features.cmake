################################################################################
# Author: John Biddiscombe (john.biddiscombe@cscs.ch)
#
# Contains some useful CMake functions.

################################################################################
# Prints a list with a title and the list length.
function(print_list title list_in)
  list(LENGTH list_in LIST_LENGTH)
  message("================ ${title} (${LIST_LENGTH}) ================ ")
  foreach (_item ${list_in})
    message("${_item}")
  endforeach()
  message("================")
endfunction()

################################################################################
# Removes all items that appear in list1 from list2.
function(remove_from_list list1 list2 out_list)
  foreach(_item ${list1})
    list(FIND list2 ${_item} _pos)
    if (NOT _pos STREQUAL "-1")
      list(REMOVE_ITEM list2 ${_item})
    endif()
  endforeach()
  set(${out_list} "${list2}" PARENT_SCOPE)
endfunction()

################################################################################
# Checks if all user supplied CXX features are supported.
# In: list of variable length
function(dca_check_cxx_features)
  set(_feature_list "${ARGN}")
  # print_list("Requested features" "${_feature_list}")

  # Get all CXX features known to CMake.
  get_property(cxx_known_features GLOBAL PROPERTY CMAKE_CXX_KNOWN_FEATURES)
  # list(SORT cxx_known_features)

  # Get all CXX features supported by the compiler.
  set(cxx_compile_features ${CMAKE_CXX_COMPILE_FEATURES})
  # list(SORT cxx_compile_features)

  # (debug) Remove supported features from known features list.
  remove_from_list("${cxx_compile_features}" "${cxx_known_features}" unsupported_features)

  # Remove supported features from requested feature list.
  remove_from_list("${cxx_compile_features}" "${_feature_list}" unsupported_requested)
  list(LENGTH unsupported_requested LENGTH_UNSUPPORTED)
  if (LENGTH_UNSUPPORTED GREATER 0)
    print_list("Requested unsupported CXX features" "${unsupported_requested}")
    print_list("All unsupported CXX features" "${unsupported_features}")
    message(WARNING "This project requires unsupported CXX features\n
                     @@@@@@@@@@@@@@@@@@@@@\n
                     Get a better compiler\n
                     @@@@@@@@@@@@@@@@@@@@@\n")
  endif()
endfunction()
