################################################################################
# Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
#

include(CMakeParseArguments)

function(dca_add_library)
  add_library(${ARGV})

  if (DCA_HAVE_CUDA)
    target_compile_definitions(${ARGV0} PRIVATE DCA_HAVE_CUDA)
    target_include_directories(${ARGV0} PUBLIC ${CUDA_TOOLKIT_INCLUDE})
  endif()
  if(DCA_HAVE_MAGMA)
    target_compile_definitions(${ARGV0} PRIVATE DCA_HAVE_MAGMA)
    target_include_directories(${ARGV0} PUBLIC ${MAGMA_INCLUDE_DIR};)
  endif()
endfunction()
