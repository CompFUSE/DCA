# test/unit/phys/dca_loop

# Disorder unit tests are only built with the disorder option, and never for a GPU build.
if(DCA_WITH_DISORDER AND NOT DCA_WITH_CUDA AND NOT DCA_WITH_HIP)
  add_subdirectory(disorder)
endif()
