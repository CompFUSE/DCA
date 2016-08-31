# Builds main_analysis.

if (DCA_BUILD_ANALYSIS)
  add_executable(main_analysis main_analysis.cpp)
  target_link_libraries(main_analysis ${DCA_LIBRARIES})
  target_include_directories(main_analysis
    PRIVATE ${DCA_LATTICE_MODELS_FOLDER}
    PRIVATE ${DCA_INCLUDES})
endif()
