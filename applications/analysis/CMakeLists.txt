# Builds main_analysis.

if (DCA_BUILD_ANALYSIS)
  add_executable(main_analysis main_analysis.cpp)
  target_include_directories(main_analysis PRIVATE ${DCA_INCLUDE_DIRS})
    if (DCA_HAVE_GPU)
      target_link_libraries(main_analysis PRIVATE ${DCA_GPU_LIBS})
    endif()

    target_link_libraries(main_analysis PUBLIC signals ${DCA_LIBS} simplex_gm_rule)

  message("DCA_MODEL_IS_MATERIAL_LATTICE ${DCA_MODEL_IS_MATERIAL_LATTICE}")
  if (DCA_LATTICE STREQUAL "square")
    #new development of exanded q chi_q_omega
    add_executable(chi_q_omega chi_q_omega.cpp)
    target_compile_definitions(chi_q_omega DCA_HAS_SOFQOMEGA)
    target_include_directories(chi_q_omega PRIVATE ${DCA_INCLUDE_DIRS})

    if (DCA_HAVE_GPU)
      target_link_libraries(chi_q_omega PRIVATE ${DCA_GPU_LIBS})
    endif()

    target_link_libraries(chi_q_omega PUBLIC signals ${DCA_LIBS} simplex_gm_rule)

  endif()

endif()
