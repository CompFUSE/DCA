set(TEST_INCLUDES ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR})
set(TEST_LIBRARIES ${DCA_LIBS};statistical_testing)

add_executable(NiO_fullDCA NiO_fullDCA.cpp)
target_link_libraries(NiO_fullDCA ${TEST_LIBRARIES})
target_include_directories(NiO_fullDCA PRIVATE ${TEST_INCLUDES})
target_compile_definitions(NiO_fullDCA PRIVATE TEST_DIRECTORY="${CMAKE_CURRENT_SOURCE_DIR}/")
