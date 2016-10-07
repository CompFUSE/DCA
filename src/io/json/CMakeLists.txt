# JSON
# TODO: Create one json library as soon as the multiple definition problem is fixed.

SET(JSON_PARSER_SOURCES
  json_parser/json_character_mapper.cpp
  json_parser/json_context.cpp
  json_parser/json_enumerations.cpp
  json_parser/json_leaf.cpp
  json_parser/json_mode_stack.cpp
  json_parser/json_operators.cpp
  json_parser/json_translation_table.cpp
  json_parser/whatever.cpp)

if (DCA_HAVE_MULDEFS)
  add_library(json_reader STATIC json_reader.cpp ${JSON_PARSER_SOURCES})
  add_library(json_writer STATIC json_writer.cpp ${JSON_PARSER_SOURCES})
else()
  add_library(json_reader SHARED json_reader.cpp ${JSON_PARSER_SOURCES})
  add_library(json_writer SHARED json_writer.cpp ${JSON_PARSER_SOURCES})
  target_link_libraries(json_reader ${LAPACK_LIBRARIES})
  target_link_libraries(json_writer ${LAPACK_LIBRARIES})
endif()

target_include_directories(json_reader PUBLIC ${FFTW_INCLUDE_DIR})
target_include_directories(json_writer PUBLIC ${FFTW_INCLUDE_DIR})
