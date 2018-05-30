# JSON

add_library(json STATIC
  json_parser/json_character_mapper.cpp
  json_parser/json_context.cpp
  json_parser/json_enumerations.cpp
  json_parser/json_leaf.cpp
  json_parser/json_mode_stack.cpp
  json_parser/json_operators.cpp
  json_parser/json_translation_table.cpp
  json_parser/whatever.cpp
  json_reader.cpp
  json_writer.cpp)
