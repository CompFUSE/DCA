add_library(json STATIC json_writer.cpp json_reader.cpp details/json_group.cpp
                        details/json_entry.cpp details/util.cpp)

if(DCA_HAVE_HIP)
  target_link_libraries(json PUBLIC hip::host)
endif()
		      
