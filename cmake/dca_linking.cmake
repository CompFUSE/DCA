
#link the correct gpu runtime library
function(dca_gpu_runtime_link target_name)
  if(DCA_HAVE_HIP)
    target_link_libraries(${target_name} PUBLIC hip::host roc::hipblas roc::hipsparse)
    message("linking target ${target_name} to hip::host")
  elseif(DCA_HAVE_CUDA)
    target_link_libraries(${target_name} PUBLIC CUDA::cudart)
  endif()
endfunction()

#link the correct gpu runtime library
function(dca_gpu_blas_link target_name)
  if(DCA_HAVE_HIP)
    target_link_libraries(${target_name} PUBLIC roc::hipblas roc::hipsparse)
    message("linking target ${target_name} to roc::hipblas")
  elseif(DCA_HAVE_CUDA)
    target_link_libraries(${target_name} PUBLIC CUDA::cublas)
  endif()
endfunction()
	
function(dca_gpu_device_link target_name)
  if(DCA_HAVE_HIP)
    set_target_properties( ${target_name} PROPERTIES LINKER_LANGUAGE "HIP")
    set_target_properties( ${target_name}
      PROPERTIES HIP_SEPARABLE_COMPILATION ON)
    set_target_properties( ${target_name}
      PROPERTIES HIP_RESOLVE_DEVICE_SYMBOLS ON)
    target_link_libraries(${target_name} PRIVATE hip::device roc::hipblas roc::hipsparse roc::rocthrust)
    get_target_property(_srcs ${target_name} SOURCES)
    get_target_property(_src_dir ${target_name} SOURCE_DIR)
    #
    # Mark all cu source files as HIP code.
    foreach(_src IN LISTS _srcs)
        #message("${_src_dir}/${_src}")
        if(_src MATCHES ".*\.cu$")
            set_source_files_properties(${_src} PROPERTIES LANGUAGE HIP)
        endif()
    endforeach()
  elseif(DCA_HAVE_CUDA)
    set_target_properties( ${target_name}
      PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    set_target_properties( ${target_name}
      PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
   # target_compile_definitions(lapack_kernels PRIVATE DCA_HAVE_CUDA)
  endif()
endfunction()
