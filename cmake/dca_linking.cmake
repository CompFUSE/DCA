

#link the correct gpu runtime library
function(dca_gpu_runtime_link target_name)
  if(DCA_HAVE_HIP)
    target_link_libraries(${target_name} PUBLIC hip::host)
    message("linking target ${target_name} to hip::host")
  elseif(DCA_HAVE_CUDA)
    target_lihk_libraries(${target_name} PUBLIC CUDA::cudart)
  endif()
endfunction()

#link the correct gpu runtime library
function(dca_gpu_blas_link target_name)
  if(DCA_HAVE_HIP)
    target_link_libraries(${target_name} PUBLIC roc::hipblas)
    message("linking target ${target_name} to roc::hipblas")
  elseif(DCA_HAVE_CUDA)
    target_lihk_libraries(${target_name} PUBLIC CUDA::cublas)
  endif()
endfunction()

function(dca_gpu_device_link target_name)
  if(DCA_HAVE_HIP)
    set_target_properties( ${target_name} PROPERTIES LINKER_LANGUAGE "HIP")
    set_target_properties( ${target_name}
      PROPERTIES HIP_SEPARABLE_COMPILATION ON)
    target_link_libraries(${target_name} PRIVATE hip::device)
  elseif(DCA_HAVE_CUDA)
    set_target_properties( ${target_name}
      PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
   # target_compile_definitions(lapack_kernels PRIVATE DCA_HAVE_CUDA)
  endif()
endfunction()
