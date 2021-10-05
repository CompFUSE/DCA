/////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright(C) 2021 Advanced Micro Devices, Inc. All rights reserved.
// Copyright(C) 2021 UT-Battelle, LLC
//
// File developed by: Jakub Kurzak, jakurzak@amd.com, Advanced Micro Devices, Inc.
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Jakub Kurzak, jakurzak@amd.com, Advanced Micro Devices, Inc.
//////////////////////////////////////////////////////////////////////////////////////

#ifndef CUDA2HIP_H
#define CUDA2HIP_H

#define CUBLAS_OP_N                     HIPBLAS_OP_N
#define CUBLAS_STATUS_ALLOC_FAILED      HIPBLAS_STATUS_ALLOC_FAILED
#define CUBLAS_STATUS_ARCH_MISMATCH     HIPBLAS_STATUS_ARCH_MISMATCH
#define CUBLAS_STATUS_EXECUTION_FAILED  HIPBLAS_STATUS_EXECUTION_FAILED
#define CUBLAS_STATUS_INTERNAL_ERROR    HIPBLAS_STATUS_INTERNAL_ERROR
#define CUBLAS_STATUS_INVALID_VALUE     HIPBLAS_STATUS_INVALID_VALUE
#define CUBLAS_STATUS_LICENSE_ERROR     HIPBLAS_STATUS_LICENSE_ERROR
#define CUBLAS_STATUS_MAPPING_ERROR     HIPBLAS_STATUS_MAPPING_ERROR
#define CUBLAS_STATUS_NOT_INITIALIZED   HIPBLAS_STATUS_NOT_INITIALIZED
#define CUBLAS_STATUS_NOT_SUPPORTED     HIPBLAS_STATUS_NOT_SUPPORTED
#define CUBLAS_STATUS_SUCCESS           HIPBLAS_STATUS_SUCCESS

#define cublasCgemmBatched      hipblasCgemmBatched
#define cublasCgetrfBatched     hipblasCgetrfBatched
#define cublasCgetriBatched     hipblasCgetriBatched
#define cublasComplex           hipblasComplex
#define cublasCreate            hipblasCreate
#define cublasDestroy           hipblasDestroy
#define cublasDgemmBatched      hipblasDgemmBatched
#define cublasDgetrfBatched     hipblasDgetrfBatched
#define cublasDgetriBatched     hipblasDgetriBatched
#define cublasDoubleComplex     hipblasDoubleComplex
#define cublasGetVersion        hipRuntimeGetVersion
#define cublasHandle_t          hipblasHandle_t
#define cublasSetStream         hipblasSetStream
#define cublasSgemmBatched      hipblasSgemmBatched
#define cublasSgetrfBatched     hipblasSgetrfBatched
#define cublasSgetriBatched     hipblasSgetriBatched
#define cublasStatus_t          hipblasStatus_t
#define cublasZgemmBatched      hipblasZgemmBatched
#define cublasZgetrfBatched     hipblasZgetrfBatched
#define cublasZgetriBatched     hipblasZgetriBatched
#define cusparseHandle_t        hipsparseHandle_t
#define cusparseCreate          hipsparseCreate
#define cusparseDestroy         hipsparseDestroy
#define magma_queue_create_from_cuda_internal  magma_queue_create_from_hip_internal

#define cuComplex                       hipComplex
#define cudaAddressModeClamp            hipAddressModeClamp
#define cudaArray                       hipArray
#define cudaBindTextureToArray          hipBindTextureToArray
#define cudaChannelFormatDesc           hipChannelFormatDesc
#define cudaChannelFormatKindFloat      hipChannelFormatKindFloat
#define cudaCreateChannelDesc           hipCreateChannelDesc
#define cudaDeviceProp                  hipDeviceProp_t
#define cudaDeviceReset                 hipDeviceReset
#define cudaDeviceSynchronize           hipDeviceSynchronize
#define cudaError_t                     hipError_t
#define cudaErrorInvalidValue           hipErrorInvalidValue
#define cudaEvent_t                     hipEvent_t
#define cudaEventCreate                 hipEventCreate
#define cudaEventCreateWithFlags        hipEventCreateWithFlags
#define cudaEventDestroy                hipEventDestroy
#define cudaEventDisableTiming          hipEventDisableTiming
#define cudaEventElapsedTime            hipEventElapsedTime
#define cudaEventQuery                  hipEventQuery
#define cudaEventRecord                 hipEventRecord
#define cudaEventSynchronize            hipEventSynchronize
#define cudaFilterModeLinear            hipFilterModeLinear
#define cudaFree                        hipFree
#define cudaFreeHost                    hipHostFree
#define cudaGetDevice                   hipGetDevice
#define cudaGetDeviceCount              hipGetDeviceCount
#define cudaGetDeviceProperties         hipGetDeviceProperties
#define cudaGetErrorString              hipGetErrorString
#define cudaGetLastError                hipGetLastError
#define cudaPeekAtLastError             hipPeekAtLastError
#define cudaHostAlloc                   hipHostMalloc
#define cudaHostAllocMapped             hipHostMallocMapped
#define cudaHostAllocDefault            hipHostMallocDefault
#define cudaIpcGetMemHandle             hipIpcGetMemHandle
#define cudaIpcMemHandle_t              hipIpcMemHandle_t
#define cudaIpcMemLazyEnablePeerAccess  hipIpcMemLazyEnablePeerAccess
#define cudaIpcOpenMemHandle            hipIpcOpenMemHandle
#define cudaMalloc                      hipMalloc
#define cudaMallocArray                 hipMallocArray
#define cudaMallocManaged               hipMallocManaged
#define cudaMemAdvise                   hipMemAdvise
#define cudaMemAdviseSetAccessedBy      hipMemAdviseSetAccessedBy
#define cudaMemAdviseSetReadMostly      hipMemAdviseSetReadMostly
#define cudaMemAttachGlobal             hipMemAttachGlobal
#define cudaMemcpy                      hipMemcpy
#define cudaMemcpyAsync                 hipMemcpyAsync
#define cudaMemcpyDeviceToDevice        hipMemcpyDeviceToDevice
#define cudaMemcpyDeviceToHost          hipMemcpyDeviceToHost
#define cudaMemcpyHostToDevice          hipMemcpyHostToDevice
#define cudaMemcpyHostToHost            hipMemcpyHostToHost
#define cudaMemset                      hipMemset
#define cudaMemsetAsync                 hipMemsetAsync
#define cudaMemcpyDefault               hipMemcpyDefault
#define cudaMemcpyToArrayAsync          hipMemcpyToArray
#define cudaMemcpyToSymbol              hipMemcpyToSymbol
#define cudaMemcpyToSymbolAsync         hipMemcpyToSymbolAsync
#define cudaMemPrefetchAsync            hipMemPrefetchAsync
#define cudaMemoryType                  hipMemoryType
#define cudaMemoryTypeDevice            hipMemoryTypeDevice
#define cudaMemoryTypeHost              hipMemoryTypeHost
#define cudaPointerAttributes           hipPointerAttribute_t
#define cudaPointerGetAttributes        hipPointerGetAttributes
#define cudaReadModeElementType         hipReadModeElementType
#define cudaSetDevice                   hipSetDevice
#define cudaStream_t                    hipStream_t
#define cudaStreamCreate                hipStreamCreate
#define cudaStreamDestroy               hipStreamDestroy
#define cudaStreamSynchronize           hipStreamSynchronize
#define cudaStreamWaitEvent             hipStreamWaitEvent
#define cudaSuccess                     hipSuccess
#define cuFloatComplex                  hipFloatComplex
#define cuDoubleComplex                 hipDoubleComplex
#define make_cuComplex                  make_hipComplex
#define make_cuDoubleComplex            make_hipDoubleComplex

#define cudaDeviceSetLimit(limit, falue) ;

#endif /* CUDA2HIP_H */
