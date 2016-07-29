//-*-C++-*-

#ifndef BASIC_CUDA_FUNCTIONS_H
#define BASIC_CUDA_FUNCTIONS_H

#include <iostream>
#include <cuda_runtime.h>

void print_device();

void print_device_info();

void initialize_magma();

#ifdef DEBUG_CUDA

bool cuda_check_for_errors(std::string function_name, std::string file_name, int line);

bool cuda_check_for_errors_bgn(std::string function_name, std::string file_name, int line);

bool cuda_check_for_errors_end(std::string function_name, std::string file_name, int line);

#endif 

cudaDeviceProp& get_device_properties();

int get_number_of_threads();
int get_number_of_blocks(int n);
int get_number_of_blocks(int n, int N_th);

void synchronize_devices();

#endif
