//-*-C++-*-

#ifndef LINALG_SWAP_GPU_H
#define LINALG_SWAP_GPU_H

namespace LIN_ALG {
namespace GPU_KERNEL {
void scale_many_rows(int Nc, int Ni, const int* r_i, const double* alpha, double* A, int LD,
                     int thread_id, int stream_id);

void many_row_copies(int N_x, int N_i, const int* i_x, const double* x, int inc_x, const int* i_y,
                     double* y, int inc_y, int thread_id, int stream_id);
void many_column_copies(int N_x, int N_i, const int* i_x, const double* x, int inc_x,
                        const int* i_y, double* y, int inc_y, int thread_id, int stream_id);
}
}

#endif
