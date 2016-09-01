//-*-C++-*-

#ifndef LINALG_SWAP_GPU_H
#define LINALG_SWAP_GPU_H

namespace LIN_ALG {
namespace GPU_KERNEL {
void scale_many_rows(int Nc, int Ni, const int * r_i, const double* alpha, double* A, int LD,
                     int thread_id, int stream_id);
}
}

#endif
