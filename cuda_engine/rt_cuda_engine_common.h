#ifndef RT_CUDA_ENGINE_COMMON_H
#define RT_CUDA_ENGINE_COMMON_H

#include <cuda.h>
#include <cuda_runtime.h>
#include "../common-lib/num_defs.h"

// helper functions declaration from rt_cuda_engine.cu source file

cudaError_t cuda_kernel_props(size_t n_elems, int *N_cuda_blocks, int *N_cuda_threads);
cudaError_t cudaMallocChunk(size_t required_Nelems, size_t *allocated_Nelems,
                            int N_real_t, bool flag_type, ...);
void cudaFreeChunk(int N, ...);
//void cudaFreeChunk(int N_real_t, bool flag_type, ...);
cudaError_t cudaMemcpyChunk(size_t Nelems, cudaMemcpyKind kind,
                            int N_real_t, bool flag_type, ... );




#endif // RT_CUDA_ENGINE_COMMON_H
