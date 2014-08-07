#include<cuda.h>
#include<cuda_runtime.h>

#include "../common-lib/num_defs.h"
#include "../common-lib/rt_engine_errors.h"

#include "rt_cuda_engine_common.h"

#include <cstdlib>
#include "iostream"

                    /********************************************
                    *                                           *
                    *     nVIDIA CUDA font-end realization      *
                    *   for TabulatedFunction class functions   *
                    *                                           *
                    ********************************************/


__global__
static void spline_init_kernel(size_t N, real_t *dev_x, real_t *dev_y, real_t *dev_dy2, real_t *dev_u)
{
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;
    real_t sig, p;

    idx += 1;
    while ( idx < (N-1) ) { // skip the first and last element computations
        sig=(dev_x[idx]-dev_x[idx-1])/(dev_x[idx+1]-dev_x[idx-1]);
        p=sig*dev_dy2[idx-1]+2.0;
        dev_dy2[idx]=(sig-1.0)/p;
        dev_u[idx]=(dev_y[idx+1]-dev_y[idx])/(dev_x[idx+1]-dev_x[idx]) - (dev_y[idx]-dev_y[idx-1])/(dev_x[idx]-dev_x[idx-1]);
        dev_u[idx]=(6.0*dev_u[idx]/(dev_x[idx+1]-dev_x[idx-1])-sig*dev_u[idx-1])/p;

        idx += blockDim.x*gridDim.x;
    }

}

// callable host function
// The function computes the second derivatives of the interpolating function at the tabulated x-points
// The algorithm is adopted from Numerical Recipts in C
__host__
RT_engine_error spline_init(size_t N, real_t *x, real_t *y, real_t *dy2)
{
    cudaError_t cuda_err;

    RT_engine_error ret_err;
    int N_cuda_blocks, N_cuda_threads;
    size_t dev_N, N_chunks, rest_N, start_elem;

    real_t *u, *dev_x, *dev_y, *dev_dy2, *dev_u;

    u = (real_t*) malloc(N*sizeof(real_t));
    if ( u == NULL ) {
        return ENGINE_ERROR_BAD_ALLOC;
    }

//    cuda_err = cuda_malloc_vectors(N,&dev_N,XYcXcY_vect_flags,&dev_x,&dev_y,&dev_dy2,&dev_u); // allocate memmory for 4 vectors
    cuda_err = cudaMallocChunk(N,&dev_N,4,false,&dev_x,&dev_y,&dev_dy2,&dev_u);
    if ( cuda_err != cudaSuccess ) {
        free(u);
        return ENGINE_ERROR_BAD_ALLOC;
    }

    cuda_err = cuda_kernel_props(dev_N-2,&N_cuda_blocks,&N_cuda_threads); // "dev_N-2" because of the first and last elements is not computed!!!
    if ( cuda_err != cudaSuccess ) {
        free(u);
//        cuda_free_mem(XYcXcY_vect_flags,dev_x,dev_y,dev_dy2,dev_u);
        cudaFreeChunk(4,dev_x,dev_y,dev_dy2,dev_u);
        return ENGINE_ERROR_FAILED;
    }



    N_chunks = N/dev_N;

    ret_err = ENGINE_ERROR_OK;

    // natural spline
    dy2[0] = 0.0;
    dy2[N-1] = 0.0;
    u[0] = 0.0;


    start_elem = 0;
    for ( size_t i = 0; i < N_chunks; ++i) {
//        cuda_err = cuda_copy_mem(dev_N,start_elem,cudaMemcpyHostToDevice,XYcXcY_vect_flags,x,dev_x,y,dev_y,dy2,dev_dy2,u,dev_u);
        cuda_err = cudaMemcpyChunk(dev_N,cudaMemcpyHostToDevice,4,false,dev_x,x+start_elem,dev_y,y+start_elem,dev_dy2,dy2+start_elem,dev_u,u+start_elem);
        if ( cuda_err != cudaSuccess ) {
            break;
        }

        spline_init_kernel<<<N_cuda_blocks,N_cuda_threads>>>(dev_N,dev_x,dev_y,dev_dy2,dev_u);

//        cuda_err = cuda_copy_mem(dev_N,start_elem,cudaMemcpyDeviceToHost,X_vect_flag | Y_vect_flag,dy2,dev_dy2,u,dev_u);
        cuda_err = cudaMemcpyChunk(dev_N,cudaMemcpyDeviceToHost,4,false,dy2+start_elem,dev_dy2,u+start_elem,dev_u);
        if ( cuda_err != cudaSuccess ) {
            break;
        }

        start_elem += dev_N;
    }

    // process rest of rays
    rest_N = N % dev_N;
    if ( rest_N && (ret_err == ENGINE_ERROR_OK) ) {
        cuda_err = cuda_kernel_props(rest_N-2,&N_cuda_blocks,&N_cuda_threads); // "rest_N-2" because of the first and last elements is not computed!!!
        if ( cuda_err != cudaSuccess ) {
            free(u);
//            cuda_free_mem(XYcXcY_vect_flags,dev_x,dev_y,dev_dy2,dev_u);
            cudaFreeChunk(4,dev_x,dev_y,dev_dy2,dev_u);
            return ENGINE_ERROR_FAILED;
        }

//        cuda_err = cuda_copy_mem(rest_N,start_elem,cudaMemcpyHostToDevice,XYcXcY_vect_flags,x,dev_x,y,dev_y,dy2,dev_dy2,u,dev_u);
        cuda_err = cudaMemcpyChunk(rest_N,cudaMemcpyHostToDevice,4,false,dev_x,x+start_elem,dev_y,y+start_elem,dev_dy2,dy2+start_elem,dev_u,u+start_elem);
        if ( cuda_err != cudaSuccess ) {
            free(u);
//            cuda_free_mem(XYcXcY_vect_flags,dev_x,dev_y,dev_dy2,dev_u);
            cudaFreeChunk(4,dev_x,dev_y,dev_dy2,dev_u);
            return ENGINE_ERROR_FAILED;
        }

        spline_init_kernel<<<N_cuda_blocks,N_cuda_threads>>>(rest_N,dev_x,dev_y,dev_dy2,dev_u);

//        cuda_err = cuda_copy_mem(rest_N,start_elem,cudaMemcpyDeviceToHost,X_vect_flag | Y_vect_flag,dy2,dev_dy2,u,dev_u);
        cuda_err = cudaMemcpyChunk(rest_N,cudaMemcpyDeviceToHost,4,false,dy2+start_elem,dev_dy2,u+start_elem,dev_u);
        if ( cuda_err != cudaSuccess ) {
            free(u);
//            cuda_free_mem(XYcXcY_vect_flags,dev_x,dev_y,dev_dy2,dev_u);
            cudaFreeChunk(4,dev_x,dev_y,dev_dy2,dev_u);
            return ENGINE_ERROR_FAILED;
        }
    }


    // final backsubstitution loop
    for (size_t k = N-2; k > 0; --k) dy2[k] = dy2[k]*dy2[k+1] + u[k];

    free(u);
//    cuda_free_mem(XYcXcY_vect_flags,dev_x,dev_y,dev_dy2,dev_u);
    cudaFreeChunk(4,dev_x,dev_y,dev_dy2,dev_u);

    return ret_err;
}


// callable host function
// The function computes interpolated value of the tabulated function at point xx (dy2 is from previous call of spline_init function)
//__host__
//RT_engine_error spline_inter(size_t N, real_t *x, real_t *y, real_t dy2, real_t xx, real_t *yy)
//{

//}
