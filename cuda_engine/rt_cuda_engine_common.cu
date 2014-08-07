#include "../common-lib/num_defs.h"
#include "../common-lib/rt_engine_errors.h"

#include "rt_cuda_engine_common.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdarg.h>


                /***********************************
                *                                  *
                *   CUDA common helper functions   *
                *                                  *
                ***********************************/

#define CUDA_GLOBAL_MEMORY_PART 3/4 // maximal part of global memory used for computations (75%)



// the function computes optimal??? number of CUDA blocks and thread using number of elements in array as input parameter
__host__
cudaError_t cuda_kernel_props(size_t n_elems, int *N_cuda_blocks, int *N_cuda_threads)
{
    cudaError_t err;

    // get device info
    cudaDeviceProp prop;

    err = cudaGetDeviceProperties(&prop,0);
    if ( err != cudaSuccess) {
        return err;
    }

    *N_cuda_threads = prop.maxThreadsPerBlock/2;

    if ( n_elems < (prop.maxThreadsPerBlock/2) ) { // trivial case
        *N_cuda_blocks = 1;
        *N_cuda_threads = ((n_elems + prop.warpSize - 1)/prop.warpSize)*prop.warpSize; // multiple to warpSize!
    } else {
        *N_cuda_blocks = (n_elems + *N_cuda_threads - 1)/(*N_cuda_threads);
        if ( *N_cuda_blocks > prop.multiProcessorCount*2 ) *N_cuda_blocks = prop.multiProcessorCount*2;
    }

    return cudaSuccess;
}



/*
 * The function allocates memory in device global heap for a number of vectors of
 * maximal possible length.
 *
 *    inputs:
 *      required_Nelems - asked number of elements per vector
 *      N_real_t - number of vectors with elements of type "real_t" (see num_defs.h)
 *      flag_type - if true memory for vector with elements of type "beam_flag_t"
 *                  is also allocated (see num_defs.h)
 *      ... - a list of variables of type "real_t**" and, if flag_type is true, variable
 *            of type "beam_flag_t**"
 *    outputs:
 *      allocated_Nelems - number of actually allocated elements per vector. It may be
 *                         less than required_Nelems
 *
 * RESTRICTION: variable of type "beam_flag_t**" must be the last element in input list!!!
 *
 * SIDE EFFECTS: If memory allocation failed for some vector then memory for all previously
 * allocated vectors will be freed!!!
 */

__host__
cudaError_t cudaMallocChunk(size_t required_Nelems, size_t *allocated_Nelems, int N_real_t, bool flag_type, ...)
{
    cudaError_t err;
    size_t vector_size, Nbytes_per_elem, free_mem, total_mem;
    int i;

    real_t **real_t_ptr;
    beam_flag_t **beam_flag_t_ptr;

    *allocated_Nelems = 0;

    if ( required_Nelems == 0) {
        return cudaErrorInvalidValue;
    }
    if ( (N_real_t == 0) && (!flag_type) ) {
        return cudaErrorInvalidValue;
    }

    // first, compute number of bytes per sum of sizes of one elements of each vector
    if ( flag_type ) {
        Nbytes_per_elem = sizeof(beam_flag_t);
    } else {
        Nbytes_per_elem = 0;
    }

    if ( N_real_t ) Nbytes_per_elem += N_real_t*sizeof(real_t);

    // compute available memory size
    err = cudaMemGetInfo(&free_mem,&total_mem);
    if ( err != cudaSuccess) {
        return err;
    }

    free_mem = free_mem*CUDA_GLOBAL_MEMORY_PART; // prevent to use of entire memory

    *allocated_Nelems = (free_mem/Nbytes_per_elem); // number of elements in available device heap
    if ( *allocated_Nelems == 0 ) return cudaErrorMemoryAllocation; // not enough memory!!!

    if (*allocated_Nelems > required_Nelems) *allocated_Nelems = required_Nelems;

    // allocate memory
    va_list args;
    va_start(args, flag_type);

    err = cudaSuccess;

    vector_size = *allocated_Nelems*sizeof(real_t);
    for ( i = 0; i < N_real_t; ++i ) {
        real_t_ptr = va_arg(args,real_t**);
        err = cudaMalloc((void**)real_t_ptr,vector_size);
        if ( err != cudaSuccess ) { // free previous allocated memory
            *allocated_Nelems = 0;
            break;
        }
    }

    if ( i < N_real_t ) { // memory allocation failed
        va_start(args, flag_type);
        for ( int j = 0; j < (i-1); ++j ) { // free previous allocated memory
            real_t_ptr = va_arg(args,real_t**);
            cudaFree(*real_t_ptr);
        }
        va_end(args);
        return err;
    }

    if ( flag_type ) {
        beam_flag_t_ptr = va_arg(args,beam_flag_t**);
        vector_size = *allocated_Nelems*sizeof(beam_flag_t);
        err = cudaMalloc((void**)beam_flag_t_ptr,vector_size);
        if ( err != cudaSuccess ) { // free previous allocated memory
            *allocated_Nelems = 0;
            va_start(args, flag_type);
            for ( i = 0; i < N_real_t; ++i) {
                real_t_ptr = va_arg(args,real_t**);
                cudaFree(*real_t_ptr);
            }
        }
    }

    va_end(args);

    return err;
}

/*
* The function frees memory for a number of vectors.
*
* inputs:
*
*
*/
__host__
//void cudaFreeChunk(int N_real_t, bool flag_type, ...)
void cudaFreeChunk(int N, ...)
{
//    real_t *real_t_ptr;
//    beam_flag_t * beam_flag_t_ptr;

//    if ( (N_real_t == 0) && (!flag_type) ) return;

    void *ptr;

    if ( !N ) return;

    va_list args;
//    va_start(args, flag_type);
    va_start(args, N);

//    for ( int i = 0; i < N_real_t; ++i ) {
//        real_t_ptr = va_arg(args,real_t*);
//        cudaFree(real_t_ptr);
//    }

//    if ( flag_type ) {
//        beam_flag_t_ptr = va_arg(args,beam_flag_t*);
//        cudaFree(real_t_ptr);
//    }

    for ( int i = 0; i < N; ++i ) {
        ptr = va_arg(args,void*);
        cudaFree(ptr);
    }

    va_end(args);
 }


/*
 *  The function copies memory from host to device or vice-versa for a number of vectors.
 *
 *  inputs:
 *    Nelems - number of elements to be copied in each input vectors
 *    kind - type of transfer (cudaMemcpyHostToDevice or cudaMemcpyDeviceToHost)
 *           see cudaMemcpy function for details
 *    N_real_t - number pairs of vectors of type "real_t*"
 *    flag_type - if it is true then the last pair of the input
 *                vectors is interpretated as of type "beam_flag_t*"
 *    ... - list of variable pairs of type "real_t*" and/or "beam_flag_t*":
 *          the first in the pair is destination vector, the second - source vector
 *          (just like in the CUDA function cudaMemcpy)
 *
 * RESTRICTION: pair of variables of type "beam_flag_t*" must be the last two elements in input list!!!
 *
 * USAGE:
 *   copy memory from two host vectors X, Y and flag to device vectors dev_X, dev_Y and dev_flag:
 *   	err = cudamemcpyChunk(100,cudaMemcpyHostToDevice,2,true,dev_X,X,dev_Y,Y,dev_flag,flag)
 *
 *   copy memory from two device vectors X and Y to host vectors dev_X and dev_Y:
 *   	err = cudamemcpyChunk(100,cudaMemcpyDeviceToHost,2,false,X,dev_X,Y,dev_Y)
 */

__host__
cudaError_t cudaMemcpyChunk(size_t Nelems, cudaMemcpyKind kind,
                            int N_real_t, bool flag_type, ... )
{
    cudaError_t err;
    real_t *real_t_dst, *real_t_src;
    beam_flag_t *beam_flag_t_dst, *beam_flag_t_src;
    size_t vector_size = Nelems*sizeof(real_t);

    if ( Nelems == 0 ) return cudaErrorInvalidValue;
    if ( (N_real_t == 0) && (!flag_type) ) return cudaErrorInvalidValue;

    va_list args;
    va_start(args, flag_type);

    err = cudaSuccess;

    for ( int i = 0; i < N_real_t; ++i ) {
        real_t_dst = va_arg(args,real_t*);
        real_t_src = va_arg(args,real_t*);
        err = cudaMemcpy(real_t_dst,real_t_src,vector_size,kind);
        if ( err != cudaSuccess ) {
            break;
        }
    }

    if ( flag_type && (err == cudaSuccess) ) {
        vector_size = Nelems*sizeof(beam_flag_t);
        beam_flag_t_dst = va_arg(args,beam_flag_t*);
        beam_flag_t_src = va_arg(args,beam_flag_t*);
        err = cudaMemcpy(beam_flag_t_dst,beam_flag_t_src,vector_size,kind);
    }

    va_end(args);

    return err;
}

/*  CUDA engine host memory allocation and freeing functions  */


__host__
void rt_engine_free_vector(real_t* data)
{
    free(data);
    data = NULL;
}

__host__
real_t* rt_engine_alocate_vector(size_t N_elem)
{
    real_t *ptr;

    ptr = (real_t*) malloc(sizeof(real_t)*N_elem);
    return ptr;
}

