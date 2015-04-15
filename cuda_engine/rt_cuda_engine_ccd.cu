#include "../common-lib/num_defs.h"
#include "../common-lib/rt_engine_errors.h"
#include "../common-lib/tabulated_function.h"
#include "rt_cuda_engine_common.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdlib>
#include <cstdio>

__constant__ real_t dev_ccd_xsize;
__constant__ real_t dev_ccd_ysize;
__constant__ real_t dev_ccd_xpix;
__constant__ real_t dev_ccd_ypix;
__constant__ size_t dev_ccd_xdim;
__constant__ size_t dev_ccd_ydim;
__constant__ real_t dev_pix_area;

__constant__ real_t hc = 1.98644560233e-12; // Planck*c in ergs*mkm


                        /****************************************
                        *                                       *
                        *   nVIDIA CUDA font-end realization    *
                        *      for CCD-related functions        *
                        *                                       *
                        ****************************************/


__global__ static
void ccd_coordinates_kernel(size_t N, real_t *dev_X, real_t *dev_Y, unsigned int *dev_ccd_X, unsigned int *dev_ccd_Y)
{
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    while ( idx < N ) {
#ifdef RT_NUM_DOUBLE
        dev_ccd_X[idx] = __double2uint_rz((dev_X[idx] + dev_ccd_xsize/2.0)/dev_ccd_xpix);
        dev_ccd_Y[idx] = __double2uint_rz((dev_Y[idx] + dev_ccd_ysize/2.0)/dev_ccd_ypix);
#else
        dev_ccd_X[idx] = __float2uint_rz((dev_X[idx] + dev_ccd_xsize/2.0)/dev_ccd_xpix);
        dev_ccd_Y[idx] = __float2uint_rz((dev_Y[idx] + dev_ccd_ysize/2.0)/dev_ccd_ypix);
#endif
        idx += blockDim.x*gridDim.x;
    }
}


// main cycle
__host__ static
RT_engine_error ccd_coordinates_cycle(size_t N, size_t start_elem,
                                      real_t *X, real_t *Y, real_t *dev_X, real_t *dev_Y,
                                      unsigned int *ccd_X, unsigned int *ccd_Y, unsigned int *dev_ccd_X, unsigned int *dev_ccd_Y,
                                      int N_cuda_blocks, int N_cuda_threads)
{
    cudaError_t cuda_err;
    size_t mem_disp, vec_len;

    mem_disp = start_elem*sizeof(real_t);
    vec_len = N*sizeof(real_t);

    cuda_err = cudaMemcpy(dev_X,X+mem_disp,vec_len,cudaMemcpyHostToDevice);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    cuda_err = cudaMemcpy(dev_Y,Y+mem_disp,vec_len,cudaMemcpyHostToDevice);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }


    ccd_coordinates_kernel<<<N_cuda_blocks,N_cuda_threads>>>(N,dev_X,dev_Y,dev_ccd_X,dev_ccd_Y);

    mem_disp = start_elem*sizeof(unsigned int);
    vec_len = N*sizeof(unsigned int);

    cuda_err = cudaMemcpy(ccd_X+mem_disp,dev_ccd_X,vec_len,cudaMemcpyDeviceToHost);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    cuda_err = cudaMemcpy(ccd_Y+mem_disp,dev_ccd_Y,vec_len,cudaMemcpyDeviceToHost);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    return ENGINE_ERROR_OK;
}

// callable host function
//
// ccd_X and ccd_Y must be allocated in calling routine (vector of Nrays length)
//
// All linear sizes (CCD light and pixel sizes ) must be in millimeter!!!
//
__host__
RT_engine_error ccd_coordinates(size_t Nrays, real_t *X, real_t *Y,
                                real_t ccd_xsize, real_t ccd_ysize,
                                real_t ccd_xpix, real_t ccd_ypix,
                                unsigned int *ccd_X, unsigned int *ccd_Y)
{
    cudaError_t cuda_err;
    size_t N_dev, N_chunks, rest_Nrays,start_elem, nbytes, free_mem, total_mem;
    int N_cuda_blocks, N_cuda_threads;

    real_t *dev_X, *dev_Y;
    unsigned int *dev_ccd_X, *dev_ccd_Y;

    RT_engine_error ret_err;

    // copy CCD parameters
    cuda_err = cudaMemcpyToSymbol(dev_ccd_xsize,&ccd_xsize,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    cuda_err = cudaMemcpyToSymbol(dev_ccd_ysize,&ccd_ysize,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    cuda_err = cudaMemcpyToSymbol(dev_ccd_xpix,&ccd_xpix,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    cuda_err = cudaMemcpyToSymbol(dev_ccd_ypix,&ccd_ypix,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    // allocate device memory
    nbytes = 2*(sizeof(unsigned int) + sizeof(real_t)); // number of bytes in the device memory per two coordinate pairs

    cuda_err = cudaMemGetInfo(&free_mem,&total_mem);
    if ( cuda_err != cudaSuccess) {
        return ENGINE_ERROR_FAILED;
    }

    free_mem = free_mem*4/5; // prevent to use of entire memory

    N_dev = free_mem/nbytes;
    if ( N_dev == 0 ) { // not enough memory
        return ENGINE_ERROR_BAD_ALLOC;
    }

    if (N_dev > Nrays) N_dev = Nrays;

    cuda_err = cuda_kernel_props(N_dev,&N_cuda_blocks,&N_cuda_threads);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    cuda_err = cudaMalloc((void**)&dev_X,N_dev*sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_BAD_ALLOC;
    }

    cuda_err = cudaMalloc((void**)&dev_Y,N_dev*sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        cudaFree(dev_X);
        return ENGINE_ERROR_BAD_ALLOC;
    }

    cuda_err = cudaMalloc((void**)&dev_ccd_X,N_dev*sizeof(unsigned int));
    if ( cuda_err != cudaSuccess ) {
        cudaFree(dev_X);
        cudaFree(dev_Y);
        return ENGINE_ERROR_BAD_ALLOC;
    }

    cuda_err = cudaMalloc((void**)&dev_ccd_Y,N_dev*sizeof(unsigned int));
    if ( cuda_err != cudaSuccess ) {
        cudaFree(dev_X);
        cudaFree(dev_Y);
        cudaFree(dev_ccd_X);
        return ENGINE_ERROR_BAD_ALLOC;
    }


    N_chunks = Nrays/N_dev;

    ret_err = ENGINE_ERROR_OK;

    start_elem = 0;
    for ( size_t i = 0; i < N_chunks; ++i ) {
        ret_err = ccd_coordinates_cycle(N_dev,start_elem,X,Y,dev_X,dev_Y,
                                        ccd_X,ccd_Y,dev_ccd_X,dev_ccd_Y,
                                        N_cuda_blocks,N_cuda_threads);
        if ( ret_err != ENGINE_ERROR_OK ) break;
        start_elem += N_dev;
    }

    rest_Nrays = Nrays % N_dev;

    if ( rest_Nrays && (ret_err == ENGINE_ERROR_OK) ) {
        cuda_err = cuda_kernel_props(rest_Nrays,&N_cuda_blocks,&N_cuda_threads);
        if ( cuda_err == cudaSuccess ) {
            ret_err = ccd_coordinates_cycle(rest_Nrays,start_elem,X,Y,dev_X,dev_Y,
                                            ccd_X,ccd_Y,dev_ccd_X,dev_ccd_Y,
                                            N_cuda_blocks,N_cuda_threads);
        }
    }

    // free memory
    cudaFree(dev_X);
    cudaFree(dev_Y);
    cudaFree(dev_ccd_X);
    cudaFree(dev_ccd_Y);

    if ( cuda_err != cudaSuccess ) ret_err = ENGINE_ERROR_FAILED;

    return ret_err;
}




__global__ static
void compute_ccd_image_kernel(size_t N, unsigned int *dev_ccd_X, unsigned int *dev_ccd_Y,
                              real_t *dev_QE, real_t *dev_spec, real_t *dev_image)
{
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;
    unsigned int index;

    while ( idx < N ) {
        index = dev_ccd_Y[idx]*dev_ccd_xdim + dev_ccd_X[idx];
        dev_image[index] += dev_spec[idx]/hc*dev_QE[idx]*dev_pix_area; // here dev_spec is in ergs/s/cm^2 (see "compute_ccd_image_cycle" function)
        idx += blockDim.x*gridDim.x;
    }

}


//
// CCD image computation main cycle
//
__host__ static
RT_engine_error compute_ccd_image_cycle(size_t N, size_t start_elem,
                                        unsigned int *ccd_X, unsigned int *ccd_Y, unsigned int *dev_ccd_X, unsigned int *dev_ccd_Y,
                                        real_t *lambda, real_t *dev_QE, real_t *dev_spec, real_t *dev_image,
                                        TabulatedFunction &QE, TabulatedFunction &incident_spec,
                                        int N_cuda_blocks, int N_cuda_threads)
{
    cudaError_t cuda_err;
    real_t *host_QE, *host_spec;
    size_t vec_len = N*sizeof(real_t);


    cuda_err = cudaMemcpy(dev_ccd_X,ccd_X+start_elem,N*sizeof(unsigned int),cudaMemcpyHostToDevice);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    cuda_err = cudaMemcpy(dev_ccd_Y,ccd_Y+start_elem,N*sizeof(unsigned int),cudaMemcpyHostToDevice);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    // allocate host memory

    host_QE = (real_t*) malloc(vec_len);
    if ( host_QE == NULL ) {
        return ENGINE_ERROR_BAD_ALLOC;
    }

    host_spec = (real_t*) malloc(vec_len);
    if ( host_spec == NULL ) {
        free(host_QE);
        return ENGINE_ERROR_BAD_ALLOC;
    }


    // compute interpolated QE and spectrum at lambda-points

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
#ifdef USING_LINUX
    for (size_t i = 0; i < N; ++i ) {
#endif
#ifdef USING_MSVC // OpenMP v.2 does not support unsigned loop counter
    for (long long i = 1; i < N; ++i ) {
#endif
        size_t j = i+start_elem;
        host_QE[i] = QE[lambda[j]];
        host_spec[i] = incident_spec[lambda[j]]*lambda[j];
    }

    cuda_err = cudaMemcpy(dev_QE,host_QE,vec_len,cudaMemcpyHostToDevice);
    if ( cuda_err != cudaSuccess ) {
        free(host_QE);
        free(host_spec);
        return ENGINE_ERROR_FAILED;
    }

    cuda_err = cudaMemcpy(dev_spec,host_spec,vec_len,cudaMemcpyHostToDevice);
    if ( cuda_err != cudaSuccess ) {
        free(host_QE);
        free(host_spec);
        return ENGINE_ERROR_FAILED;
    }


    // ran kernel

    compute_ccd_image_kernel<<<N_cuda_blocks,N_cuda_threads>>>(N,dev_ccd_X,dev_ccd_Y,dev_QE,dev_spec,dev_image);

    free(host_QE);
    free(host_spec);

    return ENGINE_ERROR_OK;
}


//
// callable host function
//
// image is assumed to be a vector of ccd_xdim*ccd_ydim elements length and continuosly formed by
// ccd_ydim rows of length ccd_xdim, i.e.:
//    image[0..(ccd_xdim-1)] - the first image row,
//    image[ccd_xdim...(2*ccd_xdim-1)] - the second image row and so on ...
//
// lambda is a vector of wavelengths in mkm
//
// QE is a function of wavelength and is in E/mkm, where E is dimensionless number and is in range of [0,1]
//
// incident_spec is a function of dimension ergs/s/cm^2/mkm
//
// All linear sizes (CCD light and pixel sizes ) must be in millimeter!!!
//
//
__host__
RT_engine_error compute_ccd_image(size_t Nrays, unsigned int *ccd_X, unsigned int *ccd_Y, real_t *lambda,
                                  TabulatedFunction &QE, TabulatedFunction &incident_spec,
                                  size_t ccd_xdim, size_t ccd_ydim, real_t ccd_xpix, real_t ccd_ypix,
                                  real_t *image)
{
    cudaError_t cuda_err;
    size_t N_dev, N_chunks, rest_Nrays, start_elem, image_size, nbytes, free_mem, total_mem;
    int N_cuda_blocks, N_cuda_threads;

    real_t *dev_QE, *dev_spec, *dev_image;
    unsigned int *dev_ccd_X, *dev_ccd_Y;

    RT_engine_error ret_err;


    // copy CCD parameters
    cuda_err = cudaMemcpyToSymbol(dev_ccd_xdim,&ccd_xdim,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    cuda_err = cudaMemcpyToSymbol(dev_ccd_ydim,&ccd_ydim,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    real_t pix_area = ccd_xpix*ccd_ypix*1.0E-2; // pixel area in cm^2

    cuda_err = cudaMemcpyToSymbol(dev_pix_area,&pix_area,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    // compute working image size, allocate it and copy image to device memory

    cuda_err = cudaMemGetInfo(&free_mem,&total_mem);
    if ( cuda_err != cudaSuccess) {
        return ENGINE_ERROR_FAILED;
    }

    free_mem = free_mem*4/5; // prevent to use of entire memory

    image_size = ccd_xdim*ccd_ydim*sizeof(real_t);

    if ( image_size > free_mem ) { // not enough memory
        return ENGINE_ERROR_BAD_ALLOC;
    }
    free_mem -= image_size;

    cuda_err = cudaMalloc((void**)&dev_image,image_size);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_BAD_ALLOC;
    }

    cuda_err = cudaMemcpy(dev_image,image,image_size,cudaMemcpyHostToDevice);
    if ( cuda_err != cudaSuccess ) {
        cudaFree(dev_image);
        return ENGINE_ERROR_BAD_ALLOC;
    }


    // allocate device memory for working vectors
    nbytes = 2*(sizeof(unsigned int) + sizeof(real_t)); // number of bytes in the device memory per coordinate pair, QE and spectrum vector

    N_dev = free_mem/nbytes;
    if ( N_dev == 0 ) { // not enough memory
        return ENGINE_ERROR_BAD_ALLOC;
    }

    if (N_dev > Nrays) N_dev = Nrays;

    cuda_err = cuda_kernel_props(N_dev,&N_cuda_blocks,&N_cuda_threads);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    cuda_err = cudaMalloc((void**)&dev_QE,N_dev*sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        cudaFree(dev_image);
        return ENGINE_ERROR_BAD_ALLOC;
    }

    cuda_err = cudaMalloc((void**)&dev_spec,N_dev*sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        cudaFree(dev_image);
        cudaFree(dev_QE);
        return ENGINE_ERROR_BAD_ALLOC;
    }

    cuda_err = cudaMalloc((void**)&dev_ccd_X,N_dev*sizeof(unsigned int));
    if ( cuda_err != cudaSuccess ) {
        cudaFree(dev_image);
        cudaFree(dev_QE);
        cudaFree(dev_spec);
        return ENGINE_ERROR_BAD_ALLOC;
    }

    cuda_err = cudaMalloc((void**)&dev_ccd_Y,N_dev*sizeof(unsigned int));
    if ( cuda_err != cudaSuccess ) {
        cudaFree(dev_image);
        cudaFree(dev_QE);
        cudaFree(dev_spec);
        cudaFree(dev_ccd_X);
        return ENGINE_ERROR_BAD_ALLOC;
    }


    N_chunks = Nrays/N_dev;

    ret_err = ENGINE_ERROR_OK;

    start_elem = 0;
    for ( size_t i = 0; i < N_chunks; ++i ) {
        ret_err = compute_ccd_image_cycle(N_dev,start_elem,ccd_X,ccd_Y,dev_ccd_X,dev_ccd_Y,
                                          lambda, dev_QE, dev_spec, dev_image, QE, incident_spec,
                                          N_cuda_blocks,N_cuda_threads);
        if ( ret_err != ENGINE_ERROR_OK ) break;
        start_elem += N_dev;
    }

    rest_Nrays = Nrays % N_dev;

    if ( rest_Nrays && (ret_err == ENGINE_ERROR_OK) ) {
        cuda_err = cuda_kernel_props(rest_Nrays,&N_cuda_blocks,&N_cuda_threads);
        if ( cuda_err == cudaSuccess ) {
            ret_err = compute_ccd_image_cycle(rest_Nrays,start_elem,ccd_X,ccd_Y,dev_ccd_X,dev_ccd_Y,
                                              lambda, dev_QE, dev_spec, dev_image, QE, incident_spec,
                                              N_cuda_blocks,N_cuda_threads);
        }
    }

    // copy imge from device to host
    if ( (cuda_err == cudaSuccess) && (ret_err == ENGINE_ERROR_OK) ) {
        cuda_err = cudaMemcpy(image,dev_image,image_size,cudaMemcpyDeviceToHost);
    }

    // free memory
    cudaFree(dev_image);
    cudaFree(dev_QE);
    cudaFree(dev_spec);
    cudaFree(dev_ccd_X);
    cudaFree(dev_ccd_Y);

    if ( cuda_err != cudaSuccess ) ret_err = ENGINE_ERROR_FAILED;

    return ret_err;

}
