#include<cstdlib>
#include<cstdio>

//#if defined(WIN32) && defined(_MSC_VER)
//    #define USING_MSVC
//#endif

//#if defined(__GNUC__) && !defined(USING_MACOSX)
//    #define USING_LINUX
//#endif


#include<cuda.h>
#include<curand.h>

#include "../common-lib/num_defs.h"
#include "../common-lib/rt_engine_errors.h"
#include "../common-lib/beam.h"

#include "rt_cuda_engine_common.h"


// headers for time stamps
#ifdef USING_LINUX // GCC
#include<time.h> // for POSIX systems
#endif

#ifdef USING_MSVC // Visual Studio
#include<Windows.h>
#include<cstdint>
#endif


                        /****************************************
                        *                                       *
                        *   nVIDIA CUDA font-end realization    *
                        *       for Beam class functions        *
                        *                                       *
                        ****************************************/

__constant__ real_t PI2 = 6.28318531; // 2*Pi


// beam parameters (see beam.cpp for details)
__constant__ real_t Par0;
__constant__ real_t Par1;
__constant__ real_t Par2;
__constant__ real_t Par3;

__constant__ real_t cenX;
__constant__ real_t cenY;
__constant__ real_t cenZ;

// translations and rotation sin and cosin
__constant__ real_t dev_dist[3];
__constant__ real_t sinA;
__constant__ real_t cosA;


// start lambda and its step (for uniform wavelength range distribution)
// for random distribution in lambda_step the maximal wavelength is stored
__constant__ real_t start_lambda;
__constant__ real_t lambda_step;

        /*  Auxialiary functions  */


// the function copies beam parameters and center coordinates to device constant memory
__host__
static RT_engine_error copy_beam_params(Beam::BeamShape shape, real_t *params, real_t* center)
{

    cudaError_t cuda_err;

    // common part of parameters
    cuda_err = cudaMemcpyToSymbol(Par0,params,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) return ENGINE_ERROR_FAILED;

    cuda_err = cudaMemcpyToSymbol(Par1,params+1,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) return ENGINE_ERROR_FAILED;

    cuda_err = cudaMemcpyToSymbol(Par2,params+2,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) return ENGINE_ERROR_FAILED;

    switch (shape) {
        case Beam::Circle: break;
        case Beam::Rectangle: {
            cuda_err = cudaMemcpyToSymbol(Par3,params+3,sizeof(real_t));
            if ( cuda_err != cudaSuccess ) return ENGINE_ERROR_FAILED;
            break;
        }
        default: return ENGINE_ERROR_FAILED;
    }

    cuda_err = cudaMemcpyToSymbol(cenX,center,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) return ENGINE_ERROR_FAILED;

    cuda_err = cudaMemcpyToSymbol(cenY,center+1,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) return ENGINE_ERROR_FAILED;

    cuda_err = cudaMemcpyToSymbol(cenZ,center+2,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) return ENGINE_ERROR_FAILED;

    return ENGINE_ERROR_OK;
}


         /*-----------------------------*
         *  Beam of random coordinates  *
         *          functions           *
         *-----------------------------*/


//__global__
//static void ss(size_t N, real_t* dev_X, real_t* dev_Y, real_t* dev_Z)
//{
//    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

//    while ( idx < N ) {
//        printf("--XYZ: %f %f %f\n",dev_X[idx],dev_Y[idx],dev_Z[idx]);
//        idx += blockDim.x*gridDim.x;
//    }
//}

// random coordinate beam kernel functions
__global__
static void random_circle_kernel(size_t N, real_t* dev_X, real_t* dev_Y, real_t* dev_Z)
{
    real_t R,phi;

    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    // random numbers are now in dev_X and dev_Y vectors
    // Par1 - inner circle radius
    // Par2 - outter circle radius

    while ( idx < N ) {
        phi = PI2*dev_Y[idx]; // polar angle
        dev_Z[idx] = cenZ;
#ifdef RT_NUM_DOUBLE
        R = (Par2-Par1)*sqrt(dev_X[idx]) + Par1;   // radius-vector
        dev_X[idx] = R*cos(phi) + cenX;
        dev_Y[idx] = R*sin(phi) + cenY;
#else
        R = (Par2-Par1)*sqrtf(dev_X[idx]) + Par1;  // radius-vector
        dev_X[idx] = R*cosf(phi) + cenX;
        dev_Y[idx] = R*sinf(phi) + cenY;
#endif
//        printf("XYZ: %f %f %f, CEN: %f %f %f\n",dev_X[idx],dev_Y[idx],dev_Z[idx],cenX,cenY,cenZ);
        idx += blockDim.x*gridDim.x;
    }
}

__global__ void random_rectangle_kernel(size_t N, real_t* dev_X, real_t* dev_Y, real_t* dev_Z)
{
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    dev_Z[idx] = cenZ;

    while ( idx < N ) {
        idx += blockDim.x*gridDim.x;
    }
}

// random coordinate beam main cycle
__host__ RT_engine_error beam_random_cycle(Beam::BeamShape shape, curandGenerator_t& gen,
                                           size_t Nrays, size_t start_elem,
                                           real_t* X, real_t* Y, real_t* Z,
                                           real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                                           int N_cuda_blocks, int N_cuda_threads)
{
    curandStatus_t curand_err;
    cudaError_t cuda_err;

#ifdef RT_NUM_DOUBLE
    curand_err = curandGenerateUniformDouble(gen,dev_X,Nrays);
#else
    curand_err = curandGenerateUniform(gen,dev_X,Nrays);
#endif
    if ( curand_err != CURAND_STATUS_SUCCESS ) {
        return ENGINE_ERROR_FAILED;
    }
#ifdef RT_NUM_DOUBLE
    curand_err = curandGenerateUniformDouble(gen,dev_Y,Nrays);
#else
    curand_err = curandGenerateUniform(gen,dev_Y,Nrays);
#endif
    if ( curand_err != CURAND_STATUS_SUCCESS ) {
        return ENGINE_ERROR_FAILED;
    }

//    ss<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_X,dev_Y,dev_Z);

    switch ( shape ) {
        case Beam::Circle: {
            random_circle_kernel<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_X,dev_Y,dev_Z);
            break;
        }
        case Beam::Rectangle: {
            random_rectangle_kernel<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_X,dev_Y,dev_Z);
            break;
        }
        default: break;
    }

    // copy results from device to host memory
//    cuda_err = cuda_copy_mem(Nrays,start_elem,cudaMemcpyDeviceToHost,XYZ_vect_flags,X,dev_X,Y,dev_Y,Z,dev_Z);
    cuda_err = cudaMemcpyChunk(Nrays,cudaMemcpyDeviceToHost,3,false,
                               X+start_elem,dev_X,Y+start_elem,dev_Y,Z+start_elem,dev_Z);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    return ENGINE_ERROR_OK;

}


// callable function
// NOTE: The CUDA function curandGenerateUniform returns random numbers in range (0,1] thus
//       it does not guarantee that points at exact maximal radius are returned in final vectors and
//       points at exact minimal radius are never returned also!!!

__host__ RT_engine_error beam_random(Beam::BeamShape shape, real_t* params, real_t* center,
                                     size_t Nrays, real_t* X, real_t* Y, real_t* Z)
{
    cudaError_t cuda_err;
    curandStatus_t curand_err;
    RT_engine_error ret_err;
    int N_cuda_blocks, N_cuda_threads;
    curandGenerator_t gen;

    size_t dev_Nrays, N_chunks, rest_Nrays, start_elem;

    real_t *dev_X, *dev_Y,  *dev_Z;

//    printf("\n\n initial addr: %p\n",dev_X);

    cuda_err = cudaMallocChunk(Nrays,&dev_Nrays,3,false,&dev_X,&dev_Y,&dev_Z);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_BAD_ALLOC;
    }

//    printf("returned addr: %p\n",dev_X);

    cuda_err = cuda_kernel_props(dev_Nrays,&N_cuda_blocks,&N_cuda_threads);
    if ( cuda_err != cudaSuccess ) {
        cudaFreeChunk(3,dev_X,dev_Y,dev_Z);
        return ENGINE_ERROR_FAILED;
    }

    N_chunks = Nrays/dev_Nrays;

    ret_err = ENGINE_ERROR_OK;

    ret_err = copy_beam_params(shape,params,center);
    if ( ret_err != ENGINE_ERROR_OK ) {
        cudaFreeChunk(3,dev_X,dev_Y,dev_Z);
        return ret_err;
    }

    // create random generator and init random seed
    curand_err = curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
    if ( curand_err != CURAND_STATUS_SUCCESS ) {
        cudaFreeChunk(3,dev_X,dev_Y,dev_Z);
        return ENGINE_ERROR_FAILED;
    }

    // generate seed
#ifdef USING_LINUX
    u_int64_t seed;
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC_RAW,&ts);
    seed = (u_int64_t)ts.tv_nsec;
#endif

#ifdef USING_MSVC
    uint64_t seed;
    LARGE_INTEGER ts;
    QueryPerformanceCounter(&ts);
    seed = (uint64_t)ts.QuadPart;
#endif

    curand_err = curandSetPseudoRandomGeneratorSeed(gen, seed);
    if ( curand_err != CURAND_STATUS_SUCCESS ) {
        cudaFreeChunk(3,dev_X,dev_Y,dev_Z);
        return ENGINE_ERROR_FAILED;
    }


    // main cycle
    start_elem = 0;
    for ( size_t i = 0; i < N_chunks; ++i) {
        ret_err = beam_random_cycle(shape,gen,dev_Nrays,start_elem,X,Y,Z,dev_X,dev_Y,dev_Z,N_cuda_blocks,N_cuda_threads);
        if ( ret_err != ENGINE_ERROR_OK ) break;
        start_elem += dev_Nrays;
    }

    // process rest of rays
    rest_Nrays = Nrays % dev_Nrays;
    if ( rest_Nrays && (ret_err == ENGINE_ERROR_OK) ) {
        cuda_err = cuda_kernel_props(rest_Nrays,&N_cuda_blocks,&N_cuda_threads);
        if ( cuda_err != cudaSuccess ) {
            cudaFreeChunk(3,dev_X,dev_Y,dev_Z);
            return ENGINE_ERROR_FAILED;
        }
        ret_err = beam_random_cycle(shape,gen,rest_Nrays,start_elem,X,Y,Z,dev_X,dev_Y,dev_Z,N_cuda_blocks,N_cuda_threads);
    }

    cudaFreeChunk(3,dev_X,dev_Y,dev_Z);
    return ret_err;
}


            /*---------------------------------------*
             *  Translation and rotation functions   *
             *---------------------------------------*/


__global__
static void translation_kernel(size_t N, real_t *dev_X, real_t *dev_Y, real_t *dev_Z)
{
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    while ( idx < N ) {
        dev_X[idx] -= dev_dist[0];
        dev_Y[idx] -= dev_dist[1];
        dev_Z[idx] -= dev_dist[2];
        idx += blockDim.x*gridDim.x;
    }
}

__global__
static void rotation_kernel(size_t N, real_t *v1, real_t *v2, real_t *v3, real_t *v4)
{
    real_t e;
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    while ( idx < N ) {
        e = v1[idx]*cosA - v2[idx]*sinA;
        v2[idx] = v1[idx]*sinA + v2[idx]*cosA;
        v1[idx] = e;

        e = v3[idx]*cosA - v4[idx]*sinA;
        v4[idx] = v3[idx]*sinA + v4[idx]*cosA;
        v3[idx] = e;

        idx += blockDim.x*gridDim.x;
    }
}

// translation main cycle
__host__
static RT_engine_error beam_translate_cycle(size_t Nrays, size_t start_elem, real_t *X, real_t *Y, real_t *Z,
                                            real_t *dev_X, real_t *dev_Y, real_t *dev_Z,
                                            int N_cuda_blocks, int N_cuda_threads)
{
    cudaError_t cuda_err;

//    cuda_err = cuda_copy_mem(Nrays,start_elem,cudaMemcpyHostToDevice,XYZ_vect_flags,
//                             X,dev_X,Y,dev_Y,Z,dev_Z);
    cuda_err = cudaMemcpyChunk(Nrays,cudaMemcpyHostToDevice,3,false,
                               dev_X,X+start_elem,dev_Y,Y+start_elem,dev_Z,Z+start_elem);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    translation_kernel<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_X,dev_Y,dev_Z);

//    cuda_err = cuda_copy_mem(Nrays,start_elem,cudaMemcpyDeviceToHost,XYZ_vect_flags,
//                             X,dev_X,Y,dev_Y,Z,dev_Z);
    cuda_err = cudaMemcpyChunk(Nrays,cudaMemcpyDeviceToHost,3,false,
                               X+start_elem,dev_X,Y+start_elem,dev_Y,Z+start_elem,dev_Z);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    return ENGINE_ERROR_OK;
}


// rotation main cycle
__host__
static RT_engine_error beam_rotate_cycle(size_t Nrays, size_t start_elem,
                                         real_t *v1, real_t *v2, real_t *cv1, real_t *cv2,
                                         real_t *dev_v1, real_t *dev_v2,
                                         real_t *dev_cv1, real_t *dev_cv2,
                                         int N_cuda_blocks, int N_cuda_threads)
{
    cudaError_t cuda_err;

//    cuda_err = cuda_copy_mem(Nrays,start_elem,cudaMemcpyHostToDevice,bit_flags,
//                             v1,dev_v1,v2,dev_v2,cv1,dev_cv1,cv2,dev_cv2);
    cuda_err = cudaMemcpyChunk(Nrays,cudaMemcpyHostToDevice,4,false,
                               dev_v1,v1+start_elem,dev_v2,v2+start_elem,
                               dev_cv1,cv1+start_elem,dev_cv2,cv2+start_elem);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    rotation_kernel<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_v1,dev_v2,dev_cv1,dev_cv2);

//    cuda_err = cuda_copy_mem(Nrays,start_elem,cudaMemcpyDeviceToHost,bit_flags,
//                             v1,dev_v1,v2,dev_v2,cv1,dev_cv1,cv2,dev_cv2);
    cuda_err = cudaMemcpyChunk(Nrays,cudaMemcpyDeviceToHost,4,false,
                               v1+start_elem,dev_v1,v2+start_elem,dev_v2,
                               cv1+start_elem,dev_cv1,cv2+start_elem,dev_cv2);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    return ENGINE_ERROR_OK;
}

// callable driver function
__host__
RT_engine_error beam_transform(size_t Nrays, real_t *angles, real_t *dist, real_t *X, real_t *Y, real_t *Z, real_t *cX, real_t *cY, real_t *cZ)
{
    cudaError_t cuda_err;
    RT_engine_error ret_err;
    int N_cuda_blocks, N_cuda_threads;

    size_t dev_Nrays, N_chunks, rest_Nrays, start_elem;

    real_t *dev_X, *dev_Y, *dev_Z;

//    printf("COO BEFORE: %f %f %f\n",X[10],Y[10],Z[10]);

    // copy translation values
    cuda_err = cudaMemcpyToSymbol(dev_dist,dist,3*sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    // first, do translation
    if ( (dist[2] != 0.0) || (dist[0] != 0.0) || (dist[1] != 0.0) ) {
//        cuda_err = cuda_malloc_vectors(Nrays,&dev_Nrays,XYZ_vect_flags,&dev_X,&dev_Y,&dev_Z);
        cuda_err = cudaMallocChunk(Nrays,&dev_Nrays,3,false,&dev_X,&dev_Y,&dev_Z);
        if ( cuda_err != cudaSuccess ) {
            return ENGINE_ERROR_BAD_ALLOC;
        }

        cuda_err = cuda_kernel_props(dev_Nrays,&N_cuda_blocks,&N_cuda_threads);
        if ( cuda_err != cudaSuccess ) {
//            cuda_free_mem(XYZ_vect_flags,dev_X,dev_Y,dev_Z);
            cudaFreeChunk(3,dev_X,dev_Y,dev_Z);
            return ENGINE_ERROR_FAILED;
        }

        N_chunks = Nrays/dev_Nrays;

        ret_err = ENGINE_ERROR_OK;

        start_elem = 0;
        for ( size_t i = 0; i < N_chunks; ++i) {
            ret_err = beam_translate_cycle(dev_Nrays,start_elem,X,Y,Z,dev_X,dev_Y,dev_Z,N_cuda_blocks,N_cuda_threads);
            if ( ret_err != ENGINE_ERROR_OK ) break;
            start_elem += dev_Nrays;
        }

        // process rest of rays
        rest_Nrays = Nrays % dev_Nrays;
        if ( rest_Nrays && (ret_err == ENGINE_ERROR_OK) ) {
            cuda_err = cuda_kernel_props(rest_Nrays,&N_cuda_blocks,&N_cuda_threads);
            if ( cuda_err != cudaSuccess ) {
                cudaFreeChunk(3,dev_X,dev_Y,dev_Z);
                return ENGINE_ERROR_FAILED;
            }
            ret_err = beam_translate_cycle(rest_Nrays,start_elem,X,Y,Z,dev_X,dev_Y,dev_Z,N_cuda_blocks,N_cuda_threads);
        }

        cudaFreeChunk(3,dev_X,dev_Y,dev_Z);
    }
    cuda_err = cudaThreadSynchronize();
    if ( cuda_err != cudaSuccess ) return ENGINE_ERROR_FAILED;

//    printf("COO AFTER: %f %f %f\n",X[10],Y[10],Z[10]);

    // rotate beam

//    bit_flag_t bf;
    real_t *v1, *v2, *dev_v1, *dev_v2, *cv1, *cv2, *dev_cv1, *dev_cv2;
    real_t sin_ang, cos_ang;

    // 0 - rotation about Y-axis
    // 1 - about X-axis
    // 2 - about Z-axis
    for ( size_t ia = 0; ia < 3; ++ia ) {
        switch ( ia ) {
            case 0: {
                v1 = X;
                v2 = Z;
                cv1 = cX;
                cv2 = cZ;
//                bf = XZcXcZ_vect_flags;
                break;
            }
            case 1: {
                v1 = Y;
                v2 = Z;
                cv1 = cY;
                cv2 = cZ;
//                bf = YZcYcZ_vect_flags;
                break;
            }
            case 2: {
                v1 = X;
                v2 = Y;
                cv1 = cX;
                cv2 = cY;
//                bf = XYcXcY_vect_flags;
                break;
            }
            default: break;
        }

        if ( angles[ia] != 0 ) {
#ifdef RT_NUM_DOUBLE
            sin_ang = sin(angles[ia]);
            cos_ang = cos(angles[ia]);
#else
            sin_ang = sinf(angles[ia]);
            cos_ang = cosf(angles[ia]);
#endif
            // copy sin and cosin to tdevice memory
            cuda_err = cudaMemcpyToSymbol(sinA,&sin_ang,sizeof(real_t));
            if ( cuda_err != cudaSuccess ) {
                return ENGINE_ERROR_BAD_ALLOC;
            }

            cuda_err = cudaMemcpyToSymbol(cosA,&cos_ang,sizeof(real_t));
            if ( cuda_err != cudaSuccess ) {
                return ENGINE_ERROR_BAD_ALLOC;
            }

//            cuda_err = cuda_malloc_vectors(Nrays,&dev_Nrays,bf,&dev_v1,&dev_v2,&dev_cv1,&dev_cv2);
            cuda_err = cudaMallocChunk(Nrays,&dev_Nrays,4,false,&dev_v1,&dev_v2,&dev_cv1,&dev_cv2);
            if ( cuda_err != cudaSuccess ) {
                return ENGINE_ERROR_BAD_ALLOC;
            }

            cuda_err = cuda_kernel_props(dev_Nrays,&N_cuda_blocks,&N_cuda_threads);
            if ( cuda_err != cudaSuccess ) {
//                cuda_free_mem(bf,dev_v1,dev_v2,dev_cv1,dev_cv2);
                cudaFreeChunk(4,dev_v1,dev_v2,dev_cv1,dev_cv2);
                return ENGINE_ERROR_FAILED;
            }

            N_chunks = Nrays/dev_Nrays;

            ret_err = ENGINE_ERROR_OK;

            start_elem = 0;
            for ( size_t i = 0; i < N_chunks; ++i) {
                ret_err = beam_rotate_cycle(dev_Nrays,start_elem,v1,v2,cv1,cv2,dev_v1,dev_v2,dev_cv1,dev_cv2,N_cuda_blocks,N_cuda_threads);
                if ( ret_err != ENGINE_ERROR_OK ) break;
                start_elem += dev_Nrays;
            }

            // process rest of rays
            rest_Nrays = Nrays % dev_Nrays;
            if ( rest_Nrays && (ret_err == ENGINE_ERROR_OK) ) {
                cuda_err = cuda_kernel_props(rest_Nrays,&N_cuda_blocks,&N_cuda_threads);
                if ( cuda_err != cudaSuccess ) {
//                    cuda_free_mem(bf,dev_v1,dev_v2,dev_cv1,dev_cv2);
                    cudaFreeChunk(4,dev_v1,dev_v2,dev_cv1,dev_cv2);
                    return ENGINE_ERROR_FAILED;
                }
                ret_err = beam_rotate_cycle(rest_Nrays,start_elem,v1,v2,cv1,cv2,dev_v1,dev_v2,dev_cv1,dev_cv2,N_cuda_blocks,N_cuda_threads);
            }

//            cuda_free_mem(bf,dev_v1,dev_v2,dev_cv1,dev_cv2);
            cudaFreeChunk(4,dev_v1,dev_v2,dev_cv1,dev_cv2);
        }
    }

    return ENGINE_ERROR_OK;
}


                /*----------------------------*
                 *  Beam cosins calculations  *
                 *----------------------------*/

__global__
static void conic_cosins_kernel(size_t N, real_t conic_angle, real_t *dev_X, real_t *dev_Y, real_t *dev_cX, real_t *dev_cY, real_t *dev_cZ)
{
    real_t sin_con_ang, cos_con_ang;
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

#ifdef RT_NUM_DOUBLE
    cos_con_ang = cos(conic_angle/2.0);
    sin_con_ang = sin(conic_angle/2.0);
#else
    cos_con_ang = cosf(conic_angle/2.0);
    sin_con_ang = sinf(conic_angle/2.0);
#endif
    while ( idx < N ) {
        dev_cZ[idx] = cos_con_ang;
        // NOTE: CUDA atan2(P1,P2) function computes arctan of P1/P2!!!
#ifdef RT_NUM_DOUBLE
        dev_cX[idx] = sin_con_ang*cos(atan2(dev_Y[idx],dev_X[idx]));
        dev_cY[idx] = sin_con_ang*sin(atan2(dev_Y[idx],dev_X[idx]));
#else
        dev_cX[idx] = sin_con_ang*cosf(atan2f(dev_Y[idx],dev_X[idx]));
        dev_cY[idx] = sin_con_ang*sinf(atan2f(dev_Y[idx],dev_X[idx]));
#endif
        idx += blockDim.x*gridDim.x;
    }
}


__host__
static RT_engine_error conic_cosins_cycle(size_t Nrays,size_t start_elem,real_t conic_angle, real_t *X, real_t *Y, real_t *cX, real_t *cY, real_t *cZ,
                                          real_t *dev_X, real_t *dev_Y, real_t *dev_cX, real_t *dev_cY, real_t *dev_cZ,
                                          int N_cuda_blocks, int N_cuda_threads)
{
    cudaError_t cuda_err;


    conic_cosins_kernel<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,conic_angle,dev_X,dev_Y,dev_cX,dev_cY,dev_cZ);

//    cuda_err = cuda_copy_mem(Nrays,start_elem,cudaMemcpyDeviceToHost,XYcXcYcZ_vect_flags,
//                             X,dev_X,Y,dev_Y,cX,dev_cX,cY,dev_cY,cZ,dev_cZ);
    cuda_err = cudaMemcpyChunk(Nrays,cudaMemcpyDeviceToHost,5,false,
                               X+start_elem,dev_X,Y+start_elem,dev_Y,
                               cX+start_elem,dev_cX,cY+start_elem,dev_cY,
                               cZ+start_elem,dev_cZ);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    return ENGINE_ERROR_OK;

}


// callable host function (conic_angle must be in radians)
__host__
RT_engine_error beam_conic_cosins(size_t Nrays, real_t conic_angle, real_t *X, real_t *Y, real_t *cX, real_t *cY, real_t *cZ)
{
    cudaError_t cuda_err;
    RT_engine_error ret_err;
    int N_cuda_blocks, N_cuda_threads;

    size_t dev_Nrays, N_chunks, rest_Nrays, start_elem;

    real_t *dev_X, *dev_Y, *dev_cX, *dev_cY, *dev_cZ;

//    cuda_err = cuda_malloc_vectors(Nrays,&dev_Nrays,XYcXcYcZ_vect_flags,
//                                   &dev_X,&dev_Y,&dev_cX,&dev_cY,&dev_cZ);
    cuda_err = cudaMallocChunk(Nrays,&dev_Nrays,5,false,
                               &dev_X,&dev_Y,&dev_cX,&dev_cY,&dev_cZ);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_BAD_ALLOC;
    }

    cuda_err = cuda_kernel_props(dev_Nrays,&N_cuda_blocks,&N_cuda_threads);
    if ( cuda_err != cudaSuccess ) {
//        cuda_free_mem(XYcXcYcZ_vect_flags,dev_X,dev_Y,dev_cX,dev_cY,dev_cZ);
        cudaFreeChunk(5,dev_X,dev_Y,dev_cX,dev_cY,dev_cZ);
        return ENGINE_ERROR_FAILED;
    }

    N_chunks = Nrays/dev_Nrays;

    ret_err = ENGINE_ERROR_OK;

    // main cycle
    start_elem = 0;
    for ( size_t i = 0; i < N_chunks; ++i) {
        ret_err = conic_cosins_cycle(dev_Nrays,start_elem,conic_angle,X,Y,cX,cY,cZ,dev_X,dev_Y,dev_cX,dev_cY,dev_cZ,N_cuda_blocks,N_cuda_threads);
        if ( ret_err != ENGINE_ERROR_OK ) break;
        start_elem += dev_Nrays;
    }

    // process rest of rays
    rest_Nrays = Nrays % dev_Nrays;
    if ( rest_Nrays && (ret_err == ENGINE_ERROR_OK) ) {
        cuda_err = cuda_kernel_props(rest_Nrays,&N_cuda_blocks,&N_cuda_threads);
        if ( cuda_err != cudaSuccess ) {
//            cuda_free_mem(XYcXcYcZ_vect_flags,dev_X,dev_Y,dev_cX,dev_cY,dev_cZ);
            cudaFreeChunk(5,dev_X,dev_Y,dev_cX,dev_cY,dev_cZ);
            return ENGINE_ERROR_FAILED;
        }
        ret_err = conic_cosins_cycle(rest_Nrays,start_elem,conic_angle,X,Y,cX,cY,cZ,dev_X,dev_Y,dev_cX,dev_cY,dev_cZ,N_cuda_blocks,N_cuda_threads);
    }

//    cuda_free_mem(XYcXcYcZ_vect_flags,dev_X,dev_Y,dev_cX,dev_cY,dev_cZ);
    cudaFreeChunk(5,dev_X,dev_Y,dev_cX,dev_cY,dev_cZ);
    return ret_err;
}



                    /*-------------------------------*
                     *  Beam wavelength calculation  *
                     *-------------------------------*/

__global__
static void uniform_range_beam_kernel(size_t N, real_t *dev_lambda)
{
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    while ( idx < N ) {
        dev_lambda[idx] = start_lambda + idx*lambda_step;
        idx += blockDim.x*gridDim.x;
    }
}


// callable host function
__host__
RT_engine_error uniform_range_beam(size_t Nrays, real_t *range, real_t *lambda)
{
    cudaError_t cuda_err;
    RT_engine_error ret_err;
    int N_cuda_blocks, N_cuda_threads;

    size_t dev_Nrays, N_chunks, rest_Nrays, start_elem;

    real_t *dev_lambda, d_l;

    // compute step for lambda vector
    d_l = (range[1]-range[0])/(Nrays-1);

    // copy start lambda and step values
    cuda_err = cudaMemcpyToSymbol(start_lambda,range,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    cuda_err = cudaMemcpyToSymbol(lambda_step,&d_l,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

//    cuda_err = cuda_malloc_vectors(Nrays,&dev_Nrays,lambda_vect_flag,&dev_lambda);
    cuda_err = cudaMallocChunk(Nrays,&dev_Nrays,1,false,&dev_lambda);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_BAD_ALLOC;
    }

    cuda_err = cuda_kernel_props(dev_Nrays,&N_cuda_blocks,&N_cuda_threads);
    if ( cuda_err != cudaSuccess ) {
//        cuda_free_mem(lambda_vect_flag,dev_lambda);
        cudaFreeChunk(1,dev_lambda);
        return ENGINE_ERROR_FAILED;
    }

    N_chunks = Nrays/dev_Nrays;

    ret_err = ENGINE_ERROR_OK;

    start_elem = 0;
    for ( size_t i = 0; i < N_chunks; ++i) {
        uniform_range_beam_kernel<<<N_cuda_blocks,N_cuda_threads>>>(dev_Nrays,dev_lambda);

//        cuda_err = cuda_copy_mem(Nrays,start_elem,cudaMemcpyDeviceToHost,
//                                 lambda_vect_flag,lambda,dev_lambda);
        cuda_err = cudaMemcpyChunk(dev_Nrays,cudaMemcpyDeviceToHost,1,false,
                                   lambda+start_elem,dev_lambda);
        if ( cuda_err != cudaSuccess ) {
            break;
        }

        start_elem += dev_Nrays;
    }

    // process rest of rays
    rest_Nrays = Nrays % dev_Nrays;
    if ( rest_Nrays && (ret_err == ENGINE_ERROR_OK) ) {
        cuda_err = cuda_kernel_props(rest_Nrays,&N_cuda_blocks,&N_cuda_threads);
        if ( cuda_err != cudaSuccess ) {
//            cuda_free_mem(lambda_vect_flag,dev_lambda);
            cudaFreeChunk(1,dev_lambda);
            return ENGINE_ERROR_FAILED;
        }

        uniform_range_beam_kernel<<<N_cuda_blocks,N_cuda_threads>>>(dev_Nrays,dev_lambda);

//        cuda_err = cuda_copy_mem(Nrays,start_elem,cudaMemcpyDeviceToHost,lambda_vect_flag,lambda,dev_lambda);
        cuda_err = cudaMemcpyChunk(rest_Nrays,cudaMemcpyDeviceToHost,1,false,
                                   lambda+start_elem,dev_lambda);
        if ( cuda_err != cudaSuccess ) {
//            cuda_free_mem(lambda_vect_flag,dev_lambda);
            cudaFreeChunk(1,dev_lambda);
            return ENGINE_ERROR_FAILED;
        }
    }

//    cuda_free_mem(lambda_vect_flag,dev_lambda);
    cudaFreeChunk(1,dev_lambda);

    return ret_err;
}



__global__
static void random_range_beam_kernel(size_t N, real_t *dev_lambda)
{
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    // random numbers in (0,1] are in dev_lambda
    while ( idx < N ) {
        dev_lambda[idx] = (lambda_step - start_lambda)*dev_lambda[idx] + start_lambda;
        idx += blockDim.x*gridDim.x;
    }

}

// callable host function
// NOTE: The CUDA function curandGenerateUniform returns random numbers in range (0,1] thus
//       it does not guarantee exact maximal wavelength in the returned vector and
//       exact minimal wavelength is never returned also!!!
__host__
RT_engine_error random_range_beam(size_t Nrays, real_t *range, real_t *lambda)
{
    cudaError_t cuda_err;
    curandStatus_t curand_err;
    RT_engine_error ret_err;
    int N_cuda_blocks, N_cuda_threads;
    curandGenerator_t gen;

    size_t dev_Nrays, N_chunks, rest_Nrays, start_elem;

    real_t *dev_lambda;

    // copy minimal and maximal lambda values
    cuda_err = cudaMemcpyToSymbol(start_lambda,range,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    cuda_err = cudaMemcpyToSymbol(lambda_step,range+1,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    // create random generator and init random seed
    curand_err = curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
    if ( curand_err != CURAND_STATUS_SUCCESS ) {
        return ENGINE_ERROR_FAILED;
    }

    // generate seed
#ifdef USING_LINUX
    u_int64_t seed;
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC_RAW,&ts);
    seed = (u_int64_t)ts.tv_nsec;
#endif

#ifdef USING_MSVC
    uint64_t seed;
    LARGE_INTEGER ts;
    QueryPerformanceCounter(&ts);
    seed = (uint64_t)ts.QuadPart;
#endif

    curand_err = curandSetPseudoRandomGeneratorSeed(gen, seed);
    if ( curand_err != CURAND_STATUS_SUCCESS ) {
        return ENGINE_ERROR_FAILED;
    }

//    cuda_err = cuda_malloc_vectors(Nrays,&dev_Nrays,lambda_vect_flag,&dev_lambda);
    cuda_err = cudaMallocChunk(Nrays,&dev_Nrays,1,false,&dev_lambda);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_BAD_ALLOC;
    }

    cuda_err = cuda_kernel_props(dev_Nrays,&N_cuda_blocks,&N_cuda_threads);
    if ( cuda_err != cudaSuccess ) {
//        cuda_free_mem(lambda_vect_flag,dev_lambda);
        cudaFreeChunk(1,dev_lambda);
        return ENGINE_ERROR_FAILED;
    }

    N_chunks = Nrays/dev_Nrays;

    ret_err = ENGINE_ERROR_OK;


    start_elem = 0;
    for ( size_t i = 0; i < N_chunks; ++i) {
#ifdef RT_NUM_DOUBLE
        curand_err = curandGenerateUniformDouble(gen,dev_lambda,dev_Nrays);
#else
        curand_err = curandGenerateUniform(gen,dev_lambda,dev_Nrays);
#endif
        if ( curand_err != CURAND_STATUS_SUCCESS ) {
            break;
        }

        random_range_beam_kernel<<<N_cuda_blocks,N_cuda_threads>>>(dev_Nrays,dev_lambda);

//        cuda_err = cuda_copy_mem(Nrays,start_elem,cudaMemcpyDeviceToHost,lambda_vect_flag,lambda,dev_lambda);
        cuda_err = cudaMemcpyChunk(dev_Nrays,cudaMemcpyDeviceToHost,1,false,
                                   lambda+start_elem,dev_lambda);
        if ( cuda_err != cudaSuccess ) {
            break;
        }

        start_elem += dev_Nrays;
    }

    // process rest of rays
    rest_Nrays = Nrays % dev_Nrays;
    if ( rest_Nrays && (ret_err == ENGINE_ERROR_OK) ) {
        cuda_err = cuda_kernel_props(rest_Nrays,&N_cuda_blocks,&N_cuda_threads);
        if ( cuda_err != cudaSuccess ) {
//            cuda_free_mem(lambda_vect_flag,dev_lambda);
            cudaFreeChunk(1,dev_lambda);
            return ENGINE_ERROR_FAILED;
        }

#ifdef RT_NUM_DOUBLE
        curand_err = curandGenerateUniformDouble(gen,dev_lambda,dev_Nrays);
#else
        curand_err = curandGenerateUniform(gen,dev_lambda,dev_Nrays);
#endif
        if ( curand_err != CURAND_STATUS_SUCCESS ) {
//            cuda_free_mem(lambda_vect_flag,dev_lambda);
            cudaFreeChunk(1,dev_lambda);
            return ENGINE_ERROR_FAILED;
        }

        random_range_beam_kernel<<<N_cuda_blocks,N_cuda_threads>>>(dev_Nrays,dev_lambda);

//        cuda_err = cuda_copy_mem(Nrays,start_elem,cudaMemcpyDeviceToHost,
//                                 lambda_vect_flag,lambda,dev_lambda);
        cuda_err = cudaMemcpyChunk(rest_Nrays,cudaMemcpyDeviceToHost,1,false,
                                   lambda+start_elem,dev_lambda);
        if ( cuda_err != cudaSuccess ) {
//            cuda_free_mem(lambda_vect_flag,dev_lambda);
            cudaFreeChunk(1,dev_lambda);
            return ENGINE_ERROR_FAILED;
        }
    }

//    cuda_free_mem(lambda_vect_flag,dev_lambda);
    cudaFreeChunk(1,dev_lambda);

    return ret_err;
}


