#include<cuda.h>
#include<cuda_runtime.h>

#include "../common-lib/num_defs.h"
#include "../common-lib/rt_engine_errors.h"
#include "../common-lib/surface.h"
#include "../common-lib/tabulated_function.h"

#include "rt_cuda_engine_common.h"

#include <stdio.h>


                    /****************************************
                    *                                       *
                    *   nVIDIA CUDA font-end realization    *
                    *     for Surface class functions       *
                    *                                       *
                    ****************************************/

// surface parameters declarations
__constant__ real_t c;        // 1/R, where R is curvature radius
__constant__ real_t k1;       // (k+1) where k = -e^2 and e is eccentricity
__constant__ real_t cR;       // 1/R, where R is radius of revolution for toric class surfaces
__constant__ real_t eps;      // tolerance, a small positive number for iterative algorithm
__constant__ size_t max_iter; // maximal number of iterations for iterative intersection algorithm
__constant__ size_t N_poly;   // number of polynom members (aspheric part of F(X,Y,Z))

__constant__ long dev_order;  // diffraction order

//static __device__ real_t *poly_coeffs = NULL;    // array of polynom coefficients (aspheric part of F(X,Y,Z))
static __device__ real_t poly_coeffs[MAX_POLY_COEFFS];    // array of polynom coefficients (aspheric part of F(X,Y,Z))

static __device__ size_t dev_N_bad;

                        /*-----------------------------------------*/
                        /*    Device routines for computing of     */
                        /*  surface functions partial derivatives  */
                        /*-----------------------------------------*/

// conic function derivatives
// NOTE: derivative of conic surface with respect to Z is 1.0 (M-component), so compute only 2 partial derivatives
__device__
static void conic_derivs(size_t idx, real_t* dFdX, real_t* dFdY, real_t* dev_X, real_t* dev_Y)
{
    real_t E, sum, r2, rr;

    r2 = dev_X[idx]*dev_X[idx] + dev_Y[idx]*dev_Y[idx];
#ifdef RT_NUM_DOUBLE
    E = c/sqrt(1.0-k1*c*c*r2);
#else
    E = c/sqrtf(1.0-k1*c*c*r2);
#endif
    sum = 0.0, rr = 1.0;
    for (size_t i = 1; i <= N_poly; ++i, rr *= r2 ) sum += i*poly_coeffs[i-1]*rr; // aspheric term
    E += 2.0*sum;
    *dFdX = -E*dev_X[idx];
    *dFdY = -E*dev_Y[idx];
}


// toric function derivatives
__device__
static void toric_derivs(size_t idx, real_t* dFdX, real_t* dFdY, real_t* dFdZ, real_t* dev_X, real_t* dev_Y, real_t* dev_Z)
{
    real_t Y2, sum, rr, Fy, dFy;

    Y2 = dev_Y[idx]*dev_Y[idx];
#ifdef RT_NUM_DOUBLE
    Fy = c*Y2/(1.0+sqrt(1.0-k1*c*c*Y2));
#else
    Fy = c*Y2/(1.0+sqrtf(1.0-k1*c*c*Y2));
#endif
    sum = 0.0, rr = Y2;
    for (size_t i = 0; i < N_poly; ++i, rr *= Y2 ) sum += poly_coeffs[i]*rr; // aspheric term
    Fy += sum;

#ifdef RT_NUM_DOUBLE
    dFy = c*dev_Y[idx]/sqrt(1.0-k1*c*c*Y2);
#else
//    dFy = c*dev_Y[idx]/sqrtf(1.0-k1*c*c*Y2);
    dFy = c*dev_Y[idx]*rsqrtf(1.0-k1*c*c*Y2);
#endif
    sum = 0.0, rr = 1.0;
    for (size_t i = 1; i <= N_poly; ++i, rr *= Y2 ) sum += i*poly_coeffs[i-1]*rr; // aspheric term
    dFy += 2.0*sum;


    *dFdX = -cR*dev_X[idx];
    *dFdY = (cR*Fy - 1.0)*dFy;
    *dFdZ = 1.0 - cR*dev_Z[idx];
}




                        /*-------------------------------------------*/
                        /*  Beam-and-surface intersection functions  */
                        /*-------------------------------------------*/


        /*  auxiliary functions  */

// the function copies surface parameters to device constant memory (and to special static array if aspheric part is presented)
__host__
static RT_engine_error copy_surf_params(Surface::SurfaceClass surf_class,real_t *params)
{
    real_t kk1;
    size_t NN;

    cudaError_t cuda_err;

    cuda_err = cudaMemcpyToSymbol(c,params,sizeof(real_t)); // 1/R, where R is curvature radius
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_BAD_ALLOC;
    }

    kk1 = params[1] + 1;
    cuda_err = cudaMemcpyToSymbol(k1,&kk1,sizeof(real_t));   // k+1, where k is conic constant
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_BAD_ALLOC;
    }

    NN = (size_t) params[2];                                 // number of polynom coefficients
    cuda_err = cudaMemcpyToSymbol(N_poly,&NN,sizeof(size_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_BAD_ALLOC;
    }


    if ( NN ) { // copy polynom coefficients to device global memory
        cuda_err = cudaMemcpyToSymbol(poly_coeffs,params+3,NN*sizeof(real_t));
        if ( cuda_err != cudaSuccess ) {
            return ENGINE_ERROR_FAILED;
        }
    }

    if ( surf_class == Surface::Toric ) {
        cuda_err = cudaMemcpyToSymbol(cR,params+(3+NN),sizeof(real_t));
        if ( cuda_err != cudaSuccess ) {
            return ENGINE_ERROR_BAD_ALLOC;
        }
    }

    return ENGINE_ERROR_OK;
}


        /*  kernel functions  */

// the function computes beam intersections with plane Z=0
__global__
static void plane_intersection(size_t N,
                               real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                               real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ,
                               beam_flag_t* dev_flag)
{
    real_t s0;
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    while ( idx < N ) {
        s0 = -dev_Z[idx]/dev_cZ[idx];
//        printf("s0 = %f\n",s0);
        // !!!!!!!!!!!! ***  NOTE ABOUT WARP DIVERGENCE! ***
        if ( s0 < 0.0 ) {   // distance must be non-negative by definition
            dev_flag[idx] = 0;
            ++dev_N_bad;
//            atomicAdd(&NN,1);
        } else {
            dev_X[idx] += dev_cX[idx]*s0;
            dev_Y[idx] += dev_cY[idx]*s0;
            dev_Z[idx] = 0.0;
        }
        idx += blockDim.x*gridDim.x;
    }
}



// the function computes beam intersection with paraboloid ( k = -1 )
__global__
static void paraboloid_intersection(size_t N,
                                    real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                                    real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ,
                                    beam_flag_t* dev_flag)
{
    real_t B,C,s;

    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    while ( idx < N ) {
        B = 2.0*(k1*dev_Z[idx]*dev_cZ[idx] + dev_X[idx]*dev_cX[idx] + dev_Y[idx]*dev_cY[idx] - dev_cZ[idx]/c);
        C = k1*dev_Z[idx]*dev_Z[idx] + dev_X[idx]*dev_X[idx] + dev_Y[idx]*dev_Y[idx] - 2.0*dev_Z[idx]/c;

//        printf("XYZ: %f %f %f, COS: %f %f %f\n",dev_X[idx],dev_Y[idx],dev_Z[idx],dev_cX[idx],dev_cY[idx],dev_cZ[idx]);
        // !!!!!!!!!!!! ***  NOTE ABOUT WARP DIVERGENCE! ***

        s = -C/B;

        if ( s < 0.0 ) { // distance must be non-negative by definition
            dev_flag[idx] = 0;
            ++dev_N_bad;
        } else {
            // compute coordinates of intersection
            dev_X[idx] += s*dev_cX[idx];
            dev_Y[idx] += s*dev_cY[idx];
            dev_Z[idx] += s*dev_cZ[idx];
        }
        idx += blockDim.x*gridDim.x;
    }
}

// intersection algorithm for conic class surfaces (without additional polynom) except paraboloid!
__global__
static void conic_intersection(size_t N,
                               real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                               real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ,
                               beam_flag_t* dev_flag)
{
    real_t A,B,C,D,s;
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    // parameters c and k1 are in device constant memory now
    while ( idx < N ) {
        // compute coefficients of quadratic equation
        A = k1*dev_cZ[idx]*dev_cZ[idx] + dev_cX[idx]*dev_cX[idx] + dev_cY[idx]*dev_cY[idx];
        B = 2.0*(k1*dev_Z[idx]*dev_cZ[idx] + dev_X[idx]*dev_cX[idx] + dev_Y[idx]*dev_cY[idx] - dev_cZ[idx]/c);
        C = k1*dev_Z[idx]*dev_Z[idx] + dev_X[idx]*dev_X[idx] + dev_Y[idx]*dev_Y[idx] - 2.0*dev_Z[idx]/c;

        // !!!!!!!!!!!! ***  NOTE ABOUT WARP DIVERGENCE! ***

        D = B*B - 4.0*A*C;

        if ( D < 0.0 ) { // it means there is no intersection
            dev_flag[idx] = 0;
            ++dev_N_bad;
        } else {
            if ( c > 0 ) { // convex surface
                s = (-B - sqrt(D))/2.0/A;
            } else {       // concave surface
                s = (-B + sqrt(D))/2.0/A;
            }
            // compute coordinates of intersection
            dev_X[idx] += s*dev_cX[idx];
            dev_Y[idx] += s*dev_cY[idx];
            dev_Z[idx] += s*dev_cZ[idx];
        }

        idx += blockDim.x*gridDim.x;
    }
}


// generic iterative intersection algorithm for conic class surfaces
__global__
static void iterative_conic_intersection(size_t N,
                                         real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                                         real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ,
                                         beam_flag_t* dev_flag)
{
    real_t s0,s1,s2,Xi,Yi,Zi,F,dF,r2,rr,E,sum;
    size_t Niter = 0;

    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    // parameters c and k1 are in device constant memory now
    while ( idx < N ) {
        // compute intersection coordinates with plane Z = 0
        s0 = - dev_Z[idx]/dev_cZ[idx];
        if ( s0 >= 0 ) { // distance must be non-negative
            s2 = 0.0;
            do {
                s1 = s2;
                Xi = dev_X[idx] + (s1+s0)*dev_cX[idx];
                Yi = dev_Y[idx] + (s1+s0)*dev_cY[idx];
                Zi = s1*dev_cZ[idx];

                // compute F(X,Y,Z)
                r2 = Xi*Xi + Yi*Yi;
#ifdef RT_NUM_DOUBLE
                F = Zi - c*r2/(1.0+sqrt(1.0-k1*c*c*r2));
#else
                F = Zi - c*r2/(1.0+sqrtf(1.0-k1*c*c*r2));
#endif
                sum = 0.0, rr = r2;
                for (size_t i = 0; i < N_poly; ++i, rr *= r2 ) sum += poly_coeffs[i]*rr; // aspheric term
                F -= sum;
                // compute dF(X,Y,Z)/ds
#ifdef RT_NUM_DOUBLE
                E = c/sqrt(1.0-k1*c*c*r2);
#else
                E = c/sqrtf(1.0-k1*c*c*r2);
#endif
                sum = 0.0, rr = 1.0;
                for (size_t i = 1; i <= N_poly; ++i, rr *= r2 ) sum += i*poly_coeffs[i-1]*rr; // aspheric term
                E += 2.0*sum;
                dF = -E*Xi*dev_cX[idx] - E*Yi*dev_cY[idx] + dev_cZ[idx];

#ifdef RT_NUM_DOUBLE
                if ( fabs(dF) < eps ) { // possibly grazing incidence
#else
                if ( fabsf(dF) < eps ) {
#endif
                    ++dev_N_bad;
                    dev_flag[idx] = 0;
                    break;
                }

                s2 = s1 - F/dF;

                ++Niter;
#ifdef RT_NUM_DOUBLE
            } while ( (fabs(s2-s1) >= eps) && (Niter <= max_iter) );
#else
            } while ( (fabsf(s2-s1) >= eps) && (Niter <= max_iter) );
#endif
            dev_X[idx] = Xi;
            dev_Y[idx] = Yi;
            dev_Z[idx] = Zi;

            if ( Niter > max_iter ) { // mark this ray as bad!!!
                ++dev_N_bad;
                dev_flag[idx] = 0;
            }
        } else {
            ++dev_N_bad;
            dev_flag[idx] = 0;
        }
        idx += blockDim.x*gridDim.x;
    }
}



// toric function computation
__device__
real_t toric_func(real_t s, real_t X0, real_t Y0, real_t Z0, real_t cX, real_t cY, real_t cZ)
{
    real_t X = X0 + s*cX;
    real_t Y = Y0 + s*cY;
    real_t Z = Z0 + s*cZ;

    real_t Y2 = Y*Y;

    real_t FY;

#ifdef RT_NUM_DOUBLE
    FY = c*Y2/(1.0+sqrt(1.0+k1*c*c*Y2));
#else
    FY = c*Y2/(1.0+sqrtf(1.0+k1*c*c*Y2));
#endif

    return Z - FY -0.5*cR*(X*X+Z*Z-FY*FY);
}

// new iterative intersection algorithm for toric class surfaces
__global__ void iter_toric(size_t N,
                           real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                           real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ,
                           beam_flag_t* dev_flag)
{
    real_t zn, xn, sn, sn1, yn, Fxn, Fyn, Fzn;

    size_t Niter;

    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    while ( idx < N ) {
        // compute intersection coordinates with plane Z = 0
        sn1 = -dev_Z[idx]/dev_cZ[idx];
        if ( sn1 >= 0 ) {
            Niter = 0;
            Fyn = 0.0;
            do {
                sn = sn1;

                Fxn = toric_func(sn,dev_X[idx],dev_Y[idx],dev_Z[idx],dev_cX[idx],dev_cY[idx],dev_cZ[idx]);
                if ( abs(Fyn-Fxn) <= eps ) break;

                zn = sn + Fxn;
                Fzn = toric_func(zn,dev_X[idx],dev_Y[idx],dev_Z[idx],dev_cX[idx],dev_cY[idx],dev_cZ[idx]);

//                yn = sn - Fxn*Fxn/(Fzn-Fxn);
//                Fyn = toric_func(yn,dev_X[idx],dev_Y[idx],dev_Z[idx],dev_cX[idx],dev_cY[idx],dev_cZ[idx]);

//                sn1 = yn - Fyn*(Fxn-Fzn)/(xn-zn)/(Fxn-Fyn)*(xn-yn)/(Fyn-Fzn)*(yn-zn);

                sn1 = sn - Fxn*Fxn/(Fzn-Fxn);
                Fyn = Fxn;
                ++Niter;
#ifdef RT_NUM_DOUBLE
            } while ( (fabs(sn1-sn) >= eps) && (Niter <= max_iter) );
#else
            } while ( (fabsf(sn1-sn) >= eps) && (Niter <= max_iter) );
#endif
            dev_X[idx] += sn1*dev_cX[idx];
            dev_Y[idx] += sn1*dev_cY[idx];
            dev_Z[idx] += sn1*dev_cZ[idx];
            printf("<><>  Fxn = %f (%d)\n",Fxn,Niter);
        } else {
            ++dev_N_bad;
            dev_flag[idx] = 0;
        }

        idx += blockDim.x*gridDim.x;
     }

}


// iterative intersection algorithm for toric class surfaces
__global__
static void iterative_toric_intersection(size_t N,
                                         real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                                         real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ,
                                         beam_flag_t* dev_flag)
{
    real_t s0,s1,s2,Xi,Yi,Zi,F,Fy,dF,rr,dFy,sum,Y2;
//    size_t Niter = 0;
    size_t Niter = 0;

    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    while ( idx < N ) {
        // compute intersection coordinates with plane Z = 0
        s0 = - dev_Z[idx]/dev_cZ[idx];
        if ( s0 >= 0 ) { // distance must be non-negative
            s2 = 0.0;
            Niter = 0;
            do {
                s1 = s2;
                Xi = dev_X[idx] + (s1+s0)*dev_cX[idx];
                Yi = dev_Y[idx] + (s1+s0)*dev_cY[idx];
                Zi = s1*dev_cZ[idx];

                // compute F(X,Y,Z)
                Y2 = Yi*Yi;
#ifdef RT_NUM_DOUBLE
                Fy = c*Y2/(1.0+sqrt(1.0-k1*c*c*Y2));
#else
                Fy = c*Y2/(1.0+sqrtf(1.0-k1*c*c*Y2));
#endif
                sum = 0.0, rr = Y2;
                for (size_t i = 0; i < N_poly; ++i, rr *= Y2 ) sum += poly_coeffs[i]*rr; // aspheric term
                Fy += sum;
                F = Zi - Fy - 0.5*cR*(Xi*Xi + Zi*Zi - Fy*Fy);
                // compute dF(X,Y,Z)/ds
#ifdef RT_NUM_DOUBLE
                dFy = c*Yi/sqrt(1.0-k1*c*c*Y2);
#else
                dFy = c*Yi/sqrtf(1.0-k1*c*c*Y2);
#endif
                sum = 0.0, rr = 1.0;
                for (size_t i = 1; i <= N_poly; ++i, rr *= Y2 ) sum += i*poly_coeffs[i-1]*rr; // aspheric term
                dFy += 2.0*sum;

                dF = -cR*Xi*dev_cX[idx] + (cR*Fy - 1.0)*dFy*dev_cY[idx] + (1.0 - cR*Zi)*dev_cZ[idx];

#ifdef RT_NUM_DOUBLE
                if ( fabs(dF) < eps ) { // possibly grazing incidence
#else
                if ( fabsf(dF) < eps ) {
#endif
                    ++dev_N_bad;
                    dev_flag[idx] = 0;
                    break;
                }

                s2 = s1 - F/dF;

                ++Niter;
#ifdef RT_NUM_DOUBLE
            } while ( (fabs(s2-s1) >= eps) && (Niter <= max_iter) );
#else
            } while ( (fabsf(s2-s1) >= eps) && (Niter <= max_iter) );
#endif
            dev_X[idx] = Xi;
            dev_Y[idx] = Yi;
            dev_Z[idx] = Zi;
            if ( Niter > max_iter ) { // mark this ray as bad!!!
                ++dev_N_bad;
                dev_flag[idx] = 0;
            }
        } else {
            ++dev_N_bad;
            dev_flag[idx] = 0;
        }
       idx += blockDim.x*gridDim.x;
    }
}


// iterative algorithm for plane Z=0 + polynom part
// coefficients and it number are in device constant memory
__global__
static void iterative_poly_intersection(size_t N,
                                        real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                                        real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ,
                                        beam_flag_t* dev_flag)
{
    real_t F, dF, r2, rr, s0, s1, s2, Xi, Yi, Zi;
    size_t Niter = 0;

    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    while ( idx < N ) {
        // compute intersection coordinates with plane Z = 0
        s0 = - dev_Z[idx]/dev_cZ[idx];
        if ( s0 > 0 ) { // distance must be non-negative
            s2 = 0.0;
            do {
                s1 = s2;
                Xi = dev_X[idx] + (s1+s0)*dev_cX[idx];
                Yi = dev_Y[idx] + (s1+s0)*dev_cY[idx];
                Zi = s1*dev_cZ[idx];

                // compute F(X,Y,Z)
                r2 = Xi*Xi + Yi*Yi;
                F = 0.0, rr = r2;
                for (size_t i = 0; i < N_poly; ++i, rr *= r2 ) F += poly_coeffs[i]*rr; // just aspheric term
                F = Zi - F;

                // compute dF(X,Y,Z)/ds
                dF = 0.0, rr = 1.0;
                for (size_t i = 1; i <= N_poly; ++i, rr *= r2 ) dF += i*poly_coeffs[i-1]*rr; // just aspheric term
                dF = -Xi*dF - Yi*dF + 1.0;

#ifdef RT_NUM_DOUBLE
                if ( fabs(dF) < eps ) { // possibly grazing incidence
#else
                if ( fabsf(dF) < eps ) {
#endif
                    ++dev_N_bad;
                    dev_flag[idx] = 0;
                    break;
                }

                s2 = s1 - F/dF;

                ++Niter;
#ifdef RT_NUM_DOUBLE
            } while ( (fabs(s2-s1) >= eps) && (Niter <= max_iter) );
#else
            } while ( (fabsf(s2-s1) >= eps) && (Niter <= max_iter) );
#endif
            dev_X[idx] = Xi;
            dev_Y[idx] = Yi;
            dev_Z[idx] = Zi;
            if ( Niter > max_iter ) { // mark this ray as bad!!!
                ++dev_N_bad;
                dev_flag[idx] = 0;
            }
        } else {
            ++dev_N_bad;
            dev_flag[idx] = 0;
        }
        idx += blockDim.x*gridDim.x;
    }
}

        /*  host functions  */

// intersection main cycle function (host-side). it rans CUDA kernel functions
__host__
static RT_engine_error engine_intersection_cycle(Surface::SurfaceClass surf_class, real_t* params,
                                                 size_t Nrays, size_t start_elem,
                                                 real_t* X, real_t* Y, real_t* Z,
                                                 real_t* cX, real_t* cY, real_t* cZ,
                                                 beam_flag_t* flag,
                                                 real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                                                 real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ,
                                                 beam_flag_t* dev_flag,
                                                 int N_cuda_blocks, int N_cuda_threads)
{
    cudaError_t cuda_err;

    size_t NN = (size_t) params[2]; // number polynom coefficients

//    cuda_err = cuda_copy_mem(Nrays,start_elem,cudaMemcpyHostToDevice,XYZcXcYcZFlag_vect_flags,
//                             X,dev_X,Y,dev_Y,Z,dev_Z,cX,dev_cX,cY,dev_cY,cZ,dev_cZ,flag,dev_flag);
    cuda_err = cudaMemcpyChunk(Nrays,cudaMemcpyHostToDevice,6,true,
                               dev_X,X+start_elem,dev_Y,Y+start_elem,
                               dev_Z,Z+start_elem,
                               dev_cX,cX+start_elem,dev_cY,cY+start_elem,
                               dev_cZ,cZ+start_elem,dev_flag,
                               flag+start_elem);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    switch (surf_class) {
        case Surface::Plane: {
            plane_intersection<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);
            break;
        }
        case Surface::Conic: {
            // first, some checks
            if ( NN == 0 ) { // no aspheric part
                if ( params[0] == 0 ) { // c == 0, it is just a plane
                    plane_intersection<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);
                    break;
                }
                if ( params[1] == -1 ) { // k = -1, it is a paraboloid (special algorithm)
                    paraboloid_intersection<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);
                    break;
                }
                conic_intersection<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);
            } else {
                if ( params[0] == 0 ) { // c == 0
                     // it is computationaly more cheaper realization
                    iterative_poly_intersection<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);
                    break;
                }
                iterative_conic_intersection<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);
            }
            break;
        }
        case Surface::Toric: {
//        iterative_toric_intersection<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);
        iter_toric<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);

//            if ( NN == 0 ) { // no aspheric part
//                if ( (params[0] == 0) &&  (params[3+NN] == 0) )
//                    plane_intersection<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag); // it is just a plane (c == 0, cR == 0)
//                else
//                    iterative_toric_intersection<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);
//            } else {
//                iterative_toric_intersection<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);
//            }
            break;
        }
        default: // it can be AUX surfaces
            break;
    }

//    cuda_err = cuda_copy_mem(Nrays,start_elem,cudaMemcpyDeviceToHost,XYZcXcYcZFlag_vect_flags,
//                             X,dev_X,Y,dev_Y,Z,dev_Z,cX,dev_cX,cY,dev_cY,cZ,dev_cZ,flag,dev_flag);
    cuda_err = cudaMemcpyChunk(Nrays,cudaMemcpyDeviceToHost,6,true,
                               X+start_elem,dev_X,Y+start_elem,dev_Y,
                               Z+start_elem,dev_Z,
                               cX+start_elem,dev_cX,cY+start_elem,dev_cY,
                               cZ+start_elem,dev_cZ,
                               flag+start_elem,dev_flag);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    return ENGINE_ERROR_OK;
}


// callable driver function
__host__
RT_engine_error surface_intersection(Surface::SurfaceClass surf_class, real_t* params,
                                     size_t Nrays, real_t* X, real_t* Y, real_t* Z,
                                     real_t* cX, real_t* cY, real_t* cZ,
                                     size_t *N_bad, beam_flag_t* flag)
{
    cudaError_t cuda_err;
    RT_engine_error ret_err;
    int N_cuda_blocks, N_cuda_threads;

    size_t dev_Nrays, N_chunks, rest_Nrays, start_elem;

    real_t iter_eps = 1.0e-7;
    size_t alg_max_iter = 100000;


    real_t *dev_X, *dev_Y, *dev_Z, *dev_cX, *dev_cY, *dev_cZ;
    beam_flag_t *dev_flag;


    cuda_err = cudaMemcpyToSymbol(eps,&iter_eps,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_BAD_ALLOC;
    }

    cuda_err = cudaMemcpyToSymbol(max_iter,&alg_max_iter,sizeof(size_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_BAD_ALLOC;
    }


//    cuda_err = cuda_malloc_vectors(Nrays,&dev_Nrays,XYZcXcYcZFlag_vect_flags,
//                                   &dev_X,&dev_Y,&dev_Z,&dev_cX,&dev_cY,&dev_cZ,&dev_flag);
    cuda_err = cudaMallocChunk(Nrays,&dev_Nrays,6,true,
                               &dev_X,&dev_Y,&dev_Z,&dev_cX,&dev_cY,&dev_cZ,&dev_flag);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_BAD_ALLOC;
    }

    cuda_err = cuda_kernel_props(dev_Nrays,&N_cuda_blocks,&N_cuda_threads);
    if ( cuda_err != cudaSuccess ) {
//        cuda_free_mem(XYZcXcYcZFlag_vect_flags,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);
        cudaFreeChunk(7,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);
        return ENGINE_ERROR_FAILED;
    }

    N_chunks = Nrays/dev_Nrays;


    ret_err = ENGINE_ERROR_OK;

    // copy surface parameters to device constant memory
    ret_err = copy_surf_params(surf_class,params);
    if ( ret_err != ENGINE_ERROR_OK ) {
//        cuda_free_mem(XYZcXcYcZFlag_vect_flags,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);
        cudaFreeChunk(7,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);
        return ret_err;
    }

    // copy number of bad rays to device memory

    cuda_err = cudaMemcpyToSymbol(dev_N_bad,N_bad,sizeof(size_t));
    if ( cuda_err != cudaSuccess ) {
//        cuda_free_mem(XYZcXcYcZFlag_vect_flags,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);
        cudaFreeChunk(7,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);
        return ENGINE_ERROR_FAILED;
    }

    // main cycle
    start_elem = 0;
    for ( size_t i = 0; i < N_chunks; ++i) {
        ret_err = engine_intersection_cycle(surf_class,params,dev_Nrays,start_elem,
                                            X,Y,Z,cX,cY,cZ,flag,
                                            dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag,
                                            N_cuda_blocks,N_cuda_threads);
        if ( ret_err != ENGINE_ERROR_OK ) break;
        start_elem += dev_Nrays;
    }

    // process rest of rays
    rest_Nrays = Nrays % dev_Nrays;
    if ( rest_Nrays && (ret_err == ENGINE_ERROR_OK) ) {
        cuda_err = cuda_kernel_props(rest_Nrays,&N_cuda_blocks,&N_cuda_threads);
        if ( cuda_err != cudaSuccess ) {
//            cuda_free_mem(XYZcXcYcZFlag_vect_flags,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);
            cudaFreeChunk(7,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);
//            cudaFree(poly_coeffs);
            return ENGINE_ERROR_FAILED;
        }
        ret_err = engine_intersection_cycle(surf_class,params,rest_Nrays,start_elem,
                                            X,Y,Z,cX,cY,cZ,flag,
                                            dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag,
                                            N_cuda_blocks,N_cuda_threads);
    }

//    cuda_free_mem(XYZcXcYcZFlag_vect_flags,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);
    cudaFreeChunk(7,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_flag);

    // copy number of bad rays from device memory
    cuda_err = cudaMemcpyFromSymbol(N_bad,dev_N_bad,sizeof(size_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    return ret_err;
}


                 /*******************************************
                 *  Reflection, refraction and diffraction  *
                 *               functions                  *
                 *******************************************/


    /*  REFLECTION FUNCTIONS */

__global__
static void plane_reflection_kernel(size_t N, real_t *dev_cZ)
{
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    // for plane normal vector has components (K,L,M) = (0,0,1)
    // so only direction cosin along Z-axis is changed: m' = m - 2*a*M, where a = m*M and finally m' = -m
    while ( idx < N ) {
        dev_cZ[idx] = -dev_cZ[idx];
        idx += blockDim.x*gridDim.x;
    }
}

__global__
static void conic_reflection_kernel(size_t N, real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                                    real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ)
{
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;
    real_t dFdX, dFdY, a2;

    while ( idx < N ) {
        // first, compute components of normal vector (K,L,M)
        // NOTE: derivative of conic surface with respect to Z is 1.0 (M-component)

        conic_derivs(idx,&dFdX,&dFdY, dev_X, dev_Y);

        a2 = 2.0*(dev_cX[idx]*dFdX + dev_cY[idx]*dFdY + dev_cZ[idx])/(dFdX*dFdX + dFdY*dFdY + 1.0); // 2*a = (k*K + l*L + m*M)/(K*K + L*L + M*M)

        dev_cX[idx] -= a2*dFdX; // k' = k - 2*a*K
        dev_cY[idx] -= a2*dFdY; // l' = l - 2*a*L
        dev_cZ[idx] -= a2;      // m' = m - 2*a*M

        idx += blockDim.x*gridDim.x;
    }
}

__global__
static void toric_reflection_kernel(size_t N, real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                                    real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ)
{
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;
    real_t dFdX, dFdY, dFdZ, a2;

    while ( idx < N ) {
        // first, compute components of normal vector (K,L,M)

        toric_derivs(idx,&dFdX,&dFdY,&dFdZ, dev_X, dev_Y, dev_Z);

        a2 = 2.0*(dev_cX[idx]*dFdX + dev_cY[idx]*dFdY + dev_cZ[idx]*dFdZ)/(dFdX*dFdX + dFdY*dFdY + dFdZ*dFdZ); // 2*a = (k*K + l*L + m*M)/(K*K + L*L + M*M)

        dev_cX[idx] -= a2*dFdX; // k' = k - 2*a*K
        dev_cY[idx] -= a2*dFdY; // l' = l - 2*a*L
        dev_cZ[idx] -= a2*dFdZ; // m' = m - 2*a*M

        idx += blockDim.x*gridDim.x;
    }
}

// reflection function main cycle
__host__
static RT_engine_error surface_reflection_cycle(Surface::SurfaceClass surf_class,
                                                size_t Nrays, size_t start_elem,
                                                real_t* X, real_t* Y, real_t* Z,
                                                real_t* cX, real_t* cY, real_t* cZ,
                                                real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                                                real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ,
                                                int N_cuda_blocks, int N_cuda_threads)
{
    cudaError_t cuda_err;

//    cuda_err = cuda_copy_mem(Nrays,start_elem,cudaMemcpyHostToDevice,XYZcXcYcZ_vect_flags,
//                             X,dev_X,Y,dev_Y,Z,dev_Z,cX,dev_cX,cY,dev_cY,cZ,dev_cZ);
    cuda_err = cudaMemcpyChunk(Nrays,cudaMemcpyHostToDevice,6,false,
                               dev_X,X+start_elem, dev_Y,Y+start_elem,
                               dev_Z,Z+start_elem,
                               dev_cX,cX+start_elem, dev_cY,cY+start_elem,
                               dev_cZ,cZ+start_elem);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    switch ( surf_class ) {
        case Surface::Plane: {
            plane_reflection_kernel<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_cZ);
            break;
        }
        case Surface::Conic: {
            conic_reflection_kernel<<<N_cuda_blocks,N_cuda_threads>>>(Nrays, dev_X, dev_Y, dev_Z, dev_cX, dev_cY, dev_cZ);
            break;
        }
        case Surface::Toric: {
            toric_reflection_kernel<<<N_cuda_blocks,N_cuda_threads>>>(Nrays, dev_X, dev_Y, dev_Z, dev_cX, dev_cY, dev_cZ);
            break;
        }
        default: break;
    }

//    cuda_err = cuda_copy_mem(Nrays,start_elem,cudaMemcpyDeviceToHost,cXcYcZ_vect_flags,
//                             cX,dev_cX,cY,dev_cY,cZ,dev_cZ);
    // copy only direction cosins
    cuda_err = cudaMemcpyChunk(Nrays,cudaMemcpyDeviceToHost,3,false,
                               cX+start_elem,dev_cX,
                               cY+start_elem,dev_cY,
                               cZ+start_elem,dev_cZ);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    return ENGINE_ERROR_OK;
}

// callable driver function
__host__
RT_engine_error surface_reflection(Surface::SurfaceClass surf_class, real_t* params, size_t Nrays,
                                   real_t* X, real_t* Y, real_t* Z,
                                   real_t* cX, real_t* cY, real_t* cZ)
{
    cudaError_t cuda_err;
    RT_engine_error ret_err;
    int N_cuda_blocks, N_cuda_threads;

    size_t dev_Nrays, N_chunks, rest_Nrays, start_elem;

    real_t *dev_X, *dev_Y, *dev_Z, *dev_cX, *dev_cY, *dev_cZ;

//    cuda_err = cuda_malloc_vectors(Nrays,&dev_Nrays,XYZcXcYcZ_vect_flags,
//                                   &dev_X,&dev_Y,&dev_Z,&dev_cX,&dev_cY,&dev_cZ);
    cuda_err = cudaMallocChunk(Nrays,&dev_Nrays,6,false,
                               &dev_X,&dev_Y,&dev_Z,&dev_cX,&dev_cY,&dev_cZ);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_BAD_ALLOC;
    }

    cuda_err = cuda_kernel_props(dev_Nrays,&N_cuda_blocks,&N_cuda_threads);
    if ( cuda_err != cudaSuccess ) {
//        cuda_free_mem(XYZcXcYcZ_vect_flags,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ);
        cudaFreeChunk(6,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ);
        return ENGINE_ERROR_FAILED;
    }

    N_chunks = Nrays/dev_Nrays;


    ret_err = ENGINE_ERROR_OK;

    // copy surface parameters to device constant memory
    ret_err = copy_surf_params(surf_class,params);
    if ( ret_err != ENGINE_ERROR_OK ) {
//        cuda_free_mem(XYZcXcYcZ_vect_flags,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ);
        cudaFreeChunk(6,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ);
        return ret_err;
    }

    // main cycle
    start_elem = 0;
    for ( size_t i = 0; i < N_chunks; ++i) {
        ret_err = surface_reflection_cycle(surf_class,dev_Nrays,start_elem,X,Y,Z,cX,cY,cZ,
                                           dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,
                                           N_cuda_blocks,N_cuda_threads);
        if ( ret_err != ENGINE_ERROR_OK ) break;
        start_elem += dev_Nrays;
    }

    // process rest of rays
    rest_Nrays = Nrays % dev_Nrays;
    if ( rest_Nrays && (ret_err == ENGINE_ERROR_OK) ) {
        cuda_err = cuda_kernel_props(rest_Nrays,&N_cuda_blocks,&N_cuda_threads);
        if ( cuda_err != cudaSuccess ) {
//            cuda_free_mem(XYZcXcYcZ_vect_flags,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ);
            cudaFreeChunk(6,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ);
            return ENGINE_ERROR_FAILED;
        }
        ret_err = surface_reflection_cycle(surf_class,rest_Nrays,start_elem,X,Y,Z,cX,cY,cZ,
                                           dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,
                                           N_cuda_blocks,N_cuda_threads);
    }

//    cuda_free_mem(XYZcXcYcZ_vect_flags,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ);
    cudaFreeChunk(6,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ);
    return ret_err;
}


    /*  REFRACTION FUNCTIONS */

__global__
static void plane_refraction_kernel(size_t N, real_t *dev_cX, real_t *dev_cY, real_t *dev_cZ,
                                    real_t *dev_mu, beam_flag_t* dev_flag)
{
    real_t a,b,a2,Gamma;
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    // for a plane a normal vector has components (K,L,M) = (0,0,1)
    // so only direction cosins are:
    //   k' = mu*k
    //   l' = mu*l
    //   m' = mu*m + Gamma, where k,l,m - direction cosins along X, Y and Z-axis of incident ray
    while ( idx < N ) {
        a = dev_mu[idx]*dev_cZ[idx];
        b = dev_mu[idx]*dev_mu[idx] - 1.0;

        a2 = a*a;
        if ( b > a2 ) { // total internal reflection!
            dev_flag[idx] = 0;
            ++dev_N_bad;
        } else {
            b = a2 - b; // discriminant (actually, b = D/4.0)
#ifdef RT_NUM_DOUBLE
            Gamma = -a + sqrt(b); // use of a greater root of quadratic equation
#else
            Gamma = -a + sqrtf(b); // use of a greater root of quadratic equation
#endif
            dev_cX[idx] *= dev_mu[idx];
            dev_cY[idx] *= dev_mu[idx];
            dev_cZ[idx] = dev_mu[idx]*dev_cZ[idx] + Gamma;
        }
        idx += blockDim.x*gridDim.x;
    }
}


__global__
static void conic_refraction_kernel(size_t N, real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                                    real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ,
                                    real_t *dev_mu, beam_flag_t* dev_flag)
{
    real_t a,b,a2,Gamma,dFdX,dFdY;

    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    while ( idx < N ) {
        // first, compute components of normal vector (K,L,M)
        // NOTE: derivative of conic surface with respect to Z is 1.0 (M-component)

        conic_derivs(idx,&dFdX,&dFdY, dev_X, dev_Y);

        Gamma = dFdX*dFdX + dFdY*dFdY + 1.0; // norm
        a = dev_mu[idx]*(dev_cX[idx]*dFdX + dev_cY[idx]*dFdY + dev_cZ[idx])/Gamma;
        b = (dev_mu[idx]*dev_mu[idx] - 1.0)/Gamma;

        a2 = a*a;
        if ( b > a2 ) { // total internal reflection!
            dev_flag[idx] = 0;
            ++dev_N_bad;
        } else {
            b = a2 - b; // discriminant (actually, b = D/4.0)
#ifdef RT_NUM_DOUBLE
            Gamma = -a + sqrt(b); // use of a greater root of quadratic equation
#else
            Gamma = -a + sqrtf(b); // use of a greater root of quadratic equation
#endif
            dev_cX[idx] = dev_mu[idx]*dev_cX[idx] + Gamma*dFdX;
            dev_cY[idx] = dev_mu[idx]*dev_cY[idx] + Gamma*dFdY;
            dev_cZ[idx] = dev_mu[idx]*dev_cZ[idx] + Gamma;
        }
        idx += blockDim.x*gridDim.x;
    }
}


__global__
static void toric_refraction_kernel(size_t N, real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                                    real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ,
                                    real_t *dev_mu, beam_flag_t* dev_flag)
{
    real_t a,b,a2,Gamma,dFdX,dFdY,dFdZ;

    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    while ( idx < N ) {
        // first, compute components of normal vector (K,L,M)

        toric_derivs(idx,&dFdX,&dFdY,&dFdZ, dev_X, dev_Y, dev_Z);

        Gamma = dFdX*dFdX + dFdY*dFdY + dFdZ*dFdZ; // norm
        a = dev_mu[idx]*(dev_cX[idx]*dFdX + dev_cY[idx]*dFdY + dev_cZ[idx]*dFdZ)/Gamma;
        b = (dev_mu[idx]*dev_mu[idx] - 1.0)/Gamma;

        a2 = a*a;
        if ( b > a2 ) { // total internal reflection!
            dev_flag[idx] = 0;
            ++dev_N_bad;
        } else {
            b = a2 - b; // discriminant (actually, b = D/4.0)
#ifdef RT_NUM_DOUBLE
            Gamma = -a + sqrt(b); // use of a greater root of quadratic equation
#else
            Gamma = -a + sqrtf(b); // use of a greater root of quadratic equation
#endif
            dev_cX[idx] = dev_mu[idx]*dev_cX[idx] + Gamma*dFdX;
            dev_cY[idx] = dev_mu[idx]*dev_cY[idx] + Gamma*dFdY;
            dev_cZ[idx] = dev_mu[idx]*dev_cZ[idx] + Gamma*dFdZ;
        }
        idx += blockDim.x*gridDim.x;
    }
}


// refraction function main cycle
__host__
static RT_engine_error surface_refraction_cycle(Surface::SurfaceClass surf_class,
                                                size_t Nrays, size_t start_elem,
                                                real_t* X, real_t* Y, real_t* Z,
                                                real_t* cX, real_t* cY, real_t* cZ, real_t *lambda,
                                                TabulatedFunction &n1, TabulatedFunction &n2, beam_flag_t* flag,
                                                real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                                                real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ, real_t *dev_mu,
                                                beam_flag_t* dev_flag,
                                                int N_cuda_blocks, int N_cuda_threads)
{
    cudaError_t cuda_err;
    real_t *host_mu;

    // first, compute refractive indices ratio mu = n1/n2

    size_t vec_len = sizeof(real_t)*Nrays;

    host_mu = (real_t*)malloc(vec_len);
    if ( host_mu == NULL ) return ENGINE_ERROR_BAD_ALLOC;

//#ifdef USE_OPENMP
//#pragma omp parallel for
//#endif
    for (size_t i = 0; i < Nrays; ++i) {
        host_mu[i] = n1[lambda[i+start_elem]];
        host_mu[i] /= n2[lambda[i+start_elem]];
    }

    // copy computed refractive indices ratio to the device memory

    cuda_err = cudaMemcpy(dev_mu,host_mu,vec_len,cudaMemcpyHostToDevice);
    if ( cuda_err != cudaSuccess ) {
        free(host_mu);
        return ENGINE_ERROR_FAILED;
    }

    free(host_mu);


//    cuda_err = cuda_copy_mem(Nrays,start_elem,cudaMemcpyHostToDevice,XYZcXcYcZFlag_vect_flags,
//                             X,dev_X,Y,dev_Y,Z,dev_Z,cX,dev_cX,cY,dev_cY,cZ,dev_cZ,flag,dev_flag);
    cuda_err = cudaMemcpyChunk(Nrays,cudaMemcpyHostToDevice,6,true,
                               dev_X,X+start_elem,dev_Y,Y+start_elem,
                               dev_Z,Z+start_elem,
                               dev_cX,cX+start_elem,dev_cY,cY+start_elem,
                               dev_cZ,cZ+start_elem,
                               dev_flag,flag+start_elem);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    switch ( surf_class ) {
        case Surface::Plane: {
            plane_refraction_kernel<<<N_cuda_blocks,N_cuda_threads>>>(Nrays, dev_cX, dev_cY, dev_cZ,dev_mu,dev_flag);
            break;
        }
        case Surface::Conic: {
            conic_refraction_kernel<<<N_cuda_blocks,N_cuda_threads>>>(Nrays, dev_X, dev_Y, dev_Z, dev_cX, dev_cY, dev_cZ, dev_mu, dev_flag);
            break;
        }
        case Surface::Toric: {
            toric_refraction_kernel<<<N_cuda_blocks,N_cuda_threads>>>(Nrays, dev_X, dev_Y, dev_Z, dev_cX, dev_cY, dev_cZ, dev_mu, dev_flag);
            break;
        }
        default: break;
    }

//    cuda_err = cuda_copy_mem(Nrays,start_elem,cudaMemcpyDeviceToHost,cXcYcZFlag_vect_flags,
//                             cX,dev_cX,cY,dev_cY,cZ,dev_cZ,flag,dev_flag);
    cuda_err = cudaMemcpyChunk(Nrays,cudaMemcpyDeviceToHost,3,true,
                               cX+start_elem,dev_cX,
                               cY+start_elem,dev_cY,
                               cZ+start_elem,dev_cZ,
                               flag+start_elem,dev_flag);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    return ENGINE_ERROR_OK;
}



// callable host function
__host__
RT_engine_error surface_refraction(Surface::SurfaceClass surf_class, real_t* params, size_t Nrays,
                                   real_t *X, real_t *Y, real_t *Z,
                                   real_t *cX, real_t *cY, real_t *cZ, real_t *lambda,
                                   TabulatedFunction &n1, TabulatedFunction &n2,
                                   size_t *N_bad, beam_flag_t *flag)
{
    cudaError_t cuda_err;
    RT_engine_error ret_err;
    int N_cuda_blocks, N_cuda_threads;

    size_t dev_Nrays, N_chunks, rest_Nrays, start_elem;

    real_t *dev_X, *dev_Y, *dev_Z, *dev_cX, *dev_cY, *dev_cZ, *dev_mu;
    beam_flag_t *dev_flag;

    // copy surface parameters to device constant memory
    ret_err = copy_surf_params(surf_class,params);
    if ( ret_err != ENGINE_ERROR_OK ) {
        return ret_err;
    }

    // copy number of bad rays to device memory
    cuda_err = cudaMemcpyToSymbol(dev_N_bad,N_bad,sizeof(size_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }


    // allocate device memory for all working and refrative indices ratio vectors

//    cuda_err = cuda_malloc_vectors(Nrays,&dev_Nrays,All_vect_flags,
//                                   &dev_X,&dev_Y,&dev_Z,&dev_cX,&dev_cY,&dev_cZ,&dev_mu,&dev_flag);
    cuda_err = cudaMallocChunk(Nrays,&dev_Nrays,7,true,
                               &dev_X,&dev_Y,&dev_Z,&dev_cX,&dev_cY,&dev_cZ,&dev_mu,&dev_flag);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_BAD_ALLOC;
    }

    cuda_err = cuda_kernel_props(dev_Nrays,&N_cuda_blocks,&N_cuda_threads);
    if ( cuda_err != cudaSuccess ) {
//        cuda_free_mem(All_vect_flags,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_mu,dev_flag);
        cudaFreeChunk(8,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_mu,dev_flag);
        return ENGINE_ERROR_FAILED;
    }

    N_chunks = Nrays/dev_Nrays;


    ret_err = ENGINE_ERROR_OK;

    // main cycle
    start_elem = 0;
    for ( size_t i = 0; i < N_chunks; ++i) {
        ret_err = surface_refraction_cycle(surf_class,dev_Nrays,start_elem,
                                           X,Y,Z,cX,cY,cZ,lambda,n1,n2,flag,
                                           dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_mu,dev_flag,
                                           N_cuda_blocks,N_cuda_threads);
        if ( ret_err != ENGINE_ERROR_OK ) break;
        start_elem += dev_Nrays;
    }

    // process rest of rays
    rest_Nrays = Nrays % dev_Nrays;
    if ( rest_Nrays && (ret_err == ENGINE_ERROR_OK) ) {
        cuda_err = cuda_kernel_props(rest_Nrays,&N_cuda_blocks,&N_cuda_threads);
        if ( cuda_err != cudaSuccess ) {
//            cuda_free_mem(All_vect_flags,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_mu,dev_flag);
            cudaFreeChunk(8,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_mu,dev_flag);
            return ENGINE_ERROR_FAILED;
        }
        ret_err = surface_refraction_cycle(surf_class,rest_Nrays,start_elem,
                                           X,Y,Z,cX,cY,cZ,lambda,n1,n2,flag,
                                           dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_mu,dev_flag,
                                           N_cuda_blocks,N_cuda_threads);
    }

//    cuda_free_mem(All_vect_flags,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_mu,dev_flag);
    cudaFreeChunk(8,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_mu,dev_flag);

    // copy number of bad rays from device memory
    cuda_err = cudaMemcpyFromSymbol(N_bad,dev_N_bad,sizeof(size_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    return ret_err;
}


            /*  DIFFRACTION FUNCTIONS  */

__global__
static void plane_diffr_trans_kernel(size_t N, real_t *dev_cX, real_t *dev_cY, real_t *dev_cZ,
                                     real_t *dev_mu, real_t *dev_L, beam_flag_t *dev_flag)
{
    real_t a,b,a2,Gamma;
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    // for a plane a normal vector has components (K,L,M) = (0,0,1)
    // so vector p is (u,v,w) = (1,0,0)
    // and finally,
    // only direction cosins along X and Z-axis are changed:
    //   k' = mu*k - L, where k - direction cosin along X-axis of incident ray
    //   m' = mu*m + Gamma, where m - direction cosin along Z-axis of incident ray
    while ( idx < N ) {
        a = dev_mu[idx]*dev_cZ[idx];
        b = dev_mu[idx]*dev_mu[idx] - 1.0 + dev_L[idx]*dev_L[idx] -
                2.0*dev_mu[idx]*dev_L[idx]*dev_cX[idx];

        a2 = a*a;
        if ( b > a2 ) { // total internal reflection!
            dev_flag[idx] = 0;
            ++dev_N_bad;
        } else {
            b = a2 - b; // discriminant (actually, b = D/4.0)
#ifdef RT_NUM_DOUBLE
            Gamma = -a + sqrt(b); // use of a greater root of quadratic equation
#else
            Gamma = -a + sqrtf(b); // use of a greater root of quadratic equation
#endif
            dev_cX[idx] = dev_mu[idx]*dev_cX[idx] - dev_L[idx];
            dev_cY[idx] *= dev_mu[idx];
            dev_cZ[idx] = dev_mu[idx]*dev_cZ[idx] + Gamma;
        }
        idx += blockDim.x*gridDim.x;
    }
}

__global__
static void plane_diffr_refl_kernel(size_t N, real_t *dev_cX, real_t *dev_cY, real_t *dev_cZ,
                                    real_t *dev_L, beam_flag_t *dev_flag)
{
    real_t a,b,a2,Gamma;
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    // for a plane a normal vector has components (K,L,M) = (0,0,1)
    // so vector p is (u,v,w) = (1,0,0)
    // and finally,
    // only direction cosins along X and Z-axis are changed:
    //   k' = mu*k - L, where k - direction cosin along X-axis of incident ray
    //   m' = mu*m + Gamma, where m - direction cosin along Z-axis of incident ray
    while ( idx < N ) {
        a = dev_cZ[idx];
        b = dev_L[idx]*dev_L[idx] - 2.0*dev_L[idx]*dev_cX[idx];

        a2 = a*a;
        if ( b > a2 ) { // total internal reflection!
            dev_flag[idx] = 0;
            ++dev_N_bad;
        } else {
            b = a2 - b; // discriminant (actually, b = D/4.0)
#ifdef RT_NUM_DOUBLE
            Gamma = -a - sqrt(b); // use of a smaller root of quadratic equation
#else
            Gamma = -a - sqrtf(b); // use of a smaller root of quadratic equation
#endif
            dev_cX[idx] -= dev_L[idx];
//            dev_cY[idx] *= dev_mu[idx];
            dev_cZ[idx] += Gamma;
        }
        idx += blockDim.x*gridDim.x;
    }
}

__global__
static void conic_diffr_trans_kernel(size_t N, real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                                     real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ,
                                     real_t *dev_mu, real_t *dev_L, beam_flag_t* dev_flag)
{
    real_t a,b,a2,Gamma,dFdX,dFdY,u,v,w;

    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    while ( idx < N ) {
        // first, compute components of normal vector (K,L,M)
        // NOTE: derivative of conic surface with respect to Z is 1.0 (M-component)

        conic_derivs(idx,&dFdX,&dFdY, dev_X, dev_Y);

        // compute components of vector p
        Gamma = dFdX/(dFdY*dFdY + 1.0); // K/(L^2+M^2)
#ifdef RT_NUM_DOUBLE
        u = 1.0/sqrt(1.0+dFdX*Gamma); // u = 1/sqrt(1+K^2/(L^2+M^2))
#else
        u = 1.0/sqrtf(1.0+dFdX*Gamma);
#endif
        v = -dFdY*u*Gamma; // v = -K*L*u/(L^2+M^2)
        w = -u*Gamma;      // v = -K*M*u/(L^2+M^2)

        Gamma = dFdX*dFdX + dFdY*dFdY + 1.0; // norm

        a = dev_mu[idx]*(dev_cX[idx]*dFdX + dev_cY[idx]*dFdY + dev_cZ[idx])/Gamma;
        b = (dev_mu[idx]*dev_mu[idx] - 1.0 + dev_L[idx]*dev_L[idx] -
             2.0*dev_L[idx]*dev_mu[idx]*(dev_cX[idx]*u + dev_cY[idx]*v + dev_cZ[idx]*w))/Gamma;

        a2 = a*a;
        if ( b > a2 ) { // total internal reflection!
            dev_flag[idx] = 0;
            ++dev_N_bad;
        } else {
            b = a2 - b; // discriminant (actually, b = D/4.0)
#ifdef RT_NUM_DOUBLE
            Gamma = -a + sqrt(b); // use of a greater root of quadratic equation
#else
            Gamma = -a + sqrtf(b); // use of a greater root of quadratic equation
#endif
            dev_cX[idx] = dev_mu[idx]*dev_cX[idx] - dev_L[idx]*u + Gamma*dFdX;
            dev_cY[idx] = dev_mu[idx]*dev_cY[idx] - dev_L[idx]*v + Gamma*dFdY;
            dev_cZ[idx] = dev_mu[idx]*dev_cZ[idx] - dev_L[idx]*w + Gamma;
        }
        idx += blockDim.x*gridDim.x;
    }
}

__global__
static void conic_diffr_refl_kernel(size_t N, real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                                    real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ,
                                    real_t *dev_L, beam_flag_t* dev_flag)
{
    real_t a,b,a2,Gamma,dFdX,dFdY,u,v,w;

    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    while ( idx < N ) {
        // first, compute components of normal vector (K,L,M)
        // NOTE: derivative of conic surface with respect to Z is 1.0 (M-component)

        conic_derivs(idx,&dFdX,&dFdY, dev_X, dev_Y);

        // compute components of vector p
        Gamma = dFdX/(dFdY*dFdY + 1.0); // K/(L^2+M^2)
#ifdef RT_NUM_DOUBLE
        u = 1.0/sqrt(1.0+dFdX*Gamma); // u = 1/sqrt(1+K^2/(L^2+M^2))
#else
        u = 1.0/sqrtf(1.0+dFdX*Gamma);
#endif
        v = -dFdY*u*Gamma; // v = -K*L*u/(L^2+M^2)
        w = -u*Gamma;      // v = -K*M*u/(L^2+M^2)

        Gamma = dFdX*dFdX + dFdY*dFdY + 1.0; // norm

        a = (dev_cX[idx]*dFdX + dev_cY[idx]*dFdY + dev_cZ[idx])/Gamma;
        b = (dev_L[idx]*dev_L[idx] -
             2.0*dev_L[idx]*(dev_cX[idx]*u + dev_cY[idx]*v + dev_cZ[idx]*w))/Gamma;

        a2 = a*a;
        if ( b > a2 ) { // total internal reflection!
            dev_flag[idx] = 0;
            ++dev_N_bad;
        } else {
            b = a2 - b; // discriminant (actually, b = D/4.0)
#ifdef RT_NUM_DOUBLE
            Gamma = -a - sqrt(b); // use of a smaller root of quadratic equation
#else
            Gamma = -a - sqrtf(b); // use of a smaller root of quadratic equation
#endif
            dev_cX[idx] = dev_cX[idx] -dev_L[idx]*u + Gamma*dFdX;
            dev_cY[idx] = dev_cY[idx] -dev_L[idx]*v + Gamma*dFdY;
            dev_cZ[idx] = dev_cZ[idx] -dev_L[idx]*w + Gamma;
//            dev_cX[idx] += -dev_L[idx]*u + Gamma*dFdX;
//            dev_cY[idx] += -dev_L[idx]*v + Gamma*dFdY;
//            dev_cZ[idx] += -dev_L[idx]*w + Gamma;
        }
        idx += blockDim.x*gridDim.x;
    }
}

__global__
static void toric_diffr_trans_kernel(size_t N, real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                                     real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ,
                                     real_t *dev_mu, real_t *dev_L, beam_flag_t* dev_flag)
{
    real_t a,b,a2,Gamma,dFdX,dFdY,dFdZ,u,v,w;

    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    while ( idx < N ) {
        // first, compute components of normal vector (K,L,M)

        toric_derivs(idx,&dFdX,&dFdY,&dFdZ, dev_X, dev_Y, dev_Z);

        // compute components of vector p
        Gamma = dFdX/(dFdY*dFdY + dFdZ+dFdZ); // K/(L^2+M^2)
#ifdef RT_NUM_DOUBLE
        u = 1.0/sqrt(1.0+dFdX*Gamma); // u = 1/sqrt(1+K^2/(L^2+M^2))
#else
        u = 1.0/sqrtf(1.0+dFdX*Gamma);
#endif
        v = -dFdY*u*Gamma; // v = -K*L*u/(L^2+M^2)
        w = -dFdZ*u*Gamma; // v = -K*M*u/(L^2+M^2)

        Gamma = dFdX*dFdX + dFdY*dFdY + dFdZ*dFdZ; // norm

        a = dev_mu[idx]*(dev_cX[idx]*dFdX + dev_cY[idx]*dFdY + dev_cZ[idx]*dFdZ)/Gamma;
        b = (dev_mu[idx]*dev_mu[idx] - 1.0 + dev_L[idx]*dev_L[idx] -
             2.0*dev_L[idx]*dev_mu[idx]*(dev_cX[idx]*u + dev_cY[idx]*v + dev_cZ[idx]*w))/Gamma;

        a2 = a*a;
        if ( b > a2 ) { // total internal reflection!
            dev_flag[idx] = 0;
            ++dev_N_bad;
        } else {
            b = a2 - b; // discriminant (actually, b = D/4.0)
#ifdef RT_NUM_DOUBLE
            Gamma = -a + sqrt(b); // use of a greater root of quadratic equation
#else
            Gamma = -a + sqrtf(b); // use of a greater root of quadratic equation
#endif
            dev_cX[idx] = dev_mu[idx]*dev_cX[idx] - dev_L[idx]*u + Gamma*dFdX;
            dev_cY[idx] = dev_mu[idx]*dev_cY[idx] - dev_L[idx]*v + Gamma*dFdY;
            dev_cZ[idx] = dev_mu[idx]*dev_cZ[idx] - dev_L[idx]*w + Gamma*dFdZ;
        }
        idx += blockDim.x*gridDim.x;
    }
}

__global__
static void toric_diffr_refl_kernel(size_t N, real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                                    real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ,
                                    real_t *dev_L, beam_flag_t* dev_flag)
{
    real_t a,b,a2,Gamma,dFdX,dFdY,dFdZ,u,v,w;

    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    while ( idx < N ) {
        // first, compute components of normal vector (K,L,M)

        toric_derivs(idx,&dFdX,&dFdY,&dFdZ, dev_X, dev_Y, dev_Z);

        // compute components of vector p
        Gamma = dFdX/(dFdY*dFdY + dFdZ+dFdZ); // K/(L^2+M^2)
#ifdef RT_NUM_DOUBLE
        u = 1.0/sqrt(1.0+dFdX*Gamma); // u = 1/sqrt(1+K^2/(L^2+M^2))
#else
//        u = 1.0/sqrtf(1.0+dFdX*Gamma);
        u = rsqrtf(1.0+dFdX*Gamma);
#endif
        v = -dFdY*u*Gamma; // v = -K*L*u/(L^2+M^2)
        w = -dFdZ*u*Gamma; // v = -K*M*u/(L^2+M^2)

        Gamma = dFdX*dFdX + dFdY*dFdY + dFdZ*dFdZ; // norm

        a = (dev_cX[idx]*dFdX + dev_cY[idx]*dFdY + dev_cZ[idx]*dFdZ)/Gamma;
        b = (dev_L[idx]*dev_L[idx] -
             2.0*dev_L[idx]*(dev_cX[idx]*u + dev_cY[idx]*v + dev_cZ[idx]*w))/Gamma;

        a2 = a*a;
        if ( b > a2 ) { // total internal reflection!
            dev_flag[idx] = 0;
            ++dev_N_bad;
        } else {
            b = a2 - b; // discriminant (actually, b = D/4.0)
#ifdef RT_NUM_DOUBLE
            Gamma = -a - sqrt(b); // use of a smaller root of quadratic equation
#else
            Gamma = -a - sqrtf(b); // use of a smaller root of quadratic equation
#endif
            dev_cX[idx] += -dev_L[idx]*u + Gamma*dFdX;
            dev_cY[idx] += -dev_L[idx]*v + Gamma*dFdY;
            dev_cZ[idx] += -dev_L[idx]*w + Gamma*dFdZ;
        }
        idx += blockDim.x*gridDim.x;
    }
}

__host__
static RT_engine_error surface_diffraction_cycle(Surface::SurfaceClass surf_class,
                                                 Grating::GratingType gr_type,
                                                 long diff_order, real_t gr_constant,
                                                 size_t Nrays, size_t start_elem,
                                                 real_t* X, real_t* Y, real_t* Z,
                                                 real_t* cX, real_t* cY, real_t* cZ, real_t *lambda,
                                                 TabulatedFunction &n1, TabulatedFunction &n2,
                                                 beam_flag_t* flag,
                                                 real_t* dev_X, real_t* dev_Y, real_t* dev_Z,
                                                 real_t* dev_cX, real_t* dev_cY, real_t* dev_cZ,
                                                 real_t *dev_mu, real_t *dev_L,
                                                 beam_flag_t* dev_flag,
                                                 int N_cuda_blocks, int N_cuda_threads)
{
    cudaError_t cuda_err;
    real_t *host_mu, *host_L;

    // first, compute refractive indices ratio mu = n1/n2

    size_t vec_len = sizeof(real_t)*Nrays;

    cuda_err = cudaMemcpyChunk(Nrays,cudaMemcpyHostToDevice,6,true,
                               dev_X,X+start_elem,dev_Y,Y+start_elem,
                               dev_Z,Z+start_elem,
                               dev_cX,cX+start_elem,dev_cY,cY+start_elem,
                               dev_cZ,cZ+start_elem,
                               dev_flag,flag+start_elem);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    switch ( gr_type ) {
        case Grating::Transparent: {
            host_mu = (real_t*)malloc(vec_len);
            if ( host_mu == NULL ) return ENGINE_ERROR_BAD_ALLOC;

            host_L = (real_t*)malloc(vec_len);
            if ( host_L == NULL ) return ENGINE_ERROR_BAD_ALLOC;

    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
            for (size_t i = 0; i < Nrays; ++i) {
                host_mu[i] = n1[lambda[i+start_elem]];
                host_mu[i] /= n2[lambda[i+start_elem]];
                host_L[i] = diff_order*lambda[i+start_elem]/n2[lambda[i+start_elem]]/gr_constant;
            }
            cuda_err = cudaMemcpy(dev_mu,host_mu,vec_len,cudaMemcpyHostToDevice);
            if ( cuda_err != cudaSuccess ) {
                free(host_mu);
                return ENGINE_ERROR_FAILED;
            }
            cuda_err = cudaMemcpy(dev_L,host_L,vec_len,cudaMemcpyHostToDevice);
            if ( cuda_err != cudaSuccess ) {
                free(host_mu);
                free(host_L);
                return ENGINE_ERROR_FAILED;
            }

            switch ( surf_class ) {
                case Surface::Plane: {
                    plane_diffr_trans_kernel<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_cX,dev_cY,dev_cZ,dev_mu,dev_L,dev_flag);
                    break;
                }
                case Surface::Conic: {
                    conic_diffr_trans_kernel<<<N_cuda_blocks,N_cuda_threads>>>(Nrays, dev_X, dev_Y, dev_Z,
                                                                               dev_cX,dev_cY,dev_cZ,dev_mu,dev_L,dev_flag);
                    break;
                }
                case Surface::Toric: {
                    toric_diffr_trans_kernel<<<N_cuda_blocks,N_cuda_threads>>>(Nrays, dev_X, dev_Y, dev_Z,
                                                                               dev_cX,dev_cY,dev_cZ,dev_mu,dev_L,dev_flag);
                    break;
                }
                default: break;
            }
            break;
        }
        case Grating::Reflective: {
            host_L = (real_t*)malloc(vec_len);
            if ( host_L == NULL ) return ENGINE_ERROR_BAD_ALLOC;
            host_mu = NULL;
    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
            for (size_t i = 0; i < Nrays; ++i) {
                host_L[i] = diff_order*lambda[i+start_elem]/n2[lambda[i+start_elem]]/gr_constant;
            }
//            printf("lambda = %f, L = %f\n",lambda[Nrays-1+start_elem],host_L[Nrays-1]);
            cuda_err = cudaMemcpy(dev_L,host_L,vec_len,cudaMemcpyHostToDevice);
            if ( cuda_err != cudaSuccess ) {
                free(host_L);
                return ENGINE_ERROR_FAILED;
            }

            switch ( surf_class ) {
                case Surface::Plane: {
                    plane_diffr_refl_kernel<<<N_cuda_blocks,N_cuda_threads>>>(Nrays,dev_cX,dev_cY,dev_cZ,dev_L,dev_flag);
                    break;
                }
                case Surface::Conic: {
                    conic_diffr_refl_kernel<<<N_cuda_blocks,N_cuda_threads>>>(Nrays, dev_X, dev_Y, dev_Z,
                                                                               dev_cX,dev_cY,dev_cZ,dev_L,dev_flag);
                    break;
                }
                case Surface::Toric: {
                    toric_diffr_refl_kernel<<<N_cuda_blocks,N_cuda_threads>>>(Nrays, dev_X, dev_Y, dev_Z,
                                                                               dev_cX,dev_cY,dev_cZ,dev_L,dev_flag);
                    break;
                }
                default: break;
            }
            break;
        }
        default: return ENGINE_ERROR_BAD_VALUE;
    }

    free(host_mu);
    free(host_L);

    cuda_err = cudaMemcpyChunk(Nrays,cudaMemcpyDeviceToHost,3,true,
                               cX+start_elem,dev_cX,
                               cY+start_elem,dev_cY,
                               cZ+start_elem,dev_cZ,
                               flag+start_elem,dev_flag);
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    return ENGINE_ERROR_OK;
}

__host__
RT_engine_error surface_diffration(Surface::SurfaceClass surf_class, Grating::GratingType gr_type,
                                   real_t* params, long diff_order, real_t gr_constant,
                                   size_t Nrays,
                                   real_t *X, real_t *Y, real_t *Z,
                                   real_t *cX, real_t *cY, real_t *cZ, real_t *lambda,
                                   TabulatedFunction &n1, TabulatedFunction &n2,
                                   size_t *N_bad, beam_flag_t *flag)
{
    cudaError_t cuda_err;
    RT_engine_error ret_err;
    int N_cuda_blocks, N_cuda_threads;

    size_t dev_Nrays, N_chunks, rest_Nrays, start_elem;

    real_t *dev_X, *dev_Y, *dev_Z, *dev_cX, *dev_cY, *dev_cZ, *dev_mu, *dev_L;
    beam_flag_t *dev_flag;

    // copy surface parameters to device constant memory
    ret_err = copy_surf_params(surf_class,params);
    if ( ret_err != ENGINE_ERROR_OK ) {
        return ret_err;
    }

    // copy number of bad rays to device memory
    cuda_err = cudaMemcpyToSymbol(dev_N_bad,N_bad,sizeof(size_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }


    // allocate device memory for all working and refrative indices ratio vectors

//    cuda_err = cuda_malloc_vectors(Nrays,&dev_Nrays,All_vect_flags,
//                                   &dev_X,&dev_Y,&dev_Z,&dev_cX,&dev_cY,&dev_cZ,&dev_mu,&dev_flag);
    switch ( gr_type ) {
        case Grating::Transparent: {
            cuda_err = cudaMallocChunk(Nrays,&dev_Nrays,8,true,
                                       &dev_X,&dev_Y,&dev_Z,
                                       &dev_cX,&dev_cY,&dev_cZ,
                                       &dev_mu,&dev_L,&dev_flag);
            break;
        }
        case Grating::Reflective: {
            cuda_err = cudaMallocChunk(Nrays,&dev_Nrays,7,true,
                                       &dev_X,&dev_Y,&dev_Z,
                                       &dev_cX,&dev_cY,&dev_cZ,
                                       &dev_L,&dev_flag);
            dev_mu = NULL;
            break;
        }
        default: return ENGINE_ERROR_BAD_VALUE;
    }
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_BAD_ALLOC;
    }

    cuda_err = cuda_kernel_props(dev_Nrays,&N_cuda_blocks,&N_cuda_threads);
    if ( cuda_err != cudaSuccess ) {
//        cuda_free_mem(All_vect_flags,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_mu,dev_flag);
        cudaFreeChunk(9,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_mu,dev_L,dev_flag);
        return ENGINE_ERROR_FAILED;
    }

    N_chunks = Nrays/dev_Nrays;


    ret_err = ENGINE_ERROR_OK;

    // main cycle
    start_elem = 0;
    for ( size_t i = 0; i < N_chunks; ++i) {
        ret_err = surface_diffraction_cycle(surf_class, gr_type,
                                            diff_order, gr_constant,
                                            dev_Nrays,start_elem,
                                            X,Y,Z,cX,cY,cZ,lambda,n1,n2,flag,
                                            dev_X,dev_Y,dev_Z,
                                            dev_cX,dev_cY,dev_cZ,
                                            dev_mu,dev_L,dev_flag,
                                            N_cuda_blocks,N_cuda_threads);
        if ( ret_err != ENGINE_ERROR_OK ) break;
        start_elem += dev_Nrays;
    }

    // process rest of rays
    rest_Nrays = Nrays % dev_Nrays;
    if ( rest_Nrays && (ret_err == ENGINE_ERROR_OK) ) {
        cuda_err = cuda_kernel_props(rest_Nrays,&N_cuda_blocks,&N_cuda_threads);
        if ( cuda_err != cudaSuccess ) {
//            cuda_free_mem(All_vect_flags,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_mu,dev_flag);
            cudaFreeChunk(9,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_mu,dev_L,dev_flag);
            return ENGINE_ERROR_FAILED;
        }
        ret_err = surface_diffraction_cycle(surf_class, gr_type,
                                            diff_order, gr_constant,
                                            rest_Nrays,start_elem,
                                            X,Y,Z,cX,cY,cZ,lambda,n1,n2,flag,
                                            dev_X,dev_Y,dev_Z,
                                            dev_cX,dev_cY,dev_cZ,
                                            dev_mu,dev_L,dev_flag,
                                            N_cuda_blocks,N_cuda_threads);
    }

//    cuda_free_mem(All_vect_flags,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_mu,dev_flag);
    cudaFreeChunk(9,dev_X,dev_Y,dev_Z,dev_cX,dev_cY,dev_cZ,dev_mu,dev_L,dev_flag);

    // copy number of bad rays from device memory
    cuda_err = cudaMemcpyFromSymbol(N_bad,dev_N_bad,sizeof(size_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    return ret_err;
}



//
//
//
__global__
static void surface_QE_kernel(size_t N, real_t *dev_QE_curve, real_t *dev_spec)
{
    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    while ( idx < N ) {
        dev_spec[idx] *= dev_QE_curve[idx];
        idx += blockDim.x*gridDim.x;
    }
}

//
// surface_QE main cycle
//
__host__
static RT_engine_error surface_QE_cycle(size_t N, size_t start_elem,
                                        real_t *lambda, real_t *spec, TabulatedFunction &QE,
                                        real_t *dev_lambda, real_t *dev_spec, real_t *dev_QE_curve,
                                        int N_cuda_blocks, int N_cuda_threads)
{
    cudaError_t cuda_err;

    real_t *QE_curve;
    size_t vec_len = N*sizeof(real_t);

    QE_curve = (real_t*)malloc(vec_len);
    if ( QE_curve == NULL ) return ENGINE_ERROR_BAD_ALLOC;

    cuda_err = cudaMemcpyChunk(N,cudaMemcpyHostToDevice,2,false,dev_lambda,lambda+start_elem,dev_spec,spec+start_elem);
    if ( cuda_err != cudaSuccess ) return ENGINE_ERROR_FAILED;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < N; ++i) {
        QE_curve[i] = QE[lambda[i+start_elem]];
    }

    cuda_err = cudaMemcpy(dev_QE_curve,QE_curve,vec_len,cudaMemcpyHostToDevice);
    if ( cuda_err != cudaSuccess ) {
        free(QE_curve);
        return ENGINE_ERROR_FAILED;
    }

    surface_QE_kernel<<<N_cuda_blocks,N_cuda_threads>>>(N,dev_QE_curve,dev_spec);

    free(QE_curve);
    cuda_err = cudaMemcpy(spec+start_elem,dev_spec,vec_len,cudaMemcpyDeviceToHost);
    if ( cuda_err != cudaSuccess ) return ENGINE_ERROR_FAILED;

    return ENGINE_ERROR_OK;
}

//
// host callable function
//
__host__
RT_engine_error surface_QE(size_t N, real_t *lambda, real_t *spec, TabulatedFunction &QE)
{
    cudaError_t cuda_err;
    RT_engine_error ret_err;
    int N_cuda_blocks, N_cuda_threads;

    size_t dev_N, N_chunks, rest_N, start_elem;

    real_t *dev_lambda, *dev_spec, *dev_QE_curve;

    cuda_err = cudaMallocChunk(N,&dev_N,3,false,&dev_lambda,&dev_spec,&dev_QE_curve);
    if ( cuda_err != cudaSuccess ) return ENGINE_ERROR_BAD_ALLOC;

    cuda_err = cuda_kernel_props(dev_N,&N_cuda_blocks,&N_cuda_threads);
    if ( cuda_err != cudaSuccess ) {
        cudaFreeChunk(3,dev_lambda,dev_spec,dev_QE_curve);
        return ENGINE_ERROR_FAILED;
    }

    N_chunks = N/dev_N;


    ret_err = ENGINE_ERROR_OK;

    // main cycle
    start_elem = 0;
    for ( size_t i = 0; i < N_chunks; ++i) {
        ret_err = surface_QE_cycle(dev_N,start_elem,lambda,spec,QE,dev_lambda,dev_spec,dev_QE_curve,N_cuda_blocks,N_cuda_threads);
        if ( ret_err != ENGINE_ERROR_OK ) break;
        start_elem += dev_N;
    }

    // process rest of rays
    rest_N = N % dev_N;
    if ( rest_N && (ret_err == ENGINE_ERROR_OK) ) {
        cuda_err = cuda_kernel_props(rest_N,&N_cuda_blocks,&N_cuda_threads);
        if ( cuda_err != cudaSuccess ) {
            cudaFreeChunk(3,dev_lambda,dev_spec,dev_QE_curve);
            return ENGINE_ERROR_FAILED;
        }
        ret_err = surface_QE_cycle(rest_N,start_elem,lambda,spec,QE,dev_lambda,dev_spec,dev_QE_curve,N_cuda_blocks,N_cuda_threads);
        if ( ret_err != ENGINE_ERROR_OK ) {
            cudaFreeChunk(3,dev_lambda,dev_spec,dev_QE_curve);
            return ENGINE_ERROR_FAILED;
        }
    }

    cudaFreeChunk(3,dev_lambda,dev_spec,dev_QE_curve);
    return ret_err;
}


// square of sinc-function
__device__
static real_t sinc2(real_t arg)
{
    real_t si;

    if ( abs(arg) < 1.0E-2*RT_PI ) { // use of approximation polynom
        real_t a2 = arg*arg;
        si = (0.00761*a2-0.16605)*a2+1.0;
    } else {
        si = sin(arg)/arg;
    }

    return si*si;
}

// I = I1*I2/I3 function
// in max_iter constant the number of grating rules is stored
__device__
static real_t I_func(real_t lambda, real_t beta)
{
    real_t I1, I2, I3, u, term;

    term = RT_PI*eps/lambda;

    u = term*cos(beta)/cos(c-beta)*(sin(c-k1)+sin(c-beta));
    I1 = sinc2(u);

//    printf("beta=%f\n",beta*180/RT_PI);
//    printf("u=%e\n",u);
//    printf("I1=%e\n",I1);

    u = term*(sin(k1)+sin(beta));
    I2 = sinc2(max_iter*u);

    I3 = sinc2(u);

    term = cos(c-k1)*cos(beta)/cos(c-beta);
    term *= term;

//    printf("I2=%e\n",I2);
//    printf("I3=%e\n",I3);

    return term*I1*I2/I3;
}


__global__
static void gr_eff_kernel(size_t N, real_t *dev_lambda, real_t *dev_energy_distr)
{
    long order_min, order_max, m;
    real_t sin_alpha, cos_gamma, I_sum, beta, I;

    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    while ( idx < N ) {
#ifdef RT_NUM_DOUBLE
        sin_alpha = sin(k1);
        cos_gamma = cos(cR);

        order_min = __double2ll_rz(eps/dev_lambda[idx]*(sin_alpha-1)*cos_gamma);
        order_max = __double2ll_rz(eps/dev_lambda[idx]*(sin_alpha+1)*cos_gamma);
#else
        sin_alpha = sinf(k1);
        cos_gamma = cosf(cR);

        order_min = __float2ll_rz(eps/dev_lambda[idx]*(sin_alpha-1)*cos_gamma);
        order_max = __float2ll_rz(eps/dev_lambda[idx]*(sin_alpha+1)*cos_gamma);
#endif

        I_sum = 0.0;

        for ( m = order_min; m <= order_max; ++m ) {
#ifdef RT_NUM_DOUBLE
            beta = asin(m*dev_lambda[idx]/eps/cos_gamma-sin_alpha);
#else
            beta = asinf(m*dev_lambda[idx]/eps/cos_gamma-sin_alpha);
#endif
            I_sum += I_func(dev_lambda[idx],beta);
        }

#ifdef RT_NUM_DOUBLE
        beta = asin(dev_order*dev_lambda[idx]/(eps*cos_gamma) - sin_alpha);
#else
        beta = asinf(dev_order*dev_lambda[idx]/(eps*cos_gamma) - sin_alpha);
#endif

        I = I_func(dev_lambda[idx],beta);

        printf("beta = %f\n",beta*180/RT_PI);
        printf("I = %e\n",I);
        printf("Isum = %e\n",I_sum);

        dev_energy_distr[idx] *= I/I_sum;

        idx += blockDim.x*gridDim.x;
    }

}


/*
//
// I1 function (exactly follows to Tarasov but without I2 and I3 computations)
//
__device__
static real_t I1_func(real_t lambda, real_t beta)
{
    real_t term, I1;

    // compute U-function
#ifdef RT_NUM_DOUBLE
    I1 = RT_PI*eps/lambda*cos(beta)/cos(c-beta)*(sin(c-k1)+sin(c-beta));
//    printf("lambda = %f, U = %f",lambda,I1);
    term = cos(beta)*cos(c-k1)/cos(c-beta);
//    I1 *= term;
    if ( abs(I1) < RT_PI/2.0) {
        I1 *= I1;
        I1 = (0.00761*I1-0.16605)*I1+1;
    } else I1 = sin(I1)/(I1);
#else
    I1 = RT_PI*eps/lambda*cosf(beta)/cosf(c-beta)*(sinf(c-k1)+sinf(c-beta));
    I1 = sinf(I1)/(I1);
    term = cosf(beta)*cosf(c-k1)/cosf(c-beta);
#endif
    I1 *= term;
    I1 *= I1;

//    printf("lambda = %f, beta = %f, I1 = %f\n",lambda,beta*180/RT_PI,I1);

//    printf("EPS = %f; k1 = %f, c = %f\n",eps,k1,c);
    //    printf("lambda = %f, BETA = %f\n",lambda,beta*180/RT_PI);
    printf("BETA = %f\n",beta*180/RT_PI);
    return I1;
}
*/

//
// I1 function (modified by Zeddo)
//
__device__
static real_t I1_func(real_t lambda, real_t beta)
{
    real_t term1, term2, I1, d, dprime;

#ifdef RT_NUM_DOUBLE
    if ( beta >= c ) {
        I1 = RT_PI/lambda*(sin(c-k1)+sin(c-beta));
        d = eps*cos(beta)/cos(c-beta);
//        dprime = 0.0;
        term1 = cos(c-k1)*sin(I1*d)/I1;
        I1 = term1*term1;
    } else {
        if ( beta > (c-RT_PI/2.0) ) {
            d = eps*cos(c);
            dprime = eps*sin(c);
            I1 = RT_PI/lambda*(sin(c-k1)+sin(c-beta));
            term1 = cos(c-k1)*sin(I1*d)/I1;
            I1 = RT_PI/lambda*(cos(c-k1)+cos(c-beta));
            term2 = sin(c-k1)*sin(I1*dprime)/I1;
            I1 = term1*term1 + term2*term2;
        } else {
//            d = 0.0;
            I1 = RT_PI/lambda*(cos(c-k1)+cos(c-beta));
            dprime = eps*cos(beta)/sin(c-beta);
            term2 = sin(c-k1)*sin(I1*dprime)/I1;
            I1 = term2*term2;
        }
    }
#else
    if ( beta >= c ) {
        I1 = RT_PI/lambda*sinf(c-k1)+sinf(c-beta);
        d = eps*cosf(beta)/cosf(c-beta);
//        dprime = 0.0;
        term1 = cosf(c-k1)*sinf(I1*d)/I1;
        I1 = term1*term1;
    } else {
        if ( beta > (c-RT_PI/2.0) ) {
            d = eps*cosf(c);
            dprime = eps*sinf(c);
            I1 = RT_PI/lambda*sinf(c-k1)+sinf(c-beta);
            term1 = cosf(c-k1)*sinf(I1*d)/I1;
            I1 = RT_PI/lambda*cosf(c-k1)+cosf(c-beta);
            term2 = sinf(c-k1)*sinf(I1*dprime)/I1;
            I1 = term1*term1 + term2*term2;
        } else {
//            d = 0.0;
            I1 = RT_PI/lambda*cosf(c-k1)+cosf(c-beta);
            dprime = eps*cosf(beta)/sinf(c-beta);
            term2 = sinf(c-k1)*sinf(I1*dprime)/I1;
            I1 = term2*term2;
        }
    }
#endif

    return I1;
}


//
// blaze_angle in c
// alpha in k1
// gamma in cR
// grating constant in eps
//
__global__
static void grating_energy_distr_kernel(size_t N, real_t *dev_lambda, real_t *dev_energy_distr)
{
//    long long order_min, order_max, m;
    long order_min, order_max, m;
    real_t sin_alpha, cos_gamma, I1_norm, beta, I1;

    size_t idx = threadIdx.x + blockIdx.x*blockDim.x;

    while ( idx < N ) {
#ifdef RT_NUM_DOUBLE
        sin_alpha = sin(k1);
        cos_gamma = cos(cR);

        order_min = __double2ll_rz(eps/dev_lambda[idx]*(sin_alpha-1)*cos_gamma);
        order_max = __double2ll_rz(eps/dev_lambda[idx]*(sin_alpha+1)*cos_gamma);
#else
        sin_alpha = sinf(k1);
        cos_gamma = cosf(cR);

        order_min = __float2ll_rz(eps/dev_lambda[idx]*(sin_alpha-1)*cos_gamma);
        order_max = __float2ll_rz(eps/dev_lambda[idx]*(sin_alpha+1)*cos_gamma);
#endif
        I1_norm = 0.0;

        for ( m = order_min; m <= order_max; ++m ) {
#ifdef RT_NUM_DOUBLE
            beta = asin(m*dev_lambda[idx]/eps/cos_gamma-sin_alpha);
#else
            beta = asinf(m*dev_lambda[idx]/eps/cos_gamma-sin_alpha);
#endif
            I1_norm += I1_func(dev_lambda[idx],beta);
//            printf("lambda = %f, beta = %f, I1_norm = %f\n",dev_lambda[idx],beta*180/RT_PI,I1_norm);
        }
//        printf("lambda = %f, min = %d\n",dev_lambda[idx],order_min);
//        printf("lambda = %f, max = %d\n",dev_lambda[idx],order_max);

#ifdef RT_NUM_DOUBLE
        beta = asin(dev_order*dev_lambda[idx]/(eps*cos_gamma) - sin_alpha);
#else
        beta = asinf(dev_order*dev_lambda[idx]/(eps*cos_gamma) - sin_alpha);
#endif
        I1 = I1_func(dev_lambda[idx],beta);
        dev_energy_distr[idx] *= I1/I1_norm;

//        dev_energy_distr[idx] *= I1_func(dev_lambda[idx],beta)/I1_norm;

        printf("lambda = %f, beta = %f, I1 = %f, I1_norm = %f\n",dev_lambda[idx],beta*180/RT_PI,I1,I1_norm);
//        printf("eff = %f\n",I1_func(dev_lambda[idx],beta)/I1_norm);

        idx += blockDim.x*gridDim.x;
    }

}

//
//  grating energy distribution main cycle
//
__host__
static RT_engine_error grating_energy_distr_cycle(size_t N, size_t start_elem, real_t *lambda, real_t *energy_distr,
                                                  real_t *dev_lambda, real_t *dev_energy_distr,
                                                  int N_cuda_blocks, int N_cuda_threads)
{
    cudaError_t cuda_err;

    size_t vec_len = N*sizeof(real_t);

    cuda_err = cudaMemcpy(dev_lambda,lambda+start_elem,vec_len,cudaMemcpyHostToDevice);
    if ( cuda_err != cudaSuccess ) return ENGINE_ERROR_FAILED;

    cuda_err = cudaMemcpy(dev_energy_distr,energy_distr+start_elem,vec_len,cudaMemcpyHostToDevice);
    if ( cuda_err != cudaSuccess ) return ENGINE_ERROR_FAILED;

//printf("\n");

    grating_energy_distr_kernel<<<N_cuda_blocks,N_cuda_threads>>>(N,dev_lambda,dev_energy_distr);
//    gr_eff_kernel<<<N_cuda_blocks,N_cuda_threads>>>(N,dev_lambda,dev_energy_distr);

    cuda_err = cudaMemcpy(energy_distr+start_elem,dev_energy_distr,vec_len,cudaMemcpyDeviceToHost);
    if ( cuda_err != cudaSuccess ) return ENGINE_ERROR_FAILED;

    return ENGINE_ERROR_OK;
}

//
// host callable function
//
// all angles must in radians
//
__host__
RT_engine_error grating_energy_distr(size_t N, real_t *lambda,
                                     long order, real_t blaze_angle, real_t alpha,
                                     real_t gamma, real_t gr_const,
                                     real_t *energy_distr)
{
    cudaError_t cuda_err;
    RT_engine_error ret_err;
    int N_cuda_blocks, N_cuda_threads;

    size_t dev_N, N_chunks, rest_N, start_elem;

    size_t Nrules = 22500;
    cuda_err = cudaMemcpyToSymbol(max_iter,&Nrules,sizeof(size_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }


    real_t *dev_lambda, *dev_energy_distr;

//    printf("CONST = %f, ORDER = %d\n",gr_const,order);

    // copy angles and order
    cuda_err = cudaMemcpyToSymbol(c,&blaze_angle,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    cuda_err = cudaMemcpyToSymbol(k1,&alpha,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    cuda_err = cudaMemcpyToSymbol(cR,&gamma,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    cuda_err = cudaMemcpyToSymbol(eps,&gr_const,sizeof(real_t));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    cuda_err = cudaMemcpyToSymbol(dev_order,&order,sizeof(long));
    if ( cuda_err != cudaSuccess ) {
        return ENGINE_ERROR_FAILED;
    }

    cuda_err = cudaMallocChunk(N,&dev_N,2,false,&dev_lambda,&dev_energy_distr);
    if ( cuda_err != cudaSuccess ) return ENGINE_ERROR_BAD_ALLOC;

    cuda_err = cuda_kernel_props(dev_N,&N_cuda_blocks,&N_cuda_threads);
    if ( cuda_err != cudaSuccess ) {
        cudaFreeChunk(2,dev_lambda,dev_energy_distr);
        return ENGINE_ERROR_FAILED;
    }

    N_chunks = N/dev_N;


    ret_err = ENGINE_ERROR_OK;

    // main cycle
    start_elem = 0;
    for ( size_t i = 0; i < N_chunks; ++i) {
        ret_err = grating_energy_distr_cycle(dev_N,start_elem,lambda,energy_distr,
                                             dev_lambda,dev_energy_distr,N_cuda_blocks,N_cuda_threads);
        if ( ret_err != ENGINE_ERROR_OK ) break;
        start_elem += dev_N;
    }

    // process rest of rays
    rest_N = N % dev_N;
    if ( rest_N && (ret_err == ENGINE_ERROR_OK) ) {
        cuda_err = cuda_kernel_props(rest_N,&N_cuda_blocks,&N_cuda_threads);
        if ( cuda_err != cudaSuccess ) {
            cudaFreeChunk(2,dev_lambda,dev_energy_distr);
            return ENGINE_ERROR_FAILED;
        }
        ret_err = grating_energy_distr_cycle(rest_N,start_elem,lambda,energy_distr,
                                             dev_lambda,dev_energy_distr,N_cuda_blocks,N_cuda_threads);
        if ( ret_err != ENGINE_ERROR_OK ) {
            cudaFreeChunk(2,dev_lambda,dev_energy_distr);
            return ENGINE_ERROR_FAILED;
        }
    }

    cudaFreeChunk(2,dev_lambda,dev_energy_distr);

    return ret_err;
}
