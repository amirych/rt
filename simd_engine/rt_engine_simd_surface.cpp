#include "../num_defs.h"
#include "simd_defs.h"
#include "../rt_engine_errors.h"
#include "../surface.h"


                        /***************************************
                        *                                      *
                        *   Intel SIMD font-end realization    *
                        *     for Surface class functions      *
                        *                                      *
                        ***************************************/





RT_engine_error surface_intersection(Surface::SurfaceClass surf_class, real_t* params,
                                     size_t Nrays, real_t* X, real_t* Y, real_t* Z,
                                     real_t* cX, real_t* cY, real_t* cZ,
                                     size_t *N_bad, beam_flag_t* flag)
{

    size_t N_simd = Nrays - (Nrays % SIMD_VEC_LEN);
    real_t zero = 0.0;

#ifdef USE_SIMD
    real_t s[SIMD_VEC_LEN];
    simd_real_t simd_s;
    simd_real_t simd_X;
    simd_real_t simd_Y;
    simd_real_t simd_Z;
    simd_real_t simd_cX;
    simd_real_t simd_cY;
    simd_real_t simd_cZ;


#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < N_simd; i += SIMD_VEC_LEN) {

        switch ( surf_class ) {
            case Surface::Plane: { // trivial case of Z=0 plane
                simd_Z = simd_load(Z+i);
                simd_cZ = simd_load(cZ+i);
                simd_s = -simd_div(simd_Z,simd_cZ);      // -Z/cZ
                simd_Z = simd_broadcast_scalar(&zero);     // Z = 0
                simd_store(s,simd_s);

                for (int j = 0; j < SIMD_VEC_LEN; ++j) {
                    if (s[j] < 0) { // distance must be non-negative!
                        flag[i+j] = 0;
                        ++(*N_bad);
                    }
                }

                simd_X = simd_load(X+i);
                simd_Y = simd_load(Y+i);
                simd_cX = simd_load(cX+i);
                simd_cY = simd_load(cY+i);

                simd_cX = simd_mul(simd_cX,simd_s); // cX*s
                simd_cY = simd_mul(simd_cY,simd_s); // cY*s

                simd_X = simd_add(simd_X,simd_cX); // X_new = X + cX*s
                simd_Y = simd_add(simd_Y,simd_cY); // Y_new = Y + cY*s

                simd_store(X+i,simd_X);
                simd_store(Y+i,simd_Y);
                simd_store(Z+i,simd_Z);
                break;
            }
        case Surface::Conic: {
            break;
        }
        }
    }
#else // generic non-SIMD realization
    real_t s;

#endif
}


RT_engine_error surface_reflection(Surface::SurfaceClass surf_class, real_t* params, size_t Nrays,
                                   real_t* X, real_t* Y, real_t* Z,
                                   real_t* cX, real_t* cY, real_t* cZ)
{
}


RT_engine_error surface_refraction(Surface::SurfaceClass surf_class, real_t* params, size_t Nrays,
                                   real_t *X, real_t *Y, real_t *Z,
                                   real_t *cX, real_t *cY, real_t *cZ, real_t *lambda,
                                   TabulatedFunction &n1, TabulatedFunction &n2,
                                   size_t *N_bad, beam_flag_t *flag)
{

}

RT_engine_error surface_diffration(Surface::SurfaceClass surf_class, Grating::GratingType gr_type,
                                   real_t* params, long diff_order, real_t gr_constant,
                                   size_t Nrays,
                                   real_t *X, real_t *Y, real_t *Z,
                                   real_t *cX, real_t *cY, real_t *cZ, real_t *lambda,
                                   TabulatedFunction &n1, TabulatedFunction &n2,
                                   size_t *N_bad, beam_flag_t *flag)
{

}
