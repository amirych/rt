#include "../num_defs.h"
#include "simd_defs.h"
#include "../rt_engine_errors.h"
#include "../beam.h"


                        /***************************************
                        *                                      *
                        *   Intel SIMD font-end realization    *
                        *      for Beam class functions        *
                        *                                      *
                        ***************************************/


RT_engine_error beam_random(Beam::BeamShape shape, real_t* params, real_t* center,
                            size_t Nrays, real_t* X, real_t* Y, real_t* Z)
{
    size_t N_chunks = Nrays - (Nrays % SIMD_VEC_LEN);



}

RT_engine_error beam_conic_cosins(size_t Nrays, real_t conic_angle,
                                  real_t *X, real_t *Y, real_t *cX, real_t *cY, real_t *cZ)
{

}


RT_engine_error uniform_range_beam(size_t Nrays, real_t *range, real_t *lambda)
{

}

RT_engine_error random_range_beam(size_t Nrays, real_t *range, real_t *lambda)
{

}

RT_engine_error beam_transform(size_t Nrays, real_t *angles, real_t *dist, real_t *X, real_t *Y, real_t *Z,
                               real_t *cX, real_t *cY, real_t *cZ)
{

}

