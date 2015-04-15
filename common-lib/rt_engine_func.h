#ifndef RT_ENGINE_FUNC_H
#define RT_ENGINE_FUNC_H

#include "num_defs.h"
#include "rt_engine_errors.h"
#include "surface.h"

#include <list>
#include <string>

using namespace std;


 /***********************************************************************
 *                                                                      *
 *       The prototypes of the ray-tracing engine functions.            *
 *              Specific realizations are in sources of                 *
 *                nVidia CUDA or Intel SIMD engines.                    *
 *                                                                      *
 ***********************************************************************/


            /*    Initialization and info functions    */

RT_engine_error rt_engine_init();
RT_engine_error rt_engine_info(list<string> &info_str);

            /*    Beam class functions    */

RT_engine_error beam_random(Beam::BeamShape shape, real_t* params, real_t* center,
                            size_t Nrays, real_t* X, real_t* Y, real_t* Z);

RT_engine_error beam_conic_cosins(size_t Nrays, real_t conic_angle, real_t *X, real_t *Y, real_t *cX, real_t *cY, real_t *cZ);

RT_engine_error beam_gauss_cosins(size_t Nrays, real_t fwhm, real_t *cX, real_t *cY, real_t *cZ);

RT_engine_error uniform_range_beam(size_t Nrays, real_t *range, real_t *lambda);

RT_engine_error random_range_beam(size_t Nrays, real_t *range, real_t *lambda);

RT_engine_error beam_transform(size_t Nrays, real_t *angles, real_t *dist, real_t *X, real_t *Y, real_t *Z,
                               real_t *cX, real_t *cY, real_t *cZ);

            /*    Surface class functions    */

RT_engine_error surface_intersection(Surface::SurfaceClass surf_class, real_t* params,
                                     size_t Nrays, real_t* X, real_t* Y, real_t* Z,
                                     real_t* cX, real_t* cY, real_t* cZ,
                                     size_t *N_bad, beam_flag_t* flag);


RT_engine_error surface_reflection(Surface::SurfaceClass surf_class, real_t* params, size_t Nrays,
                                   real_t* X, real_t* Y, real_t* Z,
                                   real_t* cX, real_t* cY, real_t* cZ);


RT_engine_error surface_refraction(Surface::SurfaceClass surf_class, real_t* params, size_t Nrays,
                                   real_t *X, real_t *Y, real_t *Z,
                                   real_t *cX, real_t *cY, real_t *cZ, real_t *lambda,
                                   TabulatedFunction &n1, TabulatedFunction &n2,
                                   size_t *N_bad, beam_flag_t *flag);


RT_engine_error surface_diffration(Surface::SurfaceClass surf_class, Grating::GratingType gr_type,
                                   real_t* params, long diff_order, real_t gr_constant,
                                   size_t Nrays,
                                   real_t *X, real_t *Y, real_t *Z,
                                   real_t *cX, real_t *cY, real_t *cZ, real_t *lambda,
                                   TabulatedFunction &n1, TabulatedFunction &n2,
                                   size_t *N_bad, beam_flag_t *flag);


RT_engine_error surface_QE(size_t N, real_t *lambda, real_t *spec, TabulatedFunction &QE);

RT_engine_error grating_energy_distr(size_t N, real_t *lambda,
                                     long order, real_t blaze_angle, real_t alpha,
                                     real_t gamma, real_t gr_const,
                                     real_t *energy_distr);



            /*      Memory driver functions     */

void rt_engine_free_vector(real_t* data);

real_t* rt_engine_alocate_vector(size_t N_elem);


            /*      Spline interpolation functions      */

RT_engine_error spline_init(size_t N, real_t *x, real_t *y, real_t *y2);


            /*      CCD-related functions       */

RT_engine_error ccd_coordinates(size_t Nrays, real_t *X, real_t *Y,
                                real_t ccd_xsize, real_t ccd_ysize,
                                real_t ccd_xpix, real_t ccd_ypix,
                                unsigned int *ccd_X, unsigned int *ccd_Y);

RT_engine_error compute_ccd_image(size_t Nrays, unsigned int *ccd_X, unsigned int *ccd_Y, real_t *lambda,
                                  TabulatedFunction &QE, TabulatedFunction &incident_spec,
                                  size_t ccd_xdim, size_t ccd_ydim, real_t ccd_xpix, real_t ccd_ypix,
                                  real_t *image);


#endif // RT_ENGINE_FUNC_H
