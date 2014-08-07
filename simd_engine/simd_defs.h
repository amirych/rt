  /*****************************
  * SIMD operation definitions *
  *****************************/

#ifndef SIMD_DEFS_H

#define SIMD_DEFS_H

#include "../num_defs.h"
#include <immintrin.h>

const real_t PI  =3.141592653589793238463;

/* alignment and vector size */

#ifdef __AVX__

  #define USE_SIMD

  #define SIMD_STR "AVX"

  #define SIMD_ALIGN 32
  #define SIMD_VEC_LEN 32/sizeof(real_t) // real_t type is defined in num_defs.h
  
  #ifdef RT_NUM_DOUBLE
	typedef __m256d simd_real_t;
	#define simd_load(a) _mm256_load_pd(a)
	#define simd_store(a,b) _mm256_store_pd(a,b)
	#define simd_add(a,b) _mm256_add_pd(a,b)
	#define simd_sub(a,b) _mm256_sub_pd(a,b)
	#define simd_mul(a,b) _mm256_mul_pd(a,b)
	#define simd_div(a,b) _mm256_div_pd(a,b)
    #define simd_broadcast_scalar(a) _mm256_broadcast_sd(a)
    #define simd_sqrt(a) _mm256_sqrt_pd(a)
    #define simd_cmpgt(a,b) _mm256_cmp_pd(a,b,_CMP_GT_OS)
  #else
	typedef __m256 simd_real_t;
	#define simd_load(a) _mm256_load_ps(a)
	#define simd_store(a,b) _mm256_store_ps(a,b)
	#define simd_add(a,b) _mm256_add_ps(a,b)
	#define simd_sub(a,b) _mm256_sub_ps(a,b)
	#define simd_mul(a,b) _mm256_mul_ps(a,b)
	#define simd_div(a,b) _mm256_div_ps(a,b)
    #define simd_broadcast_scalar(a) _mm256_broadcast_ss(a)
    #define simd_sqrt(a) _mm256_sqrt_ps(a)
    #define simd_cmpgt(a,b) _mm256_cmp_ps(a,b,_CMP_GT_OS)
  #endif

#elif defined __SSE2__  // (it is always true on x86_64 platforms)
//  #include <xmmintrin.h>

  #define USE_SIMD

  #define SIMD_STR "SSE2"

  #define SIMD_ALIGN 16
  #define SIMD_VEC_LEN 16/sizeof(real_t) // real_t type is defined in num_defs.h

  #ifdef RT_NUM_DOUBLE
	typedef __m128d simd_real_t;
	#define simd_load(a) _mm_load_pd(a)
	#define simd_store(a,b) _mm_store_pd(a,b)
	#define simd_add(a,b) _mm_add_pd(a,b)
	#define simd_sub(a,b) _mm_sub_pd(a,b)
	#define simd_mul(a,b) _mm_mul_pd(a,b)
	#define simd_div(a,b) _mm_div_pd(a,b)
    #define simd_broadcast_scalar(a) _mm_load1_pd(a)
    #define simd_sqrt(a) _mm_sqrt_pd(a)
    #define simd_cmpgt(a,b) _mm_cmpgt_pd(a,b)
  #else
	typedef __m128 simd_real_t;
	#define simd_load(a) _mm_load_ps(a)
	#define simd_store(a,b) _mm_store_ps(a,b)
	#define simd_add(a,b) _mm_add_ps(a,b)
	#define simd_sub(a,b) _mm_sub_ps(a,b)
	#define simd_mul(a,b) _mm_mul_ps(a,b)
	#define simd_div(a,b) _mm_div_ps(a,b)
    #define simd_broadcast_scalar(a) _mm_load1_ps(a)
    #define simd_sqrt(a) _mm_sqrt_ps(a)
    #define simd_cmpgt(a,b) _mm_cmpgt_ps(a,b)
  #endif
  
#else
  #define SIMD_STR "native"

  #define SIMD_ALIGN sizeof(real_t)
  #define SIMD_VEC_LEN 1
#endif

#endif // end of SIMD_DEFS_H
