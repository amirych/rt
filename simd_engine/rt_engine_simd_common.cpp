#include "../num_defs.h"

#include "simd_defs.h"

#include <cstdlib>


                /*  SIMD engine memory allocation  */

void rt_engine_free_vector(real_t* data)
{
#ifdef __linux__
    free(data);
#elif defined WIN32
    return _aligned_free(data);
#else                       // other (use valloc for page-aligned memory)
    free(data);
#endif
}

real_t* rt_engine_alocate_vector(size_t N_elem)
{
    size_t size = N_elem*sizeof(real_t);
#ifdef __linux__
    return static_cast<real_t*>(aligned_alloc(SIMD_ALIGN, size));
#elif defined WIN32
    return static_cast<real_t*>(_aligned_malloc(size, SIMD_ALIGN));
#else                       // other (use valloc for page-aligned memory?)
    return static_cast<real_t*>(valloc(size));
#endif

}
