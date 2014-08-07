#ifndef RT_ENGINE_SIMD_COMMON_H
#define RT_ENGINE_SIMD_COMMON_H

#include "num_defs.h"

class SIMD_vector
{
public:
    SIMD_vector();
    SIMD_vector(const real_t scalar);

    friend SIMD_vector
    ~SIMD_vector();
private:
    real_t *data;
};

#endif // RT_ENGINE_SIMD_COMMON_H
