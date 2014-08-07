#ifndef NUM_DEFS_H
#define NUM_DEFS_H

// what is a compiler I use
#if defined(__GNUC__) && !defined(USING_MACOSX) // GCC
    #define USING_LINUX
#elif defined(WIN32) && defined(_MSC_VER)       // Visual Studio
    #define USING_MSVC
#else
    #error "Unsuported compiler! Stop compilation!"
#endif


//typedef double real_t;
//#define RT_NUM_DOUBLE


typedef float real_t;

#define RT_PI 3.14159265359

//#ifdef RT_NUM_DOUBLE
////    static const real_t deg2rad = 3.14159265359/180.0;
//    static const real_t RT_PI = 3.14159265359;
//#else
////    static const real_t deg2rad = 3.141593f/180.0f;
//    static const real_t RT_PI = 3.14159265359f;
//#endif

//static const real_t deg2rad = RT_PI/180.0;
#define deg2rad RT_PI/180.0

typedef char beam_flag_t;

#endif // NUM_DEFS_H
