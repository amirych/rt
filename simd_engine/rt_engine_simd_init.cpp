#include<list>
#include<string>

#include "../rt_engine_func.h"
#include "../rt_engine_errors.h"

RT_engine_error rt_engine_init()
{

}

RT_engine_error rt_engine_info(list<string> &info_str)
{
    string str;
    info_str.clear();
#ifdef USE_SIMD
    info_str.push_back("INTEL SIMD engine");
    str = "Using SIMD extension: " + SIMD_STR;
    info_str.push_back(str);
#else
    info_str.push_back("Native engine realization");
    info_str.push_back("No SIMD extension is used");
#endif
}

