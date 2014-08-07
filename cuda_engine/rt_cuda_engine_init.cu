#include "../common-lib/rt_engine_errors.h"
#include "../common-lib/rt_engine_func.h"

#include <list>
#include <string>
#include <sstream>

RT_engine_error rt_engine_init()
{
    return ENGINE_ERROR_OK;
}


RT_engine_error rt_engine_info(list<string> &info_str)
{
    cudaError_t err;
    RT_engine_error ret_err;
    stringstream ss;

    // get device info
    cudaDeviceProp prop;

    err = cudaGetDeviceProperties(&prop,0);
    if ( err != cudaSuccess) {
        return ENGINE_ERROR_FAILED;
    }

    ret_err = ENGINE_ERROR_OK;

    try {
        info_str.clear();

        info_str.push_back("nVidia CUDA engine");

        ss << "Device name: " << prop.name;

        info_str.push_back(ss.str());

        ss.str("");
        ss << "Computing ability: " << prop.major << "." << prop.minor;

        info_str.push_back(ss.str());

        ss.str("");
        ss << "Global memory: " << prop.totalGlobalMem;
        info_str.push_back(ss.str());

        ss.str("");
        ss << "Number of multiprocessors: " << prop.multiProcessorCount;
        info_str.push_back(ss.str());
    } catch (bad_alloc &ex) {
        ret_err = ENGINE_ERROR_BAD_ALLOC;
    }

    return ret_err;
}
