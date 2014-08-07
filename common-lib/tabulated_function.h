#ifndef TABULATED_FUNCTION_H
#define TABULATED_FUNCTION_H


#include "num_defs.h"

#include <vector>

using namespace std;


//
//  Argument vector in the input file must be in strickly increasing order!
//  This fact will be not checked in the class realization!
//


class TabulatedFunction
{
public:
    enum TabulatedFunctionError {MemAllocError,BadValue,BadInputValues,FileFailure,NotEnoughPoints,InitFailure,BadArgument,OutOfRange};

    TabulatedFunction();
    explicit TabulatedFunction(const real_t &val);       // just one value, i.e., independant on argument function
    explicit TabulatedFunction(const char* filename);    // tabulated function is in ASCII-file
    TabulatedFunction(const vector<real_t> &x, const vector<real_t> &yval); // user values

    void Set(const real_t &val);
    void Set(const size_t Nelements, const real_t *xval, const real_t *yval);

    void Load(const char* filename);

    real_t operator [](const real_t &xval);

    class bad_tab_func: public exception
    {
    public:
        bad_tab_func(TabulatedFunctionError err);
        TabulatedFunctionError bad_tab_func_error() const;
    private:
        TabulatedFunctionError Error;
    };

private:
    bool SingleValue;

    vector<real_t> x;    // function argument
    vector<real_t> y;    // function value
    vector<real_t> dy2;  // second derivative of the function at the tabulated points


    void SplineInit();
    real_t SplineInterp(const real_t &xval);
};

#endif // TABULATED_FUNCTION_H
