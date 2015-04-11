#include "tabulated_function.h"

#include "rt_engine_errors.h"
#include "rt_engine_func.h"

#include "ascii_file.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <string>

#include <cstdlib>

            /*  TabulatedFunction::bad_tab_func class realization  */
TabulatedFunction::bad_tab_func::bad_tab_func(TabulatedFunction::TabulatedFunctionError err): Error(err)
{
}

TabulatedFunction::TabulatedFunctionError TabulatedFunction::bad_tab_func::bad_tab_func_error() const
{
    return Error;
}

            /*  TabulatedFunction class realization  */

TabulatedFunction::TabulatedFunction():
    x(vector<real_t>()), y(vector<real_t>()), dy2(vector<real_t>()), SingleValue(false)
{

}

TabulatedFunction::TabulatedFunction(const real_t &val):
    x(vector<real_t>()), y(vector<real_t>()), dy2(vector<real_t>()), SingleValue(true)
{
    Set(val);
}


TabulatedFunction::TabulatedFunction(const vector<real_t> &xval, const vector<real_t> &yval):
    x(xval), y(yval), dy2(vector<real_t>()), SingleValue(false)
{
    if ( x.size() != y.size() ) {
        x.clear();
        y.clear();
        throw bad_tab_func(TabulatedFunction::BadInputValues);
    }
}


TabulatedFunction::TabulatedFunction(const char *filename):
    x(vector<real_t>()), y(vector<real_t>()), dy2(vector<real_t>()), SingleValue(false)
{
    Load(filename);
}


void TabulatedFunction::Load(const char *filename)
{
    x.clear();
    y.clear();
    dy2.clear();
    SingleValue = false;

//    ifstream file;
    AsciiFile file;
    real_t xval, yval;

    file.open(filename);

    if ( file.fail() ) {
        throw TabulatedFunction::bad_tab_func(TabulatedFunction::FileFailure);
    }

    xval = 0.0;
    yval = 0.0;

    AsciiFile::AsciiFileFlag line_flag;
    try {
//        file >> xval >> yval;
//        line_flag = file.ReadLine(2,&xval,&yval);
//        while ( !file.eof() ) {
        while ( (line_flag = file.ReadLine(2,&xval,&yval)) != AsciiFile::Eof ) {
//            cout << "RET: " << line_flag << ": " << xval << "; " << yval << endl;
            if ( line_flag == AsciiFile::InvalidData ) throw TabulatedFunction::bad_tab_func(TabulatedFunction::BadValue);
            if ( line_flag == AsciiFile::DataString ) {
//                if ( (xval <= 0.0) || (yval <= 0.0) ) {
//                    throw TabulatedFunction::bad_tab_func(TabulatedFunction::BadValue);
//                }
                x.push_back(xval);
                y.push_back(yval);

                xval = 0.0;
                yval = 0.0;
            }
//            file >> xval >> yval;
//            line_flag = file.ReadLine(2,&xval,&yval);
        }
        if ( x.size() == 0 ) throw TabulatedFunction::bad_tab_func(TabulatedFunction::BadValue);
        if ( x.size() < 4 ) throw TabulatedFunction::bad_tab_func(TabulatedFunction::NotEnoughPoints);
    } catch (TabulatedFunction::bad_tab_func &ex) {
        file.close();
        throw;
    } catch (bad_alloc &ex) {
        file.close();
        throw;
    }

    file.close();

    try {
        dy2.resize(x.size(),0.0);
    } catch (bad_alloc &ex) {
        throw TabulatedFunction::bad_tab_func(TabulatedFunction::MemAllocError);
    }
    SplineInit();
//cout << "size X = " << x.size() << "\n";
}


void TabulatedFunction::Set(const real_t &val)
{
    x.clear();
    y.clear();
    dy2.clear();
    SingleValue = true;

    if ( val <= 0.0 ) throw TabulatedFunction::bad_tab_func(TabulatedFunction::BadValue);

    try {
        x.push_back(0.0); // special value 0
        y.push_back(val);
        dy2.push_back(0.0);
    } catch ( bad_alloc &ex ) {
        throw TabulatedFunction::bad_tab_func(TabulatedFunction::MemAllocError);
    }
}


void TabulatedFunction::Set(const size_t Nelements, const real_t *xval, const real_t *yval)
{
    if ( Nelements == 0 ) throw TabulatedFunction::bad_tab_func(TabulatedFunction::BadValue);
    if ( Nelements < 4 ) throw TabulatedFunction::bad_tab_func(TabulatedFunction::NotEnoughPoints);

    x.clear();
    y.clear();
    dy2.clear();
    SingleValue = false;

    try {
        x = vector<real_t>(xval,xval+Nelements);
        y = vector<real_t>(yval,yval+Nelements);
        dy2 = vector<real_t>(Nelements,0.0);
        SplineInit();
    } catch ( ... ) {
        throw;
    }
}


//extern RT_engine_error spline_init(size_t N, real_t *x, real_t *y, real_t *dy2);

void TabulatedFunction::SplineInit()
{
    RT_engine_error err = spline_init(x.size(),x.data(),y.data(),dy2.data());
    if ( err != ENGINE_ERROR_OK ) {
        if ( err == ENGINE_ERROR_BAD_ALLOC ) throw TabulatedFunction::bad_tab_func(TabulatedFunction::MemAllocError);
        throw TabulatedFunction::bad_tab_func(TabulatedFunction::InitFailure);
    }
}


real_t TabulatedFunction::operator [](const real_t &xval)
{
    if ( y.size() == 0 ) throw TabulatedFunction::bad_tab_func(TabulatedFunction::OutOfRange);

    if ( SingleValue ) return y[0];

    if ( (xval < x[0]) || (xval > x[x.size()-1]) ) {
        throw TabulatedFunction::bad_tab_func(TabulatedFunction::OutOfRange);
    }

    // find place of the input point in the tabulated argument vector

    size_t i_left, i_right, i;
    real_t h, a, b;

    i_left = 0;
    i_right = x.size();

    while ( (i_right - i_left) > 1 ) {
        i = (i_left + i_right) >> 1;
        if ( x[i] == xval ) return y[i];
        if ( x[i] > xval ) i_right = i; else i_left = i;
    }

    h = x[i_right] - x[i_left];
    if ( h == 0.0 ) {
        throw TabulatedFunction::bad_tab_func(TabulatedFunction::BadArgument);
    }
    a = (x[i_right] - xval)/h;
    b = (xval - x[i_left])/h;
    return a*y[i_left]+b*y[i_right]+((a*a*a-a)*dy2[i_left]+(b*b*b-b)*dy2[i_right])*(h*h)/6.0;
}
