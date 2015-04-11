#include "beam.h"

#include "rt_engine_func.h"
#include "rt_engine_errors.h"

#include "surface.h"

#include<iostream>
#include<math.h>
#include<cstdio>
#include<cstring>
#include<algorithm>

using namespace std;

        /********************************
        *                               *
        *   Beam class implementation   *
        *                               *
        ********************************/

        /*  Realization-dependant function declarations (nVidia CUDA or Intel SIMD)  */

        /*  Beam class exception implementation */

Beam::bad_beam::bad_beam(BeamError err): exception(), Error(err)
{
}

Beam::BeamError Beam::bad_beam::beam_error() const
{
    return Error;
}

        /*  Constructors and destructor  */

Beam::Beam(): Type(Beam::Parallel), Shape(Beam::Circle), Profile(Beam::Random), Range_Distr(Beam::RandomDistr),
    Params(vector<real_t>(3,0.0)), Center(vector<real_t>(3,0.0)), Nrays(0), N_good_rays(0), CompositeBeam(false),
    X(NULL), Y(NULL), Z(NULL), cX(NULL), cY(NULL), cZ(NULL), lambda(NULL), flag(NULL)
{
}


Beam::Beam(Beam::BeamType type, Beam::BeamShape shape, Beam::BeamProfile profile, Beam::BeamRangeDistr range_distr):
    Type(type), Shape(shape), Profile(profile), Range_Distr(range_distr),
    Params(vector<real_t>(3,0.0)), Center(vector<real_t>(3,0.0)),
    Nrays(0), N_good_rays(0), CompositeBeam(false), X(NULL), Y(NULL), Z(NULL), cX(NULL), cY(NULL), cZ(NULL), lambda(NULL), flag(NULL)
{
//    cerr << "Create beam!\n";
}


Beam::~Beam()
{
//    beam_free_mem(X,Y,Z,cX,cY,cZ);
//    delete[] X;
//    delete[] Y;
//    delete[] Z;

    destroy();
}


        /*  Beam class public methods   */


void Beam::SetParams(const vector<real_t> &params, bool composite_flag, size_t Nlambda)
{
    // check for input parameters
    if ( params.size() < 3 ) throw Beam::bad_beam(Beam::InvalidParams);
    if ( (Shape == Beam::Rectangle) && (params.size() < 4) ) throw Beam::bad_beam(Beam::InvalidParams);

    if ( (Type == Beam::Gauss) || Type == Beam::Conic ) {
        if ( (Shape == Beam::Circle) && (params.size() < 4) ) throw Beam::bad_beam(Beam::InvalidParams);
        if ( (Shape == Beam::Rectangle) && (params.size() < 5) ) throw Beam::bad_beam(Beam::InvalidParams);
    }

    if ( params[0] <= 0 ) throw Beam::bad_beam(Beam::InvalidParams); // number of rays or concentric shapes
    if ( (params[1] < 0) || (params[2] <= 0) ) throw Beam::bad_beam(Beam::InvalidParams); // geometric sizes must be greater than 0

    switch ( Shape ) {
        case Beam::Circle: {
            if ( (params[1] >= params[2]) ) throw Beam::bad_beam(Beam::InvalidParams); // inner circle must be less than outter
            break;
        }
        case Beam::Rectangle: {
            if ( params[3] <= 0 )  throw Beam::bad_beam(Beam::InvalidParams); // size ratio must be greater than 0
            break;
        }
        default: break;
    }

    Params = params;

    CompositeBeam = composite_flag;
    N_lambda = (Nlambda > 0) ? Nlambda : 1;
}


void Beam::SetCenter(const vector<real_t> &center)
{
    if ( center.size() < 3 ) throw Beam::bad_beam(Beam::InvalidParams);
    Center = center;
}


void Beam::SetRange(const vector<real_t> &range)
{
    if ( range.size() == 0 ) throw Beam::bad_beam(Beam::InvalidParams);
    if ( Range_Distr != Beam::Monochromatic ) {
        if ( range.size() < 2 ) throw Beam::bad_beam(Beam::InvalidParams);
        if ( (range[0] >= range[1]) || (range[0] <= 0.0) || (range[1] <= 0.0) ) {
            throw Beam::bad_beam(Beam::InvalidParams);
        }
    } else {
        if ( range[0] <= 0.0 ) throw Beam::bad_beam(Beam::InvalidParams);
    }

    Range = range;
}

Beam& Beam::operator = (const Beam &beam) // assignment operator
{

    Type = beam.Type;
    Shape = beam.Shape;
    Profile = beam.Profile;
    Range_Distr = beam.Range_Distr;

    Nrays = beam.Nrays;
    N_good_rays = beam.N_good_rays;
    N_lambda = beam.N_lambda;

    Params = beam.Params;
    Center = beam.Center;

    CompositeBeam = beam.CompositeBeam;

    destroy();

    X = rt_engine_alocate_vector(Nrays);
    if ( X == NULL ) throw Beam::bad_beam(Beam::MemAllocFailure);

    Y = rt_engine_alocate_vector(Nrays);
    if ( Y == NULL ) {
        destroy();
        throw Beam::bad_beam(Beam::MemAllocFailure);
    }

    Z = rt_engine_alocate_vector(Nrays);
    if ( Z == NULL ) {
        destroy();
        throw Beam::bad_beam(Beam::MemAllocFailure);
    }

    cX = rt_engine_alocate_vector(Nrays);
    if ( cX == NULL ) {
        destroy();
        throw Beam::bad_beam(Beam::MemAllocFailure);
    }

    cY = rt_engine_alocate_vector(Nrays);
    if ( cY == NULL ) {
        destroy();
        throw Beam::bad_beam(Beam::MemAllocFailure);
    }

    cZ = rt_engine_alocate_vector(Nrays);
    if ( cZ == NULL ) {
        destroy();
        throw Beam::bad_beam(Beam::MemAllocFailure);
    }

    lambda = rt_engine_alocate_vector(Nrays);
    if ( lambda == NULL ) {
        destroy();
        throw Beam::bad_beam(Beam::MemAllocFailure);
    }

    try {
        flag = new beam_flag_t[Nrays];
    } catch (bad_alloc &ex) {
        destroy();
        throw Beam::bad_beam(Beam::MemAllocFailure);
    }

    size_t N = Nrays*sizeof(real_t);

    memcpy(X,beam.X,N);
    memcpy(Y,beam.Y,N);
    memcpy(Z,beam.Z,N);

    memcpy(cX,beam.cX,N);
    memcpy(cY,beam.cY,N);
    memcpy(cZ,beam.cZ,N);

    memcpy(lambda,beam.lambda,N);

    N = Nrays*sizeof(beam_flag_t);
    memcpy(flag,beam.flag,N);

    return *this;
}


size_t Beam::GetTotalNumberOfRays() const
{
    return Nrays;
}

size_t Beam::GetNumberOfGoodRays() const
{
    return N_good_rays;
}

#define beam_N0 8
#define beam_dN 4

/*
    for shape == Beam::Circle:
        for profile == Beam::Concentric:

            params[0] - number of concentric circles
            params[1] - inner circle radius
            params[2] - outter circle radius

        for profile == Beam::Random:
            params[0] - number of rays
            params[1] - inner circle radius
            params[2] - outter circle radius

    for shape == Beam::Rectangle:
        for profile == Beam::Concentric:

            params[0] - number of concentric rectangles
            params[1] - inner rectangle width
            params[2] - outter rectangle width
            params[3] - ratio: width/height

        for profile == Beam::Random:
            params[0] - number of rays
            params[1] - inner rectangle width
            params[2] - outter rectangle width
            params[3] - ratio: width/height

        type = conic (shape = circle);
            params[3] - conic angle in degrees
        type = conic (shape = rectangle);
            params[4] - conic angle in degrees

        type = gauss (shape = circle);
            params[3] - Gaussian FWHM in arcsecs
        type = gauss (shape = rectangle);
            params[4] - Gaussian FWHM in arcsecs
*/

void Beam::Create(const vector<real_t> &params, const vector<real_t> &center, const vector<real_t> &range,
                  bool composite_flag, size_t Nlambda)
{
    size_t central_ray;
    real_t dR;
    RT_engine_error err;

                /* realization independed path of code ( independant on nVidia CUDA or Intel SIMD) */


    SetParams(params,composite_flag,Nlambda);
    SetCenter(center);
    SetRange(range);

    // be sure memory was not already allocated
    if ( Nrays ) destroy();

    // first, how many rays should be

    switch (Profile) {
        case Beam::Concentric: {
            // apply formula for sum of the N-first members of arithmetic progression
            // with a0 = N0, step = dN and params[0] as a number of terms
            size_t nn = static_cast<size_t>(Params[0]);
            central_ray = (Params[1] == 0) ? 1 : 0;
            Nrays = central_ray + (2*beam_N0 + (nn-1)*beam_dN)/2*nn;
            if ( central_ray ) {
                dR = Params[2]/Params[0]; // increment for circles/rectangle
            } else {
                dR = (Params[2]-Params[1])/(Params[0]-1);
            }
            break;
        }
        case Beam::Random: {
            Nrays = static_cast<size_t>(Params[0]);
            break;
        }
        default: {break;}
    }

    if ( !Nrays ) throw Beam::bad_beam(Beam::InvalidParams);

    if ( CompositeBeam ) Nrays *= N_lambda; // composite beam

    N_good_rays = Nrays; // init number of good rays to the total one

    // allocate memory


    X = rt_engine_alocate_vector(Nrays);
    if ( X == NULL ) throw Beam::bad_beam(Beam::MemAllocFailure);

    Y = rt_engine_alocate_vector(Nrays);
    if ( Y == NULL ) throw Beam::bad_beam(Beam::MemAllocFailure);

    Z = rt_engine_alocate_vector(Nrays);
    if ( Z == NULL ) throw Beam::bad_beam(Beam::MemAllocFailure);

    cX = rt_engine_alocate_vector(Nrays);
    if ( cX == NULL ) throw Beam::bad_beam(Beam::MemAllocFailure);

    cY = rt_engine_alocate_vector(Nrays);
    if ( cY == NULL ) throw Beam::bad_beam(Beam::MemAllocFailure);

    cZ = rt_engine_alocate_vector(Nrays);
    if ( cZ == NULL ) throw Beam::bad_beam(Beam::MemAllocFailure);

    lambda = rt_engine_alocate_vector(Nrays);
    if ( lambda == NULL ) throw Beam::bad_beam(Beam::MemAllocFailure);

    try {
        flag = new beam_flag_t[Nrays];
    } catch (bad_alloc &ex) {
        throw Beam::bad_beam(Beam::MemAllocFailure);
    }

    // init flag
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
#ifdef USING_LINUX
    for (size_t i = 0; i < Nrays; ++i ) flag[i] = 1;
#endif
#ifdef USING_MSVC // OpenMP v.2 does not support unsigned loop counter
    for (long long i = 0; i < Nrays; ++i ) flag[i] = 1;
#endif

    // compute coordinates of the beam
    switch ( Profile ) {
        case Beam::Random: {
            err = beam_random(Shape,Params.data(),Center.data(),Nrays,X,Y,Z);
            if ( err != ENGINE_ERROR_OK ) {
                cerr << "ERROR: " << err << endl;
                throw Beam::bad_beam(Beam::CreationFailure);
                return;
            }
            break;
        }
        case Beam::Concentric: {
            break;
        }
    }

    // compute cosins of the beam

    switch ( Type ) {
        case Beam::Parallel: {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
#ifdef USING_LINUX
        for (size_t i = 0; i < Nrays; ++i ) {
#endif
#ifdef USING_MSVC // OpenMP v.2 does not support unsigned loop counter
        for (long long i = 0; i < Nrays; ++i ) {
#endif
                cX[i] = 0;
                cY[i] = 0;
                cZ[i] = 1;
            }
            break;
        }
        case Beam::Conic: { // use of external engine function
            real_t conic_ang;
            switch ( Shape ) {
                case Beam::Circle: {
                    conic_ang = params[3]*deg2rad;
                    break;
                }
                case Beam::Rectangle: {
                    conic_ang = params[4]*deg2rad;
                    break;
                }
                default: throw Beam::bad_beam(Beam::InvalidParams);
            }
            err = beam_conic_cosins(Nrays,conic_ang,X,Y,cX,cY,cZ);
            if ( err != ENGINE_ERROR_OK ) {
                cerr << "ERROR: " << err << endl;
                throw Beam::bad_beam(Beam::CreationFailure);
                return;
            }
            break;
        }
        case Beam::Gauss: {
            real_t fwhm;
            switch ( Shape ) {
                case Beam::Circle: {
                    fwhm = params[3]; // fwhm must be given in arcsecs
                    break;
                }
                case Beam::Rectangle: {
                    fwhm = params[4]; // fwhm must be given in arcsecs
                    break;
                }
                default: throw Beam::bad_beam(Beam::InvalidParams);
            }
            err = beam_gauss_cosins(Nrays,fwhm,cX,cY,cZ);
            if ( err != ENGINE_ERROR_OK ) {
                cerr << "ERROR: " << err << endl;
                throw Beam::bad_beam(Beam::CreationFailure);
                return;
            }
            break;
        }
    }

    // compute wavelength vector of the beam

    switch ( Range_Distr ) {
        case Beam::Monochromatic: { // monochromatic beam
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
#ifdef USING_LINUX
        for (size_t i = 0; i < Nrays; ++i ) {
#endif
#ifdef USING_MSVC // OpenMP v.2 does not support unsigned loop counter
        for (long long i = 0; i < Nrays; ++i ) {
#endif
                lambda[i] = (range[1] + range[0])/2.0;
//                cout << "lambda[" << i << "] = " << lambda[i] << endl;
            }
            break;
        }
        case Beam::UniformDistr: { // uniformly distributed wavelength in increasing order
            if ( CompositeBeam ) {
                err = uniform_range_beam(N_lambda,Range.data(),lambda);
            } else {
                err = uniform_range_beam(Nrays,Range.data(),lambda);
            }
            if ( err != ENGINE_ERROR_OK ) {
//                cerr << "ERROR: " << err << endl;
                throw Beam::bad_beam(Beam::CreationFailure);
                return;
            }
            break;
        }
        case Beam::RandomDistr: { // randomly distributed wavelengths
            if ( CompositeBeam ) {
                err = random_range_beam(N_lambda,Range.data(),lambda); // generate wavelengths
            } else {
                err = random_range_beam(Nrays,Range.data(),lambda);
            }
            if ( err != ENGINE_ERROR_OK ) {
//                cerr << "ERROR: " << err << endl;
                throw Beam::bad_beam(Beam::CreationFailure);
                return;
            }
            break;
        }
        default: break;
    }

    if ( CompositeBeam ) { // "clone" lambda vector
        size_t N = static_cast<size_t>(Params[0]); // number of rays per each wavelength
        size_t Nbytes = N_lambda*sizeof(real_t);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
#ifdef USING_LINUX
        for (size_t i = 1; i < N; ++i ) {
#endif
#ifdef USING_MSVC // OpenMP v.2 does not support unsigned loop counter
        for (long long i = 1; i < N; ++i ) {
#endif
            memcpy(lambda+i*N_lambda,lambda,Nbytes);
        }
    }
}


void Beam::Recreate() // just re-create beam with already defined parameters
{
    Create(Params,Center,Range,CompositeBeam,N_lambda);
}

vector<real_t> Beam::GetBeamDim()
{
    vector<real_t> dim;

    dim.push_back(*min_element(X,X+N_good_rays-1));
    dim.push_back(*max_element(X,X+N_good_rays-1));

    dim.push_back(*min_element(Y,Y+N_good_rays-1));
    dim.push_back(*max_element(Y,Y+N_good_rays-1));

    dim.push_back(*min_element(Z,Z+N_good_rays-1));
    dim.push_back(*max_element(Z,Z+N_good_rays-1));

    return dim;
}

vector<real_t> Beam::GetBeamCosins()
{
    vector<real_t> dim;

//    dim.push_back(*min_element(cX,cX+N_good_rays-1));
//    dim.push_back(*max_element(cX,cX+N_good_rays-1));

//    dim.push_back(*min_element(cY,cY+N_good_rays-1));
//    dim.push_back(*max_element(cY,cY+N_good_rays-1));

    dim.push_back(*min_element(cZ,cZ+N_good_rays-1));
    dim.push_back(*max_element(cZ,cZ+N_good_rays-1));

    return dim;
}


// normalize diraction cosins to unity
void Beam::Normalize()
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
#ifdef USING_LINUX
        for (size_t i = 0; i < Nrays; ++i ) {
#endif
#ifdef USING_MSVC // OpenMP v.2 does not support unsigned loop counter
        for (long long i = 0; i < Nrays; ++i ) {
#endif
//            real_t acX = abs(cX[i]);
//            real_t acY = abs(cY[i]);
//            real_t acZ = abs(cZ[i]);
            real_t acX = fabs(cX[i]);
            real_t acY = fabs(cY[i]);
            real_t acZ = fabs(cZ[i]);

            real_t norm;


            if ( acX >= acY ) {
                if ( acX >= acZ ) norm = acX*sqrt(1.0 + acY*acY/acX/acX + acZ*acZ/acX/acX);
                else norm = acZ*sqrt(1.0 + acX*acX/acZ/acZ + acY*acY/acZ/acZ);
            } else {
                if ( acY >= acZ ) norm = acY*sqrt(1.0 + acX*acX/acY/acY + acZ*acZ/acY/acY);
                else norm = acZ*sqrt(1.0 + acX*acX/acZ/acZ + acY*acY/acZ/acZ);
            }

//            if ( (acX >= acY) && (acX >= acZ) ) {
//                norm = acX*sqrt(1.0 + acY*acY/acX/acX + acZ*acZ/acX/acX);
//            }
//            if ( (acY >= acX) && (acY >= acZ) ) {
//                norm = acY*sqrt(1.0 + acX*acX/acY/acY + acZ*acZ/acY/acY);
//            }
//            if ( (acZ >= acX) && (acZ >= acY) ) {
//                norm = acZ*sqrt(1.0 + acX*acX/acZ/acZ + acY*acY/acZ/acZ);
//            }
//            real_t norm = sqrt(cX[i]*cX[i] + cY[i]*cY[i] + cZ[i]*cZ[i]);

            cX[i] /= norm;
            cY[i] /= norm;
            cZ[i] /= norm;
        }
}

// rearrange coordinates, cosins and wavelength vectors in such way that bad rays are always in the end of the vectors
void Beam::Rearrange()
{
    beam_flag_t *start = flag;
    beam_flag_t *end = flag+(N_good_rays-1);
    size_t i;

//    cout << "rearr: N_good_before = " << N_good_rays << endl;
    i = 0;
    while (start <= end) {
        while (*start == 1 && start <= end) {  // Find a Zero, starting from the front
            ++start;
            ++i;
        }
        while (*end == 0 && start <= end) {  // Find a One, starting from the back
            --end;
            --N_good_rays;
        }
        if (start < end) {  // *start == Zero, and *end == One, and start is to the left of end
            *start = 1;
            *end = 0;
            swap_elements(i,N_good_rays-1);
        }
    }
//    cout << "rearr: N_good_after = " << N_good_rays << " (i = " << i << ")" << endl;
}


void Beam::Transform(const Surface &surf)
{
    vector<real_t> ang, dist;

    ang = surf.GetAngles();
    dist = surf.GetDistance();

    // transform to radians because of GetAngles return angles in degrees (see surface.cpp)
    for ( int i = 0; i < 3; ++i ) {
        ang[i] *= deg2rad;
    }

    RT_engine_error err = beam_transform(N_good_rays,ang.data(),dist.data(),
                                         X,Y,Z,cX,cY,cZ);
    if ( err != ENGINE_ERROR_OK ) throw Beam::bad_beam(Beam::TransformationFailure);
}

// angles must be in degrees
void Beam::Transform(const vector<real_t> &dist, const vector<real_t> &angs)
{
    if ( (dist.size() < 3) || (angs.size() < 3) ) throw Beam::bad_beam(Beam::InvalidParams);

    RT_engine_error err;
    vector<real_t> v = angs;

    vector<real_t> dd = dist;

    for ( int i = 0; i < 3; ++i ) v[i] *= deg2rad;

    err = beam_transform(N_good_rays,v.data(),dd.data(),X,Y,Z,cX,cY,cZ);

    if ( err != ENGINE_ERROR_OK ) throw Beam::bad_beam(Beam::TransformationFailure);
}


void Beam::Save(const string &filename)
{
    FILE *file;
    size_t count;

    file = fopen(filename.data(),"wb");

    if ( file == NULL ) {
        throw Beam::bad_beam(Beam::SaveFailure);
    }

    try {
        // first, write total and good-ray number of rays
        count = fwrite(&Nrays,sizeof(size_t),1,file);
        if ( count != 1 ) throw Beam::bad_beam(Beam::SaveFailure);

        count = fwrite(&N_good_rays,sizeof(size_t),1,file);
        if ( count != 1 ) throw Beam::bad_beam(Beam::SaveFailure);

        count = fwrite(&N_lambda,sizeof(size_t),1,file);
        if ( count != 1 ) throw Beam::bad_beam(Beam::SaveFailure);

        // first, save "good" rays

        count = fwrite(X,sizeof(real_t),N_good_rays,file);
        if ( count != N_good_rays ) throw Beam::bad_beam(Beam::SaveFailure);

        count = fwrite(Y,sizeof(real_t),N_good_rays,file);
        if ( count != N_good_rays ) throw Beam::bad_beam(Beam::SaveFailure);

        count = fwrite(Z,sizeof(real_t),N_good_rays,file);
        if ( count != N_good_rays ) throw Beam::bad_beam(Beam::SaveFailure);

        count = fwrite(cX,sizeof(real_t),N_good_rays,file);
        if ( count != N_good_rays ) throw Beam::bad_beam(Beam::SaveFailure);

        count = fwrite(cY,sizeof(real_t),N_good_rays,file);
        if ( count != N_good_rays ) throw Beam::bad_beam(Beam::SaveFailure);

        count = fwrite(cZ,sizeof(real_t),N_good_rays,file);
        if ( count != N_good_rays ) throw Beam::bad_beam(Beam::SaveFailure);

        count = fwrite(lambda,sizeof(real_t),N_good_rays,file);
        if ( count != N_good_rays ) throw Beam::bad_beam(Beam::SaveFailure);

        // save remained rays
        size_t rest_Nrays = Nrays-N_good_rays;

        count = fwrite(X+N_good_rays,sizeof(real_t),rest_Nrays,file);
        if ( count != rest_Nrays ) throw Beam::bad_beam(Beam::SaveFailure);

        count = fwrite(Y+N_good_rays,sizeof(real_t),rest_Nrays,file);
        if ( count != rest_Nrays ) throw Beam::bad_beam(Beam::SaveFailure);

        count = fwrite(Z+N_good_rays,sizeof(real_t),rest_Nrays,file);
        if ( count != rest_Nrays ) throw Beam::bad_beam(Beam::SaveFailure);

        count = fwrite(cX+N_good_rays,sizeof(real_t),rest_Nrays,file);
        if ( count != rest_Nrays ) throw Beam::bad_beam(Beam::SaveFailure);

        count = fwrite(cY+N_good_rays,sizeof(real_t),rest_Nrays,file);
        if ( count != rest_Nrays ) throw Beam::bad_beam(Beam::SaveFailure);

        count = fwrite(cZ+N_good_rays,sizeof(real_t),rest_Nrays,file);
        if ( count != rest_Nrays ) throw Beam::bad_beam(Beam::SaveFailure);

        count = fwrite(lambda+N_good_rays,sizeof(real_t),rest_Nrays,file);
        if ( count != rest_Nrays ) throw Beam::bad_beam(Beam::SaveFailure);

    } catch (Beam::bad_beam &ex) {
        fclose(file);
        throw;
    }

    fclose(file);
}


            /*  Beam class private methods  */

void Beam::destroy()
{
    if ( X != NULL ) {
        rt_engine_free_vector(X);
        X = NULL;
    }

    if ( Y != NULL ) {
        rt_engine_free_vector(Y);
        Y = NULL;
    }

    if ( Z != NULL ) {
        rt_engine_free_vector(Z);
        Z = NULL;
    }

    if ( cX != NULL ) {
        rt_engine_free_vector(cX);
        cX = NULL;
    }

    if ( cY != NULL ) {
        rt_engine_free_vector(cY);
        cY = NULL;
    }

    if ( cZ != NULL ) {
        rt_engine_free_vector(cZ);
        cZ = NULL;
    }

    if ( lambda != NULL ) {
        rt_engine_free_vector(lambda);
        lambda = NULL;
    }

    if ( flag != NULL ) {
        delete[] flag; // flag array was allocated in Beam::Create not by RT engine function!
        flag = NULL;
    }
}


void Beam::swap_elements(size_t i1, size_t i2)
{
//    real_t tmp;

#ifdef USE_OPENMP
#pragma omp sections
#endif
    {
    {
    real_t tmp = X[i1];
    X[i1] = X[i2];
    X[i2] = tmp;
    }
#ifdef USE_OPENMP
#pragma omp section
#endif
    {
    real_t tmp = Y[i1];
    Y[i1] = Y[i2];
    Y[i2] = tmp;
    }
#ifdef USE_OPENMP
#pragma omp section
#endif
    {
    real_t tmp = Z[i1];
    Z[i1] = Z[i2];
    Z[i2] = tmp;
    }
#ifdef USE_OPENMP
#pragma omp section
#endif
    {
    real_t tmp = cX[i1];
    cX[i1] = cX[i2];
    cX[i2] = tmp;
    }
#ifdef USE_OPENMP
#pragma omp section
#endif
    {
    real_t tmp = cY[i1];
    cY[i1] = cY[i2];
    cY[i2] = tmp;
    }
#ifdef USE_OPENMP
#pragma omp section
#endif
    {
    real_t tmp = cZ[i1];
    cZ[i1] = cZ[i2];
    cZ[i2] = tmp;
    }
#ifdef USE_OPENMP
#pragma omp section
#endif
    {
    real_t tmp = lambda[i1];
    lambda[i1] = lambda[i2];
    lambda[i2] = tmp;
    }
    }
}


int load_beam_data(const char *rt_filename, size_t *total_Nrays, size_t *N_good, size_t *N_lambda,
                   real_t **X, real_t **Y, real_t **Z,
                   real_t **cX, real_t **cY, real_t **cZ,
                   real_t **lambda)
{
    FILE *fid;
    size_t n, vec_len;

    real_t *ptr;

    fid = fopen(rt_filename,"rb");

    try {
        if ( fid == NULL ) {
            throw 1;
        }

        n = fread(total_Nrays,sizeof(size_t),1,fid);
        if ( n != 1 ) throw 2;

        n = fread(N_good,sizeof(size_t),1,fid);
        if ( n != 1 ) throw 2;

        vec_len = *N_good*sizeof(real_t);

        n = fread(N_lambda,sizeof(size_t),1,fid);
        if ( n != 1 ) throw 2;

        // allocate arrays
        if ( X != NULL ) {
            ptr = (real_t*) malloc(vec_len);
            if ( ptr == NULL ) throw 3;
            *X = ptr;
        }

        if ( Y != NULL ) {
            ptr = (real_t*) malloc(vec_len);
            if ( ptr == NULL ) throw 3;
            *Y = ptr;
        }

        if ( Z != NULL ) {
            ptr = (real_t*) malloc(vec_len);
            if ( ptr == NULL ) throw 3;
            *Z = ptr;
        }

        if ( cX != NULL ) {
            ptr = (real_t*) malloc(vec_len);
            if ( ptr == NULL ) throw 3;
            *cX = ptr;
        }

        if ( cY != NULL ) {
            ptr = (real_t*) malloc(vec_len);
            if ( ptr == NULL ) throw 3;
            *cY = ptr;
        }

        if ( cZ != NULL ) {
            ptr = (real_t*) malloc(vec_len);
            if ( ptr == NULL ) throw 3;
            *cZ = ptr;
        }

        if ( lambda != NULL ) {
            ptr = (real_t*) malloc(vec_len);
            if ( ptr == NULL ) throw 3;
            *lambda = ptr;
        }

        // read the data

        if ( X != NULL ) {
            n = fread(*X,sizeof(real_t),*N_good,fid);
            if ( n != *N_good ) throw 2;
        } else {
            n = fseek(fid,vec_len,SEEK_CUR);
            if ( n ) throw 2;
        }

        if ( Y != NULL ) {
            n = fread(*Y,sizeof(real_t),*N_good,fid);
            if ( n != *N_good ) throw 2;
        } else {
            n = fseek(fid,vec_len,SEEK_CUR);
            if ( n ) throw 2;
        }

        if ( Z != NULL ) {
            n = fread(*Z,sizeof(real_t),*N_good,fid);
            if ( n != *N_good ) throw 2;
        } else {
            n = fseek(fid,vec_len,SEEK_CUR);
            if ( n ) throw 2;
        }

        if ( cX != NULL ) {
            n = fread(*cX,sizeof(real_t),*N_good,fid);
            if ( n != *N_good ) throw 2;
        } else {
            n = fseek(fid,vec_len,SEEK_CUR);
            if ( n ) throw 2;
        }

        if ( cY != NULL ) {
            n = fread(*cY,sizeof(real_t),*N_good,fid);
            if ( n != *N_good ) throw 2;
        } else {
            n = fseek(fid,vec_len,SEEK_CUR);
            if ( n ) throw 2;
        }

        if ( cZ != NULL ) {
            n = fread(*cZ,sizeof(real_t),*N_good,fid);
            if ( n != *N_good ) throw 2;
        } else {
            n = fseek(fid,vec_len,SEEK_CUR);
            if ( n ) throw 2;
        }

        if ( lambda != NULL ) {
            n = fread(*lambda,sizeof(real_t),*N_good,fid);
            if ( n != *N_good ) throw 2;
        } else {
            n = fseek(fid,vec_len,SEEK_CUR);
            if ( n ) throw 2;
        }

    } catch (int i) {
        fclose(fid);
        return i;
    }

    fclose(fid);

    return 0;
}
