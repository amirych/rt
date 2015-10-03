#include "surface.h"

#include "rt_engine_func.h"
#include "rt_engine_errors.h"
#include "tabulated_function.h"

#include<iostream>
#include<fstream>



             /**************************************
             *                                     *
             *  Surface base class implementation  *
             *                                     *
             **************************************/

/*  Scheme class exceptions implementation */

Surface::bad_surface::bad_surface(SurfaceError err): exception(), Error(err)
{

}

Surface::SurfaceError Surface::bad_surface::surface_error() const
{
    return Error;
}


            /*  Base class constructor and destructor  */

Surface::Surface(const SurfaceClass sclass, const SurfaceType type, const SurfaceShape shape):
    Class(sclass), Type(type), Shape(shape),
    Params(vector<real_t>(3,0.0)), Size(vector<real_t>(2,0.0)), Center(vector<real_t>(2,0.0)),
    Distance(vector<real_t>(3,0.0)), Angles(vector<real_t>(3,0.0)),
    CurrentError(Surface::Ok), Comment(""), QE(1.0) // initialize QE by 100%
{
//    cerr << "Base class Surface is created!\n";
}

Surface::~Surface()
{
}


            /*  Base class public methods   */

/*
    Surface classes:
        CONIC:
                                                              N
            F(X,Y,Z) = Z - c*r^2/(1+SQRT(1-(1+k)*c^2*r^2)) - SUM(a_i*r^(2*i)), where r^2 = X^2+Y^2, and a_i are known coefficients
                                                             i=1

        TORIC:

            F(X,Y,Z) = Z - f(Y) - 0.5*cR*(X^2+Z^2-f(Y)^2),
                                                                   N
                       f(Y) = c*Y^2/[1 + SQRT(1-(k+1)*c^2*Y^2)] + SUM(a_i*Y^(2*i)), where a_i are known coefficients
                                                                  i=1


    For Class == Conic:
        pars = [c,k,N,...], where c = 1/R and R is vertex curvature radius (positive for convex surface and negative for concave one)
                                  k = -e^2 and e is conic surface eccentricity:
                                      k < -1 for hyperbola
                                      k = -1 for parabola
                                      -1 < k < 0 for ellipsoid
                                      k = 0 for spheroid
                                      k > 0 for oblate spheroid
                                  N - a number of coefficients in aspherical term
              SIZE(pars) = (3 + N) elements

    For Class == Toric:
        pars = [c,k,N,....,cR], where c = 1/R and R is vertex curvature radius (positive for convex surface and negative for concave one)
                                                  of function f(Y)
                                      k = -e^2 and e is conic surface eccentricity for function f(Y)
                                      N - a number of coefficients in aspherical term in function f(Y)
                                      cR = 1/Rr and Rr is distance from origin to revolution axis (positive for convex surface and negative for concave one)
              SIZE(pars) = (4 + N) elements

*/

void Surface::SetParams(const vector<real_t> &pars)
{
    CurrentError = Surface::Ok;

    if ( (pars.size() < 2) || pars.empty() ) {
        CurrentError = Surface::InvalidSurfaceParam;
        throw bad_surface(Surface::InvalidSurfaceParam);
        return;
    }

    switch ( Class ) {
        case Surface::Toric: {
            if ( pars.size() < 4 ) {
                CurrentError = Surface::InvalidSurfaceParam;
                throw bad_surface(Surface::InvalidSurfaceParam);
                return;
            }

            size_t N = static_cast<size_t>(pars[2]);

            if ( pars.size() < (4+N) ) { // c, k, N, N polynominal terms, cR
                CurrentError = Surface::InvalidSurfaceParam;
                throw bad_surface(Surface::InvalidSurfaceParam);
                return;
            }
            Params = pars;
            return;
        }
        case Surface::Conic: {
            if ( pars.size() < 3 ) {
                CurrentError = Surface::InvalidSurfaceParam;
                throw bad_surface(Surface::InvalidSurfaceParam);
                return;
            }

            size_t N = static_cast<size_t>(pars[2]);

            if ( pars.size() < (3+N) ) { // c, k, N, N polynominal terms
                CurrentError = Surface::InvalidSurfaceParam;
                throw bad_surface(Surface::InvalidSurfaceParam);
                return;
            }
            Params = pars;
            return;
        }
        case Surface::Plane: { // just ignore
            return;
        }
        default: return;
    }
}


void Surface::SetGeometry(const vector<real_t> &size, const vector<real_t> &center)
{
    CurrentError = Surface::Ok;

    switch ( Shape ) {
        case Surface::Circle: {
            if ( size.empty() ) {
                CurrentError = Surface::InvalidSurfaceSize;
                throw bad_surface(Surface::InvalidSurfaceSize);
                return;
            }
            Size[0] = size[0];
            break;
        }
        case Surface::Rectangle: {
            if ( size.empty() || (size.size() < 2) ) {
                CurrentError = Surface::InvalidSurfaceSize;
                throw bad_surface(Surface::InvalidSurfaceSize);
                return;
            }
            Size[0] = size[0];
            Size[1] = size[1];
            break;
        }
        default: break;
    }

    if ( (center.size() < 2) || center.empty() ) {
        CurrentError = Surface::InvalidSurfaceCenter;
        throw bad_surface(Surface::InvalidSurfaceCenter);
        return;
    }
    Center[0] = center[0];
    Center[1] = center[1];

//    cout << "\nSetGeometry: " << Center[0] << " " << Center[1] << endl;
}


void Surface::SetDistance(const vector<real_t> &dist)
{
    CurrentError = Surface::Ok;

    if ( dist.empty() || (dist.size() < 3) ) {
        CurrentError = Surface::InvalidSurfaceDistance;
        throw bad_surface(Surface::InvalidSurfaceDistance);
        return;
    }
    Distance.assign(dist.begin(),dist.end());
}


void Surface::SetAngles(const vector<real_t> &angs) // input angles values are in degrees (to be transformed in the function to radians)
{
    CurrentError = Surface::Ok;

    if ( angs.empty() || (angs.size() < 3) ) {
        CurrentError = Surface::InvalidSurfaceAngles;
        throw bad_surface(Surface::InvalidSurfaceAngles);
        return;
    }
    Angles.assign(angs.begin(),angs.end());
    // convert to radians
    for ( int i = 0; i < 3; ++i ) Angles[i] *= deg2rad;
}


void Surface::SetQE(const real_t qe_val)
{
    QE.Set(qe_val);
}


void Surface::SetQE(const char *qe_file)
{
    QE.Load(qe_file);
}


vector<real_t> Surface::GetDistance() const
{
    return Distance;
}

vector<real_t> Surface::GetAngles() const
{
    vector<real_t> v = Angles;

    // convert back to degrees
    for ( int i = 0; i < 3; ++i ) v[i] /= deg2rad;

    return v;
}


void Surface::SetComment(const string &comm)
{
    Comment = comm;
}


string Surface::GetComment() const
{
    return Comment;
}

void Surface::Intersection(Beam &beam)
{
    RT_engine_error err;
    size_t N_bad = 0; // it will contain a number of rays with no surface-and-beam intersection

    err = surface_intersection(Class,Params.data(),beam.N_good_rays,beam.X,beam.Y,beam.Z,
                               beam.cX,beam.cY,beam.cZ,&N_bad,beam.flag);

    if ( err != ENGINE_ERROR_OK ) throw bad_surface(Surface::IntersectionFailure);

    // rearrange coordinates and cosins vectors if there are new bad rays
    if ( N_bad ) {
        if ( N_bad == beam.N_good_rays ) throw bad_surface(Surface::NoIntersections);
//        cout << "N_bad = " << N_bad << "; N_good = " << beam.N_good_rays << endl;
        beam.Rearrange();
    }

}


void Surface::ApplyConstrains(Beam &beam)
{
    // apply size constrains
    real_t R2;
//    , ray_R2, xr, yr;

    R2 = Size[0]*Size[0]; // for circular shape

    size_t N_bad = 0;

//    cout << "\nCENTER: [" << Center[0] << ", " << Center[1] << "]\n";

    // NO PARALLEL ALGORITHM HERE!!!!!
#ifdef USE_OPENMP
//#pragma omp parallel for reduction(+:N_bad)
#pragma omp parallel for
#endif
#ifdef USING_LINUX
    for (size_t i = 0; i < beam.N_good_rays; ++i ) {
#endif
#ifdef USING_MSVC // OpenMP v.2 does not suppport unsigned loop counter
    for (long long i = 0; i < beam.N_good_rays; ++i ) {
#endif
        real_t xr = beam.X[i] - Center[0];
        real_t yr = beam.Y[i] - Center[1];
        switch ( Shape ) {
            case Surface::Circle: {
                real_t ray_R2 = xr*xr + yr*yr;
                if ( ray_R2 > R2 ) {
                    beam.flag[i] = 0;
                    ++N_bad;
                }
                break;
            }
            case Surface::Rectangle: {
//            cout << "SIZE: [" << Size[0] << ", " << Size[1] << "]\n";
//                if ( (xr < -Size[0]/2.0) || (xr > Size[0]/2.0) ||
//                     (yr < -Size[1]/2.0) || (yr > Size[1]/2.0) ) {
                if ( (beam.X[i] < (Center[0] -Size[0]/2.0)) || (beam.X[i] > (Center[0] +Size[0]/2.0)) ||
                     (beam.Y[i] < (Center[1] -Size[1]/2.0)) || (beam.Y[i] > (Center[1] + Size[1]/2.0)) ) {
//                    cout << "left edge: " << Center[0] -Size[0]/2.0 << "; right edge: " << Center[0] +Size[0]/2.0 << endl;
//                    cout << "bottom edge: " << Center[1] -Size[1]/2.0 << "; top edge: " << Center[1] +Size[1]/2.0 << endl;
//                    cout << " beam: [" << beam.X[i] << ", " << beam.Y[i] << "]\n";
                    beam.flag[i] = 0;
                    ++N_bad;
                }
                break;
            }
        }
    }

    // rearrange coordinates and cosins vectors if there are new bad rays
    if ( N_bad ) {
//        cout << " BAD RAYS: " << N_bad << endl;
        if ( N_bad == beam.N_good_rays ) throw bad_surface(Surface::NoIntersections);
        beam.Rearrange();
    }
}


void Surface::Action(Beam &beam)
{
    // nothing to do
//    cerr << "Base Surface class has no action!!!\n";
}


void Surface::ApplyQE(vector<real_t> &lambda, vector<real_t> &spec)
{
    if ( lambda.size() != spec.size() ) throw Surface::bad_surface(Surface::BadQE);

    RT_engine_error err = surface_QE(spec.size(),lambda.data(),spec.data(),QE);

    if ( err != ENGINE_ERROR_OK ) {
        throw Surface::bad_surface(Surface::BadQE);
    }
}


Surface::SurfaceClass Surface::GetClass() const
{
    return Class;
}

Surface::SurfaceType Surface::GetType() const
{
    return Type;
}

Surface::SurfaceShape Surface::GetShape() const
{
    return Shape;
}


Surface::SurfaceError Surface::GetError() const {
    return CurrentError;
}



        /*************************************
        *                                    *
        *  Inherited classes implementation  *
        *                                    *
        *************************************/


Plane::Plane(const SurfaceType type, const SurfaceShape shape): Surface(Surface::Plane,type,shape)
{
//    cerr << "Plane surface is created!\n";
}


                    /*    MIRROR CLASS    */

Mirror::Mirror(const SurfaceClass sclass, const SurfaceShape shape): Surface(sclass,Surface::Mirror,shape)
{
//    cerr << "Mirror is created!\n";
}

void Mirror::Action(Beam &beam)
{
    RT_engine_error err;

    if ( beam.N_good_rays == 0 ) throw bad_surface(Surface::NoGoodRays);

    err = surface_reflection(Class,Params.data(),beam.N_good_rays,
                             beam.X,beam.Y,beam.Z,beam.cX,beam.cY,beam.cZ);
    if ( err != ENGINE_ERROR_OK ) {
        throw bad_surface(Surface::ActionFailure);
    }
}

                    /*    LENS CLASS    */

Lens::bad_lens::bad_lens(Lens::LensError err)
{
    Error = err;
}

Lens::LensError Lens::bad_lens::bad_lens_error() const
{
    return Error;
}

Lens::Lens(const SurfaceClass sclass, const SurfaceShape shape, const real_t &ref_index1, const real_t &ref_index2):
    Surface(sclass,Surface::Lens,shape), n1(ref_index1), n2(ref_index2)
{
    if ( (ref_index1 <= 0.0) || (ref_index2 <= 0.0) ) {
        throw Lens::bad_lens(Lens::BadRefIndex);
    }
//    cerr << "Lens is created!\n";
}

Lens::Lens(const SurfaceClass sclass, const SurfaceShape shape, const char* ref_index1, const char* ref_index2)
    try :Surface(sclass,Surface::Lens,shape), n1(ref_index1), n2(ref_index2)
    {
//        cerr << "Lens is created!\n";
    } catch (TabulatedFunction::bad_tab_func &ex) {
        throw Lens::bad_lens(Lens::BadRefIndexFile);
    }

Lens::Lens(const SurfaceClass sclass, const SurfaceShape shape, const char* ref_index1, const real_t &ref_index2)
    try :Surface(sclass,Surface::Lens,shape), n1(ref_index1), n2(ref_index2)
{
    if ( ref_index2 <= 0.0 ) {
        throw Lens::bad_lens(Lens::BadRefIndex);
    }
//    cerr << "Lens is created!\n";
} catch (TabulatedFunction::bad_tab_func &ex) {
    throw Lens::bad_lens(Lens::BadRefIndexFile);
} catch (...) {
    throw;
}

Lens::Lens(const SurfaceClass sclass, const SurfaceShape shape, const real_t &ref_index1, const char* ref_index2)
    try :Surface(sclass,Surface::Lens,shape), n1(ref_index1), n2(ref_index2)
{
    if ( ref_index1 <= 0.0 ) {
        throw Lens::bad_lens(Lens::BadRefIndex);
    }
//    cerr << "Lens is created!\n";
} catch (TabulatedFunction::bad_tab_func &ex) {
    throw Lens::bad_lens(Lens::BadRefIndexFile);
} catch (...) {
    throw;
}

void Lens::Action(Beam &beam)
{
    RT_engine_error err;
    size_t N_bad = 0;


    if ( beam.N_good_rays == 0 ) throw bad_surface(Surface::NoGoodRays);

    err = surface_refraction(Class,Params.data(),beam.N_good_rays,beam.X,beam.Y,beam.Z,
                             beam.cX,beam.cY,beam.cZ,beam.lambda,n1,n2,&N_bad,beam.flag);

    if ( err != ENGINE_ERROR_OK ) {
        throw bad_surface(Surface::ActionFailure);
    }

    if ( N_bad ) {
        if ( N_bad == beam.N_good_rays ) throw bad_surface(Surface::NoGoodRays);
        beam.Rearrange();
    }
}


                    /*    GRATING CLASS    */

Grating::bad_grating::bad_grating(Grating::GratingError err)
{
    Error = err;
}

Grating::GratingError Grating::bad_grating::bad_grating_error() const
{
    return Error;
}

Grating::Grating(const SurfaceClass sclass, const SurfaceShape shape, const Grating::GratingType gr_type, const vector<long> &ord,
                 const real_t &gr_const, const real_t &ref_index1, const real_t &ref_index2)
try :Surface(sclass,Surface::Grating,shape), Grating_Type(gr_type), Order(ord), Grating_constant(gr_const), n1(ref_index1), n2(ref_index2)
{
    if ( ord.size() == 0 ) throw Grating::bad_grating(Grating::BadOrder);

    if ( ord.size() > 1 ) { // range of orders
        if ( ord[0] > ord[1] ) throw Grating::bad_grating(Grating::BadOrder);
    }

    if ( (ref_index1 <= 0.0) || (ref_index2 <= 0.0) ) {
        throw Grating::bad_grating(Grating::BadRefIndex);
    }
//    cerr << "Grating is created!\n";
} catch (...) {
    throw;
}

Grating::Grating(const SurfaceClass sclass, const SurfaceShape shape, const Grating::GratingType gr_type, const vector<long> &ord,
                 const real_t &gr_const, const char* ref_index1, const char* ref_index2)
    try :Surface(sclass,Surface::Grating,shape), Grating_Type(gr_type), Order(ord), Grating_constant(gr_const), n1(ref_index1), n2(ref_index2)
{
    if ( ord.size() == 0 ) throw Grating::bad_grating(Grating::BadOrder);

    if ( ord.size() > 1 ) { // range of orders
        if ( ord[0] > ord[1] ) throw Grating::bad_grating(Grating::BadOrder);
    }

//    cerr << "Grating is created!\n";
} catch (TabulatedFunction::bad_tab_func &ex) {
    throw Grating::bad_grating(Grating::BadRefIndexFile);
} catch (...) {
    throw;
}

Grating::Grating(const SurfaceClass sclass, const SurfaceShape shape, const Grating::GratingType gr_type, const vector<long> &ord,
                 const real_t &gr_const, const real_t &ref_index1, const char* ref_index2)
    try :Surface(sclass,Surface::Grating,shape), Grating_Type(gr_type), Order(ord), Grating_constant(gr_const), n1(ref_index1), n2(ref_index2)
{
    if ( ord.size() == 0 ) throw Grating::bad_grating(Grating::BadOrder);

    if ( ord.size() > 1 ) { // range of orders
        if ( ord[0] > ord[1] ) throw Grating::bad_grating(Grating::BadOrder);
    }

    if ( ref_index1 <= 0.0 ) {
        throw Grating::bad_grating(Grating::BadRefIndex);
    }
//    cerr << "Grating is created!\n";
} catch (TabulatedFunction::bad_tab_func &ex) {
    throw Grating::bad_grating(Grating::BadRefIndexFile);
} catch (...) {
    throw;
}

Grating::Grating(const SurfaceClass sclass, const SurfaceShape shape, const Grating::GratingType gr_type, const vector<long> &ord,
                 const real_t &gr_const, const char* ref_index1, const real_t &ref_index2)
    try :Surface(sclass,Surface::Grating,shape), Grating_Type(gr_type), Order(ord), Grating_constant(gr_const), n1(ref_index1), n2(ref_index2)
{
    if ( ord.size() == 0 ) throw Grating::bad_grating(Grating::BadOrder);

    if ( ord.size() > 1 ) { // range of orders
        if ( ord[0] > ord[1] ) throw Grating::bad_grating(Grating::BadOrder);
    }

    if ( ref_index2 <= 0.0 ) {
        throw Grating::bad_grating(Grating::BadRefIndex);
    }
//    cerr << "Grating is created!\n";
} catch (TabulatedFunction::bad_tab_func &ex) {
    throw Grating::bad_grating(Grating::BadRefIndexFile);
} catch (...) {
    throw;
}


void Grating::Action(Beam &beam)
{
    RT_engine_error err;
    size_t N_bad = 0;

    if ( beam.N_good_rays == 0 ) throw bad_surface(Surface::NoGoodRays);

    long ord = Order[0];

    err = surface_diffration(Class,Grating_Type,Params.data(),ord,Grating_constant,beam.N_good_rays,
                             beam.X,beam.Y,beam.Z,beam.cX,beam.cY,beam.cZ,beam.lambda,
                             n1,n2,&N_bad,beam.flag);

    if ( err != ENGINE_ERROR_OK ) {
        throw bad_surface(Surface::ActionFailure);
    }

    if ( N_bad ) {
        if ( N_bad == beam.N_good_rays ) throw bad_surface(Surface::NoGoodRays);
        beam.Rearrange();
    }

}


void Grating::SetOrder(const long &ord)
{
    Order[0] = ord;
}


// all the angles must be given in degrees!
void Grating::SetRuleParams(const real_t &blaze_ang, const real_t &alpha_ang, const real_t &gamma_ang,
                            const real_t &rule_start_pos, const real_t &rule_stop_pos)
{
    BlazeAngle = blaze_ang*deg2rad;
    Alpha = alpha_ang*deg2rad;
    Gamma = gamma_ang*deg2rad;

    Rule_start = rule_start_pos;
    Rule_stop = rule_stop_pos;

    if ( Rule_start < 0.0 ) Rule_start = 0.0;
    if ( Rule_start > 1.0 ) Rule_start = 1.0;
    if ( Rule_stop < 0.0 ) Rule_stop = 0.0;
    if ( Rule_stop > 1.0 ) Rule_stop = 1.0;

    if ( Rule_start >= Rule_stop ) throw Grating::bad_grating(Grating::BadRuleParams);
}


void Grating::ApplyQE(vector<real_t> &lambda, vector<real_t> &spec)
{
    if ( lambda.size() != spec.size() ) throw Surface::bad_surface(Surface::BadQE);

    RT_engine_error err = surface_QE(spec.size(),lambda.data(),spec.data(),QE);
//cout << "surface_QE err = " << err << "\n";
    if ( err != ENGINE_ERROR_OK ) {
        throw Surface::bad_surface(Surface::BadQE);
    }

    // compute grating blaze function

//    cout << "Grating angles: " << BlazeAngle << ", " << Alpha << ", " << Gamma << endl;
//    cout << "Const: " << Grating_constant << ", Order: " << Order[0] << endl;

//    err = grating_energy_distr(lambda.size(),lambda.data(),abs(Order[0]),BlazeAngle,Alpha,Gamma,Grating_constant,spec.data());
    err = grating_energy_distr(lambda.size(),lambda.data(),abs(Order[0]),BlazeAngle,Alpha,Gamma,Grating_constant,Rule_start,Rule_stop,spec.data());
//cout << "grating_QE err = " << err << "\n";

    if ( err != ENGINE_ERROR_OK ) {
        throw bad_surface(Surface::BadQE);
    }
}


            /*    STOP CLASS    */

Stop::Stop(const SurfaceClass sclass, const SurfaceShape shape): Surface(sclass,Surface::Stop,shape)
{
//    cerr << "Stop is created!\n";
}

Shade::Shade(const SurfaceClass sclass, const SurfaceShape shape): Surface(sclass,Surface::Shade,shape)
{
//    cerr << "Shade is created!\n";
}


void Shade::ApplyConstrains(Beam &beam)
{
    // apply size constrains for shade type surface
    real_t R2;
//    , ray_R2, xr, yr;

    R2 = Size[0]*Size[0]; // for circular shape

    size_t N_bad = 0;

#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:N_bad)
#endif
#ifdef USING_LINUX
    for (size_t i = 0; i < beam.N_good_rays; ++i ) {
#endif
#ifdef USING_MSVC
    for (long long i = 0; i < beam.N_good_rays; ++i ) {
#endif
        real_t xr = beam.X[i] - Center[0];
        real_t yr = beam.Y[i] - Center[1];
        switch ( Shape ) {
            case Surface::Circle: {
                real_t ray_R2 = xr*xr + yr*yr;
                if ( ray_R2 <= R2 ) {
                    beam.flag[i] = 0;
                    ++N_bad;
                }
                break;
            }
            case Surface::Rectangle: {
                if ( (xr >= -Size[0]/2.0) && (xr <= Size[0]/2.0) &&
                     (yr >= -Size[1]/2.0) && (yr <= Size[1]/2.0) ) {
                    beam.flag[i] = 0;
                    ++N_bad;
                }
                break;
            }
        }
    }

    // rearrange coordinates and cosins vectors if there are new bad rays
    if ( N_bad ) {
        if ( N_bad == beam.N_good_rays ) throw bad_surface(Surface::NoIntersections);
        beam.Rearrange();
    }

}



// just auxiliary surface (used for transformation of coordinate system)
Aux::Aux(): Surface(Surface::Plane,Surface::AuxType,Surface::Circle)
{
//    cerr << "Auxiliary surface is created!\n";
}
