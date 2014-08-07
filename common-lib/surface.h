#ifndef SURFACE_H
#define SURFACE_H

#include "num_defs.h"
#include "beam.h"
#include "tabulated_function.h"

#include<list>
#include<vector>
#include<exception>


#define MAX_POLY_COEFFS 10


using namespace std;


        /* Base class of surface */

class Surface
{
public:
    enum SurfaceClass {Plane,Conic,Toric,AuxClass};
    enum SurfaceType {Mirror,Lens,Grating,Stop,Shade,Detector,AuxType};
    enum SurfaceShape {Circle,Rectangle};

    enum SurfaceError {Ok,InvalidSurfaceParam,InvalidSurfaceSize,InvalidSurfaceCenter,
                       InvalidSurfaceDistance,InvalidSurfaceAngles,
                       MemAllocFailure,NoIntersections,IntersectionFailure,
                       ActionFailure,NoGoodRays,BadQE};

    Surface(const SurfaceClass sclass, const SurfaceType type, const SurfaceShape shape);
    ~Surface();

    SurfaceError GetError() const;

    void SetParams(const vector<real_t> &pars);
    void SetGeometry(const vector<real_t> &size, const vector<real_t> &center = vector<real_t>(2,0.0));
    void SetDistance(const vector<real_t> &dist);
    void SetAngles(const vector<real_t> &angs); // input angles values are in degrees (to be transformed in the function to radians)

    void SetQE(const real_t qe_val);
    void SetQE(const char* qe_file);

    void SetComment(const string &comm);

    vector<real_t> GetDistance() const;
    vector<real_t> GetAngles() const;
    string GetComment() const;

    Surface::SurfaceClass GetClass() const;
    Surface::SurfaceType GetType() const;
    Surface::SurfaceShape GetShape() const;

    void Intersection(Beam &beam);
    virtual void ApplyConstrains(Beam &beam);
    virtual void Action(Beam &beam);
    virtual void ApplyQE(vector<real_t> &lambda, vector<real_t> &spec);

    // Surface class exception
    class bad_surface: public exception
    {
    public:
        bad_surface(SurfaceError err);
        SurfaceError surface_error() const;
    private:
        SurfaceError Error;
    };

protected:
    SurfaceClass Class;
    SurfaceType Type;
    SurfaceShape Shape;
    string Comment;

    vector<real_t> Params;
    vector<real_t> Size;
    vector<real_t> Center;

    vector<real_t> Distance;
    vector<real_t> Angles;

    TabulatedFunction QE; // quantum efficiency of the surface

    SurfaceError CurrentError;
};


        /*  Inherited surface classes  */

class Plane: public Surface
{
public:
    Plane(const SurfaceType type, const SurfaceShape shape);
};


class Mirror: public Surface
{
public:
    Mirror(const SurfaceClass sclass, const SurfaceShape shape);
    virtual void Action(Beam &beam);
};


class Lens: public Surface
{
public:
    enum LensError {BadRefIndex,BadRefIndexFile};

    Lens(const SurfaceClass sclass, const SurfaceShape shape, const real_t &ref_index1, const real_t &ref_index2);
    Lens(const SurfaceClass sclass, const SurfaceShape shape, const char* ref_index1, const char* ref_index2);
    Lens(const SurfaceClass sclass, const SurfaceShape shape, const real_t &ref_index1, const char* ref_index2);
    Lens(const SurfaceClass sclass, const SurfaceShape shape, const char* ref_index1, const real_t &ref_index2);

    virtual void Action(Beam &beam);

    class bad_lens: public exception
    {
    public:
        bad_lens(Lens::LensError err);
        Lens::LensError bad_lens_error() const;
    private:
        Lens::LensError Error;
    };
private:
    TabulatedFunction n1, n2; // refractive index before and after the surface
};


class Grating: public Surface
{
public:
    enum GratingType {Reflective,Transparent};
    enum GratingError {BadRefIndex,BadRefIndexFile,BadOrder,BadType};

    Grating(const SurfaceClass sclass, const SurfaceShape shape, const Grating::GratingType gr_type, const vector<long> &ord,
            const real_t &gr_const, const real_t &ref_index1, const real_t &ref_index2);

    Grating(const SurfaceClass sclass, const SurfaceShape shape, const Grating::GratingType gr_type, const vector<long> &ord,
            const real_t &gr_const, const char* ref_index1, const char* ref_index2);

    Grating(const SurfaceClass sclass, const SurfaceShape shape, const Grating::GratingType gr_type, const vector<long> &ord,
            const real_t &gr_const, const real_t &ref_index1, const char* ref_index2);

    Grating(const SurfaceClass sclass, const SurfaceShape shape, const Grating::GratingType gr_type, const vector<long> &ord,
            const real_t &gr_const, const char* ref_index1, const real_t &ref_index2);

    virtual void Action(Beam &beam);
//    void Action(Beam &beam, const long &ord);

    void SetOrder(const long &ord);

    // all the angles must be given in degrees!
    void SetRuleParams(const real_t &blaze_ang, const real_t &alpha_ang, const real_t &gamma_ang);

    virtual void ApplyQE(vector<real_t> &lambda, vector<real_t> &spec);

    class bad_grating: public exception
    {
    public:
        bad_grating(Grating::GratingError err);
        Grating::GratingError bad_grating_error() const;
    private:
        Grating::GratingError Error;
    };
private:
    GratingType Grating_Type;
    vector<long> Order;       // diffraction order or a range of orders
    real_t Grating_constant;  // assumed it is given in mkm!
    TabulatedFunction n1, n2; // refractive index before and after the surface (for transparent grating)

    real_t BlazeAngle;        // blaze angle of grating rule (in radians)
    real_t Alpha;             // incident angle (in radians)
    real_t Gamma;             // incident angle along rule (in radians)
};



class Stop: public Surface
{
public:
    Stop(const SurfaceClass sclass, const SurfaceShape shape);
};


class Shade: public Surface
{
public:
    Shade(const SurfaceClass sclass, const SurfaceShape shape);
    virtual void ApplyConstrains(Beam &beam);
};


class Aux: public Surface
{
public:
    Aux(); // no parameters are needed!
};


#endif // SURFACE_H
