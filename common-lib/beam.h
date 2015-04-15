#ifndef BEAM_H
#define BEAM_H

#include "num_defs.h"

#include<vector>
#include<string>
#include<exception>
#include <stdarg.h>

using namespace std;

/*
     declaration of function for reading beam data.
        if argument is NULL-pointer then skip reading
*/
int load_beam_data(const char *rt_filename, size_t *total_Nrays, size_t *N_good, size_t *N_lambda,
                   real_t **X, real_t **Y, real_t **Z, real_t **cX, real_t **cY, real_t **cZ, real_t **lambda, bool all = false);

class Surface;

class Beam
{
    friend class Surface;
    friend class Mirror;
    friend class Lens;
    friend class Grating;
    friend class Shade;
    friend class Aux;
public:
    enum BeamType {Parallel,Conic,Gauss};
    enum BeamShape {Circle,Rectangle};
    enum BeamProfile {Random,Concentric};
    enum BeamRangeDistr {RandomDistr,UniformDistr,Monochromatic};

    enum BeamError {MemAllocFailure,InvalidParams,SaveFailure,CreationFailure,TransformationFailure};

    Beam();
    Beam(Beam::BeamType type, Beam::BeamShape shape, Beam::BeamProfile profile, BeamRangeDistr range_distr);

    Beam& operator= (const Beam &beam);

    ~Beam();

    void SetParams(const vector<real_t> &params, bool composite_flag = false, size_t Nlambda = 1);

    void SetCenter(const vector<real_t> &center);

    void SetRange(const vector<real_t> &range);

    void Create(const vector<real_t> &params, const vector<real_t> &center, const vector<real_t> &range,
                bool composite_flag = false, size_t Nlambda = 1);

    void Recreate(); // create beam with already defined parameters

    void Save(const string &filename);

    void Rearrange();

    void Transform(const vector<real_t> &dist, const vector<real_t> &angs);
    void Transform(const Surface &surf);

    vector<real_t> GetBeamDim();

    vector<real_t> GetBeamCosins();

    void Normalize();

    size_t GetTotalNumberOfRays() const;
    size_t GetNumberOfGoodRays() const;

    // Beam class exception class
    class bad_beam: public exception
    {
    public:
        bad_beam(BeamError err);
        BeamError beam_error() const;
    private:
        BeamError Error;
    };

private:
    BeamType Type;
    BeamShape Shape;
    BeamProfile Profile;
    BeamRangeDistr Range_Distr;
    bool CompositeBeam;

    vector<real_t> Params;
    vector<real_t> Center;
    vector<real_t> Range;

    size_t Nrays;       // total number of rays
    size_t N_good_rays; // number of good rays
    size_t N_lambda;    // number of wavelength. If CompositeBeam is false than Nlambda = Nrays, otherwise Nrays = Nlambda*(number of rays per each wavelength)
    real_t *X,*Y,*Z;    // coordinates of rays in beam
    real_t *cX,*cY,*cZ; // direction cosins of rays in beam
    real_t *lambda;     // wavelength of rays
    beam_flag_t *flag;  // goodness flag for each ray

    void destroy();
    void swap_elements(size_t i1, size_t i2);
};

#endif // BEAM_H
