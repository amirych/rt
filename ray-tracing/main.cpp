#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include "../common-lib/scheme.h"
#include "../common-lib/beam.h"
#include "../common-lib/surface.h"
#include "../common-lib/tabulated_function.h"
#include "../common-lib/rt_engine_errors.h"

#include <sstream>

//#include <cuda_runtime_api.h>
//#include <cuda.h>

using namespace std;

int main(int argc, char* argv[])
{

    Scheme scheme;
    ofstream log_file, rt_files;
    int exit_code = 0;
    stringstream ss;
    list<string> resultfile_list;



    if ( argc < 3 ) {
        cerr << "Usage: ray-tracing scheme_filename rt_files_filename [log_file]" << endl;
        return 1;
    }

    string file_name = argv[1];

//    log_file = new ofstream;
    if ( argc > 3 ) { // log file is given
        log_file.open(argv[3]);
    }

    try {
        cout << "Loading scheme file ... ";

        scheme.load_from_file(file_name);

        cout << "OK!\n";

        cout << "Loaded " << scheme.GetNumberOfSurfaces() << " surfaces.\n";

        cout << "Run ... ";

        if ( argc > 3 ) {
            scheme.run(resultfile_list,log_file);
        } else {
            scheme.run(resultfile_list);
        }

        cout << "OK!\n";

        rt_files.open(argv[2]);

        list<string>::iterator rt_filename = resultfile_list.begin();

        for ( size_t i = 0; i < resultfile_list.size(); ++i, ++rt_filename ) {
            rt_files << (*rt_filename) << endl;
        }

        rt_files.close();

    } catch (bad_alloc &ex) {
        ss << "Memory allocation error was detected!\n";
        exit_code = 1;
    } catch (Scheme::bad_scheme &ex) {
        Scheme::SchemeError err = ex.scheme_error();
        switch (err) {
            case Scheme::FileOpenFailure:
                ss << "Can not open scheme file!\n";
                exit_code = 1;
                break;
            case Scheme::ParseFailure: {
                ss << "A parsing error occured while loading scheme file!\n";
                exit_code = 1;
                break;
            }
            case Scheme::BadAttrValue: {
                ss << "Bad attribute in the input scheme file!\n";
                exit_code = 1;
                break;
            }
            case Scheme::MemAllocFailure: {
                ss << "Memory allocation error occured while parsing input scheme!\n";
                exit_code = 1;
                break;
            }
            case Scheme::NoSchemeDesc: {
                ss << "There was no general scheme description in the input scheme file!\n";
                exit_code = 1;
                break;
            }
            case Scheme::NoBeamDesc: {
                ss << "There was no beam description in the input scheme file!\n";
                exit_code = 1;
                break;
            }
            case Scheme::NoSurfaceDesc: {
                ss << "No one surface description was found in the input scheme file!\n";
                exit_code = 1;
                break;
            }
            case Scheme::NoRange: {
                ss << "There was neither range attrribute or range file in the input scheme file!\n";
                exit_code = 1;
                break;
            }
            case Scheme::NoRangeFile: {
                ss << "Can not find range file!\n";
                exit_code = 1;
                break;
            }
            case Scheme::EmptyScheme: {
                ss << "Try to run ray-tracing in an empty sheme! Load scheme first\n";
                exit_code = 1;
                break;
            }
            case Scheme::BadRange: {
                ss << "Ther was bad range description in the input range file!\n";
                exit_code = 1;
                break;
            }
            case Scheme::BadStartIndex: {
                ss << "Start surface index is greater than number of surfaces in the scheme!\n";
                exit_code = 1;
                break;
            }

            default: {
                ss << "Unknown scheme error!\n";
                exit_code = 1;
                break;
            }
        }
    } catch (Beam::bad_beam &ex) {
        Beam::BeamError err = ex.beam_error();
        switch (err) {
            case Beam::MemAllocFailure: {
                ss << "Memory allocation error occured during beam creation!\n";
                exit_code = 1;
                break;
            }
            case Beam::InvalidParams: {
                ss << "Beam invalid parameter!\n";
                exit_code = 1;
                break;
            }
            case Beam::SaveFailure: {
                ss << "An error occured while beam saving!\n";
                exit_code = 1;
                break;
            }
            case Beam::CreationFailure: {
                ss << "Can not create a beam!\n";
                exit_code = 1;
                break;
            }
            case Beam::TransformationFailure: {
                ss << "An error occured while beam transformation!\n";
                exit_code = 1;
                break;
            }
            default: {
                ss << "Unknown beam error!";
                exit_code = 1;
                break;
            }
        }
    } catch (Surface::bad_surface &ex) {
        Surface::SurfaceError err = ex.surface_error();
        switch (err) {
        case Surface::InvalidSurfaceAngles: {
            ss << "Invalid surface angles!\n";
            exit_code = 1;
            break;
        }
        case Surface::InvalidSurfaceDistance: {
            ss << "Invalid surface distance!\n";
            exit_code = 1;
            break;
        }
        case Surface::InvalidSurfaceParam: {
            ss << "Invalid surface parameter!\n";
            exit_code = 1;
            break;
        }
        case Surface::InvalidSurfaceSize: {
            ss << "Invalid surface size!\n";
            exit_code = 1;
            break;
        }
        case Surface::InvalidSurfaceCenter: {
            ss << "Invalid surface center!\n";
            exit_code = 1;
            break;
        }
        case Surface::NoGoodRays: {
            ss << "";
            exit_code = 1;
            break;
        }
        case Surface::NoIntersections: {
            ss << "The surface has no intersections with current beam!\n";
            exit_code = 1;
            break;
        }
        case Surface::IntersectionFailure: {
            ss << "Ray-tracing engine intesection function failed!\n";
            exit_code = 1;
            break;
        }
        case Surface::ActionFailure: {
            ss << "Ray-tracing engine action function failed!\n";
            exit_code = 1;
            break;
        }
        default: {
            ss << "Unknown surface error!\n";
            exit_code = 1;
            break;
        }
        }
    } catch (TabulatedFunction::bad_tab_func &ex) {
        TabulatedFunction::TabulatedFunctionError err = ex.bad_tab_func_error();
        switch (err) {
        case TabulatedFunction::MemAllocError: {
            ss << "Memory allocation error while compute tabulated function!\n";
            exit_code = 1;
            break;
        }
        case TabulatedFunction::BadValue: {
            ss << "Bad value of tabulated function!\n";
            exit_code = 1;
            break;
        }
        case TabulatedFunction::FileFailure: {
            ss << "Can not open tabulated function file!\n";
            exit_code = 1;
            break;
        }
        case TabulatedFunction::BadArgument: {
            ss << "Argument of tabulated function must be in increasing order!\n";
            exit_code = 1;
            break;
        }
        case TabulatedFunction::OutOfRange: {
            ss << "Interpolated point is out of range of tabulated function!\n";
            exit_code = 1;
            break;
        }
        case TabulatedFunction::InitFailure: {
            ss << "Computation of second derivatives of tabulated function was failed!\n";
            exit_code = 1;
            break;
        }
        case TabulatedFunction::NotEnoughPoints: {
            ss << "Not enough points for cubic interpolation!\n";
            exit_code = 1;
            break;
        }
        default: {
            ss << "Unknown tabulated function error!\n";
            exit_code = 1;
            break;
        }
        }
    } catch (...) {
        ss << "Something is wrong!!!\n";
        exit_code = 1;
    }

    cerr << ss.str();

    if ( argc > 2 ) {
        log_file.close();
    }

    return exit_code;
}

