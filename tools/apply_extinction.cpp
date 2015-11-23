#include<iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "../common-lib/num_defs.h"
#include "../common-lib/tabulated_function.h"

using namespace std;

int main(int argc, char* argv[])
{
    if ( argc < 4 ) {
        cout << "Usage: QE_file ext_file result_QE_file\n";
        return 1;
    }

    FILE *QE_file, *result_QE_file;

    size_t Npoints;
    real_t *QE_lambda = NULL;
    real_t *QE_val = NULL;

    TabulatedFunction ext_curve;

    QE_file = fopen(argv[1],"rb");
    if ( QE_file == NULL ) {
        cerr << "Can not open QE file!\n";
        return 2;
    }

    result_QE_file = fopen(argv[3],"wb");
    if ( result_QE_file == NULL ) {
        cerr << "Can not open result QE file!\n";
        fclose(QE_file);
        return 2;
    }

    int err_code = 0;
    size_t counts;
    size_t i_range = 1;
    try {
        ext_curve.Load(argv[2]);

        while ( !feof(QE_file) ) {
            // read QE for current range

            counts = fread(&Npoints,sizeof(size_t),1,QE_file);
            if ( counts != 1 ) {
                err_code = 10;
                throw err_code;
            }

            cout << "Reading " << i_range << "-th range " << "(" << Npoints << " points)  ...";

            QE_lambda = (real_t*) malloc(sizeof(real_t)*Npoints);
            if ( QE_lambda == NULL ) {
                err_code = 20;
                throw err_code;
            }

            QE_val = (real_t*) malloc(sizeof(real_t)*Npoints);
            if ( QE_val == NULL ) {
                err_code = 20;
                throw err_code;
            }

            counts = fread(QE_lambda,sizeof(real_t),Npoints,QE_file);
            if ( counts != Npoints ) {
                err_code = 10;
                throw err_code;
            }

            counts = fread(QE_val,sizeof(real_t),Npoints,QE_file);
            if ( counts != Npoints ) {
                err_code = 10;
                throw err_code;
            }

            cout << "   Done!\n";

            // compute corrected for extinction QE (apply extinction to QE curve)

            cout << "   Computing and writing QE curve for " << i_range << "-th range ...";
            for ( size_t i = 0; i < Npoints; ++i ) {
                QE_val[i] *= pow(10.0,-0.4*ext_curve[QE_lambda[i]]);
            }

            // write corrected QE

            counts = fwrite(&Npoints,sizeof(size_t),1,result_QE_file);
            if ( counts != 1 ) {
                err_code = 15;
                throw err_code;
            }

            counts = fwrite(QE_lambda,sizeof(real_t),Npoints,result_QE_file);
            if ( counts != Npoints ) {
                err_code = 15;
                throw err_code;
            }

            counts = fwrite(QE_val,sizeof(real_t),Npoints,result_QE_file);
            if ( counts != Npoints ) {
                err_code = 15;
                throw err_code;
            }

            cout << "   Done!\n";

            free(QE_lambda);
            free(QE_val);

            QE_lambda = NULL;
            QE_val = NULL;

            ++i_range;
        }
    } catch (TabulatedFunction::bad_tab_func &ex) {
        switch ( ex.bad_tab_func_error() ) {
        case TabulatedFunction::FileFailure: {
            cout << "Can not open extinction file!\n";
            err_code = 2;
            break;
        }
        case TabulatedFunction::OutOfRange: {
            cout << "Input QE lambda is out of range extiction curve!\n";
            err_code = 30;
            break;
        }
        default: {
            cout << "Something wrong with extinction curve!\n";
            err_code = 40;
            break;
        }
        }
    } catch (int ex) {
        switch (ex) {
        case 10: {
            cout << "Something was wrong while reading QE file!\n";
            break;
        }
        case 15: {
            cout << "Something was wrong while writing result QE file!\n";
            break;
        }
        case 20: {
            cout << "Memory allocation error!\n";
            break;
        }
        default:
            break;
        }
    }

    fclose(QE_file);
    fclose(result_QE_file);

    free(QE_lambda);
    free(QE_val);

    return err_code;
}
