#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>

#include <fitsio.h>

#include "../common-lib/num_defs.h"
#include "../common-lib/ascii_file.h"
#include "../common-lib/tabulated_function.h"
#include "../common-lib/rt_engine_func.h"
#include "../common-lib/rt_engine_errors.h"
#include "../common-lib/beam.h"


#define STR_MAXLEN 1024

using namespace std;

int main(int argc, char* argv[])
{
    AsciiFile file;
    real_t ccd_xsize, ccd_ysize; // linear size of CCD
    real_t ccd_xdim, ccd_ydim;   // number of pixels along X and Y axes
    real_t ccd_xpix, ccd_ypix;   // linear size of CCD pixel
    int status;

    AsciiFile::AsciiFileFlag file_flag;

    TabulatedFunction QE, spec;

    if ( argc < 6 ) {
        cerr << "Usage: ccd_frame rt_files_list QE_filename spec_filename ccd_spec_filename result_FITS_filename\n";
        return 1;
    }

    // read CCD specification file

    file.open(argv[4]);
    if ( file.fail() ) {
        cerr << "Can not open CCD specification file!\n";
        return 2;
    }

    while ( (file_flag = file.ReadLine(2,&ccd_xsize,&ccd_ysize)) != AsciiFile::DataString) {
        if ( (file_flag == AsciiFile::InvalidData) || (file_flag == AsciiFile::Eof) ) { // something wrong
            file.close();
            cerr << "Invalid CCD linear size specification!\n";
            return 2;
        }
    }

    while ( (file_flag = file.ReadLine(2,&ccd_xdim,&ccd_ydim)) != AsciiFile::DataString) {
        if ( (file_flag == AsciiFile::InvalidData) || (file_flag == AsciiFile::Eof) ) { // something wrong
            file.close();
            cerr << "Invalid CCD linear size specification!\n";
            return 2;
        }
    }

    while ( (file_flag = file.ReadLine(2,&ccd_xpix,&ccd_ypix)) != AsciiFile::DataString) {
        if ( (file_flag == AsciiFile::InvalidData) || (file_flag == AsciiFile::Eof) ) { // something wrong
            file.close();
            cerr << "Invalid CCD linear size specification!\n";
            return 2;
        }
    }

    cout << "CCD specifications:\n";
    cout << "Linear size: [" << ccd_xsize << ":" << ccd_ysize << "]" << endl;
    cout << "Pixel linear size: [" << ccd_ypix << ":" << ccd_ypix << "]" << endl;
    cout << "Number of pixels: [" << ccd_xdim << ":" << ccd_ydim << "]" << endl;

    file.close();

    long naxis = 2;
    long naxes[2];

    naxes[0] = static_cast<long>(ccd_xdim);
    naxes[1] = static_cast<long>(ccd_ydim);

    size_t ccd_xdim_i = static_cast<size_t>(ccd_xdim);
    size_t ccd_ydim_i = static_cast<size_t>(ccd_ydim);

    size_t image_size = naxes[0]*naxes[1];

    FILE *QE_file;


    // open QE and read input spectrum files

//    try {
//        QE.Load(argv[2]);
//    } catch (TabulatedFunction::bad_tab_func &ex) {
//        TabulatedFunction::TabulatedFunctionError err = ex.bad_tab_func_error();
//        cerr << "Can not load QE function! Error code: " << err << endl;
//        return 3;
//    }

    QE_file = fopen(argv[2],"rb");
    if ( QE_file == NULL ) {
        cerr << "Can not open QE file!\n";
        return 3;
    }

    try {
        spec.Load(argv[3]);
    } catch (TabulatedFunction::bad_tab_func &ex) {
        TabulatedFunction::TabulatedFunctionError err = ex.bad_tab_func_error();
        cerr << "Can not load input spectrum function! Error code: " << err << endl;
        return 3;
    }

    // compute image

    RT_engine_error err;
    real_t *image;
    char rt_filename[STR_MAXLEN];

    real_t *X = NULL;
    real_t *Y = NULL;
    real_t *lambda = NULL;

    unsigned int *ccd_X = NULL;
    unsigned int *ccd_Y = NULL;

    size_t Nrays,N_good,N_lambda,Npoints,n;

    real_t *QE_lambda = NULL;
    real_t *QE_val = NULL;

    int flag = 0;

    try {
        image = new real_t[image_size](); // allocate with zeroing
    } catch (bad_alloc &ex) {
        cerr << "Can not allocate memory for the result image!\n";
        return 100;
    }

    file.open(argv[1]);
    if ( file.fail() ) {
        cerr << "Can not open file with list of ray-tracing result files!\n";
        return 5;
    }

    cout << "\nSTART CCD IMAGE CREATION ...\n";

    while ( !file.eof() ) {
        file.getline(rt_filename,STR_MAXLEN);
        if ( file.eof() ) break;

        if ( file.fail() || file.bad() ) { // thomething is wrong!
            cerr << "\nAn error occured while read the file!!!\n";
            flag = 1;
            break;
        }

        cout << "   Processing file: " << rt_filename << "...  \n";

        try {
            // QE data for given range
            n = fread(&Npoints,sizeof(size_t),1,QE_file);
            if ( n != 1 ) throw 30;

            QE_lambda = (real_t*) malloc(Npoints*sizeof(real_t));
            if ( QE_lambda == NULL ) throw 50;

            QE_val = (real_t*) malloc(Npoints*sizeof(real_t));
            if ( QE_val == NULL ) throw 50;

            n = fread(QE_lambda,sizeof(real_t),Npoints,QE_file);
            if ( n != Npoints ) throw 30;

            n = fread(QE_val,sizeof(real_t),Npoints,QE_file);
            if ( n != Npoints ) throw 30;

            QE.Set(Npoints,QE_lambda,QE_val);

            status = load_beam_data(rt_filename,&Nrays,&N_good,&N_lambda,&X,&Y,NULL,NULL,NULL,NULL,&lambda);
            if ( status ) throw 1;

            ccd_X = (unsigned int*) malloc(N_good*sizeof(unsigned int));
            if ( ccd_X == NULL ) throw 50;

            ccd_Y = (unsigned int*) malloc(N_good*sizeof(unsigned int));
            if ( ccd_Y == NULL ) throw 50;

            err = ccd_coordinates(N_good,X,Y,ccd_xsize,ccd_ysize,ccd_xpix,ccd_ypix,ccd_X,ccd_Y);
            if ( err != ENGINE_ERROR_OK ) throw 10;
//            cout << "  X = " << X[10] << "; Y = " << Y[10] << endl;
//            cout << "  ccdX = " << ccd_X[10] << "; ccdY = " << ccd_Y[10] << endl;

            err = compute_ccd_image(N_good,ccd_X,ccd_Y,lambda,QE,spec,ccd_xdim_i,ccd_ydim_i,image);
            if ( err != ENGINE_ERROR_OK ) throw 20;

        } catch (int err_code) {
            flag = 1;
            switch ( err_code ) {
                case 1: {
                    cerr << "\nAn error occured while read the ray-tracing output file!\n";
                    break;
                }
                case 10: {
                    cerr << "\nAn error occured while compute CCD coordinates!\n";
                    break;
                }
                case 20: {
                    cerr << "\nAn error occured while compute CCD image!\n";
                    break;
                }
                case 30: {
                    cerr << "\nAn error occured while read QE file!\n";
                    break;
                }
                case 50: {
                    cerr << "\nMemory allocation failure!\n";
                    break;
                }
            }

        } catch (TabulatedFunction::bad_tab_func &ex) {
            flag = 1;
            cerr << "\nError occured in QE data creation! Error code: " << ex.bad_tab_func_error() << endl;
        }


        free(X);
        free(Y);
        free(lambda);
        free(ccd_X);
        free(ccd_Y);
        X = NULL;
        Y = NULL;
        ccd_X = NULL;
        ccd_Y = NULL;

        free(QE_lambda);
        free(QE_val);
        QE_lambda = NULL;
        QE_val = NULL;

        if ( flag ) break; // an error occured!
//        cout << "OK!\n";
    }

    file.close();

    fclose(QE_file);

    if ( flag ) { // error occured while read list
        delete[] image;
        return 5;
    }

    cout << "DONE!\n";

    // create and output FITS image

    cout << "\nSAVE IMAGE TO FITS FILE ...  ";

    int fits_status = 0;

    fitsfile *fptr;
    string fits_filename = argv[5];

    fits_filename = "!" + fits_filename; // add "!" to erase already existing file

    fits_create_file(&fptr,fits_filename.c_str(),&fits_status);
    if ( fits_status ) {
        cerr << "Can not create result FITS file! Error code " << fits_status << endl;
        delete[] image;
        return 4;
    }

    int bitpix,datatype;
    long fpix[2];

    fpix[0] = 1;
    fpix[1] = 1;

#ifdef RT_NUM_DOUBLE
    bitpix = DOUBLE_IMG;
    datatype = TDOUBLE;
#else
    bitpix = FLOAT_IMG;
    datatype = TFLOAT;
#endif

    fits_create_img(fptr, bitpix, naxis, naxes, &fits_status);
    if ( fits_status ) {
        cerr << "Can not create result FITS image! Error code " << fits_status << endl;
        fits_close_file(fptr,&fits_status);
        delete[] image;
        return 4;
    }

    fits_write_pix(fptr,datatype,fpix,image_size,image,&fits_status);
    if ( fits_status ) {
        cerr << "Can not create result FITS image! Error code " << fits_status << endl;
        fits_close_file(fptr,&fits_status);
        delete[] image;
        return 4;
    }

    // fill header

    fits_write_key(fptr,TSTRING,"RT_LIST",argv[1],"",&fits_status);
    fits_write_key(fptr,TSTRING,"QE_FILE",argv[2],"",&fits_status);
    fits_write_key(fptr,TSTRING,"SPECFILE",argv[3],"",&fits_status);
    fits_write_key(fptr,TSTRING,"CCDINFO",argv[4],"",&fits_status);

    fits_close_file(fptr,&fits_status);

    delete[] image;

    cout << "OK!\n";

    return 0;
}
