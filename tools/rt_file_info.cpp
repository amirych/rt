#include <iostream>
#include <algorithm>
#include "../common-lib/beam.h"
#include "../common-lib/num_defs.h"

using namespace std;

int main(int argc, char* argv[])
{
    int status = 0;
    real_t *X = NULL;
    real_t *Y = NULL;
    real_t *Z = NULL;
    real_t *cX = NULL;
    real_t *cY = NULL;
    real_t *cZ = NULL;
    real_t *lambda = NULL;

    size_t total_Nrays, good_Nrays, Nlambda;

    if ( argc < 2 ) {
        cerr << "Usage: rt_file_info rt_file\n";
    }

    status = load_beam_data(argv[1],&total_Nrays,&good_Nrays,&Nlambda,
                            &X,&Y,&Z,&cX,&cY,&cZ,&lambda);

    if ( status ) {
        cerr << "Can not read ray-tracing output file!\n";
    } else {
        cout << "Ray-tracing file info:\n";
        cout << "  Filename: " << argv[1] << endl;
        cout << "  Total number of rays: " << total_Nrays << endl;
        cout << "  Number of 'good' rays: " << good_Nrays << endl;
        cout << "  Number of wavelength per range: " << Nlambda << endl;
        cout << endl;

        cout << "  X-axis dimension: [" << *min_element(X,X+good_Nrays-1) << ", " << *max_element(X,X+good_Nrays-1) << "]" << endl;
        cout << "  Y-axis dimension: [" << *min_element(Y,Y+good_Nrays-1) << ", " << *max_element(Y,Y+good_Nrays-1) << "]" << endl;
        cout << "  Z-axis dimension: [" << *min_element(Z,Z+good_Nrays-1) << ", " << *max_element(Z,Z+good_Nrays-1) << "]" << endl;
        cout << endl;

        cout << "  X-axis cosins range: [" << *min_element(cX,cX+good_Nrays-1) << ", " << *max_element(cX,cX+good_Nrays-1) << "]" << endl;
        cout << "  Y-axis cosins range: [" << *min_element(cY,cY+good_Nrays-1) << ", " << *max_element(cY,cY+good_Nrays-1) << "]" << endl;
        cout << "  Z-axis cosins range: [" << *min_element(cZ,cZ+good_Nrays-1) << ", " << *max_element(cZ,cZ+good_Nrays-1) << "]" << endl;
        cout << endl;

        cout << "  Wavelength range: [" << *min_element(lambda,lambda+good_Nrays-1) << ", " << *max_element(lambda,lambda+good_Nrays-1) << "]" << endl;
        cout << endl;
        cout << endl;
    }

    if ( X != NULL ) free(X);
    if ( Y != NULL ) free(Y);
    if ( Z != NULL ) free(Z);
    if ( cX != NULL ) free(cX);
    if ( cY != NULL ) free(cY);
    if ( cZ != NULL ) free(cZ);
    if ( lambda != NULL ) free(lambda);

    return status;
}
