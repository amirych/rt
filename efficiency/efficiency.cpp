#include <iostream>
#include <mxml.h>
#include <cstdlib>
#include <vector>

#include "../common-lib/num_defs.h"
#include "../common-lib/scheme.h"

using namespace std;

int main(int argc, char* argv[])
{

    Scheme scheme;
    vector<real_t> QE_curve;
    size_t Nlambda;

    if ( argc < 3 ) {
        cerr << "Usage: qe_comp scheme_file number_of_wavelengths result_file\n";
        return 1;
    }

    try {
        char *endptr;
#ifdef USING_LINUX
        Nlambda = strtoll(argv[2],&endptr,10);
#else // Visual Studio
        Nlambda = _strtoi64(argv[2],&endptr,10);
#endif
        if ( endptr[0] != '\0' ) {
            throw 10;
        }

        cout << "Load scheme file ...";
        scheme.load_from_file(argv[1],true);
        cout << "  OK!\n";

        cout << "Compute QE curve ...";
        scheme.ComputeQE(Nlambda,argv[3]);
        cout << "  OK!\n";
    } catch (Scheme::bad_scheme &ex) {
        cout << "\nSomething is wrong!\n";
        cout << "Error: " << ex.scheme_error() << endl;
    }

    return 0;
}
