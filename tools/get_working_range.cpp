#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <iomanip>

#include "../common-lib/beam.h"
#include "../common-lib/num_defs.h"

using namespace std;

int main(int argc, char* argv[])
{
    if ( argc < 4 ) {
        cout << "Usage: create_range_file rt_filenames min_order cd_order [coeff]\n";
        return 100;
    }

    ifstream rt_files;
    string filename;

    size_t total_Nrays, good_Nrays, Nlambda;
    real_t *lambda = NULL;
    int st, min_order, cr_order;
    real_t coeff = 1.0;
    real_t min_lambda, max_lambda, range, cen_wave;

    if ( argc > 4 ) coeff = atof(argv[4]);

    min_order = atoi(argv[2]);
    cr_order = atoi(argv[3]);

    rt_files.open(argv[1]);

    st = 0;

    cout << std::fixed;

    while ( !rt_files.eof() ) {
        rt_files >> filename;
        if ( rt_files.fail() || rt_files.bad() ) {
            st = 200;
            break;
        }
        st = load_beam_data(filename.c_str(),&total_Nrays,&good_Nrays,&Nlambda,
                            NULL,NULL,NULL,NULL,NULL,NULL,&lambda);
        if ( st ) {
            cerr << "An error occured while reading files!\n";
            free(lambda);
            break;
        }

        min_lambda = *min_element(lambda,lambda+good_Nrays);
        max_lambda = *max_element(lambda,lambda+good_Nrays);

        range = (max_lambda-min_lambda)*coeff;
        cen_wave = (max_lambda+min_lambda)/2.0;

        min_lambda = cen_wave - range/2.0;
        max_lambda = cen_wave + range/2.0;

        cout << setprecision(7) << min_lambda << "   "
             << setprecision(7) << max_lambda << "   "
             << min_order++ << "   " << cr_order << endl;

        free(lambda);
        lambda = NULL;
    }

    rt_files.close();

    return st;
}
