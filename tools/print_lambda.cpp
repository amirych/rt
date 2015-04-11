#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstdlib>

#include "../common-lib/num_defs.h"
#include "../common-lib/beam.h"

#define STR_MAXLEN 1024

using namespace std;

int main(int argc, char* argv[])
{
    ifstream rt_files;

    size_t Nrays, N_good, N_lambda, N_uniq, j;

    char rt_filename[STR_MAXLEN+1];
    int status, flag;

    real_t *lambda,*uniq;

    if ( argc < 2 ) {
        cerr << "Usage: print_lambda rt_list_filename\n";
        return 1;
    }

    cout << std::fixed;
    cout << std::setprecision(9);

    flag = 0;
    rt_files.open(argv[1]);

    while ( !rt_files.eof() ) {
        rt_files.getline(rt_filename,STR_MAXLEN);
        if ( rt_files.eof() ) break;

        if ( rt_files.fail() || rt_files.bad() ) { // thomething is wrong!
            cerr << "\nAn error occured while read the RT-list file!!!\n";
            flag = 1;
            break;
        }

        try {
            status = load_beam_data(rt_filename,&Nrays,&N_good,&N_lambda,NULL,NULL,NULL,NULL,NULL,NULL,&lambda);
            if ( status ) throw 10;

            uniq = (real_t*) calloc(N_lambda,sizeof(real_t));
            uniq[0] = lambda[0];
            N_uniq = 1;

            cout << "# " << rt_filename << endl;


            cout << lambda[0] << endl;
            for ( size_t i = 1; i < N_good; ++i ) {
                for ( j = 0; j < N_uniq; ++j ) {
                    if ( lambda[i] == uniq[j] ) {
                        break;
                    }
                }
                if ( j == N_uniq ) {
                    cout << lambda[i] << endl;
                    uniq[N_uniq] = lambda[i];
                    ++N_uniq;
                }
                if ( N_uniq == N_lambda ) break;
            }
            cout << "\n";

            free(uniq);
            free(lambda);
        } catch (int err) {
            flag = err;
            switch ( err ) {
                case 1: {
                    cerr << "An error occured while reading RT-file!\n";
                }
            }
            break;
        }
    }

    rt_files.close();

    return flag;
}
