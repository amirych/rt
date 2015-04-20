#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <algorithm>

#include "../common-lib/num_defs.h"
#include "../common-lib/beam.h"

using namespace std;

int main(int argc, char* argv[])
{

    if (argc < 2 ) {
        cout << "Usage: losses_comp rt_files_list result_filename\n";
        return 500;
    }

    char* rt_files_list = argv[1];

    char* result_filename = argv[2];

    int ret_code;
    string filename;
    size_t total_Nrays, Ngood, Nlambda;
    real_t *lambda = NULL;
    real_t *ratio = NULL;
    real_t *ll = NULL;
    real_t *good_lambda = NULL;
    real_t *ptr;
    size_t n,nl,curr_pos;


    ifstream rt_files;
    FILE* result_file;



    try {
        rt_files.open(rt_files_list);

        result_file = fopen(result_filename,"wb");
        if ( result_file == NULL ) {
            throw 100;
        }

        while ( !rt_files.eof() ) {
            rt_files >> filename;

            cout << "Processing " << filename << " ..." << endl;

            ret_code = load_beam_data(filename.c_str(),&total_Nrays,&Ngood,&Nlambda,NULL,NULL,NULL,NULL,NULL,NULL,&lambda,true);
            if ( ret_code ) throw ret_code;

            if ( !Ngood ) {
                cerr << "No good rays in the beam! The file is " << filename << endl;
                continue;
            }

            good_lambda = (real_t*) malloc(Ngood*sizeof(real_t));
            if ( good_lambda == NULL ) throw 300;

            memcpy(good_lambda,lambda,Ngood*sizeof(real_t));

            ratio = (real_t*) calloc(Nlambda,sizeof(real_t));
            if ( ratio == NULL ) throw 300;

            ll = (real_t*) malloc(Nlambda*sizeof(real_t));
            if ( ll == NULL ) throw 300;

            nl = total_Nrays/Nlambda; // maximal number of rays with the same lambda

            // main body ...
//            sort(lambda,lambda+Ngood);
            sort(lambda,lambda+total_Nrays);
            sort(good_lambda,good_lambda+Ngood);
            for ( size_t i_lambda =  0; i_lambda < Nlambda; ++i_lambda ) ll[i_lambda] = lambda[i_lambda*nl];

            curr_pos = 0;
            for ( size_t i_lambda =  0; i_lambda < Nlambda; ++i_lambda ) {
                ptr = good_lambda+curr_pos;
                n = count(ptr,ptr+nl,ll[i_lambda]);
                ratio[i_lambda] = 1.0*n/nl;
                curr_pos += n;
                if ( curr_pos > Ngood ) break;
            }

            // save results
            n = fwrite(&Nlambda,sizeof(size_t),1,result_file);
            if ( n != 1 ) {
                throw 200;
            }

            n = fwrite(ll,sizeof(real_t),Nlambda,result_file);
            if ( n != Nlambda ) {
                throw 200;
            }

            n = fwrite(ratio,sizeof(real_t),Nlambda,result_file);
            if ( n != Nlambda ) {
                throw 200;
            }

            free(good_lambda);
            good_lambda = NULL;
            free(ratio);
            ratio = NULL;
            free(ll);
            ll = NULL;
            free(lambda);
            lambda = NULL;
        }

        fclose(result_file);
        rt_files.close();
    } catch (ios_base::failure &ex) {
        cerr << ex.what();
        return 1;
    } catch (int err) {
        cerr << "An error (code = " << err << ") occured while reading ray-tracing file " << filename << endl;
        if ( good_lambda != NULL ) free(good_lambda);
        if ( ratio != NULL ) free(ratio);
        if ( ll != NULL ) free(ll);
        if ( lambda != NULL ) free(lambda);


        rt_files.close();
        fclose(result_file);

        return 2;
    }

    return 0;
}
