#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "../common-lib/num_defs.h"

using namespace std;

int main(int argc, char* argv[])
{
    FILE *QE_file;
    size_t n, Npoints;
    real_t *QE_lambda = NULL;
    real_t *QE_val = NULL;
    int flag = 0;


    if ( argc < 2 ) {
        cerr << "Usage: print_QE QE_file\n";
        return 1;
    }


    QE_file = fopen(argv[1],"rb");
    if ( QE_file == NULL ) {
        cerr << "Can not open QE file!\n";
        return 2;
    }

    cout << "QE file:\n";

    while ( !feof(QE_file) ) {
        try {
            n = fread(&Npoints,sizeof(size_t),1,QE_file);
            if ( n != 1 ) throw 10;

            QE_lambda = (real_t*) malloc(sizeof(real_t)*Npoints);
            if ( QE_lambda == NULL ) throw 20;

            QE_val = (real_t*) malloc(sizeof(real_t)*Npoints);
            if ( QE_val == NULL ) throw 20;

            n = fread(QE_lambda,sizeof(real_t),Npoints,QE_file);
            if ( n != Npoints ) throw 10;

            n = fread(QE_val,sizeof(real_t),Npoints,QE_file);
            if ( n != Npoints ) throw 10;

            cout << endl;
            for ( size_t i = 0; i < Npoints; ++i ) {
//                cout << QE_lambda[i] << " " << QE_val[i] << endl;
                printf("%10.8f  %5.3f\n",QE_lambda[i],QE_val[i]);
            }
        } catch ( int err ) {
            flag = err;
            switch (err) {
                case 10: {
                    if ( !feof(QE_file) ) cerr << "Error while reading the file!\n";
                    break;
                }
                case 20: {
                    cerr << "Memmory allocation error!\n";
                    break;
                }
            }
        }

        free(QE_lambda);
        free(QE_val);
        QE_lambda = NULL;
        QE_val = NULL;

        if ( flag ) break;
    }

    fclose(QE_file);
    return flag;
}
