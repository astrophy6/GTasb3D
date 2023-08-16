/*****************************************************************************************************************************
 *
 *  GTasb3D:  C++ program for generalized finite differences thermal modeling
 *
 *  Developed by Yun Zhang (2022-2023)
 *
 *  Function: solve 3D transient heat conduction in a given irregular-shape body
 *  Feature:
 *
 *
 *****************************************************************************************************************************/


#include "master.hpp"
#include <ctime>


int main(int argc, char *argv[]) {

    Master S;
    int iOpt;

    if ( argc < 2 ) {
        cerr << "Usage: " << argv[0] << "[-n N] input.par\n"
             << "Options:\n"
             << "\t-n N: \t specify number of threads for parallel computing"<<endl;
        return 1;
    }
    while ((iOpt = getopt(argc, argv, "n:")) != -1) {
        switch (iOpt) {
            case 'n':
                S.iNthreads = atoi(optarg); // number of threads for parallel computing
                if ( S.iNthreads < 1 )
                    cerr << "Number of thread must be larger thant 1! (Input: " << S.iNthreads << ")" << endl;
                break;
            default:
                cerr << "Usage: " << argv[0] << "[-n N] input.par\n"
                << "Options:\n"
                << "\t-n N: \t specify number of threads for parallel computing"<<endl;
                return 1;
        }
    }
    
#ifdef USE_GAS
    cout << "USE_GAS open." << endl;
#endif
    
#ifdef USE_SELF
    cout << "USE_SELF open." << endl;
#endif
    
    omp_set_num_threads(S.iNthreads);
    cout << "GTasb3D running on " << S.iNthreads << " processors\n" << endl;

#ifdef SEETIME
    double t0,t1,t2;
    t0 = omp_get_wtime();
#endif

    
    /* Initializing the simulation */
    if ( S.InitSystem(argv[argc-1]) != 0 ) {
        cout << "Error occurred while initializing master." << endl;
        return 1;
    }
    S.Output();  /* output nodes' information */

#ifdef SEETIME
    t1 = omp_get_wtime();
    cout << "The input elapsed time: " << (double)(t1-t0) << " s." << endl;
#endif

    cout << "Simulation starts." << endl;

    while( S.CurrentStep() <= S.TotalSteps() ){
        if ( S.Run() )
            break; /* encounter running error */
        S.Output();  /* output nodes' information */
#ifdef SEETIME
        t2 = omp_get_wtime();
        cout<<"The elapsed time for "<<S.OutInterval()<<" runs: "<<(double)(t2-t1)<<" s."<<endl;
        t1 = omp_get_wtime();
#endif
    }
    cout << "Simulation ends." << std::endl;
#ifdef SEETIME
    t2 = omp_get_wtime();
    cout<<"The total elapsed time: "<<(double)(t2-t0)<<" s."<<endl;
#endif
    
    return 0;
}
