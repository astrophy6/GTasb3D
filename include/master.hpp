#ifndef master_hpp
#define master_hpp

#include "common.hpp"
#include "fdmnode.hpp"

class Master {  /* control all master information and process */
    private:
        /* user-define constant */
        bool _bReadLatest;  /* =1: read the latest file */
        unsigned short int  _nDigits;   /* output file digit mask; comparable with nSteps */
        unsigned int _iStartStep;  /* initial step numbering */
        unsigned int  _nSteps;    /* total number of steps in intervals of dDelta */
        unsigned int _nOutInterval; /* general output file interval in timesteps */
        unsigned int _nDelta_all; /* temperature update time step = nDelta_all * dDelta_gas */
        double dTimeCurrent;  /* current time, s*/
        double _dDelta_gas; /* time step of Fixed-step integral for gas dynamics */
        double _dDelta; /* time step of Fixed-step integral */
        double _dError_eq; /* allocated threshold for deciding whether gas equilibrium is achieved */
        std::string achNodeInFile;     /* Position information of nodes */
        std::string achStateInFile;     /* initial condition file for nodes' physical state */
        std::string achSfactorInFile;     /* self-shadowing S factor file */
        std::string achFfactorInFile;     /* self-heating F factor file  */
        std::string achOutName;        /* output file prefix */
        /* node information */
        unsigned long nNodes;   /* total number of nodes */
        vector3d vDimension; /* dimension of the whole system */
        vector3d vDimensionMax; /* max displacement between two nodes */
        std::vector<fdmnode> NodeVec;
        /* variables */
        unsigned int iStep;  /* current step numbering */
        unsigned int iOutSteps; /* output steps in intervals of dDelta */
        unsigned int iNplayer; /* number of layer nodes */
        bodyphysics bPara;

    public:
        Master(){}
        ~Master(){}
        unsigned int iNthreads = 1; /* number of threads for parallel computing */
        unsigned int TotalSteps() const {return _nSteps;}
        unsigned int OutInterval() const {return _nOutInterval;}
        unsigned int CurrentStep() const {return iStep;}

        /* Parameter initialization */
        bool ParaInit(double df_input, double dr1, double dr2, const double (&df_ice0)[NVOLATILE] );
    
        /* simulation controller */
        bool InitSystem( const string & Filename );
        bool Run(void);
        void Output(void);
};

#endif /* master_hpp */
