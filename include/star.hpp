#ifndef star_hpp
#define star_hpp

#include "common.hpp"
#include "vector3d.hpp"

class star{ /* contain the whole road used in one simulation */
    friend class fdmnode;
public:
    void Weighting_Cal(const vector3d (&vDis)[NNODE_STAR], double (&dWeight)[NNODE_STAR] ) { /* calculate the weighting function for each node in the star */
        double dDis_norm;
        double dDis_min;
        unsigned int i;
        
        dDis_min = 1.0e8;
        for ( i=0; i<iNnode_s; i++ ) {
            dDis_norm = Norm3d(vDis[i]);
            if ( dDis_norm < dDis_min )
                dDis_min = dDis_norm;
        }
        for ( i=0; i<iNnode_s; i++ ) {
            dDis_norm = Norm3d(vDis[i])/dDis_min;
            dWeight[i] = 1.0/pow(dDis_norm,3.0); /* Eq.(26) in Gavete et al. 2016 */
        }
        
        return;
    }
    private:
        unsigned int _nIDs[NNODE_STAR]; /* neighbour node ID number list */
        unsigned int iNnode_s;  /* actual number of node in the star */
        arma::mat E = arma::Mat<double>(6, NNODE_STAR+1,arma::fill::zeros);
};


#endif /* star_hpp */
