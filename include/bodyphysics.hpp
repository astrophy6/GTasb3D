#include "common.hpp"
#include "vector3d.hpp"

#ifndef bodyphysics_hpp
#define bodyphysics_hpp

class Master;

#ifdef USE_SELF
class S_factor {
    friend class fdmnode;
    friend class Master;
    friend class orbit;
private:
    vector3d vSunVector;
    std::vector<double> S_partial_fraction;  /* self-shadowing partial fraction for given solar direction */
};
#endif

class bodyphysics;

class orbit { /* orbital elements and function */
    friend class fdmnode;
    friend class Master;
private:
    double da;  /* semi-major axis, m */
    double de;  /* eccentricity */
    double di;  /* inclination, rad */
    double dOmega;  /* ascending node, rad */
    double domega;  /* argument of periapsis, rad */
    double df0;  /* initial true anomaly, rad */
    double dM0;  /* initial mean anomaly, rad */
    double df;  /* true anomaly, rad */
    double drt_sq_AU;  /* distance at a given time, AU */
    vector3d vrt_norm;  /* normalized position at a given time */
#ifdef USE_SELF
        unsigned int iCloseSolarDirection[3]; /* index of the 3 closest solar direction */
        double dCloseSolarDirection[3]; /* scaling fraction of the 3 closest solar direction */
#endif
public:
    double f() {return df;}
    void f2M() {
        double dvar1, dvar2;
        double dE0=0.0;  /* eccentric anomaly, rad */
        
        dvar1 = floor(0.25*df0/M_PI);
        dvar2 = 0.5 * df0 - 2.0 * dvar1 * M_PI;

        if ( dvar2 >= 0.0 && dvar2 < 0.5*M_PI )
            dE0 = 2.0*acos(sqrt((((de+cos(df0))/(1.0 + de*cos(df0)))+1.0)*0.5)) + 4.0*dvar1*M_PI;
        else if ( dvar2 >= 0.5*M_PI && dvar2 < M_PI )
            dE0 = 2.0*acos(-sqrt((((de+cos(df0))/(1.0 + de*cos(df0)))+1.0)*0.5)) + 4.0*dvar1*M_PI;
        else if ( dvar2 >= M_PI && dvar2 < 1.5*M_PI)
            dE0 = 4.0*M_PI-2.0*acos(-sqrt((((de+cos(df))/(1.0 + de*cos(df0)))+1.0)*0.5)) + 4.0*dvar1*M_PI;
        else if ( dvar2 >= 1.5*M_PI && dvar2 < 2.0*M_PI)
            dE0 = 4.0*M_PI-2.0*acos(sqrt((((de+cos(df))/(1.0 + de*cos(df0)))+1.0)*0.5)) + 4.0*dvar1*M_PI;
        dM0 = dE0 - de * sin(dE0);
    }
    void E2f(const double & dEt) {
        double dvar1, dvar2;
        
        dvar1 = floor(0.25*dEt/M_PI);
        dvar2 = 0.5 * dEt - 2.0 * dvar1 * M_PI;

        if ( dvar2 >= 0.0 && dvar2 < 0.5*M_PI )
            df = 2.0*acos(sqrt((((cos(dEt)-de)/(1.0-de*cos(dEt)))+1.0)*0.5)) + 4.0*dvar1*M_PI;
        else if ( dvar2 >= 0.5*M_PI && dvar2 < M_PI )
            df = 2.0*acos(-sqrt((((cos(dEt)-de)/(1.0-de*cos(dEt)))+1.0)*0.5)) + 4.0*dvar1*M_PI;
        else if ( dvar2 >= M_PI && dvar2 < 1.5*M_PI)
            df = 4.0*M_PI-2.0*acos(-sqrt((((cos(dEt)-de)/(1.0-de*cos(dEt)))+1.0)*0.5)) + 4.0*dvar1*M_PI;
        else if ( dvar2 >= 1.5*M_PI && dvar2 < 2.0*M_PI)
            df = 4.0*M_PI-2.0*acos(sqrt((((cos(dEt)-de)/(1.0-de*cos(dEt)))+1.0)*0.5)) + 4.0*dvar1*M_PI;
    }
    double E2M(const double & dEt) {
        double dMt;
        dMt = dEt - de * sin(dEt);
        return dMt;
    }
    double M2E(const double & dMt) {
        double dEt, fe, dfe;
        dEt = dMt;
        if ( abs(de) < 1e-13 )
            return dEt;
        fe = 1.0;
        while ( abs(fe) > 1e-13 ) {
            fe = dEt - de * sin(dEt) - dMt;
            dfe = 1.0 - de * cos(dEt);
            dEt = dEt - fe/dfe;
        }
        return dEt;
    }
    void coe2rv(const double & dtime, const bodyphysics & bPara);
};


class bodyphysics {  /* body physical parameters */
    friend class fdmnode;
    friend class Master;
    friend class orbit;
    private:
        double dP_rot;  /* rotational period, s */
        double dOmega_rot;  /* rotation rate, rad/s */
        double dObliquity;  /* obliquity of spin axis, rad */
        double dSolstice;  /* solstice true anomaly, rad */
        double df_dust; /* porosity */
        double dRho_ice[3]; /* compact ice density, kg/m^3 */
        double dRho_dust; /* compact dust density, kg/m^3 */
        double dMass_dust; /* dust density, kg/m^3 */
        double dA_B;  /* Bond albedo */
        double dEpsilon; /* bolometric emissivity */
        double dC_gas[NVOLATILE]; /* specific heat capacity (h2o, co, co2) */
        double dTortuosity; /* channel tortuosity */
        double dL_pore;  /* pore length, m */
        double dR_pore;  /* pore radius, m */
        double dR_dust;  /* dust particle radius, m */
        double dBackflux; /* surface backflux fraction */
        double dAlpha_s;  /* sublimation coefficient */
        double dAlpha_c;  /* condensation coefficient */
        double dD_kn_eff[NVOLATILE]; /* Knudsen diffusion coefficient (h2o, co, co2) */
        double dCoeff_kbmsq[3];  /* coefficient */
        double dCoeff_mkb_sat[3];  /* coefficient */
        orbit oEle; /* orbital elements */
        vector3d mRotate[3];  /* rotation matrix */
        std::vector<double> dArea_1Dlayer;  /* total area of each layer */
#ifdef USE_SELF
        unsigned int iS_number;  /* total solar direction number dataset */
        std::vector<S_factor> S_factor_all;  /* self-shadowing partial fraction */
#endif
    public:
        double CalHertz(const double & dPoro) const; /* Shoshany, Prialnik & Podolak, 2002 */
        double CalCapaDust(const double & dTemp, const double & dTemp_sq) const; /* Robie et al. 1982 */
        double CalConductDust(const double & dTemp) const; /* Opeil et al 2010 */
    
        double CalCapaIce(const double & dTemp, const unsigned int nVol) const; /* Opeil et al 2010 */
        double CalLatentIce(const double & dTemp, const unsigned int nVol) const; /* Orosei et al. 1995 */
        double CalSatRhoIce(const double & dTemp, const unsigned int nVol) const; /* Fanale & Salvail 1984 */
        double CalConductIce(const double & dTemp, const unsigned int nVol) const; /* Seiferlin et al 1995 */
        double CalGasProduct(const double & dRho_g, const double & df_i, const double & dTemp, const double & dTemp_sq, const unsigned int nVol) const;
    
};
#endif /* bodyphysics_hpp */
