#ifndef fdmnode_hpp
#define fdmnode_hpp

#include "common.hpp"
#include "star.hpp"
#include "bodyphysics.hpp"

class fdmnode{ /* contain the whole node used in one simulation */
    friend istream & operator >> (istream & is, fdmnode & fn);
    private:
        unsigned long _nID; /* node's serial ID number (constant)  */
        unsigned int _nTYPE; /* node type: 0, boundary; 1, interior 1D; Others, interior 3D */
        double dq_gas[NVOLATILE]; /* mean gas producation rate during nDelta_all */
        double dFlux_gas_all[NVOLATILE]; /* mean gas diffusion mass flux rate during nDelta_all (1D) */
        double dFlux_gas_t[NVOLATILE]; /* mean gas diffusion mass flux rate during nDelta_all (1D) */
        double dT_current; /* temperature at the current time step, K */
        double dT_current_sq; /* square root of temperature at the current time step */
        double dT_next; /* temperature at the next time step */
        double dRho_gas_current[NVOLATILE]; /* gas mass density at the current time step, kg/m^3 */
        double dRho_gas_next[NVOLATILE]; /* gas mass density at the next time step, kg/m^3 */
        double df_ice_current[NVOLATILE]; /* ice fraction at the current time step */
        double dPorosity; /* local porosity of the node */
        double dKappa; /* local heat conductivity */
        double df_ice_delta_sum[NVOLATILE]; /* accumulative ice consumption within a time period to ensure bf_ice_delta_sum >1e-12f_i for numerical accuracy  */
        double dr_upper, dr_lower;  /* spacing between the 1D node */
        vector3d vPos; /* position */
        vector3d norm_surf0; /* surface normal direction for 1D node */
        vector3d vFlux_gas_all[NVOLATILE]; /* mean gas diffusion mass flux rate during nDelta_all (3D) */
        vector3d vFlux_gas_t[NVOLATILE];
        star sNode;  /* star information */
    public:
#ifdef USE_SELF
        std::vector<unsigned int> iF_view_index;  /* self-heating corresponding node */
        std::vector<double> dF_view_fraction;  /* self-heating view fraction */
#endif
        void UpdateTemp(const std::vector<fdmnode> & fnvec, const bodyphysics & bPara, const unsigned int & iNplayer);
        void SetID(unsigned long i) { _nID = i; }
        unsigned long nID() const { return _nID; }
        bool InitNode(const std::vector<fdmnode> & fnvec, const bodyphysics & bPara, const double & dDelta, const unsigned int & iNplayer);
        bool AdvanceInterTemp(const std::vector<fdmnode> & fnvec, const bodyphysics & bPara, const double & dt, const unsigned int & iNplayer, const unsigned int & nDelta_all);
        bool AdvanceBoundaryTemp(const std::vector<fdmnode> & fnvec, const bodyphysics & bPara, const double & dt_current, const unsigned int & iNplayer, const unsigned int & nDelta_all);
#ifdef USE_GAS
        bool AdvanceInterGas(const std::vector<fdmnode> & fnvec, const bodyphysics & bPara, const double & dt, const double & error_eq, const unsigned int & j, const unsigned int & iNplayer, bool & flag_eq);
        bool AdvanceBoundaryGas(const std::vector<fdmnode> & fnvec, const bodyphysics & bPara, const double & dt, const double & error_eq, const unsigned int & j, const unsigned int & iNplayer, bool & flag_eq);
        void UpdateGas(const bodyphysics & bPara);
        void InitGas() {
            unsigned int i_vol;
            for ( i_vol=0; i_vol<NVOLATILE; i_vol++ ) {
                dq_gas[i_vol] = 0.0;
                vFlux_gas_all[i_vol].SetVec(0.0, 0.0, 0.0);
                dFlux_gas_all[i_vol] = 0.0;
            }
        }
        void AdvanceEquil( const bodyphysics & bPara, const double & dDelta_gas, const unsigned int & j_eq, const unsigned int & _nDelta_all ) {
            double dTimeExtrap, dFactor, dtemp, dq_g_t;
            double dOverDelta = 1.0 / _nDelta_all;
            unsigned int i_vol;
            
            for ( i_vol=0; i_vol<NVOLATILE; i_vol++ ) {
                dtemp = df_ice_current[i_vol];
                dTimeExtrap = ( _nDelta_all - j_eq - 1 ) * dDelta_gas;
                dFactor = 3.0 / bPara.dR_dust * ( bPara.CalSatRhoIce(dT_current,i_vol) * bPara.dAlpha_s - dRho_gas_current[i_vol] * bPara.dAlpha_c ) * bPara.dCoeff_kbmsq[i_vol] * dT_current_sq;
                df_ice_current[i_vol] *= exp( - dFactor / bPara.dRho_ice[i_vol] * dTimeExtrap );
                dq_g_t = 0.5 * ( df_ice_current[i_vol] + dtemp ) * dFactor;
                dq_gas[i_vol] += dq_g_t * ( _nDelta_all - j_eq - 1 );
                dq_gas[i_vol] *= dOverDelta;
                if ( _nTYPE > 1 ) { /* 3D interior */
                    vFlux_gas_all[i_vol] += vFlux_gas_t[i_vol] * ( _nDelta_all - j_eq - 1 );
                    vFlux_gas_all[i_vol] *= dOverDelta;
                }
                else { /* 1D interior */
                    dFlux_gas_all[i_vol] += dFlux_gas_t[i_vol] * ( _nDelta_all - j_eq - 1 );
                    dFlux_gas_all[i_vol] *= dOverDelta;
                }
            }
        }
#endif
        vector3d abs3d( const fdmnode & fn) {return Abs3d(vPos,fn.vPos);}
        double Dis3d( const fdmnode & fn) {return Distance3d(vPos,fn.vPos);}
        unsigned int & type() { return _nTYPE; }
        unsigned int type() const { return _nTYPE; }
        double & x() {return vPos.x();}
        double x() const {return vPos.x();}
        double & y() {return vPos.y();}
        double y() const {return vPos.y();}
        double & z() {return vPos.z();}
        double z() const {return vPos.z();}
        double & Temp() {return dT_current;}
        double Temp() const {return dT_current;}
        double & Rho_gas(unsigned int i) {return dRho_gas_current[i];}
        double Rho_gas(unsigned int i) const {return dRho_gas_current[i];}
        double & f_ice(unsigned int i) {return df_ice_current[i];}
        double f_ice(unsigned int i) const {return df_ice_current[i];}
        double & q_gas(unsigned int i) {return dq_gas[i];}
        double q_gas(unsigned int i) const {return dq_gas[i];}
        double & Flux_gas1D(unsigned int i) {return dFlux_gas_all[i];}
        double Flux_gas1D(unsigned int i) const {return dFlux_gas_all[i];}
        double & Flux_gas3Dx(unsigned int i) {return vFlux_gas_all[i].x();}
        double Flux_gas3Dx(unsigned int i) const {return vFlux_gas_all[i].x();}
        double & Flux_gas3Dy(unsigned int i) {return vFlux_gas_all[i].y();}
        double Flux_gas3Dy(unsigned int i) const {return vFlux_gas_all[i].y();}
        double & Flux_gas3Dz(unsigned int i) {return vFlux_gas_all[i].z();}
        double Flux_gas3Dz(unsigned int i) const {return vFlux_gas_all[i].z();}
};

#endif /* fdmnode_hpp */
