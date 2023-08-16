#include "fdmnode.hpp"

istream & operator >> (istream & is, fdmnode & fn) {
    is  >> fn._nTYPE >> fn.vPos;
  return is;
}

/* Initialization of nodes and stars */
bool fdmnode::InitNode(const std::vector<fdmnode> & fnvec, const bodyphysics & bPara, const double & dDelta, const unsigned int & iNplayer) {

    unsigned int i, j, k;
    int iOct;
    double dDis, dvar;
    vector3d vDis;
    std::deque<unsigned int> iOctdeq[8];
    std::deque<double> dOctdeq[8];
    vector3d vDis_all[NNODE_STAR]; /* distance of the center node to each neighbour node */
    double dWeight_all[NNODE_STAR]; /* weighting factor */
    double dMass_all;
    
    /* initialize porosity */
    dT_current_sq = sqrt(dT_current);
    dPorosity = 1.0 - bPara.df_dust;
    dMass_all = bPara.dMass_dust;
    for ( k = 0; k < NVOLATILE; k++ ) {
        df_ice_delta_sum[k] = 0.0;
        dPorosity -= df_ice_current[k];
        dMass_all += df_ice_current[k]*bPara.dRho_ice[k];
    }
    dKappa = bPara.CalConductDust(dT_current)*bPara.dMass_dust/dMass_all;
    for ( k = 0; k < NVOLATILE; k++ ) {
        dKappa += bPara.CalConductIce(dT_current,k)*df_ice_current[k]*bPara.dRho_ice[k]/dMass_all;
    }
    dKappa *= bPara.CalHertz(dPorosity);
    dKappa += 4.0 * sigma * bPara.dR_pore * bPara.dEpsilon * pow(dT_current, 3.0);
    
    
    if ( _nTYPE < 2 ) {
        /* find surface normal direction at current time step */
        norm_surf0 = vPos - fnvec[_nID+iNplayer].vPos;
        dr_lower = Distance3d(fnvec[_nID+iNplayer].vPos, vPos);
        norm_surf0 *= 1.0/dr_lower;
        if ( _nID > iNplayer - 1 )
            dr_upper = Distance3d(fnvec[_nID-iNplayer].vPos, vPos);
        
        return 0;
    }
    
    /* find neighbour for the selection of star */
    for ( i=0; i<fnvec.size(); i++ ) {
        if ( fnvec[i].nID() != _nID ) {
            vDis = fnvec[i].vPos - vPos;
            dDis = Norm3d(vDis);
            iOct = vDis.octant();
            if ( iOctdeq[iOct].size() < NNODE ) {
                iOctdeq[iOct].push_front(i);
                dOctdeq[iOct].push_front(dDis);
            }
            else {
                dvar = 0.0;
                for ( j=0; j<NNODE; j++ ) { /* find the most distance node to replace */
                    if ( dvar < dOctdeq[iOct][j] ) {
                        dvar = dOctdeq[iOct][j];
                        k = j;
                    }
                }
                if ( dDis < dvar ) {
                    iOctdeq[iOct][k] = i;
                    dOctdeq[iOct][k] = dDis;
                }
            }
        }
    }
    
    /* initialize star parameter */
    sNode.iNnode_s= 0;
    for ( i=0; i<NOCTANT; i++ ) {
        for ( j=0; j<iOctdeq[i].size(); j++) {
            sNode._nIDs[sNode.iNnode_s] = iOctdeq[i][j];
            vDis_all[sNode.iNnode_s] = fnvec[iOctdeq[i][j]].vPos - vPos;
            sNode.iNnode_s++;
        }
    }
    sNode.Weighting_Cal( vDis_all, dWeight_all );
    
    /* derive the explicit difference formulae */
    arma::Mat<double> A(9,9,arma::fill::zeros), A_inv;
    arma::Mat<double> II(9,9,arma::fill::ones);
    arma::Mat<double> E_k(9,9,arma::fill::zeros);
    arma::Mat<double> B(9,NNODE_STAR+1,arma::fill::zeros);
    
    for ( k=0; k<sNode.iNnode_s; k++ ) {
        E_k(0,0) = dWeight_all[k]*vDis_all[k].x();
        E_k(1,1) = dWeight_all[k]*vDis_all[k].y();
        E_k(2,2) = dWeight_all[k]*vDis_all[k].z();
        E_k(3,3) = 0.5 * dWeight_all[k]*pow(vDis_all[k].x(),2.0);
        E_k(4,4) = 0.5 * dWeight_all[k]*pow(vDis_all[k].y(),2.0);
        E_k(5,5) = 0.5 * dWeight_all[k]*pow(vDis_all[k].z(),2.0);
        E_k(6,6) = dWeight_all[k]*vDis_all[k].x()*vDis_all[k].y();
        E_k(7,7) = dWeight_all[k]*vDis_all[k].x()*vDis_all[k].z();
        E_k(8,8) = dWeight_all[k]*vDis_all[k].y()*vDis_all[k].z();
        A += E_k * II * E_k;

        for ( i=0; i<9; i++ ) {
            B(i,0) -= dWeight_all[k]*E_k(i,i);
            B(i,k+1) = dWeight_all[k]*E_k(i,i);
        }
    }
    
    cout.precision(15);
    cout.setf(ios::fixed);
    //A.raw_print("A:");/*DEBUG*/
    A_inv = arma::inv(A);
    sNode.E = A_inv.rows(0,5) * B;
    
    /* examine the convergency of the explicit scheme (Eq.60 in Urena et al. 2019) */
    dvar = 0.5 * (bPara.df_dust*bPara.dRho_dust*bPara.CalCapaDust(dT_current, dT_current_sq))/( abs(sNode.E(3,0)+sNode.E(4,0)+sNode.E(5,0)) * dKappa ) ;
    if ( dvar < dDelta ) {
        cout << "ERROR: dDelta (" << dDelta << " s dose not ensure stability (need <" << dvar << " s)!" <<endl;
        return 1;
    }
    if ( abs(A_inv(0,0)) > 1.0e6 ) {
        A_inv.raw_print("A_inv:");
        cout << "ERROR: A matrix is close to singular!" <<endl;
        return 1;
    }

    return 0;
}


/******************************************/
/* asteroid temperature evolution */
/******************************************/

#ifdef USE_GAS

void fdmnode::UpdateGas(const bodyphysics & bPara) {
    unsigned short i_vol;
    dPorosity = 1.0 - bPara.df_dust;
    for ( i_vol=0; i_vol<NVOLATILE; i_vol++ ) {
        dRho_gas_current[i_vol] = dRho_gas_next[i_vol];
        dPorosity -= df_ice_current[i_vol];
    }
}

bool fdmnode::AdvanceInterGas(const std::vector<fdmnode> & fnvec, const bodyphysics & bPara, const double & dt, const double & error_eq, const unsigned int & j_loop, const unsigned int & iNplayer, bool & flag_eq) {
    
    double dM, dM_left, dq_g_t;
    double dRho_sat;
    unsigned int i_vol;
    
    for ( i_vol=0; i_vol<NVOLATILE; i_vol++ ) {
        dRho_sat = bPara.CalSatRhoIce(dT_current, i_vol);
        dM = 3.0 * df_ice_current[i_vol] / bPara.dR_dust * dT_current_sq * bPara.dCoeff_kbmsq[i_vol];
        dM_left = 1.0 + bPara.dAlpha_c * dM * dt / ( 2.0 * dPorosity );
        
        if ( _nTYPE > 1 ) { /* 3D GFDM node */
            arma::Mat<double> mT(NNODE_STAR+1,1,arma::fill::zeros);
            arma::Mat<double> mD(6,1);
            unsigned int i;
            
            /* update neighbor variable list */
            mT(0) = dRho_gas_current[i_vol] * dT_current_sq;
            for ( i=0; i<sNode.iNnode_s; i++ )
                mT(i+1) = fnvec[sNode._nIDs[i]].dRho_gas_current[i_vol]*fnvec[sNode._nIDs[i]].dT_current_sq;
            mD = sNode.E.rows(0,5) * mT;
            
            /* mass flux of gas at current time step */
            vFlux_gas_t[i_vol].SetVec(mD(0),mD(1),mD(2));
            vFlux_gas_t[i_vol] *= bPara.dD_kn_eff[i_vol] * dPorosity;
            vFlux_gas_all[i_vol] += vFlux_gas_t[i_vol];
            
            /* gas mass conservation equation */
            dRho_gas_next[i_vol] = ( dRho_gas_current[i_vol] + dt/dPorosity * ( dM * (bPara.dAlpha_s * dRho_sat - 0.5*bPara.dAlpha_c*dRho_gas_current[i_vol]) - bPara.dD_kn_eff[i_vol] * dPorosity * (mD(3)+mD(4)+mD(5)) ) )/dM_left;
            if ( dRho_gas_next[i_vol] < 0.0) {
                cout << "ERROR: Gas density Rho_gas = " << dRho_gas_next[i_vol] << " ("<< i_vol <<")  is negative (should use smaller dDelta_gas)!" <<endl;
                return 1;
            }
        }
        else { /* 1D FDM node */
            double dFluxrr, dPorosity_lower, dPorosity_upper, dA_self, dA_upper, dA_lower;
            unsigned int iNumLayer;
            
            dPorosity_upper = 0.5 * ( dPorosity + fnvec[_nID-iNplayer].dPorosity );
            dPorosity_lower = 0.5 * ( dPorosity + fnvec[_nID+iNplayer].dPorosity );
            iNumLayer = _nID/iNplayer;
            dA_self = bPara.dArea_1Dlayer[iNumLayer];
            dA_upper = 0.5 * ( dA_self + bPara.dArea_1Dlayer[iNumLayer-1] );
            dA_lower = 0.5 * ( dA_self + bPara.dArea_1Dlayer[iNumLayer+1] );
            dFlux_gas_t[i_vol] = bPara.dD_kn_eff[i_vol] * dPorosity * ( fnvec[_nID+iNplayer].dRho_gas_current[i_vol]*fnvec[_nID+iNplayer].dT_current_sq
                      - fnvec[_nID-iNplayer].dRho_gas_current[i_vol]*fnvec[_nID-iNplayer].dT_current_sq ) / ( dr_upper + dr_lower );
            dFlux_gas_all[i_vol] += dFlux_gas_t[i_vol];
            
            /* gas conservation equation */
            dFluxrr = 2.0 * bPara.dD_kn_eff[i_vol] * ( dPorosity_lower*dA_lower*dr_upper*( fnvec[_nID+iNplayer].dRho_gas_current[i_vol]*fnvec[_nID+iNplayer].dT_current_sq - dRho_gas_current[i_vol]*dT_current_sq )
                                 - dPorosity_upper*dA_upper*dr_lower*( dRho_gas_current[i_vol]*dT_current_sq - fnvec[_nID-iNplayer].dRho_gas_current[i_vol]*fnvec[_nID-iNplayer].dT_current_sq ) ) / ( dA_self * dr_upper * dr_lower * ( dr_upper + dr_lower ) );
            dRho_gas_next[i_vol] = ( dRho_gas_current[i_vol] + dt/dPorosity * ( dM * (bPara.dAlpha_s * dRho_sat - 0.5*bPara.dAlpha_c*dRho_gas_current[i_vol]) - dFluxrr ) )/dM_left;
            if ( dRho_gas_next[i_vol] < 0.0 ) {
                cout << "ERROR: Gas density Rho_gas = " << dRho_gas_next[i_vol] << " ("<< i_vol <<")  is negative (should use smaller dDelta_gas)!" <<endl;
                return 1;
            }
        }

        if ( abs( dRho_gas_next[i_vol] - dRho_gas_current[i_vol] ) > error_eq * abs(dRho_gas_current[i_vol]) )
            flag_eq = false;
        
        /* volume gas mass production/condensation rate at current time step */
        dq_g_t = dM * ( dRho_sat * bPara.dAlpha_s - 0.5*(dRho_gas_current[i_vol]+dRho_gas_next[i_vol]) * bPara.dAlpha_c );
        dq_gas[i_vol] += dq_g_t;
        
        /* ice conservation equation */
        df_ice_delta_sum[i_vol] += dt / bPara.dRho_ice[i_vol] * dq_g_t;
        if ( abs(df_ice_delta_sum[i_vol]) > 1.0e-13 * df_ice_current[i_vol] ) { /* node's ice consumption at one timestep is larger than 1e-13f_i */
            df_ice_current[i_vol] -= df_ice_delta_sum[i_vol];
            df_ice_delta_sum[i_vol] = 0.0;
            if ( df_ice_current[i_vol] < 0.0 )
                df_ice_current[i_vol] = 0.0;
        }
    }
    
    return 0;
}


bool fdmnode::AdvanceBoundaryGas(const std::vector<fdmnode> & fnvec, const bodyphysics & bPara, const double & dt, const double & error_eq, const unsigned int & j_loop, const unsigned int & iNplayer, bool & flag_eq) {
    
    double dM, dQ_sur;
    double dRho_sat;
    unsigned int i_vol;
    
    for ( i_vol=0; i_vol<NVOLATILE; i_vol++ ) {
        dRho_sat = bPara.CalSatRhoIce(dT_current, i_vol);
        dM = dT_current_sq * bPara.dCoeff_kbmsq[i_vol];
        dQ_sur = df_ice_current[i_vol] * dM * ( dRho_sat * bPara.dAlpha_s - dRho_gas_current[i_vol] * bPara.dAlpha_c );
        
        if ( df_ice_current[i_vol] > 1.0e-30 && dRho_gas_current[i_vol] > 0.5 * dRho_sat ) { /* surface ice present */
            dFlux_gas_t[i_vol] =  bPara.dD_kn_eff[i_vol] * 0.5 * ( dPorosity + fnvec[_nID+iNplayer].dPorosity ) * ( fnvec[_nID+iNplayer].dRho_gas_current[i_vol]*fnvec[_nID+iNplayer].dT_current_sq - 0.5 * dRho_sat * dT_current_sq ) / dr_lower;
        }
        else {
            dFlux_gas_t[i_vol] =  bPara.dD_kn_eff[i_vol] * 0.5 * ( dPorosity + fnvec[_nID+iNplayer].dPorosity ) * ( fnvec[_nID+iNplayer].dRho_gas_current[i_vol]*fnvec[_nID+iNplayer].dT_current_sq ) / dr_lower;
        }
        dFlux_gas_all[i_vol] += dFlux_gas_t[i_vol];
        dRho_gas_next[i_vol] = ( dQ_sur - dFlux_gas_t[i_vol] ) / ( ( 1.0 - bPara.dBackflux ) * dM );
        if ( dRho_gas_next[i_vol] < 0.0 ) {
            cout << "WARNING: Gas density Rho_gas = " << dRho_gas_next[i_vol] <<" is negative (the initial boundary condition may be inappropriate)!" <<endl;
            dRho_gas_next[i_vol] = 0.0;
        }
        
        /* ice conservation equation */
        dq_gas[i_vol] += 3.0 * dQ_sur / bPara.dR_dust;
        
        df_ice_delta_sum[i_vol] += dt / bPara.dRho_ice[i_vol] * 3.0 * dQ_sur / bPara.dR_dust;
        if ( abs(df_ice_delta_sum[i_vol]) > 1.0e-13 * df_ice_current[i_vol] ) { /* node's ice consumption at one timestep is larger than 1e-13f_i */
            df_ice_current[i_vol] -= df_ice_delta_sum[i_vol];
            df_ice_delta_sum[i_vol] = 0.0;
            if ( df_ice_current[i_vol] < 0.0 )
                df_ice_current[i_vol] = 0.0;
        }
        
        if ( abs( dRho_gas_next[i_vol] - dRho_gas_current[i_vol] ) > error_eq * abs(dRho_gas_current[i_vol]) )
            flag_eq = false;
    }
    
    return 0;
}

#endif


void fdmnode::UpdateTemp( const std::vector<fdmnode> & fnvec, const bodyphysics & bPara, const unsigned int & iNplayer ) {
    
    double dMass_all = bPara.dMass_dust;
    unsigned short i_vol;
    
    dT_current = dT_next;
    dT_current_sq = sqrt(dT_current);
    
    /* calculate the total mass of the node */
    for ( i_vol = 0; i_vol < NVOLATILE; i_vol++ )
        dMass_all += df_ice_current[i_vol]*bPara.dRho_ice[i_vol];
    
    /* calculate the thermal conductivity */
    dKappa = bPara.CalConductDust(dT_current)*bPara.dMass_dust/dMass_all;
    for ( i_vol = 0; i_vol < NVOLATILE; i_vol++ )
        dKappa += bPara.CalConductIce(dT_current,i_vol)*df_ice_current[i_vol]*bPara.dRho_ice[i_vol]/dMass_all;
    dKappa *= bPara.CalHertz(dPorosity);
    dKappa += 4.0 * sigma * bPara.dR_pore * bPara.dEpsilon * pow(dT_current, 3.0);
    
}


bool fdmnode::AdvanceInterTemp(const std::vector<fdmnode> & fnvec, const bodyphysics & bPara, const double & dt, const unsigned int & iNplayer, const unsigned int & nDelta_all ) {
    
    /* update temperature-dependent properties */
    double dFactor;
    unsigned int i_vol;
    
    dFactor = bPara.dRho_dust*bPara.CalCapaDust(dT_current, dT_current_sq)*bPara.df_dust;
    for ( i_vol=0; i_vol<NVOLATILE; i_vol++ ) {
        dFactor += bPara.dRho_ice[i_vol]*bPara.CalCapaIce(dT_current, i_vol)*df_ice_current[i_vol];
    }
    
    if ( _nTYPE > 1 ) { /* 3D GFDM node */
        unsigned int i;
        arma::Mat<double> mT(NNODE_STAR+1,1,arma::fill::zeros);
        arma::Mat<double> mKappa(NNODE_STAR+1,1,arma::fill::zeros);
        arma::Mat<double> mD(6,1);
        arma::Mat<double> mDK(3,1);
        vector3d vGrad_T, vGrad_K, vGrad_T2, vGrad_TT;
        
        /* update neighbor temperature list */
        mT(0) = dT_current;
        for ( i=0; i<sNode.iNnode_s; i++ ) {
            mT(i+1) = fnvec[sNode._nIDs[i]].dT_current;
            mKappa(i+1) = fnvec[sNode._nIDs[i]].dKappa;
        }
        mD = sNode.E.rows(0,5) * mT;
        mDK = sNode.E.rows(0,2) * mKappa;
        vGrad_T.SetVec(mD(0),mD(1),mD(2));
        vGrad_T2.SetVec(mD(3),mD(4),mD(5));
        vGrad_K.SetVec(mDK(0),mDK(1),mDK(2));
        
        /* energy conservation equation */
        dT_next = dT_current + dt/dFactor * ( dKappa * ( mD(3) + mD(4) + mD(5) + ScalProd3d(vGrad_K,vGrad_T) ) );
        
#ifdef USE_GAS
        for ( i_vol=0; i_vol<NVOLATILE; i_vol++ ) {
            /* energy conservation equation */
            dT_next -= dt/dFactor * ( bPara.dC_gas[i_vol] * ScalProd3d( vFlux_gas_all[i_vol], vGrad_T ) + dq_gas[i_vol] * bPara.CalLatentIce(dT_current,i_vol) );
        }
#endif
    }
    else { /* 1D FDM node */
        double dTrr, dKappa_upper, dKappa_lower, dA_upper, dA_lower, dA_self;
        unsigned int iNumLayer;
        
        dKappa_upper = 0.5 * ( dKappa + fnvec[_nID-iNplayer].dKappa );
        dKappa_lower = 0.5 * ( dKappa + fnvec[_nID+iNplayer].dKappa );
        iNumLayer = _nID/iNplayer;
        dA_self = bPara.dArea_1Dlayer[iNumLayer];
        dA_upper = 0.5 * ( dA_self + bPara.dArea_1Dlayer[iNumLayer-1] );
        dA_lower = 0.5 * ( dA_self + bPara.dArea_1Dlayer[iNumLayer+1] );
        
        dTrr = 2.0 * (dr_upper*dKappa_lower*dA_lower*( fnvec[_nID+iNplayer].dT_current - dT_current ) - dr_lower*dKappa_upper*dA_upper*( dT_current - fnvec[_nID-iNplayer].dT_current)) / (dA_self*dr_upper*dr_lower*(dr_upper+dr_lower));
        
        /* energy conservation equation */
        dT_next = dT_current + dt/dFactor * dTrr;
        
#ifdef USE_GAS
        double dTr;
        dTr = ( fnvec[_nID+iNplayer].dT_current - fnvec[_nID-iNplayer].dT_current ) / ( dr_upper + dr_lower );
        for ( i_vol=0; i_vol<NVOLATILE; i_vol++ ) {
            dT_next -= dt/dFactor * ( bPara.dC_gas[i_vol] * dFlux_gas_all[i_vol] * dTr + dq_gas[i_vol] * bPara.CalLatentIce(dT_current,i_vol) );
        }
#endif
    }
    
    if ( dT_next < 0.0) {
        cout << "ERROR: Temperature = " << dT_next <<" is negative (node No."<< _nID <<")!" <<endl;
        return 1;
    }
    
    return 0;
}


bool fdmnode::AdvanceBoundaryTemp(const std::vector<fdmnode> & fnvec, const bodyphysics & bPara, const double & dt_current, const unsigned int & iNplayer, const unsigned int & nDelta_all) {
    unsigned int j;
    double cos_illumination0, cos_illumination, dT_solution, dKappa_mean;
    arma::vec P(5);
    
    dKappa_mean = 0.5 * ( dKappa + fnvec[_nID+iNplayer].dKappa );
    
    cos_illumination0 = ScalProd3d(norm_surf0, bPara.oEle.vrt_norm);
    if ( cos_illumination0 > 0.0 ) {
#ifdef USE_SELF
        double dS_shadow; /* find self-shadowing partial factor */
        dS_shadow = bPara.oEle.dCloseSolarDirection[0] * bPara.S_factor_all[bPara.oEle.iCloseSolarDirection[0]].S_partial_fraction[_nID] +
                    bPara.oEle.dCloseSolarDirection[1] * bPara.S_factor_all[bPara.oEle.iCloseSolarDirection[1]].S_partial_fraction[_nID] +
                    bPara.oEle.dCloseSolarDirection[2] * bPara.S_factor_all[bPara.oEle.iCloseSolarDirection[2]].S_partial_fraction[_nID];
        cos_illumination = cos_illumination0 * dS_shadow;
#else
        cos_illumination = cos_illumination0;
#endif
    }
    else
        cos_illumination = 0.0;
    
    /* solve boundary condition */
    P(0) = bPara.dEpsilon * sigma;
    P(1) = 0.0; P(2) = 0.0;
    P(3) = dKappa_mean / dr_lower;
    
#ifdef USE_SELF
    if ( iF_view_index.size() > 0 ) {
        double dF_heating_illu, dF_heating_irra, dS_shadow_j;
        dF_heating_illu = 0.0;
        dF_heating_irra = 0.0;
        for ( j=0; j<iF_view_index.size(); j++ ) {
            cos_illumination0 = ScalProd3d(fnvec[iF_view_index[j]].norm_surf0, bPara.oEle.vrt_norm);
            if ( cos_illumination0 > 0.0 ) {
                dS_shadow_j = bPara.oEle.dCloseSolarDirection[0] * bPara.S_factor_all[bPara.oEle.iCloseSolarDirection[0]].S_partial_fraction[iF_view_index[j]] +
                              bPara.oEle.dCloseSolarDirection[1] * bPara.S_factor_all[bPara.oEle.iCloseSolarDirection[1]].S_partial_fraction[iF_view_index[j]] +
                              bPara.oEle.dCloseSolarDirection[2] * bPara.S_factor_all[bPara.oEle.iCloseSolarDirection[2]].S_partial_fraction[iF_view_index[j]];
                dF_heating_illu += dF_view_fraction[j] * cos_illumination0 * dS_shadow_j;
            }
            dF_heating_irra += dF_view_fraction[j] * pow( fnvec[iF_view_index[j]].dT_current, 4.0 );
        }
        P(4) = - dKappa_mean / dr_lower * fnvec[_nID+iNplayer].dT_current - (1.0 - bPara.dA_B) * F_sun / bPara.oEle.drt_sq_AU * ( cos_illumination + bPara.dA_B * dF_heating_illu ) -
                bPara.dEpsilon * sigma * (1.0 - bPara.dA_B) * dF_heating_irra;
    }
    else {
        P(4) = - dKappa_mean / dr_lower * fnvec[_nID+iNplayer].dT_current - (1.0 - bPara.dA_B) * F_sun / bPara.oEle.drt_sq_AU * cos_illumination;
    }
#else
    P(4) = - dKappa_mean / dr_lower * fnvec[_nID+iNplayer].dT_current - (1.0 - bPara.dA_B) * F_sun / bPara.oEle.drt_sq_AU * cos_illumination;
#endif
    
#ifdef USE_GAS
    unsigned int i_vol;
    for ( i_vol=0; i_vol<NVOLATILE; i_vol++ ) {
        P(4) += dq_gas[i_vol] * bPara.dR_dust / 3.0 *  bPara.CalLatentIce(dT_current,i_vol);
    }
#endif
    
    arma::cx_vec poly_roots = roots(P);
    dT_solution = 1e6;
    for ( j=0; j<poly_roots.n_rows; j++ ){ // choose the real and closest value
        if ( abs(std::imag(poly_roots(j))) < 1e-15 )
            if ( abs( std::real(poly_roots(j)) - dT_current ) < dT_solution && std::real(poly_roots(j)) > 0.0 ) {
                dT_solution = abs( std::real(poly_roots(j)) - dT_current );
                dT_next = std::real(poly_roots(j));
            }
    }
    
    if ( dT_next < 0.0) {
        cout << "ERROR: Temperature = " << dT_next <<" is negative (node No."<< _nID <<")!" <<endl;
        return 1;
    }
    
    return 0;
}



