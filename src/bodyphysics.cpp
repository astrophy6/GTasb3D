#include "bodyphysics.hpp"


void orbit::coe2rv( const double & dtime, const bodyphysics & bPara ) {
    double dMt, dEt, dW;
    vector3d vp, vq, vrt_norm_inertia, vrt_norm_body_rot;
    
    dMt = sqrt( MU_sun /pow(da,3)) * dtime + dM0;

    vp.SetVec(cos(dOmega)*cos(domega)-sin(dOmega)*sin(domega)*cos(di),sin(dOmega)*cos(domega)+cos(dOmega)*sin(domega)*cos(di), sin(domega)*sin(di));

    vq.SetVec(-cos(dOmega)*sin(domega)-sin(dOmega)*cos(domega)*cos(di),-sin(dOmega)*sin(domega)+cos(dOmega)*cos(domega)*cos(di), cos(domega)*sin(di));

    if ( de < 1.0 ) {
        dEt = M2E(dMt);
        vrt_norm_inertia = da*(cos(dEt) - de) * vp + da*sqrt(1.0-de*de) * sin(dEt) * vq;
    }
    else
        cout << "Cann't process hyperbolic orbit for now!" << endl;
    drt_sq_AU = MagSq3d(vrt_norm_inertia);
    vrt_norm_inertia *= -1.0/sqrt(drt_sq_AU);
    drt_sq_AU /= AU*AU;
    E2f(dEt); /* update the true anomaly */
    
    /* transfer solar direction to the body-fixed frame */
    dW = bPara.dOmega_rot * dtime;
    vrt_norm_body_rot.SetVec(ScalProd3d(vrt_norm_inertia, bPara.mRotate[0]), ScalProd3d(vrt_norm_inertia, bPara.mRotate[1]), ScalProd3d(vrt_norm_inertia, bPara.mRotate[2]));
    vrt_norm.SetVec(vrt_norm_body_rot.x()*cos(dW)+vrt_norm_body_rot.y()*sin(dW), -vrt_norm_body_rot.x()*sin(dW)+vrt_norm_body_rot.y()*cos(dW), vrt_norm_body_rot.z());
    
#ifdef USE_SELF
    /* find the solar direction for the S factor matrix */
    unsigned int ii;
    double dDiff, dDiff_min1, dDiff_min2, dDiff_min3;
    
    dDiff_min1 = 100.0; dDiff_min2 = 100.0; dDiff_min3 = 100.0;
    for ( ii=0; ii<bPara.iS_number; ii++ ) {
        dDiff = MagSq3d( bPara.S_factor_all[ii].vSunVector + vrt_norm);
        if ( dDiff < dDiff_min1 ) {
            dDiff_min3 = dDiff_min2;
            iCloseSolarDirection[2] = iCloseSolarDirection[1];
            dDiff_min2 = dDiff_min1;
            iCloseSolarDirection[1] = iCloseSolarDirection[0];
            dDiff_min1 = dDiff;
            iCloseSolarDirection[0] = ii;
        }
        else if ( dDiff < dDiff_min2 ) {
            dDiff_min3 = dDiff_min2;
            iCloseSolarDirection[2] = iCloseSolarDirection[1];
            dDiff_min2 = dDiff;
            iCloseSolarDirection[1] = ii;
        }
        else if ( dDiff < dDiff_min3 ) {
            dDiff_min3 = dDiff;
            iCloseSolarDirection[2] = ii;
        }
        dDiff = 2.0 * ( dDiff_min1 + dDiff_min2 + dDiff_min3 );
        dCloseSolarDirection[0] = ( dDiff_min2 + dDiff_min3 ) / dDiff;
        dCloseSolarDirection[1] = ( dDiff_min1 + dDiff_min3 ) / dDiff;
        dCloseSolarDirection[2] = ( dDiff_min1 + dDiff_min2 ) / dDiff;
    }
#endif
}
double bodyphysics::CalHertz(const double & dPoro) const { /* Shoshany, Prialnik & Podolak, 2002 */
    double dp0 = 1.0 - sqrt( 1.0 - dPoro );
    return pow( 1.0 - dp0/0.7, pow( 4.1*dp0 + 0.22, 2.0 ) );
}

double bodyphysics::CalCapaDust(const double & dTemp, const double & dTemp_sq) const { /* Robie et al. 1982 */
    if ( dTemp < 50.0 )
        return 0.05922*dTemp*dTemp - 6.173 * dTemp +41.41 * dTemp_sq + 92.77/dTemp_sq - 103.1344;
    else if ( dTemp < 380.0 )
        return -6.361*dTemp + 316.7*dTemp_sq + 10160.0/dTemp_sq - 3316.0;
    else //if ( dTemp < 1800 )
        return -1.5945e-04 * dTemp * dTemp + 0.6213 * dTemp + 6.0130e+03/dTemp_sq - 2.6366e+07 / (dTemp*dTemp) + 622.6967;
}
double bodyphysics::CalConductDust(const double & dTemp) const {
    return 0.5;
}
double bodyphysics::CalCapaIce(const double & dTemp, const unsigned int nVol) const {
    switch(nVol) {
        case 0: /* H2O */
            return 7.49 * dTemp + 90.0;
            break;
        case 1: /* CO */
            return 35.7 * dTemp - 187.0;
            break;
    }
    return 0;
} /* Opeil et al 2010 */
double bodyphysics::CalLatentIce(const double & dTemp, const unsigned int nVol) const {
    switch(nVol) {
        case 0: /* H2O */
            return 2.867127e6 - 1.10807e3*dTemp;
            break;
        case 1: /* CO */
            return 0.227e6;
            break;
    }
    return 0;
} /* Orosei et al. 1995 */
double bodyphysics::CalSatRhoIce(const double & dTemp, const unsigned int nVol) const {
    switch(nVol) {
        case 0: /* H2O */
            return dCoeff_mkb_sat[0] * exp(-6141.667/dTemp) / dTemp;
            break;
        case 1: /* CO */
            return dCoeff_mkb_sat[1] * exp(-764.16/dTemp) / dTemp;;
            break;
    }
    return 0;
} /* Fanale & Salvail 1984 */
double bodyphysics::CalConductIce(const double & dTemp, const unsigned int nVol) const {
    switch(nVol) {
        case 0: /* H2O */
            return 567.0 / dTemp;
            break;
        case 1: /* CO */
            return 567.0 / dTemp;
            break;
    }
    return 0;
} /* Seiferlin et al 1995 */
double bodyphysics::CalGasProduct(const double & dRho_g, const double & df_i, const double & dTemp, const double & dTemp_sq, const unsigned int nVol) const { /* Davidsson 2021 */
    double dRho_sat = CalSatRhoIce(dTemp,nVol);
    return 3.0*df_i/dR_pore * ( dRho_sat*dAlpha_s - dRho_g*dAlpha_c ) * dTemp_sq * dCoeff_kbmsq[nVol];
}
