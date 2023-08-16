#ifndef _COMMON_H
#define _COMMON_H

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <deque>
#include <algorithm>
#include <armadillo>
#include <omp.h>

using namespace std;
using namespace arma;

/* define constants */
namespace myconstants
{
    const unsigned int NVOLATILE = 1;  /* number of volatile spieces considered in the simulation, 1: H2O, 2: CO 3: CO2*/
    const unsigned int NNODE = 3; /* number of nodes in each octant (eight octant criterion) */
    const unsigned int NOCTANT = 8; /* number of octant (3D: 8) */
    const unsigned int NNODE_STAR = NNODE*NOCTANT;   /* number of nodes in each star (eight octant criterion) */
    const double G = 6.67430e-11; /* gravitational constant, m */
    const double AU = 1.495978707e11;  /* astronomical unit, m */
    const double SID_YR = 31558149.7635000; /* one sidereal year in seconds (2000.0) */
    const double M_Sun = 1.9884e30;  /* solar mass, kg */
    const double MU_sun = G*M_Sun; /* gravitaional parameter */
    const double F_sun = 1370.0;  /* integrated stellar flux at 1 au, W/m^2 */
    const double sigma = 5.670374419e-8;  /* Stefan-Boltzmann constant, W/(m^2.K^4) */
    const double k_B = 1.380649e-23;  /* Boltzmann constant, J/K */
    const double m_h = 1.6738e-27;  /* H molecule mass, kg */
    const double m_h2o = 2.9915e-26;  /* H2O molecule mass, kg */
    const double m_co = 4.6514e-26;  /* CO molecule mass, kg */
    const double m_co2 = 7.3082e-26;  /* CO2 molecule mass, kg */
    const double m_gas[3] = {m_h2o, m_co, m_co2};  /* H molecule mass, kg */
    const double sigma_J = 1.4e-18;  /* H2O-H2O collision cross-section, m (Mammoli, et al. 2016, Scientific reports) */
    const double gamma = 4.0/3.0;  /* specific heat ratio */
    const double zeta = 3.0;  /* specific heat ratio */
    const double Tfree = 20.0;  /* free temperature */
}

using namespace myconstants;


/* define compile setup */
#define USE_GAS  /* turn on gas and ice activity */
//#define USE_SELF /* include self-shadow and self-heating */
#define SEETIME      /* set the parameter when need to estimate the computation time */

#endif
