#ifndef _SIMULATION_H
#define _SIMULATION_H

#include "master.hpp"
#include <map>
#include <dirent.h>


bool is_number(const std::string& str) {
    for (char const &c : str)
        if ( std::isdigit(c) == 0 )
            return false;
    return true;
}

bool isDouble(const std::string& str) {
    std::istringstream iss(str);
    double d;
    
    iss >> noskipws >> d;
    
    return iss.eof() && !iss.fail();
}


/* Input Data & Initialization */
bool Master::InitSystem(const string & Filename){

    std::map< string, int > TokenMap;
    ifstream finput;
    string sWord, achDataSubPath;
    double df_input;
    unsigned long ii, jj;
    bool flag_error = false;

    /* initializations for tokens */
    TokenMap["#"] = 0;
    TokenMap["!"] = 0;
    TokenMap["nDigits"] = 1;
    TokenMap["nSteps"] = 2;
    TokenMap["iStartStep"] = 3;
    TokenMap["dDelta_gas"] = 4;
    TokenMap["nDelta_all"] = 5;
    TokenMap["achDataSubPath"] = 6;
    TokenMap["achNodeInFile"] = 7;
    TokenMap["achStateInFile"] = 8;
    TokenMap["achOutName"] = 9;
    TokenMap["nOutInterval"] = 11;
    TokenMap["dP_rot"] = 20;
    TokenMap["dObliquity"] = 21;
    TokenMap["dSolstice"] = 22;
    TokenMap["df_dust"] = 23;
    TokenMap["dRho_h2o_ice"] = 24;
    TokenMap["dRho_co_ice"] = 25;
    TokenMap["dRho_co2_ice"] = 26;
    TokenMap["dRho_dust"] = 29;
    TokenMap["dA_B"] = 30;
    TokenMap["dEpsilon"] = 31;
    TokenMap["dTortuosity"] = 40;
    TokenMap["dL_pore"] = 41;
    TokenMap["dR_pore"] = 42;
    TokenMap["dR_dust"] = 43;
    TokenMap["dBackflux"] = 45;
    TokenMap["dAlpha_s"] = 46;
    TokenMap["dAlpha_c"] = 47;
    TokenMap["da"] = 50;
    TokenMap["de"] = 51;
    TokenMap["di"] = 52;
    TokenMap["dOmega"] = 53;
    TokenMap["domega"] = 54;
    TokenMap["df0"] = 55;
    TokenMap["achSfactorInFile"] = 60;
    TokenMap["achFfactorInFile"] = 61;
    TokenMap["bReadLatest"] = 70;
    TokenMap["dError_eq"] = 80;


    /********************/
    /**   Input Data   **/
    /********************/
    /* general information */
    finput.open( Filename.c_str(), ios::in );
    cout << "Read simulation general parameters from " << Filename << endl;
    if ( finput.is_open() ) {
        while ( !finput.eof() ) {
            if ( finput.peek() == '#' || finput.peek() == '!') /* skip comment */
                finput.ignore( 1000, '\n' );
            else {
                finput >> sWord;
                switch ( TokenMap[sWord] ) {
                    case 0:
                        finput.ignore( 1000, '\n' ); /* skip comment */
                        break;
                    case 70:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( is_number(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> _bReadLatest;
                            cout << "bReadLatest: " << _bReadLatest << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for bReadLatest!" << endl;
                            return 1;
                        }
                        break;
                    case 1:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( is_number(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> _nDigits;
                            cout << "nDigits: " << _nDigits << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for nDigits!" << endl;
                            return 1;
                        }
                        break;
                    case 2:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( is_number(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> _nSteps;
                            cout << "nSteps: " << _nSteps << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for nSteps!" << endl;
                            return 1;
                        }
                        break;
                    case 3:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( is_number(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> _iStartStep;
                            cout << "iStartStep: " << _iStartStep << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for iStartStep!" << endl;
                            return 1;
                        }
                        break;
                    case 4:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> _dDelta_gas;
                            cout << "dDelta_gas: " << _dDelta_gas << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dDelta_gas!" << endl;
                            return 1;
                        }
                        break;
                    case 5:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( is_number(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> _nDelta_all;
                            cout << "nDelta_all: " << _nDelta_all << endl;
                            _dDelta = _dDelta_gas * _nDelta_all;
                            cout << "dDelta: " << _dDelta << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for nDelta_all!" << endl;
                            return 1;
                        }
                        break;
                    case 80:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> _dError_eq;
                            cout << "dError_eq: " << _dError_eq << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dError_eq!" << endl;
                            return 1;
                        }
                        break;
                    case 6:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        achDataSubPath = sWord;
                        cout << "achDataSubPath: " << achDataSubPath << endl;
                        break;
                    case 7:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        {
                            stringstream sstream;
                            sstream << achDataSubPath << "/"  << sWord;
                            sstream >> achNodeInFile;
                        }
                        cout << "achNodeInFile: " << achNodeInFile << endl;
                        break;
                    case 8:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        {
                            stringstream sstream;
                            sstream << achDataSubPath << "/"  << sWord;
                            sstream >> achStateInFile;
                        }
                        cout << "achStateInFile: " << achStateInFile << endl;
                        break;
                    case 9:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> achOutName;
                        }
                        cout << "achOutName: " << achOutName << endl;
                        break;
                    case 11:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( is_number(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> _nOutInterval;
                            cout << "nOutInterval: " << _nOutInterval << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for nOutInterval!" << endl;
                            return 1;
                        }
                        break;
                    case 20:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.dP_rot;
                            cout << "dP_rot: " <<  bPara.dP_rot << " hr" << endl;
                            if ( bPara.dP_rot < 0.0 ) {
                                cout << Filename << ": Unvalid argument set for dP_rot!" << endl;
                                return 1;
                            }
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dP_rot!" << endl;
                            return 1;
                        }
                        break;
                    case 21:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.dObliquity;
                            cout << "dObliquity: " <<  bPara.dObliquity << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dObliquity!" << endl;
                            return 1;
                        }
                        break;
                    case 22:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.dSolstice;
                            cout << "dSolstice: " <<  bPara.dSolstice << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dSolstice!" << endl;
                            return 1;
                        }
                        break;
                    case 23:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.df_dust;
                            cout << "df_dust: " <<  bPara.df_dust << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for df_dust!" << endl;
                            return 1;
                        }
                        break;
                    case 24:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.dRho_ice[0];
                            cout << "dRho_h2o_ice: " <<  bPara.dRho_ice[0] << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dRho_h2o_ice!" << endl;
                            return 1;
                        }
                        break;
                    case 25:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.dRho_ice[1];
                            cout << "dRho_co_ice: " <<  bPara.dRho_ice[1] << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dRho_co_ice!" << endl;
                            return 1;
                        }
                        break;
                    case 26:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.dRho_ice[2];
                            cout << "dRho_co2_ice: " <<  bPara.dRho_ice[2] << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dRho_co2_ice!" << endl;
                            return 1;
                        }
                        break;
                    case 29:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.dRho_dust;
                            cout << "dRho_dust: " <<  bPara.dRho_dust << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dRho_dust!" << endl;
                            return 1;
                        }
                        break;
                    case 30:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.dA_B;
                            cout << "dA_B: " <<  bPara.dA_B << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dA_B!" << endl;
                            return 1;
                        }
                        break;
                    case 31:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.dEpsilon;
                            cout << "dEpsilon: " <<  bPara.dEpsilon << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dEpsilon!" << endl;
                            return 1;
                        }
                        break;
                    case 40:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.dTortuosity;
                            cout << "dTortuosity: " <<  bPara.dTortuosity << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dTortuosity!" << endl;
                            return 1;
                        }
                        break;
                    case 41:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.dL_pore;
                            cout << "dL_pore: " <<  bPara.dL_pore << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dL_pore!" << endl;
                            return 1;
                        }
                        break;
                    case 42:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.dR_pore;
                            cout << "dR_pore: " <<  bPara.dR_pore << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dR_pore!" << endl;
                            return 1;
                        }
                        break;
                    case 43:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.dR_dust;
                            cout << "dR_dust: " <<  bPara.dR_dust << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dR_dust!" << endl;
                            return 1;
                        }
                        break;
                    case 45:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.dBackflux;
                            cout << "dBackflux: " <<  bPara.dBackflux << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dBackflux!" << endl;
                            return 1;
                        }
                        break;
                    case 46:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.dAlpha_s;
                            cout << "dAlpha_s: " <<  bPara.dAlpha_s << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dAlpha_s!" << endl;
                            return 1;
                        }
                        break;
                    case 47:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.dAlpha_c;
                            cout << "dAlpha_c: " <<  bPara.dAlpha_c << endl;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dAlpha_c!" << endl;
                            return 1;
                        }
                        break;
                    case 50:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.oEle.da;
                            cout << "da: " <<  bPara.oEle.da << endl;
                            if ( bPara.oEle.da < 0.0 ) {
                                cout << Filename << ": Unvalid argument set for da!" << endl;
                                return 1;
                            }
                            bPara.oEle.da *= AU;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for da!" << endl;
                            return 1;
                        }
                        break;
                    case 51:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.oEle.de;
                            cout << "de: " <<  bPara.oEle.de << endl;
                            if ( bPara.oEle.de < 0.0 ) {
                                cout << Filename << ": Unvalid argument set for de!" << endl;
                                return 1;
                            }
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for de!" << endl;
                            return 1;
                        }
                        break;
                    case 52:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.oEle.di;
                            cout << "di: " <<  bPara.oEle.di << endl;
                            if ( bPara.oEle.di < 0.0 || bPara.oEle.di > M_PI ) {
                                cout << Filename << ": Unvalid argument set for di!" << endl;
                                return 1;
                            }
                            bPara.oEle.di *= M_PI/180.0;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for di!" << endl;
                            return 1;
                        }
                        break;
                    case 53:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.oEle.dOmega;
                            cout << "dOmega: " <<  bPara.oEle.dOmega << endl;
                            bPara.oEle.dOmega *= M_PI/180.0;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for dOmega!" << endl;
                            return 1;
                        }
                        break;
                    case 54:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.oEle.domega;
                            cout << "domega: " <<  bPara.oEle.domega << endl;
                            bPara.oEle.domega *= M_PI/180.0;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for domega!" << endl;
                            return 1;
                        }
                        break;
                    case 55:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        if( isDouble(sWord) ) {
                            stringstream sstream;
                            sstream << sWord;
                            sstream >> bPara.oEle.df0;
                            cout << "df0: " <<  bPara.oEle.df0 << endl;
                            if ( bPara.oEle.de*cos(bPara.oEle.df0) < -1.0) {
                                cout << Filename << ": Unvalid argument set for df0!" << endl;
                                return 1;
                            }
                            bPara.oEle.df0 *= M_PI/180.0;
                        }
                        else {
                            cout << Filename << ": Unvalid argument set for df0!" << endl;
                            return 1;
                        }
                        break;
                    case 60:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        {
                            stringstream sstream;
                            sstream << achDataSubPath << "/"  << sWord;
                            sstream >> achSfactorInFile;
                        }
                        cout << "achSfactorInFile: " << achSfactorInFile << endl;
                        break;
                    case 61:
                        finput >> sWord;
                        if( sWord == "=" ) finput >> sWord;
                        {
                            stringstream sstream;
                            sstream << achDataSubPath << "/"  << sWord;
                            sstream >> achFfactorInFile;
                        }
                        cout << "achFfactorInFile: " << achFfactorInFile << endl;
                        break;
                    default:
                        cout << Filename << ": Unrecognized argument " << sWord << endl;
                        return 1;
                }
            }
        }
        finput.close();
        finput.clear();
    }
    else {
        cout << Filename << ": Unable to open file." << endl;
        return 1;
    }
    iOutSteps = _iStartStep;
    iStep = _iStartStep;
    cout << "Simulation general parameters have been successfully initialized.\n" << endl;

    /* nodes position information */
    std::cout << "Read node position information from " << achNodeInFile << "." << endl;
    finput.open(achNodeInFile.c_str(),ios::in | ios::binary );
    iNplayer = 0;
    if ( finput.is_open() ) {
        double dtemp;
        unsigned short nLayer1D;
        
        finput.read( (char*) &dtemp, sizeof(double) );
        nNodes = (unsigned long) dtemp; /* total number of nodes */
        
        finput.read( (char*) &dtemp, sizeof(double) );
        nLayer1D = (unsigned short) dtemp; /* number of 1D layers */
        for ( ii=0; ii<nLayer1D; ii++ ) { /* read 1D layer area */
            finput.read( (char*) &dtemp, sizeof(double) );
            bPara.dArea_1Dlayer.push_back(dtemp);
        }
        
        for ( ii=0; ii<nNodes; ii++ ) { /* read node type and position */
            fdmnode fn;
            fn.SetID(ii);
            finput.read( (char*) &dtemp, sizeof(double) );
            fn.type() = (unsigned int ) dtemp;
            finput.read( (char*) &fn.x(), sizeof(double) );
            finput.read( (char*) &fn.y(), sizeof(double) );
            finput.read( (char*) &fn.z(), sizeof(double) );
            NodeVec.push_back(fn);
            if ( NodeVec[ii].type() == 0 )
                iNplayer++;
            if ( !finput ) {
                cout << "Error occur while reading "<< achStateInFile << "." << endl;
                return 1;
            }
        }
        finput.close();
        finput.clear();
        std::cout << nNodes << " nodes in total (with " << iNplayer << " surface nodes).\n" << endl;
    }
    else {
        cout << achNodeInFile << ": Unable to open file." << endl;
        return 1;
    }
    
    /* Nodes state information */
    if ( _bReadLatest ) { /* read latest file from the result directory */
        DIR *dir;
        struct dirent *entry;
        int filelen;
        unsigned int n_file;
        string dir_str = achDataSubPath + "/result/";
        filelen = _nDigits + 5;
        _iStartStep = 0;
        dir = opendir(dir_str.c_str());
        if ( dir != NULL ) {
            while ( ( entry = readdir(dir)) != NULL ) {
                if ( strlen(entry->d_name) == filelen ) {
                    stringstream sstream;
                    sstream << entry->d_name;
                    sstream >> sWord;
                    n_file = std::stoi(sWord.substr( 5, filelen));
                    if ( n_file > _iStartStep ) {
                        _iStartStep = n_file;
                        achStateInFile = dir_str + entry->d_name;
                    }
                }
            }
            closedir(dir);
            iOutSteps = _iStartStep;
            iStep = _iStartStep;
        }
        else {
            cout << "Error: No result directory." << endl;
            return 1;
        }
    }
    std::cout << "Read node state information from " << achStateInFile << "." << endl;
    finput.open(achStateInFile.c_str(),ios::in | ios::binary);
    if ( finput.is_open() ) {
        unsigned int i_vol;
        finput.read( (char*) &dTimeCurrent, sizeof(double) );
        finput.read( (char*) &df_input, sizeof(double) );
        for ( ii=0; ii<nNodes; ii++ ) {
            finput.read( (char*) &NodeVec[ii].Temp(), sizeof(double) );
            for ( i_vol=0; i_vol<NVOLATILE; i_vol++ ) {
                finput.read( (char*) &NodeVec[ii].Rho_gas(i_vol), sizeof(double) );
                finput.read( (char*) &NodeVec[ii].f_ice(i_vol), sizeof(double) );
                if ( !finput ) {
                    cout << "Error occur while reading "<< achStateInFile << "." << endl;
                    return 1;
                }
            }
        }
        finput.close();
        if ( abs( dTimeCurrent - _iStartStep * _dDelta ) > 1.0e-6 )
            cout << "WARNING: Starting time ("<< dTimeCurrent << " s) of "<< achStateInFile << " is different from the one inferred from the initial time step ("<< _iStartStep * _dDelta << " s).\n" << endl;
    }
    else {
        cout << achStateInFile << ": Unable to open file." << endl;
        return 1;
    }

    
#ifdef USE_SELF
    std::cout << "Read self-shadowing information from " << achSfactorInFile << "." << endl;
    finput.open(achSfactorInFile.c_str(),ios::in | ios::binary);
    if ( finput.is_open() ) {
        double dtemp;
        finput.read( (char*) &dtemp, sizeof(double) );
        bPara.iS_number = (unsigned int) dtemp;
        for ( ii=0; ii<bPara.iS_number; ii++ ) {
            S_factor Sf;
            finput.read( (char*) &dtemp, sizeof(double) );
            Sf.vSunVector.Setx(dtemp);
            finput.read( (char*) &dtemp, sizeof(double) );
            Sf.vSunVector.Sety(dtemp);
            finput.read( (char*) &dtemp, sizeof(double) );
            Sf.vSunVector.Setz(dtemp);
            for ( jj=0; jj<iNplayer; jj++ ) {
                finput.read( (char*) &dtemp, sizeof(double) );
                Sf.S_partial_fraction.push_back(dtemp);
                if ( !finput ) {
                    cout << "Error occur while reading "<< achSfactorInFile << "." << endl;
                    return 1;
                }
            }
            bPara.S_factor_all.push_back(Sf);
        }
        finput.close();
        finput.clear();
    }
    else {
        cout << achSfactorInFile << ": Unable to open file." << endl;
        return 1;
    }
    
    std::cout << "Read self-heating information from " << achFfactorInFile << "." << endl;
    finput.open(achFfactorInFile.c_str(),ios::in | ios::binary);
    if ( finput.is_open() ) {
        double dtemp;
        unsigned int iFsize;
        for ( ii=0; ii<iNplayer; ii++ ) {
            finput.read( (char*) &dtemp, sizeof(double) );
            iFsize = (unsigned int) dtemp;
            for ( jj=0; jj<iFsize; jj++ ) {
                unsigned int iViewFacet;
                double dViewFactor;
                finput.read( (char*) &dtemp, sizeof(double) );
                iViewFacet = (unsigned int) dtemp;
                finput.read( (char*) &dtemp, sizeof(double) );
                dViewFactor = dtemp;
                NodeVec[ii].iF_view_index.push_back(iViewFacet);
                NodeVec[ii].dF_view_fraction.push_back(dViewFactor);
                if ( !finput ) {
                    cout << "Error occur while reading "<< achFfactorInFile << "." << endl;
                    return 1;
                }
            }
        }
        finput.close();
        finput.clear();
    }
    else {
        cout << achFfactorInFile << ": Unable to open file." << endl;
        return 1;
    }
#endif

    /************************/
    /**   Initialization   **/
    /************************/
    double dr1, dr2, f_ice0[NVOLATILE];
    unsigned int nChunk = 100;//ceil(nNodes/iNthreads)+1;
    dr1 = NodeVec[0].Dis3d(NodeVec[iNplayer]);
    dr2 = NodeVec[iNplayer].Dis3d(NodeVec[iNplayer+iNplayer]);
    for ( ii = 0; ii < NVOLATILE; ii++  )
        f_ice0[ii] = NodeVec[0].f_ice(ii);
    
    if ( ParaInit( df_input, dr1, dr2, f_ice0 ) )
        return 1;
    
#pragma omp parallel for shared(flag_error), schedule(static,nChunk)
    for ( ii = 0; ii < nNodes; ii++ ) {
        if ( NodeVec[ii].InitNode( NodeVec, bPara, _dDelta, iNplayer ) != 0 ) {
            cout << "Error occurred while initializing node No." << NodeVec[ii].nID() << endl;
            flag_error = true;
        }
    }
    if (flag_error)
        return 1;

    return 0;
}


bool Master::ParaInit( double df_input, double dr1, double dr2, const double (&df_ice0)[NVOLATILE] ) {
    vector3d v1,v2;
    double dRatio = bPara.dL_pore/bPara.dR_pore;
    double dT_assume, dKappa, dFactor, dDelta_recom, dTorbit;
    double dPorosity, dMass_all;
    
    dT_assume = 180.0;
    
    bPara.dP_rot *= 3600.0; /* hr to s*/
    dTorbit = 2.0 * M_PI * sqrt(pow(bPara.oEle.da,3.0)/MU_sun);
    bPara.dOmega_rot = 2.0*M_PI/bPara.dP_rot;
    bPara.dObliquity *= M_PI/180.0;
    bPara.dSolstice *= M_PI/180.0;
    
    bPara.mRotate[0].SetVec(-cos(bPara.dObliquity)*cos(bPara.dSolstice), -cos(bPara.dObliquity)*sin(bPara.dSolstice), sin(bPara.dObliquity));
    bPara.mRotate[1].SetVec(sin(bPara.dSolstice), -cos(bPara.dSolstice), 0.0);
    bPara.mRotate[2].SetVec(cos(bPara.dSolstice)*sin(bPara.dObliquity), sin(bPara.dObliquity)*sin(bPara.dSolstice), cos(bPara.dObliquity));
    
    cout << "Initial ";
    dPorosity = 1.0 - bPara.df_dust;
    bPara.dMass_dust = bPara.dRho_dust * bPara.df_dust;
    dMass_all = bPara.dMass_dust;
    for (unsigned short k = 0; k < NVOLATILE; k++ ) {
        cout << "df_ice[" << k << "]: "<< df_ice0[k] << ", ";
        dPorosity -= df_ice0[k];
        dMass_all += df_ice0[k]*bPara.dRho_ice[k];
        /* heat capacity of gas */
        bPara.dC_gas[k] = 3.0*k_B/m_gas[k];
        /* coefficient */
        bPara.dCoeff_kbmsq[k] = sqrt(k_B/(2.0*M_PI*m_gas[k]));
        /* coefficient */
        bPara.dCoeff_mkb_sat[k] = 3.56e12 * m_gas[k] / k_B;
        /* porosity is taken into account at each node */
        bPara.dD_kn_eff[k] = - bPara.dL_pore * ( 20.0 + 8.0*dRatio ) / ( 20.0 + 19.0*dRatio + 3.0*dRatio*dRatio ) / pow(bPara.dTortuosity,2.0) * sqrt(k_B/(2.0*M_PI*m_gas[k]));
    }
    cout << "df_dust:" << bPara.df_dust << endl;
    cout << "dD_kn_eff: " <<  bPara.dD_kn_eff[0] * dPorosity << endl;
    cout << "dHertz: " <<  bPara.CalHertz(dPorosity) << endl;
    dKappa = bPara.CalConductDust(dT_assume)*bPara.dMass_dust/dMass_all;
    dFactor = bPara.dRho_dust*bPara.CalCapaDust(dT_assume, sqrt(dT_assume))*bPara.df_dust;
    for (unsigned short k = 0; k < NVOLATILE; k++ ) {
        dKappa += bPara.CalConductIce(dT_assume,k)*df_ice0[k]*bPara.dRho_ice[k]/dMass_all;
        dFactor += bPara.dRho_ice[k]*bPara.CalCapaIce(dT_assume,k)*df_ice0[k];
    }
    dKappa *= bPara.CalHertz(dPorosity);
    dKappa += 4.0 * sigma * bPara.dR_pore * bPara.dEpsilon * pow(dT_assume, 3.0);
    cout << "Conductivity at T = "<< dT_assume <<" K: " << dKappa << " W/(m.K)" << endl;
    cout << "Dust heat capacity at T = "<< dT_assume <<" K: " << bPara.CalCapaDust(dT_assume,sqrt(dT_assume)) << " J/(kg.K)" << endl;
    cout << "Thermal inertia at T = "<< dT_assume <<" K: " << sqrt(dKappa*dFactor) << " W/(m.K)" << endl;
    cout << "Diurnal thermal skin depth (assuming T = "<< dT_assume <<" K): " << sqrt(bPara.dP_rot*dKappa/(2.0*M_PI*dFactor)) << " m" << endl;
    cout << "Annual thermal skin depth (assuming T = "<< dT_assume <<" K): " << sqrt(2.0*M_PI*pow(bPara.oEle.da,1.5)/sqrt(MU_sun)*dKappa/(2.0*M_PI*dFactor))  << " m" << endl;
    
    bPara.oEle.f2M();
    bPara.oEle.coe2rv( dTimeCurrent, bPara );
    cout << "df at initial time ("<< dTimeCurrent/SID_YR << " yr): " <<  bPara.oEle.f()*180.0/M_PI << " deg" << endl;
    if ( abs( df_input - bPara.oEle.f() ) > 1e-6 ) {
        cout << "Error: df ("<< df_input*180.0/M_PI << " deg) of "<< achStateInFile << " is different from the calculated one" << bPara.oEle.f() << endl;
        return 1;
    }
    
    /* check time step */
    dDelta_recom = 0.5 * dr1*dr2 * dFactor / dKappa;
    cout << "Possible maximum dDelta: " << dDelta_recom << " s." <<endl;
    if ( _dDelta > dDelta_recom ) {
        cout << "ERROR: dDelta (" << _dDelta << " s dose not ensure stability (need <" << dDelta_recom << " s)" <<endl;
        return 1;
    }
    
    return 0;
}



/* Run simulation */
bool Master::Run(void) {

    unsigned long i;
    unsigned int j;
    unsigned nChunk;
    bool flag_error = false;
    bool flag_eq; /* check if gas state equilibrium is estabilished */

    iOutSteps += _nOutInterval;

    while ( iStep < iOutSteps ) {
        
        iStep++; /* advance the step counter */
        dTimeCurrent += _dDelta;
        bPara.oEle.coe2rv( dTimeCurrent, bPara );
        nChunk = 100;//ceil(nNodes/iNthreads)+1;
        
#ifdef USE_GAS
#pragma omp parallel for schedule(static,nChunk)
        for ( i=0; i<nNodes; i++ ) {
            NodeVec[i].InitGas();
        }
        /* advance gas density and ice fraction */
        for ( j=0; j<_nDelta_all; j++ ) {
            flag_eq = true;
#pragma omp parallel for shared(flag_error,flag_eq), schedule(static,nChunk)
            for ( i=0; i<nNodes; i++ ) {
                if ( NodeVec[i].type() > 0 ) { /* interior node */
                    if ( NodeVec[i].AdvanceInterGas( NodeVec, bPara, _dDelta_gas, _dError_eq, j, iNplayer, flag_eq ) ) {
                        cout << "Error occurred while integrate the gas dynamics of internal node No." << NodeVec[i].nID() << endl;
                        flag_error = true;
                    }
                }
                else { /* boundary node */
                    if ( NodeVec[i].AdvanceBoundaryGas( NodeVec, bPara, _dDelta_gas, _dError_eq, j, iNplayer, flag_eq ) ) {
                        cout << "Error occurred while integrate the gas dynamics of boundary node No." << NodeVec[i].nID() << endl;
                        flag_error = true;
                    }
                }
            }
#pragma omp parallel for schedule(static,nChunk)
            for ( i=0; i<nNodes; i++ ) {
                NodeVec[i].UpdateGas(bPara);
            }
            
            /* check if equilibrium is estabilished */
            if ( flag_eq ) { /* advance the thermal physical state to next T integration point and jump out of the small step integration */
#pragma omp parallel for schedule(static,nChunk)
                for ( i=0; i<nNodes; i++ ) {
                    NodeVec[i].AdvanceEquil( bPara, _dDelta_gas, j, _nDelta_all );
                }
                break;
            }
            
        }
#endif
        /* advance temperature */
#pragma omp parallel for shared(flag_error), schedule(static,nChunk)
        for ( i=0; i<nNodes; i++ ) {
            if ( NodeVec[i].type() > 0 ) { /* interior node */
                if ( NodeVec[i].AdvanceInterTemp( NodeVec, bPara, _dDelta, iNplayer, _nDelta_all ) ) {
                    cout << "Error occurred while integrate the temperature of internal node No." << NodeVec[i].nID() << endl;
                    flag_error = true;
                }
            }
            else { /* boundary node */
                if ( NodeVec[i].AdvanceBoundaryTemp( NodeVec, bPara, dTimeCurrent, iNplayer, _nDelta_all ) ) {
                    cout << "Error occurred while integrate the temperature of boundary node No." << NodeVec[i].nID() << endl;
                    flag_error = true;
                }
            }
        }
        if (flag_error)
            return 1;
#pragma omp parallel for schedule(static,nChunk)
        for ( i=0; i<nNodes; i++ ) {
            NodeVec[i].UpdateTemp( NodeVec, bPara, iNplayer );
        }
        
    }
    
    return 0;
}



/* Output Data for Analyses and Visualization */
void Master::Output() {

    ofstream foutput;
    string filename;
    stringstream ss;
    double ftemp;
    unsigned long i;
    unsigned int i_vol;

    ss << std::setw(_nDigits) << std::setfill('0') << iStep;
    filename = achOutName + "." + ss.str();
    filename.insert(0,"result/");
    foutput.open( filename.c_str(), ios::out | ios::binary );
    ftemp = (double) dTimeCurrent;
    foutput.write( (char*) &ftemp, sizeof(double) );
    ftemp = (double) bPara.oEle.df;
    foutput.write( (char*) &ftemp, sizeof(double) );
    for ( i = 0; i<nNodes; i++ ) {
        ftemp = (double) NodeVec[i].Temp();
        foutput.write( (char*) &ftemp, sizeof(double) );
        for ( i_vol=0; i_vol<NVOLATILE; i_vol++ ) {
            ftemp = (double) NodeVec[i].Rho_gas(i_vol);
            foutput.write( (char*) &ftemp, sizeof(double) );
            ftemp = (double) NodeVec[i].f_ice(i_vol);
            foutput.write( (char*) &ftemp, sizeof(double) );
        }
    }
    foutput.close();
    
#ifdef USE_GAS
    if ( iStep > _iStartStep ) {
        filename.append(".gas");
        foutput.open( filename.c_str(), ios::out | ios::binary );
        for ( i = 0; i<nNodes; i++ ) {
            for ( i_vol=0; i_vol<NVOLATILE; i_vol++ ) {
                ftemp = (double) NodeVec[i].q_gas(i_vol);
                foutput.write( (char*) &ftemp, sizeof(double) );
                if ( NodeVec[i].type() < 2 ) {
                    ftemp = (double) NodeVec[i].Flux_gas1D(i_vol);
                    foutput.write( (char*) &ftemp, sizeof(double) );
                }
                else {
                    ftemp = (double) NodeVec[i].Flux_gas3Dx(i_vol);
                    foutput.write( (char*) &ftemp, sizeof(double) );
                    ftemp = (double) NodeVec[i].Flux_gas3Dy(i_vol);
                    foutput.write( (char*) &ftemp, sizeof(double) );
                    ftemp = (double) NodeVec[i].Flux_gas3Dz(i_vol);
                    foutput.write( (char*) &ftemp, sizeof(double) );
                }
            }
        }
        foutput.close();
    }
#endif
}

#endif
