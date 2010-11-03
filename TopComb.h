// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project TopComb
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#ifndef __TOPCOMB__H
#define __TOPCOMB__H

#include <BAT/BCModel.h>
#include <BAT/BCMath.h>
#include <BAT/BCModelOutput.h>
#include <BAT/BCEngineMCMC.h>
#include <BAT/BCDataPoint.h>
#include <BAT/BCDataSet.h>
#include <BAT/BCH1D.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include "InputData.h"
#include <vector.h>

// This is a TopComb header file.
// Model source code is located in file TopComb/TopComb.cxx

// ---------------------------------------------------------
class TopComb : public BCModel
{
	public:

		// Constructors and destructor
		TopComb();
		TopComb(const char * name);
		~TopComb();

		// Methods to overload, see file TopComb.cxx
                void ReadCards(const char*);
		void DefineParameters();
		double LogAPrioriProbability(std::vector <double> parameters);
		double LogLikelihood(std::vector <double> parameters);
                InputData* GetCards();
                void SetNobs();
                TH2D *hnbgqcd_we, *hnbgqcd_wm;
                TH1D *hnsig, *hnbg, *hnbgw, *hnbgwnob, *hnbgqcd, *hnbgst, *hnbgz;
                TH1D *hnsige, *hnbge, *hnbgwe, *hnbgwnobe, *hnbgqcde, *hnbgste, *hnbgze;
                TH1D *hnsigm, *hnbgm, *hnbgwm, *hnbgwnobm, *hnbgqcdm, *hnbgstm, *hnbgzm;
                TH1D *hQCDSYSTE, *hQCDSTATE, *hWBGSTATE, *hWBGSYST, *hW2TAGSTAT, *hBTAGEFF, *hCTAGEFF, *hLTAGEFF, *hWHFRAC, *hWCFRAC, *hZXSEC, *hLUMI, *hJES;
                TH1D *hJER, *hJERECO, *hSTXS, *hEID, *hEES, *hEER, *hIFSR, *hNLOGEN, *hPSGEN, *hPDF, *hOTX, *hQCDSYSTM, *hQCDSTATM, *hWBGSTATM, *hMID;
                TH1D *hMES, *hMER, *hxsec, *heff_e, *heff_m, *like_xsec, *hpvalue;

                vector<TH1D*> Gethistos();
                vector<TH2D*> Gethistos2();
                void LikeXsec();

        private:
                int nneg;
                double Nsig_tot, Nbg_tot, Nbgw_tot, Nbgwnob_tot, Nbgqcd_tot, Nbgst_tot, Nbgz_tot;
                double Nsig_e, Nbg_e, Nbgw_e, Nbgwnob_e, Nbgqcd_e, Nbgst_e, Nbgz_e;
                double Nsig_m, Nbg_m, Nbgw_m, Nbgwnob_m, Nbgqcd_m, Nbgst_m, Nbgz_m, testg, testg2, eff_e, eff_m;
                vector<double> sys;
                vector<int> likelihood;
                vector< map<string,double> > m_parcv;
                vector< map<string,map<string,vector<double> > > > m_parsyst;
                map<string, vector<double> > m_systintlim;
                map<string, vector<double> > m_gammapar;
                map<string, vector<int> > m_systtyp;
                map<string, vector<int> > m_meastyp;
                vector< map<string, vector<double> > > m_data;
                map<string,vector<double> > m_measpar, m_measlim;
                void MCMCUserIterationInterface();
                void MCMCrun();
                vector<TH1D*> hvec;
                vector<TH2D*> hvec2;
                void parabola(vector<double>,double*, double*);  
                double Gamma(double, double, double);
		InputData* mydata;
                bool CheckforNegativeValues();
};
// ---------------------------------------------------------

#endif

