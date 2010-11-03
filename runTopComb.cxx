// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project TopComb
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCDataSet.h>
#include <iostream>
#include <fstream>

#include "TopComb.h"
#include <TRandom3.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

int main()
{

	// set nice style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
//	BCLog::OpenLog("log.txt");
	BCLog::SetLogLevel(BCLog::detail);

	// create new TopComb object
	// choose input cards
	TopComb * m = new TopComb("ljets_dilep.cards");
	//	TopComb * m = new TopComb("ljets.cards");
	//	TopComb * m = new TopComb("dilep.cards");
         
	// run settings
        m->MCMCSetNIterationsBurnIn(10000);
        m->MCMCSetNIterationsRun(100000);

        map<string, vector<double> >::iterator it;
        InputData* mydata = m->GetCards();
        vector< map<string, vector<double> > >* m_obs = mydata->GetData();
        vector<double> nobs, nobs_pois;
	BCDataSet dataset;
	BCDataPoint* DNobs = new BCDataPoint(5);
        for (int iv=0; iv<mydata->GetNDataSet(); iv++){
	  for(it=m_obs->at(iv).begin(); it != m_obs->at(iv).end(); it++){
	    cout << " Nobs  " << it->second.at(0) << endl;
            nobs.push_back(it->second.at(0));
	    DNobs->SetValue(iv,it->second.at(0));
	  }
	}
	dataset.AddDataPoint(DNobs);
	m->SetDataSet(&dataset);
        m->SetDataBoundaries(0,0.,100.);
        m->SetDataBoundaries(1,0.,100.);
        TRandom3 rnd;
        for (int iter=0; iter<1; iter++){
	  cout << " Iteration number ==============>   " << iter << endl;
  	  BCLog::OutSummary("Test model created");
          for(int i=0; i<nobs.size(); i++){
            nobs_pois.push_back((double)rnd.Poisson(nobs.at(i)));
	  }
	  //	mydata->SetNobs(nobs_pois);
	  //        m->SetNobs();

        
          for (int iv=0; iv<mydata->GetNDataSet(); iv++){
	    for(it=mydata->GetData()->at(iv).begin(); it != mydata->GetData()->at(iv).end(); it++){
	      cout << " Nobs  Poiss " << it->second.at(0) << endl;
	    }
	  }

	  // perform your analysis here

	  // normalize the posterior, i.e. integrate posterior
	  // over the full parameter space
          //	m -> Normalize();

   	  // run MCMC and marginalize posterior wrt. all parameters
	  // and all combinations of two parameters
          m->MCMCGetTRandom3()->SetSeed(12345);
	  //          m->Normalize();
	  m -> MarginalizeAll();

          BCParameter * a = m->GetParameter(0);
          double mean = m->GetMarginalized(a)->GetMean();
          double median = m->GetMarginalized(a)->GetMedian();
	  //          double max = m->GetMarginalized(a)->GetMaximum();
          double limm =  m->GetMarginalized(a)->GetMedian() - m->GetMarginalized(a)->GetQuantile(0.16);
          double limp =  m->GetMarginalized(a)->GetQuantile(0.84) - m->GetMarginalized(a)->GetMedian();

          int niter = m->MCMCGetNIterationsConvergenceGlobal();

	  nobs_pois.clear();
            
	  //	  std::cout << " maximum ---> " << max  << std::endl;
	  std::cout << " mean   ----> " << mean << std::endl;
	  std::cout << " median ----> " << median << " +" << limp << " -" << limm << std::endl;
          std::cout << " 64% C.L. ----> " << m->GetMarginalized(a)->GetQuantile(0.16) << " " << m->GetMarginalized(a)->GetQuantile(0.84) << std::endl;

	  // run mode finding; by default using Minuit
//	  m -> FindMode();

  	  // if MCMC was run before (MarginalizeAll()) it is
	  // possible to use the mode found by MCMC as
	  // starting point of Minuit minimization
//	  m -> FindMode( m -> GetBestFitParameters() );

          //	  m->FindMode();
          //          m->LikeXsec();
	  // draw all marginalized distributions into a PostScript file
	  m -> PrintAllMarginalized("TopComb_plots.ps");

  	  // save histograms in root file
          //  TH2D* h01;
	  TH1D* h0=m->GetMarginalized(m->GetParameter(1))->GetHistogram();
	  TH1D* h1=m->GetMarginalized(m->GetParameter(2))->GetHistogram();
          //TH2D* h01=m->GetMarginalized(m->GetParameter(0),m->GetParameter(1))->GetHistogram();
	  h0->SetName("h0");
	  h0->SetTitle("MMSTATE");
	  h1->SetName("h1");
	  h1->SetTitle("MMSTATM");
          //h01->SetName("h01");
          //h01->SetTitle("xsec vs Jes");
	  // calculate p-value
//	  m -> CalculatePValue( m -> GetBestFitParameters() );

	  // print results of the analysis into a text file
//	  m -> PrintResults("TopComb_results.txt");

	  // close log file
//	  BCLog::CloseLog();

	}
	m->PrintSummary();
	delete m;
	BCLog::OutSummary("Test program ran successfully");
	BCLog::OutSummary("Exiting");
	return 0;
}

