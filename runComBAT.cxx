// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project ComBAT
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCDataSet.h>
#include <BAT/BCDataPoint.h>
#include <BAT/BCH1D.h>

#include "ComBAT.h"

#include "commons.h"

int main( int argc, char *argv[] )
{
  const char* inFile = ( argc > 1 ) ? argv[1]  : "analysis_10TeV.txt" ;
  const unsigned int NParamsToRead = 6;

  // set nice style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
//	BCLog::OpenLog("log.txt");
	BCLog::SetLogLevel(BCLog::detail);

	// create new ComBAT object
	ComBAT * m = new ComBAT();

	// ele channel
	//	BCDataPoint* ch_ele_0 = new BCDataPoint(5);
	// ch_ele_0->SetValue( VAL::Nobs, 3260 );
// 	ch_ele_0->SetValue( VAL::Nbkg,  1239 );
// 	ch_ele_0->SetValue( VAL::eff_sig, 1703./21495. );
// 	ch_ele_0->SetValue( VAL::eff_bkg, 1239./5261417. );
// 	ch_ele_0->SetValue( VAL::iLumi, 99. );

	// mu channel
	//	BCDataPoint* ch_mu_0 = new BCDataPoint(5);
	// ch_mu_0->SetValue( VAL::Nobs, 3374 );
// 	ch_mu_0->SetValue( VAL::Nbkg,  1322 );
// 	ch_mu_0->SetValue( VAL::eff_sig, 1831. / 21495. );
// 	ch_mu_0->SetValue( VAL::eff_bkg, 1322. / 5261417. );
// 	ch_mu_0->SetValue( VAL::iLumi, 99. );

	// JES up	
	// BCDataPoint* ch_ele_JES_up = new BCDataPoint(5);
// 	ch_ele_0->SetValue( VAL::Nobs, 1746 );
// 	ch_ele_0->SetValue( VAL::Nbkg,  820 );
// 	ch_ele_0->SetValue( VAL::eff_sig, 0.99 );
// 	ch_ele_0->SetValue( VAL::eff_bkg, 0.055 );
// 	ch_ele_0->SetValue( VAL::iLumi, 99. );

	// Add the data point to the data set
        BCDataSet* dataSet = new BCDataSet(); 
	dataSet->ReadDataFromFile( inFile, NParamsToRead );

	// dataSet->AddDataPoint( ch_ele_0 );
	//dataSet->AddDataPoint( ch_mu_0 );
	 
	m->SetDataSet( dataSet );

	BCLog::OutSummary("Test model created");

	// perform your analysis here

	// normalize the posterior, i.e. integrate posterior
	// over the full parameter space
	//	m -> Normalize();

	// run MCMC and marginalize posterior wrt. all parameters
	// and all combinations of two parameters
	m -> MarginalizeAll();

	// run mode finding; by default using Minuit
//	m -> FindMode();

	// if MCMC was run before (MarginalizeAll()) it is
	// possible to use the mode found by MCMC as
	// starting point of Minuit minimization
	m -> FindMode( m -> GetBestFitParameters() );

	// draw all marginalized distributions into a PostScript file
		m -> PrintAllMarginalized("ComBAT_plots.ps");

	// calculate p-value
	m -> CalculatePValue( m -> GetBestFitParameters() );

	int npar = m->GetNParameters();
        for (int i=0; i<npar; i++){
          BCParameter * a = m->GetParameter(i);
	  std::cout << " 95% C.L. ----> " << m->GetMarginalized(a)->GetQuantile(0.95) << std::endl;
	}

	// print results of the analysis into a text file
		m -> PrintResults("ComBAT_results.txt");

	// close log file
	BCLog::CloseLog();

	
	//delete ch_ele;
	//delete ch_mu;
	//delete m;

	BCLog::OutSummary("Test program ran successfully");
	BCLog::OutSummary("Exiting");

	return 0;

}

