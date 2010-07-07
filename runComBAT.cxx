// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project ComBAT
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include <string>

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCDataSet.h>
#include <BAT/BCDataPoint.h>
#include <BAT/BCH1D.h>
#include <BAT/BCModelOutput.h>

#include "commons.h"
#include "ComBAT.h"

int main( int argc, char *argv[] )
{
  //const string dataFileName   = ( argc > 1 ) ? argv[1]  : "analysis_10TeV.txt" ;
  //const string configFileName = ( argc > 2 ) ? argv[2]  : "combat.ini" ;
  //const unsigned int NParamsToRead = 6;
  const string configFileName = ( argc > 1 ) ? argv[1]  : "combat.bat" ;

  // set nice style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
//	BCLog::OpenLog("log.txt");
	BCLog::SetLogLevel(BCLog::detail);

	// create new ComBAT object
	ComBAT * m = new ComBAT();

	m->SetCardFileName( configFileName );
	m->DefineParameters();

	//BCModelOutput * o = new BCModelOutput(m,"ComBAT_MarkovChains.root");
	//o->WriteMarkovChain(true);

	// Add the data point to the data set
	//cout << "Reading data points from file" << endl;
        //BCDataSet* dataSet = new BCDataSet(); 
	//dataSet->ReadDataFromFile( dataFileName.c_str(), NParamsToRead );

	// dataSet->AddDataPoint( ch_ele_0 );
	//dataSet->AddDataPoint( ch_mu_0 );
	 
	//m->SetDataSet( dataSet );

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

	//o->WriteMarginalizedDistributions();
	//o->Close();

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

	//delete o;
	delete m;
	//delete ch_ele;
	//delete ch_mu;
	//delete m;

	BCLog::OutSummary("Test program ran successfully");
	BCLog::OutSummary("Exiting");

	return 0;

}

