// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project ComBAT
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include <string>
#include <sstream>
#include <vector>

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCDataSet.h>
#include <BAT/BCDataPoint.h>
#include <BAT/BCH1D.h>
#include <BAT/BCModelOutput.h>

#include "commons.h"
#include "ComBAT.h"

std::vector<std::string> &split(const std::string &s, const char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


int main( int argc, char *argv[] )
{
  //const string dataFileName   = ( argc > 1 ) ? argv[1]  : "analysis_10TeV.txt" ;
  //const string configFileName = ( argc > 2 ) ? argv[2]  : "combat.ini" ;
  //const unsigned int NParamsToRead = 6;
  const string configFileName = ( argc > 1 ) ? argv[1]  : "combat.bat" ;
  std::vector<std::string> tokens;
  split( configFileName, '.', tokens);

  // set nice style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
//	BCLog::OpenLog("log.txt");
	BCLog::SetLogLevel(BCLog::detail);

	// create new ComBAT object
	ComBAT * m = new ComBAT( tokens[0].c_str() );

	m->SetCardFileName( configFileName );
	m->DefineParameters();

        string rootFileName = tokens[0] + "_MarkovChains.root";
	BCModelOutput * o = new BCModelOutput(m, rootFileName.c_str() );
	o->WriteMarkovChain(true);

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
	const string psFileName = tokens[0] + "_plots.ps";
	m -> PrintAllMarginalized( psFileName.c_str() );


	// calculate p-value
	m -> CalculatePValue( m -> GetBestFitParameters() );

	int npar = m->GetNParameters();
        for (int i=0; i<npar; i++){
          BCParameter * a = m->GetParameter(i);
	  std::cout << " 95% C.L. ----> " << m->GetMarginalized(a)->GetQuantile(0.95) << std::endl;
	}

	// print results of the analysis into a text file
	const string txtFileName = tokens[0] + "_results.txt";
	m -> PrintResults( txtFileName.c_str() );

	o->WriteMarginalizedDistributions();
	o->Close();

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

