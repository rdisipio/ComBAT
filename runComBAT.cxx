#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include "CrossSectionFinder.h"

int main()
{	
// set nice style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
//	BCLog::OpenLog("log.txt");
	BCLog::SetLogLevel(BCLog::detail);

	// create new CountingExp object
	CrossSectionFinder * m = new CrossSectionFinder();

	BCLog::OutSummary("Test model created");

	// perform your analysis here

	// normalize the posterior, i.e. integrate posterior
	// over the full parameter space
//	m -> Normalize();

	// run MCMC and marginalize posterior wrt. all parameters
	// and all combinations of two parameters
//	m -> MarginalizeAll();

	// run mode finding; by default using Minuit
//	m -> FindMode();

	// if MCMC was run before (MarginalizeAll()) it is
	// possible to use the mode found by MCMC as
	// starting point of Minuit minimization
//	m -> FindMode( m -> GetBestFitParameters() );

	// draw all marginalized distributions into a PostScript file
//	m -> PrintAllMarginalized("CountingExp_plots.ps");

	// calculate p-value
//	m -> CalculatePValue( m -> GetBestFitParameters() );

	// print results of the analysis into a text file
//	m -> PrintResults("CountingExp_results.txt");

	// close log file
//	BCLog::CloseLog();

	delete m;

	BCLog::OutSummary("Test program ran successfully");
	BCLog::OutSummary("Exiting");

	return 0;
}
