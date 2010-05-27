#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCDataSet.h>
#include <BAT/BCDataPoint.h>

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


	 //create background measurement data point (T,R_B)=(100, 100) 
	BCDataPoint * backgroundMeasurement = new BCDataPoint(2);
	backgroundMeasurement->SetValue(0, 100); 
	backgroundMeasurement->SetValue(1, 100); 

	// add the single single measurement to the data set 
	BCDataSet * dataSet = new BCDataSet(); 
	dataSet->AddDataPoint(backgroundMeasurement); // register the data set with the model
	m->SetDataSet(dataSet);


	BCLog::OutSummary("Test model created");


	// perform your analysis here

	// normalize the posterior, i.e. integrate posterior
	// over the full parameter space
	m -> Normalize();

	// run MCMC and marginalize posterior wrt. all parameters
	// and all combinations of two parameters
	//	m -> MarginalizeAll();

	// run mode finding; by default using Minuit
//	m -> FindMode();

	// if MCMC was run before (MarginalizeAll()) it is
	// possible to use the mode found by MCMC as
	// starting point of Minuit minimization
	m -> FindMode( m -> GetBestFitParameters() );

	// draw all marginalized distributions into a PostScript file
	m -> PrintAllMarginalized("combination.ps");

	// calculate p-value
//	m -> CalculatePValue( m -> GetBestFitParameters() );

	// print results of the analysis into a text file
	m -> PrintResults("combination_results.txt");

	// close log file
	BCLog::CloseLog();

	delete backgroundMeasurement;
	delete dataSet;
	delete m;

	BCLog::OutSummary("Test program ran successfully");
	BCLog::OutSummary("Exiting");

	return 0;
}
