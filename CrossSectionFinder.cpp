// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project CrossSectionFinder
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include "CrossSectionFinder.h"

// ---------------------------------------------------------
CrossSectionFinder::CrossSectionFinder() : BCModel()
{  // default constructor
  DefineParameters();
};

// ---------------------------------------------------------
CrossSectionFinder::CrossSectionFinder(const char * name) : BCModel(name)
{  // constructor
  DefineParameters();
};

// ---------------------------------------------------------
CrossSectionFinder::~CrossSectionFinder()
{};  // default destructor

// ---------------------------------------------------------
void CrossSectionFinder::DefineParameters()
{
	// Add parameters to your model here.
	// You can then use them in the methods below by calling the
	// parameters.at(i) or parameters[i], where i is the index
	// of the parameter. The indices increase from 0 according to the
	// order of adding the parameters.

//	this -> AddParameter("par1", 0.0, 1.0);   // index 0
//	this -> AddParameter("par2", -7.0, 28.0); // index 1
  AddParameter("R_B", 0.0, 5.0);

}

// ---------------------------------------------------------
double CrossSectionFinder::LogLikelihood(std::vector <double> parameters)
{
	// This methods returns the logarithm of the conditional probability
	// p(data|parameters). This is where you have to define your model.

	double logprob = 0.;

	// get background measurement
	double T = GetDataPoint(0)->GetValue(0);
	double N1 = GetDataPoint(0)->GetValue(1);

	// extract value of background rate
	double R_B = parameters.at(0);

	// calculate expected counts given background rate
	double n_B = R_B * T;

	// update likelihood
	logprob += -n_B + N1 * log(n_B) - BCMath::LogFact(N1);

	return logprob;
}

// ---------------------------------------------------------
double CrossSectionFinder::LogAPrioriProbability(std::vector <double> parameters)
{
	// This method returns the logarithm of the prior probability for the
	// parameters p(parameters).

	double logprob = 0.;

	// For flat prior it's very easy.
	for(unsigned int i=0; i < this -> GetNParameters(); i++)
		logprob -= log(this -> GetParameter(i) -> GetRangeWidth());

	return logprob;
}
// ---------------------------------------------------------



