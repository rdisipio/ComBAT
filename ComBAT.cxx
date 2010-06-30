// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project ComBAT
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include <iostream>
#include <BAT/BCMath.h>

#include "ComBAT.h"

#include "commons.h"

// ---------------------------------------------------------
ComBAT::ComBAT() : BCModel()
{  // default constructor
	DefineParameters();
};

// ---------------------------------------------------------
ComBAT::ComBAT(const char * name) : BCModel(name)
{  // constructor
	DefineParameters();
};

// ---------------------------------------------------------
ComBAT::~ComBAT()
{
  printf("Goodbye");
}  // default destructor

// ---------------------------------------------------------
void ComBAT::DefineParameters()
{
	// Add parameters to your model here.
	// You can then use them in the methods below by calling the
	// parameters.at(i) or parameters[i], where i is the index
	// of the parameter. The indices increase from 0 according to the
	// order of adding the parameters.

//	this -> AddParameter("par1", 0.0, 1.0);   // index 0
//	this -> AddParameter("par2", -7.0, 28.0); // index 1

  AddParameter( "xsec", 0., 600. ); //0

  AddParameter( "eff_sig_ele", -3, 3.  ); //1
  AddParameter( "eff_sig_mu", -3, 3.  ); //2

  // systematics
  AddParameter( "JES",   -5, 5. ); //3
  AddParameter( "Wnorm", -5, 5.); //4
  AddParameter( "iLumi",  -5, 5. ); //5
  
}

// ---------------------------------------------------------
double ComBAT::LogLikelihood(std::vector <double> parameters)
{
	// This methods returns the logarithm of the conditional probability
	// p(data|parameters). This is where you have to define your model.

	double logprob  = 0.;

	param_t br      = 0.545;

	// Get parameters
	param_t xsec_0        = parameters[PARAM::xsec];
	
	param_t eff_sig_ele_0 = parameters[PARAM::eff_sig_ele];
	param_t eff_sig_mu_0  = parameters[PARAM::eff_sig_mu];

	//systematics
	syst_t d_JES        = parameters[SYST::JES];
	syst_t d_Wnorm      = parameters[SYST::Wnorm];
	syst_t d_iLumi      = parameters[SYST::iLumi];

	// define systematic effects (1sigma)
	double systMatrix[2][3];
	//syst_matrix_t systMatrix = new 
	
	systMatrix[0][0] = 0.1;
	systMatrix[0][1] = 0.3;
	systMatrix[0][2] = 0.02;
	
	systMatrix[0][0] = 0.14;
	systMatrix[1][1] = 0.245;
	systMatrix[2][2] = 0.01;


	/////////////////////////////
	// Likelihood
	
	const double Nobs_ele    = GetDataPoint(DATAPOINT::ELE_0)->GetValue( VAL::Nobs );
	const double Nbkg_ele    = GetDataPoint(DATAPOINT::ELE_0)->GetValue( VAL::Nbkg );
	const double Nsig_ele    = Nobs_ele - Nbkg_ele;
	const double Nsig_ele_MC = GetDataPoint(DATAPOINT::ELE_0)->GetValue( VAL::NsigMC );
	const double eff_sig_ele = GetDataPoint(DATAPOINT::ELE_0)->GetValue( VAL::eff_sig ) * 
	  ( 1 +
	    systMatrix[0][0] * d_JES   +
	    systMatrix[0][1] * d_Wnorm +
	    systMatrix[0][2] * d_iLumi
	    );
	const double iLumi_ele_0 = GetDataPoint(DATAPOINT::ELE_0)->GetValue( VAL::iLumi ) * ( 1. + systMatrix[2][0] * d_iLumi );
	double xsec_ele          = Nsig_ele / ( br * iLumi_ele_0 * eff_sig_ele ) ;


	const double Nobs_mu    = GetDataPoint(DATAPOINT::MU_0)->GetValue( VAL::Nobs );
	const double Nbkg_mu    = GetDataPoint(DATAPOINT::MU_0)->GetValue( VAL::Nbkg );
	const double Nsig_mu    = Nobs_mu - Nbkg_mu;
	const double Nsig_mu_MC = GetDataPoint(DATAPOINT::MU_0)->GetValue( VAL::NsigMC );
	const double eff_sig_mu = GetDataPoint(DATAPOINT::MU_0)->GetValue( VAL::eff_sig ) * 
	  ( 1 +
	    systMatrix[1][0] * d_JES   +
	    systMatrix[2][1] * d_Wnorm +
	    systMatrix[3][2] * d_iLumi
	    );
	const double iLumi_mu_0 = GetDataPoint(DATAPOINT::MU_0)->GetValue( VAL::iLumi ) * ( 1. + systMatrix[2][1] * d_iLumi );
	double xsec_mu          = Nsig_mu / ( br * iLumi_mu_0 * eff_sig_mu );
	

	//logprob += log( Nsig_ele_MC + Nbkg_ele ) * Nobs_ele - ( Nsig_ele + Nbkg_ele ) - BCMath::ApproxLogFact( int(Nobs_ele) );
	//logprob += log( Nsig_mu_MC + Nbkg_ele ) * Nobs_mu - ( Nsig_mu + Nbkg_mu ) - BCMath::ApproxLogFact( int(Nobs_mu) );

	//logprob += log( Nsig_mu )  * Nobs_mu  - Nsig_ele - BCMath::LogFact( int(Nobs_mu) );
	//logprob += BCMath::LogGaus( Nsig_ele );
	//logprob += BCMath::LogGaus( Nsig_mu );

	logprob += BCMath::LogGaus( xsec_ele, xsec_0, sqrt( Nsig_ele_MC + Nbkg_ele) ) + BCMath::LogGaus( xsec_mu, xsec_0, sqrt( Nsig_ele_MC + Nbkg_ele)  );
	//logprob += BCMath::LogPoisson( xsec_ele, xsec_0 ) + BCMath::LogPoisson( xsec_mu, xsec_0 );
	//logprob += BCMath::LogPoisson( Nobs_ele,  
	return logprob;
}

// ---------------------------------------------------------
double ComBAT::LogAPrioriProbability(std::vector <double> parameters)
{
	// This method returns the logarithm of the prior probability for the
	// parameters p(parameters).

	double logprob = 0.;

	// For flat prior it's very easy.
//	for(unsigned int i=0; i < this -> GetNParameters(); i++)
	logprob -= log( GetParameter("xsec")->GetRangeWidth() ); // flat a priori on xs
	//logprob -= log( GetParameter("lumi")->GetRangeWidth() );
	logprob += BCMath::LogGaus( parameters[SYST::JES] );
	logprob += BCMath::LogGaus( parameters[SYST::Wnorm] );
	logprob += BCMath::LogGaus( parameters[SYST::iLumi] );
	logprob += BCMath::LogGaus( parameters[PARAM::eff_sig_ele] );
	logprob += BCMath::LogGaus( parameters[PARAM::eff_sig_mu] );

	return logprob;
}
// ---------------------------------------------------------
