// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project ComBAT
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include <iostream>
#include <fstream>

#include <BAT/BCMath.h>
#include <BAT/BCDataSet.h>
#include <BAT/BCDataPoint.h>

#include "ComBAT.h"

#include "commons.h"

////////////////////////////////////////////

string TypeToString( const Systematic& s )
{
  string type = "";
  switch( s.type ) {
  case Unspecified:
    type = "Unspecified";
    break;
  case Symmetric:
    type = "Symmetric";
    break;
  case AsymmetricParabolic:
    type = "Asymmetric with parabolic interpolation";
    break;
  case Asymmetric2HalfGaussian:
    type = "Asymmetric with 2 half-gaussian interpolation";
    break;
  case AsymmetricDiscrete:
    type = "Asymmetric with descrete values";
    break;
  default:
    type = "Error decoding type of systematic";
    break;
  }
  return type;
}

////////////////////////////////////////////

static void FitParabola(const double y1, const double y2, const double y3, double * a, double * b, double * c)
{
  const double x1 = -1.;
  const double x2 = 0.;
  const double x3 = 1.;

  const double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
  *a  = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
  *b  = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
  *c  = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;
}


////////////////////////////////////////////


std::ostream& operator<<( std::ostream& os, const Systematic& s )
{
  // save to os
  os << s.name << " " << s.lowerBound << " " << s.upperBound << endl;
  return os;
}

std::istream& operator>>( std::istream& is, Systematic& s )
{
  // load from is
  is >> s.name >>  s.lowerBound >> s.upperBound;// >> type;
  return is;
}


double ComBAT::LogGamma( double x, double k = 1, double theta = 1)
{
  // gamma(x,k,theta) = pow(x, k-1) * exp( -1.* x / theta ) / ( pow(theta, k) * );
  return (k-1)*x - x/theta - k*log(theta) - BCMath::ApproxLogFact(k);
}

// ---------------------------------------------------------
ComBAT::ComBAT() : BCModel()
{  // default constructor
  //DefineParameters();
};

// ---------------------------------------------------------
ComBAT::ComBAT(const char * name) : BCModel(name)
{  // constructor
  const string iniFileName = string(name) + ".ini";
  //DefineParameters();
};

// ---------------------------------------------------------
ComBAT::~ComBAT()
{
  printf("Goodbye");
}  // default destructor



// ---------------------------------------------------------

void ComBAT::AddSystematic( const Systematic& param )
{
  m_systematics.push_back( param );

  AddParameter( param.name.c_str(), param.lowerBound, param.upperBound ); 
  m_paramIndex[ param.name ] = m_N_systematics + m_N_params;

  cout << "Added systematic " << m_paramIndex[ param.name ] << " " << param.name << " between " << param.lowerBound << " / " << param.upperBound <<
    " of type " << TypeToString(param) << endl;

  ++m_N_systematics; 
}


void ComBAT::AddParam( const Systematic& param  )    
{
  AddParameter( param.name.c_str(), param.lowerBound, param.upperBound ); 
  m_paramIndex[ param.name ] = m_N_systematics + m_N_params;
  
  cout << "Added parameter " << m_paramIndex[ param.name ] << " " << param.name << " between " << param.lowerBound << " / " << param.upperBound << endl;
  
  ++m_N_params; 
}

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

  BCDataSet* dataSet = new BCDataSet();

  cout << "Opening cardfile " << m_cardFileName << endl;
  ifstream cardfile( m_cardFileName.c_str(), ifstream::in  );
  if( cardfile.is_open() ) {
    string token;
    while( cardfile >> token ) {
      //cout << "token " << token << endl;
      if( token == "param" ) {
	Systematic param;
	cardfile >> param;
	AddParam( param );
      }
      else if( token == "syst" ) {
	Systematic syst;
	cardfile >> syst;

	int type = 0;
	cardfile >> type;
	syst.type = (TypeOfSystematic)type;
	
	AddSystematic( syst );
      }
      else if( token == "delta" ) {
	string whichParam = "unsetParam";
	cardfile >> whichParam;

	cout << "Sigma for param " << whichParam << ": ";

	vector< double > delta;
	for( unsigned int s = 0 ; s < 2*m_N_systematics ; ++s ) {
	  double syst = 0.0;
	  cardfile >> syst;

	  delta.push_back( syst );
	  cout << " " << syst;
	}
	cout << endl;
	m_sigma[whichParam] = delta;
      }
      else if( token == "data" ) {
	BCDataPoint * dp = new BCDataPoint( 6 );
	
	double Nobs, NsigMC, NbkgMC, eff_sig, eff_bkg, ilumi = 0.0;
	cardfile >> Nobs >> NsigMC >> NbkgMC >> eff_sig >> eff_bkg >> ilumi;

	dp->SetValue( 0, Nobs );
	dp->SetValue( 1, NsigMC );
	dp->SetValue( 2, NbkgMC );
	dp->SetValue( 3, eff_sig );
	dp->SetValue( 4, eff_bkg );
	dp->SetValue( 5, ilumi );
	dataSet->AddDataPoint( dp );
	cout << "Data point " << Nobs << " " << NsigMC << " " << NbkgMC << " " << eff_sig << " " << eff_bkg << " " << ilumi << endl;
      }
      else if( token == "end" ) {
	
	break;
      }
      else {
	cout << "Impossible to parse token " << token << endl;
	continue;
      }
    } // next token

  } // file opened

  SetDataSet( dataSet );

  cout << "End of configuration" << endl;
  cardfile.close();

}

// ---------------------------------------------------------
double ComBAT::LogLikelihood(std::vector <double> parameters)
{
	// This methods returns the logarithm of the conditional probability
	// p(data|parameters). This is where you have to define your model.

	double logprob  = 0.;

	param_t br      = 0.545;//no full had



	/////////////////////////////
	// Likelihood
	
	const double xsec_0      = parameters[ 0 ]; //m_paramIndex["xsec"] ];

	const double Nobs_ele      = GetDataPoint(DATAPOINT::ELE_0)->GetValue( VAL::Nobs );
	const double Nbkg_ele_0    = GetDataPoint(DATAPOINT::ELE_0)->GetValue( VAL::Nbkg );
	const double Nbkg_ele      = Nbkg_ele_0 * ( 1 + CalculateTotalVariation("Nbkg_ele", parameters, Nbkg_ele_0) );
	//	const double Nsig_ele_MC = GetDataPoint(DATAPOINT::ELE_0)->GetValue( VAL::NsigMC );
	const double eff_sig_ele_0 = GetDataPoint(DATAPOINT::ELE_0)->GetValue( VAL::eff_sig );
	const double eff_sig_ele   = eff_sig_ele_0 * ( 1 + CalculateTotalVariation("eff_sig_ele", parameters, eff_sig_ele_0) );
	const double iLumi_ele     = GetDataPoint(DATAPOINT::ELE_0)->GetValue( VAL::iLumi );// * ( 1 + CalculateTotalVariation("iLumi", parameters) );
	const double Nsig_ele      = xsec_0 * iLumi_ele * br * eff_sig_ele;

	const double Nobs_mu      = GetDataPoint(DATAPOINT::MU_0)->GetValue( VAL::Nobs );
	const double Nbkg_mu_0    = GetDataPoint(DATAPOINT::MU_0)->GetValue( VAL::Nbkg ); 
	const double Nbkg_mu      = Nbkg_mu_0 * ( 1 + CalculateTotalVariation("Nbkg_mu", parameters, Nbkg_mu_0) );
	//	const double Nsig_mu_MC = GetDataPoint(DATAPOINT::MU_0)->GetValue( VAL::NsigMC );
	const double eff_sig_mu_0 = GetDataPoint(DATAPOINT::MU_0)->GetValue( VAL::eff_sig );
	const double eff_sig_mu   = eff_sig_mu_0 * ( 1 + CalculateTotalVariation("eff_sig_mu", parameters, eff_sig_mu_0) );
	const double iLumi_mu     = GetDataPoint(DATAPOINT::MU_0)->GetValue( VAL::iLumi );//  * ( 1 + CalculateTotalVariation("iLumi", parameters) );
	const double Nsig_mu      = xsec_0 * iLumi_mu * br * eff_sig_mu;


	// Poisson likelihood
	const double mu_ele     = Nsig_ele + Nbkg_ele;
	const double mu_mu      = Nsig_mu  + Nbkg_mu;
	//logprob += BCMath::LogPoisson( Nobs_ele, Nsig_ele + Nbkg_ele );
	//logprob += BCMath::LogPoisson( Nobs_mu,  Nsig_mu + Nbkg_mu );

	logprob += log( mu_ele ) * Nobs_ele - mu_ele - BCMath::ApproxLogFact( int(Nobs_ele) );
	logprob += log( mu_mu )  * Nobs_mu  - mu_mu  - BCMath::ApproxLogFact( int(Nobs_mu) );

	//printf("Nobs_ele=%5.0f Nsig_ele=%5.0f Nbkg_ele=%5.0f mu_ele=%5.0f\n", Nobs_ele, Nsig_ele, Nbkg_ele, mu_ele );
	//printf("Nobs_mu=%5.0f  Nsig_mu=%5.0f  Nbkg_mu=%5.0f  mu_mu=%5.0f logprob=%5.2f\n\n", Nobs_mu, Nsig_mu, Nbkg_mu, mu_mu, logprob );

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

	for( unsigned int s = 1 ; s <= parameters.size() ; ++s ) {
	  const TypeOfSystematic type = m_systematics[s-1].type;

	  if( type == Symmetric || type == AsymmetricParabolic || type == Asymmetric2HalfGaussian ) {
	    logprob += BCMath::LogGaus( parameters[s] );
	  }
	  else {
	    logprob -= log( GetParameter(s-1)->GetRangeWidth() );
	  }
	  //logprob += LogGamma( parameters[s] );
	}


	return logprob;
}
// ---------------------------------------------------------


double ComBAT::CalculateTotalVariation( const string& param, std::vector <double>& parameters, const double& val_0 )
{
  double v = 0.;// parameters[0];  //parameters[0] = sigma
  //cout << param << ": ";
  m_N_systematics = m_sigma.size();
  for( unsigned int p = 0 ; p < 1 + m_N_systematics; ++p ) {
    const unsigned int s = 2 * p;

    double delta = 0.0;
    double a, b, c;
    double val_m, val_p;

    const double sigma_p = m_sigma[param][s];
    const double sigma_m = m_sigma[param][s+1];
    //cout << "s_p=" << sigma_p << "  s_m=" << sigma_m << endl;
    switch( m_systematics[p].type ) {

    case Unspecified:
      delta = 0.0;
      break;
    case Symmetric:
      delta = sigma_p * parameters[ p + 1 ];
      break;
    case AsymmetricParabolic:
      val_m = (1. + sigma_m) * val_0;
      val_p = (1. + sigma_p) * val_0;
      FitParabola( val_m, val_0, val_p, &a, &b, &c);
      //cout << a << " " << b << " " << c << " " << endl;
      delta = c + b*parameters[p+1] + a*pow(parameters[p+1],2);
      break;
    case Asymmetric2HalfGaussian:
      if( parameters[ p + 1 ] > 0. ) {
	delta = sigma_p * parameters[ p + 1 ] ;
      }
      else {
	delta = -sigma_m * parameters[ p + 1 ] ;
      }
      break;
    case AsymmetricDiscrete:
      delta = sigma_p * parameters[ p + 1 ];
      break;
    default:
      delta = sigma_p * parameters[ p + 1 ];
      break;
    }

    v += delta;
  }
  //cout << endl;
  return v;
}
