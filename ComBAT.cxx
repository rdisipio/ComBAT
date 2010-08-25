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

string PriorToString( const Systematic& s )
{
  string type = "";
  switch( s.type ) {
  case FlatPrior:
    type = "Flat";
    break;
  case GaussPrior:
    type = "Gaussian";
    break;
  case JeffreyPrior:
    type = "Jeffrey's";
    break;
  default:
    type = "Unknown";
    break;
  }
 
  return type;
};

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
  m_h_Nbkg->Write();
  m_h_Nbkg_conv_poisson->Write();
  m_histoFile->Close();
  printf("User-defined histograms saved");

  //delete m_h_Nbkg;
  //delete m_h_Nbkg_conv_poisson;
  //delete m_histoFile;

  printf("Goodbye");
}  // default destructor



// ---------------------------------------------------------

void ComBAT::AddSystematic( const Systematic& param )
{
  m_systematics.push_back( param );

  AddParameter( param.name.c_str(), param.lowerBound, param.upperBound ); 
  m_paramIndex[ param.name ] = m_N_systematics + m_N_params;

  cout << "Added systematic " << m_paramIndex[ param.name ] << " " << param.name << " range " << param.lowerBound << " / " << param.upperBound <<
    " of type " << TypeToString(param) << endl;

  ++m_N_systematics; 
}


void ComBAT::AddParam( const Systematic& param  )    
{
  AddParameter( param.name.c_str(), param.lowerBound, param.upperBound ); 
  m_paramIndex[ param.name ] = m_N_systematics + m_N_params;
  
  cout << "Added parameter " << m_paramIndex[ param.name ] << " " << param.name << " range " << param.lowerBound << " / " << param.upperBound
       << " with prior " << PriorToString(param) << endl;
  
  ++m_N_params; 
  m_xsec = param;
}


void ComBAT::AddSigma( const string& delta, const Systematic& param  )    
{
  m_sigma[delta].push_back( param.upperBound );
  m_sigma[delta].push_back( param.lowerBound );

  const unsigned int size = m_sigma[delta].size();
  //cout << size << endl;
  cout << "   sigma_" << param.name << "_up= " << m_sigma[delta][size-2] << " / sigma_" << param.name << "_down = " << m_sigma[delta][size-1] << endl;
}


// ---------------------------------------------------------



void ComBAT::DefineParameters()
{
  /*
	// Add parameters to your model here.
	// You can then use them in the methods below by calling the
	// parameters.at(i) or parameters[i], where i is the index
	// of the parameter. The indices increase from 0 according to the
	// order of adding the parameters.

//	this -> AddParameter("par1", 0.0, 1.0);   // index 0
//	this -> AddParameter("par2", -7.0, 28.0); // index 1


*/
  cout << "Before starting calculation:" << endl;
  cout << "No. of systematics: " << m_systematics.size() << endl;// m_N_systematics << endl;
  cout << "No. of data points: " << GetNDataPoints() << endl;

  
}


// ---------------------------------------------------------


void ComBAT::SetUserDefinedHistogramFile( string filename )
{
  m_histoFile = new TFile( filename.c_str(), "recreate" );

  printf("User-defined histograms will be saved in file %s\n", filename.c_str() );

  m_h_Nbkg = new TH1D( "Nbkg", "Number of background events", 31, -0.5, 30 );
  m_h_Nbkg_conv_poisson = new TH1D( "Nbkg_conv_poisson", "Number of background events convoluted w/ Poisson", 31, -0.5, 30 );
}


// ---------------------------------------------------------


double ComBAT::LogLikelihood(std::vector <double> parameters)
{
	// This methods returns the logarithm of the conditional probability
	// p(data|parameters). This is where you have to define your model.

	double logprob  = 0.;

	/////////////////////////////
	// Likelihood

	// Poisson likelihood
	logprob += MarginalLikelihood( DATAPOINT::ELE_0, parameters );
	logprob += MarginalLikelihood( DATAPOINT::MU_0, parameters );
	
	return logprob;
}


// ---------------------------------------------------------


double ComBAT::MarginalLikelihood( nparam_t datapoint, const std::vector <double>& parameters )
{
  m_Nbkg = 0.;
  string ch   = ( datapoint == DATAPOINT::ELE_0 ) ? "ele" : "mu";  
  BCDataPoint* data = GetDataPoint(datapoint);

  const double xsec_0        = parameters[ 0 ]; //m_paramIndex["xsec"] ];
  const double br            = 0.545;

  const double Nobs          = data->GetValue( VAL::Nobs );

  const double f_b           = data->GetValue( VAL::W_b_frac )      * ( 1 + CalculateTotalVariation( "f_b", parameters) );
  const double f_c           = data->GetValue( VAL::W_c_frac )      * ( 1 + CalculateTotalVariation( "f_c", parameters) );
  const double eff_b         = data->GetValue( VAL::eff_btag )      * ( 1 + CalculateTotalVariation( "eff_b", parameters) );
  const double eff_mistag_c  = data->GetValue( VAL::eff_mistag_c )  * ( 1 + CalculateTotalVariation( "eff_mistag_c", parameters) );
  const double eff_mistag_lq = data->GetValue( VAL::eff_mistag_lq ) * ( 1 + CalculateTotalVariation( "eff_mistag_lq", parameters) );
  const double eff_btag      = ( f_b * + eff_b ) +  ( f_c * eff_mistag_c ) + ( (1 - f_c - f_b ) * eff_mistag_lq ); 
	                            
  const double NW_0      = data->GetValue( VAL::NbkgW_datadriven);
  const double NW_b      = NW_0 * eff_btag;
	
  const double NbkgQCD    = data->GetValue( VAL::NbkgQCD_datadriven) * ( 1 + CalculateTotalVariation( "NbkgQCD_" + ch, parameters) );
  const double NbkgST     = data->GetValue( VAL::NbkgST)             * ( 1 + CalculateTotalVariation( "NbkgST_" + ch, parameters) );
  const double NbkgZ      = data->GetValue( VAL::NbkgZ)              * ( 1 + CalculateTotalVariation( "NbkgZ_" + ch, parameters) );
  const double NbkgOthers = data->GetValue( VAL::NbkgOthers)         * ( 1 + CalculateTotalVariation( "NbkgOthers_" + ch, parameters) );
  m_Nbkg       = NW_b + NbkgQCD + NbkgST + NbkgZ + NbkgOthers;
  const double Nbkg = m_Nbkg;

  const double eff_sig   = data->GetValue( VAL::eff_sig ) * ( 1 + CalculateTotalVariation( "eff_sig_" + ch , parameters ) );
  const double iLumi     = data->GetValue( VAL::iLumi ) * ( 1 + CalculateTotalVariation("ilumi", parameters) );
  const double Nsig      = xsec_0 * iLumi * br * eff_sig;

  //logprob += log( mu_mu )  * Nobs_mu  - mu_mu  - BCMath::ApproxLogFact( int(Nobs_mu) );
  //printf("Nobs_ele=%5.0f Nsig_ele=%5.0f Nbkg_ele=%5.0f mu_ele=%5.0f\n", Nobs_ele, Nsig_ele, Nbkg_ele, mu_ele );
  return BCMath::LogPoisson( Nobs, Nsig + Nbkg );
}

// ---------------------------------------------------------


double ComBAT::LogAPrioriProbability(std::vector <double> parameters)
{
	// This method returns the logarithm of the prior probability for the
	// parameters p(parameters).

	double logprob = 0.;

	// For flat prior it's very easy.
//	for(unsigned int i=0; i < this -> GetNParameters(); i++)

	if( m_xsec.type == FlatPrior ) {
	  logprob -= log( GetParameter("xsec")->GetRangeWidth() ); // flat a priori on xs
	}
	else if( m_xsec.type == GaussPrior ) { 
	  logprob += BCMath::LogGaus( parameters[0], GetParameter("xsec")->GetRangeWidth()/2., GetParameter("xsec")->GetRangeWidth()/4. );
	}
	else if( m_xsec.type == JeffreyPrior ) {
	  logprob += 0;
	}
	
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


double ComBAT::CalculateTotalVariation( string param, const std::vector <double>& parameters ) 
{
  double v = 0.;// parameters[0];  //parameters[0] = sigma
  //cout << param << ": ";
  //m_N_systematics = m_sigma.size();
  for( unsigned int p = 0 ; p < 1 + m_systematics.size() /* m_N_systematics */; ++p ) {
    const unsigned int s = 2 * p;

    double delta = 0.0;
    double a, b, c;
    //double val_m, val_p;

    const double sigma_p = m_sigma[param][s];
    const double sigma_m = m_sigma[param][s+1];
    //cout << m_systematics[p].name << " s_p=" << sigma_p << "  s_m=" << sigma_m << endl;
    switch( m_systematics[p].type ) {

    case Unspecified:
      delta = 0.0;
      break;
    case Symmetric:
      delta = sigma_p * parameters[ p + 1 ];
      break;
    case AsymmetricParabolic:
      //val_m = (1. + sigma_m) * val_0;
      //val_p = (1. + sigma_p) * val_0;
      FitParabola( sigma_m, 0, sigma_p, &a, &b, &c);
      //cout << a << " " << b << " " << c << " " << endl;
      delta = c + b*parameters[p+1] + a*pow(parameters[p+1],2);
      //cout << "param " << param << ": a=" << a << " b=" << b << " c=" << c << endl;

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


// ---------------------------------------------------------

void ComBAT::MCMCUserIterationInterface()
{
  m_h_Nbkg->Fill( m_Nbkg );

}
