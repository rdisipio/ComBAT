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
  switch( s.prior ) {
  case FlatPrior:
    type = "Flat";
    break;
  case GaussPrior:
    type = "Gaussian";
    break;
  case JeffreysPrior:
    type = "Jeffrey's";
    break;
  case GammaPrior:
    type = "Gamma";
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
  switch( s.shape ) {
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


// ---------------------------------------------------------


// double ComBAT::LogGamma( double x, double k, double theta )
// {
//  return (k -1 ) * log(x) - ( x / theta ) - k * log(theta) - GammaLanczos(k); //BCMath::ApproxLogFact(k);
// }


//---------------------------------------------------------

double ComBAT::LogGamma( double mu, double sigma, double x )
{
  //  return -pow((x-mu)/sig,2)/2;
  const double k     = pow(mu/sigma,2);
  const double theta = pow(sigma,2)/mu;

  return log( TMath::GammaDist( x, k, 0, theta ));
}

//---------------------------------------------------------


double ComBAT::GammaLanczos( double z )
{
  static double q[7] = { 75122.6331530, 80916.6278952, 36308.2951477, 8687.24529705, 1168.92649479, 83.8676043424, 2.50662827511 };

  double f = pow( z + 5.5, z + 0.5) * exp( -1. * (z + 5.5) );

  double Num = 0;
  double Den = 1;
  for( unsigned int n = 0 ; n < 7 ; ++n ) {
    Num += q[n] * pow( z, n );
    Den *=  ( z + n );
  }
  
  return ( Num * f / Den );
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
  //delete m_h_Nbkg;
  //delete m_h_Nbkg_conv_poisson;
  //delete m_histoFile;

  printf("Goodbye");
}  // default destructor



// ---------------------------------------------------------

void ComBAT::AddSystematic( Systematic& param )
{
  m_systematics[param.name] = param;

  AddParameter( param.name.c_str(), param.lowerBound, param.upperBound ); 
  param.index =  m_N_systematics + m_N_params;
  m_paramIndex[ param.name ] = param.index;

  const double sigma = param.params["sigma"]; //param.params["sigma"];

  cout << "Added systematic " << m_paramIndex[ param.name ] << " " << param.name << " range " << param.lowerBound << " / " << param.upperBound <<
    " with prior " << PriorToString(param) << " shape " << TypeToString(param) << " sigma=" << sigma << endl;

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
  m_sigma[delta][param.name] = param;
  
  cout << "   sigma_" << param.name << "_up= " << m_sigma[delta][param.name].upperBound 
       << " / sigma_" << param.name << "_down = " << m_sigma[delta][param.name].lowerBound << endl;
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

  const int bmax  = 61;
  const int nbins = 10 * bmax;
  
  m_h_Nbkg        = new TH1D( "Nbkg",        "Number of background events (All)",    nbins, -0.5, bmax - 0.5 );
  m_h_Nbkg_W      = new TH1D( "Nbkg_W",      "Number of background events (W+jets)", nbins, -0.5, bmax - 0.5 );
  m_h_Nbkg_ST     = new TH1D( "Nbkg_ST",     "Number of background events (Single top)", nbins, -0.5, bmax - 0.5 );
  m_h_Nbkg_Z      = new TH1D( "Nbkg_Z",      "Number of background events (Z+jets)", nbins, -0.5, bmax - 0.5 );
  m_h_Nbkg_QCD    = new TH1D( "Nbkg_QCD",    "Number of background events (QCD)",    nbins, -0.5, bmax - 0.5 );
  m_h_Nbkg_Others = new TH1D( "Nbkg_Others", "Number of background events (Others)", nbins, -0.5, bmax - 0.5 );
  m_h_Nsig        = new TH1D( "Nsig",        "Number of signal events",              nbins, -0.5, bmax - 0.5 );
  m_h_Nsig_data   = new TH1D( "Nsig_data",   "Number of signal events (from data)",  nbins, -0.5, bmax - 0.5 );
  m_h_xsec_data   = new TH1D( "xsec",        "XS from data",                         500, 0, 500 );
  
  //m_h_Nbkg_conv_poisson = new TH1D( "Nbkg_conv_poisson", "Number of background events convoluted w/ Poisson", 31, -0.5, 30 );
}


// ---------------------------------------------------------


double ComBAT::LogLikelihood(std::vector <double> parameters)
{
	// This methods returns the logarithm of the conditional probability
	// p(data|parameters). This is where you have to define your model.

	double logprob  = 0.;

	// reset global variables
	m_Nbkg       = 0;
	m_NbkgW      = 0;
	m_NbkgST     = 0;
	m_NbkgZ      = 0;
	m_NbkgQCD    = 0;
	m_NbkgOthers = 0;
	m_Nsig       = 0;
	m_Nsig_data  = 0;
	

	/////////////////////////////
	// Likelihood
	/*	cout << "parameters: ";
	for( unsigned int n = 0; n < parameters.size() ; ++n ) cout << parameters[n] << " ";
	cout << endl;
	*/

	// Poisson likelihood
	logprob += MarginalLikelihood( DATAPOINT::ELE_0, parameters );
	logprob += MarginalLikelihood( DATAPOINT::MU_0,  parameters );

	return logprob;
}


// ---------------------------------------------------------


double ComBAT::MarginalLikelihood( nparam_t datapoint, const std::vector <double>& parameters )
{
  //cout << "Marginal likelihood" << endl;

  string ch   = ( datapoint == DATAPOINT::ELE_0 ) ? "ele" : "mu";  
  BCDataPoint* data = GetDataPoint(datapoint);

  const double xsec_0        = parameters[ 0 ]; //m_paramIndex["xsec"] ];
  const double br            = 0.543;
  const double NsigMC        = data->GetValue( VAL::NsigMC );
  const double Nobs          = data->GetValue( VAL::Nobs );

  //const double f_b           = data->GetValue( VAL::W_b_frac ); //      * ( 1 + CalculateTotalVariation( "f_b", parameters) );
  //const double f_c           = data->GetValue( VAL::W_c_frac ); //      * ( 1 + CalculateTotalVariation( "f_c", parameters) );
  //const double eff_b         = data->GetValue( VAL::eff_btag ); //      * ( 1 + CalculateTotalVariation( "eff_b", parameters) );
  //const double eff_mistag_c  = data->GetValue( VAL::eff_mistag_c ); //  * ( 1 + CalculateTotalVariation( "eff_mistag_c", parameters) );
  //const double eff_mistag_lq = data->GetValue( VAL::eff_mistag_lq ); // * ( 1 + CalculateTotalVariation( "eff_mistag_lq", parameters) );
  //const double eff_btag      = ( f_b * eff_b ) +  ( f_c * eff_mistag_c ) + ( (1 - f_c - f_b ) * eff_mistag_lq ); 
	                            
  const double NbkgW_b    = data->GetValue( VAL::NbkgW_datadriven) * ( 1 + CalculateTotalVariation( "NbkgW_" + ch, parameters) );
	
  const double NbkgQCD    = data->GetValue( VAL::NbkgQCD_datadriven) * ( 1 + CalculateTotalVariation( "NbkgQCD_" + ch, parameters) );
  const double NbkgST     = data->GetValue( VAL::NbkgST)             * ( 1 + CalculateTotalVariation( "NbkgST_" + ch, parameters) );
  const double NbkgZ      = data->GetValue( VAL::NbkgZ)              * ( 1 + CalculateTotalVariation( "NbkgZ_" + ch, parameters) );
  const double NbkgOthers = data->GetValue( VAL::NbkgOthers);//         * ( 1 + CalculateTotalVariation( "NbkgOthers_" + ch, parameters) );
  const double Nbkg       = NbkgW_b + NbkgQCD + NbkgST + NbkgZ + NbkgOthers;

  //const double eff_sig   = data->GetValue( VAL::eff_sig ) * ( 1 + CalculateTotalVariation( "eff_sig_" + ch , parameters ) );
  
  const double iLumi     = data->GetValue( VAL::iLumi );// * ( 1 + CalculateTotalVariation("ilumi", parameters) );
  const double Ngen      = 160 * br * iLumi;
  const double eff_sig   = NsigMC / Ngen;
  const double Nsig      = xsec_0 * iLumi * br * eff_sig;
  const double Nsig_data = Nobs - Nbkg;
  //const double xsec_data = Nsig_data / ( eff_sig * iLumi * br );
  //cout << Nobs << " " << Nbkg << " " << NsigMC << " " << eff_sig << " " << Nsig_data << endl;
  // update global variables
  m_Nbkg        += Nbkg;
  m_NbkgW       += NbkgW_b;
  m_NbkgST      += NbkgST;
  m_NbkgZ       += NbkgZ;
  m_NbkgQCD     += NbkgQCD;
  m_NbkgOthers  += NbkgOthers;
  m_Nsig        += Nsig;
  m_Nsig_data   += Nsig_data;

  //logprob += log( mu_mu )  * Nobs_mu  - mu_mu  - BCMath::ApproxLogFact( int(Nobs_mu) );
  //printf("Nobs_ele=%5.0f Nsig_ele=%5.0f Nbkg_ele=%5.0f mu_ele=%5.0f\n", Nobs_ele, Nsig_ele, Nbkg_ele, mu_ele );
  //cout << "channel " << ch << " logLikelihood=" << BCMath::LogPoisson( Nobs, Nsig + Nbkg ) << endl;
  return BCMath::LogPoisson( Nobs, Nsig + Nbkg );
}

// ---------------------------------------------------------


double ComBAT::LogAPrioriProbability(std::vector <double> parameters)
{
	// This method returns the logarithm of the prior probability for the
	// parameters p(parameters).
  //cout << "log prior" << endl;

	double logprob = 0.;

	// For flat prior it's very easy.
//	for(unsigned int i=0; i < this -> GetNParameters(); i++)

	if( m_xsec.prior == FlatPrior ) {
	  logprob -= log( GetParameter("xsec")->GetRangeWidth() ); // flat a priori on xs
	}
	else if( m_xsec.prior == GaussPrior ) { 
	  logprob += BCMath::LogGaus( parameters[0], GetParameter("xsec")->GetRangeWidth()/2., GetParameter("xsec")->GetRangeWidth()/4. );
	}
	else if( m_xsec.prior == GammaPrior ) {
	  logprob += 0;
	}


	//for( unsigned int s = 1 ; s <= parameters.size() ; ++s ) {
	for( ArrayOfSystematics::iterator itr = m_systematics.begin() ; itr != m_systematics.end() ; ++itr ) { 
	  Systematic * p_syst = &itr->second; //&m_systematics[s-1];
	  unsigned int s =  m_paramIndex[ p_syst->name ];
	  const TypeOfPrior prior = p_syst->prior;
	  //cout << "prior for " << p_syst->name << endl;

	  if( prior == GaussPrior ) {
	    //cout << "gauss" << endl;
	    //const double sigma = p_syst->params["sigma"];
	    logprob += BCMath::LogGaus( parameters[s], 0., 1. );
	  }
	  else if( prior == GammaPrior ) { 
	    //cout << "gamma" << endl;
	    const double x = parameters[s];
	    const double sigma = p_syst->params["sigma"];
	
	    logprob += LogGamma( 1, sigma, x );
	  }
	  else {
	    logprob -= log( GetParameter(s)->GetRangeWidth() );
	  }
	}


	return logprob;
}


//---------------------------------------------------------


double ComBAT::CalculateTotalVariation( string param, const std::vector <double>& parameters ) 
{
  //cout << "variation for " << param << endl;
  double v = 0.;
  //cout << param << ": ";
  
  for( ArrayOfSystematics::iterator itr = m_sigma[param].begin() ; itr != m_sigma[param].end() ; ++itr ) { 
    Systematic * p_syst = &itr->second; 
    const unsigned int p = m_paramIndex[ p_syst->name ]; 
    //cout << "delta " << p_syst->name << " " << endl;

    double delta = 0.0;
    double a, b, c;

    const double sigma_p = p_syst->upperBound;
    const double sigma_m = p_syst->lowerBound;
    //cout << p_syst->name << " s_p=" << sigma_p << "  s_m=" << sigma_m << " * " << parameters[p] << endl;

    if( m_systematics[param].prior == GammaPrior ) {
      delta =  parameters[p] - 1;
    }
    else {
      switch( p_syst->shape ) 
	{
	case Unspecified:
	  cout << "Unspecified" << endl;
	  delta = 0.0;
	  break;
	case Symmetric:
	  delta = sigma_p * parameters[p];
	  break;
	case AsymmetricParabolic:
	  FitParabola( sigma_m, 0, sigma_p, &a, &b, &c);
	  cout << a << " " << b << " " << c << " " << endl;
	  delta = c + b*parameters[p] + a*pow(parameters[p],2);
	  //cout << "param " << param << ": a=" << a << " b=" << b << " c=" << c << endl;
	  break;
	case Asymmetric2HalfGaussian:
	  if( parameters[p] > 0. ) {
	    delta = sigma_p * parameters[p] ;
	  }
	  else {
	    delta = -sigma_m * parameters[p] ;
	  }
	  break;
	case AsymmetricDiscrete:
	  delta = sigma_p * parameters[p];
	  break;
	default:
	  delta = sigma_p * parameters[p];
	  break;
	}
    }
   

    v += delta;
  }
  //cout << endl;
  //cout << param << " v=" << v << endl << endl;;
  return v;
}


// ---------------------------------------------------------

void ComBAT::MCMCUserIterationInterface()
{
  //cout << m_Nbkg << " " << m_NbkgW << endl;

  m_h_Nbkg->Fill( m_Nbkg );
  m_h_Nbkg_W->Fill( m_NbkgW );
  m_h_Nbkg_ST->Fill( m_NbkgST );
  m_h_Nbkg_Z->Fill( m_NbkgZ );
  m_h_Nbkg_QCD->Fill( m_NbkgQCD );
  m_h_Nbkg_Others->Fill( m_NbkgOthers );
  m_h_Nsig->Fill( m_Nsig ); 
  m_h_Nsig_data->Fill( m_Nsig_data );

  m_h_xsec_data->Fill(MCMCGetx().at(0));
}

// ---------------------------------------------------------


void ComBAT::SaveUserDefinedHistograms()
{
  m_h_Nbkg->Write();
  m_h_Nbkg_W->Write();
  m_h_Nbkg_ST->Write();
  m_h_Nbkg_Z->Write();
  m_h_Nbkg_QCD->Write();
  m_h_Nbkg_Others->Write();
  m_h_Nsig->Write();
  m_h_Nsig_data->Write();
  m_h_xsec_data->Write();

  m_histoFile->Close();
  printf("User-defined histograms saved");
}
