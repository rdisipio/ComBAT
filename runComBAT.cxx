// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project ComBAT
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include <string>
#include <sstream>
#include <vector>
#include <TRandom3.h>

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCDataSet.h>
#include <BAT/BCDataPoint.h>
#include <BAT/BCH1D.h>
#include <BAT/BCModelOutput.h>

#include "TinyXML/tinyxml.h"

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

	//m->SetCardFileName( configFileName );

	// read configuration
	TiXmlDocument * cardFile = new TiXmlDocument( configFileName.c_str() );
	if( !cardFile->LoadFile() ) {
	  cout << "Impossibile to read cardfile " << configFileName << endl;
	  exit(1);
	}

	TiXmlHandle h_cardFile( cardFile );
	TiXmlHandle analysis( h_cardFile.FirstChildElement("analysis").Element() );
	TiXmlElement * xsec = analysis.FirstChildElement("param").Element() ;
	
	Systematic param;
	param.name = xsec->Attribute("name");

	const string prior = xsec->Attribute("prior");

	if( prior == "gauss" ) {
	  param.type = (TypeOfSystematic)GaussPrior;
	  
	}
	else if( prior == "jeffrey" ) {
	  param.type = (TypeOfSystematic)JeffreyPrior;
	}
	else if( prior == "flat" ) {
	  param.type = (TypeOfSystematic)FlatPrior;
	}
	else {
	  cout << "Wrong prior " << prior << endl;
	  exit(1);
	}
	 xsec->Attribute( "min", &param.lowerBound );
	 xsec->Attribute( "max", &param.upperBound );
	//cout << "Found param " << param.name << " " << param.lowerBound << " / " << param.upperBound << " with prior " << prior << endl;
	m->AddParam( param );


	// get systematics and their parametrization
	for( TiXmlElement * systXML = xsec->FirstChildElement("syst") ; systXML ; systXML = systXML->NextSiblingElement() ) {
	  Systematic syst;
	  
	  syst.name = systXML->Attribute("name");
	  systXML->Attribute("min", &syst.lowerBound );
	  systXML->Attribute("max", &syst.upperBound );
	  
	  int tos = 0;
	  systXML->QueryIntAttribute( "prior", &tos );
	  syst.type = (TypeOfSystematic)tos;

	  m->AddSystematic( syst );
	} // end loop over xsec/syst


	// Now look for variations due to systematic errors
	for( TiXmlElement * deltaXML = analysis.FirstChildElement("delta").Element() ; deltaXML ; deltaXML = deltaXML->NextSiblingElement() ) {
	  if( !deltaXML) break;
	  const char * delta_name = deltaXML->Attribute("name");
	  if( !delta_name ) continue;
	  cout << "delta " << delta_name << endl;

	  for( TiXmlElement * sourceXML = deltaXML->FirstChildElement("sigma") ; sourceXML ; sourceXML = sourceXML->NextSiblingElement() ) {
	    Systematic source;

	    source.name = sourceXML->Attribute("source");
	    sourceXML->Attribute( "up",     &source.upperBound );
	    sourceXML->Attribute( "down",   &source.lowerBound );

	    m->AddSigma( delta_name, source );
	    //cout << source.name << endl;
	  } // end loop over sources of systematics
	} //end loop over delta


	// read data points
	BCDataSet* dataSet = new BCDataSet();
	for( TiXmlElement * dataXML = analysis.FirstChildElement("data").Element() ; dataXML ; dataXML = dataXML->NextSiblingElement() ) {

	  //determine the number of parameters
	  TiXmlNode * child = 0;
	  unsigned int nParams = 1;
	  while( child = ( dataXML->IterateChildren( child )) ) {
	    ++nParams;
	  }

	  BCDataPoint * dp = new BCDataPoint( nParams );
	  
	  double Nobs, NsigMC, NbkgW, NbkgQCD, NbkgST, NbkgZ, NbkgOthers, W_c_frac, W_b_frac, eff_btag, eff_mistag_c, eff_mistag_lq, eff_sig, ilumi = 0.0;

	  Nobs               = atof( dataXML->FirstChildElement("Nobs")->GetText() ); //Attribute( "Nobs",    &Nobs    );
	  NsigMC             = atof( dataXML->FirstChildElement("NsigMC")->GetText() );
	  NbkgW              = atof( dataXML->FirstChildElement("NbkgW")->GetText() );
	  NbkgQCD            = atof( dataXML->FirstChildElement("NbkgQCD")->GetText() );
	  NbkgST             = atof( dataXML->FirstChildElement("NbkgST")->GetText() );
	  NbkgZ              = atof( dataXML->FirstChildElement("NbkgZ")->GetText() );
	  NbkgOthers         = atof( dataXML->FirstChildElement("NbkgOthers")->GetText() );
	  W_c_frac           = atof( dataXML->FirstChildElement("W_c_frac")->GetText() );
	  W_b_frac           = atof( dataXML->FirstChildElement("W_b_frac")->GetText() );
	  eff_btag           = atof( dataXML->FirstChildElement("eff_btag")->GetText() );
	  eff_mistag_c       = atof( dataXML->FirstChildElement("eff_mistag_c")->GetText() );
	  eff_mistag_lq      = atof( dataXML->FirstChildElement("eff_mistag_lq")->GetText() );
	  eff_sig            = atof( dataXML->FirstChildElement("eff_sig")->GetText() );
	  ilumi              = atof( dataXML->FirstChildElement("ilumi")->GetText() );


	  dp->SetValue( VAL::Nobs, Nobs );
	  dp->SetValue( VAL::NsigMC, NsigMC );
	  dp->SetValue( VAL::NbkgW_datadriven, NbkgW );
	  dp->SetValue( VAL::NbkgQCD_datadriven, NbkgQCD );
	  dp->SetValue( VAL::NbkgST, NbkgST );
	  dp->SetValue( VAL::NbkgZ, NbkgZ );
	  dp->SetValue( VAL::NbkgOthers, NbkgOthers ); 
	  dp->SetValue( VAL::W_c_frac, W_c_frac );
	  dp->SetValue( VAL::W_b_frac, W_b_frac );
	  dp->SetValue( VAL::eff_btag, eff_btag );
	  dp->SetValue( VAL::eff_mistag_c, eff_mistag_c );
	  dp->SetValue( VAL::eff_mistag_lq, eff_mistag_lq );
	  dp->SetValue( VAL::eff_sig, eff_sig );
	  dp->SetValue( VAL::iLumi, ilumi );

	  dataSet->AddDataPoint( dp );

	  const double calc_eff_btag = ( W_b_frac * + eff_btag ) +  ( W_c_frac * eff_mistag_c ) + ( (1 - W_c_frac - W_b_frac ) * eff_mistag_lq );
	  const double NbkgW_btag = NbkgW * calc_eff_btag;

	  const double Nbkg = NbkgW_btag + NbkgST + NbkgZ + NbkgQCD + NbkgOthers;
	  cout << "\nData point:  Nobs=" << Nobs << " NsigMC=" << NsigMC << " NbkgW=" << NbkgW 
	       << " Nbkg=" << Nbkg << " eff_btag=" << eff_btag
	       << " eff_s=" << eff_sig << endl;
	     
	
	  cout << "Calculated xs from MC: " << NsigMC /( eff_sig * ilumi * 0.545 ) << " pb" << endl;
	  cout << "Calculated xs from data: " << (Nobs - Nbkg) / ( eff_sig * ilumi * 0.545 ) << " pb" << endl << endl;
	}
	m->SetDataSet( dataSet );
	m->DefineParameters();

	string UDHistoFileName = tokens[0] + "_UDHistograms.root";
	m->SetUserDefinedHistogramFile( UDHistoFileName );

        

	BCLog::OutSummary("Test model created");

	// perform your analysis here

	// normalize the posterior, i.e. integrate posterior
	// over the full parameter space
	//	m -> Normalize();

	// run MCMC and marginalize posterior wrt. all parameters
	// and all combinations of two parameters
	m -> MCMCSetNIterationsBurnIn(10000);
	m->MCMCGetTRandom3()->SetSeed(21340);
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

	string rootFileName = tokens[0] + "_MarkovChains.root";
	BCModelOutput * o = new BCModelOutput(m, rootFileName.c_str() );

	int npar = m->GetNParameters();
        for (int i=0; i<npar; i++){
          BCParameter * a = m->GetParameter(i);
	  std::cout << " 95% C.L. ----> " << m->GetMarginalized(a)->GetQuantile(0.95) << std::endl;
	}

	// print results of the analysis into a text file
	const string txtFileName = tokens[0] + "_results.txt";
	m -> PrintResults( txtFileName.c_str() );

	o->WriteMarkovChain(true);
	o->WriteMarginalizedDistributions();
	o->Close();

	// close log file
	BCLog::CloseLog();

	//delete cardFile;
	//delete o;
	//delete m;
	//delete ch_ele;
	//delete ch_mu;
	//delete m;

	BCLog::OutSummary("Test program ran successfully");
	BCLog::OutSummary("Exiting");

	return 0;

}

