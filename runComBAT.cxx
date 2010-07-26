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

	}
	else if( prior == "jeffrey" ) {

	}
	else if( prior == "flat" ) {
	  xsec->Attribute( "min", &param.lowerBound );
	  xsec->Attribute( "max", &param.upperBound );
	}
	else {
	  cout << "Wrong prior " << prior << endl;
	  exit(1);
	}

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
	  BCDataPoint * dp = new BCDataPoint( 6 );
	  
	  double Nobs, NsigMC, NbkgMC, eff_sig, eff_bkg, ilumi = 0.0;
	  dataXML->Attribute( "Nobs",    &Nobs    );
	  dataXML->Attribute( "NsigMC",  &NsigMC  );
	  dataXML->Attribute( "NbkgMC",  &NbkgMC  );
	  dataXML->Attribute( "eff_sig", &eff_sig );
	  dataXML->Attribute( "eff_bkg", &eff_bkg );
	  dataXML->Attribute( "ilumi",   &ilumi   );

	  dp->SetValue( 0, Nobs );
	  dp->SetValue( 1, NsigMC );
	  dp->SetValue( 2, NbkgMC );
	  dp->SetValue( 3, eff_sig );
	  dp->SetValue( 4, eff_bkg );
	  dp->SetValue( 5, ilumi );

	  dataSet->AddDataPoint( dp );

	  cout << "\nData point:  Nobs=" << Nobs << " NsigMC=" << NsigMC << " NbkgMC=" << NbkgMC << " eff_s=" << eff_sig 
	     << " eff_b=" << eff_bkg << " iLumi=" << ilumi << endl;
	
	  cout << "Calculated xs from MC: " << NsigMC /( eff_sig * ilumi * 0.545 ) << " pb" << endl;
	  cout << "Calculated xs from data: " << (Nobs - NbkgMC) / ( eff_sig * ilumi * 0.545 ) << " pb" << endl << endl;
	}
	m->SetDataSet( dataSet );
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

	delete cardFile;
	//delete o;
	delete m;
	//delete ch_ele;
	//delete ch_mu;
	//delete m;

	BCLog::OutSummary("Test program ran successfully");
	BCLog::OutSummary("Exiting");

	return 0;

}

