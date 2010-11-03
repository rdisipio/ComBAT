// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project TopComb
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************
#include "TopComb.h"
// ---------------------------------------------------------
TopComb::TopComb() : BCModel()
{  // default constructor
	DefineParameters();
};

// ---------------------------------------------------------
TopComb::TopComb(const char * name) : BCModel(name)
{  // constructor
        ReadCards(name);
	DefineParameters();
        hxsec = new TH1D("hxsec","xsec",500,0.,500.);
        heff_e = new TH1D("heff_e","eff e",500,0.,1.);
        heff_m = new TH1D("heff_m","eff mu",500,0.,1.);
        hnsig = new TH1D("hnsig","Nsig",600,0.,60.);
        hnbg = new TH1D("hnbg","Nbg",1000,0.,100.);
        hnbgw = new TH1D("hnbgw","Nbg W",600,0.,60.);
        hnbgwnob = new TH1D("hnbgwnob","Nbg W nob",600,0.,200.);
        hnbgqcd = new TH1D("hnbgqcd","Nbg QCD",600,0.,60.);
        hnbgst = new TH1D("hnbgst","Nbg single top",600,0.,60.);
        hnbgz = new TH1D("hnbgz","Nbg Z",600,0.,60.);
        hnsige = new TH1D("hnsige","Nsig e",600,0.,60.);
        hnbge = new TH1D("hnbge","Nbg e",600,0.,60.);
        hnbgwe = new TH1D("hnbgwe","Nbg W e",600,0.,60.);
        hnbgwnobe = new TH1D("hnbgwnobe","Nbg W nob e",600,0.,200.);
        hnbgqcde = new TH1D("hnbgqcde","Nbg QCD e",600,0.,60.);
        hnbgqcd_we = new TH2D("hnbgqcd_we","Nbg QCD vs W e",600,0.,60.,600,0.,60.);
        hnbgste = new TH1D("hnbgste","Nbg single top e",600,0.,60.);
        hnbgze = new TH1D("hnbgze","Nbg Z e",600,0.,60.);
        hnsigm = new TH1D("hnsigm","Nsig mu",600,0.,60.);
        hnbgm = new TH1D("hnbgm","Nbg mu",600,0.,60.);
        hnbgwm = new TH1D("hnbgwm","Nbg W mu",600,0.,60.);
        hnbgwnobm = new TH1D("hnbgwnobm","Nbg W nob mu",600,0.,200.);
        hnbgqcdm = new TH1D("hnbgqcdm","Nbg QCD mu",600,0.,60.);
        hnbgqcd_wm = new TH2D("hnbgqcd_wm","Nbg QCD vs W m",600,0.,60.,600,0.,60.);
        hnbgstm = new TH1D("hnbgstm","Nbg single top mu",600,0.,60.);
        hnbgzm = new TH1D("hnbgzm","Nbg Z mu",600,0.,60.);
        hQCDSYSTE =  new TH1D("hQCDSYSTE","QCDSYSTE",200,-10.,10.);
        hQCDSTATE = new TH1D("hQCDSTATE","QCDSTATE",200,0.,20.);
        hWBGSTATE  = new TH1D("hWBGSTATE","WBGSTATE",200,0.,20.);
        hWBGSYST = new TH1D("hWBGSYST","WBGSYST",200,-10.,10.);
        hW2TAGSTAT = new TH1D("hW2TAGSTAT","W2TAGSTAT",200,0.,20.);
        hBTAGEFF = new TH1D("hBTAGEFF","BTAGEFF",200,-10.,10.);
        hCTAGEFF = new TH1D("hCTAGEFF","BCAGEFF",200,-10.,10.);
        hLTAGEFF = new TH1D("hLTAGEFF","BLAGEFF",200,-10.,10.);
        hWHFRAC = new TH1D("hWHFRAC","WHFRAC",200,-10.,10.);
        hWCFRAC = new TH1D("hWCFRAC","WCFRAC",200,-10.,10.);
        hZXSEC  = new TH1D("hZXSEC","ZXSEC",200,-10.,10.);
        hLUMI  = new TH1D("hLUMI","LUMI",200,-10.,10.);
        hJES  = new TH1D("hJES","JES",200,-10.,10.);
        hJER  = new TH1D("hJER","JER",200,-10.,10.);
        hJERECO  = new TH1D("hJERECO","JERECO",200,-10.,10.);
        hSTXS  = new TH1D("hSTXS","STXS",200,-10.,10.);
        hEID   = new TH1D("hEID","EID",200,-10.,10.);
        hEES   = new TH1D("hEES","EES",200,-10.,10.);
        hEER   = new TH1D("hEER","EER",200,-10.,10.);
        hIFSR   = new TH1D("hIFSR","IFSR",200,-10.,10.);
        hNLOGEN  = new TH1D("hNLOGEN","NLOGEN",200,-10.,10.);
        hPSGEN  = new TH1D("hPSGEN","PSGEN",200,-10.,10.);
        hPDF  = new TH1D("hPDF","PDF",200,-10.,10.);
        hOTX  = new TH1D("hOTX","OTX",200,-10.,10.);
        hQCDSYSTM = new TH1D("hQCDSYSTM","QCDSYSTM",200,-10.,10.);
        hQCDSTATM = new TH1D("hQCDSTATM","QCDSTATM",200,0.,20.);
        hWBGSTATM = new TH1D("hWBGSTATM","WBGSTATM",200,0.,20.);
        hMID   = new TH1D("hMID","MID",200,-10.,10.);
        hMES   = new TH1D("hMES","MES",200,-10.,10.);
        hMER   = new TH1D("hMER","MER",200,-10.,10.);                
        like_xsec = new TH1D("like_xsec","likelihood vs xsec",500,-0.5,499.5);
	// data points
	/*
	map<string,vector<double> >::iterator itnobs;
        DNobs = new BCDataPoint(2);
        for(int ids=0; ids<m_parsyst.size(); ids++){
          itnobs=m_data.at(ids).begin();
          DNobs->SetValue(ids,itnobs->second.at(0));
	}
        SetDataBoundaries(0,0.,100.);
        SetDataBoundaries(1,0.,100.);
	*/
        nneg=0;
};

// ---------------------------------------------------------
TopComb::~TopComb()
{

};  // default destructor

// ---------------------------------------------------------
void TopComb::LikeXsec(){
  std::vector<double> bestfit=GetBestFitParameters();
  std::vector<double> vv;
  std::vector<double>::iterator it;
  for(int i=1; i<bestfit.size(); i++){
    vv.push_back(bestfit.at(i));
  }
  it=vv.begin();
  cout << " Minimum of the likelihood... " << LogLikelihood(bestfit) << endl;
  for(int i=0; i<500; i++){
    double x=i;
    vv.insert(it,x);
    vv.erase(it+1);
    double w=LogLikelihood(bestfit)-LogLikelihood(vv);
    if(i % 100 == 0)
    cout << x << " " << w << endl;
    like_xsec->Fill(x,w);   
  }  
  hpvalue = CalculatePValue(bestfit,true)->GetHistogram();
  hpvalue->SetName("hpvalue");
  hpvalue->SetTitle("P value");
  cout << " P value......." << GetPValue() << endl;
}
void TopComb::ReadCards(const char* file)
{
// read control cards
  map<string,double>::iterator it;
  map<string,map<string,vector<double> > >::iterator it2;
  map<string,vector<double> >::iterator it1;
  mydata = new InputData(file);
  if(mydata->Read()) {
    m_meastyp = *mydata->GetMeasTyp();
    m_data = *mydata->GetData();
    m_measpar = *mydata->GetMeasPar();
    m_measlim = *mydata->GetMeasLim();
    m_parcv = *mydata->GetparCV();
    m_parsyst = *mydata->GetparSyst();
    m_systintlim =  *mydata->GetSystIntLim();
    m_systtyp =  *mydata->GetSystTyp();
    m_gammapar = *mydata->GetGammaPar();

    likelihood = *mydata->GetLikelihood();
  }
  mydata->PrintparSyst();
  mydata->PrintparCV();
  mydata->PrintSystIntLim();
  mydata->PrintSystTyp();
};  

void TopComb::SetNobs()
{
    m_data = *mydata->GetData();
}

InputData* TopComb::GetCards(){
  return mydata;
}

// ---------------------------------------------------------
void TopComb::DefineParameters()
{
	// Add parameters to your model here.
	// You can then use them in the methods below by calling the
	// parameters.at(i) or parameters[i], where i is the index
	// of the parameter. The indices increase from 0 according to the
	// order of adding the parameters.

//	this -> AddParameter("par1", 0.0, 1.0);   // index 0
//	this -> AddParameter("par2", -7.0, 28.0); // index 1

  map<string,vector<double> >::iterator itmeas;
  int id=0;
  itmeas=m_measlim.begin();

  cout << "Added parameter " << id << " " << itmeas->first << " " << itmeas->second.at(0) << " " << itmeas->second.at(1) << endl;
  AddParameter(itmeas->first.c_str(), itmeas->second.at(0), itmeas->second.at(1)); // index 0
  map<string, vector<int> >::iterator it;
  id++;
  double x1,x2;
  while(id <= m_systtyp.size()){
    for(it=m_systtyp.begin(); it!=m_systtyp.end(); it++){
      if(it->second.at(0) == id){
        x1=m_systintlim.find(it->first)->second.at(0);
        x2=m_systintlim.find(it->first)->second.at(1);
        AddParameter(it->first.c_str(), x1, x2);
        cout << "Added parameter " << id << " " << it->first << " " << x1 << " " << x2 << endl;
      }
    }
    id++;
  }
}

// ---------------------------------------------------------
double TopComb::LogLikelihood(std::vector <double> parameters)
{
	// This methods returns the logarithm of the conditional probability
	// p(data|parameters). This is where you have to define your model.

        const double pig=3.14159;
	double logprob = 0.;
        double par0, systdev;
        double Nsig, mu, Nobs;
        vector<double> parx;
        map<string,map<string,vector<double> > >::iterator it1;
	map<string,vector<double> >::iterator it2, itnobs;
	map<string,double>::iterator it3;
	map<string,vector<int> >::iterator it4;
        double a,b;
	bool gamma=false;
        Nsig_tot=0;
        Nbg_tot=0;
        Nbgw_tot=0;
        Nbgwnob_tot=0;
        Nbgqcd_tot=0;
        Nbgst_tot=0;
        Nbgz_tot=0;
        Nsig_e=0;
        Nbg_e=0;
        Nbgw_e=0;
        Nbgwnob_e=0;
        Nbgqcd_e=0;
        Nbgst_e=0;
        Nbgz_e=0;
        Nsig_m=0;
        Nbg_m=0;
        Nbgw_m=0;
        Nbgwnob_m=0;
        Nbgqcd_m=0;
        Nbgst_m=0;
        Nbgz_m=0;
        testg=-1;
        testg2=-1;
        sys.clear();
	// loop over different datasets
        for (int ids=0; ids<m_parsyst.size(); ids++){
    	  // loop on parameters
          for(it3=m_parcv.at(ids).begin(); it3!=m_parcv.at(ids).end(); it3++){
            par0=it3->second;
            it1=m_parsyst.at(ids).find(it3->first);
            systdev=0;
	    //            gamma=false;
	    // loop on systematics
            for(it2=it1->second.begin(); it2!=it1->second.end(); it2++)
	    {
              if(it2->second.at(0) != 0){
                it4=m_systtyp.find(it2->first);
                if(it4 != m_systtyp.end()){
		  //		  cout << "DEBUG    " << it3->first << " " << it2->first << " " << it2->second.at(0) << " " << it4->first << " " << it4->second.at(0) << " " << it4->second.at(1) << " " << par0 << endl;
		  // Flat or Gaussian Syst. Prior
                  if(it4->second.at(1) < 2){
                    if(it2->second.at(1) == 0){
                      systdev += it2->second.at(0)*parameters[it4->second.at(0)];
		    }
                    else{
		      parabola(it2->second,&a,&b);
		      systdev += a*pow(parameters[it4->second.at(0)],2)+b*parameters[it4->second.at(0)];
		    }
		  }
		  else if(it4->second.at(1) == 2){
		    // Gamma Syst. Prior
		    systdev += parameters[it4->second.at(0)]-1;
                    if(ids==0){
		      testg=(parameters[it4->second.at(0)]);
		      //		      if(parameters[it4->second.at(0)]>3)cout << "ids/par0   " << ids << " " << it4->second.at(0) << " " << parameters[it4->second.at(0)] << " " << par0 << endl;
		    }
		    if(ids==1){
		      testg2=(parameters[it4->second.at(0)]);
		      //		      if(parameters[it4->second.at(0)]>3)cout << "ids/par0   " << ids << " " << it4->second.at(0) << " " << parameters[it4->second.at(0)] << " " << par0 << endl;
		    }		    //                    gamma=true;
		  }
		  /*
                  else if(it4->second.at(1) == 3){
                    systdev += pow(it2->second.at(0),parameters[it4->second.at(0)]);
		  }
		  */
                }
	      }
	    }
            if(!gamma) {
	      if(par0*(1+systdev)>0)parx.push_back(par0*(1+systdev));
	      else parx.push_back(0.);
	    }
	    else {
	      if((par0*systdev)>0)parx.push_back(par0*(systdev-1));
	      else parx.push_back(0.);
	    }
            gamma=false;
            systdev=0;
	  }
	  //          for(int id=1; id<GetNParameters(); id++){
	  //	    sys.push_back(parameters[id]);            
	  //	  }
          itnobs=m_data.at(ids).begin();
	  Nsig=parameters[0]*parx.at(0)*parx.at(1)*parx.at(2);
	  //	  Nsig=0.;
          mu=Nsig+parx.at(3)+parx.at(4)+parx.at(5)+parx.at(6);
	  //          Nobs=itnobs->second.at(0);
	  //          Nobs=DNobs->GetValue(ids);
          Nobs=GetDataSet()->GetDataPoint(0)->GetValue(ids);
	  //	  cout << " ===================== " << Nobs << endl;
	  //	  cout << ids << " " << mu << " " << Nsig << " " << Nobs << endl;
	  //	  cout << parameters[0] << " " << parx.at(0) << " " << parx.at(1) << " " << parx.at(2) << endl;
	  //	  cout << parx.at(3) << " " << parx.at(4) << " " << parx.at(5) << " " << parx.at(6) << " " << parx.at(7) << " " << parx.at(8) << endl;
          Nsig_tot+=Nobs-(parx.at(3)+parx.at(4)+parx.at(5)+parx.at(6));
          Nbg_tot+=parx.at(3)+parx.at(4)+parx.at(5)+parx.at(6);
          Nbgw_tot+=parx.at(5);
	  if(parx.at(3) <0 || parx.at(4)<0 || parx.at(5)<0)
            cout << "parx " << parx.at(3)<< " " << parx.at(4) << " " << parx.at(5) << endl;
          Nbgwnob_tot+=0;
          Nbgqcd_tot+=parx.at(3);
          Nbgst_tot+=parx.at(4);
          Nbgz_tot+=parx.at(6);
          if(ids==0){
            Nbg_e=parx.at(3)+parx.at(4)+parx.at(5)+parx.at(6); 
            Nsig_e=Nobs-(parx.at(3)+parx.at(4)+parx.at(5)+parx.at(6));
            Nbgw_e=parx.at(5);
            Nbgwnob_e=0.;
            Nbgqcd_e=parx.at(3);
            Nbgst_e=parx.at(4);
            Nbgz_e=parx.at(6);
          }
          if(ids==1){
            Nbg_m=parx.at(3)+parx.at(4)+parx.at(5)+parx.at(6); 
            Nsig_m=Nobs-(parx.at(3)+parx.at(4)+parx.at(5)+parx.at(6));
            Nbgw_m=parx.at(5);
            Nbgwnob_m=0.;
            Nbgqcd_m=parx.at(3);
            Nbgst_m=parx.at(4);
            Nbgz_m=parx.at(6);
          }
	  if(likelihood.at(0) == 1){
 	    logprob +=  log(mu)*Nobs -mu -BCMath::ApproxLogFact(int(Nobs));
	    //	    cout << ids << " " << Nsig << " " << parx.at(3)+parx.at(4)+parx.at(5)+parx.at(6) << " " << Nobs << " " << log(mu)*Nobs -mu -BCMath::ApproxLogFact(int(Nobs)) << " " << BCMath::LogPoisson(Nobs,mu) << endl;
	  }
	  else if(likelihood.at(0) == 2 )
	    logprob +=  -log(sqrt(2*pig*mu)) -pow((mu-Nobs),2)/(2*mu);
	  parx.clear();
	}
	//	if(CheckforNegativeValues())
	return logprob;
	//        else return -1.e999;
	//	return 1.;
}

// ---------------------------------------------------------
double TopComb::LogAPrioriProbability(std::vector <double> parameters)
{
	// This method returns the logarithm of the prior probability for the
	// parameters p(parameters).

	double logprob = 0.;

	// Prior cross section
        map<string, vector<int> >::iterator itm;
        map<string, vector<double> >::iterator itmp;
        itm=m_meastyp.begin();
        itmp=m_measpar.begin();
	//	cout << "Debug xsec prior " << itm->first << " " <<  itm->second.at(0) << " " <<  GetParameter(0) -> GetRangeWidth() << endl;
        if(itm->second.at(0)==0)
        logprob -= log(GetParameter(0) -> GetRangeWidth());
        if(itm->second.at(0)==1 )
	  logprob -= pow((parameters.at(0)-itmp->second.at(0)),2)/(2*pow(itmp->second.at(1),2));

        map<string, vector<int> >::iterator it;
        map<string, vector<double> >::iterator itg;
        for(it=m_systtyp.begin(); it!=m_systtyp.end(); it++){
	  // Flat prior
          if(it->second.at(1)==0)
          logprob -= log(GetParameter(it->second.at(0)) -> GetRangeWidth());
	  // Gaussian Prior
          if(it->second.at(1)==1)
	    logprob -= pow(parameters.at(it->second.at(0)),2)/2;
	  // Gamma prior
          if(it->second.at(1)==2){
	    itg=m_gammapar.find(it->first);
	    //            cout << "test test test " << it->first << " " << itg->second.at(0) << " " << it->second.at(0) << " " << parameters.at(it->second.at(0)) << endl;
	    logprob += Gamma(1.,itg->second.at(0),parameters.at(it->second.at(0)));
	    //	    logprob += Gamma(itg->second.at(0),itg->second.at(1),parameters.at(it->second.at(0)));
	  }
	}
	// For flat prior it's very easy.
//	for(unsigned int i=0; i < this -> GetNParameters(); i++)
//		logprob -= log(this -> GetParameter(i) -> GetRangeWidth());

	return logprob;
}
// ---------------------------------------------------------

double TopComb::Gamma(double mu, double sig, double x)
{
  //  return -pow((x-mu)/sig,2)/2;
  //  double k=pow(mu/sig,2);
  //  double theta=pow(sig,2)/mu;
  double k=1+(pow(mu,2)+mu*sqrt(pow(mu,2)+4*pow(sig,2)))/2/pow(sig,2);
  double theta=mu/(k-1);
  return log(TMath::GammaDist(x,k,0,theta));
}

void TopComb::parabola(std::vector <double> sig, double* a, double* b)
{
  *a=(sig.at(0)+sig.at(1))/2.;
  *b=(sig.at(0)-sig.at(1))/2.;
}

bool TopComb::CheckforNegativeValues()
{  
  if(Nbgqcd_e < 0 || Nbgqcd_m < 0 || Nbgst_e < 0 || Nbgst_m < 0 || Nbgz_e < 0 || Nbgz_m < 0 || Nbgw_e < 0 || Nbgw_m < 0 || Nsig_m < 0 ||  Nsig_e < 0){ 
    //     nneg++; 
     //     cout << " !!!! Negative Value !!!! " << nneg << endl; 
     return false; 
  }
  return true;
}

void TopComb::MCMCUserIterationInterface()
{
  MCMCrun();
  if(!CheckforNegativeValues()){
    //   cout << " negative!!! " << endl; 
   return;
  }
  hnsig->Fill(Nsig_tot);
  hnbg->Fill(Nbg_tot);
  if(Nbg_tot < 0.0001){
    cout << Nbg_tot << endl;
    cout << Nbgw_tot << " " << Nbgqcd_tot << " " << Nbgst_tot << " " << Nbgz_tot << endl;
  }
  heff_e->Fill(eff_e);
  heff_m->Fill(eff_m);
  hnbgw->Fill(Nbgw_tot);
  hnbgwnob->Fill(Nbgwnob_tot);
  hnbgqcd->Fill(Nbgqcd_tot);
  hnbgst->Fill(Nbgst_tot);
  hnbgz->Fill(Nbgz_tot);
  hnsige->Fill(Nsig_e);
  hnbge->Fill(Nbg_e);
  hnbgwe->Fill(Nbgw_e);
  hnbgwnobe->Fill(Nbgwnob_e);
  hnbgqcde->Fill(Nbgqcd_e);
  hnbgqcd_we->Fill(Nbgqcd_e,Nbgw_e);
  hnbgste->Fill(Nbgst_e);
  hnbgze->Fill(Nbgz_e);
  hnsigm->Fill(Nsig_m);
  hnbgm->Fill(Nbg_m);
  hnbgwm->Fill(Nbgw_m);
  hnbgwnobm->Fill(Nbgwnob_m);
  hnbgqcdm->Fill(Nbgqcd_m);
  hnbgqcd_wm->Fill(Nbgqcd_m,Nbgw_m);
  hnbgstm->Fill(Nbgst_m);
  hnbgzm->Fill(Nbgz_m);
  hxsec->Fill(MCMCGetx().at(0));
  hQCDSYSTE->Fill(MCMCGetx().at(1));
  hQCDSTATE->Fill(MCMCGetx().at(2));
  hWBGSTATE->Fill(MCMCGetx().at(3));
  hWBGSYST->Fill(MCMCGetx().at(4));
  hW2TAGSTAT->Fill(MCMCGetx().at(5));
  hBTAGEFF->Fill(MCMCGetx().at(6));
  hCTAGEFF->Fill(MCMCGetx().at(7));
  hLTAGEFF->Fill(MCMCGetx().at(8));
  hWHFRAC->Fill(MCMCGetx().at(9));
  hWCFRAC->Fill(MCMCGetx().at(10));
  hZXSEC->Fill(MCMCGetx().at(11));
  hLUMI->Fill(MCMCGetx().at(12));
  hJES->Fill(MCMCGetx().at(13));
  hJER->Fill(MCMCGetx().at(14));
  hJERECO->Fill(MCMCGetx().at(15));
  hSTXS->Fill(MCMCGetx().at(16));
  hEID->Fill(MCMCGetx().at(17));
  hEES->Fill(MCMCGetx().at(18));
  hEER->Fill(MCMCGetx().at(19));
  hIFSR->Fill(MCMCGetx().at(20));
  hNLOGEN->Fill(MCMCGetx().at(21));
  hPSGEN->Fill(MCMCGetx().at(22));
  hPDF->Fill(MCMCGetx().at(23));
  hOTX->Fill(MCMCGetx().at(24));
  hQCDSYSTM->Fill(MCMCGetx().at(25));
  hQCDSTATM->Fill(MCMCGetx().at(26));
  hWBGSTATM->Fill(MCMCGetx().at(27));
  hMID->Fill(MCMCGetx().at(28));
  hMES->Fill(MCMCGetx().at(29));
  hMER->Fill(MCMCGetx().at(30));
  //  if(testg > 5 || testg2 > 5)cout << " user " << testg << " " << testg2 << endl;
}
vector<TH2D*> TopComb::Gethistos2()
{
  hvec2.push_back(hnbgqcd_we);
  hvec2.push_back(hnbgqcd_wm);
  return hvec2;
}
vector<TH1D*> TopComb::Gethistos()
{
  hvec.push_back(hnsig);
  hvec.push_back(hnbg);
  hvec.push_back(hnbgw);
  hvec.push_back(hnbgwnob);
  hvec.push_back(hnbgqcd);
  hvec.push_back(hnbgst);
  hvec.push_back(hnbgz);
  hvec.push_back(hnsige);
  hvec.push_back(hnbge);
  hvec.push_back(hnbgwe);
  hvec.push_back(hnbgwnobe);
  hvec.push_back(hnbgqcde);
  hvec.push_back(hnbgste);
  hvec.push_back(hnbgze);
  hvec.push_back(hnsigm);
  hvec.push_back(hnbgm);
  hvec.push_back(hnbgwm);
  hvec.push_back(hnbgwnobm);
  hvec.push_back(hnbgqcdm);
  hvec.push_back(hnbgstm);
  hvec.push_back(hxsec);
  hvec.push_back(hQCDSYSTE);
  hvec.push_back(hQCDSTATE);
  hvec.push_back(hWBGSTATE);
  hvec.push_back(hWBGSYST);
  hvec.push_back(hW2TAGSTAT);
  hvec.push_back(hBTAGEFF);
  hvec.push_back(hCTAGEFF);
  hvec.push_back(hLTAGEFF);
  hvec.push_back(hWHFRAC);
  hvec.push_back(hWCFRAC);
  hvec.push_back(hZXSEC);
  hvec.push_back(hLUMI);
  hvec.push_back(hJES);
  hvec.push_back(hJER);
  hvec.push_back(hJERECO);
  hvec.push_back(hSTXS);
  hvec.push_back(hEID);
  hvec.push_back(hEES);
  hvec.push_back(hEER);
  hvec.push_back(hIFSR);
  hvec.push_back(hNLOGEN);
  hvec.push_back(hPSGEN);
  hvec.push_back(hPDF);
  hvec.push_back(hOTX);
  hvec.push_back(hQCDSYSTM);
  hvec.push_back(hQCDSTATM);
  hvec.push_back(hWBGSTATM);
  hvec.push_back(hMID);
  hvec.push_back(hMES);
  hvec.push_back(hMER);
  hvec.push_back(heff_e);
  hvec.push_back(heff_m);
  hvec.push_back(like_xsec);
  hvec.push_back(hpvalue);
  return hvec;
}

void TopComb::MCMCrun()
{
        const double pig=3.14159;
        double par0, systdev;
        vector<double> parxx;
        map<string,map<string,vector<double> > >::iterator it1;
	map<string,vector<double> >::iterator it2;
	map<string,double>::iterator it3;
	map<string,vector<int> >::iterator it4;
        double a,b;
        Nsig_tot=0;
        Nbg_tot=0;
        Nbgw_tot=0;
        Nbgwnob_tot=0;
        Nbgqcd_tot=0;
        Nbgst_tot=0;
        Nbgz_tot=0;
        Nsig_e=0;
        Nbg_e=0;
        Nbgw_e=0;
        Nbgwnob_e=0;
        Nbgqcd_e=0;
        Nbgst_e=0;
        Nbgz_e=0;
        Nsig_m=0;
        Nbg_m=0;
        Nbgw_m=0;
        Nbgwnob_m=0;
        Nbgqcd_m=0;
        Nbgst_m=0;
        Nbgz_m=0;
	// loop over different datasets
        for (int ids=0; ids<m_parsyst.size(); ids++){
    	  // loop on parameters
          for(it3=m_parcv.at(ids).begin(); it3!=m_parcv.at(ids).end(); it3++){
            par0=it3->second;
            it1=m_parsyst.at(ids).find(it3->first);
            systdev=0;
	    //            gamma=false;
	    // loop on systematics
            for(it2=it1->second.begin(); it2!=it1->second.end(); it2++)
	    {
              if(it2->second.at(0) != 0){
                it4=m_systtyp.find(it2->first);
                if(it4 != m_systtyp.end()){
		  //		  cout << "DEBUG    " << it3->first << " " << it2->first << " " << it2->second.at(0) << " " << it4->first << " " << it4->second.at(0) << " " << it4->second.at(1) << " " << par0 << endl;
		  // Flat or Gaussian Syst. Prior
                  if(it4->second.at(1) < 2){
                    if(it2->second.at(1) == 0){
          	      systdev += it2->second.at(0)*MCMCGetx().at(it4->second.at(0));
		    }
                    else{
		      parabola(it2->second,&a,&b);
		      systdev += a*pow(MCMCGetx().at(it4->second.at(0)),2)+b*MCMCGetx().at(it4->second.at(0));
		    }
		  }
		  else{
		    // Gamma Syst. Prior
		    systdev += MCMCGetx().at(it4->second.at(0))-1;
		  }
                }
	      }
	    }
	    parxx.push_back(par0*(1+systdev));
            systdev=0;
	  }
          for(int id=1; id<GetNParameters(); id++){
	    sys.push_back(MCMCGetx().at(id));            
	  }
	  //          Nsig_tot+=Nobs-(parxx.at(3)+parxx.at(4)+parxx.at(5)+parxx.at(6)+parxx.at(7)+parxx.at(8));
          Nbg_tot+=parxx.at(3)+parxx.at(4)+parxx.at(5)+parxx.at(6);
          Nbgw_tot+=parxx.at(5);
          Nbgwnob_tot+=0;
          Nbgqcd_tot+=parxx.at(3);
          Nbgst_tot+=parxx.at(4);
          Nbgz_tot+=parxx.at(6);
          if(ids==0){
            Nbg_e=parxx.at(3)+parxx.at(4)+parxx.at(5)+parxx.at(6); 
	    //            Nsig_e=Nobs-(parxx.at(3)+parxx.at(4)+parxx.at(5)+parxx.at(6)+parxx.at(7)+parxx.at(8));
            Nbgw_e=parxx.at(5);
            Nbgwnob_e=0.;
            Nbgqcd_e=parxx.at(3);
            Nbgst_e=parxx.at(4);
            Nbgz_e=parxx.at(6);
            eff_e=parxx.at(1);            
          }
          if(ids==1){
            Nbg_m=parxx.at(3)+parxx.at(4)+parxx.at(5)+parxx.at(6); 
	    //            Nsig_m=Nobs-(parxx.at(3)+parxx.at(4)+parxx.at(5)+parxx.at(6)+parxx.at(7)+parxx.at(8));
            Nbgw_m=parxx.at(5);
            Nbgwnob_m=0.;
            Nbgqcd_m=parxx.at(3);
            Nbgst_m=parxx.at(4);
            Nbgz_m=parxx.at(6);
            eff_m=parxx.at(1);
          }
	  parxx.clear();
	}
}
