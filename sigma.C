void sigma(double Z)
{
  TH1D* hh = new TH1D("hh","gauss",2000,-10.,10.);
  TF1 *gaus = new TF1("gaus","TMath::Gaus(x,0,1,1)",-10.,10.);
  int Nbins=hh->GetNbinsX();
  Double_t x[0], sum, x0, x1;
  for (int ib=0; ib<Nbins; ib++){
    x[0]=-10.+ib*0.01+0.005;
    hh->Fill(x[0],gaus->EvalPar(x)*0.01);
  }
  if(Z>1)cout << "Z....." << Z << "-----> Prob....." << hh->Integral(1000+Z*100,2000) << endl;
  else{
    sum=0;
    for (int ib=0; ib<Nbins; ib++){
      x0=-10.+ib*0.01+0.005;
      x1=-10.+(ib+1)*0.01+0.005;
      sum+=hh->Integral(1,ib);
      if(Z<(1.-hh->Integral(1,ib)) && Z>(1-hh->Integral(1,ib+1))){
	cout << "Prob....." << Z << " -----> Z....." << (x1+x0)/2 << endl;
        break;
      }
    }
  }
}


////////////////////////////////////////////////////////


Double_t fitf(Double_t *v, Double_t *par)
{
  Double_t x = v[0];
  Double_t val = par[0] + 
    par[1] * x + 
    par[2] * ( pow(x,2) - 1 )  +
    par[3] * ( 4 * pow(x, 3) - 3 * x )  +
    par[4] * ( 8 * pow(x, 4) - 8 * pow(x,2) + 1 );

  return val;
}


////////////////////////////////////////////////////////


void MakeInvertedHistogram( const TH1D* source )
{
  const Double_t max = source->GetMaximum();
  const Int_t binmax = source->GetMaximumBin();
  const Double_t xs  = source->GetXaxis()->GetBinCenter( binmax );

  const Double_t xmax = source->GetXaxis()->GetXmax();
  cout << "X range: 0 - " << xmax << " xs=" << xs << " pb" << endl;

  TH1D * target = new TH1D("inverted", source->GetTitle(), source->GetNbinsX(), 0, xmax/160 );

  for( unsigned int bin = 1 ; bin < source->GetNbinsX() + 1 ; ++bin ) {
    if( max <= 0 || source->GetBinContent(bin) <= 0 ) continue;
    const Double_t val = log(max) - log( source->GetBinContent(bin) );
    target->SetBinContent( bin, val );
  }


  target->Draw();

  TF1 * fit = new TF1( "fit", fitf, 0, 2, 5);
  fit->SetParameters( 1, 1, 1, 1, 1 );
  fit->SetLineColor(kRed);
  fit->SetLineStyle(kDashed);
  target->Fit( "fit" );

  //const Double_t L0 = target->GetBinContent(2);
  const Double_t L0 = fit->Eval( 0 );
  cout << "Value at x=0: " << L0 << endl;

  cout << "Chisq = " << fit->GetChisquare() << " Chisq/Ndof = " << fit->GetChisquare() / fit->GetNDF() << endl;
  
  const Double_t significance = sqrt( 2 * L0 );
  cout << "Significance: " << significance << endl;
}


////////////////////////////////////////////////////////




void doit()
{
  TFile * f050   = TFile::Open("analysis_baseline_rediscovery_0_50_UDHistograms.root");
  TFile * f50500 = TFile::Open("analysis_baseline_rediscovery_UDHistograms.root");

  TH1D * xs1 = (TH1D*)f050->Get("xsec");
  TH1D * xs2 = (TH1D*)f50500->Get("xsec");

  xs1->Scale( 1. / xs1->Integral() );
  xs2->Scale( 1. / xs2->Integral() );
  
  const Double_t n1 = xs1->Integral( 50, 80 );
  const Double_t n2 = xs2->Integral( 50, 80 );
  const Double_t sf = n1/n2;
  cout << "Integral(50,80) of xs1 = " << n1 << ", of xs2 = " << n2 << ", ratio = " << sf << endl;
  xs2->Scale( sf );


  TH1D * hsum = new TH1D( "hsum", "Cross Section PDF", 500, 0, 500 );

  for( unsigned int bin = 1 ; bin <= 51 ; ++bin ) {
    hsum->SetBinContent( bin, xs1->GetBinContent(bin) );
  }

  for( unsigned int bin = 51 ; bin <= 501 ; ++bin ) {
    hsum->SetBinContent( bin, xs2->GetBinContent(bin) );
  }

  //hsum->SetLineStyle(kDashed);
  xs1->SetLineColor(kRed);
  xs2->SetLineColor(kBlue);

  hsum->Draw();
  //xs1->Draw("same");
  //xs2->Draw("same");

  MakeInvertedHistogram( hsum );
}
