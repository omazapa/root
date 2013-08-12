// From gain calibration of the DESY CMS Pixel Upgrade group
// Copyright (C) 2013 A. Burgmeier, D. Pitzl (DESY)

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TClass.h"
#include "TKey.h"
#include "TObject.h"
#include "TMath.h"
#include "TSystem.h"
#include <Fit/FitConfig.h>
#include <TFitResult.h>

#include <string>
#include <iostream>

Double_t dWeibdp( Double_t *x, Double_t *par, Double_t *grd )
{
  // Weibull derivative w.r.t. par
  double t = par[0] + 1E-5 * x[0] / par[1]; // t near 1, decorrelate p0
  double p = par[2];

  double l = log(t); // small, ln(1) = 0
  double a = exp( p*l ); // t^p

  double c = exp(-a); // not 1-exp(-a), decorrelates p4
  double f = par[4] + par[3]*c;

  if(grd)
  {
    double b = exp( (p-1)*l ); // t^(p-1)

    grd[4] = 1;
    grd[3] = c;
    grd[2] =-par[3] * c * a * l; // dWeib/dp2
    grd[1] = par[3] * c * p * b * 1E-5 * x[0]/par[1]/par[1]; // dWeib/dp1
    grd[0] =-par[3] * c * p * b; // dWeib/dp0
  }

  return f; // Weib
}

Double_t Weib( Double_t *x, Double_t *par)
{
  return dWeibdp(x, par, NULL);
}

int main()
{
  TFile* f = new TFile("gainmap-Ia34-trim26.root", "READ");
  if(!f || f->IsZombie())
  {
    std::cerr << "Could not read input file gainmap-Ia34-trim26.root" << std::endl;
    return 1;
  }

  // Set minimizer:
  ROOT::Fit::FitConfig::SetDefaultMinimizer("LVMini", "");
  const double err = sqrt(1.0/0.3);

  // loop over histos:
  TIter nextkey( gDirectory->GetListOfKeys() );

  int nh = 0;

  while( TKey* key = (TKey*)nextkey() ) {

    TObject * obj = key->ReadObj();

    if( obj->IsA()->InheritsFrom("TH1") ) {

      TH1 *h = (TH1*)obj;

      std::string hn = h->GetName(); // PHvsVcal_c05r07_C0
      if( hn.substr( 0, 8 ) != "PHvsVcal" ) continue;

      int hl = hn.length();
      if( hl != 18 ) continue;

      std::string hc = hn.substr( 10, 2 ); // col
      int icol = atoi( hc.c_str() );

      std::string hr = hn.substr( 13, 2 ); // row
      int irow = atoi( hr.c_str() );

      std::cout << std::endl << nh << ": col" << icol << ", row " << irow << std::endl;

      int imax = h->GetMaximumBin();
      double amax = h->GetBinContent(imax);

      if( amax > 9 ) {
	// set bin error:

	int nb = h->GetNbinsX();
	for( int ii = 1; ii <= nb; ++ii ){
	  h->SetBinError( ii, err ); // large Vcal
	}

	// fit range:

	int ib1 = 1;// zero is underflow
	int ib9 = nb;// last bin

	for( int ii = ib9; ii > 0; --ii ) { // scan from right to left

	  if( h->GetBinContent(ii) > 0 )
	    ib1 = ii; // overwritten
	  else
	    break; // first zero bin from right
	}

        h->GetXaxis()->SetRange(ib1, ib9); // for restricting fit range
        double x1 = h->GetBinLowEdge(ib1);
        double x9 = h->GetBinLowEdge(ib9)+h->GetBinWidth(ib9);

        // set start values:
	TF1 *fWeib = new TF1( "Weib", Weib, x1, x9, 5 );
        fWeib->SetGradientFunction(dWeibdp);
        fWeib->SetParName(0, "horizontal offset");
        fWeib->SetParameter(0, 1.);
        fWeib->SetParName(1, "width");
        fWeib->SetParameter(1, 1.5); // x scaled
        fWeib->SetParName(2, "power");
        fWeib->SetParameter(2, 2794.);
        fWeib->SetParameter(3, -amax-22.);
        fWeib->SetParName(4, "offset");
        fWeib->SetParameter(4, amax);
        fWeib->CheckGradientFunction(1e-6); // decrease epsilon for better precision

        // fit
        TFitResultPtr result = h->Fit(fWeib, "SQG");
        for(int i = 0; i < 5; ++i)
          std::cout << "Parameter " << i << ": " << result->Parameter(i) << "  +/-  " << result->ParError(i) << std::endl;
      } // amax > 9

    }//TH1
  }//while

  std::cout << "These were " << nh << " histos\n" << std::endl;
  return 0;
}
