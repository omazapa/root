#include <stdio.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <iostream>//for cout
#include <iomanip> //for setprecision
#include <string>
#include <sstream>
#include <TROOT.h>
#include <TChain.h>
#include <TObject.h>
#include <TFile.h>
#include <math.h>
#include <TF1NormSum.h>
#include <TF1.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TStopwatch.h>
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "Math/MinimizerOptions.h"

using namespace RooFit ;
using namespace std;


void fitNormSum()
{
   
   //compare fitting of a normalized sum of several functions with ROOT or ROOFIT
   
   // ROOT -------------------------------------------------------------------------------------------------------
  
   //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   
   
   Int_t NEvents = 1e6;
   Int_t NBins   = 1e4;
   
   // PARAMETERS ...........................
   
   TF1 *f_ExpGauss  = new TF1("f_ExpGauss","expo(0)+gaus(2)+gaus(5)",-5.,5.);

   //*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*
   Double_t Coeff0  = 1.;
   Double_t Alpha   = 0.5;
   Double_t Coeff1  = 10.;
   Double_t Mean1   = 2.;
   Double_t Sigma1  = 0.5;
   Double_t Coeff2  = 10.;
   Double_t Mean2   = 0.;
   Double_t Sigma2  = 0.3;
   
   double LowerLimits[]     = {100.,       100.,       100.,       0.,        -10.,      0.,         -10.,      0.};
   double UpperLimits[]     = {1.e9,       1.e9,       1.e9,       10.,       10.,       100.,       10.,       100.};
   double ParametersArray[] = {1e5,        1.e5,      1.e5, Alpha, Mean1, Sigma1, Mean2, Sigma2};
   string StringArray[]     = {"nbkg",     "nsig1",    "nsig2",    "alpha",   "mean1",   "sigma1",   "mean2",   "sigma2"};
   //*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*
    
   f_ExpGauss -> SetParameters(Coeff0,Alpha,Coeff1,Mean1,Sigma1,Coeff2,Mean2,Sigma2);//exact function that we want to fit
   // HISTOGRAM TO FIT ....................

   TH1F *h_ExpGauss = new TH1F("h_ExpGauss","Exponential Bkg + 2 gaussian peaks signal",NBins,-5.,5.);
   h_ExpGauss -> FillRandom("f_ExpGauss", NEvents);
   TH1F *h_copie    = new TH1F(*h_ExpGauss);
   h_ExpGauss -> Sumw2();
   h_ExpGauss -> Scale(1.,"width");
   // h_copie = h_ExpGauss;
   //CONSTRUCTION OF THE TF1NORMSUM ........
    
   TF1 *f_exp    = new TF1("Exponential","expo",-5.,5.);
   TF1 *f_gauss  = new TF1("Gaussian",   "gaus",-5.,5.);
   TF1 *f_gauss2 = new TF1("Gaussian2",  "gaus",-5.,5.);

   TF1NormSum *fnorm_exp_gauss = new TF1NormSum("expo + gaus + gaus");
   Int_t NofParams = fnorm_exp_gauss -> GetNpar();

   TF1   * fsum = new TF1("fsum",*fnorm_exp_gauss, -5., 5., NofParams);

   fsum -> SetParameters(1e5, 1.e5, 1.e5, Alpha, Mean1, Sigma1, Mean2, Sigma2);
   for (int i=0; i<NofParams; i++)
   {
      fsum -> SetParLimits(i, LowerLimits[i],UpperLimits[i]);
   }
   
   // FIT ......................................
    
   new TCanvas("c_gaussexp","c_gaussexp",800,1000);
   TStopwatch tw;
   tw.Start();
   for (int i=0;i<10;i++)
   {
     fsum -> SetParameters(1e5, 1.e5, 1.e5, Alpha, Mean1, Sigma1, Mean2, Sigma2);
      h_ExpGauss -> Fit("fsum","L");
   }
   cout << "**************************************************" << endl;
   tw.Print(); //time to fit
   cout << "**************************************************" << endl;
   h_ExpGauss -> Draw();

   cout << "Integral: " << fsum -> Integral(-5.,5.) << endl;
    
   
   // --------------------------------------------------------------------------------------------------------------
   // ROOFIT -------------------------------------------------------------------------------------------------------
   // --------------------------------------------------------------------------------------------------------------
    
   cout << endl;
   cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
   cout << "%%%%%%%%%%%%%%%%%%%%%%%%%  ROOFIT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
   cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
   cout << endl;
    
   //creates workspace and pdf ....
    
   RooWorkspace w("w",kTRUE);//kTRUE creates a C++ namespace
   w.factory("Exponential:: Expo  (x[-5.,5.], alpha[0.6,0.,10.])");
   w.factory("Gaussian   :: Gauss (x[-5.,5.], mean1[2.1, -10.,10.], sigma1[0.4,0., 100.])");
   w.factory("Gaussian   :: Gauss2(x[-5.,5.], mean2[0.3,-10.,10.],  sigma2[0.2,0.,100.])");
   w.factory("SUM:model( nsig1[1.e5,100,1.E9]*Gauss, nsig2[1.e5,100.,1.E9]*Gauss2, nbkg[1.e5,100.,1.E9]*Expo )");  // for extended model
   
   // set paramters and makes data .....
    
   RooRealVar *x      = w.var("x");
   for (int i=0;i<NofParams; i++)
   {
      w.var(StringArray[i].c_str())   -> setVal(ParametersArray[i]);
   }
   RooDataHist data("data","mydata",*x, h_copie);
    
   //fit ...........
    
   TStopwatch tw2;
   tw2.Start();
    for (int i=0;i<10;i++)
    {
     for (int i=0;i<NofParams; i++)
      {
      w.var(StringArray[i].c_str())   -> setVal(ParametersArray[i]);
      }
        w.pdf("model")->fitTo(data);
    }
   cout << "**************************************************" << endl;
   tw2.Print(); // time to fit
   cout << "**************************************************" << endl;
   
   //draw .........
    
   RooPlot* frame = x->frame();
   data.plotOn(frame);
   w.pdf("model")->plotOn(frame);
   new TCanvas("Roofit gauss","Roofit gauss",800,1000);
   frame->Draw();
}


