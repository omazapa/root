//
//  TestConv.cpp
//  
//
//  Created by Aurélie Flandi on 28.08.14.
//
//

#include <stdio.h>
#include <TMath.h>
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
#include <TRandom.h>
#include <TFile.h>
#include <math.h>
#include <TF1Convolution.h>
#include <TF1.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStopwatch.h>
#include "Math/MinimizerOptions.h"
#include <Math/PdfFuncMathCore.h>

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"


using namespace RooFit ;
using namespace std;


void TestConv()
{
   
   
   TF1* f_try = new TF1("MyGaussian", "gaus(0)",-10.,10.);
   f_try -> SetParameters(1.,0.,1.);

   TH1F *h_Gauss = new TH1F("h_Gauss","gaussian",100,-10,10.);
   h_Gauss->FillRandom("MyGaussian",1e5);
   new TCanvas("try","try",800,1000);
   h_Gauss -> Draw();
   
   
   return;
  // gErrorIgnoreLevel=3001;
   TH1F *h_ExpGauss = new TH1F("h_ExpGauss","Exponential convoluted by gaussian",100,0.,5.);
   
   for (int i=0;i<1e6;i++)
   {
      Double_t x = gRandom->Exp(1./0.3);//gives a alpha of -0.3 in the
      x += gRandom->Gaus(0.,1.);
      h_ExpGauss->Fill(x);
   }
   new TCanvas("h","h",800,1000);
   h_ExpGauss->Draw();
   gPad->Update();//return;
   
   
   //convolution
   //TF1 *f_exp    = new TF1("Exponential","expo",0.,10.);
   //TF1 *f_gauss  = new TF1("Gaussian",   "gaus",-10.,10.);
   TF1Convolution *f_conv = new TF1Convolution("expo","gaus");
   f_conv->SetRange(-1.,6.);
   f_conv->SetNofPointsFFT(1000);
  // f_conv->SetNumConv(true);
   TF1   *f = new TF1("f",*f_conv, 0., 5., f_conv->GetNpar());
   f->SetParameters(1.,-0.3,0.,1.);
   new TCanvas("f","f",800,1000);
   f->Draw();
  

  
   //histogram to fit
  // TH1F *h_ExpGauss = new TH1F("h_ExpGauss","Exponential convoluted by gaussian",100,0.,5.);
  // h_ExpGauss -> FillRandom("f", 1e6);
  
   
   //fit
   new TCanvas("c","c",800,1000);
   TStopwatch tw1;
   tw1.Start();
   for (int i=0;i<100;i++)
   {
       f->SetParameters(1.,-0.3,0.,1.);
      h_ExpGauss -> Fit("f");
   }
   cout << "**************************************************" << endl;
   tw1.Print();
   cout << "**************************************************" << endl;
   h_ExpGauss->Draw();
   //f->Draw();
   //
   TF1 *f_cb = new TF1("crystal ball","ROOT::Math::crystalball_pdf(x,[0],[1],[2],[3])",-5.,5.);
   f_cb ->SetParameters(1,2,2.5,0.8);
   new TCanvas("crystal ball","crystal ball",800,1000);
   f_cb -> Draw();
   
   
   
   // --------------------------------------------------------------------------------------------------------------
   // ROOFIT -------------------------------------------------------------------------------------------------------
   // --------------------------------------------------------------------------------------------------------------
   
   cout << endl;
   cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
   cout << "%%%%%%%%%%%%%%%%%%%%%%%%%  ROOFIT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
   cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
   cout << endl;
   
   RooWorkspace w("w",kTRUE);//kTRUE creates a C++ namespace
   w.factory("Exponential::MyExpo (x[0.,5.], alpha[-0.3,-5.,0.])");
   w.factory("Gaussian::MyGauss(x, mean1[0, -5.,5.], sigma1[1.0,0., 10.])");
    RooRealVar *x = w.var("x");
   x->setBins(1e4,"cache");
   w.factory("FCONV::N(x,MyExpo,MyGauss)");
  
   RooFFTConvPdf* pdf = (RooFFTConvPdf*)w.pdf("N");
   pdf->setBufferFraction(0.5);
   

  /* RooRealVar t("t","t",0.,5.);

   RooAbsPdf *r_expo =  w.pdf("MyExpo");
   RooAbsPdf *r_gaus =  w.pdf("MyGauss");
   RooNumConvPdf("Exp conv.by gaus", "f_conv", *x, *r_expo, *r_gaus);*/
   
   RooDataHist data("data","mydata",*x, h_ExpGauss);
   TStopwatch tw2;
   tw2.Start();
   w.pdf("N")->fitTo(data);
   cout << "**************************************************" << endl;
   tw2.Print();
   cout << "**************************************************" << endl;
   //draw .........
   
   RooPlot* frame = x->frame();
   data.plotOn(frame);
   w.pdf("N")->plotOn(frame);
   new TCanvas("Roofit conv","Roofit conv",800,1000);
   frame->Draw();return;

   
   
   
}
