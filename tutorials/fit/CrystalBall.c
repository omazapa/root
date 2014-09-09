//
//  CrystalBall.c
//  
//
//  Created by Aurélie Flandi on 09.09.14.
//
//

#include <stdio.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TChain.h>
#include <TObject.h>
#include <TRandom.h>
#include <TFile.h>
#include <math.h>
#include <TF1NormSum.h>
#include <TF1.h>
#include <TH1F.h>
#include <TGraph.h>
#include <Math/PdfFuncMathCore.h>
#include <Math/IntegratorOptions.h>

using namespace std;


void CrystalBall()
{

   // tutorial for normalized sum of two functions
   // here: a background exponential and a crystalball function
   
   
   // parameters can be set:
   // I. with the TF1 object before adding the function
   // II. with the TF1NormSum object (first two are the coefficients, then the parameters without the constant coefficients
   // III. with the TF1 object after adding the function
   
   // sum can be constructed by:
   // 1) by a string containing the names of the functions and/or the coefficient in front
   // 2) by a string containg formulas like expo, gaus...
   // 3) by the list of functions and coefficients (are 1 by default)
   // 4) by a std::vector for functions and coefficients
   
   
   Int_t NEvents = 1e5;
   Int_t NBins   = 1e2;
   Double_t Epsilon = 0.1;
   
   ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("GAUSS");
   ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1e-6);
   ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-6);

   
   TF1 *f_cb  = new TF1("MyCrystalBall","ROOT::Math::crystalball_pdf(x,[0],[1],[2],[3])",-5.,5.);
   TF1 *f_exp = new TF1("MyExponential","expo",-5.,5.);
   TF1 *f_gauss  = new TF1("MyGaussian",   "gaus",-5.,5.);
   // I.:
   //f_exp-> SetParameters(0.,-0.3);
   //f_cb -> SetParameters(1,2,3,0.3);
   // 1) :
   TF1NormSum *fnorm_exp_cb = new TF1NormSum("0.2*MyExponential + MyCrystalBall");
   // 2) :
   //TF1NormSum *fnorm_exp_cb = new TF1NormSum("0.2*expo + MyCrystalBall");
   // 3) :
   //TF1NormSum *fnorm_exp_cb = new TF1NormSum(f_exp,f_cb,0.2,1.);
   // 4) :
   //const std::vector < TF1*     > functions  = {f_exp, f_cb};
   //const std::vector < Double_t > coeffs     = {0.2,1};
   //TF1NormSum *fnorm_exp_cb = new TF1NormSum(functions,coeffs);

   // II. :
   //fnorm_exp_cb -> SetParameters(1.,1.,-0.3,1,2,3,0.3);
   TF1   * f_sum = new TF1("fsum",*fnorm_exp_cb, -5., 5., fnorm_exp_cb->GetNpar());
   // III.:
   f_sum -> SetParameters(3.,1.,-0.5,1.,2.,3,0.3);
  
   //histogram to fit
   TH1F *h_sum = new TH1F("h_ExpCB","Exponential Bkg + CrystalBall function",NBins,-5.,5.);
   h_sum -> FillRandom("fsum", NEvents);
   h_sum -> Sumw2();
   h_sum -> Scale(1.,"width");

   //fit
   f_sum -> SetParameters(3.+Epsilon,1.-Epsilon,-0.5+Epsilon,1.+Epsilon, 2.+Epsilon,3+Epsilon,0.3+Epsilon);
   new TCanvas("crystal ball","crystal ball",800,1000);
   //h_sum -> Fit("fsum");
   h_sum -> Draw();
   
   
   
}
