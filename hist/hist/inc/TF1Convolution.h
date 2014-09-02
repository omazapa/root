//
//  TF1Convolution.h
//  
//
//  Created by Aurélie Flandi on 27.08.14.
//
//

#ifndef ____TF1Convolution__
#define ____TF1Convolution__

#include <iostream>
#include "TF1.h"
#include <memory>

class TF1Convolution
{
   
   
   protected:
   
   std::shared_ptr <TF1> fFunction1;
   std::shared_ptr <TF1> fFunction2;
   std::vector < Double_t > fParams1;
   std::vector < Double_t > fParams2;
   Double_t fXmin;
   Double_t fXmax;
   Int_t fNofParams1;
   Int_t fNofParams2;
   
   
   public:
   TF1Convolution(TF1* f, TF1* g);
 //  TF1Convolution(TF1* f, TF1* g, Double_t xmin, Double_t xmax);
   Double_t MakeConvolution(Double_t x);
   Double_t operator()(Double_t* t, Double_t* p);
   void     SetParameters(Double_t* p);
   void     SetParameters(Double_t p0, Double_t p1, Double_t p2=0., Double_t p3=0.,
                          Double_t p4=0., Double_t p5=0., Double_t p6=0., Double_t p7=0.);
   
   Int_t    GetNPar() const {return (fNofParams1+fNofParams2);}
   Double_t GetXmin() const {return fXmin;}
   Double_t GetXmax() const {return fXmax;}

   
   //ClassDef(TF1Convolution,1)
};


#endif