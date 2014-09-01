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


//class TF1Convolution;
//class TF1Convolution_EvalWrapper;



class TF1Convolution
{
   
   
   protected:
   
   std::shared_ptr <TF1> fFunction1;
   std::shared_ptr <TF1> fFunction2;
   Double_t fXmin;
   Double_t fXmax;
   
   public:
   TF1Convolution(TF1* f, TF1* g);
   TF1Convolution(TF1* f, TF1* g, Double_t xmin, Double_t xmax);
   Double_t MakeConvolution(Double_t x);
   Double_t operator()(Double_t* t, Double_t* p);
   void SetParameters(Double_t* p);
   void SetParameters(Double_t p0, Double_t p1, Double_t p2=0., Double_t p3=0.,
                      Double_t p4=0., Double_t p5=0., Double_t p6=0., Double_t p7=0.);
   
   //ClassDef(TF1Convolution,1)
};


#endif