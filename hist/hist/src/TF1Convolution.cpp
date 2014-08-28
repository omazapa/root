//
//  TF1Convolution.cpp
//  
//
//  Created by Aurélie Flandi on 27.08.14.
//
//

#include "TF1Convolution.h"


// class wrapping evaluation of TF1(t) * TF1(x-t)
class TF1Convolution_EvalWrapper{
   const TF1* fFunction1;
   const TF1* fFunction2;
   const double fT0;
public:
   TF1Convolution_EvalWrapper(const TF1* function1 , const TF1* function2, double t)
   :fFunction1(function1), fFunction2(function2), fT0(t){}
   double operator()(double x) const
   {
      return fFunction1->Eval(x)*fFunction2->Eval(x-fT0);
   }
   
};
/*
TF1Convolution::TF1Convolution(TF1* function1, TF1* function2, Double_t xmin, Double_t xmax)
{

   

}
Double_t TF1Convolution::MakeConvolution(Double_t x)
{
 
 ROOT::Math::IntegratorOneDim func(TF1Convolution_EvalWrapper(function1, function2, x), ROOT::Math::IntegratorOneDimOptions::DefaultIntegratorType(), epsabs, epsrel);
 return func.Integral(xmin, xmax);
 
}
*/