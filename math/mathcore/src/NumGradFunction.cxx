// @(#)root/mathcore:$Id$
// Authors: L. Moneta, J.T. Offermann    08/2013 
/**********************************************************************
 *                                                                    *
 * Copyright (c) 2013 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/
// implementation file for NumGradFunction

#include "Math/IFunction.h"
#include "Math/NumGradFunction.h"
#include "Math/NumericalDerivator.h"
#include <cmath>
#include <algorithm>
#include <Math/IFunction.h>
#include <TMath.h>
#include <iostream>

namespace ROOT {
namespace Math {


NumGradFunction::NumGradFunction(const ROOT::Math::IBaseFunctionMultiDim &f):
fFunction(&f),
fGradientCounter(0),
fDoEvalCounter(0)
{
    fDerivator = new NumericalDerivator(f, 0.003, 0.0005, 3);
}

NumGradFunction::~NumGradFunction() { 
    
}

void NumGradFunction::Gradient(const double * x, double * g) const // const double * x
{
    const double* deriv = fDerivator->Differentiate(x);
    std::copy (deriv, deriv+fFunction->NDim(), g);
    fGradientCounter++;
    
//    for (unsigned int i = 0; i<fFunction->NDim(); i++) {
//        std::cout << "g[" << i<< "] = " << g[i] << std::endl;
//    }
    
}

void NumGradFunction::Gradient (const double * x, double * g, double * g2) const {
    Gradient(x, g);
    const double * tmp = fDerivator->GetG2();
    std::copy (tmp, tmp+NDim(), g2);
}

double NumGradFunction::DoEval (const double * x) const {
    fDoEvalCounter ++;
  return (*fFunction)(x);
}

double NumGradFunction::DoDerivative (const double * x, unsigned int i) const {
    const double * derivs = fDerivator->Differentiate(x);
    return derivs[i-1];
}

unsigned int NumGradFunction::NDim() const {
    return fFunction->NDim();
}

ROOT::Math::IMultiGenFunction * NumGradFunction::Clone() const {
    return new NumGradFunction(*fFunction);
}

void NumGradFunction::FdF (const double * x, double &f, double * df) const {
    const double* deriv = fDerivator->Differentiate(x);
    std::copy (deriv, deriv+fFunction->NDim(), df);
    f = fDerivator->GetFValue();
    fGradientCounter ++;
//    for (unsigned int i = 0; i<fFunction->NDim(); i++) {
//        std::cout << "df[" << i << "] = " << df[i] << std::endl;
//    }
}

void NumGradFunction::SetStrategy (int i) {
    if (i == 0) {
        fDerivator->SetStepTolerance(0.5);
        fDerivator->SetGradTolerance(0.1);
        fDerivator->SetNCycles(2);
    }
    if (i == 1) {
        fDerivator->SetStepTolerance(0.3);
        fDerivator->SetGradTolerance(0.05);
        fDerivator->SetNCycles(3);
    }
    else
    {
        fDerivator->SetStepTolerance(0.1);
        fDerivator->SetGradTolerance(0.02);
        fDerivator->SetNCycles(5);
    }
}

void NumGradFunction::SetInitialGradient( const double * s) {
    fDerivator->SetInitialGradient(s);
}

} // namespace Math
} // namespace ROOT

