// @(#)root/mathcore:$Id$
// Authors: L. Moneta, J.T. Offermann    08/2013 
/**********************************************************************
 *                                                                    *
 * Copyright (c) 2013 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/
/*
 * NumericalDerivator.cpp
 *
 *  Created on: Aug 14, 2013
 *      Authors: L. Moneta, J. T. Offermann
 *
 *      This is essentially a slightly modified copy of code
 *      written by M. Winkler, F. James, L. Moneta, and A. Zsenei for Minuit2,
 *      Copyright (c) 2005 LCG ROOT Math team, CERN/PH-SFT.
 */

#include "Math/NumericalDerivator.h"
#include <cmath>
#include <algorithm>
#include <Math/IFunction.h>
#include <iostream>
#include <TMath.h>
#include <cassert>



namespace ROOT {
namespace Math {

NumericalDerivator::NumericalDerivator() :
   fFunction(0), 
   fStepTolerance(0.5), 
   fGradTolerance(0.1),
   fNCycles(2), 
   fVal(0), 
   fN(0)
{}



NumericalDerivator::NumericalDerivator(const ROOT::Math::IBaseFunctionMultiDim &f, double step_tolerance, double grad_tolerance, unsigned int ncycles):
    fFunction(&f),
    fStepTolerance(step_tolerance),
    fGradTolerance(grad_tolerance),
    fNCycles(ncycles)
{
	// constructor with function, and tolerances (coordinates must be specified for differentiate function, not constructor)
//    fStepTolerance=step_tolerance;
//    fGradTolerance=grad_tolerance;
//    fFunction=&f;
    
    fN = f.NDim(); //number of dimensions, will look at vector size
    fGrd.resize(fN);
    for (unsigned int i = 0; i<fN; i++) {
        fGrd[i]=0.1;
    }
	fG2.resize(fN);
    for (unsigned int i = 0; i<fN; i++) {
        fG2[i]=0.1;
    }
	fGstep.resize(fN);
    for (unsigned int i = 0; i<fN; i++) {
        fGstep[i]=0.001;
    }
    fVal = 0;
}

void NumericalDerivator::SetStepTolerance(double value) {
    fStepTolerance = value;
}

void NumericalDerivator::SetGradTolerance(double value) {
    fGradTolerance = value;
}

void NumericalDerivator::SetNCycles(int value) {
    fNCycles = value;
}

NumericalDerivator::~NumericalDerivator() {
	// TODO Auto-generated destructor stub
}

void NumericalDerivator::SetInitialValues(const double* g, const double* g2, const double* s) {
    for (unsigned int i = 0; i<fN; i++) {
        fGrd[i]=g[i];
        fG2[i]=g2[i];
        fGstep[i]=s[i];
    }
}

const double* NumericalDerivator::Differentiate(const double* cx) {
    // std::cout <<"Start:" << std::endl;
    // for (unsigned int i = 0; i<fN; i++) {
    //     std::cout << "fGrd[" << i <<"] = " << fGrd[i] << std::endl;
    //     //std::cout << "fG2[" << i <<"] = " << fG2[i] << std::endl;
    //     std::cout << "fGstep[" << i <<"] = " << fGstep[i] << std::endl;
    //  }
    assert(fFunction != 0);
    std::vector<double> vx(fFunction->NDim());
    assert (vx.size() > 0);
    double *x = &vx[0];
    std::copy (cx, cx+fFunction->NDim(), x);
    double step_tolerance = fStepTolerance;
    double grad_tolerance = fGradTolerance;
    const ROOT::Math::IBaseFunctionMultiDim &f = *fFunction;
    fVal = f(x); //value of function at given points
    double eps = std::numeric_limits<double>::epsilon();
    double eps2 = std::sqrt(eps);//1.e-8; //sqrt(epsilon)
    const double Up = 1;
    double dfmin = double(8.*eps2*(std::abs(fVal)+Up)); //had to cast to double, otherwise "statement has no effect"
    double vrysml = 8.*eps*eps;
    unsigned int ncycle = fNCycles;

    for (int i = 0; i < int(fN); i++) {

       double xtf=x[i]; //current value of coordinate x(i) (looping on i)
       double epspri = eps2 + std::abs(double(fGrd[i]*eps2)); //had to cast to double because I am using std::abs instead of fabs
       double step_old = 0.;
       for (unsigned int j = 0; j < ncycle; ++ j) {
          
          double optstp = std::sqrt(dfmin/(std::abs(fG2[i])+epspri));
          double step = std::max(optstp, std::abs(0.1*fGstep[i])); //had to cast to double again

          double stpmax = 10.*std::abs(fGstep[i]);
          if (step > stpmax) step = stpmax;
          
          double stpmin = std::max(vrysml, 8.*std::abs(eps2*x[i])); //8.*std::abs(double(eps2*x[i]))
          if (step < stpmin) step = stpmin;
          if (std::abs((step-step_old)/step) < step_tolerance) {
             //answer = fGrd[i];
             break;
          }
          fGstep[i] = step;
          step_old = step;
          // std::cout << "step = " << step << std::endl;
          x[i] = xtf + step;
          //std::cout << "x[" << i << "] = " << x[i] <<std::endl;
          double fs1 = f(x);
          //std::cout << "xtf + step = " << x[i] << ", fs1 = " << fs1 << std::endl;
          x[i] = xtf - step;
          double fs2 = f(x);
          //std::cout << "xtf - step = " << x[i] << ", fs2 = " << fs2 << std::endl;
          x[i] = xtf;
          
          
          double fGrd_old = fGrd[i];
          fGrd[i] = 0.5*(fs1-fs2)/step;
//            std::cout << "int i = " << i << std::endl;
//            std::cout << "fs1 = " << fs1 << std::endl;
//            std::cout << "fs2 = " << fs2 << std::endl;
//            std::cout << "fVal = " << fVal << std::endl;
//            std::cout << "step^2 = " << (step*step) << std::endl;
//            std::cout << std::endl;
          fG2[i] = (fs1 + fs2 -2.*fVal)/step/step;
          
          if (std::abs(fGrd_old-fGrd[i])/(std::abs(fGrd[i]+dfmin/step)) < grad_tolerance)
          {
             //answer = fGrd[i];
             break;
          }
       }
    }
    
    // std::cout <<"End:" << std::endl;
    // for (unsigned int i = 0; i<fN; i++) {
    //     std::cout << "fGrd[" << i <<"] = " << fGrd[i] << std::endl;
    //     //std::cout << "fG2[" << i <<"] = " << fG2[i] << std::endl;
    //     std::cout << "fGstep[" << i <<"] = " << fGstep[i] << std::endl;
    // }
    
    return &fGrd[0];
}


void NumericalDerivator::SetInitialGradient( const double * s) {
   // set an initial gradient using some given steps 
   // (used in the first iteration)
 
   for (unsigned int i = 0; i < fN; ++i)  {
        //double eps2 = TMath::Sqrt(fEpsilon);
        //double gsmin = 8.*eps2*(fabs(x[i])) + eps2;
      double dirin = s[i];
      double g2 = 2.0 /(dirin*dirin);
      double gstep = 0.1*dirin;
      double grd = g2*dirin;
      
      fGrd[i] = grd;
      fG2[i] = g2;
      fGstep[i] = gstep;
      
   }
    
}


} // namespace Math
} // namespace ROOT



