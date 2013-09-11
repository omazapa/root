// @(#)root/mathcore:$Id$
// Authors: L. Moneta, J.T. Offermann    08/2013 
/**********************************************************************
 *                                                                    *
 * Copyright (c) 2013 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/
/*
 * NumericalDerivator.h
 *
 *  Created on: Aug 14, 2013
 *      Authors: L. Moneta, J. T. Offermann
 */

#ifndef ROOT_Math_NumericalDerivator
#define ROOT_Math_NumericalDerivator

#ifndef ROOT_Math_IFunctionfwd
#include <Math/IFunctionfwd.h>
#endif

#include <vector>


namespace ROOT {
namespace Math {


class NumericalDerivator {
public:

   NumericalDerivator();
   NumericalDerivator(const ROOT::Math::IBaseFunctionMultiDim &f, double step_tolerance, double grad_tolerance, unsigned int ncycles);
   virtual ~NumericalDerivator();
   const double* Differentiate(const double* x);
   double GetFValue() const {
      return fVal;
   }
   const double * GetG2() {
      return &fG2[0];
   }
   void SetStepTolerance(double value);
   void SetGradTolerance(double value);
   void SetNCycles(int value);
   
   void SetInitialValues(const double* g, const double* g2, const double* s);
       
   void SetInitialGradient (const double * s);

private:

    std::vector <double> fGrd;
    std::vector <double> fG2;
    std::vector <double> fGstep;
    const ROOT::Math::IBaseFunctionMultiDim* fFunction;
    double fStepTolerance;
    double fGradTolerance;
    double fNCycles;
    double fVal;
    unsigned int fN;

};


} // namespace Math
} // namespace ROOT

#endif /* NumericalDerivator_H_ */
