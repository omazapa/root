// @(#)root/mathcore:$Id$
// Authors: L. Moneta, J.T. Offermann    08/2013 
/**********************************************************************
 *                                                                    *
 * Copyright (c) 2013 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/
/*
 * NumGradFunction.h
 *
 *  Created on: Aug 15, 2013
 *      Author: L. Moneta, J.T. Offermann
 */

#ifndef ROOT_Math_NumGradFunction
#define ROOT_Math_NumGradFunction

#ifndef ROOT_Math_IFunctionfwd
#include <Math/IFunctionfwd.h>
#endif

#include <vector>

namespace ROOT {
namespace Math {

class NumericalDerivator; 

class NumGradFunction : public ROOT::Math::IGradientFunctionMultiDim {
    
public:
    
    NumGradFunction(const ROOT::Math::IBaseFunctionMultiDim &f);
    virtual ~NumGradFunction();
    ROOT::Math::IMultiGenFunction * Clone() const;
    
    unsigned int NDim() const;
    
    void Gradient (const double * x, double * g) const;
    
    void Gradient (const double * x, double * g, double * g2) const;
    
    void setfGradientCounter (int val)
    {
        this->fGradientCounter = val;
    }
    
    virtual void FdF (const double * x, double &f, double * df) const;
    
    void SetStrategy (int i);
    
    void SetInitialGradient( const double * s);
    
private:
    
    double DoEval( const double * x) const;
    
    double DoDerivative (const double * x, unsigned int i) const;
//    {
//        std::vector<double> g(fDim);
//        Gradient(x,&g[0]);
//        return  g[i];
//    }
    
    
    const ROOT::Math::IBaseFunctionMultiDim* fFunction;
    NumericalDerivator* fDerivator;
    mutable int fGradientCounter;
    mutable int fDoEvalCounter;
};


} // namespace Math
} // namespace ROOT


#endif /* NumGradFunction_H_ */
