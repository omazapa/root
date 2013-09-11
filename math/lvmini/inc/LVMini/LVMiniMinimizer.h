// @(#)root/lvmini:$Id$
// Author: A. Burgmeier

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2013  DESY                                           *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class Minuit2Minimizer

#ifndef ROOT_LVMini_LVMiniMinimizer
#define ROOT_LVMini_LVMiniMinimizer

#ifndef ROOT_Math_Minimizer
#include "Math/Minimizer.h"
#endif

#ifndef ROOT_Math_IFunctionfwd
#include "Math/IFunctionfwd.h"
#endif

#ifndef ROOT_Math_IParamFunctionfwd
#include "Math/IParamFunctionfwd.h"
#endif

namespace ROOT { 

   namespace LVMini { 

//______________________________________________________________________________
/** 
   LVMiniMinimizer class implementing the ROOT::Math::Minimizer interface for
   the LVMini algorithm.
   In ROOT it can be instantiated using the plug-in manager (plug-in "LVMini")
*/
class LVMiniMinimizer : public ROOT::Math::Minimizer {

public: 

   /** 
      Default constructor
   */ 
   LVMiniMinimizer (bool calcerrors = true, float eps = 1e-4, float wlf1 = 1e-4, float wlf2 = 0.9);

   /** 
      Constructor with a char (used by PM) 
   */ 
   LVMiniMinimizer (const char *  type); 

   /** 
      Destructor (no operations)
   */ 
   virtual ~LVMiniMinimizer ();

private:
   // usually copying is non trivial, so we make this unaccessible

   /** 
      Copy constructor
   */ 
   LVMiniMinimizer(const LVMiniMinimizer &); 

   /** 
      Assignment operator
   */ 
   LVMiniMinimizer & operator = (const LVMiniMinimizer & rhs); 

public: 
   /* Implementation of the ROOT::Math::Minimizer interface */
   virtual void Clear();
   virtual void SetFunction(const ROOT::Math::IMultiGenFunction& func);
   virtual void SetFunction(const ROOT::Math::IMultiGradFunction& func);
   virtual bool SetVariable(unsigned int var, const std::string& varname, double start, double step);
   virtual bool Minimize();
   virtual double MinValue() const;
   virtual double Edm() const;
   virtual const double* X() const;
   virtual const double* MinGradient() const;
   virtual unsigned int NCalls() const;
   virtual unsigned int NDim() const;
   virtual unsigned int NFree() const;
   virtual bool ProvidesError() const;
   virtual const double* Errors() const;
   virtual double CovMatrix(unsigned int i, unsigned int j) const;
   virtual bool GetCovMatrix(double * cov) const;
   virtual int CovMatrixStatus() const;
   virtual double GlobalCC(unsigned int i) const;
   virtual void PrintResults();
protected: 
   
private: 
   void common_init();

   const bool fCalcErrors;

   ROOT::Math::IMultiGradFunction* fFunc;

   std::vector<double> fVariables;
   std::vector<double> fSteps;
   std::vector<std::string> fVariableNames;

   // Filled after minimization
   double* fAux;
   double fMin;
   mutable std::vector<double> fMinGradient;
   unsigned int fIterations;
}; 

   } // end namespace LVMini

} // end namespace ROOT



#endif /* ROOT_LVMini_Minuit2Minimizer */
