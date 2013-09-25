// @(#)root/lvmini:$Id$
// Author: A. Burgmeier

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2013  DESY                                           *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Implementation file for class LVMiniMinimizer

#include "LVMini/LVMiniMinimizer.h"

#include "Math/IFunction.h"
#include "Math/IOptions.h"
#include "Math/NumGradFunction.h"
#include "Math/MinimTransformFunction.h"
#include <cassert>
#include <iostream>

// Declare the fortran symbols
extern "C" {
   extern void lvmeps_(float* eps, float* wlf1, float* wlf2);
   extern int lvmdim_(int* npar, int* mvec);
   extern void lvmini_(int* npar, int* mvec, int* nfcn, double* aux);
   extern void lvmfun_(double* x, double* f, int* iret, double* aux);
   extern int lvmind_(int* iarg);
   extern void lvmprt_(int* lup, double* aux, int* iarg);
}

namespace {
   double getCovMatrixEntry(double aux[], int i, int j) {
      int iarg = 5;
      int ind = lvmind_(&iarg);
      if(ind == 0) return -1.;

      // Matrix is symmetric
      if(j > i) std::swap(i, j);

      const unsigned int ijInd = ind - 1 + (j + 1) + ((i + 1) * (i + 1) - (i + 1)) / 2;
      //std::cout << "Index: " << ijInd << ", entry=" << aux[ijInd] << std::endl;
      return aux[ijInd];
   }
}

namespace ROOT { 

namespace LVMini {

LVMiniMinimizer::LVMiniMinimizer(bool calcerrors, float eps, float wlf1, float wlf2):
   fCalcErrors(calcerrors), fFunc(NULL), fAux(NULL), fMin(0.), fIterations(0)
{
   lvmeps_(&eps, &wlf1, &wlf2);
}

LVMiniMinimizer::LVMiniMinimizer(const char* type):
   fCalcErrors(strcmp(type, "NoErrors") != 0), fFunc(NULL), fAux(NULL), fMin(0.), fIterations(0)
{
}

LVMiniMinimizer::~LVMiniMinimizer()
{
   delete[] fAux;
   if (fFunc)
   {
      delete fFunc;
   }
}

void LVMiniMinimizer::Clear()
{
   fVariables.clear();
   fVariableNames.clear();

   delete[] fAux;
   if (fFunc)
   {
      delete fFunc;
      fFunc = 0;
   }
   fAux = 0;
}

void LVMiniMinimizer::SetFunction(const ROOT::Math::IMultiGenFunction& func)
{
//   std::cerr << "LVMiniMinimizer::SetFunction: LVMini needs gradient, use LVMiniMinimizer::SetFunction(const IMultiGradFunction&) instead!" << std::endl;
//   std::cerr << "If you are using TH1::Fit, use the fitting option \"G\", and, ideally, provide an analytic gradient calculation with TF1::SetGradientFunction" << std::endl;
   ROOT::Math::NumGradFunction* gradf = new ROOT::Math::NumGradFunction(func);
   // need to delete (memory leak)
   if (fFunc)
   {
      delete fFunc;
   }
   fFunc = gradf; //new
}

void LVMiniMinimizer::SetFunction(const ROOT::Math::IMultiGradFunction& func)
{
   if (fFunc)
   {
      delete fFunc;
   }
   fFunc = dynamic_cast< ROOT::Math::IMultiGradFunction*>(func.Clone());
   assert(fFunc);
   //fFunc = const_cast<ROOT::Math::IMultiGradFunction*>(&func);
}

bool LVMiniMinimizer::SetVariable(unsigned int var, const std::string& varname, double start, double step)
{
   if(var != fVariables.size())
   {
      std::cerr << "LVMiniMinimizer::SetVariable: Wrong index! Variables need to be added consecutively, use index=" << fVariables.size() << std::endl;
      return false;
   }
   else
   {
      fVariables.push_back(start);
      fSteps.push_back(step);
      fVariableNames.push_back(varname);
      fVarTypes.push_back(ROOT::Math::kDefault);
      return true;
   }
}

bool LVMiniMinimizer::SetLowerLimitedVariable(unsigned int ivar, const std::string & name, double val, double step, double lower) {
   bool ret = SetVariable(ivar, name, val, step);
   if (!ret) return false;
   fBounds[ivar] = std::make_pair( lower, lower);
   fVarTypes[ivar] = ROOT::Math::kLowBound;
   return true;
}

bool LVMiniMinimizer::SetUpperLimitedVariable(unsigned int ivar, const std::string & name, double val, double step, double upper ) {
   bool ret = SetVariable(ivar, name, val, step);
   if (!ret) return false;
   fBounds[ivar] = std::make_pair( upper, upper);
   fVarTypes[ivar] = ROOT::Math::kUpBound;
   return true;
}

bool LVMiniMinimizer::SetLimitedVariable(unsigned int ivar, const std::string & name, double val, double step, double lower, double upper) {
   bool ret = SetVariable(ivar, name, val, step);
   if (!ret) return false;
   fBounds[ivar] = std::make_pair( lower, upper);
   fVarTypes[ivar] = ROOT::Math::kBounds;
   return true;
}

bool LVMiniMinimizer::SetFixedVariable(unsigned int ivar , const std::string & name , double val ) {
   bool ret = SetVariable(ivar, name, val, 0.);
   if (!ret) return false;
   fVarTypes[ivar] = ROOT::Math::kFix;
   return true;
}

bool LVMiniMinimizer::Minimize()
{
   if(fVariables.empty())
   {
      std::cerr << "LVMiniMinimizer::Minimize: No free variables!" << std::endl;
      return false;
   }

   if(!fFunc)
   {
      std::cerr << "LVMiniMinimizer::Minimize: No function to minimize provided!" << std::endl;
      return false;
   }
   
   // check if using the numerical gradient calculation using the 
   // NumGradFunction
   ROOT::Math::NumGradFunction* gradf = dynamic_cast< ROOT::Math::NumGradFunction*>(fFunc);
   if (gradf) {
      // assert(fVariables.size() == gradf->NDim() );
      // std::cout<< "f = " << (*gradf)(&fVariables[0]) << std::endl;
      gradf->SetStrategy(Strategy());
      gradf->SetInitialGradient( &fSteps[0] );
   }

   // Transformation of input variables
   ROOT::Math::MinimTransformFunction* trFunc = 0;
   ROOT::Math::IMultiGradFunction* func = fFunc;
   bool doTransform = false;
   for(unsigned int ivar = 0; !doTransform && ivar < fVarTypes.size(); ++ivar)
      doTransform = (fVarTypes[ivar] != ROOT::Math::kDefault);

   int npar = fVariables.size();
   std::vector<double> x(fVariables);
   if(doTransform)
   {
      trFunc = new ROOT::Math::MinimTransformFunction(fFunc, fVarTypes, fVariables, fBounds);
      trFunc->InvTransformation(&fVariables[0], &x[0]);

      npar = trFunc->NDim();
      x.resize(npar);

      func = trFunc;
      //N.B. the Minim transform function manages the given function pointer 
      fFunc = trFunc; 
   }

   // I haven't seen a clear guide on how many vector pairs
   // should be used. This is a rough estimate from table 1 of
   // Blobel's manual and the recommendation to use values between
   // 6 and 29.
   int mvec = std::min(std::max(6, npar / 5), 29);
   if(fCalcErrors) mvec = -mvec;
   int npar_sign = npar;
   if(fDebug) npar_sign = -npar;

   int mdim = lvmdim_(&npar, &mvec);
    
#if 0
    float tol = Tolerance(); // edited by L. Moneta
    float wlf1 = 0.001;
    float wlf2 = 0.9;
    lvmeps_(&tol, &wlf1, &wlf2); //end of edits
#endif

   delete[] fAux;
   fAux = new double[mdim];

   int nfcn = fMaxCalls;

   lvmini_(&npar_sign, &mvec, &nfcn, fAux);

   unsigned int fMaxIterations = fMaxIter;
   if(fMaxIterations == 0) fMaxIterations = 10000;

   int iret = -1;
   fIterations = 0;

   // check if the second derivatives are computed
   std::vector<double> g2(npar);
   do
   {
      func->FdF(&x[0], fMin, &fAux[0], &g2[0]);
      if (g2[0] != std::numeric_limits<double>::quiet_NaN() )
         std::copy(g2.begin(), g2.end(), &fAux[npar]);
      lvmfun_(&x[0], &fMin, &iret, &fAux[0]);
   } while(++fIterations < fMaxIterations && iret < 0);

   fStatus = iret;

   // Remember minimum location
   const double* xout = &x[0];
   if(doTransform)
      xout = trFunc->Transformation(&x[0]);
   std::copy(&xout[0], &xout[npar], fVariables.begin());

   if(iret == -1)
   {
      // Max number of iterations reached
      fValidError = false;
      std::cerr << "LVMiniMinimizer::Minimize: Max iterations reached" << std::endl;
      return false;
   }
   else if(iret > 0)
   {
      fValidError = false;
      std::cerr << "LVMiniMinimizer::Minimize: Minimization failed, error code=" << iret << std::endl;
      return false;
   }
   else
   {
      // Read and transform the errors and covariance matrix.
      // First, read (internal) covariance matrix into an array
      double* covMatrix;
      if(doTransform)
      {
         covMatrix = new double[npar * npar];
      }
      else
      {
         fCovMatrix.resize(npar * npar);
         covMatrix = &fCovMatrix[0];
      }

      if(fCalcErrors)
      {
         fValidError = true;
         // TODO: Could be optimized, by only obtaining the index once
         // and then iterating.
         for(unsigned int i = 0; i < static_cast<unsigned int>(npar); ++i)
         {
            for(unsigned int j = 0; j < i; ++j)
            {
               covMatrix[i * npar + j] = covMatrix[j * npar + i] = getCovMatrixEntry(fAux, i, j) * 2 * fUp;
            }

            covMatrix[i * npar + i] = getCovMatrixEntry(fAux, i, i) * 2 * fUp;
         }

#if 0
         int iarg = 3;
         int ind = lvmind_(&iarg);
         for(unsigned int i = 0; i < static_cast<unsigned int>(npar); ++i)
            std::cout << "Error " << i << ": " << fAux[ind + i] << std::endl;
#endif
      }
      else
      {
         fValidError = false;

         // Read approximate errors, and fill a diagonal covariance matrix
         int iarg = 2;
         int ind = lvmind_(&iarg);
         if(ind != 0)
         {
            for(unsigned int i = 0; i < static_cast<unsigned int>(npar); ++i)
            {
               for(unsigned int j = 0; j < i; ++j)
                  covMatrix[i * npar + j] = covMatrix[j + npar + i] = 0.;
               covMatrix[i * npar + i] = fAux[ind + i] * 2 * fUp;
            }
         }
      }

      // Next, convert internal to external covariance matrix
      if(doTransform)
      {
         fCovMatrix.resize(fVariables.size() * fVariables.size());
         trFunc->MatrixTransformation(&x[0], covMatrix, &fCovMatrix[0]);
         delete[] covMatrix;
      }

      // Finally, fill external errors from diagonal elements of covariance matrix
      fErrors.resize(fVariables.size());
      for(unsigned int i = 0; i < fVariables.size(); ++i)
         fErrors[i] = sqrt(fCovMatrix[i * fVariables.size() + i]);

      // Minimization complete
      return true;
   }
}

double LVMiniMinimizer::MinValue() const
{
   int iarg = 0;
   int ind = lvmind_(&iarg);
   return fAux[ind - 1];
}

double LVMiniMinimizer::Edm() const
{
   // Return magnitude of gradient vector at minimum
   const double* grad = MinGradient();
   double gradSqr = 0.;
   for(unsigned int i = 0; i < fVariables.size(); ++i)
      if(fVarTypes[i] != ROOT::Math::kFix)
         gradSqr += grad[i]*grad[i]; // TODO: *2.0*fUp?
   return sqrt(gradSqr);
}

const double* LVMiniMinimizer::X() const
{
   return &fVariables[0];
}

const double* LVMiniMinimizer::MinGradient() const
{
   fMinGradient.resize(fVariables.size());
   fFunc->Gradient(X(), &fMinGradient[0]);
   return &fMinGradient[0];
}

unsigned int LVMiniMinimizer::NCalls() const
{
   int iarg = -2;
   return lvmind_(&iarg);
}

unsigned int LVMiniMinimizer::NDim() const
{
   return fVariables.size();
}

unsigned int LVMiniMinimizer::NFree() const
{
   unsigned int nFree = 0;
   for(unsigned int i = 0; i < fVarTypes.size(); ++i)
      if(fVarTypes[i] != ROOT::Math::kFix)
         ++nFree;
   return nFree;
}

bool LVMiniMinimizer::ProvidesError() const
{
   return fCalcErrors;
}

const double* LVMiniMinimizer::Errors() const
{
   return &fErrors[0];
}

double LVMiniMinimizer::CovMatrix(unsigned int i, unsigned int j) const
{
   return fCovMatrix[i * fVariables.size() + j];
}

bool LVMiniMinimizer::GetCovMatrix(double* cov) const
{
   if(!fCalcErrors)
      return false;

   std::copy(fCovMatrix.begin(), fCovMatrix.end(), cov);
   return true;
}

int LVMiniMinimizer::CovMatrixStatus() const
{
   int iarg = 5;
   int ind = lvmind_(&iarg);
   if(ind == 0) return -1;
   return 3;
}

double LVMiniMinimizer::GlobalCC(unsigned int i) const
{
   // TODO: Recalculate this from the covariance matrix
/*
   int iarg = 4;
   int ind = lvmind_(&iarg);
   if(ind == 0) return -1.;

   return fAux[ind + i];*/
   return -1.;
}

void LVMiniMinimizer::PrintResults()
{
   int lup = 0;
   int iarg = 6;
   lvmprt_(&lup, fAux, &iarg);
}

} // end namespace LVMini

} // end namespace ROOT

