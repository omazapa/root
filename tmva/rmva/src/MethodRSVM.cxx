// @(#)root/rmva $Id$
// Author: Omar Zapata 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodRSVM-                                                           *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *       Support Vector Machines                                                  *
 *                                                                                *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 *                                                                                *
 **********************************************************************************/

#include <iomanip>

#include "TMath.h"
#include "Riostream.h"
#include "TMatrix.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include "TMVA/VariableTransformBase.h"
#include "TMVA/MethodRSVM.h"
#include "TMVA/Tools.h"
#include "TMVA/Ranking.h"
#include "TMVA/Types.h"
#include "TMVA/PDF.h"
#include "TMVA/ClassifierFactory.h"

#include "TMVA/Results.h"

using namespace TMVA;

REGISTER_METHOD(RSVM)

ClassImp(MethodRSVM)


//_______________________________________________________________________
MethodRSVM::MethodRSVM(const TString &jobName,
                       const TString &methodTitle,
                       DataSetInfo &dsi,
                       const TString &theOption,
                       TDirectory *theTargetDir) :
   RMethodBase(jobName, Types::kRSVM, methodTitle, dsi, theOption, theTargetDir), fMvaCounter(0)
{
   // standard constructor for the RSVM
   //Booking options
   fScale = kTRUE;
   fType = "C-classification";
   fKernel = "radial";
   fDegree = 3;

   fGamma = (fDfTrain.GetNcols() == 1) ? 1.0 : (1.0 / fDfTrain.GetNcols());
   fCoef0 = 0;
   fCost = 1;
   fNu = 0.5;
   fCacheSize = 40;
   fTolerance = 0.001;
   fEpsilon = 0.1;
   fShrinking = kTRUE;
   fCross = 0;
   fProbability = kTRUE;
   fFitted = kTRUE;

}

//_______________________________________________________________________
MethodRSVM::MethodRSVM(DataSetInfo &theData, const TString &theWeightFile, TDirectory *theTargetDir)
   : RMethodBase(Types::kRSVM, theData, theWeightFile, theTargetDir), fMvaCounter(0)
{
   // standard constructor for the RSVM
   //Booking options
   fScale = kTRUE;
   fType = "C-classification";
   fKernel = "radial";
   fDegree = 3;

   fGamma = (fDfTrain.GetNcols() == 1) ? 1.0 : (1.0 / fDfTrain.GetNcols());
   fCoef0 = 0;
   fCost = 1;
   fNu = 0.5;
   fCacheSize = 40;
   fTolerance = 0.001;
   fEpsilon = 0.1;
   fShrinking = kTRUE;
   fCross = 0;
   fProbability = kTRUE;
   fFitted = kTRUE;

}


//_______________________________________________________________________
MethodRSVM::~MethodRSVM(void)
{
}

//_______________________________________________________________________
Bool_t MethodRSVM::HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets)
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   return kFALSE;
}


//_______________________________________________________________________
void     MethodRSVM::Init()
{
   if (!r.IsInstalled("e1071")) {
      Error("Init", "R's package e1071 is not installed.");
      Log() << kFATAL << " R's package e1071 is not installed."
            << Endl;
      return;
   }

   if (!r.Require("e1071")) {
      Error("Init", "R's package e1071 can not be loaded.");
      Log() << kFATAL << " R's package e1071 can not be loaded."
            << Endl;
      return;
   }
   //Paassing Data to R's environment
   //NOTE:need improved names in R's environment using JobName of TMVA
   r["RMVA.RSVM.fDfTrain"] = fDfTrain;
   
   r["RMVA.RSVM.fWeightTrain"] = fWeightTrain;
   r << "write.table(RMVA.RSVM.fDfTrain,file='fDfTrain.txt')";

   r["RMVA.RSVM.fDfTest"] = fDfTest;
   r["RMVA.RSVM.fWeightTest"] = fWeightTest;
   r << "write.table(RMVA.RSVM.fDfTest,file='fDfTest.txt')";

   //factors creations
   r["RMVA.RSVM.fFactorTrain"] = fFactorTrain;
   r << "RMVA.RSVM.fFactorTrain<-factor(RMVA.RSVM.fFactorTrain)";
   r["RMVA.RSVM.fFactorTest"] = fFactorTest;
   r << "RMVA.RSVM.fFactorTest<-factor(RMVA.RSVM.fFactorTest)";

   //SVM require a named vector
   ROOT::R::TRDataFrame ClassWeightsTrain;
   ClassWeightsTrain["background"]=Data()->GetNEvtBkgdTrain();
   ClassWeightsTrain["signal"]=Data()->GetNEvtSigTrain();
   r["RMVA.RSVM.ClassWeightsTrain"]=ClassWeightsTrain;   

   //Spectator creation
   r["RMVA.RSVM.fDfSpectators"] = fDfSpectators;


}

void MethodRSVM::Train()
{
   if (Data()->GetNTrainingEvents() == 0) Log() << kFATAL << "<Train> Data() has zero events" << Endl;

   r << "RMVA.RSVM.Model<-svm( x             = RMVA.RSVM.fDfTrain, \
                               y             = RMVA.RSVM.fFactorTrain,\
                               scale         = RMVA.RSVM.Scale,\
                               type          = RMVA.RSVM.Type,\
                               kernel        = RMVA.RSVM.Kernel, \
                               degree        = RMVA.RSVM.Degree, \
                               gamma         = RMVA.RSVM.Gamma,\
                               coef0         = RMVA.RSVM.Coef0,\
                               cost          = RMVA.RSVM.Cost,\
                               nu            = RMVA.RSVM.Nu,\
                               class.weights = RMVA.RSVM.ClassWeightsTrain,\
                               cachesize     = RMVA.RSVM.CacheSize,\
                               tolerance     = RMVA.RSVM.Tolerance,\
                               epsilon       = RMVA.RSVM.Epsilon,\
                               shrinking     = RMVA.RSVM.Shrinking,\
                               cross         = RMVA.RSVM.Cross,\
                               probability   = RMVA.RSVM.Probability,\
                               fitted        = RMVA.RSVM.Fitted)";
   r.SetVerbose(1);
   r << "summary(RMVA.RSVM.Model)";
   r.SetVerbose(0);
}

//_______________________________________________________________________
void MethodRSVM::DeclareOptions()
{
   DeclareOptionRef(fScale, "Scale", "A logical vector indicating the variables to be scaled. If\
                                       ‘scale’ is of length 1, the value is recycled as many times \
                                       as needed.  Per default, data are scaled internally (both ‘x’\
                                       and ‘y’ variables) to zero mean and unit variance. The center \
                                       and scale values are returned and used for later predictions.");
   DeclareOptionRef(fType, "Type", "‘svm’ can be used as a classification machine, as a \
                                     regression machine, or for novelty detection.  Depending of\
                                     whether ‘y’ is a factor or not, the default setting for\
                                     ‘type’ is ‘C-classification’ or ‘eps-regression’,\
                                     respectively, but may be overwritten by setting an explicit value.\
                                     Valid options are:\
                                      - ‘C-classification’\
                                      - ‘nu-classification’\
                                      - ‘one-classification’ (for novelty detection)\
                                      - ‘eps-regression’\
                                      - ‘nu-regression’");
   DeclareOptionRef(fKernel, "Kernel", "the kernel used in training and predicting. You might\
                                        consider changing some of the following parameters, depending on the kernel type.\
                                        linear: u'*v\
                                        polynomial: (gamma*u'*v + coef0)^degree\
                                        radial basis: exp(-gamma*|u-v|^2)\
                                        sigmoid: tanh(gamma*u'*v + coef0)");
   DeclareOptionRef(fDegree, "Degree", "parameter needed for kernel of type ‘polynomial’ (default: 3)");
   DeclareOptionRef(fGamma, "Gamma", "parameter needed for all kernels except ‘linear’ (default:1/(data dimension))");
   DeclareOptionRef(fCoef0, "Coef0", "parameter needed for kernels of type ‘polynomial’ and ‘sigmoid’ (default: 0)");
   DeclareOptionRef(fCost, "Cost", "cost of constraints violation (default: 1)-it is the ‘C’-constant of the regularization term in the Lagrange formulation.");
   DeclareOptionRef(fNu, "Nu", "parameter needed for ‘nu-classification’, ‘nu-regression’,and ‘one-classification’");
   DeclareOptionRef(fCacheSize, "CacheSize", "cache memory in MB (default 40)");
   DeclareOptionRef(fTolerance, "Tolerance", "tolerance of termination criterion (default: 0.001)");
   DeclareOptionRef(fEpsilon, "Epsilon", "epsilon in the insensitive-loss function (default: 0.1)");
   DeclareOptionRef(fShrinking, "Shrinking", "option whether to use the shrinking-heuristics (default:‘TRUE’)");
   DeclareOptionRef(fCross, "Cross", "if a integer value k>0 is specified, a k-fold cross validation on the training data is performed to assess the\
                                       quality of the model: the accuracy rate for classification and the Mean Squared Error for regression");
   DeclareOptionRef(fProbability, "Probability", "logical indicating whether the model should allow for probability predictions.");
   DeclareOptionRef(fFitted, "Fitted", "logical indicating whether the fitted values should be computed and included in the model or not (default: ‘TRUE’)");

}

//_______________________________________________________________________
void MethodRSVM::ProcessOptions()
{
   r["RMVA.RSVM.Scale"] = fScale;
   r["RMVA.RSVM.Type"] = fType;
   r["RMVA.RSVM.Kernel"] = fKernel;
   r["RMVA.RSVM.Degree"] = fDegree;
   r["RMVA.RSVM.Gamma"] = fGamma;
   r["RMVA.RSVM.Coef0"] = fCoef0;
   r["RMVA.RSVM.Cost"] = fCost;
   r["RMVA.RSVM.Nu"] = fNu;
   r["RMVA.RSVM.CacheSize"] = fCacheSize;
   r["RMVA.RSVM.Tolerance"] = fTolerance;
   r["RMVA.RSVM.Epsilon"] = fEpsilon;
   r["RMVA.RSVM.Shrinking"] = fShrinking;
   r["RMVA.RSVM.Cross"] = fCross;
   r["RMVA.RSVM.Probability"] = fProbability;
   r["RMVA.RSVM.Fitted"] = fFitted;

}

//_______________________________________________________________________
void MethodRSVM::TestClassification()
{
   Log() << kINFO << "Testing Classification RSVM METHOD  " << Endl;

   MethodBase::TestClassification();
}


//_______________________________________________________________________
Double_t MethodRSVM::GetMvaValue(Double_t *errLower, Double_t *errUpper)
{
   Double_t mvaValue;
   if (Data()->GetCurrentType() == Types::kTraining) {
      if (fProbResultForTrainSig.size() == 0) {
         r << "RMVA.RSVM.Predictor.Train.Prob<-predict(RMVA.RSVM.Model,RMVA.RSVM.fDfTrain,type='prob', decision.values = T, probability = T)"; //pridiction type prob use for ROC curves
         r["as.vector(attributes(RMVA.RSVM.Predictor.Train.Prob)$decision.values)"] >> fProbResultForTrainSig;
      }
      mvaValue = fProbResultForTrainSig[fMvaCounter];

      if (fMvaCounter < Data()->GetNTrainingEvents() - 1) fMvaCounter++;
      else fMvaCounter = 0;
   } else {
      if (fProbResultForTestSig.size() == 0) {
         r << "RMVA.RSVM.Predictor.Test.Prob <-predict(RMVA.RSVM.Model,RMVA.RSVM.fDfTest,type='prob', decision.values = T, probability = T)";
         r["as.vector(attributes(RMVA.RSVM.Predictor.Test.Prob)$decision.values)"] >> fProbResultForTestSig;
      }
      mvaValue = fProbResultForTestSig[fMvaCounter];
      if (fMvaCounter < Data()->GetNTestEvents() - 1) fMvaCounter++;
      else fMvaCounter = 0;
   }
   return mvaValue;
}

//_______________________________________________________________________
void MethodRSVM::GetHelpMessage() const
{
// get help message text
//
// typical length of text line:
//         "|--------------------------------------------------------------|"
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Short description:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "Decision Trees and Rule-Based Models " << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance optimisation:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance tuning via configuration options:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "<None>" << Endl;
}

