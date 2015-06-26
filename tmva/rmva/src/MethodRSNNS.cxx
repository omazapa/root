// @(#)root/rmva $Id$
// Author: Omar Zapata 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodRSNNS                                                           *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Neural Networks in R using the Stuttgart Neural Network Simulator         *
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
#include "TMVA/MethodRSNNS.h"
#include "TMVA/Tools.h"
#include "TMVA/Ranking.h"
#include "TMVA/Types.h"
#include "TMVA/PDF.h"
#include "TMVA/ClassifierFactory.h"

#include "TMVA/Results.h"

using namespace TMVA;

REGISTER_METHOD(RSNNS)

ClassImp(MethodRSNNS)


//_______________________________________________________________________
MethodRSNNS::MethodRSNNS(const TString &jobName,
                         const TString &methodTitle,
                         DataSetInfo &dsi,
                         const TString &theOption,
                         TDirectory *theTargetDir) :
   RMethodBase(jobName, Types::kRSNNS, methodTitle, dsi, theOption, theTargetDir), fMvaCounter(0)
{
   fNetType = methodTitle;
   if (fNetType != "RMLP") {
      Log() << kFATAL << " Unknow Method" + fNetType
            << Endl;
      return;
   }

   // standard constructor for the RSNNS
   //RSNNS Options for all NN methods
   fSize = "c(5)";
   fMaxit = 100;

   fInitFunc = "Randomize_Weights";
   fInitFuncParams = "c(-0.3,0.3)"; //the maximun number of pacameter is 5 see RSNNS::getSnnsRFunctionTable() type 6

   fLearnFunc = "Std_Backpropagation"; //
   fLearnFuncParams = "c(0.2,0)";

   fUpdateFunc = "Topological_Order";
   fUpdateFuncParams = "c(0)";

   fHiddenActFunc = "Act_Logistic";
   fShufflePatterns = kTRUE;
   fLinOut = kFALSE;
   fPruneFunc = "NULL";
   fPruneFuncParams = "NULL";
}

//_______________________________________________________________________
MethodRSNNS::MethodRSNNS(DataSetInfo &theData, const TString &theWeightFile, TDirectory *theTargetDir)
   : RMethodBase(Types::kRSNNS, theData, theWeightFile, theTargetDir), fMvaCounter(0)
{
   fNetType = "RMLP"; //GetMethodName();//GetMethodName() is not returning RMLP is reting MethodBase why?
   if (fNetType != "RMLP") {
      Log() << kFATAL << " Unknow Method = " + fNetType
            << Endl;
      return;
   }

   // standard constructor for the RSNNS
   //RSNNS Options for all NN methods
   fSize = "c(5)";
   fMaxit = 100;

   fInitFunc = "Randomize_Weights";
   fInitFuncParams = "c(-0.3,0.3)"; //the maximun number of pacameter is 5 see RSNNS::getSnnsRFunctionTable() type 6

   fLearnFunc = "Std_Backpropagation"; //
   fLearnFuncParams = "c(0.2,0)";

   fUpdateFunc = "Topological_Order";
   fUpdateFuncParams = "c(0)";

   fHiddenActFunc = "Act_Logistic";
   fShufflePatterns = kTRUE;
   fLinOut = kFALSE;
   fPruneFunc = "NULL";
   fPruneFuncParams = "NULL";

}


//_______________________________________________________________________
MethodRSNNS::~MethodRSNNS(void)
{
}

//_______________________________________________________________________
Bool_t MethodRSNNS::HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets)
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   return kFALSE;
}


//_______________________________________________________________________
void     MethodRSNNS::Init()
{
   if (!r.Require("Rcpp")) {
      Error("Init", "R's package Rcpp can not be loaded.");
      Log() << kFATAL << " R's package Rcpp can not be loaded."
            << Endl;
      return;
   }
   if (!r.IsInstalled("RSNNS")) {
      Error("Init", "R's package RSNNS is not installed.");
      Log() << kFATAL << " R's package RSNNS is not installed."
            << Endl;
      return;
   }

   if (!r.Require("RSNNS")) {
      Error("Init", "R's package RSNNS can not be loaded.");
      Log() << kFATAL << " R's package RSNNS can not be loaded."
            << Endl;
      return;
   }

   if (!r.Require("caret")) {
      Error("Init", "R's package caret can not be loaded.");
      Log() << kFATAL << " R's package caret can not be loaded."
            << Endl;
      return;
   }
   //Paassing Data to R's environment
   //NOTE:need improved names in R's environment using JobName of TMVA
   if (fNetType == "RMLP") {
      r["RMVA.RSNNSRMLP.fDfTrain"] = fDfTrain;
      r["RMVA.RSNNSRMLP.fWeightTrain"] = fWeightTrain;

      r["RMVA.RSNNSRMLP.fDfTest"] = fDfTest;
      r["RMVA.RSNNSRMLP.fWeightTest"] = fWeightTest;
   }
   //factors creations
   //RSNNS mlp require a numeric factor then background=0 signal=1 from fFactorTrain/fFactorTest
   UInt_t size = fFactorTrain.size();
   std::vector<UInt_t>  fFactorNumeric(size);

   for (UInt_t i = 0; i < size; i++) {
      if (fFactorTrain[i] == "signal") fFactorNumeric[i] = 1;
      else fFactorNumeric[i] = 0;
   }
   if (fNetType == "RMLP") {
      r["RMVA.RSNNSRMLP.fFactorTrain"] = fFactorNumeric;
   }
   fFactorNumeric.clear();
   size = fFactorTest.size();
   fFactorNumeric.resize(size);
   for (UInt_t i = 0; i < size; i++) {
      if (fFactorTest[i] == "signal") fFactorNumeric[i] = 1;
      else fFactorNumeric[i] = 0;
   }
   if (fNetType == "RMLP") {
      r["RMVA.RSNNSRMLP.fFactorTest"] = fFactorNumeric;

      //Spectator creation
      r["RMVA.RSNNSRMLP.fDfSpectators"] = fDfSpectators;
   }
}

void MethodRSNNS::Train()
{
   if (Data()->GetNTrainingEvents() == 0) Log() << kFATAL << "<Train> Data() has zero events" << Endl;
   if (fNetType == "RMLP") {
      r << "RMVA.RSNNSRMLP.Model<-mlp(x=RMVA.RSNNSRMLP.fDfTrain,y=RMVA.RSNNSRMLP.fFactorTrain,size = RMVA.RSNNSRMLP.Size,maxit = RMVA.RSNNSRMLP.Maxit,\
                                  initFunc = RMVA.RSNNSRMLP.InitFunc,initFuncParams = RMVA.RSNNSRMLP.InitFuncParams,\
                                  learnFunc = RMVA.RSNNSRMLP.LearnFunc,learnFuncParams = RMVA.RSNNSRMLP.LearnFuncParams,\
                                  updateFunc = RMVA.RSNNSRMLP.UpdateFunc,updateFuncParams = RMVA.RSNNSRMLP.UpdateFuncParams,\
                                  hiddenActFunc = RMVA.RSNNSRMLP.HiddenActFunc, shufflePatterns = RMVA.RSNNSRMLP.ShufflePatterns, linOut =RMVA.RSNNSRMLP.LinOut,\
                                  pruneFunc = RMVA.RSNNSRMLP.PruneFunc, pruneFuncParams = RMVA.RSNNSRMLP.PruneFuncParams)";
      r.SetVerbose(1);
      r<<"RMVA.RSNNSRMLP.Model";
      r.SetVerbose(0);
   }
}

//_______________________________________________________________________
void MethodRSNNS::DeclareOptions()
{
   //RSNNS Options for all NN methods
//       TVectorF  fSize;//number of units in the hidden layer(s)
   DeclareOptionRef(fSize, "Size", "number of units in the hidden layer(s)");
   DeclareOptionRef(fMaxit, "Maxit", "Maximum of iterations to learn");

   DeclareOptionRef(fInitFunc, "InitFunc", "the initialization function to use");
   DeclareOptionRef(fInitFuncParams, "InitFuncParams", "the parameters for the initialization function");

   DeclareOptionRef(fLearnFunc, "LearnFunc", "the learning function to use");
   DeclareOptionRef(fLearnFuncParams, "LearnFuncParams", "the parameters for the learning function");

   DeclareOptionRef(fUpdateFunc, "UpdateFunc", "the update function to use");
   DeclareOptionRef(fUpdateFuncParams, "UpdateFuncParams", "the parameters for the update function");

   DeclareOptionRef(fHiddenActFunc, "HiddenActFunc", "the activation function of all hidden units");
   DeclareOptionRef(fShufflePatterns, "ShufflePatterns", "should the patterns be shuffled?");
   DeclareOptionRef(fLinOut, "LinOut", "sets the activation function of the output units to linear or logistic");

   DeclareOptionRef(fPruneFunc, "PruneFunc", "the prune function to use");
   DeclareOptionRef(fPruneFuncParams, "PruneFuncParams", "the parameters for the pruning function. Unlike the\
                                                     other functions, these have to be given in a named list. See\
                                                     the pruning demos for further explanation.the update function to use");

}

//_______________________________________________________________________
void MethodRSNNS::ProcessOptions()
{
   if (fMaxit <= 0) {
      Log() << kERROR << " fMaxit <=0... that does not work !! "
            << " I set it to 50 .. just so that the program does not crash"
            << Endl;
      fMaxit = 1;
   }
   // standard constructor for the RSNNS
   //RSNNS Options for all NN methods
   if (fNetType == "RMLP") {
      r << "RMVA.RSNNSRMLP.Size<-" + fSize;
      r["RMVA.RSNNSRMLP.Maxit"] = fMaxit;

      r["RMVA.RSNNSRMLP.InitFunc"] = fInitFunc;
      r << "RMVA.RSNNSRMLP.InitFuncParams<-" + fInitFuncParams;

      r["RMVA.RSNNSRMLP.LearnFunc"] = fLearnFunc;
      r << "RMVA.RSNNSRMLP.LearnFuncParams<-" + fLearnFuncParams;

      r["RMVA.RSNNSRMLP.UpdateFunc"] = fUpdateFunc;
      r << "RMVA.RSNNSRMLP.UpdateFuncParams<-" + fUpdateFuncParams;

      r["RMVA.RSNNSRMLP.HiddenActFunc"] = fHiddenActFunc;
      r["RMVA.RSNNSRMLP.ShufflePatterns"] = fShufflePatterns;
      r["RMVA.RSNNSRMLP.LinOut"] = fLinOut;

      r << "RMVA.RSNNSRMLP.PruneFunc<-" + fPruneFunc;
      r << "RMVA.RSNNSRMLP.PruneFuncParams<-" + fPruneFuncParams;
   }


}

//_______________________________________________________________________
void MethodRSNNS::TestClassification()
{
   Log() << kINFO << "Testing Classification " << fNetType << " METHOD  " << Endl;

   MethodBase::TestClassification();
}


//_______________________________________________________________________
Double_t MethodRSNNS::GetMvaValue(Double_t *errLower, Double_t *errUpper)
{
   Double_t mvaValue;
   if (Data()->GetCurrentType() == Types::kTraining) {
      if (fProbResultForTrainSig.size() == 0) {
         if (fNetType == "RMLP") {
            r << "RMVA.RSNNSRMLP.Predictor.Train.Prob<-predict(RMVA.RSNNSRMLP.Model,RMVA.RSNNSRMLP.fDfTrain,type='prob')";
            r["as.vector(RMVA.RSNNSRMLP.Predictor.Train.Prob[,1])"] >> fProbResultForTrainSig;
         }
      }
      mvaValue = fProbResultForTrainSig[fMvaCounter];
      if (fMvaCounter < Data()->GetNTrainingEvents() - 1) fMvaCounter++;
      else fMvaCounter = 0;
   } else {
      if (fProbResultForTestSig.size() == 0) {
         if (fNetType == "RMLP") {
            r << "RMVA.RSNNSRMLP.Predictor.Test.Prob<-predict(RMVA.RSNNSRMLP.Model,RMVA.RSNNSRMLP.fDfTest,type='prob')";
            r["as.vector(RMVA.RSNNSRMLP.Predictor.Test.Prob[,1])"] >> fProbResultForTestSig;
         }
      }

      mvaValue = fProbResultForTestSig[fMvaCounter];

      if (fMvaCounter < Data()->GetNTestEvents() - 1) fMvaCounter++;
      else fMvaCounter = 0;
   }
   return mvaValue;
}


//_______________________________________________________________________
void MethodRSNNS::GetHelpMessage() const
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

