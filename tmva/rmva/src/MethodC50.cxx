// @(#)root/rmva $Id$
// Author: Omar Zapata 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodC50                                                             *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Decision Trees and Rule-Based Models                                      *
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
#include "TMVA/MethodC50.h"
#include "TMVA/Tools.h"
#include "TMVA/Ranking.h"
#include "TMVA/Types.h"
#include "TMVA/PDF.h"
#include "TMVA/ClassifierFactory.h"

#include "TMVA/Results.h"

using namespace TMVA;

REGISTER_METHOD(C50)

ClassImp(MethodC50)

//creating an Instance
Bool_t MethodC50::IsModuleLoaded=ROOT::R::TRInterface::Instance().Require("C50");

//_______________________________________________________________________
MethodC50::MethodC50(const TString &jobName,
                     const TString &methodTitle,
                     DataSetInfo &dsi,
                     const TString &theOption,
                     TDirectory *theTargetDir) :RMethodBase(jobName, Types::kC50, methodTitle, dsi, theOption, theTargetDir),
   fNTrials(1), 
   fRules(kFALSE), 
   fMvaCounter(0),
   predict("predict.C5.0"),
   C50("C5.0"),
   C50Control("C5.0Control"),
   asfactor("as.factor")
{
   // standard constructor for the C50

   //C5.0Control options
   fControlSubset = kTRUE;
   fControlBands = 0;
   fControlWinnow = kFALSE;
   fControlNoGlobalPruning = kFALSE;
   fControlCF = 0.25;
   fControlMinCases = 2;
   fControlFuzzyThreshold = kFALSE;
   fControlSample = 0;
   r["sample.int(4096, size = 1) - 1L"] >> fControlSeed;
   fControlEarlyStopping = kTRUE;

}

//_______________________________________________________________________
MethodC50::MethodC50(DataSetInfo &theData, const TString &theWeightFile, TDirectory *theTargetDir)
   : RMethodBase(Types::kC50, theData, theWeightFile, theTargetDir),
   fNTrials(1), 
   fRules(kFALSE), 
   fMvaCounter(0),
   predict("predict.C5.0"),
   C50("C5.0"),
   C50Control("C5.0Control"),
   asfactor("as.factor")
{

   // constructor from weight file
   fControlSubset = kTRUE;
   fControlBands = 0;
   fControlWinnow = kFALSE;
   fControlNoGlobalPruning = kFALSE;
   fControlCF = 0.25;
   fControlMinCases = 2;
   fControlFuzzyThreshold = kFALSE;
   fControlSample = 0;
   r["sample.int(4096, size = 1) - 1L"] >> fControlSeed;
   fControlEarlyStopping = kTRUE;

}


//_______________________________________________________________________
MethodC50::~MethodC50(void)
{
}

//_______________________________________________________________________
Bool_t MethodC50::HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets)
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   return kFALSE;
}


//_______________________________________________________________________
void     MethodC50::Init()
{

   if (!IsModuleLoaded) {
      Error("Init", "R's package C50 can not be loaded.");
      Log() << kFATAL << " R's package C50 can not be loaded."
            << Endl;
      return;
   }

    
//    if (!r.Require("partykit")) {
//       Error("Init", "R's package partykit can not be loaded.");
//       Log() << kFATAL << " R's package partykit can not be loaded."
//             << Endl;
//       return;
//    }
//    if (!r.IsInstalled("C50")) {
//       Error("Init", "R's package C50 is not installed.");
//       Log() << kFATAL << " R's package C50 is not installed."
//             << Endl;
//       return;
//    }
// 
//    if (!r.Require("C50")) {
//       Error("Init", "R's package C50 can not be loaded.");
//       Log() << kFATAL << " R's package C50 can not be loaded."
//             << Endl;
//       return;
//    }
//    //Paassing Data to R's environment
//    //NOTE:need improved names in R's environment using JobName of TMVA
//    r["RMVA.C50.fDfTrain"] = fDfTrain;
//    r["RMVA.C50.fWeightTrain"] = fWeightTrain;
//    r << "write.table(RMVA.C50.fDfTrain,file='fDfTrain.txt')";
// 
//    r["RMVA.C50.fDfTest"] = fDfTest;
//    r["RMVA.C50.fWeightTest"] = fWeightTest;
//    r << "write.table(RMVA.C50.fDfTest,file='fDfTest.txt')";
// 
//    //factors creations
//    r["RMVA.C50.fFactorTrain"] = fFactorTrain;
//    r << "RMVA.C50.fFactorTrain<-factor(RMVA.C50.fFactorTrain)";
//    r["RMVA.C50.fFactorTest"] = fFactorTest;
//    r << "RMVA.C50.fFactorTest<-factor(RMVA.C50.fFactorTest)";
// 
//    //Spectator creation
//    r["RMVA.C50.fDfSpectators"] = fDfSpectators;
// 
//    r["RMVA.C50.fCounter"] = 0;



}

void MethodC50::Train()
{
   if (Data()->GetNTrainingEvents() == 0) Log() << kFATAL << "<Train> Data() has zero events" << Endl;

   fModel=C50(ROOT::R::Label["x"]=fDfTrain, \
              ROOT::R::Label["y"]=asfactor(fFactorTrain), \
              ROOT::R::Label["trials"]=fNTrials, \
              ROOT::R::Label["rules"]=fRules, \
              ROOT::R::Label["weights"]=fWeightTrain, \
              ROOT::R::Label["control"]=fModelControl);
}

//_______________________________________________________________________
void MethodC50::DeclareOptions()
{
   //
   DeclareOptionRef(fNTrials, "NTrials", "An integer specifying the number of boosting iterations");
   DeclareOptionRef(fRules, "Rules", "A logical: should the tree be decomposed into a rule-basedmodel?");

   //C5.0Control Options
   DeclareOptionRef(fControlSubset, "ControlSubset", "A logical: should the model evaluate groups of discrete \
                                      predictors for splits? Note: the C5.0 command line version defaults this \
                                      parameter to ‘FALSE’, meaning no attempted gropings will be evaluated \
                                      during the tree growing stage.");
   DeclareOptionRef(fControlBands, "ControlBands", "An integer between 2 and 1000. If ‘TRUE’, the model orders \
                                     the rules by their affect on the error rate and groups the \
                                     rules into the specified number of bands. This modifies the \
                                     output so that the effect on the error rate can be seen for \
                                     the groups of rules within a band. If this options is \
                                     selected and ‘rules = kFALSE’, a warning is issued and ‘rules’ \
                                     is changed to ‘kTRUE’.");
   DeclareOptionRef(fControlWinnow, "ControlWinnow", "A logical: should predictor winnowing (i.e feature selection) be used?");
   DeclareOptionRef(fControlNoGlobalPruning, "ControlNoGlobalPruning", "A logical to toggle whether the final, global pruning \
                                                                         step to simplify the tree.");
   DeclareOptionRef(fControlCF, "ControlCF", "A number in (0, 1) for the confidence factor.");
   DeclareOptionRef(fControlMinCases, "ControlMinCases", "an integer for the smallest number of samples that must be \
                                                           put in at least two of the splits.");

   DeclareOptionRef(fControlFuzzyThreshold, "ControlFuzzyThreshold", "A logical toggle to evaluate possible advanced splits \
                                                                      of the data. See Quinlan (1993) for details and examples.");
   DeclareOptionRef(fControlSample, "ControlSample", "A value between (0, .999) that specifies the random \
                                                       proportion of the data should be used to train the model. By \
                                                       default, all the samples are used for model training. Samples \
                                                       not used for training are used to evaluate the accuracy of \
                                                       the model in the printed output.");
   DeclareOptionRef(fControlSeed, "ControlSeed", " An integer for the random number seed within the C code.");
   DeclareOptionRef(fControlEarlyStopping, "ControlEarlyStopping", " A logical to toggle whether the internal method for \
                                                                      stopping boosting should be used.");


}

//_______________________________________________________________________
void MethodC50::ProcessOptions()
{
   if (fNTrials <= 0) {
      Log() << kERROR << " fNTrials <=0... that does not work !! "
            << " I set it to 1 .. just so that the program does not crash"
            << Endl;
      fNTrials = 1;
   }
   fModelControl=C50Control(ROOT::R::Label["subset"]=fControlSubset, \
                            ROOT::R::Label["bands"]=fControlBands, \
                            ROOT::R::Label["winnow"]=fControlWinnow, \
                            ROOT::R::Label["noGlobalPruning"]=fControlNoGlobalPruning, \
                            ROOT::R::Label["CF"]=fControlCF, \
                            ROOT::R::Label["minCases"]=fControlMinCases, \
                            ROOT::R::Label["fuzzyThreshold"]=fControlFuzzyThreshold, \
                            ROOT::R::Label["sample"]=fControlSample, \
                            ROOT::R::Label["seed"]=fControlSeed, \
                            ROOT::R::Label["earlyStopping"]=fControlEarlyStopping);
}

//_______________________________________________________________________
void MethodC50::TestClassification()
{
   Log() << kINFO << "Testing Classification C50 METHOD  " << Endl;
//    r << "RMVA.C50.Predictor.Test.Class <-predict.C5.0(RMVA.C50.Model,RMVA.C50.fDfTest,type='class')";

   gSystem->MakeDirectory("C50");
   gSystem->MakeDirectory("C50/plots");

//    if (r.IsInstalled("ROCR")) {
//       if (r.Require("ROCR")) {
//          //calculation ROC curves https://ifordata.wordpress.com/category/predictions-in-r/
//          r << "RMVA.C50.ROCRPredictionSig<-ROCR::prediction(predictions = RMVA.C50.Predictor.Test.Prob[,2], labels=RMVA.C50.fFactorTest)";
//          //at the moment I am using default performance method for ROCR but it mush be an option to parse from booking
// 
//          //plots TPR (True Positive Rate) and FPR (False Positive Rate)
//          r << "RMVA.C50.ROCRPerformanceSig<-ROCR::performance(RMVA.C50.ROCRPredictionSig, measure='tpr', x.measure='fpr')";
// 
//          r << "setEPS()";
//          r << "postscript('C50/plots/C50ROCSigTPR-FPR.eps')";
//          r << "plot(RMVA.C50.ROCRPerformanceSig, main='ROC curve for Signal', col='darkmagenta', lwd=3)";
//          r << "abline(a=0, b=1, lwd=2, lty=2)";
//          r << "dev.off()";
// 
//          //Getting AUC  (Area Under the Curve)
//          r << "RMVA.C50.ROCRPerformanceSig<-ROCR::performance(RMVA.C50.ROCRPredictionSig, measure='auc')";
//          Log() << "----------------------------------" << Endl;
//          Float_t ROC_AUC;
//          r["RMVA.C50.ROCRPerformanceSig@y.values[[1]]"] >> ROC_AUC;
//          r.SetVerbose(1);
//          Log() << gTools().Color("bold") << "Area under the ROC curve " << gTools().Color("reset") << "= " << ROC_AUC << " (see ranking below) " << Endl;
//          Log() << "0.9-1.0 -> A (perfect)" << Endl;
//          Log() << "0.8-0.9 -> B (excellent)" << Endl;
//          Log() << "0.7-0.8 -> C (fair)" << Endl;
//          Log() << "0.6-0.7 -> D (poor)" << Endl;
//          Log() << "0.5-0.6 -> F (no value)" << Endl;
//          Log() << "----------------------------------" << Endl;
//          r.SetVerbose(0);
//       }
//    }
//    if (r.IsInstalled("caret")) {
//       if (r.Require("caret")) {
//          //performing confusion matrix with the analysis of tests
//          r.SetVerbose(1);
//          r << "RMVA.C50.TestConfusionMatrix<-caret::confusionMatrix(RMVA.C50.Predictor.Test.Class,RMVA.C50.fFactorTest,positive='signal')";
//          r.SetVerbose(0);
//       }
//    }
   MethodBase::TestClassification();
}


//_______________________________________________________________________
Double_t MethodC50::GetMvaValue(Double_t *errLower, Double_t *errUpper)
{
   Double_t mvaValue;
   const TMVA::Event *ev=GetEvent();
   const UInt_t nvar = DataInfo().GetNVariables();
   ROOT::R::TRDataFrame fDfEvent;
   for(UInt_t i=0;i<nvar;i++)
   {
      fDfEvent[GetInputLabel( i ).Data()]=ev->GetValues()[i];
   }
   TVectorD result=predict(fModel,fDfEvent,ROOT::R::Label["type"]="prob");
   mvaValue=result[1];
   return mvaValue;
}

//_______________________________________________________________________
void MethodC50::GetHelpMessage() const
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

