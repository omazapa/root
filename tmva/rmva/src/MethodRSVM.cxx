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
MethodRSVM::MethodRSVM( const TString& jobName,
                      const TString& methodTitle,
                      DataSetInfo& dsi,
                      const TString& theOption,
                      TDirectory* theTargetDir ) :
    RMethodBase( jobName, Types::kRSVM, methodTitle, dsi, theOption, theTargetDir ),fMvaCounter(0)
{
    // standard constructor for the RSVM


}

//_______________________________________________________________________
MethodRSVM::MethodRSVM( DataSetInfo& theData, const TString& theWeightFile, TDirectory* theTargetDir )
    : RMethodBase( Types::kRSVM, theData, theWeightFile, theTargetDir ),fMvaCounter(0)
{

}


//_______________________________________________________________________
MethodRSVM::~MethodRSVM( void )
{
}

//_______________________________________________________________________
Bool_t MethodRSVM::HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets )
{
    if (type == Types::kClassification && numberClasses == 2) return kTRUE;
    return kFALSE;
}


//_______________________________________________________________________
void     MethodRSVM::Init()
{
    if(!r.IsInstalled("e1071"))
    {
        Error( "Init","R's package e1071 is not installed.");
        Log() << kFATAL << " R's package e1071 is not installed."
              << Endl;
        return;
    }

    if(!r.Require("e1071"))
    {
        Error("Init","R's package e1071 can not be loaded.");
        Log() << kFATAL << " R's package e1071 can not be loaded."
              << Endl;
        return;
    }
    //Paassing Data to R's environment
    //NOTE:need improved names in R's environment using JobName of TMVA
    r["RMVA.RSVM.fDfTrain"]=fDfTrain;
    r["RMVA.RSVM.fWeightTrain"]=fWeightTrain;
    r<<"write.table(RMVA.RSVM.fDfTrain,file='fDfTrain.txt')";

    r["RMVA.RSVM.fDfTest"]=fDfTest;
    r["RMVA.RSVM.fWeightTest"]=fWeightTest;
    r<<"write.table(RMVA.RSVM.fDfTest,file='fDfTest.txt')";

    //factors creations
    r["RMVA.RSVM.fFactorTrain"]=fFactorTrain;
    r<<"RMVA.RSVM.fFactorTrain<-factor(RMVA.RSVM.fFactorTrain)";
    r["RMVA.RSVM.fFactorTest"]=fFactorTest;
    r<<"RMVA.RSVM.fFactorTest<-factor(RMVA.RSVM.fFactorTest)";

    //Spectator creation
    r["RMVA.RSVM.fDfSpectators"]=fDfSpectators;

    r["RMVA.RSVM.fCounter"]=0;



}

void MethodRSVM::Train()
{
    if (Data()->GetNTrainingEvents()==0) Log() << kFATAL << "<Train> Data() has zero events" << Endl;

    r<<"RMVA.RSVM.Model<-svm( x             = RMVA.RSVM.fDfTrain, \
                              y             = RMVA.RSVM.fFactorTrain)";
//                              class.weights = RMVA.RSVM.fWeightTrain)";
//    r.SetVerbose(1);
    r<<"summary(RMVA.RSVM.Model)";
   
//    Log() << kWARNING << " RMVA.RSVM.Predictor.ClassResultForTest SIze. = "<<fClassResultForTest.size()<< Endl;        
//    r.SetVerbose(0);
}

//_______________________________________________________________________
void MethodRSVM::DeclareOptions()
{

}

//_______________________________________________________________________
void MethodRSVM::ProcessOptions()
{

}

//_______________________________________________________________________
void MethodRSVM::TestClassification()
{
    Log()<<kINFO<<"Testing Classification RSVM METHOD  "<<Endl;
    
    gSystem->MakeDirectory("RSVM");
    gSystem->MakeDirectory("RSVM/plots");

    MethodBase::TestClassification();
}


//_______________________________________________________________________
Double_t MethodRSVM::GetMvaValue( Double_t* errLower, Double_t* errUpper)
{
    Double_t mvaValue;
    if(Data()->GetCurrentType()==Types::kTraining) 
    {
        if(fProbResultForTrainSig.size()==0)
        {
           r<<"RMVA.RSVM.Predictor.Train.Prob<-predict(RMVA.RSVM.Model,RMVA.RSVM.fDfTrain,type='prob', decision.values = T, probability = T)";//pridiction type prob use for ROC curves
           r["as.vector(attributes(RMVA.RSVM.Predictor.Train.Prob)$decision.values)"]>>fProbResultForTrainSig;
        }
          mvaValue=fProbResultForTrainSig[fMvaCounter];
       
       if(fMvaCounter < Data()->GetNTrainingEvents()-1) fMvaCounter++;
       else fMvaCounter=0;
    }else
    {
        if(fProbResultForTestSig.size()==0)
        {
        r<<"RMVA.RSVM.Predictor.Test.Prob <-predict(RMVA.RSVM.Model,RMVA.RSVM.fDfTest,type='prob', decision.values = T, probability = T)";
        r["as.vector(attributes(RMVA.RSVM.Predictor.Test.Prob)$decision.values)"]>>fProbResultForTestSig;
        }
       mvaValue=fProbResultForTestSig[fMvaCounter];
       if(fMvaCounter < Data()->GetNTestEvents()-1) fMvaCounter++;
       else fMvaCounter=0;
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

