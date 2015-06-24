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
MethodRSNNS::MethodRSNNS( const TString& jobName,
                      const TString& methodTitle,
                      DataSetInfo& dsi,
                      const TString& theOption,
                      TDirectory* theTargetDir ) :
    RMethodBase( jobName, Types::kRSNNS, methodTitle, dsi, theOption, theTargetDir ),fMvaCounter(0)
{
    fNetType=this->GetMethodName();
    if(fNetType!="RMLP")
    {
        Log() << kFATAL << " Unknow Method"+fNetType
              << Endl;
        return;        
    }

    // standard constructor for the RSNNS
       //RSNNS Options for all NN methods
       fSize.ResizeTo(1);
       fSize[0]=5;
       fMaxit=100;
       
       fInitFunc="Randomize_Weights";
       fInitFuncParams[0]=-0.3;
       fInitFuncParams[1]=0.3;
       
       fLearnFunc="Std_Backpropagation";
       fLearnFuncParams[0]=0.2;
       fLearnFuncParams[1]=0;
       
       fUpdateFunc="Topological_Order";
       fUpdateFuncParams.ResizeTo(1);
       fUpdateFuncParams[0]=0;

       fHiddenActFunc="Act_Logistic";
       fShufflePatterns=kTRUE;
       fLinOut=kFALSE;
}

//_______________________________________________________________________
MethodRSNNS::MethodRSNNS( DataSetInfo& theData, const TString& theWeightFile, TDirectory* theTargetDir )
    : RMethodBase( Types::kRSNNS, theData, theWeightFile, theTargetDir ),fMvaCounter(0)
{
    fNetType=this->GetMethodName();
    if(fNetType!="RMLP")
    {
        Log() << kFATAL << " Unknow Method"+fNetType
              << Endl;
        return;        
    }

    // standard constructor for the RSNNS
       //RSNNS Options for all NN methods
       fSize.ResizeTo(1);
       fSize[0]=5;
       fMaxit=100;
       
       fInitFunc="Randomize_Weights";
       fInitFuncParams[0]=-0.3;
       fInitFuncParams[1]=0.3;
       
       fLearnFunc="Std_Backpropagation";
       fLearnFuncParams[0]=0.2;
       fLearnFuncParams[1]=0;
       
       fUpdateFunc="Topological_Order";
       fUpdateFuncParams.ResizeTo(1);
       fUpdateFuncParams[0]=0;

       fHiddenActFunc="Act_Logistic";
       fShufflePatterns=kTRUE;
       fLinOut=kFALSE;

}


//_______________________________________________________________________
MethodRSNNS::~MethodRSNNS( void )
{
}

//_______________________________________________________________________
Bool_t MethodRSNNS::HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets )
{
    if (type == Types::kClassification && numberClasses == 2) return kTRUE;
    return kFALSE;
}


//_______________________________________________________________________
void     MethodRSNNS::Init()
{
    if(!r.Require("Rcpp"))
    {
        Error("Init","R's package Rcpp can not be loaded.");
        Log() << kFATAL << " R's package Rcpp can not be loaded."
              << Endl;
        return;
    }
    if(!r.IsInstalled("RSNNS"))
    {
        Error( "Init","R's package RSNNS is not installed.");
        Log() << kFATAL << " R's package RSNNS is not installed."
              << Endl;
        return;
    }

    if(!r.Require("RSNNS"))
    {
        Error("Init","R's package RSNNS can not be loaded.");
        Log() << kFATAL << " R's package RSNNS can not be loaded."
              << Endl;
        return;
    }
    
    if(!r.Require("caret"))
    {
        Error("Init","R's package caret can not be loaded.");
        Log() << kFATAL << " R's package caret can not be loaded."
              << Endl;
        return;
    }
    //Paassing Data to R's environment
    //NOTE:need improved names in R's environment using JobName of TMVA
    r["RMVA.RSNNS.fDfTrain"]=fDfTrain;
    r["RMVA.RSNNS.fWeightTrain"]=fWeightTrain;
    r<<"write.table(RMVA.RSNNS.fDfTrain,file='fDfTrain.txt')";

    r["RMVA.RSNNS.fDfTest"]=fDfTest;
    r["RMVA.RSNNS.fWeightTest"]=fWeightTest;
    r<<"write.table(RMVA.RSNNS.fDfTest,file='fDfTest.txt')";

    //factors creations
    //RSNNS mlp require a numeric factor then background=0 signal=1 from fFactorTrain/fFactorTest
    UInt_t size=fFactorTrain.size();
    std::vector<UInt_t>  fFactorNumeric(size);
    
    for(UInt_t i=0;i<size;i++)
    {
        if(fFactorTrain[i]=="signal") fFactorNumeric[i]=1;
        else fFactorNumeric[i]=0;
    }
    r["RMVA.RSNNS.fFactorTrain"]=fFactorNumeric;
    fFactorNumeric.clear();    
    size=fFactorTest.size();
    fFactorNumeric.resize(size);
    for(UInt_t i=0;i<size;i++)
    {
        if(fFactorTest[i]=="signal") fFactorNumeric[i]=1;
        else fFactorNumeric[i]=0;
    }    
    r["RMVA.RSNNS.fFactorTest"]=fFactorNumeric;

    //Spectator creation
    r["RMVA.RSNNS.fDfSpectators"]=fDfSpectators;

    r["RMVA.RSNNS.fCounter"]=0;
        
}

void MethodRSNNS::Train()
{
    if (Data()->GetNTrainingEvents()==0) Log() << kFATAL << "<Train> Data() has zero events" << Endl;
    r<<"RMVA.RSNNS.Model<-mlp(x=RMVA.RSNNS.fDfTrain,y=RMVA.RSNNS.fFactorTrain,size = c(3,5),maxit = 200)";
    r.SetVerbose(1);
    r<<"RMVA.RSNNS.Model";
    r.SetVerbose(0);
//    r<<"RMVA.RSNNS.Predictor.Train.Prob<-predict(RMVA.RSNNS.Model,RMVA.RSNNS.fDfTrain,type='prob')";
}

//_______________________________________________________________________
void MethodRSNNS::DeclareOptions()
{
       //RSNNS Options for all NN methods
//       TVectorF  fSize;//number of units in the hidden layer(s)
    DeclareOptionRef(fMaxit, "Maxit", "Maximum of iterations to learn");
    
    DeclareOptionRef(fInitFunc, "InitFunc", "the initialization function to use");

    fInitFuncParams=new Float_t[5];//the maximun number of pacameter is 5 see RSNNS::getSnnsRFunctionTable() type 6
    DeclareOptionRef(fInitFuncParams,5, "InitFuncParams", "the parameters for the initialization function"); 
       
    DeclareOptionRef(fLearnFunc, "LearnFunc", "the learning function to use");
    fLearnFuncParams=new Float_t[5];
    DeclareOptionRef(fLearnFuncParams,5, "LearnFuncParams", "the parameters for the learning function"); 

    DeclareOptionRef(fUpdateFunc, "UpdateFunc", "the update function to use");  
//    TVectorF fUpdateFuncParams;//the parameters for the update function

    DeclareOptionRef(fHiddenActFunc, "HiddenActFunc", "the activation function of all hidden units");  
    DeclareOptionRef(fShufflePatterns, "ShufflePatterns", "should the patterns be shuffled?");  
    DeclareOptionRef(fLinOut, "LinOut", "sets the activation function of the output units to linear or logistic");

}

//_______________________________________________________________________
void MethodRSNNS::ProcessOptions()
{

}

//_______________________________________________________________________
void MethodRSNNS::TestClassification()
{
    Log()<<kINFO<<"Testing Classification RSNNS METHOD  "<<Endl;
    gSystem->MakeDirectory("RSNNS");
    gSystem->MakeDirectory("RSNNS/plots");

    MethodBase::TestClassification();
}


//_______________________________________________________________________
Double_t MethodRSNNS::GetMvaValue( Double_t* errLower, Double_t* errUpper)
{
    Double_t mvaValue;
    if(Data()->GetCurrentType()==Types::kTraining) 
    {
        if(fProbResultForTrainSig.size()==0)
        {
           r<<"RMVA.RSNNS.Predictor.Train.Prob<-predict(RMVA.RSNNS.Model,RMVA.RSNNS.fDfTrain,type='prob')";
           r["as.vector(RMVA.RSNNS.Predictor.Train.Prob[,1])"]>>fProbResultForTrainSig;

        }
        mvaValue=fProbResultForTrainSig[fMvaCounter];
       if(fMvaCounter < Data()->GetNTrainingEvents()-1) fMvaCounter++;
       else fMvaCounter=0;
    }else
    {
        if(fProbResultForTestSig.size()==0)
        {
        r<<"RMVA.RSNNS.Predictor.Test.Prob<-predict(RMVA.RSNNS.Model,RMVA.RSNNS.fDfTest,type='prob')";
        r["as.vector(RMVA.RSNNS.Predictor.Test.Prob[,1])"]>>fProbResultForTestSig;
        }
        
        mvaValue=fProbResultForTestSig[fMvaCounter];
        
       if(fMvaCounter < Data()->GetNTestEvents()-1) fMvaCounter++;
       else fMvaCounter=0;
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

