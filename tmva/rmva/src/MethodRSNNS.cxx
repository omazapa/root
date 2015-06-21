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
    RMethodBase( jobName, Types::kRSNNS, methodTitle, dsi, theOption, theTargetDir )
{
    // standard constructor for the RSNNS

}

//_______________________________________________________________________
MethodRSNNS::MethodRSNNS( DataSetInfo& theData, const TString& theWeightFile, TDirectory* theTargetDir )
    : RMethodBase( Types::kRSNNS, theData, theWeightFile, theTargetDir )
{

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
    //Paassing Data to R's environment
    //NOTE:need improved names in R's environment using JobName of TMVA
    r["RMVA.RSNNS.fDfTrain"]=fDfTrain;
    r["RMVA.RSNNS.fWeightTrain"]=fWeightTrain;
    r<<"write.table(RMVA.RSNNS.fDfTrain,file='fDfTrain.txt')";

    r["RMVA.RSNNS.fDfTest"]=fDfTest;
    r["RMVA.RSNNS.fWeightTest"]=fWeightTest;
    r<<"write.table(RMVA.RSNNS.fDfTest,file='fDfTest.txt')";

    //factors creations
    r["RMVA.RSNNS.fFactorTrain"]=fFactorTrain;
    r<<"RMVA.RSNNS.fFactorTrain<-factor(RMVA.RSNNS.fFactorTrain)";
    r["RMVA.RSNNS.fFactorTest"]=fFactorTest;
    r<<"RMVA.RSNNS.fFactorTest<-factor(RMVA.RSNNS.fFactorTest)";

    //Spectator creation
    r["RMVA.RSNNS.fDfSpectators"]=fDfSpectators;

    r["RMVA.RSNNS.fCounter"]=0;
}

void MethodRSNNS::Train()
{
    if (Data()->GetNTrainingEvents()==0) Log() << kFATAL << "<Train> Data() has zero events" << Endl;

//    r.SetVerbose(1);
   //    r.SetVerbose(0);
}

//_______________________________________________________________________
void MethodRSNNS::DeclareOptions()
{
    //

}

//_______________________________________________________________________
void MethodRSNNS::ProcessOptions()
{

}

//_______________________________________________________________________
void MethodRSNNS::TestClassification()
{
    Log()<<kINFO<<"Testing Classification RSNNS METHOD  "<<Endl;
    
//    r.SetVerbose(1);
//    r.SetVerbose(0);
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
       return mvaValue;
    }else
    {
       return mvaValue;
    }
}


Double_t MethodRSNNS::GetMvaValue( const TMVA::Event* const ev, Double_t* errLower, Double_t* errUpper )
{
         // cannot determine error
         NoErrorCalc(errLower, errUpper);
         const UInt_t nvar = DataInfo().GetNVariables();
         ROOT::R::TRDataFrame fDfEvent;
         for(UInt_t i=0;i<nvar;i++)
     {
         fDfEvent[GetInputLabel( i ).Data()]=ev->GetValues()[i];
     }
         r["RMVA.RSNNS.fDfEvent"]<<fDfEvent;
         //   r<<"print(RMVA.RSNNS.Event)";

         TString type;
         if(type=="signal") return Types::kSignal;
         else return Types::kBackground;

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

