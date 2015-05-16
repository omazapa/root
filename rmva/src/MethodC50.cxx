// @(#)root/rmva $Id$
// Author: Omar Zapata 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodC50                                                              *
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

#include "TMVA/VariableTransformBase.h"
#include "TMVA/MethodC50.h"
#include "TMVA/Tools.h"
#include "TMVA/Ranking.h"
#include "TMVA/Types.h"
#include "TMVA/PDF.h"
#include "TMVA/ClassifierFactory.h"

#ifndef ROOT_R_TRInterface
#include<TRInterface.h>
#endif


using namespace TMVA;
REGISTER_METHOD(C50)


ClassImp(MethodC50)


//_______________________________________________________________________
MethodC50::MethodC50( const TString& jobName,
                          const TString& methodTitle,
                          DataSetInfo& dsi,
                          const TString& theOption,
                          TDirectory* theTargetDir ) :
   RMethodBase( jobName, Types::kC50, methodTitle, dsi, theOption, theTargetDir )
{
   // standard constructor for the C50
        
}

//_______________________________________________________________________
MethodC50::MethodC50( DataSetInfo& theData, const TString& theWeightFile, TDirectory* theTargetDir )
   : RMethodBase( Types::kC50, theData, theWeightFile, theTargetDir )
{
   // constructor from weight file
}


//_______________________________________________________________________
MethodC50::~MethodC50( void )
{
}

//_______________________________________________________________________
Bool_t MethodC50::HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t /*numberTargets*/ )
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   return kFALSE;
}


//_______________________________________________________________________
void     MethodC50::Init()
{
    if(!r.IsInstalled("C50"))
    {
        Error("Init","R's package C50 is not installed.");
        return;
    }
    
    if(!r.Require("C50"))
    {
        Error("Init","R's package C50 can not be loaded.");
        return;
    }
    const UInt_t nvar = DataInfo().GetNVariables();
    
    const UInt_t ntrains=Data()->GetNEvtSigTest()+Data()->GetNEvtBkgdTest();
    
    std::vector<std::vector<Float_t> > fArrayTrain(nvar);
    
    for(UInt_t j=0;j<ntrains;j++)
    {  
        for(UInt_t i=0;i<nvar;i++)
        {  
            fArrayTrain[i].push_back( Data()->GetEvent(  j, Types::ETreeType::kTraining )->GetValues()[i]);
        }    
        
    }    
    for(UInt_t i=0;i<nvar;i++)
    {
        fDfTrain[GetInputLabel( i ).Data()]=fArrayTrain[i];
    }
    
    TString jobName=GetJobName();
    
//    r["fDfTrain"]=fDfTrain;
//    r<<"print(fDfTrain)";

//    const UInt_t trainsize = Data()->GetNTrainingEvents();
//    const UInt_t testsize = Data()->GetNTestEvents();
//    
//    Log() << Endl;
//    Log() << gTools().Color("bold") << "--- Nvars:"<<nvar << gTools().Color("reset") << Endl;
//    Log() << gTools().Color("bold") << "--- Train Size:"<<trainsize << gTools().Color("reset") << Endl;
//    Log() << gTools().Color("bold") << "--- Test Size:"<<testsize << gTools().Color("reset") << Endl;
//    
//    Log() << gTools().Color("bold") << "--- NEvtSigTest:"<<Data()->GetNEvtSigTest() << gTools().Color("reset") << Endl;
//    Log() << gTools().Color("bold") << "--- NEvtBkgdTest:"<<Data()->GetNEvtBkgdTest() << gTools().Color("reset") << Endl;
//    Log() << gTools().Color("bold") << "--- NEvtSigTrain:"<<Data()->GetNEvtSigTrain() << gTools().Color("reset") << Endl;
//    Log() << gTools().Color("bold") << "--- NEvtBkgdTrain:"<<Data()->GetNEvtBkgdTrain() << gTools().Color("reset") << Endl;
//    
//    
  
//    ROOT::R::TRDataFrame df_testing;
//    for(int i=0;i<nvar;i++)
//    {
//    Log() << Endl;
//    Log() << gTools().Color("bold") << "--- InputVar:"<< GetInputVar( i ) << gTools().Color("reset") << Endl;
//    Log() << gTools().Color("bold") << "--- InputLabel:"<< GetInputLabel( i ) << gTools().Color("reset") << Endl;
//    Log() << gTools().Color("bold") << "--- GetInputTitle:"<< GetInputTitle( i ) << gTools().Color("reset") << Endl;
//    Log() << gTools().Color("bold") << "--- Size Training:"<< Data()->GetEvent(  i, Types::ETreeType::kTraining )->GetValues().size() << gTools().Color("reset") << Endl;
//    Log() << gTools().Color("bold") << "--- Size Testing:"<< Data()->GetEvent(  i, Types::ETreeType::kTesting )->GetValues().size() << gTools().Color("reset") << Endl;
    //df_testing[ GetInputLabel( i ).Data()] =  Data()->GetEvent(  i, Types::ETreeType::kTraining )->GetValues();  
//    df_testing[ GetInputLabel( i ).Data()] =  Data()->GetEvent(  i, Types::ETreeType::kTraining )->GetValues();  
//    }
//    r["df_testing"]=df_testing;
//    r<<"print(df_testing)";
//   
    
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

