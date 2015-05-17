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
#include "TVectorD.h"

#include "TMVA/VariableTransformBase.h"
#include "TMVA/MethodC50.h"
#include "TMVA/Tools.h"
#include "TMVA/Ranking.h"
#include "TMVA/Types.h"
#include "TMVA/PDF.h"
#include "TMVA/ClassifierFactory.h"



using namespace TMVA;
REGISTER_METHOD(C50)


ClassImp(MethodC50)


//_______________________________________________________________________
MethodC50::MethodC50( const TString& jobName,
                          const TString& methodTitle,
                          DataSetInfo& dsi,
                          const TString& theOption,
                          TDirectory* theTargetDir ) :
   RMethodBase( jobName, Types::kC50, methodTitle, dsi, theOption, theTargetDir ),fNTrials(1),fRules(kFALSE)
{
   // standard constructor for the C50
        
}

//_______________________________________________________________________
MethodC50::MethodC50( DataSetInfo& theData, const TString& theWeightFile, TDirectory* theTargetDir )
   : RMethodBase( Types::kC50, theData, theWeightFile, theTargetDir ),fNTrials(1),fRules(kFALSE)
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
        Error( "Init","R's package C50 is not installed.");
        Log() << kERROR << " R's package C50 is not installed."
              << Endl;    
        return;
    }
        
    if(!r.Require("C50"))
    {
        Error("Init","R's package C50 can not be loaded.");
        Log() << kERROR << " R's package C50 can not be loaded."
              << Endl;
        return;
    }
    const UInt_t nvar = DataInfo().GetNVariables();
    
    const UInt_t ntrains=Data()->GetNEvtBkgdTrain()+Data()->GetNEvtSigTrain();
    
    //array of columns for every var to create a dataframe for training
    std::vector<std::vector<Float_t> > fArrayTrain(nvar);
    
    for(UInt_t j=0;j<ntrains;j++)
    {  
        //creating array with class type(signal or background) for factor required
        if(Data()->GetEvent(  j, Types::ETreeType::kTraining )->GetClass()==fSignalClass) fFactorTrain.push_back("signal");
        else fFactorTrain.push_back("background");
        
        //filling vector of columns for training
        for(UInt_t i=0;i<nvar;i++)
        {  
            fArrayTrain[i].push_back( Data()->GetEvent(  j, Types::ETreeType::kTraining )->GetValues()[i]);
        }    
        
    }    
    for(UInt_t i=0;i<nvar;i++)
    {
        fDfTrain[GetInputLabel( i ).Data()]=fArrayTrain[i];
    }
    
    //NOTE:need improved names in R's environment using JobName of TMVA
    r["RMVA.C50.fDfTrain"]=fDfTrain;
    r<<"write.table(RMVA.C50.fDfTrain,file='fDfTrain.txt')";

    
    const UInt_t ntests = Data()->GetNEvtSigTest()+Data()->GetNEvtBkgdTest();

    //array of columns for every var to create a dataframe for testing
    std::vector<std::vector<Float_t> > fArrayTest(nvar);
    
    for(UInt_t j=0;j<ntests;j++)
    {  
        //creating array with class type(signal or background) for factor required
        if(Data()->GetEvent(  j, Types::ETreeType::kTesting )->GetClass()==fSignalClass) fFactorTest.push_back("signal");
        else fFactorTest.push_back("background");
        
        for(UInt_t i=0;i<nvar;i++)
        {  
            fArrayTest[i].push_back( Data()->GetEvent(  j, Types::ETreeType::kTesting )->GetValues()[i]);
        }    
        
    }    
    for(UInt_t i=0;i<nvar;i++)
    {
        fDfTest[GetInputLabel( i ).Data()]=fArrayTest[i];
    }
    r["RMVA.C50.fDfTest"]=fDfTest;
    r<<"write.table(RMVA.C50.fDfTest,file='fDfTest.txt')";
   
    //factors creations
    r["RMVA.C50.FactorTrain"]=fFactorTrain;
    r<<"RMVA.C50.FactorTrain<-factor(RMVA.C50.FactorTrain)";
    
}

void MethodC50::Train()
{
    r<<"RMVA.C50.Model<-C5.0(RMVA.C50.fDfTrain,RMVA.C50.FactorTrain,RMVA.C50.NTrials,rules=RMVA.C50.Rules)";
    r.SetVerbose(1);
    r<<"summary(RMVA.C50.Model)";
    r.SetVerbose(0);
}

//_______________________________________________________________________
void MethodC50::DeclareOptions()
{
    //
    DeclareOptionRef(fNTrials, "NTrials", "An integer specifying the number of boosting iterations");
    DeclareOptionRef(fRules, "Rules", "A logical: should the tree be decomposed into a rule-basedmodel?");
}

//_______________________________________________________________________
void MethodC50::ProcessOptions()
{
    if (fNTrials<=0){
      Log() << kERROR << " fNTrials <=0... that does not work !! "
            << " I set it to 1 .. just so that the program does not crash"
            << Endl;
      fNTrials = 1;
   }
    r["RMVA.C50.NTrials"]=fNTrials;
    Log()<<"NTrials  "<<fNTrials<<Endl;
    r["RMVA.C50.Rules"]=fRules;
    Log()<<"Rules  "<<fRules<<Endl;
 
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

