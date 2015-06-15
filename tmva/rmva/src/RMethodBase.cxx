// @(#)root/tmva $Id$
// Author: Omar Zapata 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : RMethodBase                                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Virtual base class for all MVA method based on ROOTR                      *
 *                                                                                *
 **********************************************************************************/

#include<TMVA/RMethodBase.h>
#include<TApplication.h>
using namespace TMVA;

ClassImp(RMethodBase)
        
//_______________________________________________________________________
RMethodBase::RMethodBase( const TString& jobName,
                  Types::EMVA methodType,
                  const TString& methodTitle,
                  DataSetInfo& dsi,
                  const TString& theOption ,
                  TDirectory* theBaseDir ,ROOT::R::TRInterface &_r):MethodBase(jobName,methodType, methodTitle,dsi,theOption,theBaseDir),r(_r)
{
    LoadData();
    InitWrap();
}

//_______________________________________________________________________
RMethodBase::RMethodBase( Types::EMVA methodType,
                  DataSetInfo& dsi,
                  const TString& weightFile,
                  TDirectory* theBaseDir,ROOT::R::TRInterface &_r ):MethodBase(methodType,dsi,weightFile,theBaseDir),r(_r)
{
    LoadData();
    InitWrap();
}

//_______________________________________________________________________
void RMethodBase::InitWrap()
{
      gApplication->ProcessLine("#include<TMVA/EventWrapper.h>");
      gApplication->ProcessLine("LOAD_ROOTR_MODULE(RMVA_EventWrapper);");
      r.Execute("Event       <- .GlobalEnv$.__C__Rcpp_Event");
}

//_______________________________________________________________________
void RMethodBase::LoadData()
{
    ///////////////////////////
    //Loading Training Data  //
    ///////////////////////////
   const UInt_t nvar = DataInfo().GetNVariables();
    
    const UInt_t ntrains=Data()->GetNEvtBkgdTrain()+Data()->GetNEvtSigTrain();
    
    //array of columns for every var to create a dataframe for training
    std::vector<std::vector<Float_t> > fArrayTrain(nvar);
    fWeightTrain.ResizeTo(ntrains);
    for(UInt_t j=0;j<ntrains;j++)
    {  
        const Event *ev=Data()->GetEvent(  j, Types::ETreeType::kTraining );
        //creating array with class type(signal or background) for factor required
        if(ev->GetClass()==fSignalClass) fFactorTrain.push_back("signal");
        else fFactorTrain.push_back("background");
        
        fWeightTrain[j]=ev->GetOriginalWeight();
        
        //filling vector of columns for training
        for(UInt_t i=0;i<nvar;i++)
        {  
            fArrayTrain[i].push_back(ev->GetValues()[i]);
        }    
        
    }    
    for(UInt_t i=0;i<nvar;i++)
    {
        fDfTrain[GetInputLabel( i ).Data()]=fArrayTrain[i];
    }
 
    ////////////////////////
    //Loading Test  Data  //
    ////////////////////////
    
    const UInt_t ntests = Data()->GetNEvtSigTest()+Data()->GetNEvtBkgdTest();
    const UInt_t nspectators = DataInfo().GetNSpectators(kTRUE);

    //array of columns for every var to create a dataframe for testing
    std::vector<std::vector<Float_t> > fArrayTest(nvar);
    //array of columns for every spectator to create a dataframe for testing
    std::vector<std::vector<Float_t> > fArraySpectators(nvar);
    fWeightTest.ResizeTo(ntests);
    for(UInt_t j=0;j<ntests;j++)
    {  
        const Event *ev=Data()->GetEvent(  j, Types::ETreeType::kTesting );
        //creating array with class type(signal or background) for factor required
        if(ev->GetClass()==fSignalClass) fFactorTest.push_back("signal");
        else fFactorTest.push_back("background");
        
        fWeightTest[j]=ev->GetOriginalWeight();
        
        for(UInt_t i=0;i<nvar;i++)
        {  
            fArrayTest[i].push_back( ev->GetValues()[i]);
        }    
        for(UInt_t i=0;i<nspectators;i++)
        {  
            fArraySpectators[i].push_back( ev->GetSpectator(i));
        }    
    }    
    for(UInt_t i=0;i<nvar;i++)
    {
        fDfTest[GetInputLabel( i ).Data()]=fArrayTest[i];
    }
    for(UInt_t i=0;i<nspectators;i++)
    {
        fDfSpectators[DataInfo().GetSpectatorInfo(i).GetLabel().Data()]=fArraySpectators[i];
//        Log() <<kERROR<< " Spectator Label "<<i<<" "<<DataInfo().GetSpectatorInfo(i).GetLabel()
//        << Endl;    
//        Log() <<kERROR<< " Spectator "<<i<<" "<<DataInfo().GetSpectatorInfo(i).GetExpression()
//        << Endl;    
    }

}