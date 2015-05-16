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
using namespace TMVA;

ClassImp(RMethodBase)
        
// default constructur
RMethodBase::RMethodBase( const TString& jobName,
                  Types::EMVA methodType,
                  const TString& methodTitle,
                  DataSetInfo& dsi,
                  const TString& theOption ,
                  TDirectory* theBaseDir ,ROOT::R::TRInterface &_r):MethodBase(jobName,methodType, methodTitle,dsi,theOption,theBaseDir),r(_r)
{

}

// constructor used for Testing + Application of the MVA, only (no training),
// using given weight file
RMethodBase::RMethodBase( Types::EMVA methodType,
                  DataSetInfo& dsi,
                  const TString& weightFile,
                  TDirectory* theBaseDir,ROOT::R::TRInterface &_r ):MethodBase(methodType,dsi,weightFile,theBaseDir),r(_r)
{
}
