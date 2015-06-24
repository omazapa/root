// @(#)root/rmva $Id$
// Author: Omar Zapata 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : RMethodRSVM                                                            *
 *                                                                                *
 * Description:                                                                   *
 *      R´s Package RSVM  method based on ROOTR                                    *
 *                                                                                *
 **********************************************************************************/

#ifndef ROOT_TMVA_RMethodRSVM
#define ROOT_TMVA_RMethodRSVM

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// RMethodRSVM                                                          //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TMVA_RMethodBase
#include "TMVA/RMethodBase.h"
#endif

namespace TMVA {

   class Factory;  // DSMTEST
   class Reader;   // DSMTEST
   class DataSetManager;  // DSMTEST
   class Types;
   class MethodRSVM : public RMethodBase {

   public :

      // constructors
      MethodRSVM( const TString& jobName,
                   const TString& methodTitle,
                   DataSetInfo& theData,
                   const TString& theOption = "",
                   TDirectory* theTargetDir = NULL );

      MethodRSVM( DataSetInfo& dsi,
                   const TString& theWeightFile,
                   TDirectory* theTargetDir = NULL );


      ~MethodRSVM( void );
      void     Train();
      // options treatment
      void     Init();
      void     DeclareOptions();
      void     ProcessOptions();
      // create ranking
      const Ranking* CreateRanking(){return NULL;}// = 0;
      
      
      Bool_t HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets );
     
      // performs classifier testing
      virtual void     TestClassification();
      
      
      Double_t GetMvaValue( Double_t* errLower = 0, Double_t* errUpper = 0);
      
      using MethodBase::ReadWeightsFromStream;
      // the actual "weights"
      virtual void AddWeightsXMLTo      ( void* parent ) const {}// = 0;
      virtual void ReadWeightsFromXML   ( void* wghtnode ){}// = 0;
      virtual void ReadWeightsFromStream( std::istream& ) {}//= 0;       // backward compatibility

   private :
      DataSetManager*    fDataSetManager;     // DSMTEST
      friend class Factory;                   // DSMTEST
      friend class Reader;                    // DSMTEST      
   protected:
       
       UInt_t fMvaCounter;       
       std::vector<Float_t> fProbResultForTrainSig;
       std::vector<Float_t> fProbResultForTestSig;
      // get help message text
      void GetHelpMessage() const;

      ClassDef(MethodRSVM,0)
   };
} // namespace TMVA
#endif
