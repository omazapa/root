// @(#)root/rmva $Id$
// Author: Omar Zapata 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : RMethodRSNNS                                                          *
 *                                                                                *
 * Description:                                                                   *
 *      R´s Package RSNNS  method based on ROOTR                                  *
 *                                                                                *
 **********************************************************************************/

#ifndef ROOT_TMVA_RMethodRSNNS
#define ROOT_TMVA_RMethodRSNNS

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// RMethodRSNNS                                                         //
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
   class MethodRSNNS : public RMethodBase {

   public :

      // constructors
      MethodRSNNS( const TString& jobName,
                   const TString& methodTitle,
                   DataSetInfo& theData,
                   const TString& theOption = "",
                   TDirectory* theTargetDir = NULL );

      MethodRSNNS( DataSetInfo& dsi,
                   const TString& theWeightFile,
                   TDirectory* theTargetDir = NULL );


      ~MethodRSNNS( void );
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
      // signal/background classification response
      Double_t GetMvaValue( const TMVA::Event* const ev, Double_t* err = 0, Double_t* errUpper = 0 );
      
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
       
      // get help message text
      void GetHelpMessage() const;

      ClassDef(MethodRSNNS,0)
   };
} // namespace TMVA
#endif
