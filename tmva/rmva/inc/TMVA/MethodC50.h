// @(#)root/rmva $Id$
// Author: Omar Zapata 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : RMethodC50                                                            *
 *                                                                                *
 * Description:                                                                   *
 *      R´s Package C50  method based on ROOTR                                    *
 *                                                                                *
 **********************************************************************************/

#ifndef ROOT_TMVA_RMethodC50
#define ROOT_TMVA_RMethodC50

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// RMethodC50                                                          //
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
   class MethodC50 : public RMethodBase {

   public :

      // constructors
      MethodC50( const TString& jobName,
                   const TString& methodTitle,
                   DataSetInfo& theData,
                   const TString& theOption = "",
                   TDirectory* theTargetDir = NULL );

      MethodC50( DataSetInfo& dsi,
                   const TString& theWeightFile,
                   TDirectory* theTargetDir = NULL );


      virtual ~MethodC50( void );
      virtual void     Train();
      // options treatment
      virtual void     Init();
      virtual void     DeclareOptions();
      virtual void     ProcessOptions();
      // create ranking
      virtual const Ranking* CreateRanking(){return NULL;}// = 0;
      
      virtual Double_t GetMvaValue( Double_t* errLower = 0, Double_t* errUpper = 0) {return 0;}//= 0;

     Bool_t HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets );
     

      
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
       UInt_t fNTrials;//number of trials with boost enabled
       Bool_t fRules;//A logical: should the tree be decomposed into a rule-based model?

      // get help message text
      void GetHelpMessage() const;

      ClassDef(MethodC50,0)
   };
} // namespace TMVA
#endif
