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

       TString fNetType;//default RMPL
       //RSNNS Options for all NN methods
       TVectorF  fSize;//number of units in the hidden layer(s)
       UInt_t fMaxit;//maximum of iterations to learn
       
       TString fInitFunc;//the initialization function to use
       Float_t *fInitFuncParams;//the parameters for the initialization function (type 6 see getSnnsRFunctionTable() in RSNNS package)
       
       TString fLearnFunc;//the learning function to use
       Float_t *fLearnFuncParams;//the parameters for the learning function
       
       TString fUpdateFunc;//the update function to use
       TVectorF fUpdateFuncParams;//the parameters for the update function

       TString fHiddenActFunc;//the activation function of all hidden units
       Bool_t fShufflePatterns;//should the patterns be shuffled?
       Bool_t fLinOut;//sets the activation function of the output units to linear or logistic
       
       //The next function will not implemented yet because is not well documented.
       TString fPruneFunc;//the pruning function to use
       TVectorF fPruneFuncParams;//the parameters for the pruning function. Unlike the
                                 //other functions, these have to be given in a named list. See
                                 //the pruning demos for further explanation.
       
      // get help message text
      void GetHelpMessage() const;

      ClassDef(MethodRSNNS,0)
   };
} // namespace TMVA
#endif
