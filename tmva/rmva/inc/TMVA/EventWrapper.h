// @(#)root/tmva/rmva $Id$
// Author: Omar Zapata 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : EventWrapper                                                          *
 *                                                                                *
 * Description:                                                                   *
 *      TMVA::Event wrapper to support pass objects to R's environment            *
 *                                                                                *
 **********************************************************************************/
#ifndef ROOT_TMVA_EventWrapper
#define ROOT_TMVA_EventWrapper

#ifndef ROOT_TMVA_Event
#include<TMVA/Event.h>
#endif

#ifndef ROOT_R_RExports
#include<RExports.h>
#endif

//________________________________________________________________________________________________________
/**
   This is TMVA::Event's wrapper class for R


   @ingroup RMVA
*/
namespace ROOT {
   namespace R {
      class EventWrapper: public TMVA::Event {
      public:
      // constructors
      EventWrapper():TMVA::Event(){}
      EventWrapper( const TMVA::Event& e):TMVA::Event(e){}
      
      // accessors
      Bool_t  IsDynamic()         const {return TMVA::Event::IsDynamic(); }

      //      Double_t GetWeight()         const { return fWeight*fBoostWeight; }
      Double_t GetWeight()         const{return TMVA::Event::GetWeight();}
      Double_t GetOriginalWeight() const { return TMVA::Event::GetOriginalWeight(); }
      Double_t GetBoostWeight()    const { return TMVA::Event::GetBoostWeight(); }
      UInt_t   GetClass()          const { return TMVA::Event::GetClass(); }  

      UInt_t   GetNVariables()        const{return TMVA::Event::GetNVariables();};
      UInt_t   GetNTargets()          const{return TMVA::Event::GetNTargets();}
      UInt_t   GetNSpectators()       const{return TMVA::Event::GetNSpectators();}

      Float_t  GetValue( UInt_t ivar) const{return TMVA::Event::GetValue(ivar);}
      std::vector<Float_t>& GetValues() { return TMVA::Event::GetValues();}
      void     Print() {TMVA::Event::Print(std::cout);}
//      const std::vector<Float_t>& GetValues() const{return ::GetValues();}

      Float_t  GetTarget( UInt_t itgt ) const { return TMVA::Event::GetTarget(itgt); }
      std::vector<Float_t>& GetTargets()  { return TMVA::Event::GetTargets(); }
//      const std::vector<Float_t> GetTargets() const { return ::GetTargets(); }

//      Float_t  GetSpectator( UInt_t ivar) const{return ::GetSpectator(ivar);}
//      std::vector<Float_t>& GetSpectators()  { return ::GetSpectators(); }
//      const std::vector<Float_t>& GetSpectators() const { return ::GetSpectators(); }
      };
   }
}

namespace Rcpp {
   template<> SEXP wrap(const TMVA::Event &e)
   {
      return Rcpp::wrap(ROOT::R::EventWrapper(e));
   }
   template<> ROOT::R::EventWrapper as(SEXP e)
   {
      return Rcpp::as<ROOT::R::EventWrapper>(e);
   }
}

ROOTR_EXPOSED_CLASS_INTERNAL(EventWrapper)



ROOTR_MODULE(RMVA_EventWrapper)
{

   ROOT::R::class_<ROOT::R::EventWrapper>("Event", "Event container")
   .constructor()
   .method("IsDynamic", (Bool_t (ROOT::R::EventWrapper::*)())&ROOT::R::EventWrapper::IsDynamic)
   .method("GetWeight", (Double_t (ROOT::R::EventWrapper::*)())&ROOT::R::EventWrapper::GetWeight)
   .method("GetClass", (UInt_t (ROOT::R::EventWrapper::*)())&ROOT::R::EventWrapper::GetClass)
   .method("GetNVariables", (UInt_t (ROOT::R::EventWrapper::*)())&ROOT::R::EventWrapper::GetNVariables)
   .method("GetNTargets", (UInt_t (ROOT::R::EventWrapper::*)())&ROOT::R::EventWrapper::GetNTargets)
   .method("GetNSpectators", (UInt_t (ROOT::R::EventWrapper::*)())&ROOT::R::EventWrapper::GetNSpectators)
   .method("GetValue", (Float_t (ROOT::R::EventWrapper::*)(UInt_t))&ROOT::R::EventWrapper::GetValue)
   .method("GetValues", (std::vector<Float_t>& (ROOT::R::EventWrapper::*)())&ROOT::R::EventWrapper::GetValues)
   .method("GetTarget", (Float_t (ROOT::R::EventWrapper::*)(UInt_t))&ROOT::R::EventWrapper::GetTarget)
   .method("GetTargets", (std::vector<Float_t>& (ROOT::R::EventWrapper::*)())&ROOT::R::EventWrapper::GetTargets)
   .method("Print", (void (ROOT::R::EventWrapper::*)())&ROOT::R::EventWrapper::Print)
   ;
}
#endif
