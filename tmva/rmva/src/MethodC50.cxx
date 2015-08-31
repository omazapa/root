// @(#)root/rmva $Id$
// Author: Omar Zapata 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodC50                                                             *
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
#include "TMVA/Config.h"
#include "TMVA/Ranking.h"
#include "TMVA/Types.h"
#include "TMVA/PDF.h"
#include "TMVA/ClassifierFactory.h"

#include "TMVA/Results.h"

using namespace TMVA;

REGISTER_METHOD(C50)

ClassImp(MethodC50)

//creating an Instance
Bool_t MethodC50::IsModuleLoaded=ROOT::R::TRInterface::Instance().Require("C50");

//_______________________________________________________________________
MethodC50::MethodC50(const TString &jobName,
                     const TString &methodTitle,
                     DataSetInfo &dsi,
                     const TString &theOption,
                     TDirectory *theTargetDir) :RMethodBase(jobName, Types::kC50, methodTitle, dsi, theOption, theTargetDir),
   fNTrials(1), 
   fRules(kFALSE), 
   fMvaCounter(0),
   predict("predict.C5.0"),
   C50("C5.0"),
   C50Control("C5.0Control"),
   asfactor("as.factor"),
   fModel(NULL)   
{
   // standard constructor for the C50

   //C5.0Control options
   fControlSubset = kTRUE;
   fControlBands = 0;
   fControlWinnow = kFALSE;
   fControlNoGlobalPruning = kFALSE;
   fControlCF = 0.25;
   fControlMinCases = 2;
   fControlFuzzyThreshold = kFALSE;
   fControlSample = 0;
   r["sample.int(4096, size = 1) - 1L"] >> fControlSeed;
   fControlEarlyStopping = kTRUE;

   ListOfVariables=DataInfo().GetListOfVariables();
// default extension for weight files
   SetWeightFileDir( gConfig().GetIONames().fWeightFileDir );
}

//_______________________________________________________________________
MethodC50::MethodC50(DataSetInfo &theData, const TString &theWeightFile, TDirectory *theTargetDir)
   : RMethodBase(Types::kC50, theData, theWeightFile, theTargetDir),
   fNTrials(1), 
   fRules(kFALSE), 
   fMvaCounter(0),
   predict("predict.C5.0"),
   C50("C5.0"),
   C50Control("C5.0Control"),
   asfactor("as.factor"),
   fModel(NULL)   
{

   // constructor from weight file
   fControlSubset = kTRUE;
   fControlBands = 0;
   fControlWinnow = kFALSE;
   fControlNoGlobalPruning = kFALSE;
   fControlCF = 0.25;
   fControlMinCases = 2;
   fControlFuzzyThreshold = kFALSE;
   fControlSample = 0;
   r["sample.int(4096, size = 1) - 1L"] >> fControlSeed;
   fControlEarlyStopping = kTRUE;
// default extension for weight files
   SetWeightFileDir( gConfig().GetIONames().fWeightFileDir );
}


//_______________________________________________________________________
MethodC50::~MethodC50(void)
{
    if(fModel) delete fModel;
}

//_______________________________________________________________________
Bool_t MethodC50::HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets)
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   return kFALSE;
}


//_______________________________________________________________________
void     MethodC50::Init()
{

   if (!IsModuleLoaded) {
      Error("Init", "R's package C50 can not be loaded.");
      Log() << kFATAL << " R's package C50 can not be loaded."
            << Endl;
      return;
   }
}

void MethodC50::Train()
{
   if (Data()->GetNTrainingEvents() == 0) Log() << kFATAL << "<Train> Data() has zero events" << Endl;
   SEXP Model=C50(ROOT::R::Label["x"]=fDfTrain, \
              ROOT::R::Label["y"]=asfactor(fFactorTrain), \
              ROOT::R::Label["trials"]=fNTrials, \
              ROOT::R::Label["rules"]=fRules, \
              ROOT::R::Label["weights"]=fWeightTrain, \
              ROOT::R::Label["control"]=fModelControl);
   fModel=new ROOT::R::TRObject(Model);
   TString path=GetWeightFileDir()+"/C50Model.RData";
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Saving State File In:" << gTools().Color("reset")<<path<< Endl;
   Log() << Endl;
   r["C50Model"]<<Model;
   r<<"save(C50Model,file='"+path+"')";
}

//_______________________________________________________________________
void MethodC50::DeclareOptions()
{
   //
   DeclareOptionRef(fNTrials, "NTrials", "An integer specifying the number of boosting iterations");
   DeclareOptionRef(fRules, "Rules", "A logical: should the tree be decomposed into a rule-basedmodel?");

   //C5.0Control Options
   DeclareOptionRef(fControlSubset, "ControlSubset", "A logical: should the model evaluate groups of discrete \
                                      predictors for splits? Note: the C5.0 command line version defaults this \
                                      parameter to ‘FALSE’, meaning no attempted gropings will be evaluated \
                                      during the tree growing stage.");
   DeclareOptionRef(fControlBands, "ControlBands", "An integer between 2 and 1000. If ‘TRUE’, the model orders \
                                     the rules by their affect on the error rate and groups the \
                                     rules into the specified number of bands. This modifies the \
                                     output so that the effect on the error rate can be seen for \
                                     the groups of rules within a band. If this options is \
                                     selected and ‘rules = kFALSE’, a warning is issued and ‘rules’ \
                                     is changed to ‘kTRUE’.");
   DeclareOptionRef(fControlWinnow, "ControlWinnow", "A logical: should predictor winnowing (i.e feature selection) be used?");
   DeclareOptionRef(fControlNoGlobalPruning, "ControlNoGlobalPruning", "A logical to toggle whether the final, global pruning \
                                                                         step to simplify the tree.");
   DeclareOptionRef(fControlCF, "ControlCF", "A number in (0, 1) for the confidence factor.");
   DeclareOptionRef(fControlMinCases, "ControlMinCases", "an integer for the smallest number of samples that must be \
                                                           put in at least two of the splits.");

   DeclareOptionRef(fControlFuzzyThreshold, "ControlFuzzyThreshold", "A logical toggle to evaluate possible advanced splits \
                                                                      of the data. See Quinlan (1993) for details and examples.");
   DeclareOptionRef(fControlSample, "ControlSample", "A value between (0, .999) that specifies the random \
                                                       proportion of the data should be used to train the model. By \
                                                       default, all the samples are used for model training. Samples \
                                                       not used for training are used to evaluate the accuracy of \
                                                       the model in the printed output.");
   DeclareOptionRef(fControlSeed, "ControlSeed", " An integer for the random number seed within the C code.");
   DeclareOptionRef(fControlEarlyStopping, "ControlEarlyStopping", " A logical to toggle whether the internal method for \
                                                                      stopping boosting should be used.");


}

//_______________________________________________________________________
void MethodC50::ProcessOptions()
{
   if (fNTrials <= 0) {
      Log() << kERROR << " fNTrials <=0... that does not work !! "
            << " I set it to 1 .. just so that the program does not crash"
            << Endl;
      fNTrials = 1;
   }
   fModelControl=C50Control(ROOT::R::Label["subset"]=fControlSubset, \
                            ROOT::R::Label["bands"]=fControlBands, \
                            ROOT::R::Label["winnow"]=fControlWinnow, \
                            ROOT::R::Label["noGlobalPruning"]=fControlNoGlobalPruning, \
                            ROOT::R::Label["CF"]=fControlCF, \
                            ROOT::R::Label["minCases"]=fControlMinCases, \
                            ROOT::R::Label["fuzzyThreshold"]=fControlFuzzyThreshold, \
                            ROOT::R::Label["sample"]=fControlSample, \
                            ROOT::R::Label["seed"]=fControlSeed, \
                            ROOT::R::Label["earlyStopping"]=fControlEarlyStopping);
}

//_______________________________________________________________________
void MethodC50::TestClassification()
{
   Log() << kINFO << "Testing Classification C50 METHOD  " << Endl;
   MethodBase::TestClassification();
}


//_______________________________________________________________________
Double_t MethodC50::GetMvaValue(Double_t *errLower, Double_t *errUpper)
{
   NoErrorCalc(errLower,errUpper);
   Double_t mvaValue;
   const TMVA::Event *ev=GetEvent();
   const UInt_t nvar = DataInfo().GetNVariables();
   ROOT::R::TRDataFrame fDfEvent;
   for(UInt_t i=0;i<nvar;i++)
   {
      fDfEvent[DataInfo().GetListOfVariables()[i].Data()]=ev->GetValues()[i];
   }
   //if using persistence model
   if(!fModel)
   {
       ReadStateFromFile();
   }
   TVectorD result=predict(*fModel,fDfEvent,ROOT::R::Label["type"]="prob");
   mvaValue=result[1];//returning signal prob
   return mvaValue;
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

//_______________________________________________________________________
void TMVA::MethodC50::ReadStateFromFile()
{
   ROOT::R::TRInterface::Instance().Require("C50");
   TString path=GetWeightFileDir()+"/C50Model.RData";
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Loading State File From:" << gTools().Color("reset")<<path<< Endl;
   Log() << Endl;
   r<<"load('"+path+"')"; 
   SEXP Model;
   r["C50Model"]>>Model;
   fModel=new ROOT::R::TRObject(Model);
   
}

//_______________________________________________________________________
void TMVA::MethodC50::MakeClass( const TString& theClassFileName ) const
{
   // the default consists of
   TString classFileName = "";
   if (theClassFileName == "")
      classFileName = GetWeightFileDir() + "/" + GetJobName() + "_" + GetMethodName() + ".class.C";
   else
      classFileName = theClassFileName;

   TString className = TString("Read") + GetMethodName();

   TString tfname( classFileName );
   Log() << kINFO << "Creating standalone response class: "
         << gTools().Color("lightblue") << classFileName << gTools().Color("reset") << Endl;

   std::ofstream fout( classFileName );
   if (!fout.good()) { // file could not be opened --> Error
      Log() << kFATAL << "<MakeClass> Unable to open file: " << classFileName << Endl;
   }

   // now create the class
   // preamble
   fout << "// Class: " << className << std::endl;
   fout << "// Automatically generated by MethodBase::MakeClass" << std::endl << "//" << std::endl;

   // print general information and configuration state
   fout << std::endl;
   fout << "/* configuration options =====================================================" << std::endl << std::endl;
   WriteStateToStream( fout );
   fout << std::endl;
   fout << "============================================================================ */" << std::endl;

   // generate the class
   fout << "" << std::endl;
   fout << "#include <vector>" << std::endl;
   fout << "#include <cmath>" << std::endl;
   fout << "#include <string>" << std::endl;
   fout << "#include <iostream>" << std::endl;
   fout << "#include <TRInterface.h>" << std::endl;
   fout << "" << std::endl;
   // now if the classifier needs to write some addicional classes for its response implementation
   // this code goes here: (at least the header declarations need to come before the main class
   this->MakeClassSpecificHeader( fout, className );

   fout << "#ifndef IClassifierReader__def" << std::endl;
   fout << "#define IClassifierReader__def" << std::endl;
   fout << std::endl;
   fout << "class IClassifierReader {" << std::endl;
   fout << std::endl;
   fout << " public:" << std::endl;
   fout << std::endl;
   fout << "   // constructor" << std::endl;
   fout << "   IClassifierReader() : fStatusIsClean( true ) {}" << std::endl;
   fout << "   virtual ~IClassifierReader() {}" << std::endl;
   fout << std::endl;
   fout << "   // return classifier response" << std::endl;
   fout << "   virtual double GetMvaValue( const std::vector<double>& inputValues ) const = 0;" << std::endl;
   fout << std::endl;
   fout << "   // returns classifier status" << std::endl;
   fout << "   bool IsStatusClean() const { return fStatusIsClean; }" << std::endl;
   fout << std::endl;
   fout << " protected:" << std::endl;
   fout << std::endl;
   fout << "   bool fStatusIsClean;" << std::endl;
   fout << "};" << std::endl;
   fout << std::endl;
   fout << "#endif" << std::endl;
   fout << std::endl;
   fout << "static Bool_t IsModuleLoaded=ROOT::R:TRInterface::Instance().Require(\"C50\");"<<std::endl;
   fout << "class " << className << " : public IClassifierReader {" << std::endl;
   fout << std::endl;
   fout << "ROOT::R::TRFunctionImport predict;"<<std::endl;
   fout << "ROOT::R::TRObject fModel;"<<std::endl;
   fout << " public:" << std::endl;
   fout << std::endl;
   fout << "   // constructor" << std::endl;
   fout << "   " << className << "( std::vector<std::string>& theInputVars ) " << std::endl;
   fout << "      : IClassifierReader()," << std::endl;
   fout << "        fClassName( \"" << className << "\" )," << std::endl;
   fout << "        fNvars( " << GetNvar() << " )," << std::endl;
   fout << "        fIsNormalised( " << (IsNormalised() ? "true" : "false") << " )," << std::endl;
   fout << "        predict( \"predict.C5.0\" )" << std::endl;
   fout << "   {      " << std::endl;
   fout << "   //ROOT R Stuff//" << std::endl;
   fout << "      ROOT::R:TRInterface &r=ROOT::R:TRInterface::Instance();      " << std::endl;
   fout << "      if (!IsModuleLoaded) {" << std::endl;
   fout << "      Error(\"Init\", \"R's package C50 can not be loaded.\");" << std::endl;
   fout << "      Log() << kFATAL << \" R's package C50 can not be loaded.\"" << std::endl;
   fout << "      << Endl;" << std::endl;
   fout << "      return;" << std::endl;
   fout << "      }   " << std::endl;
   fout << "      r<<\"load(\\\"C50Model.RData\\\")\";" << std::endl;
   fout << "      r[\"C50Model\"]>>fModel;" << std::endl;   
   fout << "   //ROOT R Stuff//" << std::endl;
   fout << "      // the training input variables" << std::endl;
   fout << "      const char* inputVars[] = { ";
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      fout << "\"" << GetOriginalVarName(ivar) << "\"";
      if (ivar<GetNvar()-1) fout << ", ";
   }
   fout << " };" << std::endl;
   fout << std::endl;
   fout << "      // sanity checks" << std::endl;
   fout << "      if (theInputVars.size() <= 0) {" << std::endl;
   fout << "         std::cout << \"Problem in class \\\"\" << fClassName << \"\\\": empty input vector\" << std::endl;" << std::endl;
   fout << "         fStatusIsClean = false;" << std::endl;
   fout << "      }" << std::endl;
   fout << std::endl;
   fout << "      if (theInputVars.size() != fNvars) {" << std::endl;
   fout << "         std::cout << \"Problem in class \\\"\" << fClassName << \"\\\": mismatch in number of input values: \"" << std::endl;
   fout << "                   << theInputVars.size() << \" != \" << fNvars << std::endl;" << std::endl;
   fout << "         fStatusIsClean = false;" << std::endl;
   fout << "      }" << std::endl;
   fout << std::endl;
   fout << "      // validate input variables" << std::endl;
   fout << "      for (size_t ivar = 0; ivar < theInputVars.size(); ivar++) {" << std::endl;
   fout << "         if (theInputVars[ivar] != inputVars[ivar]) {" << std::endl;
   fout << "            std::cout << \"Problem in class \\\"\" << fClassName << \"\\\": mismatch in input variable names\" << std::endl" << std::endl;
   fout << "                      << \" for variable [\" << ivar << \"]: \" << theInputVars[ivar].c_str() << \" != \" << inputVars[ivar] << std::endl;" << std::endl;
   fout << "            fStatusIsClean = false;" << std::endl;
   fout << "         }" << std::endl;
   fout << "      }" << std::endl;
   fout << std::endl;
   fout << "      // initialize min and max vectors (for normalisation)" << std::endl;
   for (UInt_t ivar = 0; ivar < GetNvar(); ivar++) {
      fout << "      fVmin[" << ivar << "] = " << std::setprecision(15) << GetXmin( ivar ) << ";" << std::endl;
      fout << "      fVmax[" << ivar << "] = " << std::setprecision(15) << GetXmax( ivar ) << ";" << std::endl;
   }
   fout << std::endl;
   fout << "      // initialize input variable types" << std::endl;
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      fout << "      fType[" << ivar << "] = \'" << DataInfo().GetVariableInfo(ivar).GetVarType() << "\';" << std::endl;
   }
   fout << std::endl;
   fout << "      // initialize constants" << std::endl;
   fout << "      Initialize();" << std::endl;
   fout << std::endl;
   if (GetTransformationHandler().GetTransformationList().GetSize() != 0) {
      fout << "      // initialize transformation" << std::endl;
      fout << "      InitTransform();" << std::endl;
   }
   fout << "   }" << std::endl;
   fout << std::endl;
   fout << "   // destructor" << std::endl;
   fout << "   virtual ~" << className << "() {" << std::endl;
   fout << "      Clear(); // method-specific" << std::endl;
   fout << "   }" << std::endl;
   fout << std::endl;
   fout << "   // the classifier response" << std::endl;
   fout << "   // \"inputValues\" is a vector of input values in the same order as the " << std::endl;
   fout << "   // variables given to the constructor" << std::endl;
   fout << "   double GetMvaValue( const std::vector<double>& inputValues ) const;" << std::endl;
   fout << std::endl;
   fout << " private:" << std::endl;
   fout << std::endl;
   fout << "   // method-specific destructor" << std::endl;
   fout << "   void Clear();" << std::endl;
   fout << std::endl;
   if (GetTransformationHandler().GetTransformationList().GetSize()!=0) {
      fout << "   // input variable transformation" << std::endl;
      GetTransformationHandler().MakeFunction(fout, className,1);
      fout << "   void InitTransform();" << std::endl;
      fout << "   void Transform( std::vector<double> & iv, int sigOrBgd ) const;" << std::endl;
      fout << std::endl;
   }
   fout << "   // common member variables" << std::endl;
   fout << "   const char* fClassName;" << std::endl;
   fout << std::endl;
   fout << "   const size_t fNvars;" << std::endl;
   fout << "   size_t GetNvar()           const { return fNvars; }" << std::endl;
   fout << "   char   GetType( int ivar ) const { return fType[ivar]; }" << std::endl;
   fout << std::endl;
   fout << "   // normalisation of input variables" << std::endl;
   fout << "   const bool fIsNormalised;" << std::endl;
   fout << "   bool IsNormalised() const { return fIsNormalised; }" << std::endl;
   fout << "   double fVmin[" << GetNvar() << "];" << std::endl;
   fout << "   double fVmax[" << GetNvar() << "];" << std::endl;
   fout << "   double NormVariable( double x, double xmin, double xmax ) const {" << std::endl;
   fout << "      // normalise to output range: [-1, 1]" << std::endl;
   fout << "      return 2*(x - xmin)/(xmax - xmin) - 1.0;" << std::endl;
   fout << "   }" << std::endl;
   fout << std::endl;
   fout << "   // type of input variable: 'F' or 'I'" << std::endl;
   fout << "   char   fType[" << GetNvar() << "];" << std::endl;
   fout << std::endl;
   fout << "   // initialize internal variables" << std::endl;
   fout << "   void Initialize(){}" << std::endl;

   fout << "   double GetMvaValue__( const std::vector<double>& inputValues ) const" << std::endl;
   fout << "   {" << std::endl;
   fout << "      // classifier response value" << std::endl;
   fout << "      double retval = 0;" << std::endl;
   fout << "      ROOT::R::TRDataFrame fDfEvent;"<< std::endl;
   fout << "      for(UInt_t i=0;i<fNvars;i++)"<< std::endl;
   fout << "      {"<< std::endl;
   fout << "          fDfEvent[inputVars[i]]=inputValues[i];"<< std::endl;
   fout << "      }"<< std::endl;
   fout << "      TVectorD result=predict(fModel,fDfEvent,ROOT::R::Label[\"type\"]=\"prob\");"<< std::endl;
   fout << "      retval=result[1];"<< std::endl;
   fout << std::endl;
   fout << "      return retval;" << std::endl;
   fout << "   }" << std::endl;
   fout << "" << std::endl;
   fout << "   // private members (method specific)" << std::endl;

   fout << "   double GetMvaValue( const std::vector<double>& inputValues ) const" << std::endl;
   fout << "   {" << std::endl;
   fout << "      // classifier response value" << std::endl;
   fout << "      double retval = 0;" << std::endl;
   fout << std::endl;
   fout << "      // classifier response, sanity check first" << std::endl;
   fout << "      if (!IsStatusClean()) {" << std::endl;
   fout << "         std::cout << \"Problem in class \\\"\" << fClassName << \"\\\": cannot return classifier response\"" << std::endl;
   fout << "                   << \" because status is dirty\" << std::endl;" << std::endl;
   fout << "         retval = 0;" << std::endl;
   fout << "      }" << std::endl;
   fout << "      else {" << std::endl;
   fout << "         if (IsNormalised()) {" << std::endl;
   fout << "            // normalise variables" << std::endl;
   fout << "            std::vector<double> iV;" << std::endl;
   fout << "            iV.reserve(inputValues.size());" << std::endl;
   fout << "            int ivar = 0;" << std::endl;
   fout << "            for (std::vector<double>::const_iterator varIt = inputValues.begin();" << std::endl;
   fout << "                 varIt != inputValues.end(); varIt++, ivar++) {" << std::endl;
   fout << "               iV.push_back(NormVariable( *varIt, fVmin[ivar], fVmax[ivar] ));" << std::endl;
   fout << "            }" << std::endl;
      // call the classifier specific output (the classifier must close the class !)
   MakeClassSpecific( fout, className );

   if (GetTransformationHandler().GetTransformationList().GetSize()!=0 &&
       GetMethodType() != Types::kLikelihood &&
       GetMethodType() != Types::kHMatrix) {
      fout << "            Transform( iV, -1 );" << std::endl;
   }
   fout << "            retval = GetMvaValue__( iV );" << std::endl;
   fout << "         }" << std::endl;
   fout << "         else {" << std::endl;
   if (GetTransformationHandler().GetTransformationList().GetSize()!=0 &&
       GetMethodType() != Types::kLikelihood &&
       GetMethodType() != Types::kHMatrix) {
      fout << "            std::vector<double> iV;" << std::endl;
      fout << "            int ivar = 0;" << std::endl;
      fout << "            for (std::vector<double>::const_iterator varIt = inputValues.begin();" << std::endl;
      fout << "                 varIt != inputValues.end(); varIt++, ivar++) {" << std::endl;
      fout << "               iV.push_back(*varIt);" << std::endl;
      fout << "            }" << std::endl;
      fout << "            Transform( iV, -1 );" << std::endl;
      fout << "            retval = GetMvaValue__( iV );" << std::endl;
   }
   else {
      fout << "            retval = GetMvaValue__( inputValues );" << std::endl;
   }
   fout << "         }" << std::endl;
   fout << "      }" << std::endl;
   fout << std::endl;
   fout << "      return retval;" << std::endl;
   fout << "   }" << std::endl;

   
   // create output for transformation - if any
   if (GetTransformationHandler().GetTransformationList().GetSize()!=0)
      GetTransformationHandler().MakeFunction(fout, className,2);

   // close the file
   fout.close();
}
