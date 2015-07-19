/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAMulticlass                                                     *
 *                                                                                *
 * This macro provides a simple example for the training and testing of the TMVA  *
 * multiclass classification                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"


#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/TMVAMultiClassGui.h"


using namespace TMVA;

void Multiclass()
{
   
   // This loads the library
   TMVA::Tools::Instance();

   // to get access to the GUI and all tmva macros
   // TString tmva_dir(TString(gRootDir) + "/tmva");
   // if(gSystem->Getenv("TMVASYS"))
   //    tmva_dir = TString(gSystem->Getenv("TMVASYS"));
   // gROOT->SetMacroPath(tmva_dir + "/test/:" + gROOT->GetMacroPath() );
   // gROOT->ProcessLine(".L TMVAMultiClassGui.C");

   
   //---------------------------------------------------------------
   // default MVA methods to be trained + tested
   std::map<std::string,int> Use;
   Use["MLP"]             = 1;
   Use["BDTG"]            = 1;
   Use["FDA_GA"]          = 0;
   Use["PDEFoam"]         = 0;
   //---------------------------------------------------------------
   
   std::cout << std::endl;
   std::cout << "==> Start TMVAMulticlass" << std::endl;
   
   // Create a new root output file.
   TString outfileName = "TMVAMulticlass.root";
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   
   TMVA::Factory *factory = new TMVA::Factory( "TMVAMulticlass", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );
   TMVA::DataLoader  *loader1=new TMVA::DataLoader("dataset1");
   TMVA::DataLoader  *loader2=new TMVA::DataLoader("dataset2");
   
   loader1->AddVariable( "var1", 'F' );
   loader1->AddVariable( "var2", "Variable 2", "", 'F' );
   loader1->AddVariable( "var3", "Variable 3", "units", 'F' );
   loader1->AddVariable( "var4", "Variable 4", "units", 'F' );

   loader2->AddVariable( "var1", 'F' );
   loader2->AddVariable( "var2", "Variable 2", "", 'F' );
   loader2->AddVariable( "var3", "Variable 3", "units", 'F' );
   loader2->AddVariable( "var4", "Variable 4", "units", 'F' );
   TFile *input(0);
   TString fname = "./tmva_example_multiple_background.root";
   if (!gSystem->AccessPathName( fname )) {
      // first we try to find the file in the local directory
      std::cout << "--- TMVAMulticlass   : Accessing " << fname << std::endl;
      input = TFile::Open( fname );
   }
   else {
      std::cout << "Creating testdata...." << std::endl;
      gROOT->ProcessLine(".L createData.C");
      gROOT->ProcessLine("create_MultipleBackground(2000)");
      std::cout << " created tmva_example_multiple_background.root for tests of the multiclass features"<<std::endl;
      input = TFile::Open( fname );
   }
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }

   TTree *signal      = (TTree*)input->Get("TreeS");
   TTree *background0 = (TTree*)input->Get("TreeB0");
   TTree *background1 = (TTree*)input->Get("TreeB1");
   TTree *background2 = (TTree*)input->Get("TreeB2");
   
//    gROOT->cd( outfileName+TString(":/") );
   loader1->AddTree    (signal,"Signal");
   loader1->AddTree    (background0,"bg0");
   loader1->AddTree    (background1,"bg1");
   loader1->AddTree    (background2,"bg2");
   
   loader2->AddTree    (signal,"Signal");
   loader2->AddTree    (background0,"bg0");
   loader2->AddTree    (background1,"bg1");
   loader2->AddTree    (background2,"bg2");

   loader1->PrepareTrainingAndTestTree( "","SplitMode=Random:NormMode=NumEvents:!V" );
   
   loader2->PrepareTrainingAndTestTree( "","SplitMode=Random:NormMode=NumEvents:!V" );

   factory->BookMethod( loader1,TMVA::Types::kBDT, "BDTG", "!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=20:MaxDepth=2");
   factory->BookMethod( loader1,TMVA::Types::kMLP, "MLP", "!V:NeuronType=tanh:NCycles=1000:HiddenLayers=N+5,5:TestRate=5:EstimatorType=MSE");
   factory->BookMethod( loader1,TMVA::Types::kFDA, "FDA_GA", "!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );
   factory->BookMethod( loader1,TMVA::Types::kPDEFoam, "PDEFoam", "!V:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );
   
   factory->BookMethod( loader2,TMVA::Types::kBDT, "BDTG", "!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=20:MaxDepth=2");
   factory->BookMethod( loader2,TMVA::Types::kMLP, "MLP", "!V:NeuronType=tanh:NCycles=1000:HiddenLayers=N+5,5:TestRate=5:EstimatorType=MSE");
   factory->BookMethod( loader2,TMVA::Types::kFDA, "FDA_GA", "!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );
   factory->BookMethod( loader2,TMVA::Types::kPDEFoam, "PDEFoam", "!V:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );
  // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------
   
   // Save the output
   outputFile->Close();
   
   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;
   
   delete factory;
   delete loader1;
   delete loader2;
   
}
