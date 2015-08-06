#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"


#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/DataLoader.h"
#include "TMVA/MethodC50.h"
#include "TRInterface.h"

void Classification()
{
   // This loads the library
   TMVA::Tools::Instance();
   ROOT::R::TRInterface &r=ROOT::R::TRInterface::Instance();

    // --------------------------------------------------------------------------------------------------

   // --- Here the preparation phase begins

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "TMVA.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory is 
   // the only TMVA object you have to interact with
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:Silent:ROC=1:Color:DrawProgressBar:Correlations=kFALSE:AnalysisType=Classification" );
   
   TMVA::DataLoader *loader1=new TMVA::DataLoader("dataset1");
   TMVA::DataLoader *loader2=new TMVA::DataLoader("dataset2");
    // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
   loader1->AddVariable( "myvar1 := var1+var2", 'F' );
   loader1->AddVariable( "myvar2 := var1-var2", "Expression 2", "", 'F' );
   loader1->AddVariable( "var3",                "Variable 3", "units", 'F' );
   loader1->AddVariable( "var4",                "Variable 4", "units", 'F' );

   loader2->AddVariable( "myvar1 := var1+var2", 'F' );
   loader2->AddVariable( "myvar2 := var1-var2", "Expression 2", "", 'F' );
   loader2->AddVariable( "var3",                "Variable 3", "units", 'F' );
   loader2->AddVariable( "var4",                "Variable 4", "units", 'F' );

   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
   loader1->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
   loader1->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );

   loader2->AddSpectator( "spec1 := var1*3",  "Spectator 1", "units", 'F' );
   loader2->AddSpectator( "spec2 := var1*4",  "Spectator 2", "units", 'F' );
     // Read training and test data
   // (it is also possible to use ASCII format as input -> see TMVA Users Guide)
   TString fname = "./tmva_class_example.root";
   
   if (gSystem->AccessPathName( fname ))  // file does not exist in local directory
      gSystem->Exec("curl -O http://root.cern.ch/files/tmva_class_example.root");
   
   TFile *input = TFile::Open( fname );
   
   std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;
   
   // --- Register the training and test trees

   TTree *tsignal     = (TTree*)input->Get("TreeS");
   TTree *tbackground = (TTree*)input->Get("TreeB");
   
   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;
   
   // You can add an arbitrary number of signal or background trees
   loader1->AddSignalTree    ( tsignal,     signalWeight     );
   loader1->AddBackgroundTree( tbackground, backgroundWeight );
 
   // You can add an arbitrary number of signal or background trees
   loader2->AddSignalTree    ( tsignal,     signalWeight     );
   loader2->AddBackgroundTree( tbackground, backgroundWeight );
   
    // Set individual event weights (the variables must exist in the original TTree)
   //    for signal    : factory->SetSignalWeightExpression    ("weight1*weight2");
   //    for background: factory->SetBackgroundWeightExpression("weight1*weight2");
   loader1->SetBackgroundWeightExpression( "weight" );
   
   loader2->SetBackgroundWeightExpression( "weight" );
   
      // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

      // Tell the factory how to use the training and testing events
   //
   // If no numbers of events are given, half of the events in the tree are used 
   // for training, and the other half for testing:
   //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
   // To also specify the number of testing events, use:
   //    factory->PrepareTrainingAndTestTree( mycut,
   //                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
   loader1->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

   loader2->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=1000:nTrain_Background=1000:nTest_Signal=1000:nTest_Background=1000:SplitMode=Random:NormMode=NumEvents:!V" );
   // Boosted Decision Trees
   factory->BookMethod( loader1, TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

//    factory->BookMethod( loader1, TMVA::Types::kBDT, "BDT",
//                            "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

//    factory->BookMethod( loader1, TMVA::Types::kBDT, "BDTB",
//                            "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

   // Boosted Decision Trees
   factory->BookMethod( loader2, TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );


   factory->BookMethod( loader2, TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

   
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

//     delete factory;
    delete loader1;
    delete loader2;
}
