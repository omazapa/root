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
#include "TMVA/MethodC50.h"
#include "TMVA/MethodRSVM.h"


#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include<TRInterface.h>


//Global Objects to use from ROOT prompt for debug and testing
ROOT::R::TRInterface &r=ROOT::R::TRInterface::Instance();

void higgs()
{
   // This loads the library
   TMVA::Tools::Instance();
   
   TString outfileName( "TMVA-Higgs.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   TMVA::Factory *factory = new TMVA::Factory( "RMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
   /////////////////////////
   //Adding Training Data //
   /////////////////////////
   ROOT::R::TRDataFrame TrainData;
   r<<"TrainData=read.csv('training.csv')";
   
   //TrainData is adataframe without the first colunm Event ID is not needed//
   r<<"TrainData=TrainData[2:33]";
   r["TrainData"]>>TrainData;
   UInt_t ncols=TrainData.GetNcols();
   UInt_t nrows=TrainData.GetNrows();
   
   //Matrix with data without col 1= EventId , col 32=Weight and col 33=Label
   TMatrixF TrainDataMat;
   TrainDataMat.ResizeTo(nrows,30);//30 cols less 1,32 and 33
   r["as.matrix(TrainData[1:30])"]>>TrainDataMat;
   Rcpp::StringVector fColNames(30);
   r["names(TrainData[1:30])"]>>fColNames;//30 cols less 1,32 and 33
   
   Rcpp::StringVector TrainLabels(nrows);//ifis signal/background
   TrainData["Label"]>>TrainLabels;
   std::vector<Float_t> TrainWeights(nrows);
   TrainData["Weight"]>>TrainWeights;
//
//   //ncols 30 because the last two cols are weights and class and the first is EventID
   for(int i=0;i<30;i++)    factory->AddVariable( std::string(fColNames[i]).c_str(), 'F' );

   
   ////////////////////////////////////////////
   //Adding training Data every event at time//
   ////////////////////////////////////////////
   for(UInt_t i=0;i<nrows;i++)
   {
       
       std::vector<Double_t> ev(30);//ncols -3 because last two are Eventid, weights and class(s/b))
       for(UInt_t j=0;j<30;j++) ev[j]=TrainDataMat[i][j];
       if(TrainLabels[i]=="s") factory->AddSignalTrainingEvent( ev, TrainWeights[i] ); //adding signals 
       else factory->AddBackgroundTrainingEvent( ev,TrainWeights[i]); //adding background
   }


   ////////////////////////
   //Adding Testing Data //
   ///////////////////////
   ROOT::R::TRDataFrame TestData;
   r<<"TestData=read.csv('test.csv')";
   
   //TrainData is adataframe without the first colunm Event ID is not needed//
   r<<"TestData=TestData[2:31]";
   r["TestData"]>>TestData;
   ncols=TestData.GetNcols();
   nrows=TestData.GetNrows();
   
   //Matrix with data without col 1= EventId , col 32=Weight and col 33=Label
   TMatrixF TestDataMat;
   TestDataMat.ResizeTo(nrows,30);//30 cols less 1,32 and 33
   r["as.matrix(TestData[1:30])"]>>TestDataMat;
   
   
   //weights no given  and labels given in random_submission.csv
   r<<"RandomSubmission<-read.csv('random_submission.csv')";
   ROOT::R::TRDataFrame Class;
   r["RandomSubmission['Class']"]>>Class;
   
   Rcpp::StringVector TestLabels(nrows);//ifis signal/background
   Class["Class"]>>TestLabels;
   
   ////////////////////////////////////////////
   //Adding training Data every event at time//
   ////////////////////////////////////////////
   for(UInt_t i=0;i<nrows;i++)
   {
       
       std::vector<Double_t> ev(30);//ncols -3 because last two are Eventid, weights and class(s/b))
       for(UInt_t j=0;j<30;j++) ev[j]=TestDataMat[i][j];
       
       if(TestLabels[i]=="s") factory->AddSignalTestEvent( ev ); //adding signals 
       else factory->AddBackgroundTestEvent( ev); //adding background
   }


   
   factory->SetBackgroundWeightExpression( "weight" );
   
   
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
   factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   // Boosted Decision Trees

   factory->BookMethod( TMVA::Types::kC50, "C50", "!H:NTrials=50:Rules=kFALSE:ControlSubSet=kFALSE:ControlBands=0:ControlWinnow=kFALSE:ControlNoGlobalPruning=kTRUE:ControlCF=0.25:ControlMinCases=2:ControlFuzzyThreshold=kTRUE:ControlSample=0:ControlEarlyStopping=kTRUE:!V" );
   
   factory->BookMethod( TMVA::Types::kRSNNS, "RMLP","!H:VarTransform=N:Size=c(5):Maxit=800:InitFunc=Randomize_Weights:LearnFunc=Std_Backpropagation:LearnFuncParams=c(0.2,0):!V" );
    
   factory->BookMethod( TMVA::Types::kRSVM, "RSVM","!H:Kernel=linear:Type=C-classification:VarTransform=Norm:Probability=kTRUE:Tolerance=0.001:!V" );

   
      // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

//    ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();
   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << "TMVA-Higgs.root" << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;
    r.SetVerbose(1);

}
