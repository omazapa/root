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

void ClassificationAll()
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
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
   
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
                                        "nTrain_Signal=1000:nTrain_Background=1000:nTest_Signal=1000:nTest_Background=1000:SplitMode=Random:NormMode=NumEvents:!V" );

   loader2->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=500:nTrain_Background=500:nTest_Signal=500:nTest_Background=500:SplitMode=Random:NormMode=NumEvents:!V" );
      factory->BookMethod( loader1,TMVA::Types::kCuts, "Cuts",
                           "!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

      factory->BookMethod( loader1,TMVA::Types::kCuts, "CutsD",
                           "!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

      factory->BookMethod( loader1,TMVA::Types::kCuts, "CutsPCA",
                           "!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

      factory->BookMethod( loader1,TMVA::Types::kCuts, "CutsGA",
                           "!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

      factory->BookMethod( loader1,TMVA::Types::kCuts, "CutsSA",
                           "!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   // Likelihood ("naive Bayes estimator")
      factory->BookMethod( loader1,TMVA::Types::kLikelihood, "Likelihood",
                           "!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

   // Decorrelated likelihood
      factory->BookMethod( loader1,TMVA::Types::kLikelihood, "LikelihoodD",
                           "!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );

   // PCA-transformed likelihood
      factory->BookMethod( loader1,TMVA::Types::kLikelihood, "LikelihoodPCA",
                           "!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 

   // Use a kernel density estimator to approximate the PDFs
      factory->BookMethod( loader1,TMVA::Types::kLikelihood, "LikelihoodKDE",
                           "!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 

   // Use a variable-dependent mix of splines and kernel density estimator
      factory->BookMethod( loader1,TMVA::Types::kLikelihood, "LikelihoodMIX",
                           "!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

   // Test the multi-dimensional probability density estimator
   // here are the options strings for the MinMax and RMS methods, respectively:
   //      "!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
   //      "!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
      factory->BookMethod( loader1,TMVA::Types::kPDERS, "PDERS",
                           "!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

      factory->BookMethod( loader1,TMVA::Types::kPDERS, "PDERSD",
                           "!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

      factory->BookMethod( loader1,TMVA::Types::kPDERS, "PDERSPCA",
                           "!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

//    Multi-dimensional likelihood estimator using self-adapting phase-space binning
      factory->BookMethod( loader1,TMVA::Types::kPDEFoam, "PDEFoam",
                           "!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

      factory->BookMethod( loader1,TMVA::Types::kPDEFoam, "PDEFoamBoost",
                           "!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

   // K-Nearest Neighbour classifier (KNN)
      factory->BookMethod( loader1,TMVA::Types::kKNN, "KNN",
                           "nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

   // H-Matrix (chi2-squared) method
      factory->BookMethod( loader1,TMVA::Types::kHMatrix, "HMatrix", "!V:VarTransform=None" );

   // Linear discriminant (same as Fisher discriminant)
      factory->BookMethod( loader1,TMVA::Types::kLD, "LD", "!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher discriminant (same as LD)
      factory->BookMethod( loader1,TMVA::Types::kFisher, "Fisher", "!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher with Gauss-transformed input variables
      factory->BookMethod( loader1,TMVA::Types::kFisher, "FisherG", "!V:VarTransform=Gauss" );

   // Composite classifier: ensemble (tree) of boosted Fisher classifiers
      factory->BookMethod( loader1,TMVA::Types::kFisher, "BoostedFisher", 
                           "!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );

   // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
      factory->BookMethod( loader1,TMVA::Types::kFDA, "FDA_MC",
                           "!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

   // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( loader1,TMVA::Types::kFDA, "FDA_GA",
                           "!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

      // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( loader1,TMVA::Types::kFDA, "FDA_SA",
                           "!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

      factory->BookMethod( loader1,TMVA::Types::kFDA, "FDA_MT",
                           "!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

      factory->BookMethod( loader1,TMVA::Types::kFDA, "FDA_GAMT",
                           "!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

      factory->BookMethod( loader1,TMVA::Types::kFDA, "FDA_MCMT",
                           "!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

   // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
      factory->BookMethod( loader1,TMVA::Types::kMLP, "MLP", "!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

      factory->BookMethod( loader1,TMVA::Types::kMLP, "MLPBFGS", "!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

      factory->BookMethod( loader1,TMVA::Types::kMLP, "MLPBNN", "!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

//    CF(Clermont-Ferrand)ANN
      factory->BookMethod( loader1,TMVA::Types::kCFMlpANN, "CFMlpANN", "!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  

   // Tmlp(Root)ANN
      factory->BookMethod( loader1,TMVA::Types::kTMlpANN, "TMlpANN", "!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

   // Support Vector Machine
      factory->BookMethod( loader1,TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

   // Boosted Decision Trees
   // Gradient Boost
      factory->BookMethod( loader1,TMVA::Types::kBDT, "BDTG",
                           "!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

   // Adaptive Boost
      factory->BookMethod( loader1,TMVA::Types::kBDT, "BDT",
                           "!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

      // Bagging
      factory->BookMethod( loader1,TMVA::Types::kBDT, "BDTB",
                           "!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

      // Decorrelation + Adaptive Boost
      factory->BookMethod( loader1,TMVA::Types::kBDT, "BDTD",
                           "!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );

      // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
      factory->BookMethod( loader1,TMVA::Types::kBDT, "BDTMitFisher",
                           "!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

   // RuleFit -- TMVA implementation of Friedman's method
      factory->BookMethod( loader1,TMVA::Types::kRuleFit, "RuleFit",
                           "!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );
// 
// 
//       
//       
//       
//       
//       
//       
//       
//       
//        factory->BookMethod( loader2,TMVA::Types::kCuts, "Cuts",
//                            "!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );
// 
//       factory->BookMethod( loader2,TMVA::Types::kCuts, "CutsD",
//                            "!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );
// 
//       factory->BookMethod( loader2,TMVA::Types::kCuts, "CutsPCA",
//                            "!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );
// 
//       factory->BookMethod( loader2,TMVA::Types::kCuts, "CutsGA",
//                            "!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );
// 
//       factory->BookMethod( loader2,TMVA::Types::kCuts, "CutsSA",
//                            "!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );
// 
//    // Likelihood ("naive Bayes estimator")
//       factory->BookMethod( loader2,TMVA::Types::kLikelihood, "Likelihood",
//                            "!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );
// 
//    // Decorrelated likelihood
//       factory->BookMethod( loader2,TMVA::Types::kLikelihood, "LikelihoodD",
//                            "!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );
// 
//    // PCA-transformed likelihood
//       factory->BookMethod( loader2,TMVA::Types::kLikelihood, "LikelihoodPCA",
//                            "!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 
// 
//    // Use a kernel density estimator to approximate the PDFs
//       factory->BookMethod( loader2,TMVA::Types::kLikelihood, "LikelihoodKDE",
//                            "!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 
// 
//    // Use a variable-dependent mix of splines and kernel density estimator
//       factory->BookMethod( loader2,TMVA::Types::kLikelihood, "LikelihoodMIX",
//                            "!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 
// 
//    // Test the multi-dimensional probability density estimator
//    // here are the options strings for the MinMax and RMS methods, respectively:
//    //      "!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
//    //      "!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
//       factory->BookMethod( loader2,TMVA::Types::kPDERS, "PDERS",
//                            "!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );
// 
//       factory->BookMethod( loader2,TMVA::Types::kPDERS, "PDERSD",
//                            "!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );
// 
//       factory->BookMethod( loader2,TMVA::Types::kPDERS, "PDERSPCA",
//                            "!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );
// 
// //    Multi-dimensional likelihood estimator using self-adapting phase-space binning
//       factory->BookMethod( loader2,TMVA::Types::kPDEFoam, "PDEFoam",
//                            "!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );
// 
//       factory->BookMethod( loader2,TMVA::Types::kPDEFoam, "PDEFoamBoost",
//                            "!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );
// 
//    // K-Nearest Neighbour classifier (KNN)
//       factory->BookMethod( loader2,TMVA::Types::kKNN, "KNN",
//                            "nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );
// 
//    // H-Matrix (chi2-squared) method
//       factory->BookMethod( loader2,TMVA::Types::kHMatrix, "HMatrix", "!V:VarTransform=None" );
// 
//    // Linear discriminant (same as Fisher discriminant)
//       factory->BookMethod( loader2,TMVA::Types::kLD, "LD", "!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );
// 
//    // Fisher discriminant (same as LD)
//       factory->BookMethod( loader2,TMVA::Types::kFisher, "Fisher", "!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );
// 
//    // Fisher with Gauss-transformed input variables
//       factory->BookMethod( loader2,TMVA::Types::kFisher, "FisherG", "!V:VarTransform=Gauss" );
// 
//    // Composite classifier: ensemble (tree) of boosted Fisher classifiers
//       factory->BookMethod( loader2,TMVA::Types::kFisher, "BoostedFisher", 
//                            "!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );
// 
//    // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
//       factory->BookMethod( loader2,TMVA::Types::kFDA, "FDA_MC",
//                            "!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );
// 
//    // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
//       factory->BookMethod( loader2,TMVA::Types::kFDA, "FDA_GA",
//                            "!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );
// 
//       // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
//       factory->BookMethod( loader2,TMVA::Types::kFDA, "FDA_SA",
//                            "!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );
// 
//       factory->BookMethod( loader2,TMVA::Types::kFDA, "FDA_MT",
//                            "!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );
// 
//       factory->BookMethod( loader2,TMVA::Types::kFDA, "FDA_GAMT",
//                            "!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );
// 
//       factory->BookMethod( loader2,TMVA::Types::kFDA, "FDA_MCMT",
//                            "!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

//    // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
//       factory->BookMethod( loader2,TMVA::Types::kMLP, "MLP", "!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

//       factory->BookMethod( loader2,TMVA::Types::kMLP, "MLPBFGS", "!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

//       factory->BookMethod( loader2,TMVA::Types::kMLP, "MLPBNN", "!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

//    // CF(Clermont-Ferrand)ANN
//       factory->BookMethod( loader2,TMVA::Types::kCFMlpANN, "CFMlpANN", "!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  

   // Tmlp(Root)ANN
//       factory->BookMethod( loader2,TMVA::Types::kTMlpANN, "TMlpANN", "!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

//    // Support Vector Machine
//       factory->BookMethod( loader2,TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );
// 
//    // Boosted Decision Trees
//    // Gradient Boost
//       factory->BookMethod( loader2,TMVA::Types::kBDT, "BDTG",
//                            "!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
// 
//    // Adaptive Boost
//       factory->BookMethod( loader2,TMVA::Types::kBDT, "BDT",
//                            "!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
// 
//       // Bagging
//       factory->BookMethod( loader2,TMVA::Types::kBDT, "BDTB",
//                            "!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );
// 
//       // Decorrelation + Adaptive Boost
//       factory->BookMethod( loader2,TMVA::Types::kBDT, "BDTD",
//                            "!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );
// 
//       // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
//       factory->BookMethod( loader2,TMVA::Types::kBDT, "BDTMitFisher",
//                            "!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );
// 
//    // RuleFit -- TMVA implementation of Friedman's method
      factory->BookMethod( loader2,TMVA::Types::kRuleFit, "RuleFit",
                           "!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );
//      
      
      
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
