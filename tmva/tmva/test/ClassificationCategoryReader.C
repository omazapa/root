#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

// two types of category methods are implemented
Bool_t UseOffsetMethod = kTRUE;

void ClassificationCategoryReader()
{
   // ---------------------------------------------------------------

   std::cout << std::endl
             << "==> Start TMVAClassificationCategoryApplication" << std::endl;

   // --- Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and spectators and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   Float_t var1, var2, var3, var4, eta;
   reader->AddVariable( "var1", &var1 );
   reader->AddVariable( "var2", &var2 );
   reader->AddVariable( "var3", &var3 );
   reader->AddVariable( "var4", &var4 );

   reader->AddSpectator( "eta", &eta );

   // --- Book the MVA methods

   TString weightfile = "dataset1/weights/TMVAClassificationCategory_FisherCat.weights.xml";
   reader->BookMVA( "FisherCat method", weightfile ); 
   
   UInt_t nbin = 100;
   TH1*    hist= new TH1F( "MVA_FisherCat",       "MVA_FisherCat",     nbin, -4, 4 );

   TString fname = "data/";
   TFile *input(0);
   if (UseOffsetMethod) fname += "toy_sigbkg_categ_offset.root";
   else                 fname += "toy_sigbkg_categ_varoff.root";
   
   if (!gSystem->AccessPathName( fname )) {
      gSystem->mkdir("data");
      // file does not exist in local directory
      gSystem->Exec("cd data;curl -O http://files.oproject.org/root/tmva/data/toy_sigbkg_categ_offset.root");
      gSystem->Exec("cd data;curl -O http://files.oproject.org/root/tmva/data/toy_sigbkg_categ_varoff.root");
   
      // first we try to find tmva_example.root in the local directory
   }
   std::cout << "--- TMVAClassificationCategory: Accessing " << fname << std::endl;
   input = TFile::Open( fname );

   
   if (!input) {
      std::cout << "ERROR: could not open data file: " << fname << std::endl;
      exit(1);
   }

   // --- Event loop

   // Prepare the tree
   // - here the variable names have to corresponds to your tree
   // - you can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   TTree* theTree = (TTree*)input->Get("TreeS");
   std::cout << "--- Use signal sample for evalution" << std::endl;
   theTree->SetBranchAddress( "var1", &var1 );
   theTree->SetBranchAddress( "var2", &var2 );
   theTree->SetBranchAddress( "var3", &var3 );
   theTree->SetBranchAddress( "var4", &var4 );

   theTree->SetBranchAddress( "eta",  &eta ); // spectator

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      theTree->GetEntry(ievt);

      // --- Return the MVA outputs and fill into histograms
         hist->Fill( reader->EvaluateMVA( "FisherCat method" ) );         

   }
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   // --- Write histograms

   TFile *target  = new TFile( "TMVApp.root","RECREATE" );
   hist->Write();

   target->Close();
   std::cout << "--- Created root file: \"TMVApp.root\" containing the MVA output histograms" << std::endl;

   delete reader;
   std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
}
