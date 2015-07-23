//Based on example to copy file
//Autor: Omar Zapata, Base in Rene Brun Code
#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TString.h"


//dataset name and TMVA output file
void TMVAMultiClassGui(TString dataset,const char *file="TMVA.root");



void CopyDir(const char *dataset,TDirectory *source) {
   //copy all objects and subdirs of directory source as a subdir of the current directory
      TDirectory *savdir = gDirectory;
      TDirectory *adir;
      if(dataset==TString(""))
      {
	  adir = savdir->mkdir(source->GetName());
      }else
      {     
	  adir = savdir;
      }
      adir->cd();
      //loop on all entries of this directory
      TKey *key;
      TIter nextkey(source->GetListOfKeys());          
      while ((key = (TKey*)nextkey())) {
	  const char *classname = key->GetClassName();
	  TClass *cl = gROOT->GetClass(classname);
	  if (!cl) continue;
	  if (cl->InheritsFrom(TDirectory::Class())) {
	    source->cd(key->GetName());
	    TDirectory *subdir = gDirectory;
	    adir->cd();
	    CopyDir("",subdir);
	    adir->cd();
	  } else if (cl->InheritsFrom(TTree::Class())) {
	    TTree *T = (TTree*)source->Get(key->GetName());
	    adir->cd();
	    TTree *newT = T->CloneTree(-1,"fast");
	    newT->Write();
	  } else {
	    source->cd();
	    TObject *obj = key->ReadObj();
	    adir->cd();
	    obj->Write();
	    delete obj;
	}
      }
      adir->SaveSelf(kTRUE);
      savdir->cd();
}
void CopyFile(const char* dataset,const char *fname) {
   //Copy all objects and subdirs of file fname as a subdir of the current directory
   TDirectory *target = gDirectory;

   TFile *f = TFile::Open(fname);
   if (!f || f->IsZombie()) {
      printf("Cannot copy file: %s\n",fname);
      target->cd();
      return;
   }
   target->cd();
   CopyDir(dataset,f->GetDirectory(dataset));
   delete f;
   target->cd();
}

void TMVAMultiClassGui(TString dataset,const char *file="TMVA.root")
{
  TString datasetfile=dataset+".root";
  TFile *f = new TFile(datasetfile.Data(),"recreate");
  
  CopyFile(dataset.Data(),file);
  f->Close();
  TMVA::TMVAMultiClassGui(datasetfile.Data());
}