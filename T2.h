//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec 24 17:02:17 2020 by ROOT version 5.34/36
// from TTree T2/INODIGI
// found on file: corsika76300_FLUKA_SIBYLL_Pethuraj_20201203_5M_1_digi.root
//////////////////////////////////////////////////////////

#ifndef T2_h
#define T2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class T2 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          irun;
   UInt_t          ievt;
   UInt_t          ngent;
   Int_t           pidin[10];   //[ngent]
   Float_t         ievt_wt;
   Int_t           intxn_id;
   Float_t         momin[10];   //[ngent]
   Float_t         thein[10];   //[ngent]
   Float_t         phiin[10];   //[ngent]
   Float_t         posxin[10];   //[ngent]
   Float_t         posyin[10];   //[ngent]
   Float_t         poszin[10];   //[ngent]
   UInt_t          ndigiht;
   Int_t           trigx;
   Int_t           trigy;
   UInt_t          ngenerated;
   UInt_t          naperture;
   UInt_t          triggeracceptance;
   UInt_t          stripid[200];   //[ndigiht]
   Int_t           digipdgid[200];   //[ndigiht]
   Int_t           digitime[200];   //[ndigiht]
   Int_t           digitruetime[200];   //[ndigiht]
   Float_t         digienr[200];   //[ndigiht]
   Float_t         digivx[200];   //[ndigiht]
   Float_t         digivy[200];   //[ndigiht]
   Float_t         digivz[200];   //[ndigiht]
   Float_t         digipx[200];   //[ndigiht]
   Float_t         digipy[200];   //[ndigiht]
   Float_t         digipz[200];   //[ndigiht]

   // List of branches
   TBranch        *b_irun;   //!
   TBranch        *b_ievt;   //!
   TBranch        *b_ngent;   //!
   TBranch        *b_pidin;   //!
   TBranch        *b_ievt_wt;   //!
   TBranch        *b_intxn_id;   //!
   TBranch        *b_momin;   //!
   TBranch        *b_thein;   //!
   TBranch        *b_phiin;   //!
   TBranch        *b_posxin;   //!
   TBranch        *b_posyin;   //!
   TBranch        *b_poszin;   //!
   TBranch        *b_ndigiht;   //!
   TBranch        *b_trigx;   //!
   TBranch        *b_trigy;   //!
   TBranch        *b_ngenerated;   //!
   TBranch        *b_naperture;   //!
   TBranch        *b_triggeracceptance;   //!
   TBranch        *b_stripid;   //!
   TBranch        *b_digipdgid;   //!
   TBranch        *b_digitime;   //!
   TBranch        *b_digitruetime;   //!
   TBranch        *b_digienr;   //!
   TBranch        *b_digivx;   //!
   TBranch        *b_digivy;   //!
   TBranch        *b_digivz;   //!
   TBranch        *b_digipx;   //!
   TBranch        *b_digipy;   //!
   TBranch        *b_digipz;   //!

   T2(TTree *tree=0);
   virtual ~T2();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef T2_cxx
T2::T2(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("corsika76300_FLUKA_SIBYLL_Pethuraj_20201203_5M_1_digi.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("corsika76300_FLUKA_SIBYLL_Pethuraj_20201203_5M_1_digi.root");
      }
      f->GetObject("T2",tree);

   }
   Init(tree);
}

T2::~T2()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t T2::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t T2::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void T2::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("irun", &irun, &b_irun);
   fChain->SetBranchAddress("ievt", &ievt, &b_ievt);
   fChain->SetBranchAddress("ngent", &ngent, &b_ngent);
   fChain->SetBranchAddress("pidin", pidin, &b_pidin);
   fChain->SetBranchAddress("ievt_wt", &ievt_wt, &b_ievt_wt);
   fChain->SetBranchAddress("intxn_id", &intxn_id, &b_intxn_id);
   fChain->SetBranchAddress("momin", momin, &b_momin);
   fChain->SetBranchAddress("thein", thein, &b_thein);
   fChain->SetBranchAddress("phiin", phiin, &b_phiin);
   fChain->SetBranchAddress("posxin", posxin, &b_posxin);
   fChain->SetBranchAddress("posyin", posyin, &b_posyin);
   fChain->SetBranchAddress("poszin", poszin, &b_poszin);
   fChain->SetBranchAddress("ndigiht", &ndigiht, &b_ndigiht);
   fChain->SetBranchAddress("trigx", &trigx, &b_trigx);
   fChain->SetBranchAddress("trigy", &trigy, &b_trigy);
   fChain->SetBranchAddress("ngenerated", &ngenerated, &b_ngenerated);
   fChain->SetBranchAddress("naperture", &naperture, &b_naperture);
   /* fChain->SetBranchAddress("triggeracceptance", &triggeracceptance, &b_triggeracceptance); */
   fChain->SetBranchAddress("stripid", stripid, &b_stripid);
   fChain->SetBranchAddress("digipdgid", digipdgid, &b_digipdgid);
   fChain->SetBranchAddress("digitime", digitime, &b_digitime);
   fChain->SetBranchAddress("digitruetime", digitruetime, &b_digitruetime);
   fChain->SetBranchAddress("digienr", digienr, &b_digienr);
   fChain->SetBranchAddress("digivx", digivx, &b_digivx);
   fChain->SetBranchAddress("digivy", digivy, &b_digivy);
   fChain->SetBranchAddress("digivz", digivz, &b_digivz);
   fChain->SetBranchAddress("digipx", digipx, &b_digipx);
   fChain->SetBranchAddress("digipy", digipy, &b_digipy);
   fChain->SetBranchAddress("digipz", digipz, &b_digipz);
   Notify();
}

Bool_t T2::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void T2::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t T2::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef T2_cxx
