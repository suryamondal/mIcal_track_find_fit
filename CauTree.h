//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 22 16:50:05 2017 by ROOT version 5.34/26
// from TTree cautree/cau offset
// found on file: RPC_evtraw-r212230.dat
//////////////////////////////////////////////////////////

#ifndef CauTree_h
#define CauTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TTimeStamp.h>
#include <vector>

#define NO_TDC_CHNL 16      //added by Esha on 22/06/2017



// Fixed size dimensions of array or collections stored in the TTree if any.

class CauTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           CauNum;
   TTimeStamp      *Cautime;
   Int_t           select_line;
   ULong64_t       raw_trig_cnt;
   ULong64_t       acpt_trig_cnt;
   Int_t           cau_status;
   ULong64_t       cau_ref1_l;
   ULong64_t       cau_ref2_l;
   ULong64_t       cau_ref1_t;
   ULong64_t       cau_ref2_t;
   std::vector<unsigned int> *cau_tdc_l_0;
   std::vector<unsigned int> *cau_tdc_ref1_l_0;
   std::vector<unsigned int> *cau_tdc_t_0;
   std::vector<unsigned int> *cau_tdc_ref1_t_0;
   std::vector<unsigned int> *cau_tdc_l_1;
   std::vector<unsigned int> *cau_tdc_ref1_l_1;
   std::vector<unsigned int> *cau_tdc_t_1;
   std::vector<unsigned int> *cau_tdc_ref1_t_1;
   std::vector<unsigned int> *cau_tdc_l_2;
   std::vector<unsigned int> *cau_tdc_ref1_l_2;
   std::vector<unsigned int> *cau_tdc_t_2;
   std::vector<unsigned int> *cau_tdc_ref1_t_2;
   std::vector<unsigned int> *cau_tdc_l_3;
   std::vector<unsigned int> *cau_tdc_ref1_l_3;
   std::vector<unsigned int> *cau_tdc_t_3;
   std::vector<unsigned int> *cau_tdc_ref1_t_3;
   std::vector<unsigned int> *cau_tdc_l_4;
   std::vector<unsigned int> *cau_tdc_ref1_l_4;
   std::vector<unsigned int> *cau_tdc_t_4;
   std::vector<unsigned int> *cau_tdc_ref1_t_4;
   std::vector<unsigned int> *cau_tdc_l_5;
   std::vector<unsigned int> *cau_tdc_ref1_l_5;
   std::vector<unsigned int> *cau_tdc_t_5;
   std::vector<unsigned int> *cau_tdc_ref1_t_5;
   std::vector<unsigned int> *cau_tdc_l_6;
   std::vector<unsigned int> *cau_tdc_ref1_l_6;
   std::vector<unsigned int> *cau_tdc_t_6;
   std::vector<unsigned int> *cau_tdc_ref1_t_6;
   std::vector<unsigned int> *cau_tdc_l_7;
   std::vector<unsigned int> *cau_tdc_ref1_l_7;
   std::vector<unsigned int> *cau_tdc_t_7;
   std::vector<unsigned int> *cau_tdc_ref1_t_7;
   std::vector<unsigned int> *cau_tdc_l_8;
   std::vector<unsigned int> *cau_tdc_ref1_l_8;
   std::vector<unsigned int> *cau_tdc_t_8;
   std::vector<unsigned int> *cau_tdc_ref1_t_8;
   std::vector<unsigned int> *cau_tdc_l_9;
   std::vector<unsigned int> *cau_tdc_ref1_l_9;
   std::vector<unsigned int> *cau_tdc_t_9;
   std::vector<unsigned int> *cau_tdc_ref1_t_9;
   std::vector<unsigned int> *cau_tdc_l_10;
   std::vector<unsigned int> *cau_tdc_ref1_l_10;
   std::vector<unsigned int> *cau_tdc_t_10;
   std::vector<unsigned int> *cau_tdc_ref1_t_10;
   std::vector<unsigned int> *cau_tdc_l_11;
   std::vector<unsigned int> *cau_tdc_ref1_l_11;
   std::vector<unsigned int> *cau_tdc_t_11;
   std::vector<unsigned int> *cau_tdc_ref1_t_11;
   std::vector<unsigned int> *cau_tdc_l_12;
   std::vector<unsigned int> *cau_tdc_ref1_l_12;
   std::vector<unsigned int> *cau_tdc_t_12;
   std::vector<unsigned int> *cau_tdc_ref1_t_12;
   std::vector<unsigned int> *cau_tdc_l_13;
   std::vector<unsigned int> *cau_tdc_ref1_l_13;
   std::vector<unsigned int> *cau_tdc_t_13;
   std::vector<unsigned int> *cau_tdc_ref1_t_13;
   std::vector<unsigned int> *cau_tdc_l_14;
   std::vector<unsigned int> *cau_tdc_ref1_l_14;
   std::vector<unsigned int> *cau_tdc_t_14;
   std::vector<unsigned int> *cau_tdc_ref1_t_14;
   std::vector<unsigned int> *cau_tdc_l_15;
   std::vector<unsigned int> *cau_tdc_ref1_l_15;
   std::vector<unsigned int> *cau_tdc_t_15;
   std::vector<unsigned int> *cau_tdc_ref1_t_15;

   // List of branches
   TBranch        *b_CauNum;   //!
   TBranch        *b_Cautime;   //!
   TBranch        *b_select_line;   //!
   TBranch        *b_raw_trig_cnt;   //!
   TBranch        *b_acpt_trig_cnt;   //!
   TBranch        *b_cau_status;   //!
   TBranch        *b_cau_ref1_l;   //!
   TBranch        *b_cau_ref2_l;   //!
   TBranch        *b_cau_ref1_t;   //!
   TBranch        *b_cau_ref2_t;   //!
   TBranch        *b_cau_tdc_l_0;   //!
   TBranch        *b_cau_tdc_ref1_l_0;   //!
   TBranch        *b_cau_tdc_t_0;   //!
   TBranch        *b_cau_tdc_ref1_t_0;   //!
   TBranch        *b_cau_tdc_l_1;   //!
   TBranch        *b_cau_tdc_ref1_l_1;   //!
   TBranch        *b_cau_tdc_t_1;   //!
   TBranch        *b_cau_tdc_ref1_t_1;   //!
   TBranch        *b_cau_tdc_l_2;   //!
   TBranch        *b_cau_tdc_ref1_l_2;   //!
   TBranch        *b_cau_tdc_t_2;   //!
   TBranch        *b_cau_tdc_ref1_t_2;   //!
   TBranch        *b_cau_tdc_l_3;   //!
   TBranch        *b_cau_tdc_ref1_l_3;   //!
   TBranch        *b_cau_tdc_t_3;   //!
   TBranch        *b_cau_tdc_ref1_t_3;   //!
   TBranch        *b_cau_tdc_l_4;   //!
   TBranch        *b_cau_tdc_ref1_l_4;   //!
   TBranch        *b_cau_tdc_t_4;   //!
   TBranch        *b_cau_tdc_ref1_t_4;   //!
   TBranch        *b_cau_tdc_l_5;   //!
   TBranch        *b_cau_tdc_ref1_l_5;   //!
   TBranch        *b_cau_tdc_t_5;   //!
   TBranch        *b_cau_tdc_ref1_t_5;   //!
   TBranch        *b_cau_tdc_l_6;   //!
   TBranch        *b_cau_tdc_ref1_l_6;   //!
   TBranch        *b_cau_tdc_t_6;   //!
   TBranch        *b_cau_tdc_ref1_t_6;   //!
   TBranch        *b_cau_tdc_l_7;   //!
   TBranch        *b_cau_tdc_ref1_l_7;   //!
   TBranch        *b_cau_tdc_t_7;   //!
   TBranch        *b_cau_tdc_ref1_t_7;   //!
   TBranch        *b_cau_tdc_l_8;   //!
   TBranch        *b_cau_tdc_ref1_l_8;   //!
   TBranch        *b_cau_tdc_t_8;   //!
   TBranch        *b_cau_tdc_ref1_t_8;   //!
   TBranch        *b_cau_tdc_l_9;   //!
   TBranch        *b_cau_tdc_ref1_l_9;   //!
   TBranch        *b_cau_tdc_t_9;   //!
   TBranch        *b_cau_tdc_ref1_t_9;   //!
   TBranch        *b_cau_tdc_l_10;   //!
   TBranch        *b_cau_tdc_ref1_l_10;   //!
   TBranch        *b_cau_tdc_t_10;   //!
   TBranch        *b_cau_tdc_ref1_t_10;   //!
   TBranch        *b_cau_tdc_l_11;   //!
   TBranch        *b_cau_tdc_ref1_l_11;   //!
   TBranch        *b_cau_tdc_t_11;   //!
   TBranch        *b_cau_tdc_ref1_t_11;   //!
   TBranch        *b_cau_tdc_l_12;   //!
   TBranch        *b_cau_tdc_ref1_l_12;   //!
   TBranch        *b_cau_tdc_t_12;   //!
   TBranch        *b_cau_tdc_ref1_t_12;   //!
   TBranch        *b_cau_tdc_l_13;   //!
   TBranch        *b_cau_tdc_ref1_l_13;   //!
   TBranch        *b_cau_tdc_t_13;   //!
   TBranch        *b_cau_tdc_ref1_t_13;   //!
   TBranch        *b_cau_tdc_l_14;   //!
   TBranch        *b_cau_tdc_ref1_l_14;   //!
   TBranch        *b_cau_tdc_t_14;   //!
   TBranch        *b_cau_tdc_ref1_t_14;   //!
   TBranch        *b_cau_tdc_l_15;   //!
   TBranch        *b_cau_tdc_ref1_l_15;   //!
   TBranch        *b_cau_tdc_t_15;   //!
   TBranch        *b_cau_tdc_ref1_t_15;   //!

   CauTree(TTree *tree=0);
   virtual ~CauTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   std::vector<unsigned int> *cau_tdc_l[NO_TDC_CHNL];
   std::vector<unsigned int> *cau_tdc_ref1_l[NO_TDC_CHNL];

   std::vector<unsigned int> *cau_tdc_t[NO_TDC_CHNL];
   std::vector<unsigned int> *cau_tdc_ref1_t[NO_TDC_CHNL];
};

#endif

#ifdef CauTree_cxx
CauTree::CauTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RPC_evtraw-r212230.dat");
      if (!f || !f->IsOpen()) {
         f = new TFile("RPC_evtraw-r212230.dat");
      }
      f->GetObject("cautree",tree);

   }
   Init(tree);
}

CauTree::~CauTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CauTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t CauTree::LoadTree(Long64_t entry)
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

void CauTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Cautime = 0;
   cau_tdc_l_0 = 0;
   cau_tdc_ref1_l_0 = 0;
   cau_tdc_t_0 = 0;
   cau_tdc_ref1_t_0 = 0;
   cau_tdc_l_1 = 0;
   cau_tdc_ref1_l_1 = 0;
   cau_tdc_t_1 = 0;
   cau_tdc_ref1_t_1 = 0;
   cau_tdc_l_2 = 0;
   cau_tdc_ref1_l_2 = 0;
   cau_tdc_t_2 = 0;
   cau_tdc_ref1_t_2 = 0;
   cau_tdc_l_3 = 0;
   cau_tdc_ref1_l_3 = 0;
   cau_tdc_t_3 = 0;
   cau_tdc_ref1_t_3 = 0;
   cau_tdc_l_4 = 0;
   cau_tdc_ref1_l_4 = 0;
   cau_tdc_t_4 = 0;
   cau_tdc_ref1_t_4 = 0;
   cau_tdc_l_5 = 0;
   cau_tdc_ref1_l_5 = 0;
   cau_tdc_t_5 = 0;
   cau_tdc_ref1_t_5 = 0;
   cau_tdc_l_6 = 0;
   cau_tdc_ref1_l_6 = 0;
   cau_tdc_t_6 = 0;
   cau_tdc_ref1_t_6 = 0;
   cau_tdc_l_7 = 0;
   cau_tdc_ref1_l_7 = 0;
   cau_tdc_t_7 = 0;
   cau_tdc_ref1_t_7 = 0;
   cau_tdc_l_8 = 0;
   cau_tdc_ref1_l_8 = 0;
   cau_tdc_t_8 = 0;
   cau_tdc_ref1_t_8 = 0;
   cau_tdc_l_9 = 0;
   cau_tdc_ref1_l_9 = 0;
   cau_tdc_t_9 = 0;
   cau_tdc_ref1_t_9 = 0;
   cau_tdc_l_10 = 0;
   cau_tdc_ref1_l_10 = 0;
   cau_tdc_t_10 = 0;
   cau_tdc_ref1_t_10 = 0;
   cau_tdc_l_11 = 0;
   cau_tdc_ref1_l_11 = 0;
   cau_tdc_t_11 = 0;
   cau_tdc_ref1_t_11 = 0;
   cau_tdc_l_12 = 0;
   cau_tdc_ref1_l_12 = 0;
   cau_tdc_t_12 = 0;
   cau_tdc_ref1_t_12 = 0;
   cau_tdc_l_13 = 0;
   cau_tdc_ref1_l_13 = 0;
   cau_tdc_t_13 = 0;
   cau_tdc_ref1_t_13 = 0;
   cau_tdc_l_14 = 0;
   cau_tdc_ref1_l_14 = 0;
   cau_tdc_t_14 = 0;
   cau_tdc_ref1_t_14 = 0;
   cau_tdc_l_15 = 0;
   cau_tdc_ref1_l_15 = 0;
   cau_tdc_t_15 = 0;
   cau_tdc_ref1_t_15 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("CauNum", &CauNum, &b_CauNum);
   fChain->SetBranchAddress("Cautime", &Cautime, &b_Cautime);
   fChain->SetBranchAddress("select_line", &select_line, &b_select_line);
   fChain->SetBranchAddress("raw_trig_cnt", &raw_trig_cnt, &b_raw_trig_cnt);
   fChain->SetBranchAddress("acpt_trig_cnt", &acpt_trig_cnt, &b_acpt_trig_cnt);
   fChain->SetBranchAddress("cau_status", &cau_status, &b_cau_status);
   fChain->SetBranchAddress("cau_ref1_l", &cau_ref1_l, &b_cau_ref1_l);
   fChain->SetBranchAddress("cau_ref2_l", &cau_ref2_l, &b_cau_ref2_l);
   fChain->SetBranchAddress("cau_ref1_t", &cau_ref1_t, &b_cau_ref1_t);
   fChain->SetBranchAddress("cau_ref2_t", &cau_ref2_t, &b_cau_ref2_t);
   fChain->SetBranchAddress("cau_tdc_l_0", &cau_tdc_l_0, &b_cau_tdc_l_0);
   fChain->SetBranchAddress("cau_tdc_ref1_l_0", &cau_tdc_ref1_l_0, &b_cau_tdc_ref1_l_0);
   fChain->SetBranchAddress("cau_tdc_t_0", &cau_tdc_t_0, &b_cau_tdc_t_0);
   fChain->SetBranchAddress("cau_tdc_ref1_t_0", &cau_tdc_ref1_t_0, &b_cau_tdc_ref1_t_0);
   fChain->SetBranchAddress("cau_tdc_l_1", &cau_tdc_l_1, &b_cau_tdc_l_1);
   fChain->SetBranchAddress("cau_tdc_ref1_l_1", &cau_tdc_ref1_l_1, &b_cau_tdc_ref1_l_1);
   fChain->SetBranchAddress("cau_tdc_t_1", &cau_tdc_t_1, &b_cau_tdc_t_1);
   fChain->SetBranchAddress("cau_tdc_ref1_t_1", &cau_tdc_ref1_t_1, &b_cau_tdc_ref1_t_1);
   fChain->SetBranchAddress("cau_tdc_l_2", &cau_tdc_l_2, &b_cau_tdc_l_2);
   fChain->SetBranchAddress("cau_tdc_ref1_l_2", &cau_tdc_ref1_l_2, &b_cau_tdc_ref1_l_2);
   fChain->SetBranchAddress("cau_tdc_t_2", &cau_tdc_t_2, &b_cau_tdc_t_2);
   fChain->SetBranchAddress("cau_tdc_ref1_t_2", &cau_tdc_ref1_t_2, &b_cau_tdc_ref1_t_2);
   fChain->SetBranchAddress("cau_tdc_l_3", &cau_tdc_l_3, &b_cau_tdc_l_3);
   fChain->SetBranchAddress("cau_tdc_ref1_l_3", &cau_tdc_ref1_l_3, &b_cau_tdc_ref1_l_3);
   fChain->SetBranchAddress("cau_tdc_t_3", &cau_tdc_t_3, &b_cau_tdc_t_3);
   fChain->SetBranchAddress("cau_tdc_ref1_t_3", &cau_tdc_ref1_t_3, &b_cau_tdc_ref1_t_3);
   fChain->SetBranchAddress("cau_tdc_l_4", &cau_tdc_l_4, &b_cau_tdc_l_4);
   fChain->SetBranchAddress("cau_tdc_ref1_l_4", &cau_tdc_ref1_l_4, &b_cau_tdc_ref1_l_4);
   fChain->SetBranchAddress("cau_tdc_t_4", &cau_tdc_t_4, &b_cau_tdc_t_4);
   fChain->SetBranchAddress("cau_tdc_ref1_t_4", &cau_tdc_ref1_t_4, &b_cau_tdc_ref1_t_4);
   fChain->SetBranchAddress("cau_tdc_l_5", &cau_tdc_l_5, &b_cau_tdc_l_5);
   fChain->SetBranchAddress("cau_tdc_ref1_l_5", &cau_tdc_ref1_l_5, &b_cau_tdc_ref1_l_5);
   fChain->SetBranchAddress("cau_tdc_t_5", &cau_tdc_t_5, &b_cau_tdc_t_5);
   fChain->SetBranchAddress("cau_tdc_ref1_t_5", &cau_tdc_ref1_t_5, &b_cau_tdc_ref1_t_5);
   fChain->SetBranchAddress("cau_tdc_l_6", &cau_tdc_l_6, &b_cau_tdc_l_6);
   fChain->SetBranchAddress("cau_tdc_ref1_l_6", &cau_tdc_ref1_l_6, &b_cau_tdc_ref1_l_6);
   fChain->SetBranchAddress("cau_tdc_t_6", &cau_tdc_t_6, &b_cau_tdc_t_6);
   fChain->SetBranchAddress("cau_tdc_ref1_t_6", &cau_tdc_ref1_t_6, &b_cau_tdc_ref1_t_6);
   fChain->SetBranchAddress("cau_tdc_l_7", &cau_tdc_l_7, &b_cau_tdc_l_7);
   fChain->SetBranchAddress("cau_tdc_ref1_l_7", &cau_tdc_ref1_l_7, &b_cau_tdc_ref1_l_7);
   fChain->SetBranchAddress("cau_tdc_t_7", &cau_tdc_t_7, &b_cau_tdc_t_7);
   fChain->SetBranchAddress("cau_tdc_ref1_t_7", &cau_tdc_ref1_t_7, &b_cau_tdc_ref1_t_7);
   fChain->SetBranchAddress("cau_tdc_l_8", &cau_tdc_l_8, &b_cau_tdc_l_8);
   fChain->SetBranchAddress("cau_tdc_ref1_l_8", &cau_tdc_ref1_l_8, &b_cau_tdc_ref1_l_8);
   fChain->SetBranchAddress("cau_tdc_t_8", &cau_tdc_t_8, &b_cau_tdc_t_8);
   fChain->SetBranchAddress("cau_tdc_ref1_t_8", &cau_tdc_ref1_t_8, &b_cau_tdc_ref1_t_8);
   fChain->SetBranchAddress("cau_tdc_l_9", &cau_tdc_l_9, &b_cau_tdc_l_9);
   fChain->SetBranchAddress("cau_tdc_ref1_l_9", &cau_tdc_ref1_l_9, &b_cau_tdc_ref1_l_9);
   fChain->SetBranchAddress("cau_tdc_t_9", &cau_tdc_t_9, &b_cau_tdc_t_9);
   fChain->SetBranchAddress("cau_tdc_ref1_t_9", &cau_tdc_ref1_t_9, &b_cau_tdc_ref1_t_9);
   fChain->SetBranchAddress("cau_tdc_l_10", &cau_tdc_l_10, &b_cau_tdc_l_10);
   fChain->SetBranchAddress("cau_tdc_ref1_l_10", &cau_tdc_ref1_l_10, &b_cau_tdc_ref1_l_10);
   fChain->SetBranchAddress("cau_tdc_t_10", &cau_tdc_t_10, &b_cau_tdc_t_10);
   fChain->SetBranchAddress("cau_tdc_ref1_t_10", &cau_tdc_ref1_t_10, &b_cau_tdc_ref1_t_10);
   fChain->SetBranchAddress("cau_tdc_l_11", &cau_tdc_l_11, &b_cau_tdc_l_11);
   fChain->SetBranchAddress("cau_tdc_ref1_l_11", &cau_tdc_ref1_l_11, &b_cau_tdc_ref1_l_11);
   fChain->SetBranchAddress("cau_tdc_t_11", &cau_tdc_t_11, &b_cau_tdc_t_11);
   fChain->SetBranchAddress("cau_tdc_ref1_t_11", &cau_tdc_ref1_t_11, &b_cau_tdc_ref1_t_11);
   fChain->SetBranchAddress("cau_tdc_l_12", &cau_tdc_l_12, &b_cau_tdc_l_12);
   fChain->SetBranchAddress("cau_tdc_ref1_l_12", &cau_tdc_ref1_l_12, &b_cau_tdc_ref1_l_12);
   fChain->SetBranchAddress("cau_tdc_t_12", &cau_tdc_t_12, &b_cau_tdc_t_12);
   fChain->SetBranchAddress("cau_tdc_ref1_t_12", &cau_tdc_ref1_t_12, &b_cau_tdc_ref1_t_12);
   fChain->SetBranchAddress("cau_tdc_l_13", &cau_tdc_l_13, &b_cau_tdc_l_13);
   fChain->SetBranchAddress("cau_tdc_ref1_l_13", &cau_tdc_ref1_l_13, &b_cau_tdc_ref1_l_13);
   fChain->SetBranchAddress("cau_tdc_t_13", &cau_tdc_t_13, &b_cau_tdc_t_13);
   fChain->SetBranchAddress("cau_tdc_ref1_t_13", &cau_tdc_ref1_t_13, &b_cau_tdc_ref1_t_13);
   fChain->SetBranchAddress("cau_tdc_l_14", &cau_tdc_l_14, &b_cau_tdc_l_14);
   fChain->SetBranchAddress("cau_tdc_ref1_l_14", &cau_tdc_ref1_l_14, &b_cau_tdc_ref1_l_14);
   fChain->SetBranchAddress("cau_tdc_t_14", &cau_tdc_t_14, &b_cau_tdc_t_14);
   fChain->SetBranchAddress("cau_tdc_ref1_t_14", &cau_tdc_ref1_t_14, &b_cau_tdc_ref1_t_14);
   fChain->SetBranchAddress("cau_tdc_l_15", &cau_tdc_l_15, &b_cau_tdc_l_15);
   fChain->SetBranchAddress("cau_tdc_ref1_l_15", &cau_tdc_ref1_l_15, &b_cau_tdc_ref1_l_15);
   fChain->SetBranchAddress("cau_tdc_t_15", &cau_tdc_t_15, &b_cau_tdc_t_15);
   fChain->SetBranchAddress("cau_tdc_ref1_t_15", &cau_tdc_ref1_t_15, &b_cau_tdc_ref1_t_15);
   Notify();
}

Bool_t CauTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void CauTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t CauTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef CauTree_cxx
