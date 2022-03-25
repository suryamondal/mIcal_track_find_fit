#define CauTree_cxx
#include "CauTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void CauTree::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L CauTree.C
//      Root > CauTree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
    cau_tdc_l[0] = cau_tdc_l_0;
    cau_tdc_l[1] = cau_tdc_l_1;
    cau_tdc_l[2] = cau_tdc_l_2;
    cau_tdc_l[3] = cau_tdc_l_3;
    cau_tdc_l[4] = cau_tdc_l_4;
    cau_tdc_l[5] = cau_tdc_l_5;
    cau_tdc_l[6] = cau_tdc_l_6;
    cau_tdc_l[7] = cau_tdc_l_7;
    cau_tdc_l[8] = cau_tdc_l_8;
    cau_tdc_l[9] = cau_tdc_l_9;
    cau_tdc_l[10] = cau_tdc_l_10;
    cau_tdc_l[11] = cau_tdc_l_11;
    cau_tdc_l[12] = cau_tdc_l_12;
    cau_tdc_l[13] = cau_tdc_l_13;
    cau_tdc_l[14] = cau_tdc_l_14;
    cau_tdc_l[15] = cau_tdc_l_15;

    //trailing
    cau_tdc_t[0] = cau_tdc_t_0;
    cau_tdc_t[1] = cau_tdc_t_1;
    cau_tdc_t[2] = cau_tdc_t_2;
    cau_tdc_t[3] = cau_tdc_t_3;
    cau_tdc_t[4] = cau_tdc_t_4;
    cau_tdc_t[5] = cau_tdc_t_5;
    cau_tdc_t[6] = cau_tdc_t_6;
    cau_tdc_t[7] = cau_tdc_t_7;
    cau_tdc_t[8] = cau_tdc_t_8;
    cau_tdc_t[9] = cau_tdc_t_9;
    cau_tdc_t[10] = cau_tdc_t_10;
    cau_tdc_t[11] = cau_tdc_t_11;
    cau_tdc_t[12] = cau_tdc_t_12;
    cau_tdc_t[13] = cau_tdc_t_13;
    cau_tdc_t[14] = cau_tdc_t_14;
    cau_tdc_t[15] = cau_tdc_t_15;

    cau_tdc_ref1_l[0] = cau_tdc_ref1_l_0;
    cau_tdc_ref1_l[1] = cau_tdc_ref1_l_1;
    cau_tdc_ref1_l[2] = cau_tdc_ref1_l_2;
    cau_tdc_ref1_l[3] = cau_tdc_ref1_l_3;
    cau_tdc_ref1_l[4] = cau_tdc_ref1_l_4;
    cau_tdc_ref1_l[5] = cau_tdc_ref1_l_5;
    cau_tdc_ref1_l[6] = cau_tdc_ref1_l_6;
    cau_tdc_ref1_l[7] = cau_tdc_ref1_l_7;
    cau_tdc_ref1_l[8] = cau_tdc_ref1_l_8;
    cau_tdc_ref1_l[9] = cau_tdc_ref1_l_9;
    cau_tdc_ref1_l[10] = cau_tdc_ref1_l_10;
    cau_tdc_ref1_l[11] = cau_tdc_ref1_l_11;
    cau_tdc_ref1_l[12] = cau_tdc_ref1_l_12;
    cau_tdc_ref1_l[13] = cau_tdc_ref1_l_13;
    cau_tdc_ref1_l[14] = cau_tdc_ref1_l_14;
    cau_tdc_ref1_l[15] = cau_tdc_ref1_l_15;

    cau_tdc_ref1_t[0] = cau_tdc_ref1_t_0;
    cau_tdc_ref1_t[1] = cau_tdc_ref1_t_1;
    cau_tdc_ref1_t[2] = cau_tdc_ref1_t_2;
    cau_tdc_ref1_t[3] = cau_tdc_ref1_t_3;
    cau_tdc_ref1_t[4] = cau_tdc_ref1_t_4;
    cau_tdc_ref1_t[5] = cau_tdc_ref1_t_5;
    cau_tdc_ref1_t[6] = cau_tdc_ref1_t_6;
    cau_tdc_ref1_t[7] = cau_tdc_ref1_t_7;
    cau_tdc_ref1_t[8] = cau_tdc_ref1_t_8;
    cau_tdc_ref1_t[9] = cau_tdc_ref1_t_9;
    cau_tdc_ref1_t[10] = cau_tdc_ref1_t_10;
    cau_tdc_ref1_t[11] = cau_tdc_ref1_t_11;
    cau_tdc_ref1_t[12] = cau_tdc_ref1_t_12;
    cau_tdc_ref1_t[13] = cau_tdc_ref1_t_13;
    cau_tdc_ref1_t[14] = cau_tdc_ref1_t_14;
    cau_tdc_ref1_t[15] = cau_tdc_ref1_t_15;

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   
   // Long64_t nbytes = 0, nb = 0;
   // for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //    Long64_t ientry = LoadTree(jentry);
   //    if (ientry < 0) break;
   //    nb = fChain->GetEntry(jentry);   nbytes += nb;
   //    // if (Cut(ientry) < 0) continue;
   // }
}
