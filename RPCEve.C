#define RPCEve_cxx
#include "RPCEve.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void RPCEve::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L RPCEve.C
//      Root > RPCEve t
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

    EveTS[0]=Evetime_0;
    EveTS[1]=Evetime_1;
    EveTS[2]=Evetime_2;
    EveTS[3]=Evetime_3;
    EveTS[4]=Evetime_4;
    EveTS[5]=Evetime_5;
    EveTS[6]=Evetime_6;
    EveTS[7]=Evetime_7;
    EveTS[8]=Evetime_8;
    EveTS[9]=Evetime_9;
    EveTS[10]=Evetime_10;
    EveTS[11]=Evetime_11;


    //start adding here
    xLayer[0]=xstriphitsL0;
    yLayer[0]=ystriphitsL0;
    xLayer[1]=xstriphitsL1;
    yLayer[1]=ystriphitsL1;
    xLayer[2]=xstriphitsL2;
    yLayer[2]=ystriphitsL2;
    xLayer[3]=xstriphitsL3;
    yLayer[3]=ystriphitsL3;
    xLayer[4]=xstriphitsL4;
    yLayer[4]=ystriphitsL4;
    xLayer[5]=xstriphitsL5;
    yLayer[5]=ystriphitsL5;
    xLayer[6]=xstriphitsL6;
    yLayer[6]=ystriphitsL6;
    xLayer[7]=xstriphitsL7;
    yLayer[7]=ystriphitsL7;
    xLayer[8]=xstriphitsL8;
    yLayer[8]=ystriphitsL8;
    xLayer[9]=xstriphitsL9;
    yLayer[9]=ystriphitsL9;
    xLayer[10]=xstriphitsL10;
    yLayer[10]=ystriphitsL10;
    xLayer[11]=xstriphitsL11;
    yLayer[11]=ystriphitsL11;

    //x-side
    vxtdc_l[0][0] = xtdc_l_0_0;     //layer-0 x-side Leading edge
    vxtdc_l[0][1] = xtdc_l_0_1;
    vxtdc_l[0][2] = xtdc_l_0_2;
    vxtdc_l[0][3] = xtdc_l_0_3;
    vxtdc_l[0][4] = xtdc_l_0_4;
    vxtdc_l[0][5] = xtdc_l_0_5;
    vxtdc_l[0][6] = xtdc_l_0_6;
    vxtdc_l[0][7] = xtdc_l_0_7;

    vxtdc_t[0][0] = xtdc_t_0_0;     //layer-0 x-side trailing edge
    vxtdc_t[0][1] = xtdc_t_0_1;
    vxtdc_t[0][2] = xtdc_t_0_2;
    vxtdc_t[0][3] = xtdc_t_0_3;
    vxtdc_t[0][4] = xtdc_t_0_4;
    vxtdc_t[0][5] = xtdc_t_0_5;
    vxtdc_t[0][6] = xtdc_t_0_6;
    vxtdc_t[0][7] = xtdc_t_0_7;

    vxtdc_l[1][0] = xtdc_l_1_0;     //layer-1 x-side Leading edge
    vxtdc_l[1][1] = xtdc_l_1_1;
    vxtdc_l[1][2] = xtdc_l_1_2;
    vxtdc_l[1][3] = xtdc_l_1_3;
    vxtdc_l[1][4] = xtdc_l_1_4;
    vxtdc_l[1][5] = xtdc_l_1_5;
    vxtdc_l[1][6] = xtdc_l_1_6;
    vxtdc_l[1][7] = xtdc_l_1_7;

    vxtdc_t[1][0] = xtdc_t_1_0;     //layer-1 x-side trailing edge
    vxtdc_t[1][1] = xtdc_t_1_1;
    vxtdc_t[1][2] = xtdc_t_1_2;
    vxtdc_t[1][3] = xtdc_t_1_3;
    vxtdc_t[1][4] = xtdc_t_1_4;
    vxtdc_t[1][5] = xtdc_t_1_5;
    vxtdc_t[1][6] = xtdc_t_1_6;
    vxtdc_t[1][7] = xtdc_t_1_7;


    vxtdc_l[2][0] = xtdc_l_2_0;     //layer-2 x-side Leading edge
    vxtdc_l[2][1] = xtdc_l_2_1;
    vxtdc_l[2][2] = xtdc_l_2_2;
    vxtdc_l[2][3] = xtdc_l_2_3;
    vxtdc_l[2][4] = xtdc_l_2_4;
    vxtdc_l[2][5] = xtdc_l_2_5;
    vxtdc_l[2][6] = xtdc_l_2_6;
    vxtdc_l[2][7] = xtdc_l_2_7;

    vxtdc_t[2][0] = xtdc_t_2_0;     //layer-2 x-side trailing edge
    vxtdc_t[2][1] = xtdc_t_2_1;
    vxtdc_t[2][2] = xtdc_t_2_2;
    vxtdc_t[2][3] = xtdc_t_2_3;
    vxtdc_t[2][4] = xtdc_t_2_4;
    vxtdc_t[2][5] = xtdc_t_2_5;
    vxtdc_t[2][6] = xtdc_t_2_6;
    vxtdc_t[2][7] = xtdc_t_2_7;


    vxtdc_l[3][0] = xtdc_l_3_0;     //layer-3 x-side Leading edge
    vxtdc_l[3][1] = xtdc_l_3_1;
    vxtdc_l[3][2] = xtdc_l_3_2;
    vxtdc_l[3][3] = xtdc_l_3_3;
    vxtdc_l[3][4] = xtdc_l_3_4;
    vxtdc_l[3][5] = xtdc_l_3_5;
    vxtdc_l[3][6] = xtdc_l_3_6;
    vxtdc_l[3][7] = xtdc_l_3_7;

    vxtdc_t[3][0] = xtdc_t_3_0;     //layer-3 x-side trailing edge
    vxtdc_t[3][1] = xtdc_t_3_1;
    vxtdc_t[3][2] = xtdc_t_3_2;
    vxtdc_t[3][3] = xtdc_t_3_3;
    vxtdc_t[3][4] = xtdc_t_3_4;
    vxtdc_t[3][5] = xtdc_t_3_5;
    vxtdc_t[3][6] = xtdc_t_3_6;
    vxtdc_t[3][7] = xtdc_t_3_7;


    vxtdc_l[4][0] = xtdc_l_4_0;     //layer-4 x-side Leading edge
    vxtdc_l[4][1] = xtdc_l_4_1;
    vxtdc_l[4][2] = xtdc_l_4_2;
    vxtdc_l[4][3] = xtdc_l_4_3;
    vxtdc_l[4][4] = xtdc_l_4_4;
    vxtdc_l[4][5] = xtdc_l_4_5;
    vxtdc_l[4][6] = xtdc_l_4_6;
    vxtdc_l[4][7] = xtdc_l_4_7;

    vxtdc_t[4][0] = xtdc_t_4_0;     //layer-4 x-side trailing edge
    vxtdc_t[4][1] = xtdc_t_4_1;
    vxtdc_t[4][2] = xtdc_t_4_2;
    vxtdc_t[4][3] = xtdc_t_4_3;
    vxtdc_t[4][4] = xtdc_t_4_4;
    vxtdc_t[4][5] = xtdc_t_4_5;
    vxtdc_t[4][6] = xtdc_t_4_6;
    vxtdc_t[4][7] = xtdc_t_4_7;


    vxtdc_l[5][0] = xtdc_l_5_0;     //layer-5 x-side Leading edge
    vxtdc_l[5][1] = xtdc_l_5_1;
    vxtdc_l[5][2] = xtdc_l_5_2;
    vxtdc_l[5][3] = xtdc_l_5_3;
    vxtdc_l[5][4] = xtdc_l_5_4;
    vxtdc_l[5][5] = xtdc_l_5_5;
    vxtdc_l[5][6] = xtdc_l_5_6;
    vxtdc_l[5][7] = xtdc_l_5_7;

    vxtdc_t[5][0] = xtdc_t_5_0;     //layer-5 x-side trailing edge
    vxtdc_t[5][1] = xtdc_t_5_1;
    vxtdc_t[5][2] = xtdc_t_5_2;
    vxtdc_t[5][3] = xtdc_t_5_3;
    vxtdc_t[5][4] = xtdc_t_5_4;
    vxtdc_t[5][5] = xtdc_t_5_5;
    vxtdc_t[5][6] = xtdc_t_5_6;
    vxtdc_t[5][7] = xtdc_t_5_7;


    vxtdc_l[6][0] = xtdc_l_6_0;     //layer-6 x-side Leading edge
    vxtdc_l[6][1] = xtdc_l_6_1;
    vxtdc_l[6][2] = xtdc_l_6_2;
    vxtdc_l[6][3] = xtdc_l_6_3;
    vxtdc_l[6][4] = xtdc_l_6_4;
    vxtdc_l[6][5] = xtdc_l_6_5;
    vxtdc_l[6][6] = xtdc_l_6_6;
    vxtdc_l[6][7] = xtdc_l_6_7;

    vxtdc_t[6][0] = xtdc_t_6_0;     //layer-6 x-side trailing edge
    vxtdc_t[6][1] = xtdc_t_6_1;
    vxtdc_t[6][2] = xtdc_t_6_2;
    vxtdc_t[6][3] = xtdc_t_6_3;
    vxtdc_t[6][4] = xtdc_t_6_4;
    vxtdc_t[6][5] = xtdc_t_6_5;
    vxtdc_t[6][6] = xtdc_t_6_6;
    vxtdc_t[6][7] = xtdc_t_6_7;


    vxtdc_l[7][0] = xtdc_l_7_0;     //layer-7 x-side Leading edge
    vxtdc_l[7][1] = xtdc_l_7_1;
    vxtdc_l[7][2] = xtdc_l_7_2;
    vxtdc_l[7][3] = xtdc_l_7_3;
    vxtdc_l[7][4] = xtdc_l_7_4;
    vxtdc_l[7][5] = xtdc_l_7_5;
    vxtdc_l[7][6] = xtdc_l_7_6;
    vxtdc_l[7][7] = xtdc_l_7_7;

    vxtdc_t[7][0] = xtdc_t_7_0;     //layer-7 x-side trailing edge
    vxtdc_t[7][1] = xtdc_t_7_1;
    vxtdc_t[7][2] = xtdc_t_7_2;
    vxtdc_t[7][3] = xtdc_t_7_3;
    vxtdc_t[7][4] = xtdc_t_7_4;
    vxtdc_t[7][5] = xtdc_t_7_5;
    vxtdc_t[7][6] = xtdc_t_7_6;
    vxtdc_t[7][7] = xtdc_t_7_7;


    vxtdc_l[8][0] = xtdc_l_8_0;     //layer-8 x-side Leading edge
    vxtdc_l[8][1] = xtdc_l_8_1;
    vxtdc_l[8][2] = xtdc_l_8_2;
    vxtdc_l[8][3] = xtdc_l_8_3;
    vxtdc_l[8][4] = xtdc_l_8_4;
    vxtdc_l[8][5] = xtdc_l_8_5;
    vxtdc_l[8][6] = xtdc_l_8_6;
    vxtdc_l[8][7] = xtdc_l_8_7;

    vxtdc_t[8][0] = xtdc_t_8_0;     //layer-8 x-side trailing edge
    vxtdc_t[8][1] = xtdc_t_8_1;
    vxtdc_t[8][2] = xtdc_t_8_2;
    vxtdc_t[8][3] = xtdc_t_8_3;
    vxtdc_t[8][4] = xtdc_t_8_4;
    vxtdc_t[8][5] = xtdc_t_8_5;
    vxtdc_t[8][6] = xtdc_t_8_6;
    vxtdc_t[8][7] = xtdc_t_8_7;


    vxtdc_l[9][0] = xtdc_l_9_0;     //layer-9 x-side Leading edge
    vxtdc_l[9][1] = xtdc_l_9_1;
    vxtdc_l[9][2] = xtdc_l_9_2;
    vxtdc_l[9][3] = xtdc_l_9_3;
    vxtdc_l[9][4] = xtdc_l_9_4;
    vxtdc_l[9][5] = xtdc_l_9_5;
    vxtdc_l[9][6] = xtdc_l_9_6;
    vxtdc_l[9][7] = xtdc_l_9_7;

    vxtdc_t[9][0] = xtdc_t_9_0;     //layer-9 x-side trailing edge
    vxtdc_t[9][1] = xtdc_t_9_1;
    vxtdc_t[9][2] = xtdc_t_9_2;
    vxtdc_t[9][3] = xtdc_t_9_3;
    vxtdc_t[9][4] = xtdc_t_9_4;
    vxtdc_t[9][5] = xtdc_t_9_5;
    vxtdc_t[9][6] = xtdc_t_9_6;
    vxtdc_t[9][7] = xtdc_t_9_7;


    vxtdc_l[10][0] = xtdc_l_10_0;     //layer-10 x-side Leading edge
    vxtdc_l[10][1] = xtdc_l_10_1;
    vxtdc_l[10][2] = xtdc_l_10_2;
    vxtdc_l[10][3] = xtdc_l_10_3;
    vxtdc_l[10][4] = xtdc_l_10_4;
    vxtdc_l[10][5] = xtdc_l_10_5;
    vxtdc_l[10][6] = xtdc_l_10_6;
    vxtdc_l[10][7] = xtdc_l_10_7;

    vxtdc_t[10][0] = xtdc_t_10_0;     //layer-10 x-side trailing edge
    vxtdc_t[10][1] = xtdc_t_10_1;
    vxtdc_t[10][2] = xtdc_t_10_2;
    vxtdc_t[10][3] = xtdc_t_10_3;
    vxtdc_t[10][4] = xtdc_t_10_4;
    vxtdc_t[10][5] = xtdc_t_10_5;
    vxtdc_t[10][6] = xtdc_t_10_6;
    vxtdc_t[10][7] = xtdc_t_10_7;

    vxtdc_l[11][0] = xtdc_l_11_0;     //layer-11 x-side Leading edge
    vxtdc_l[11][1] = xtdc_l_11_1;
    vxtdc_l[11][2] = xtdc_l_11_2;
    vxtdc_l[11][3] = xtdc_l_11_3;
    vxtdc_l[11][4] = xtdc_l_11_4;
    vxtdc_l[11][5] = xtdc_l_11_5;
    vxtdc_l[11][6] = xtdc_l_11_6;
    vxtdc_l[11][7] = xtdc_l_11_7;

    vxtdc_t[11][0] = xtdc_t_11_0;     //layer-11 x-side trailing edge
    vxtdc_t[11][1] = xtdc_t_11_1;
    vxtdc_t[11][2] = xtdc_t_11_2;
    vxtdc_t[11][3] = xtdc_t_11_3;
    vxtdc_t[11][4] = xtdc_t_11_4;
    vxtdc_t[11][5] = xtdc_t_11_5;
    vxtdc_t[11][6] = xtdc_t_11_6;
    vxtdc_t[11][7] = xtdc_t_11_7;

    //Y-side

    vytdc_l[0][0] = ytdc_l_0_0;     //layer-0 x-side Leading edge
    vytdc_l[0][1] = ytdc_l_0_1;
    vytdc_l[0][2] = ytdc_l_0_2;
    vytdc_l[0][3] = ytdc_l_0_3;
    vytdc_l[0][4] = ytdc_l_0_4;
    vytdc_l[0][5] = ytdc_l_0_5;
    vytdc_l[0][6] = ytdc_l_0_6;
    vytdc_l[0][7] = ytdc_l_0_7;

    vytdc_t[0][0] = ytdc_t_0_0;     //layer-0 x-side trailing edge
    vytdc_t[0][1] = ytdc_t_0_1;
    vytdc_t[0][2] = ytdc_t_0_2;
    vytdc_t[0][3] = ytdc_t_0_3;
    vytdc_t[0][4] = ytdc_t_0_4;
    vytdc_t[0][5] = ytdc_t_0_5;
    vytdc_t[0][6] = ytdc_t_0_6;
    vytdc_t[0][7] = ytdc_t_0_7;

    vytdc_l[1][0] = ytdc_l_1_0;     //layer-1 x-side Leading edge
    vytdc_l[1][1] = ytdc_l_1_1;
    vytdc_l[1][2] = ytdc_l_1_2;
    vytdc_l[1][3] = ytdc_l_1_3;
    vytdc_l[1][4] = ytdc_l_1_4;
    vytdc_l[1][5] = ytdc_l_1_5;
    vytdc_l[1][6] = ytdc_l_1_6;
    vytdc_l[1][7] = ytdc_l_1_7;

    vytdc_t[1][0] = ytdc_t_1_0;     //layer-1 x-side trailing edge
    vytdc_t[1][1] = ytdc_t_1_1;
    vytdc_t[1][2] = ytdc_t_1_2;
    vytdc_t[1][3] = ytdc_t_1_3;
    vytdc_t[1][4] = ytdc_t_1_4;
    vytdc_t[1][5] = ytdc_t_1_5;
    vytdc_t[1][6] = ytdc_t_1_6;
    vytdc_t[1][7] = ytdc_t_1_7;


    vytdc_l[2][0] = ytdc_l_2_0;     //layer-2 x-side Leading edge
    vytdc_l[2][1] = ytdc_l_2_1;
    vytdc_l[2][2] = ytdc_l_2_2;
    vytdc_l[2][3] = ytdc_l_2_3;
    vytdc_l[2][4] = ytdc_l_2_4;
    vytdc_l[2][5] = ytdc_l_2_5;
    vytdc_l[2][6] = ytdc_l_2_6;
    vytdc_l[2][7] = ytdc_l_2_7;

    vytdc_t[2][0] = ytdc_t_2_0;     //layer-2 x-side trailing edge
    vytdc_t[2][1] = ytdc_t_2_1;
    vytdc_t[2][2] = ytdc_t_2_2;
    vytdc_t[2][3] = ytdc_t_2_3;
    vytdc_t[2][4] = ytdc_t_2_4;
    vytdc_t[2][5] = ytdc_t_2_5;
    vytdc_t[2][6] = ytdc_t_2_6;
    vytdc_t[2][7] = ytdc_t_2_7;


    vytdc_l[3][0] = ytdc_l_3_0;     //layer-3 x-side Leading edge
    vytdc_l[3][1] = ytdc_l_3_1;
    vytdc_l[3][2] = ytdc_l_3_2;
    vytdc_l[3][3] = ytdc_l_3_3;
    vytdc_l[3][4] = ytdc_l_3_4;
    vytdc_l[3][5] = ytdc_l_3_5;
    vytdc_l[3][6] = ytdc_l_3_6;
    vytdc_l[3][7] = ytdc_l_3_7;

    vytdc_t[3][0] = ytdc_t_3_0;     //layer-3 x-side trailing edge
    vytdc_t[3][1] = ytdc_t_3_1;
    vytdc_t[3][2] = ytdc_t_3_2;
    vytdc_t[3][3] = ytdc_t_3_3;
    vytdc_t[3][4] = ytdc_t_3_4;
    vytdc_t[3][5] = ytdc_t_3_5;
    vytdc_t[3][6] = ytdc_t_3_6;
    vytdc_t[3][7] = ytdc_t_3_7;


    vytdc_l[4][0] = ytdc_l_4_0;     //layer-4 x-side Leading edge
    vytdc_l[4][1] = ytdc_l_4_1;
    vytdc_l[4][2] = ytdc_l_4_2;
    vytdc_l[4][3] = ytdc_l_4_3;
    vytdc_l[4][4] = ytdc_l_4_4;
    vytdc_l[4][5] = ytdc_l_4_5;
    vytdc_l[4][6] = ytdc_l_4_6;
    vytdc_l[4][7] = ytdc_l_4_7;

    vytdc_t[4][0] = ytdc_t_4_0;     //layer-4 x-side trailing edge
    vytdc_t[4][1] = ytdc_t_4_1;
    vytdc_t[4][2] = ytdc_t_4_2;
    vytdc_t[4][3] = ytdc_t_4_3;
    vytdc_t[4][4] = ytdc_t_4_4;
    vytdc_t[4][5] = ytdc_t_4_5;
    vytdc_t[4][6] = ytdc_t_4_6;
    vytdc_t[4][7] = ytdc_t_4_7;


    vytdc_l[5][0] = ytdc_l_5_0;     //layer-5 x-side Leading edge
    vytdc_l[5][1] = ytdc_l_5_1;
    vytdc_l[5][2] = ytdc_l_5_2;
    vytdc_l[5][3] = ytdc_l_5_3;
    vytdc_l[5][4] = ytdc_l_5_4;
    vytdc_l[5][5] = ytdc_l_5_5;
    vytdc_l[5][6] = ytdc_l_5_6;
    vytdc_l[5][7] = ytdc_l_5_7;

    vytdc_t[5][0] = ytdc_t_5_0;     //layer-5 x-side trailing edge
    vytdc_t[5][1] = ytdc_t_5_1;
    vytdc_t[5][2] = ytdc_t_5_2;
    vytdc_t[5][3] = ytdc_t_5_3;
    vytdc_t[5][4] = ytdc_t_5_4;
    vytdc_t[5][5] = ytdc_t_5_5;
    vytdc_t[5][6] = ytdc_t_5_6;
    vytdc_t[5][7] = ytdc_t_5_7;


    vytdc_l[6][0] = ytdc_l_6_0;     //layer-6 x-side Leading edge
    vytdc_l[6][1] = ytdc_l_6_1;
    vytdc_l[6][2] = ytdc_l_6_2;
    vytdc_l[6][3] = ytdc_l_6_3;
    vytdc_l[6][4] = ytdc_l_6_4;
    vytdc_l[6][5] = ytdc_l_6_5;
    vytdc_l[6][6] = ytdc_l_6_6;
    vytdc_l[6][7] = ytdc_l_6_7;

    vytdc_t[6][0] = ytdc_t_6_0;     //layer-6 x-side trailing edge
    vytdc_t[6][1] = ytdc_t_6_1;
    vytdc_t[6][2] = ytdc_t_6_2;
    vytdc_t[6][3] = ytdc_t_6_3;
    vytdc_t[6][4] = ytdc_t_6_4;
    vytdc_t[6][5] = ytdc_t_6_5;
    vytdc_t[6][6] = ytdc_t_6_6;
    vytdc_t[6][7] = ytdc_t_6_7;


    vytdc_l[7][0] = ytdc_l_7_0;     //layer-7 x-side Leading edge
    vytdc_l[7][1] = ytdc_l_7_1;
    vytdc_l[7][2] = ytdc_l_7_2;
    vytdc_l[7][3] = ytdc_l_7_3;
    vytdc_l[7][4] = ytdc_l_7_4;
    vytdc_l[7][5] = ytdc_l_7_5;
    vytdc_l[7][6] = ytdc_l_7_6;
    vytdc_l[7][7] = ytdc_l_7_7;

    vytdc_t[7][0] = ytdc_t_7_0;     //layer-7 x-side trailing edge
    vytdc_t[7][1] = ytdc_t_7_1;
    vytdc_t[7][2] = ytdc_t_7_2;
    vytdc_t[7][3] = ytdc_t_7_3;
    vytdc_t[7][4] = ytdc_t_7_4;
    vytdc_t[7][5] = ytdc_t_7_5;
    vytdc_t[7][6] = ytdc_t_7_6;
    vytdc_t[7][7] = ytdc_t_7_7;


    vytdc_l[8][0] = ytdc_l_8_0;     //layer-8 x-side Leading edge
    vytdc_l[8][1] = ytdc_l_8_1;
    vytdc_l[8][2] = ytdc_l_8_2;
    vytdc_l[8][3] = ytdc_l_8_3;
    vytdc_l[8][4] = ytdc_l_8_4;
    vytdc_l[8][5] = ytdc_l_8_5;
    vytdc_l[8][6] = ytdc_l_8_6;
    vytdc_l[8][7] = ytdc_l_8_7;

    vytdc_t[8][0] = ytdc_t_8_0;     //layer-8 x-side trailing edge
    vytdc_t[8][1] = ytdc_t_8_1;
    vytdc_t[8][2] = ytdc_t_8_2;
    vytdc_t[8][3] = ytdc_t_8_3;
    vytdc_t[8][4] = ytdc_t_8_4;
    vytdc_t[8][5] = ytdc_t_8_5;
    vytdc_t[8][6] = ytdc_t_8_6;
    vytdc_t[8][7] = ytdc_t_8_7;


    vytdc_l[9][0] = ytdc_l_9_0;     //layer-9 x-side Leading edge
    vytdc_l[9][1] = ytdc_l_9_1;
    vytdc_l[9][2] = ytdc_l_9_2;
    vytdc_l[9][3] = ytdc_l_9_3;
    vytdc_l[9][4] = ytdc_l_9_4;
    vytdc_l[9][5] = ytdc_l_9_5;
    vytdc_l[9][6] = ytdc_l_9_6;
    vytdc_l[9][7] = ytdc_l_9_7;

    vytdc_t[9][0] = ytdc_t_9_0;     //layer-9 x-side trailing edge
    vytdc_t[9][1] = ytdc_t_9_1;
    vytdc_t[9][2] = ytdc_t_9_2;
    vytdc_t[9][3] = ytdc_t_9_3;
    vytdc_t[9][4] = ytdc_t_9_4;
    vytdc_t[9][5] = ytdc_t_9_5;
    vytdc_t[9][6] = ytdc_t_9_6;
    vytdc_t[9][7] = ytdc_t_9_7;


    vytdc_l[10][0] = ytdc_l_10_0;     //layer-10 x-side Leading edge
    vytdc_l[10][1] = ytdc_l_10_1;
    vytdc_l[10][2] = ytdc_l_10_2;
    vytdc_l[10][3] = ytdc_l_10_3;
    vytdc_l[10][4] = ytdc_l_10_4;
    vytdc_l[10][5] = ytdc_l_10_5;
    vytdc_l[10][6] = ytdc_l_10_6;
    vytdc_l[10][7] = ytdc_l_10_7;

    vytdc_t[10][0] = ytdc_t_10_0;     //layer-10 x-side trailing edge
    vytdc_t[10][1] = ytdc_t_10_1;
    vytdc_t[10][2] = ytdc_t_10_2;
    vytdc_t[10][3] = ytdc_t_10_3;
    vytdc_t[10][4] = ytdc_t_10_4;
    vytdc_t[10][5] = ytdc_t_10_5;
    vytdc_t[10][6] = ytdc_t_10_6;
    vytdc_t[10][7] = ytdc_t_10_7;

    vytdc_l[11][0] = ytdc_l_11_0;     //layer-11 x-side Leading edge
    vytdc_l[11][1] = ytdc_l_11_1;
    vytdc_l[11][2] = ytdc_l_11_2;
    vytdc_l[11][3] = ytdc_l_11_3;
    vytdc_l[11][4] = ytdc_l_11_4;
    vytdc_l[11][5] = ytdc_l_11_5;
    vytdc_l[11][6] = ytdc_l_11_6;
    vytdc_l[11][7] = ytdc_l_11_7;

    vytdc_t[11][0] = ytdc_t_11_0;     //layer-11 x-side trailing edge
    vytdc_t[11][1] = ytdc_t_11_1;
    vytdc_t[11][2] = ytdc_t_11_2;
    vytdc_t[11][3] = ytdc_t_11_3;
    vytdc_t[11][4] = ytdc_t_11_4;
    vytdc_t[11][5] = ytdc_t_11_5;
    vytdc_t[11][6] = ytdc_t_11_6;
    vytdc_t[11][7] = ytdc_t_11_7;




    //End of adding
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
