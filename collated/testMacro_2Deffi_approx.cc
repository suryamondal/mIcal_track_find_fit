



/* 
   compile with:
   g++ `root-config --cflags` -o testMacro_2Deffi_approx testMacro_2Deffi_approx.cc `root-config --glibs` -lMinuit

   20220105aa: Making cor inefficiency 1 if more than 0.3
   : Making strp efficiency 1 if more than 0.6


   copyHisto() function is to decrease size of the histogram after generating.
   So use c++ for main function, then root prompt to execute the copyHisto function.

   
*/


#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <ctime>
#include <bitset>

#include "TTimeStamp.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TObject.h"
#include "TRandom.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TMinuit.h"
#include "TF1.h"
#include "TF2.h"
// #include "TSpectrum2D.h"


using namespace std;


const int nstrip = 64;
const int nlayer = 10;
const int nside = 2;



int main(int argc, char** argv) {

  /* 
     argv[0] : main
     argv[1] : layer number
     argv[2] : x or y
     argv[3] : 
     argv[4] : 
     argv[5] : 
     argv[6] : 
     argv[7] : 
     argv[8] : 
     argv[9] :
     argv[10] :
  */

  int whichLayer = stoi(argv[1]);
  int whichSide  = stoi(argv[2]);
  
  cout<<" layer "<<whichLayer<<" side "<<whichSide<<endl;

  TFile *f_sim = new TFile("outputrootfiles/RPCv4t_evtraw_201811_20211209as_trg5of8_corr.root","READ");
  TFile *f_sim1 = new TFile("outputrootfiles/RPCv4t_evtraw_201812_20211209as_trg6789_corr.root","READ");
  
  TFile *fileC = new TFile(TString::Format("outputrootfiles/effi_iter/Collated_evtraw_201811_20220105av_trg5of8_%s%i_20220105aa.root",whichSide==0?"x":"y",whichLayer),"recreate");
    
  TH2S *h1[nside][nlayer][nstrip];     // unc
  TH2S *h1Effi[nside][nlayer][nstrip]; // 
  TH2S *h1Effi1[nside][nlayer][nstrip]; //
  TH2S *h1Effi2, *h1Effi3;
  TH2S *h2[nside][nlayer];	       // total fine
  TH2S *h2ut[nside][nlayer];	       // unc fine
  TH2S *h2ct[nside][nlayer];	       // cor fine
  TH2D *h2ueffi[nside][nlayer], *h2ueff1[nside][nlayer]; // unc fine
  TH2D *h2ceffi[nside][nlayer], *h2ceff1[nside][nlayer]; // unc fine
  TH2D *h3[nlayer];			// tot
  TH2D *h3ut[nside][nlayer];		// unc

  bool intialised[2] = {false};
  int maxbins[2][2] = {0};	// [nside][x or y max bin]
  TH2D *imp_response[2];
  
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  for(int nl=0;nl<nlayer;nl++) {
    
    // if(whichLayer!=nl) {continue;}
    
    for(int nj=0;nj<nside;nj++) {
      
      cout<<" nl "<<nl<<" nj "<<nj<<endl;
      
      if(nl>=6) {
	f_sim->cd();
	if(nj==0) {
	  h3[nl] = (TH2D*)f_sim->Get(TString::Format("inefficiency_tot_m0_xr0_yr0_l%i",nl));}
	h3ut[nj][nl] = (TH2D*)f_sim->Get(TString::Format("inefficiency_unc_tot_m0_xr0_yr0_%s%i",nj==0?"x":"y",nl));
	h2[nj][nl] = (TH2S*)f_sim->Get(TString::Format("laymul_2DmTotFine_m0_xr0_yr0_%s%i",nj==0?"x":"y",nl));
	h2ut[nj][nl] = (TH2S*)f_sim->Get(TString::Format("laymul_2DuTotFine_m0_xr0_yr0_%s%i",nj==0?"x":"y",nl));
	h2ct[nj][nl] = (TH2S*)f_sim->Get(TString::Format("laymul_2DcTotFine_m0_xr0_yr0_%s%i",nj==0?"x":"y",nl));
	for(int ns=0;ns<nstrip;ns++) {
	  h1[nj][nl][ns] = (TH2S*)f_sim->Get(TString::Format("laymul_2Dstrp_m0_xr0_yr0_%s%i_s%i",nj==0?"x":"y",nl,ns));}
      } else {
	f_sim1->cd();
	if(nj==0) {
	  h3[nl] = (TH2D*)f_sim1->Get(TString::Format("inefficiency_tot_m0_xr0_yr0_l%i",nl));}
	h3ut[nj][nl] = (TH2D*)f_sim1->Get(TString::Format("inefficiency_unc_tot_m0_xr0_yr0_%s%i",nj==0?"x":"y",nl));
	h2[nj][nl] = (TH2S*)f_sim1->Get(TString::Format("laymul_2DmTotFine_m0_xr0_yr0_%s%i",nj==0?"x":"y",nl));
	h2ut[nj][nl] = (TH2S*)f_sim1->Get(TString::Format("laymul_2DuTotFine_m0_xr0_yr0_%s%i",nj==0?"x":"y",nl));
	h2ct[nj][nl] = (TH2S*)f_sim1->Get(TString::Format("laymul_2DcTotFine_m0_xr0_yr0_%s%i",nj==0?"x":"y",nl));
	for(int ns=0;ns<nstrip;ns++) {
	  h1[nj][nl][ns] = (TH2S*)f_sim1->Get(TString::Format("laymul_2Dstrp_m0_xr0_yr0_%s%i_s%i",nj==0?"x":"y",nl,ns));}
      }
      
      cout<<"got histos"<<endl;
          
      if(nj==0) {
	double maxval = h3[nl]->GetMaximum();
	for(int ijx=0;ijx<h3[nl]->GetNbinsX();ijx++) {
	  for(int ijy=0;ijy<h3[nl]->GetNbinsY();ijy++) {
	    if(h3[nl]->GetBinContent(ijx+1,ijy+1)/maxval<0.1) {
	      h3[nl]->SetBinContent(ijx+1,ijy+1,0);}}}
      }
      
      // if(whichSide!=nj) {continue;}
      cout<<" after continue nl "<<nl<<" nj "<<nj<<endl;
      
      h2ueffi[nj][nl] = new TH2D(TString::Format("laymul_2DuEffiFine_m0_xr0_yr0_%s%i",nj==0?"x":"y",nl),
				 TString::Format("laymul_2DuEffiFine_m0_xr0_yr0_%s%i",nj==0?"x":"y",nl),
				 h2[nj][nl]->GetNbinsX(),
				 h2[nj][nl]->GetXaxis()->GetBinLowEdge(1),
				 h2[nj][nl]->GetXaxis()->GetBinUpEdge(h2[nj][nl]->GetNbinsX()),
				 h2[nj][nl]->GetNbinsY(),
				 h2[nj][nl]->GetYaxis()->GetBinLowEdge(1),
				 h2[nj][nl]->GetYaxis()->GetBinUpEdge(h2[nj][nl]->GetNbinsY())
				 );
      h2ueff1[nj][nl] = new TH2D(TString::Format("laymul_2DuEff1Fine_m0_xr0_yr0_%s%i",nj==0?"x":"y",nl),
				 TString::Format("laymul_2DuEff1Fine_m0_xr0_yr0_%s%i",nj==0?"x":"y",nl),
				 h2[nj][nl]->GetNbinsX(),
				 h2[nj][nl]->GetXaxis()->GetBinLowEdge(1),
				 h2[nj][nl]->GetXaxis()->GetBinUpEdge(h2[nj][nl]->GetNbinsX()),
				 h2[nj][nl]->GetNbinsY(),
				 h2[nj][nl]->GetYaxis()->GetBinLowEdge(1),
				 h2[nj][nl]->GetYaxis()->GetBinUpEdge(h2[nj][nl]->GetNbinsY())
				 );
      
      h2ceffi[nj][nl] = new TH2D(TString::Format("laymul_2DcEffiFine_m0_xr0_yr0_%s%i",nj==0?"x":"y",nl),
				 TString::Format("laymul_2DcEffiFine_m0_xr0_yr0_%s%i",nj==0?"x":"y",nl),
				 h2[nj][nl]->GetNbinsX(),
				 h2[nj][nl]->GetXaxis()->GetBinLowEdge(1),
				 h2[nj][nl]->GetXaxis()->GetBinUpEdge(h2[nj][nl]->GetNbinsX()),
				 h2[nj][nl]->GetNbinsY(),
				 h2[nj][nl]->GetYaxis()->GetBinLowEdge(1),
				 h2[nj][nl]->GetYaxis()->GetBinUpEdge(h2[nj][nl]->GetNbinsY())
				 );
      h2ceff1[nj][nl] = new TH2D(TString::Format("laymul_2DcEff1Fine_m0_xr0_yr0_%s%i",nj==0?"x":"y",nl),
				 TString::Format("laymul_2DcEff1Fine_m0_xr0_yr0_%s%i",nj==0?"x":"y",nl),
				 h2[nj][nl]->GetNbinsX(),
				 h2[nj][nl]->GetXaxis()->GetBinLowEdge(1),
				 h2[nj][nl]->GetXaxis()->GetBinUpEdge(h2[nj][nl]->GetNbinsX()),
				 h2[nj][nl]->GetNbinsY(),
				 h2[nj][nl]->GetYaxis()->GetBinLowEdge(1),
				 h2[nj][nl]->GetYaxis()->GetBinUpEdge(h2[nj][nl]->GetNbinsY())
				 );
      
      {
	double maxval = h2[nj][nl]->GetMaximum();
	for(int ijx=0;ijx<h2[nj][nl]->GetNbinsX();ijx++) {
	  for(int ijy=0;ijy<h2[nj][nl]->GetNbinsY();ijy++) {
	    double tbinc = 1.;
	    if(h2[nj][nl]->GetBinContent(ijx+1,ijy+1)>0) {
	      tbinc = 1.0*h2ut[nj][nl]->GetBinContent(ijx+1,ijy+1)/h2[nj][nl]->GetBinContent(ijx+1,ijy+1);
	    }
	    h2ueff1[nj][nl]->SetBinContent(ijx+1,ijy+1,1.-tbinc);
	    if(tbinc>0.3) {tbinc = 1.;}
	    else {tbinc = 0.;}
	    h2ueffi[nj][nl]->SetBinContent(ijx+1,ijy+1,1.-tbinc);
	    
	    tbinc = 1.;
	    if(h2[nj][nl]->GetBinContent(ijx+1,ijy+1)>0) {
	      tbinc = 1.0*h2ct[nj][nl]->GetBinContent(ijx+1,ijy+1)/h2[nj][nl]->GetBinContent(ijx+1,ijy+1);
	    }
	    h2ceff1[nj][nl]->SetBinContent(ijx+1,ijy+1,1.-tbinc);
	    if(tbinc>0.3) {tbinc = 1.;}
	    else {tbinc = 0.;}
	    h2ceffi[nj][nl]->SetBinContent(ijx+1,ijy+1,1.-tbinc);
	  }}
	fileC->cd();
	h2ueff1[nj][nl]->Write();
	h2ueffi[nj][nl]->Write();
	h2ceff1[nj][nl]->Write();
	h2ceffi[nj][nl]->Write();
      }
      
      h3ut[nj][nl]->Divide(h3[nl]);
      for(int ijx=0;ijx<h3[nl]->GetNbinsX();ijx++) {
	for(int ijy=0;ijy<h3[nl]->GetNbinsY();ijy++) {
	  if(h3ut[nj][nl]->GetBinContent(ijx+1,ijy+1)==0) {
	    h3ut[nj][nl]->SetBinContent(ijx+1,ijy+1,1.);}}}
      // fileC->cd();
      // h3ut[nj][nl]->Write();
      
      // if(nl!=7 || nj!=1) {continue;}
      // cout<<" nl "<<nl<<" nj "<<nj<<endl;
      
      for(int ns=0;ns<nstrip;ns++) {
	
	// if(ns!=32) {continue;}
	// if(ns>15) {continue;}
	cout<<"\t ns "<<ns<<endl;
	
	int binx = h1[nj][nl][ns]->GetNbinsX();
	int biny = h1[nj][nl][ns]->GetNbinsY();
	
	h1Effi[nj][nl][ns] = (TH2S*)h1[nj][nl][ns]->Clone("xxx");
	h1Effi[nj][nl][ns]->SetNameTitle(TString::Format("%s_Effi",h1[nj][nl][ns]->GetName()),
					 TString::Format("%s_Effi",h1[nj][nl][ns]->GetName()));
	h1Effi[nj][nl][ns]->Scale(0.);
	
	h1Effi1[nj][nl][ns] = (TH2S*)h1[nj][nl][ns]->Clone("xxx1");
	h1Effi1[nj][nl][ns]->SetNameTitle(TString::Format("%s_Effi1",h1[nj][nl][ns]->GetName()),
					  TString::Format("%s_Effi1",h1[nj][nl][ns]->GetName()));
	h1Effi1[nj][nl][ns]->Scale(0.);
	
	double maxval = h2[nj][nl]->GetMaximum();
	for(int ijx=0;ijx<binx;ijx++) {
	  for(int ijy=0;ijy<biny;ijy++) {
	    double tbinc0 = h1[nj][nl][ns]->GetBinContent(ijx+1,ijy+1);
	    double xpos = h1[nj][nl][ns]->GetXaxis()->GetBinCenter(ijx+1);
	    double ypos = h1[nj][nl][ns]->GetYaxis()->GetBinCenter(ijy+1);
	    int tbin = h2[nj][nl]->FindBin(xpos,ypos);
	    double tbinc1 = h2[nj][nl]->GetBinContent(tbin);
	    double tbinc2 = h2ceffi[nj][nl]->GetBinContent(tbin);
	    double tbinc3 = 0.;
	    if(tbinc1>0) {tbinc3 = tbinc0/tbinc1;}
	    h1Effi[nj][nl][ns]->SetBinContent(ijx+1,ijy+1,10000.*tbinc3);
	    
	    tbinc3 *= tbinc2;
	    if(tbinc3>0.6) {tbinc3 = 1.;}
	    if(tbinc1/maxval<0.1) {tbinc3 = 0.;}
	    h1Effi1[nj][nl][ns]->SetBinContent(ijx+1,ijy+1,10000.*tbinc3);
	  }}			// for(int ijx=0;ijx<binx;ijx++) {
	
	fileC->cd();
	h1Effi[nj][nl][ns]->Write();
	h1Effi1[nj][nl][ns]->Write();
	
      }	// for(int ns=0;ns<nstrip;ns++) {

    }	// for(int nl=0;nl<nlayer;nl++) {
  }	// for(int nj=0;nj<side;nj++) {
  
  fileC->Purge();
  fileC->Close();

  return 0;
  
}




void copyHisto() {
  
  TFile *f_sim = new TFile("outputrootfiles/Collated_evtraw_201811_20220105av_trg5of8_20220105aa.root","READ");
  
  TFile *fileC = new TFile("outputrootfiles/Collated_evtraw_201811_20220105av_trg5of8_20220105aa_1.root","recreate");
  
  TH2S *h1[nside][nlayer][nstrip];     // unc
  TH2S *h2[nside][nlayer][nstrip];     // unc
  TH2S *h1Effi[nside][nlayer][nstrip]; // 
  TH2S *h1Effi1[nside][nlayer][nstrip]; // 
  
  double xrange[2] = {-2.,2.};
  double yrange1[2] = {-1.,59.};
  double yrange2[2] = {-1.,62.};
  
  for(int nl=0;nl<nlayer;nl++) {
    
    for(int nj=0;nj<nside;nj++) {
      
      // if(nl!=5) {continue;}
      // if(nl!=5 || nj!=1) {continue;}
      cout<<" nl "<<nl<<" nj "<<nj<<endl;
      
      for(int ns=0;ns<nstrip;ns++) {
	
	f_sim->cd();
	
	h2[nj][nl][ns] = (TH2S*)f_sim->Get(TString::Format("laymul_2Dstrp_m0_xr0_yr0_%s%i_s%i_Effi",nj==0?"x":"y",nl,ns));
	h1[nj][nl][ns] = (TH2S*)f_sim->Get(TString::Format("laymul_2Dstrp_m0_xr0_yr0_%s%i_s%i_Effi1",nj==0?"x":"y",nl,ns));
	
	int binx = h1[nj][nl][ns]->GetNbinsX();
	int biny = h1[nj][nl][ns]->GetNbinsY();
	double binwx = h1[nj][nl][ns]->GetXaxis()->GetBinWidth(1);
	double binwy = h1[nj][nl][ns]->GetYaxis()->GetBinWidth(1);
	cout<<" "<<ns<<" "<<binx<<" "<<binwx<<" "<<biny<<" "<<binwy<<endl;
	double meanx = ((TH1D*)h1[nj][nl][ns]->ProjectionX())->GetMean();
	double meany = ((TH1D*)h1[nj][nl][ns]->ProjectionY())->GetMean();
	if(nj==0 && meanx==0) {
	  meanx = ns+0.5;
	} else if(nj==1 && meany==0) {
	  meany = ns+0.5;}
	cout<<"\t mean "<<meanx<<" "<<meany<<endl;
	
	double mnx, mxx, mny, mxy;
	int binx1, biny1;
	binx1 = (nj==0?(xrange[1]-xrange[0]):(yrange1[1]-yrange1[0]))/binwx;
	biny1 = (nj==1?(xrange[1]-xrange[0]):(yrange2[1]-yrange2[0]))/binwy;
	mnx = (nj==0?xrange[0] + meanx:yrange1[0]);
	mny = (nj==1?xrange[0] + meany:yrange2[0]);
	mxx = (nj==0?xrange[1] + meanx:yrange1[1]);
	mxy = (nj==1?xrange[1] + meany:yrange2[1]);
	cout<<"\t "<<binx1<<" "<<mnx<<" "<<mxx<<" "<<biny1<<" "<<mny<<" "<<mxy<<endl;
	
	h1Effi[nj][nl][ns] = new TH2S(TString::Format("laymul_2Dstrp_m0_xr0_yr0_%s%i_s%i_Main",nj==0?"x":"y",nl,ns),
				      TString::Format("laymul_2Dstrp_m0_xr0_yr0_%s%i_s%i_Main",nj==0?"x":"y",nl,ns),
				      binx1,mnx,mxx,
				      biny1,mny,mxy);
	h1Effi1[nj][nl][ns] = new TH2S(TString::Format("laymul_2Dstrp_m0_xr0_yr0_%s%i_s%i",nj==0?"x":"y",nl,ns),
				       TString::Format("laymul_2Dstrp_m0_xr0_yr0_%s%i_s%i",nj==0?"x":"y",nl,ns),
				       binx1,mnx,mxx,
				       biny1,mny,mxy);
	
	int startx = h1[nj][nl][ns]->GetXaxis()->FindBin(h1Effi1[nj][nl][ns]->GetXaxis()->GetBinCenter(1));
	int starty = h1[nj][nl][ns]->GetYaxis()->FindBin(h1Effi1[nj][nl][ns]->GetYaxis()->GetBinCenter(1));
	cout<<"\t "<<startx<<" "<<starty<<endl;
	for(int ijx=0;ijx<binx1;ijx++) {
	  for(int ijy=0;ijy<biny1;ijy++) {
	    h1Effi1[nj][nl][ns]->SetBinContent(ijx+1,ijy+1,
					       h1[nj][nl][ns]->GetBinContent(ijx+1+startx,ijy+1+starty));
	    h1Effi[nj][nl][ns]->SetBinContent(ijx+1,ijy+1,
					      h2[nj][nl][ns]->GetBinContent(ijx+1+startx,ijy+1+starty));
	  }}
	
	fileC->cd();
	h1Effi[nj][nl][ns]->Write();
	h1Effi1[nj][nl][ns]->Write();
      }}}
  fileC->Purge();
  fileC->Close();
  
}

