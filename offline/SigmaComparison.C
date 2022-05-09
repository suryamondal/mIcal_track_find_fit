/** ROOT macro for plotting efficiency plots together.
 *
 */
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TLegend.h>
#include <TColor.h>
#include <TStyle.h>
#include <vector>
#include <stdexcept>
using namespace std;

#define STRNG(x) #x
#define STRNG2(x) STRNG(x)
#define CHECK(assertion) if(!(assertion)) throw runtime_error("At " STRNG2(__FILE__) ":" STRNG2(__LINE__))

const Color_t Palette[] = {kBlack, kRed};
const Int_t NPalette = sizeof(Palette) / sizeof(Color_t);
const auto MyBlue = TColor::GetColor("#348ABD");
const Int_t Font = 63;
const Double_t TSize = 19; // In pixels

vector<TString> Ls(TDirectory* d)
{
  vector<TString> res;
  for (const auto& obj : *(d->GetListOfKeys()))
    res.push_back(obj->GetName());
  return res;
}

void setStyle()
{
  gStyle->SetTextFont(Font);
  gStyle->SetTextSize(TSize);
  gStyle->SetLabelFont(Font, "xyz");
  gStyle->SetTitleFont(Font, "xyz");
  gStyle->SetLabelSize(TSize, "xyz");
  gStyle->SetTitleSize(TSize, "xyz");
  gStyle->SetTitleX(0.5); gStyle->SetTitleY(0.97);
  gStyle->SetTitleW(0.7); gStyle->SetTitleH(0.1);
}

void setFont(TH1* h)
{
  h->SetTitleFont(Font, "XYZ");
  h->SetTitleSize(TSize, "XYZ");
  h->SetLabelFont(Font, "XYZ");
  h->SetLabelSize(TSize, "XYZ");
  h->SetTitleOffset(3, "X");
  h->SetTitleOffset(1.1, "Y");
}

void SigmaComparison(TString outPDF, const vector<TString>& filesNames,
		     const vector<TString>& titles)
{

  setStyle();
  TCanvas c("c", "c", 640, 480);
  c.Print(outPDF + "[");
  // Subpads
  TPad topPad("topPad", "topPad", 0, 0.33, 1, 1, kWhite);
  TPad btmPad("btmPad", "btmPad", 0, 0, 1, 0.33, kWhite);
  topPad.SetMargin(0.1, 0.1, 0, 0.149);
  btmPad.SetMargin(0.1, 0.1, 0.303, 0);
  topPad.SetNumber(1);
  btmPad.SetNumber(2);
  topPad.Draw();
  btmPad.Draw();

  double sigmaval = 68.;
  TString nhits_cut = "nhits_finder > 5";
  TString chi2_cut = "chi2/(nhits_finder - 5.) < 10.";

  TString cuts =
    nhits_cut + " && " + chi2_cut;
  
  // Open files
  vector<TH2D*> mom_cors;
  vector<TH1D*> sigmas, means;
  for (int ij=0; ij<int(filesNames.size());ij++) {
    TFile* f = TFile::Open(filesNames[ij]);
    CHECK(f);
    TString datatypes;
    TTree *t = f->Get<TTree>("T1");
    if(t == nullptr) {t = f->Get<TTree>("MomTree");
      datatypes = "snm";}
    else {datatypes = "kalman";}
    CHECK(t);

    ROOT::RDataFrame rdf(*t);

    auto ddd = rdf.Define("testClmn","2+3");

    if (datatypes == "kalman") {
      ddd = ddd.Alias("chi2","chisq")
	       .Alias("momout","trkmm");
    }
    
    auto ddf = ddd.Filter(cuts.Data());

    auto mom_cor =
      ddd.Histo2D(ROOT::RDF::TH2DModel(TString(titles[ij]).ReplaceAll(" ","_").Data(),
				       titles[ij].Data(), 60, -6., 6., 1200, -6., 6.),
		  "momin", "momout");
    mom_cors.push_back((TH2D*)&(*mom_cor));

    TH1D *mom_sigma = new TH1D((TString(titles[ij]).ReplaceAll(" ","_") + "_sigma_{68}").Data(),
			       (titles[ij] + " sigma_{68}").Data(), 60, -6, 6);
    TH1D *mom_mean  = new TH1D((TString(titles[ij]).ReplaceAll(" ","_") + "_mean").Data(),
			       (titles[ij] + " sigma_{68} mean").Data(), 60, -6, 6);
    
    const int nbinx = mom_cors.back()->GetNbinsX();

    for(int nx=0; nx<nbinx;nx++) {

      TH1D *h = (TH1D*) mom_cors.back()->ProjectionY("",nx+1,nx+1);
      
      const double samples = h->GetEntries();
      const double remainder = (100.0 - sigmaval) / 200.0;
      int binLow, binUp;
      double accuLow, accuUp;
      for (accuLow = 0.0, binLow = 0; binLow < h->GetNbinsX() + 1; binLow++)
	if ((accuLow += h->GetBinContent(binLow)) / samples >= remainder)
	  break;
      for (accuUp = 0.0, binUp = h->GetNbinsX() + 1; binUp > binLow; binUp--)
	if ((accuUp += h->GetBinContent(binUp)) / samples >= remainder)
	  break;
      // Interval -> from up edge of binLow to low edge of binUp
      const double xLow = h->GetBinLowEdge(binLow + 1);
      const double xUp = h->GetBinLowEdge(binUp);
      const double xWidth = (xUp - xLow) / 2.0;
      const double xCenter = (xUp + xLow) / 2.0;
      const double xErr = h->GetBinWidth(1);

      mom_mean->SetBinContent(nx+1,xCenter);
      mom_sigma->SetBinContent(nx+1,xWidth);
      mom_sigma->SetBinError(nx+1,xErr);

    } // for(int nx=0; nx<nbinx;nx++) {

    sigmas.push_back(mom_sigma);
    means.push_back(mom_mean);
    
  }   // for (int ij=0; ij<int(filesNames.size());ij++) {


  TFile *outfile = new TFile(outPDF.ReplaceAll(".pdf",".root"),"reacreate");
  outfile->cd();
  for (int ij=0; ij<int(sigmas.size());ij++) {
    mom_cors[ij]->Write();
    sigmas[ij]->Write();
    means[ij]->Write();
  } // for (int ij=0; ij<int(sigmas.size());ij++) {
  outfile->Close();


  // for (int ij=0; ij<int(sigmas.size());ij++) {

  //   // Print
  //   topPad.cd();
      
  //   TLegend leg(0.65, 0.85 - 0.06 * (hists.size() + 1), 0.93, 0.85);

  //   int icol = 0;
  //   for (int i = 0; i < hists.size(); i++) {
  //     icol ++; if(icol == 5) icol++;
  //     if(!i) {
  // 	TString titlename = hists[i]->GetYaxis()->GetTitle();
  // 	hMC->GetYaxis()->SetTitle(titlename);
  //     }
  //     hists[i]->SetLineColor(icol);
  //     // hists[i]->SetLineWidth(i == 0 ? 2 : 1);
  //     hists[i]->SetMarkerColor(icol);
  //     leg.AddEntry(hists[i], titles.at(i), "LE");
  //     hists[i]->Draw("same");
  //   }
  //   leg.Draw();
  //   topPad.SetGrid();

  //   btmPad.cd();
  //   TH1* hRatio = (TH1*)hists.at(0)->Clone(hists.at(0)->GetName() + TString("_ratio"));
  //   hRatio->Divide(hists.at(1));
  //   // hRatio->GetYaxis()->SetTitle(titles.at(0) + " / " + titles.at(1));
  //   hRatio->GetYaxis()->SetTitle("black / red");
  //   hRatio->SetTitle("");
  //   // hRatio->SetMinimum(0); hRatio->SetMaximum(3.5);
  //   hRatio->Draw();
  //   setFont(hRatio);
  //   hRatio->GetYaxis()->SetNdivisions(505);
  //   btmPad.SetGrid(1, 2);

  //   c.Print(outPDF, TString("Title:") + hists[0]->GetTitle());
    
    
  // } // for (int ij=0; ij<int(mom_cors.size());ij++) {
  

  // // Loop over (known) directories
  // for (const TString& dirName : {"KpiCuts", "K3piCuts"}) {
  //   vector<TDirectory*> dirs;
  //   for (const auto& f : files) {
  //     TDirectory* dir = f->GetDirectory(dirName);
  //     CHECK(dir);
  //     dirs.push_back(dir);
  //   }

  //   // Loop over hist names
  //   for (const auto& hName : Ls(dirs.front())) {
  //     if (!hName.Contains("_eff_") &&
  // 	  !hName.Contains("_purity_") &&
  // 	  !hName.Contains("_sigma")
  // 	  ) continue;
  //     vector<TH1*> hists;
  //     for (const auto& d : dirs) {
  //       TH1* h = d->Get<TH1>(hName);
  //       CHECK(h);
  //       hists.push_back(h);
  //     }

  //     TString hNameMC = hName;
  //     hNameMC = hNameMC.ReplaceAll("_eff_", "_MC_");
  //     hNameMC = hNameMC.ReplaceAll("_purity_", "_All_");
  //     hNameMC = hNameMC.ReplaceAll("_sigma", "_MC");
  //     cout<<hNameMC<<endl;
  //     TH1 *hMC = dirs.front()->Get<TH1>(hNameMC);
  //     CHECK(hMC);

  //     double getMaximum = 0.; 
  //     for (int i = 0; i < hists.size(); i++) {
  // 	double ttmp = hists[i]->GetMaximum();
  // 	if(getMaximum<ttmp) {getMaximum = ttmp;}
  //     }

  //     // Print
  //     topPad.cd();
      
  //     hMC->SetLineWidth(0);
  //     hMC->SetFillColorAlpha(MyBlue, 0.4);
  //     hMC->SetFillStyle(1001);
  //     hMC->Scale(getMaximum * 0.9 / hMC->GetBinContent(hMC->GetMaximumBin()));
  //     hMC->SetMinimum(0.0001); hMC->SetMaximum(getMaximum*1.1);
  //     hMC->Draw("hist");
  //     setFont(hMC);
      
  //     TLegend leg(0.65, 0.85 - 0.06 * (hists.size() + 1), 0.93, 0.85);
  //     leg.AddEntry(hMC, "MC", "F");

  //     int icol = 0;
  //     for (int i = 0; i < hists.size(); i++) {
  // 	icol ++; if(icol == 5) icol++;
  // 	if(!i) {
  // 	  TString titlename = hists[i]->GetYaxis()->GetTitle();
  // 	  hMC->GetYaxis()->SetTitle(titlename);
  // 	}
  //       hists[i]->SetLineColor(icol);
  //       // hists[i]->SetLineWidth(i == 0 ? 2 : 1);
  //       hists[i]->SetMarkerColor(icol);
  //       leg.AddEntry(hists[i], titles.at(i), "LE");
  //       hists[i]->Draw("same");
  //     }
  //     leg.Draw();
  //     topPad.SetGrid();

  //     btmPad.cd();
  //     TH1* hRatio = (TH1*)hists.at(0)->Clone(hists.at(0)->GetName() + TString("_ratio"));
  //     hRatio->Divide(hists.at(1));
  //     // hRatio->GetYaxis()->SetTitle(titles.at(0) + " / " + titles.at(1));
  //     hRatio->GetYaxis()->SetTitle("black / red");
  //     hRatio->SetTitle("");
  //     // hRatio->SetMinimum(0); hRatio->SetMaximum(3.5);
  //     hRatio->Draw();
  //     setFont(hRatio);
  //     hRatio->GetYaxis()->SetNdivisions(505);
  //     btmPad.SetGrid(1, 2);

  //     c.Print(outPDF, TString("Title:") + hists[0]->GetTitle());
  //   }
  // }

  c.Print(outPDF + "]");
}
