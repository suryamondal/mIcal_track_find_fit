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
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
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

  TCanvas c1("c1", "c1", 640, 480);
  // c1.Print(outPDF + "[");
  TPad onePad("onePad", "onePad", 0, 0, 1, 1, kWhite);
  onePad.Range(-1.257796,-45.36986,1.349272,404.3836);
  onePad.SetFillColor(0);
  onePad.SetBorderMode(0);
  onePad.SetBorderSize(2);
  onePad.SetRightMargin(0.1339713);
  onePad.SetFrameBorderMode(0);
  onePad.SetFrameBorderMode(0);
  onePad.SetNumber(1);
  onePad.Draw();

  double sigmaval = 68.;
  TString nhits_cut = "n_nhits > 5";
  TString chi2_cut = "n_chi2/(n_nhits - 5.) < 10.";
  TString mom_cut = "n_momin * n_momout > 0.";

  TString cuts =
    nhits_cut + " && " + chi2_cut;
  
  // Open files
  vector<ROOT::RDF::RResultPtr<TH2D>> mom_cors;
  vector<ROOT::RDF::RResultPtr<TH2D>> mom_response;
  vector<TH1D*> sigmas, means;
  for (int ij=0; ij<int(filesNames.size());ij++) {
    
    std::cout<<" processing file "<<filesNames[ij].Data()<<std::endl;
    
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
      ddd = ddd.Define("n_chi2", [](const ROOT::RVec<Float_t> &chisq) { Float_t tmp = 0.; if(int(chisq.size())) {tmp = chisq[0];} return tmp; }, {"chisq"})
	       .Define("n_momout", [](const ROOT::RVec<Float_t> &trkmm) { Float_t tmp = 0.; if(int(trkmm.size())) {tmp = trkmm[0];} return tmp; }, {"trkmm"})
	       .Define("n_momin", [](const ROOT::RVec<Float_t> &momin) { Float_t tmp = 0.; if(int(momin.size())) {tmp = momin[0];} return tmp; }, {"momin"})
	       .Define("n_nhits", [](const ROOT::RVec<Int_t> &nhits_finder) { Int_t tmp = 0.; if(int(nhits_finder.size())) {tmp = nhits_finder[0];} return tmp; }, {"nhits_finder"});
    } else if (datatypes == "snm") {
      ddd = ddd.Alias("n_chi2","chi2")
	       .Alias("n_momout","momout")
	       .Alias("n_momin","momin")
	       .Alias("n_nhits","nhits_finder");
    }

    if(datatypes == "kalman") {
      cuts = "ntrkt==1 && " + cuts;
    }
    auto ddf = ddd.Filter(cuts.Data(),"chi2ndf_cut");
    
    // auto ddf = ddd.Filter("true","All");
    // if(datatypes == "kalman") {
    //   ddf = ddf.Filter("ntrkt==1","trk cut");
    // }
    // ddf = ddd.Filter(cuts.Data(),"chi2ndf_cut");

    int nbinx   =  40;
    double xmin = -4.;
    double xmax =  4.;
    int nbiny   = 500;
    double ymin = -1.;
    double ymax =  4.;

    auto mom_cor =
      ddf.Filter(mom_cut.Data())
         // .Define("n_p_residual","n_momin - n_momout")
         .Define("n_p_residual", [](Float_t in, Float_t out) {Float_t res = in - out; if(in<0) {res *= -1;} return res;}, {"n_momin", "n_momout"})
         .Histo2D(ROOT::RDF::TH2DModel(TString(titles[ij]).ReplaceAll(" ","_").Data(),
				       titles[ij].Data(), nbinx, xmin, xmax, nbiny, ymin, ymax),
		  "n_momin", "n_p_residual");
    mom_cors.push_back(mom_cor);

    auto mom_res =
      ddf.Histo2D(ROOT::RDF::TH2DModel((TString(titles[ij]).ReplaceAll(" ","_") + "_response").Data(),
				       (titles[ij] + " response").Data(), 80, -4, 4, 80, -4, 4),
		  "n_momout", "n_momin");
    mom_response.push_back(mom_res);

    TH1D *mom_sigma = new TH1D((TString(titles[ij]).ReplaceAll(" ","_") +
				TString::Format("_sigma%.0f",sigmaval)).Data(),
			       (titles[ij] + TString::Format(" #sigma_{%.0f}",sigmaval)).Data(),
			       nbinx, xmin, xmax);
    std::cout<<" name "<<mom_sigma->GetName()<<" "<<mom_sigma->GetTitle()<<std::endl;
    TH1D *mom_mean  = new TH1D((TString(titles[ij]).ReplaceAll(" ","_") + "_mean").Data(),
			       (titles[ij] + TString::Format(" #sigma_{%.0f} mean",sigmaval)).Data(),
			       nbinx, xmin, xmax);
    std::cout<<" name "<<mom_mean->GetName()<<" "<<mom_mean->GetTitle()<<std::endl;
    
    for(int nx=0; nx<nbinx;nx++) {

      TH1D *h = (TH1D*) mom_cors.back()->ProjectionY("",nx+1,nx+1);

      h->SetBinContent(0,0);
      h->SetBinContent(nbiny,0);
      // double hmaxval = h->GetMaximum();
      // for(int ny=0; ny<nbiny;ny++) {
      // 	if(h->GetBinContent(ny+1)<0.01*hmaxval) {
      // 	  h->SetBinContent(ny+1,0);
      // 	}
      // }
      
      // double samples = h->GetEntries();
      double samples = h->GetSumOfWeights();
      if(samples < 50) {continue;}

      double remainder = (100.0 - sigmaval) / 200.0;
      int binLow, binUp;
      double accuLow, accuUp;
      for (accuLow = 0.0, binLow = 0; binLow < h->GetNbinsX() + 1; binLow++)
	if ((accuLow += h->GetBinContent(binLow)) / samples >= remainder)
	  break;
      for (accuUp = 0.0, binUp = h->GetNbinsX() + 1; binUp > binLow; binUp--)
	if ((accuUp += h->GetBinContent(binUp)) / samples >= remainder)
	  break;
      // Interval -> from up edge of binLow to low edge of binUp
      double xLow = h->GetBinLowEdge(binLow + 1);
      double xUp = h->GetBinLowEdge(binUp);
      double xWidth = (xUp - xLow) / 2.0;
      double xCenter = (xUp + xLow) / 2.0;
      double xErr = h->GetBinWidth(1);
      // std::cout<<" nx "<<nx<<" xWidth "<<xWidth<<" xCenter "<<xCenter<<std::endl;

      mom_mean->SetBinContent(nx+1,xCenter);
      mom_sigma->SetBinContent(nx+1,xWidth);
      mom_sigma->SetBinError(nx+1,xErr);

    } // for(int nx=0; nx<nbinx;nx++) {

    sigmas.push_back(mom_sigma);
    means.push_back(mom_mean);

    // if(f) delete f;

  }   // for (int ij=0; ij<int(filesNames.size());ij++) {

  // TFile *outfile = TFile::Open(outPDF.ReplaceAll(".pdf",".root").Data(),"RECREATE");
  // outfile->cd();
  // for (int ij=0; ij<int(sigmas.size());ij++) {
  //   mom_cors[ij]->Write();
  //   sigmas.at(ij)->Write();
  //   means.at(ij)->Write();
  // } // for (int ij=0; ij<int(sigmas.size());ij++) {
  // outfile->Close();


  // Print
  topPad.cd();

  TLegend leg(0.65, 0.85 - 0.06 * (sigmas.size() + 1), 0.93, 0.85);
  leg.SetFillStyle(0);

  int icol = 0;
  for (int i = 0; i < sigmas.size(); i++) {
    icol ++; if(icol == 5) icol++;
    sigmas[i]->SetLineColor(icol);
    // sigmas[i]->SetLineWidth(i == 0 ? 2 : 1);
    sigmas[i]->SetMarkerColor(icol);
    leg.AddEntry(sigmas[i], titles.at(i), "LE");
    sigmas[i]->Draw(!i ? "" : "same");
  }
  leg.Draw();
  topPad.SetGrid();

  btmPad.cd();
  TH1D* hRatio = (TH1D*)sigmas.at(0)->Clone(sigmas.at(0)->GetName() + TString("_ratio"));
  hRatio->Divide(sigmas.at(1));
  // hRatio->GetYaxis()->SetTitle(titles.at(0) + " / " + titles.at(1));
  hRatio->GetYaxis()->SetTitle("black / red");
  hRatio->SetTitle("");
  hRatio->SetMinimum(0.6); hRatio->SetMaximum(1.4);
  hRatio->Draw();
  setFont(hRatio);
  hRatio->GetYaxis()->SetNdivisions(505);
  btmPad.SetGrid();

  c.Print(outPDF, TString("Title:") + sigmas[0]->GetTitle());

  topPad.cd();

  TLegend leg1(0.65, 0.85 - 0.06 * (sigmas.size() + 1), 0.93, 0.85);
  leg1.SetFillStyle(0);

  icol = 0;
  for (int i = 0; i < means.size(); i++) {
    icol ++; if(icol == 5) icol++;
    means[i]->SetLineColor(icol);
    // means[i]->SetLineWidth(i == 0 ? 2 : 1);
    means[i]->SetMarkerColor(icol);
    leg1.AddEntry(means[i], titles.at(i), "LE");
    means[i]->Draw(!i ? "" : "same");
  }
  leg1.Draw();
  topPad.SetGrid();

  btmPad.cd();
  TH1D* hRatio1 = (TH1D*)means.at(0)->Clone(means.at(0)->GetName() + TString("_ratio"));
  hRatio1->Add(hRatio1,means.at(1),1,-1);
  // hRatio1->GetYaxis()->SetTitle(titles.at(0) + " / " + titles.at(1));
  hRatio1->GetYaxis()->SetTitle("black - red");
  hRatio1->SetTitle("");
  // hRatio1->SetMinimum(0.6); hRatio1->SetMaximum(1.4);
  hRatio1->Draw();
  setFont(hRatio1);
  hRatio1->GetYaxis()->SetNdivisions(505);
  btmPad.SetGrid();

  c.Print(outPDF, TString("Title:") + means[0]->GetTitle());


  for (int i = 0; i < mom_cors.size(); i++) {
    onePad.cd();
    mom_cors[i]->Draw("COLZ");
    onePad.SetGrid();
    c1.Print(outPDF, TString("Title:") + mom_cors[i]->GetTitle());
  }

  for (int i = 0; i < mom_response.size(); i++) {
    onePad.cd();
    mom_response[i]->Draw("COLZ");
    onePad.SetGrid();
    c1.Print(outPDF, TString("Title:") + mom_response[i]->GetTitle());
  }



    
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
