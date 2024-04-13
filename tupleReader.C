#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TTree.h>

void setPlottingOptions(TH1F *hist, TString title = "; cos#theta; counts") {

  hist->GetYaxis()->SetNdivisions(505, 4);
  hist->GetXaxis()->SetNdivisions(505, 4);
  hist->GetXaxis()->SetLabelSize(0.045);
  hist->GetYaxis()->SetLabelSize(0.043);
  hist->GetXaxis()->SetTitleSize(0.04);
  hist->GetYaxis()->SetTitleSize(0.04);
  hist->GetYaxis()->SetLabelOffset(0.005);
  hist->GetXaxis()->SetLabelOffset(0.005);
  hist->GetXaxis()->CenterTitle(true);
  hist->GetYaxis()->CenterTitle(true);
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->SetTitle(title);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetTitleFontSize();
  gStyle->SetOptStat(0);
}

TH1F *getAngleDistribution(TTree *tuple,
                           TString condition = "wall > 100 && nring == 1",
                           bool weight = true) {
  // used to calculate cosine of zenith angle and oscillation weights
  Float_t dir[3];
  Float_t oscweight3f;
  // Float_t cp1p2oscweight3f;
  // Float_t cp3p2oscweight3f;
  // Float_t invoscweight3f;
  tuple->SetBranchAddress("dir", dir);
  tuple->SetBranchAddress("oscweight3f", &oscweight3f);
  // tuple->SetBranchAddress("cp1p2oscweight3f", &cp1p2oscweight3f);
  // tuple->SetBranchAddress("cp3p2oscweight3f", &cp3p2oscweight3f);
  // tuple->SetBranchAddress("invoscweight3f", &invoscweight3f);

  // used to set cuts
  // Int_t nring;
  // Int_t ip;
  // Int_t ipnu;
  // Int_t itype;
  // Int_t mode;
  // Float_t wall;
  // Float_t amom;
  // tuple->SetBranchAddress("nring", &nring);
  // tuple->SetBranchAddress("ip", &ip);
  // tuple->SetBranchAddress("ipnu", &ipnu);
  // tuple->SetBranchAddress("itype", &itype);
  // tuple->SetBranchAddress("mode", &mode);
  // tuple->SetBranchAddress("wall", &wall);
  // tuple->SetBranchAddress("amom", &amom);

  TH1F *cosZenith = new TH1F("", condition + "; cos#theta; counts", 10, -1, 1);
  // set cuts on cosine of zenith angle distribution
  Long_t nEntries = tuple->GetEntries(condition);
  for (Long_t i = 0; i < nEntries; i++) {
    tuple->GetEntry(i);
    if (weight) {
      cosZenith->Fill(-dir[2] / std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] +
                                          dir[2] * dir[2]),
                      oscweight3f);
    } else {
      cosZenith->Fill(-dir[2] / std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] +
                                          dir[2] * dir[2]),
                      1);
    }
  }

  return cosZenith;
}

void tupleReader() {

  TString mcFileName = "../expanded_fv/MC/fcmc.sk4.19b.mrbdt-2020.all.root";
  TString dataFileName = "../expanded_fv/data/fcdt.sk4.19b.mrbdt-2020.all.root";

  TFile *mcFile = new TFile(mcFileName, "READ");
  TFile *dataFile = new TFile(dataFileName, "READ");

  TTree *mcTuple, *dataTuple;
  mcFile->GetObject("osc_tuple", mcTuple);
  dataFile->GetObject("osc_tuple", dataTuple);

  ///////////////////////////////////////
  TString eDataCondition = "amom > 400 && wall > 200 && nring == 1 && ip == 2";
  TH1F *eDatahist = getAngleDistribution(dataTuple, eDataCondition);
  setPlottingOptions(eDatahist, "e data; cos#theta; counts"); //" + eDataCondition + "

  TString eMCCondition = 
      "pnu > 0.4 && wall > 200 && nring == 1 && (mode == 1 || mode == -1) && "
      "ip == 2 && (ipnu == 12 || ipnu == -12)"; //&& itype == 59
  TH1F *eMChistWeight = getAngleDistribution(mcTuple, eMCCondition);
  setPlottingOptions(eMChistWeight, "e-like " + eDataCondition + "; cos#theta; counts");
  eMChistWeight->Scale(eDatahist->Integral() / eMChistWeight->Integral());

  TH1F *eMChist = getAngleDistribution(mcTuple, eMCCondition, false);
  setPlottingOptions(eMChist, "e-like; cos#theta; counts");
  eMChist->Scale(eDatahist->Integral() / eMChist->Integral());

  ////////////////////////////////////////
  TString muDataCondition = "amom > 400 && wall > 100 && nring == 1 && ip == 3";
  TH1F *muDatahist = getAngleDistribution(dataTuple, muDataCondition);
  setPlottingOptions(muDatahist, "#mu-like; cos#theta; counts"); //" + muDataCondition + "

  TString muMCCondition =
      "pnu > 0.4 && wall > 100 && nring == 1 && (mode == 1 || mode == -1) && "
      "ip == 3 && (ipnu == 14 || ipnu == -14)"; //&& itype == 63
  TH1F *muMChistWeight = getAngleDistribution(mcTuple, muMCCondition);
  setPlottingOptions(muMChistWeight, "#mu-like " + muDataCondition + "; cos#theta; counts"); //" + muMCCondition +"
  muMChistWeight->Scale(muDatahist->Integral() / muMChistWeight->Integral());

  TH1F *muMChist = getAngleDistribution(mcTuple, muMCCondition, false);
  setPlottingOptions(muMChist, "#mu-like; cos#theta; counts");
  muMChist->Scale(muDatahist->Integral() / muMChist->Integral());
  ////////////////

  TH1F* mcRatio = (TH1F*)muMChist->Clone();
  mcRatio->Divide(eMChist);

  TH1F* dataRatio = (TH1F*)muDatahist->Clone();
  dataRatio->Divide(eDatahist);

  TH1F* doubleRatio = (TH1F*)dataRatio->Clone();
  doubleRatio->Divide(mcRatio);
  setPlottingOptions(doubleRatio, "R = #frac{(N_{#mu}/N_{e})_{data}}{(N_{#mu}/N_{e})_{MC}} ; cos#theta; R");

  TH1F* muRatio = (TH1F*)muDatahist->Clone();
  muRatio->Divide(muMChist);
  setPlottingOptions(muRatio, "#mu -like ; cos#theta; data/MC");

  TH1F* eRatio = (TH1F*)eDatahist->Clone();
  eRatio->Divide(eMChist);
  setPlottingOptions(eRatio, "e-like; cos#theta; data/MC");

  ///////////////////////////////
  TCanvas *canvas = new TCanvas("canvas1", "", 800, 600);
  canvas->Divide(2, 2);

  canvas->cd(1);
  eMChistWeight->SetLineColor(kGreen + 2);
  eMChistWeight->Draw("hist");
  eMChist->SetLineColor(kRed);
  eMChist->Draw("hist, same");
  eDatahist->SetLineColor(kBlack);
  eDatahist->SetMarkerStyle(20);
  eDatahist->SetMarkerSize(0.8);
  eDatahist->Draw("PE, same");

  canvas->cd(2);
  muMChistWeight->SetLineColor(kGreen + 2);
  muMChistWeight->Draw("hist");
  muMChist->SetLineColor(kRed);
  muMChist->Draw("hist, same");
  muDatahist->SetLineColor(kBlack);
  muDatahist->SetMarkerStyle(20);
  muDatahist->SetMarkerSize(0.8);
  muDatahist->Draw("PE, same");

  canvas->cd(3);
  eRatio->Draw("pe");

  canvas->cd(4);
  muRatio->Draw("pe");


}
