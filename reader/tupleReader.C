#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTree.h>
#include <TTreeFormula.h>
#include <TLegend.h>
#include <TLatex.h>
// #include <iostream>

void setPlottingOptions(TH1F *hist, TString title = "; cos#theta; events") {

  hist->GetYaxis()->SetNdivisions(505, 4);
  hist->GetXaxis()->SetNdivisions(505, 4);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitleSize(0.2);
  hist->GetYaxis()->SetTitleSize(0.2);
  hist->GetYaxis()->SetLabelOffset(0.003);
  hist->GetXaxis()->SetLabelOffset(0.004);
  hist->GetXaxis()->CenterTitle(true);
  hist->GetYaxis()->CenterTitle(true);
  hist->GetXaxis()->SetTitleOffset(0.4);
  hist->GetYaxis()->SetTitleOffset(0.8);
  hist->SetTitle(title);
  hist->SetTitleOffset(2);
  hist->SetMarkerStyle(4);
  hist->SetMarkerSize(0.8);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetTitleFontSize();
  gStyle->SetOptStat(0);
}

TH1F *getAngleDistribution(TTree *tuple,
                           TString condition = "wall > 200 && nring == 1",
                           bool weight = true) {
  // use to calculate cosine of zenith angle and oscillation weights
  Float_t dir[3];
  Float_t oscweight3f;
  Float_t weightx;
  tuple->SetBranchAddress("dir", dir);
  tuple->SetBranchAddress("oscweight3f", &oscweight3f);
  tuple->SetBranchAddress("weightx", &weightx);

  TH1F *cosZenith = new TH1F("", condition + "; cos#theta; events", 10, -1, 1);

  // set cuts on cosine of zenith angle distribution
  TTreeFormula *formula = new TTreeFormula("", condition, tuple);
  Long64_t nEntries = tuple->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    tuple->GetEntry(i);
    //  evaluate the cut condition using the TTreeFormula
    if (formula->EvalInstance()) {
      if (weight) {
        cosZenith->Fill(-dir[2] / std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] +
                                            dir[2] * dir[2]),
                        weightx * oscweight3f);
      } else {
        cosZenith->Fill(-dir[2] / std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] +
                                            dir[2] * dir[2]),
                        1);
      }
    }
  }

  return cosZenith;
}

void plotOnePad(TH1F *histData, TH1F *histMC, TH1F *histMCWeight, TString momentumString = "single - ring events") {

  TPad *subpad = new TPad("subpad", "subpad", 0, 0.2, 1, 1);
  subpad->Draw();
  subpad->cd();
  subpad->SetTopMargin(0.055);
  subpad->SetRightMargin(0.1);
  subpad->SetLeftMargin(0.1);
  subpad->SetBottomMargin(0);

  histMCWeight->SetMinimum(histData->GetMinimum() * 0.9);
  histMCWeight->SetMaximum(histData->GetMaximum() * 1.24);
  histMCWeight->SetLineColor(kGreen + 2);
  histMCWeight->Draw("hist");

  histMC->SetMaximum(histData->GetMaximum() * 1.12);
  histMC->SetLineColor(kRed);
  histMC->Draw("hist, same");

  histData->SetLineColor(kBlack);
  histData->SetMarkerStyle(20);
  histData->Draw("PE, same");

  TLatex latex ;
  latex.SetNDC ();
  latex.SetTextSize (0.06);
  latex.SetTextColor(kBlue + 4);
  latex.DrawText (0.6 ,0.83 , momentumString );
}

void plotOneRatio(TH1F *histData, TH1F *histMC, TH1F *histMCWeight,
                  TString title = "; cos#theta ; data/MC") { 

  TPad *subpadRatio = new TPad("subpadRatio", "subpadRatio", 0, 0, 1, 0.2);
  subpadRatio->Draw();
  subpadRatio->cd();
  subpadRatio->SetTopMargin(0);
  subpadRatio->SetRightMargin(0.1);
  subpadRatio->SetLeftMargin(0.1);
  subpadRatio->SetBottomMargin(0.3);

  TH1F *ratio = (TH1F *)histData->Clone();
  ratio->Divide(histMC);
  setPlottingOptions(ratio, title);

  TH1F *ratioWeight = (TH1F *)histData->Clone();
  ratioWeight->Divide(histMCWeight);
  setPlottingOptions(ratioWeight, title);

  ratio->SetLineColor(kRed);
  ratio->GetYaxis()->SetRangeUser(0.5, 1.28);
  ratio->GetXaxis()->SetLabelSize(0.2);
  ratio->GetYaxis()->SetLabelSize(0.2);
  ratio->GetYaxis()->SetTitleSize(0.1);
  ratio->GetYaxis()->SetTitleOffset(0.4);
  ratio->GetXaxis()->SetTitleOffset(0.65);
  ratio->GetYaxis()->SetNdivisions(606, 3);
  ratio->Draw("pe");

  ratioWeight->SetLineColor(kGreen + 2);
  ratioWeight->SetMarkerColor(kGreen + 2);
  ratioWeight->GetYaxis()->SetRangeUser(0.5, 1.28);
  ratioWeight->GetXaxis()->SetLabelSize(0.1);
  ratioWeight->GetYaxis()->SetLabelSize(0.1);
  ratioWeight->GetYaxis()->SetTitleSize(0.1);
  ratioWeight->GetYaxis()->SetTitleOffset(0.4);
  ratioWeight->GetXaxis()->SetTitleOffset(0.6);
  ratioWeight->GetYaxis()->SetNdivisions(505, 2);
  ratioWeight->Draw("hist, same");
}

void tupleReader() {

  TString mcFileName = "../expanded_fv/MC/fcmc.sk4.19b.mrbdt-2020.all.root";
  TString dataFileName = "../expanded_fv/data/fcdt.sk4.19b.mrbdt-2020.all.root";

  TFile *mcFile = new TFile(mcFileName, "READ");
  TFile *dataFile = new TFile(dataFileName, "READ");

  TTree *mcTuple, *dataTuple;
  mcFile->GetObject("osc_tuple", mcTuple);
  dataFile->GetObject("osc_tuple", dataTuple);

  double mcLivetime = 182625.0; // eq 500years sk4 // 1851250
  double dataLivetime = 3244.4; // sk4 fcdata // 30621
  double livetimeRatio = (double)1.07 * dataLivetime / mcLivetime;

  TString eCondition = "wall > 200 && nring == 1 && ip == 2 ";
  TString muCondition = "wall > 200 && nring == 1 && ip == 3 ";

  ///////////////////////////////////////

  TH1F *eDatahist =
      getAngleDistribution(dataTuple, eCondition + "&& amom < 400", false);
      setPlottingOptions(eDatahist, "; ; ");

  TH1F *eMChistWeight =
      getAngleDistribution(mcTuple, eCondition + "&& amom < 400");
      setPlottingOptions(eMChistWeight, "e-like; ; events");
      eMChistWeight->GetYaxis()->SetTitleSize(0.05);
  eMChistWeight->Scale(livetimeRatio);

  TH1F *eMChist =
      getAngleDistribution(mcTuple, eCondition + "&& amom < 400", false);
      setPlottingOptions(eMChist, "; ; ");
  eMChist->Scale(livetimeRatio);

  ////////////////////////////////////////
  TH1F *muDatahist =
      getAngleDistribution(dataTuple, muCondition + "&& amom < 400", false);
      setPlottingOptions(muDatahist, "; ; ");

  TH1F *muMChistWeight =
      getAngleDistribution(mcTuple, muCondition + "&& amom < 400");
      setPlottingOptions(muMChistWeight, "#mu-like; ; events");
      muMChistWeight->GetYaxis()->SetTitleSize(0.05);
  muMChistWeight->Scale(livetimeRatio);

  TH1F *muMChist =
      getAngleDistribution(mcTuple, muCondition + "&& amom < 400", false);
      setPlottingOptions(muMChist, "; ; ");
  muMChist->Scale(livetimeRatio * 0.89);

  ///////////////////////////////
  //       OVER 400 MeV/c      //
  ///////////////////////////////

  TH1F *eDatahistOver400MeV =
      getAngleDistribution(dataTuple, eCondition + "&& amom > 400", false);
      setPlottingOptions(eDatahistOver400MeV, "; cos#theta; ");

  TH1F *eMChistWeightOver400MeV =
      getAngleDistribution(mcTuple, eCondition + "&& amom > 400");
      setPlottingOptions(eMChistWeightOver400MeV, "; cos#theta; events");
      eMChistWeightOver400MeV->GetYaxis()->SetTitleSize(0.04);
  eMChistWeightOver400MeV->Scale(livetimeRatio);

  TH1F *eMChistOver400MeV =
      getAngleDistribution(mcTuple, eCondition + "&& amom > 400", false);
      setPlottingOptions(eMChistOver400MeV, "; cos#theta; ");
  eMChistOver400MeV->Scale(livetimeRatio);

  ////////////////////////////////////////

  TH1F *muDatahistOver400MeV =
      getAngleDistribution(dataTuple, muCondition + "&& amom > 400", false);
      setPlottingOptions(muDatahistOver400MeV, "; cos#theta; ");

  TH1F *muMChistWeightOver400MeV =
      getAngleDistribution(mcTuple, muCondition + "&& amom > 400");
      setPlottingOptions(muMChistWeightOver400MeV, "; cos#theta; events");
      muMChistWeightOver400MeV->GetYaxis()->SetTitleSize(0.04);
  muMChistWeightOver400MeV->Scale(livetimeRatio);

  TH1F *muMChistOver400MeV =
      getAngleDistribution(mcTuple, muCondition + "&& amom > 400", false);
      setPlottingOptions(muMChistOver400MeV, "; cos#theta; ");
  muMChistOver400MeV->Scale(livetimeRatio); 

  ///////////////////////////////

  TCanvas *canvas = new TCanvas("canvas", "", 800, 700);
  canvas->Divide(2, 2, 0, 0);

  ///////////////////////////////

  canvas->cd(1);
  plotOnePad(eDatahist, eMChist, eMChistWeight, "p < 400 MeV/c");

  canvas->cd(2);
  plotOnePad(muDatahist, muMChist, muMChistWeight, "p < 400 MeV/c");

  canvas->cd(1);
  plotOneRatio(eDatahist, eMChist, eMChistWeight);

  canvas->cd(2);
  plotOneRatio(muDatahist, muMChist, muMChistWeight);


  ///////////////////////////////
  //       OVER 400 MeV/c      //
  ///////////////////////////////

  canvas->cd(3);
  plotOnePad(eDatahistOver400MeV, eMChistOver400MeV, eMChistWeightOver400MeV, "p > 400 MeV/c");

  canvas->cd(4);
  plotOnePad(muDatahistOver400MeV, muMChistOver400MeV, muMChistWeightOver400MeV, "p > 400 MeV/c");

  canvas->cd(3);
  plotOneRatio(eDatahistOver400MeV, eMChistOver400MeV, eMChistWeightOver400MeV);

  canvas->cd(4);
  plotOneRatio(muDatahistOver400MeV, muMChistOver400MeV, muMChistWeightOver400MeV);

}
