#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TTree.h>

void setPlottingOptions(TH1F *hist) {

  hist->GetYaxis()->SetNdivisions(505, 4);
  hist->GetXaxis()->SetNdivisions(505, 4);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetLabelOffset(0.005);
  hist->GetXaxis()->SetLabelOffset(0.005);
  hist->GetXaxis()->CenterTitle(true);
  hist->GetYaxis()->CenterTitle(true);
  hist->GetXaxis()->SetTitleOffset(1.2);
  gStyle->SetOptStat(0);
}

void compareTuple() {
  TFile *MCFile = new TFile("/home/mk/Documents/SK/analysis/expanded_fv/MC/"
                            "fcmc.sk4.19b.mrbdt-2020.all.root",
                            "READ");
  TFile *dataFile = new TFile("/home/mk/Documents/SK/analysis/expanded_fv/data/"
                              "fcdt.sk4.19b.mrbdt-2020.all.root",
                              "READ");
  TTree *oscillationTupleMC;
  TTree *oscillationTupleData;
  MCFile->GetObject("osc_tuple", oscillationTupleMC);
  // oscillationTupleMC->Print();
  dataFile->GetObject("osc_tuple", oscillationTupleData);
  // oscillationTupleData->Print();

  // oscillationTuple->Draw("dirnu[2]/TMath::Sqrt(dirnu[0] * dirnu[0] + dirnu[1]
  // * dirnu[1] + dirnu[2] * dirnu[2])", "nring == 1 && itype == 60 && ip == 2
  // && (mode == 1 || mode == -1)", ""); oscillationTuple->Draw("amom>>h1(100,
  // 180, 1300)", "nring == 1 && itype == 63 && ip == 3 && (mode == 1 || mode ==
  // -1)", ""); oscillationTuple->Draw("itype", " mode == 1", "");
  // oscillationTuple->Draw("dirnu[2]/TMath::Sqrt(dirnu[0] * dirnu[0] + dirnu[1]
  // * dirnu[1] + dirnu[2] * dirnu[2])", "nring == 1 && itype == 63 && ip == 3
  // && (mode == 1 || mode == -1)", "");

  TCanvas *canvas = new TCanvas("canvas1", "", 800, 600);
  canvas->Divide(2, 2);
  TString zenithAngle =
      "-dir[2] / TMath::Sqrt(dir[0] * dir[0] + dir[1] * "
      "dir[1] + dir[2] * dir[2])";
  TString muCut = "amom > 400 && wall > 1 && nring == 1 && ip == 3 &&  "
                  "(mode == 1 || mode == -1)";    //+ " && (ipnu == 14 || ipnu == -14) "
  TString eCut =
      "amom > 400 && wall > 1 && nring == 1 && ip == 2  && "
      "(mode == 1 || mode == -1)";    //+ "&& (ipnu == 12 || ipnu == -12)"

  canvas->cd(1);            
  oscillationTupleMC->Draw(zenithAngle, eCut); // + " * oscweight3f"
  TH1F *histEMC = (TH1F *)gPad->GetPrimitive("htemp");
  histEMC->SetLineColor(kRed);
  histEMC->SetMarkerColor(kRed);
  histEMC->SetTitle("e MC; cos#theta; counts");
  setPlottingOptions(histEMC);

  canvas->cd(2);
  oscillationTupleMC->Draw(zenithAngle, muCut, ""); //+ " * oscweight3f"
  TH1F *histMuMC = (TH1F *)gPad->GetPrimitive("htemp");
  histMuMC->SetLineColor(kBlue);
  histMuMC->SetMarkerColor(kBlue);
  histMuMC->SetTitle("#mu MC; cos#theta; counts");
  setPlottingOptions(histMuMC);

  canvas->cd(3);
  oscillationTupleData->Draw(zenithAngle, "ip == 2", "");
  TH1F *histEdata = (TH1F *)gPad->GetPrimitive("htemp");
  histEdata->SetLineColor(kMagenta);
  histEdata->SetMarkerColor(kMagenta);
  histEdata->SetTitle("e data; cos#theta; counts");
  setPlottingOptions(histEdata);

  canvas->cd(4);
  oscillationTupleData->Draw(zenithAngle, "ip == 3", "");
  TH1F *histMuData = (TH1F *)gPad->GetPrimitive("htemp");
  histMuData->SetLineColor(kGreen);
  histMuData->SetMarkerColor(kGreen);
  histMuData->SetTitle("#mu data; cos#theta; counts");
  setPlottingOptions(histMuData);
}