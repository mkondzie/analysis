//  calibrationSi.c
//  Created by Giulia Colucci on 28/02/24
//  Modified by Maria Kondzielska on 22/04/24

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TVirtualFitter.h"
#include <RtypesCore.h>
#include <TFile.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>

template <typename T> void setPlottingOptions(T *plot) {

  plot->SetMarkerSize(0.5);
  plot->SetMarkerColor(kBlack);
  plot->SetMarkerStyle(21);
  plot->SetLineColor(kAzure + 2);
  plot->SetLineWidth(1);

  plot->GetYaxis()->SetNdivisions(505, 4);
  plot->GetXaxis()->SetNdivisions(505, 4);

  plot->GetXaxis()->CenterTitle(true);
  plot->GetYaxis()->CenterTitle(true);

  gStyle->SetCanvasColor(kWhite);
  gStyle->SetOptStat(0);
}

template <typename T, typename P> T functionFormula(T *x, P *p) {
  return p[0] + p[1] * (*x) + p[2] * TMath::Power(*x, 2) +
         p[3] * TMath::Power(*x, 2.5) + p[4] * TMath::Power(*x, 3);
}

void calibrationSi() {

  std::vector<int> workingDetectors = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15};
  const int nSpectra = workingDetectors.size();
  // spectra
  TH1F *h[nSpectra];

  // alpha variables
  TH1F *halpha[nSpectra];

  TSpectrum *sa[nSpectra];
  Int_t nFoundAlpha[nSpectra];
  Double_t *xpeaksa[nSpectra];

  // pulser variables
  TH1F *hp[nSpectra];

  TSpectrum *s[nSpectra];
  Int_t nFound[nSpectra];
  Double_t *xpeaks[nSpectra];

  TF1 *fChVeryk[nSpectra];
  TGraph *gr[nSpectra];

  Double_t voltage[] = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
                        0.4, 0.45, 0.5, 0.55, 0.6, 0.65};

  Double_t A[nSpectra], B[nSpectra], C[nSpectra], D[nSpectra], E[nSpectra],
      K[nSpectra];

  // auto functionFormula = [](double *x, double *p)
  // {
  //     return p[0] + p[1] * (*x) + p[2] * TMath::Power(*x, 2) + p[3] *
  //     TMath::Power(*x, 2.5) + p[4] * TMath::Power(*x, 3);
  // };

  // double (*functionFormulaPointer)(double *x, double *p) = functionFormula;
  //////////////////

  UShort_t energyChannel;
  UChar_t detectorN;

  TFile *file = new TFile("cal_Si_02.root");
  TTree *cc21 = (TTree *)file->Get("cc21");

  cc21->SetBranchAddress("det", &detectorN);
  cc21->SetBranchAddress("e", &energyChannel);
  /////////////////

  UShort_t lowerBound = 250;
  UShort_t lowerBoundAlpha = 310;
  UShort_t upperBoundAlpha = 385;

  for (const auto &l : workingDetectors) {

    Int_t nEntries = (Int_t)cc21->GetEntries();

    //-----------------Create and fill spectra
    h[l] = new TH1F(Form("h%i", l), Form("Energy spectrum Silicon %i", l), 4092,
                    0, 4092);
    halpha[l] =
        new TH1F(Form("halpha%i", l),
                 Form("Energy alpha spectrum Silicon %i", l), 4092, 0, 4092);
    hp[l] = new TH1F(Form("hp%i", l), Form("Pulser spectrum Silicon %i", l),
                     4092, 0, 4092);

    for (Long_t i = 0; i < nEntries; i++) {

      cc21->GetEntry(i);
      if (detectorN == l) {

        h[l]->Fill(energyChannel);

        if (energyChannel > lowerBound) {

          if (energyChannel > lowerBoundAlpha &&
              energyChannel < upperBoundAlpha) {

            halpha[l]->Fill(energyChannel);

          }

          else if (energyChannel < lowerBoundAlpha ||
                   energyChannel > upperBoundAlpha) {

            hp[l]->Fill(energyChannel);
          }
        }
      }
    }

    //-----------------------------------------

    //---------Find alpha-----------
    // Use TSpectrum to find the peak candidates

    sa[l] = new TSpectrum();

    nFoundAlpha[l] = sa[l]->Search(halpha[l], 2, "nodraw", 0.50);
    std::cout << "----Silicon detector : " << l << std::endl;
    printf("Found %d candidate alpha peaks to fit\n", nFoundAlpha[l]);
    xpeaksa[l] = sa[l]->GetPositionX();

    long long sizesa = nFoundAlpha[l];
    long long indsa[22];
    for (int k = 0; k < nFoundAlpha[l]; k++)
      indsa[k] = 0;
    TMath::Sort(sizesa, xpeaksa[l], indsa, kFALSE);
    for (int j = 0; j < nFoundAlpha[l]; j++) {

      std::cout << xpeaksa[l][indsa[j]] << std::endl;
    }

    //-------------------------

    //---------Find pulser-----------

    // Use TSpectrum to find the peak candidates

    s[l] = new TSpectrum();
    nFound[l] = s[l]->Search(hp[l], 8, "nodraw", 0.10);
    printf("Found %d candidate PULSER peaks to fit\n", nFound[l]);

    xpeaks[l] = s[l]->GetPositionX();

    Double_t nf;
    nf = nFound[l];
    Double_t pos[nFound[l]];

    long long sizes = nFound[l];
    long long inds[22];
    for (int k = 0; k < nFound[l]; k++)
      inds[k] = 0;
    TMath::Sort(sizes, xpeaks[l], inds, kFALSE);

    std::cout << "sorted values " << std::endl;
    for (int j = 0; j < nFound[l]; j++) {

      std::cout << xpeaks[l][inds[j]] << std::endl;
      pos[j] = xpeaks[l][inds[j]];
    }
    std::cout << "---------------" << std::endl;
    //-----------------------

    // fChVeryk[l] = new TF1(Form("fChVeryk%i", l),
    //                       "[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,"
    //                       "2.5)+[4]*TMath::Power(x,3)",
    //                       500, 4000);
    fChVeryk[l] = new TF1(Form("fChVeryk%i", l), functionFormula, 500, 4000, 5);

    ///////////////////////////
    // TString filename =  "kal_bk_det_run_cal_si_02.dat";
    // std::ifstream file(filename.Data());

    // if (!file.is_open()) {
    // std::cerr << "Error opening file: " << filename << std::endl;
    fChVeryk[l]->SetParameter(0, -0.025);
    fChVeryk[l]->SetParameter(1, 0.00018);
    // }

    // TString line;
    // TString format("%d %lf %lf %lf %lf %lf %lf");

    // while (line.ReadLine(file)) {

    //     Int_t lInitial;
    //     Double_t initialA, initialB, initialC, initialD, initialE, initialK;

    //     if (sscanf(line.Data(), format.Data(), &lInitial, &initialA,
    //     &initialB, &initialC, &initialD, &initialE, &initialK) != 7) {
    //         if (file.tellg() != 0) {  // not the beginning of the file
    //             std::cerr << "Error reading line: " << line << std::endl;
    //         }
    //         continue;
    //     }

    //  if(lInitial == l){
    //     fChVeryk[l]->SetParameter(0, initialA);
    //     fChVeryk[l]->SetParameter(1, initialB);
    //     fChVeryk[l]->SetParameter(2, initialC);
    //     fChVeryk[l]->SetParameter(3, initialD);
    //     fChVeryk[l]->SetParameter(4, initialE);}

    // }

    // file.close();

    gr[l] = new TGraph(nf, pos, voltage);
    gr[l]->Fit(fChVeryk[l], "RMEN", "", 500, 3900); // 4000
    // gr[l]->Fit(fChVeryk[l], "0");

    A[l] = fChVeryk[l]->GetParameter(0);
    B[l] = fChVeryk[l]->GetParameter(1);
    C[l] = fChVeryk[l]->GetParameter(2);
    D[l] = fChVeryk[l]->GetParameter(3);
    E[l] = fChVeryk[l]->GetParameter(4);

    double calibrationParameters[5] = {A[l], B[l], C[l], D[l], E[l]};
    Double_t xAlphaPeakChannel;

    // calibrationParameters[0] = A[l];
    // calibrationParameters[1] = B[l];
    // calibrationParameters[2] = C[l];
    // calibrationParameters[3] = D[l];
    // calibrationParameters[4] = E[l];

    xAlphaPeakChannel = xpeaksa[l][0];
    double energyAlphaMeV = 5.4395;

    K[l] = energyAlphaMeV /
           functionFormula(&xAlphaPeakChannel, calibrationParameters);
  }

  TString to_open1 = "kal_bk_det_run_cal_si_02.dat";

  FILE *out1 = fopen(std::string(to_open1).c_str(), "w");

  fprintf(out1, "%-8s %-15s %-15s %-15s %-15s %-15s %-10s\n", "l", "A[l]",
          "B[l]", "C[l]", "D[l]", "E[l]", "K[l]");

  //////
  std::vector<double> NeAuChannel = {
      0.0,    0.0,    0.0,    0.0,    2170.0, 2022.0, 2027.0, 0.0,
      2072.0, 2017.0, 2075.0, 2010.0, 2046.0, 2170.0, 0.0,    2095.0};

  TString to_open2 = "NeAuTargetcal_si_02.dat";

  FILE *out2 = fopen(std::string(to_open2).c_str(), "w");

  fprintf(out2, "%-8s %-15s %-15s\n", "l", "NeAuChannel[l]",
          "NeAuEnergyMeV[l]");
  std::cout << "l\t\t"
            << "NeAuChannel[l]\t"
            << "NeAuEnergyMeV[l]" << std::endl;
  //////
  for (const auto &l : workingDetectors) {

    fprintf(out1, "%-8d %-15g %-15g %-15g %-15g %-15g %-10g\n", l, A[l], B[l],
            C[l], D[l], E[l], K[l]);

    double calibrationParameters[5] = {A[l], B[l], C[l], D[l], E[l]};

    if (l != 7) {
      Double_t NeAuEnergyMeV =
          K[l] * functionFormula(&NeAuChannel.at(l), calibrationParameters);

      std::cout << l << "\t\t" << NeAuChannel.at(l) << "\t\t" << NeAuEnergyMeV
                << std::endl;

      fprintf(out2, "%-8d %-15g %-15g\n", l, NeAuChannel.at(l), NeAuEnergyMeV);
    }
    //////
  }

  fclose(out1);
  fclose(out2);

  /////////////////////////////////////////////////////////////////
  TFile *outFile = new TFile("cc21CalibrationResult.root", "RECREATE");
  TTree *cc21Calibration = cc21->CloneTree(0);
  UShort_t energyMeV;
  cc21Calibration->Branch("eMeV", &energyMeV, "eMeV/s");

  for (auto &l : workingDetectors) {
    double calibrationParameters[5] = {A[l], B[l], C[l], D[l], E[l]};

    for (Long64_t i = 0; i < cc21->GetEntries(); i++) {
      cc21->GetEntry(i);
      if (detectorN == l) {
        energyMeV =
            K[l] * functionFormula(&energyChannel, calibrationParameters);

        cc21Calibration->Fill();
      }
    }
  }
  cc21Calibration->Write();
  ///////////////////////////////////////////////////////////////
  TCanvas *canvas = new TCanvas("canvas ", "canvas ", 800, 600);
  canvas->Divide(5, 3);
  // TCanvas *canvasSpectra =
  // new TCanvas("canvasSpectra ", "canvasSpectra ", 800, 600);
  // canvasSpectra->Divide(5, 3);

  for (const auto &l : workingDetectors) {

    canvas->cd(l);
    gr[l]->GetXaxis()->SetTitle("Channel (V)");
    gr[l]->GetYaxis()->SetTitle("Voltage (V)");
    gr[l]->SetTitle(Form("Voltage(Channel) - Fit Eryk -  Silicon %i", l));
    setPlottingOptions(gr[l]);
    gr[l]->Draw("AP");
    setPlottingOptions(fChVeryk[l]);
    fChVeryk[l]->Draw("same");

    // canvasSpectra->cd(l);
    // h[l]->Draw("Hist");
    // h[l]->GetXaxis()->SetTitle("Channels");
    // h[l]->SetTitle(Form("Spectra -  Silicon %i", l));
    // setPlottingOptions(h[l]);
  }

  /*
     TCanvas *c1spectra = new TCanvas("c1spectra ","c1spectra ",10,10,1000,900);
     c1spectra->Divide(5,3);
    //  c1spectra->cd(1);
    //  h[1]->Draw("Hist");
    //  h[1]->GetXaxis()->SetTitle("Channels");
    //  h[1]->SetTitle("Spectra -  Silicon 4");
    //  h[1]->GetXaxis()->SetRange(1,4092);
    //  c1spectra->cd(2);
    //  h[2]->Draw("Hist");
    //  h[2]->GetXaxis()->SetTitle("Channels");
    //  h[2]->SetTitle("Spectra -  Silicon 5");
    //  h[2]->GetXaxis()->SetRange(1,4092);
    //  c1spectra->cd(3);
    //  h[3]->Draw("Hist");
    //  h[3]->GetXaxis()->SetTitle("Channels");
    //  h[3]->SetTitle("Spectra -  Silicon 6");
    //  h[3]->GetXaxis()->SetRange(1,4092);
     c1spectra->cd(4);
     h[4]->Draw("Hist");
     h[4]->GetXaxis()->SetTitle("Channels");
     h[4]->SetTitle("Spectra -  Silicon 7");
     h[4]->GetXaxis()->SetRange(1,4092);
     c1spectra->cd(6);
     h[5]->Draw("Hist");
     h[5]->GetXaxis()->SetTitle("Channels");
     h[5]->SetTitle("Spectra -  Silicon 8");
     h[5]->GetXaxis()->SetRange(1,4092);
     c1spectra->cd(7);
     h[6]->Draw("Hist");
     h[6]->GetXaxis()->SetTitle("Channels");
     h[6]->SetTitle("Spectra -  Silicon 9");
     h[6]->GetXaxis()->SetRange(1,4092);
     c1spectra->cd(8);
     h[7]->Draw("Hist");
     h[7]->GetXaxis()->SetTitle("Channels");
     h[7]->SetTitle("Spectra -  Silicon 10");
     h[7]->GetXaxis()->SetRange(1,4092);
     c1spectra->cd(9);
     h[8]->Draw("Hist");
     h[8]->GetXaxis()->SetTitle("Channels");
     h[8]->SetTitle("Spectra -  Silicon 11");
     h[8]->GetXaxis()->SetRange(1,4092);
     c1spectra->cd(10);
     h[9]->Draw("Hist");
     h[9]->GetXaxis()->SetTitle("Channels");
     h[9]->SetTitle("Spectra -  Silicon 12");
     h[9]->GetXaxis()->SetRange(1,4092);
     c1spectra->cd(11);
     h[10]->Draw("Hist");
     h[10]->GetXaxis()->SetTitle("Channels");
     h[10]->SetTitle("Spectra -  Silicon 13");
     h[10]->GetXaxis()->SetRange(1,4092);
     c1spectra->cd(12);
     h[11]->Draw("Hist");
     h[11]->GetXaxis()->SetTitle("Channels");
     h[11]->SetTitle("Spectra -  Silicon 14");
     h[11]->GetXaxis()->SetRange(1,4092);
     c1spectra->cd(13);
     h[12]->Draw("Hist");
     h[12]->GetXaxis()->SetTitle("Channels");
     h[12]->SetTitle("Spectra -  Silicon 15");
     h[12]->GetXaxis()->SetRange(1,4092);
     c1spectra->cd(14);
     h[13]->Draw("Hist");
     h[13]->GetXaxis()->SetTitle("Channels");
     h[13]->SetTitle("Spectra -  Silicon 16");
     h[13]->GetXaxis()->SetRange(1,4092);

     TCanvas *c2spectra = new TCanvas("c2spectra ","c2spectra ",10,10,1000,900);
     c2spectra->Divide(5,3);
     c2spectra->cd(1);
     halpha[1]->Draw("Hist");
     halpha[1]->GetXaxis()->SetTitle("Channels");
     halpha[1]->SetTitle("Spectra -  Silicon 4");
     halpha[1]->GetXaxis()->SetRange(310,380);
     c2spectra->cd(2);
     halpha[2]->Draw("Hist");
     halpha[2]->GetXaxis()->SetTitle("Channels");
     halpha[2]->SetTitle("Spectra -  Silicon 5");
     halpha[2]->GetXaxis()->SetRange(310,380);
     c2spectra->cd(3);
     halpha[3]->Draw("Hist");
     halpha[3]->GetXaxis()->SetTitle("Channels");
     halpha[3]->SetTitle("Spectra -  Silicon 6");
     halpha[3]->GetXaxis()->SetRange(310,380);
     c2spectra->cd(4);
     halpha[4]->Draw("Hist");
     halpha[4]->GetXaxis()->SetTitle("Channels");
     halpha[4]->SetTitle("Spectra -  Silicon 7");
     halpha[4]->GetXaxis()->SetRange(310,380);
     c2spectra->cd(6);
     halpha[5]->Draw("Hist");
     halpha[5]->GetXaxis()->SetTitle("Channels");
     halpha[5]->SetTitle("Spectra -  Silicon 8");
     halpha[5]->GetXaxis()->SetRange(310,380);
     c2spectra->cd(7);
     halpha[6]->Draw("Hist");
     halpha[6]->GetXaxis()->SetTitle("Channels");
     halpha[6]->SetTitle("Spectra -  Silicon 9");
     halpha[6]->GetXaxis()->SetRange(310,380);
     c2spectra->cd(8);
     halpha[7]->Draw("Hist");
     halpha[7]->GetXaxis()->SetTitle("Channels");
     halpha[7]->SetTitle("Spectra -  Silicon 10");
     halpha[7]->GetXaxis()->SetRange(310,380);
     c2spectra->cd(9);
     halpha[8]->Draw("Hist");
     halpha[8]->GetXaxis()->SetTitle("Channels");
     halpha[8]->SetTitle("Spectra -  Silicon 11");
     halpha[8]->GetXaxis()->SetRange(310,380);
     c2spectra->cd(10);
     halpha[9]->Draw("Hist");
     halpha[9]->GetXaxis()->SetTitle("Channels");
     halpha[9]->SetTitle("Spectra -  Silicon 12");
     halpha[9]->GetXaxis()->SetRange(310,380);
     c2spectra->cd(11);
     halpha[10]->Draw("Hist");
     halpha[10]->GetXaxis()->SetTitle("Channels");
     halpha[10]->SetTitle("Spectra -  Silicon 13");
     halpha[10]->GetXaxis()->SetRange(310,380);
     c2spectra->cd(12);
     halpha[11]->Draw("Hist");
     halpha[11]->GetXaxis()->SetTitle("Channels");
     halpha[11]->SetTitle("Spectra -  Silicon 14");
     halpha[11]->GetXaxis()->SetRange(310,380);
     c2spectra->cd(13);
     halpha[12]->Draw("Hist");
     halpha[12]->GetXaxis()->SetTitle("Channels");
     halpha[12]->SetTitle("Spectra -  Silicon 15");
     halpha[12]->GetXaxis()->SetRange(310,380);
     c2spectra->cd(14);
     halpha[13]->Draw("Hist");
     halpha[13]->GetXaxis()->SetTitle("Channels");
     halpha[13]->SetTitle("Spectra -  Silicon 16");
     halpha[13]->GetXaxis()->SetRange(310,380);

     TCanvas *c3spectra = new TCanvas("c3spectra ","c3spectra ",10,10,1000,900);
     c3spectra->Divide(5,3);
     c3spectra->cd(1);
     hp[1]->Draw("Hist");
     hp[1]->GetXaxis()->SetTitle("Channels");
     hp[1]->SetTitle("Spectra Pulser -  Silicon 4");
     hp[1]->GetXaxis()->SetRange(1,4092);
     c3spectra->cd(2);
     hp[2]->Draw("Hist");
     hp[2]->GetXaxis()->SetTitle("Channels");
     hp[2]->SetTitle("Spectra Pulser -  Silicon 5");
     hp[2]->GetXaxis()->SetRange(1,4092);
     c3spectra->cd(3);
     hp[3]->Draw("Hist");
     hp[3]->GetXaxis()->SetTitle("Channels");
     hp[3]->SetTitle("Spectra Pulser -  Silicon 6");
     hp[3]->GetXaxis()->SetRange(1,4092);
     c3spectra->cd(4);
     hp[4]->Draw("Hist");
     hp[4]->GetXaxis()->SetTitle("Channels");
     hp[4]->SetTitle("Spectra Pulser -  Silicon 7");
     hp[4]->GetXaxis()->SetRange(1,4092);
     c3spectra->cd(6);
     hp[5]->Draw("Hist");
     hp[5]->GetXaxis()->SetTitle("Channels");
     hp[5]->SetTitle("Spectra Pulser -  Silicon 8");
     hp[5]->GetXaxis()->SetRange(1,4092);
     c3spectra->cd(7);
     hp[6]->Draw("Hist");
     hp[6]->GetXaxis()->SetTitle("Channels");
     hp[6]->SetTitle("Spectra Pulser -  Silicon 9");
     hp[6]->GetXaxis()->SetRange(1,4092);
     c3spectra->cd(8);
     hp[7]->Draw("Hist");
     hp[7]->GetXaxis()->SetTitle("Channels");
     hp[7]->SetTitle("Spectra Pulser -  Silicon 10");
     hp[7]->GetXaxis()->SetRange(1,4092);
     c3spectra->cd(9);
     hp[8]->Draw("Hist");
     hp[8]->GetXaxis()->SetTitle("Channels");
     hp[8]->SetTitle("Spectra Pulser -  Silicon 11");
     hp[8]->GetXaxis()->SetRange(1,4092);
     c3spectra->cd(10);
     hp[9]->Draw("Hist");
     hp[9]->GetXaxis()->SetTitle("Channels");
     hp[9]->SetTitle("Spectra Pulser -  Silicon 12");
     hp[9]->GetXaxis()->SetRange(1,4092);
     c3spectra->cd(11);
     hp[10]->Draw("Hist");
     hp[10]->GetXaxis()->SetTitle("Channels");
     hp[10]->SetTitle("Spectra Pulser -  Silicon 13");
     hp[10]->GetXaxis()->SetRange(1,4092);
     c3spectra->cd(12);
     hp[11]->Draw("Hist");
     hp[11]->GetXaxis()->SetTitle("Channels");
     hp[11]->SetTitle("Spectra Pulser -  Silicon 14");
     hp[11]->GetXaxis()->SetRange(1,4092);
     c3spectra->cd(13);
     hp[12]->Draw("Hist");
     hp[12]->GetXaxis()->SetTitle("Channels");
     hp[12]->SetTitle("Spectra Pulser -  Silicon 15");
     hp[12]->GetXaxis()->SetRange(1,4092);
     c3spectra->cd(14);
     hp[13]->Draw("Hist");
     hp[13]->GetXaxis()->SetTitle("Channels");
     hp[13]->SetTitle("Spectra Pulser -  Silicon 16");
     hp[13]->GetXaxis()->SetRange(1,4092);
  //

     // //---------BWD 1-15 - pulser - ---------
     TCanvas *c1fiteryk = new TCanvas("c1fiteryk","c1fiteryk",10,10,1000,900);
     c1fiteryk->Divide(5,3);
     c1fiteryk->cd(1);
     gr[1]->Draw("AP");
     fChVeryk[1]->Draw("same");
     // gr[1]->GetXaxis()->SetTitle("Channel (V)");
     // gr[1]->GetYaxis()->SetTitle("Voltage (V)");
     // gr[1]->SetTitle("Voltage(Channel) - Fit Eryk -  Silicon 4");
     // gr[1]->GetXaxis()->SetRange(1,4092);
     // gr[1]->GetXaxis()->SetRange(-1000,3092);
     // gr[1]->GetYaxis()->SetRange(-1,1);
     c1fiteryk->cd(2);
     gr[2]->Draw("AP");
     fChVeryk[2]->Draw("same");
     // gr[2]->GetXaxis()->SetTitle("Channel (V)");
     // gr[1]->GetYaxis()->SetTitle("Voltage (V)");
     // gr[2]->SetTitle("Voltage(Channel) - Fit Eryk -  Silicon 5");
     // gr[2]->GetXaxis()->SetRange(1,4092);
     c1fiteryk->cd(3);
     gr[3]->Draw("AP");
     fChVeryk[3]->Draw("same");
     gr[3]->GetXaxis()->SetTitle("Channel (V)");
     gr[3]->GetYaxis()->SetTitle("Voltage (V)");
     gr[3]->SetTitle("Voltage(Channel) - Fit Eryk -  Silicon 6");
     gr[3]->GetXaxis()->SetRange(1,4092);
     c1fiteryk->cd(4);
     gr[4]->Draw("AP");
     fChVeryk[4]->Draw("same");
     gr[4]->GetXaxis()->SetTitle("Channel (V)");
     gr[4]->GetYaxis()->SetTitle("Voltage (V)");
     gr[4]->SetTitle("Voltage(Channel) - Fit Eryk -  Silicon 7");
     gr[4]->GetXaxis()->SetRange(1,4092);
     c1fiteryk->cd(6);
     gr[5]->Draw("AP");
     fChVeryk[5]->Draw("same");
     gr[5]->GetXaxis()->SetTitle("Channel (V)");
     gr[5]->GetYaxis()->SetTitle("Voltage (V)");
     gr[5]->SetTitle("Voltage(Channel) - Fit Eryk -  Silicon 8");
     gr[5]->GetXaxis()->SetRange(1,4092);
     c1fiteryk->cd(7);
     gr[6]->Draw("AP");
     fChVeryk[6]->Draw("same");
     gr[6]->GetXaxis()->SetTitle("Channel (V)");
     gr[6]->GetYaxis()->SetTitle("Voltage (V)");
     gr[6]->SetTitle("Voltage(Channel) - Fit Eryk -  Silicon 9");
     gr[6]->GetXaxis()->SetRange(1,4092);
     c1fiteryk->cd(8);
     gr[7]->Draw("AP");
     fChVeryk[7]->Draw("same");
     gr[7]->GetXaxis()->SetTitle("Channel (V)");
     gr[7]->GetYaxis()->SetTitle("Voltage (V)");
     gr[7]->SetTitle("Voltage(Channel) - Fit Eryk -  Silicon 10");
     gr[7]->GetXaxis()->SetRange(1,4092);
     c1fiteryk->cd(9);
     gr[8]->Draw("AP");
     fChVeryk[8]->Draw("same");
     gr[8]->GetXaxis()->SetTitle("Channel (V)");
     gr[8]->GetYaxis()->SetTitle("Voltage (V)");
     gr[8]->SetTitle("Voltage(Channel) - Fit Eryk -  Silicon 11");
     gr[8]->GetXaxis()->SetRange(1,4092);
     c1fiteryk->cd(10);
     gr[9]->Draw("AP");
     fChVeryk[9]->Draw("same");
     gr[9]->GetXaxis()->SetTitle("Channel (V)");
     gr[9]->GetYaxis()->SetTitle("Voltage (V)");
     gr[9]->SetTitle("Voltage(Channel) - Fit Eryk -  Silicon 12");
     gr[9]->GetXaxis()->SetRange(1,4092);
     c1fiteryk->cd(11);
     gr[10]->Draw("AP");
     fChVeryk[10]->Draw("same");
     gr[10]->GetXaxis()->SetTitle("Channel (V)");
     gr[10]->GetYaxis()->SetTitle("Voltage (V)");
     gr[10]->SetTitle("Voltage(Channel) - Fit Eryk -  Silicon 13");
     gr[10]->GetXaxis()->SetRange(1,4092);
     c1fiteryk->cd(12);
     gr[11]->Draw("AP");
     fChVeryk[11]->Draw("same");
     gr[11]->GetXaxis()->SetTitle("Channel (V)");
     gr[11]->GetYaxis()->SetTitle("Voltage (V)");
     gr[11]->SetTitle("Voltage(Channel) - Fit Eryk -  Silicon 14");
     gr[11]->GetXaxis()->SetRange(1,4092);
     c1fiteryk->cd(13);
     gr[12]->Draw("AP");
     fChVeryk[12]->Draw("same");
     gr[12]->GetXaxis()->SetTitle("Channel (V)");
     gr[12]->GetYaxis()->SetTitle("Voltage (V)");
     gr[12]->SetTitle("Voltage(Channel) - Fit Eryk -  Silicon 15");
     gr[12]->GetXaxis()->SetRange(1,4092);
     c1fiteryk->cd(14);
     gr[13]->Draw("AP");
     fChVeryk[13]->Draw("same");
     gr[13]->GetXaxis()->SetTitle("Channel (V)");
     gr[13]->GetYaxis()->SetTitle("Voltage (V)");
     gr[13]->SetTitle("Voltage(Channel) - Fit Eryk -  Silicon 16");
     gr[13]->GetXaxis()->SetRange(1,4092);


     */

  delete file;
  delete outFile;
}
