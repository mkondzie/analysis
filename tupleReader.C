#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>

void setPlottingOptions(TH1F *hist){

	hist->GetYaxis()->SetNdivisions(505, 4);
	hist->GetXaxis()->SetLabelSize(0.04);
	hist->GetYaxis()->SetLabelSize(0.04);
	hist->GetYaxis()->SetLabelOffset(0.005);
	hist->GetXaxis()->SetLabelOffset(0.005);
	hist->GetXaxis()->CenterTitle(true);
	hist->GetYaxis()->CenterTitle(true);
	hist->GetXaxis()->SetTitleOffset(1.2);

}

void plotMomentum(TTree *tuple){
  float momentum;
  tuple->SetBranchAddress("amom", &momentum);

      TH1F *momentumHist = new TH1F("", "", 100, 100, 1300);
      Long_t nEntries = tuple->GetEntries();
      for(Long_t i = 0; i < nEntries; i++){
        tuple->GetEntry(i);
        momentumHist->Fill(momentum);
      }

  momentumHist->GetXaxis()->SetNdivisions(505, 4);
	momentumHist->GetXaxis()->SetTitle("p (MeV/c)");
	momentumHist->GetYaxis()->SetTitle("counts");
	momentumHist->SetTitle("SK4 MC");
  setPlottingOptions(momentumHist);
	momentumHist->SetLineColor(kAzure + 2);
	momentumHist->SetFillColor(kAzure + 2);
	momentumHist->SetFillStyle(3002);

	gStyle->SetCanvasColor(kWhite);
	gStyle->SetOptStat(0);

  TCanvas *canvas = new TCanvas("canvas", "Momentum - MC", 15,15,900,900);
    
    canvas->SetLeftMargin(0.12);
    canvas->SetRightMargin(0.05);
    canvas->SetTopMargin(0.08);
    canvas->SetBottomMargin(0.1);
    canvas->Draw();
    canvas->cd();
	  momentumHist->Draw("");
}

void tupleReader(TString fileName = "/home/mk/Documents/SK/analysis/expanded_fv/MC/fcmc.sk4.19b.mrbdt-2020.9.root"){
      TFile *rootFile = new TFile(fileName, "READ");
      TTree *oscillationTuple;
      rootFile->GetObject("osc_tuple", oscillationTuple);
      oscillationTuple->Print();
      //plotMomentum(oscillationTuple);
      //oscillationTuple->Draw("dirnu[2]/TMath::Sqrt(dirnu[0] * dirnu[0] + dirnu[1] * dirnu[1] + dirnu[2] * dirnu[2])", "nring == 1 && itype == 60 && ip == 2 && (mode == 1 || mode == -1)", "");
      //oscillationTuple->Draw("amom>>h1(100, 180, 1300)", "nring == 1 && itype == 63 && ip == 3 && (mode == 1 || mode == -1)", "");
      //oscillationTuple->Draw("itype", " mode == 1", "");
      // oscillationTuple->Draw("dirnu[2]/TMath::Sqrt(dirnu[0] * dirnu[0] + dirnu[1] * dirnu[1] + dirnu[2] * dirnu[2])", "nring == 1 && itype == 63 && ip == 3 && (mode == 1 || mode == -1)", "");

    TCanvas *canvas = new TCanvas("canvas1", "", 800, 600);
    canvas->Divide(1,2);

    canvas->cd(1);                                                                                                              //
    oscillationTuple->Draw("dirnu[2]/TMath::Sqrt(dirnu[0] * dirnu[0] + dirnu[1] * dirnu[1] + dirnu[2] * dirnu[2])", "nring == 1 && ip == 3 && itype == 63 && (mode == 1 || mode == -1)");
    TH1F *histMu = (TH1F*)gPad->GetPrimitive("htemp"); // Get the histogram object
    histMu->SetLineColor(kBlue); // Set line color
    histMu->SetMarkerColor(kBlue); // Set marker color
    histMu->SetTitle("mu; cos#theta; counts");
    setPlottingOptions(histMu);

    canvas->cd(2);                                                                                                            //
    oscillationTuple->Draw("dirnu[2]/TMath::Sqrt(dirnu[0]*dirnu[0] + dirnu[1]*dirnu[1] + dirnu[2]*dirnu[2])", "nring == 1 && ip == 2 && itype == 60 && (mode == 1 || mode == -1)");
    TH1F *histE = (TH1F*)gPad->GetPrimitive("htemp"); // Get the second histogram object
    histE->SetLineColor(kRed); // Set line color
    histE->SetMarkerColor(kRed); // Set marker color
    histE->SetTitle("e; cos#theta; counts");
    setPlottingOptions(histE);

}