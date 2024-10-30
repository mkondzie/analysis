#include "TError.h"
#include "TROOT.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTree.h>
#include <TTreeFormula.h>
#include <map>
#include <vector>

// Mapping of variable and sample names to more legible labels
std::map<TString, TString> variable_labels = {
    {"nring", "Number of rings"},
    {"merprmslg", "Most energetic ring PID"},
    {"nmue", "Number of decay electrons"},
    {"dpose", "L_{decay e} (cm)"},
    {"evis", "E_{Vis} (MeV)"},
    {"transmom", "T_{mom}"},
    {"fracmom", "F_{mom}"}};

std::map<TString, TString> sample_labels = {
    {"CCnue", "#nu_{e} CC"},
    {"CCnuebar", "#bar{#nu}_{e} CC"},
    {"CCnue_nuebar", "#nu_{e} & #bar{#nu}_{e} CC"},
    {"CCnumu", "#nu_{#mu} CC"},
    {"CCnumubar", "#bar{#nu}_{#mu} CC"},
    {"CCnumu_numubar", "#nu_{#mu} & #bar{#nu}_{#mu} CC"},
    {"CCnutau", "#nu_{#tau} CC"},
    {"NC", "NC"}};

// Function to set plotting options for each histogram
void set_plotting_options(TH1F *hist, TString variable_name) {
  // Use mapped label for axis titles
  TString title = ";" + variable_labels[variable_name] + "; events";
  hist->GetYaxis()->SetNdivisions(505, 4);
  hist->GetXaxis()->SetNdivisions(505, 4);
  hist->GetXaxis()->CenterTitle(true);
  hist->GetYaxis()->CenterTitle(true);
  hist->GetYaxis()->SetTitleOffset(1.7);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetLabelSize(0.035);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.042);
  hist->GetYaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleOffset(0.95);
  hist->SetTitle(title);
  hist->SetMarkerStyle(4);
  hist->SetMarkerSize(0.8);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetOptStat(0);
}

// Function to find min and max of a given branch variable
std::pair<float, float> find_min_max(TTree *tree, TString variable_name) {
  tree->Draw(variable_name);
  TH1F *htemp = (TH1F *)gPad->GetPrimitive("htemp");
  float max = htemp->GetXaxis()->GetBinUpEdge(htemp->GetNbinsX());
  float min = htemp->GetXaxis()->GetBinLowEdge(1);
  delete htemp;
  gROOT->cd();
  delete gROOT->FindObject("c1");
  return std::make_pair(min, max);
}

// Function to find n bins of a given branch variable hist
Int_t find_n_bins(TTree *tree, TString variable_name) {
  tree->Draw(variable_name);
  TH1F *htemp = (TH1F *)gPad->GetPrimitive("htemp");
  Int_t n_bins = htemp->GetNbinsX();
  delete htemp;
  gROOT->cd();
  delete gROOT->FindObject("c1");
  return n_bins;
}

// Function to fill histogram for each variable
void fill_histogram(TTree *tree, TTreeFormula *formula, TH1F *hist,
                    TString variable_name, void *var_ptr, bool is_float) {
  Long64_t n_entries = tree->GetEntries();

  // Loop over entries
  for (Long64_t i = 0; i < n_entries; ++i) {
    tree->GetEntry(i);
    if (formula->EvalInstance()) {
      // Depending on whether the variable is float or int, cast it accordingly
      if (is_float) {
        hist->Fill(*(Float_t *)var_ptr); // For Float_t variables
      } else {
        hist->Fill(*(Int_t *)var_ptr); // For Int_t variables
      }
    }
  }
}

void plot_one_tree(TString sample = "CCnuebar") {

  gErrorIgnoreLevel = kError;

  TFile *file = new TFile("outfileAppMulti_sk4.root", "READ");
  TTree *tree;
  file->GetObject(sample, tree);
  // Branch variables
  Float_t merprmslg, dpose, evis, transmom, fracmom;
  Int_t nring, nmue;

  tree->SetBranchAddress("nring", &nring);
  tree->SetBranchAddress("merprmslg", &merprmslg);
  tree->SetBranchAddress("nmue", &nmue);
  tree->SetBranchAddress("dpose", &dpose);
  tree->SetBranchAddress("evis", &evis);
  tree->SetBranchAddress("transmom", &transmom);
  tree->SetBranchAddress("fracmom", &fracmom);

  // Variables to plot
  std::vector<TString> variables = {"nring", "merprmslg", "nmue",   "dpose",
                                    "evis",  "transmom",  "fracmom"};
  std::vector<void *> var_ptrs = {&nring, &merprmslg, &nmue,   &dpose,
                                  &evis,  &transmom,  &fracmom};
  std::vector<bool> is_float = {
      false, true, false, true,
      true,  true, true}; // Indicate if the variable is float or int

  // Create histograms for each variable
  std::vector<TH1F *> histograms;

  // FC events cut condition
  TString condition = "wall > 200 && evis > 30 && nhitac < 16";
  TTreeFormula *formula = new TTreeFormula("FC", condition, tree);

  // Initialize histograms
  for (size_t i = 0; i < variables.size(); ++i) {
    TString variable_name = variables[i];
    auto min_max_pair = find_min_max(tree, variable_name);

    // Create histogram and store in vector
    if (variables.at(i).BeginsWith("evis")) {
      min_max_pair.second = 12e3;
    } else if (variables.at(i).BeginsWith("merprmslg")) {

      min_max_pair.first = -3.5e3;
      min_max_pair.second = 0.5e3;
    } else if (variables.at(i).BeginsWith("dpose")) {

      min_max_pair.second = 2.1e3;
    }

    TH1F *hist = new TH1F(variable_name + "_hist",
                          "" + sample + ";" + variable_labels[variable_name] +
                              "; events",
                          find_n_bins(tree, variable_name), min_max_pair.first,
                          min_max_pair.second);
    set_plotting_options(hist, variable_name);
    histograms.push_back(hist);
  }

  // Fill histograms
  for (size_t i = 0; i < variables.size(); ++i) {
    fill_histogram(tree, formula, histograms[i], variables[i], var_ptrs[i],
                   is_float[i]);
  }

  // Create canvas with 8 pads
  TCanvas *canvas = new TCanvas("canvas", sample, 10, 10, 1400, 800);
  canvas->Divide(4, 2, 0.015, 0.045);

  // Draw histograms on canvas
  for (size_t i = 0; i < histograms.size(); ++i) {
    canvas->cd(i + 1);
    gPad->SetTicks(1, 1);
    histograms[i]->Draw();
    if (variables.at(i).BeginsWith("dpose")) {
      gPad->SetLogy();
    } else {
      gPad->SetLogy(0);
    }
  }
  // Add last pad for legend
  canvas->cd(8);
  TLatex latex;
  TString legend_title = "SK IV MC multi-ring events";
  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.SetTextColor(kBlue + 4);
  latex.DrawText(0.1, 0.9, legend_title);
  TLegend *legend = new TLegend(0.1, 0.6, 0.48, 0.7);
  // legend->SetHeader (legend_title);
  legend->AddEntry(histograms.at(0), sample_labels[sample], "l");
  legend->SetBorderSize(0);
  legend->Draw();
  // Update the canvas to show all histograms
  canvas->Update();
  delete file;
  delete formula;
}
