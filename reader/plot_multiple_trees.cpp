#include "TError.h"
#include "TROOT.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
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
template <typename T>
void set_plotting_options(T *hist, TString variable_name = "", int color = 0) {
  TString title = ";" + variable_labels[variable_name] + "; events";
  hist->SetTitle(title);
  hist->SetLineColor(color);
  hist->SetLineWidth(1);
  gStyle->SetOptStat(0);
}

// Function to set plotting options for each stack
template <> void set_plotting_options<THStack>(THStack *stack, TString, int) {
  stack->GetXaxis()->CenterTitle(true);
  stack->GetYaxis()->CenterTitle(true);
  stack->GetYaxis()->SetNdivisions(505, 4);
  stack->GetXaxis()->SetNdivisions(505, 4);
  stack->GetXaxis()->CenterTitle(true);
  stack->GetYaxis()->CenterTitle(true);
  stack->GetYaxis()->SetTitleOffset(1.7);
  stack->GetXaxis()->SetTitleSize(0.05);
  stack->GetYaxis()->SetTitleSize(0.042);
  stack->GetYaxis()->SetTitleOffset(1.3);
  stack->GetXaxis()->SetTitleOffset(0.95);
  stack->GetXaxis()->SetLabelSize(0.04);
  stack->GetYaxis()->SetLabelSize(0.035);
}

std::pair<float, float> find_min_max_X(TTree *tree, TString variable_name) {
  tree->Draw(variable_name);
  TH1F *htemp = (TH1F *)gPad->GetPrimitive("htemp");
  float max = htemp->GetXaxis()->GetBinUpEdge(htemp->GetNbinsX());
  float min = htemp->GetXaxis()->GetBinLowEdge(1);
  delete htemp;
  gROOT->cd();
  delete gROOT->FindObject("c1");
  return std::make_pair(min, max);
}

Int_t find_n_bins(TTree *tree, TString variable_name) {
  tree->Draw(variable_name);
  TH1F *htemp = (TH1F *)gPad->GetPrimitive("htemp");
  Int_t n_bins = htemp->GetNbinsX();
  delete htemp;
  gROOT->cd();
  delete gROOT->FindObject("c1");
  return n_bins;
}

void fill_histogram(TTree *tree, TTreeFormula *formula, TH1F *hist,
                    TString variable_name, void *var_ptr, bool is_float) {
  Long64_t n_entries = tree->GetEntries();
  for (Long64_t i = 0; i < n_entries; ++i) {
    tree->GetEntry(i);
    if (formula->EvalInstance()) {
      if (is_float) {
        hist->Fill(*(Float_t *)var_ptr);
      } else {
        hist->Fill(*(Int_t *)var_ptr);
      }
    }
  }
}

void plot_multiple_trees() {
  gErrorIgnoreLevel = kError;
  TFile *file = new TFile("outfileAppMulti_sk4.root", "READ");

  std::vector<TString> samples = {"CCnue", "CCnuebar"};
  std::vector<TString> variables = {"nring", "merprmslg", "nmue",   "dpose",
                                    "evis",  "transmom",  "fracmom"};
  std::vector<bool> is_float = {false, true, false, true, true, true, true};
  std::vector<TH1F *> histograms;

  TCanvas *canvas =
      new TCanvas("canvas", "SKIVMC", 10, 10, 1400, 800);
  canvas->Divide(4, 2, 0.02, 0.06);

  // FC selectrion cut && evis > 30
  TString condition = "wall > 200 && evis > 1330 && nhitac < 16"; //for multi-GeV
  std::vector<int> colors = {kRed + 1, kRed - 8};

  for (size_t var_idx = 0; var_idx < variables.size(); ++var_idx) {
    TString variable_name = variables[var_idx];
    canvas->cd(var_idx + 1);
    gPad->SetTicks(1, 1);
    if (variable_name.BeginsWith("dpose")) {
      gPad->SetLogy();
    } else {
      gPad->SetLogy(0);
    }

    THStack *stack =
        new THStack(variable_name + "_stack",
                    ";" + variable_labels[variable_name] + "; events");

    for (size_t sample_idx = 0; sample_idx < samples.size(); ++sample_idx) {

      TString sample = samples[sample_idx];
      TTree *tree;
      file->GetObject(sample, tree);

      Float_t merprmslg, dpose, evis, transmom, fracmom;
      Int_t nring, nmue;

      tree->SetBranchAddress("nring", &nring);
      tree->SetBranchAddress("merprmslg", &merprmslg);
      tree->SetBranchAddress("nmue", &nmue);
      tree->SetBranchAddress("dpose", &dpose);
      tree->SetBranchAddress("evis", &evis);
      tree->SetBranchAddress("transmom", &transmom);
      tree->SetBranchAddress("fracmom", &fracmom);
      std::vector<void *> var_ptrs = {&nring, &merprmslg, &nmue,   &dpose,
                                      &evis,  &transmom,  &fracmom};

      auto min_max_pair = find_min_max_X(tree, variable_name);
      if (variable_name.BeginsWith("evis"))
        min_max_pair.second = 12e3;
      if (variable_name.BeginsWith("merprmslg"))
        min_max_pair = {-3.5e3, 0.5e3};
      if (variable_name.BeginsWith("dpose"))
        min_max_pair.second = 2.2e3;

      auto n_bins = find_n_bins(tree, variable_name);
      if(variable_name.BeginsWith("fracmom")){
        n_bins = 18;
      }
      else if(variable_name.BeginsWith("transmom")){
        n_bins = 22;
      }
      else if(variable_name.BeginsWith("evis")){
        n_bins = 90;
      }
      TH1F *hist = new TH1F(variable_name + "_hist_" + sample, "",
                            n_bins,
                            min_max_pair.first, min_max_pair.second);

      TTreeFormula *formula = new TTreeFormula("FC", condition, tree);
      fill_histogram(tree, formula, hist, variable_name, var_ptrs[var_idx],
                     is_float[var_idx]);
      hist->Scale(1 / hist->Integral());
      histograms.push_back(hist);

      stack->Add(hist);
      set_plotting_options(hist, variables[var_idx], colors[sample_idx]);

      delete formula;
    }

    stack->Draw("NOSTACK HIST");
    set_plotting_options(stack);
    gPad->Update();
  }

  canvas->cd(8);
  TString legend_title = "SK IV MC multi-ring events";
  TLegend *legend = new TLegend(0.01, 0.6, 0.48, 0.99);
  legend->SetHeader(legend_title);
  legend->SetTextSize(0.07);
  legend->SetBorderSize(0);
  for (Int_t sample_idx = 0; sample_idx < samples.size(); ++sample_idx) {
    legend->AddEntry(histograms.at(sample_idx),
                     sample_labels[samples.at(sample_idx)], "l");
  }
  legend->Draw();
  canvas->Update();
  delete file;
}
