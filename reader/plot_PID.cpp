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
#include <THStack.h>
#include <map>
#include <vector>

// Mapping of mode enumeration to more legible interaction labels
std::map<Int_t, TString> mode_labels = {
    // Charged-current
    {1, "CCQE"},
    {-1, "CCQE"},
    {2, "2p2h"},
    {-2, "2p2h"},
    {11, "CCRes1#pi^{+}"},
    {13, "CCRes1#pi^{+}"},
    {12, "CCRes1#pi^{0}"},
    {-12, "CCRes1#pi^{0}"},
    {-11, "CCRes1#pi^{−}"},
    {-13, "CCRes1#pi^{−}"},
    {15, "CCDif1#pi^{+}"},
    {-15, "CCDif1#pi^{−}"},
    {16, "CCCoh1#pi^{+}"},
    {-16, "CCCoh1#pi^{−}"},
    {17, "CCRes1#gamma"},
    {-17, "CCRes1#gamma"},
    {21, "CCN#pi"},
    {-21, "CCN#pi"},
    {22, "CCRes1#eta^{0}"},
    {-22, "CCRes1#eta^{0}"},
    {23, "CCRes1K^{0}"},
    {-23, "CCRes1K^{+}"},
    {26, "CCDIS"},
    {-26, "CCDIS"},

    // Neutral-current
    {31, "NCRes1#pi^{0}"},
    {-31, "NCRes1#pi^{0}"},
    {32, "NCRes1#pi^{0}"},
    {-32, "NCRes1#pi^{0}"},
    {33, "NCRes1#pi^{−}"},
    {-33, "NCRes1#pi^{−}"},
    {34, "NCRes1#pi^{+}"},
    {-34, "NCRes1#pi^{+}"},
    {35, "NCDif1#pi^{0}"},
    {-35, "NCDif1#pi^{0}"},
    {36, "NCCoh1#pi^{+}"},
    {-36, "NCCoh1#pi^{+}"},
    {38, "NCRes1#gamma"},
    {-38, "NCRes1#gamma"},
    {39, "NCRes1#gamma"},
    {-39, "NCRes1#gamma"},
    {41, "NCN#pi"},
    {-41, "NCN#pi"},
    {42, "NCRes1#eta^{0}"},
    {-42, "NCRes1#eta^{0}"},
    {43, "NCRes1#eta^{0}"},
    {-43, "NCRes1#eta^{0}"},
    {44, "NCRes1K^{0}"},
    {-44, "NCRes1K^{0}"},
    {45, "NCRes1K^{+}"},
    {-45, "NCRes1K^{+}"},
    {46, "NCDIS"},
    {-46, "NCDIS"},
    {51, "NCEL (1p1h)"},
    {-51, "NCEL (1p1h)"},
    {52, "NCEL (1p1h)"},
    {-52, "NCEL (1p1h)"}};

std::map<TString, TString> sample_labels = {
    {"CCnue", "#nu_{e} CC"},
    {"CCnuebar", "#bar{#nu}_{e} CC"},
    {"CCnue_nuebar", "#nu_{e} & #bar{#nu}_{e} CC"},
    {"CCnumu", "#nu_{#mu} CC"},
    {"CCnumubar", "#bar{#nu}_{#mu} CC"},
    {"CCnumu_numubar", "#nu_{#mu} & #bar{#nu}_{#mu} CC"},
    {"CCnutau", "#nu_{#tau} CC"},
    {"NC", "NC"}};

std::map<TString, TString> variable_labels = {
    {"merprmslg", "Most energetic ring PID"}};

// std::vector<TString> samples = {"NC", "CCnutau", "CCnumu_numubar", "CCnuebar", "CCnue"};

// Function to set plotting options for each histogram
template <typename T>
void set_plotting_options(T *hist, TString variable_name = "", int color = 0)
{
  TString title = ";" + variable_labels[variable_name] + "; events";
  hist->SetTitle(title);
  hist->SetLineColor(color);
  hist->SetLineWidth(1);
  gStyle->SetOptStat(0);
}

// Function to set plotting options for each stack
template <>
void set_plotting_options<THStack>(THStack *stack, TString, int)
{
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
  stack->GetXaxis()->SetLabelSize(0.035);
  stack->GetYaxis()->SetLabelSize(0.035);
}

std::pair<float, float> find_min_max_X(TTree *tree, TString variable_name)
{
  tree->Draw(variable_name);
  TH1F *htemp = (TH1F *)gPad->GetPrimitive("htemp");
  float max = htemp->GetXaxis()->GetBinUpEdge(htemp->GetNbinsX());
  float min = htemp->GetXaxis()->GetBinLowEdge(1);
  delete htemp;
  gROOT->cd();
  delete gROOT->FindObject("c1");
  return std::make_pair(min, max);
}

Int_t find_n_bins(TTree *tree, TString variable_name)
{
  tree->Draw(variable_name);
  TH1F *htemp = (TH1F *)gPad->GetPrimitive("htemp");
  Int_t n_bins = htemp->GetNbinsX();
  delete htemp;
  gROOT->cd();
  delete gROOT->FindObject("c1");
  return n_bins;
}

void plot_PID(TString sample = "CCnue")
{
  TFile *file = new TFile("outfileAppMulti_sk4_ip.root", "READ");
  TTree *tree;
  file->GetObject(sample, tree);
  Float_t merprmslg;
  Int_t mode;
  tree->SetBranchAddress("merprmslg", &merprmslg);
  tree->SetBranchAddress("mode", &mode);

  tree->Draw("mode");
  TH1F *htemp = (TH1F *)gPad->GetPrimitive("htemp");

  std::vector<std::pair<Int_t, Float_t>> mode_frequencies;
  for (int i = 1; i <= htemp->GetNbinsX(); ++i)
{
    Int_t mode_value = htemp->GetBinCenter(i);  
    if (mode_value <= 0)
    {
        mode_value -= 1; // Adjust mode value for negative modes
    }
    Float_t frequency = htemp->GetBinContent(i); 
    {if(frequency > 0)
        mode_frequencies.push_back(std::make_pair(mode_value, frequency)); 
    }
}

  // Sort the modes by frequency in descending order
  std::sort(mode_frequencies.begin(), mode_frequencies.end(),
            [](const auto &a, const auto &b)
            { return a.second > b.second; });

  // Get the top 5 modes
  std::vector<Int_t> top_modes;
  for (size_t i = 0; i < 5 && i < mode_frequencies.size(); ++i)
  {
    top_modes.push_back(mode_frequencies[i].first);
  }

  std::vector<std::pair<Int_t, Float_t>> mode_percentages;
  for (size_t i = 0; i < mode_frequencies.size(); ++i)
  {
    mode_percentages.push_back(std::make_pair(mode_frequencies.at(i).first, mode_frequencies.at(i).second / htemp->GetEntries()));
   std::cout << "Mode: " << mode_frequencies[i].first
              << ", Frequency: " << mode_frequencies[i].second << ", Percentage: " << mode_percentages.at(i).second * 100.0 << std::endl;
  }

  // Double_t dominant_mode = tree->GetMaximum("mode");
  // tree->Draw("merprmslg",Form("mode==%f", dominant_mode),"");
}