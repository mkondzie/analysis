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

std::map<Int_t, TString> mode_labels = {
    // Charged-current
    {1, "CCQE"},
    {-1, "CCQE"},
    {2, "2p2h"},
    {-2, "2p2h"},
    {11, "CCRes1#it{#pi^{+}}"},
    {13, "CCRes1#it{#pi^{+}}"},
    {12, "CCRes1#it{#pi^{0}}"},
    {-12, "CCRes1#it{#pi^{0}}"},
    {-11, "CCRes1#it{#pi^{-}}"},
    {-13, "CCRes1#it{#pi^{-}}"},
    {15, "CCDif1#it{#pi^{+}}"},
    {-15, "CCDif1#it{#pi^{-}}"},
    {16, "CCCoh1#it{#pi^{+}}"},
    {-16, "CCCoh1#it{#pi^{-}}"},
    {17, "CCRes1#it{#gamma}"},
    {-17, "CCRes1#it{#gamma}"},
    {21, "CCN#it{#pi}"},
    {-21, "CCN#it{#pi}"},
    {22, "CCRes1#it{#eta^{0}}"},
    {-22, "CCRes1#it{#eta^{0}}"},
    {23, "CCRes1#it{K^{0}}"},
    {-23, "CCRes1#it{K^{+}}"},
    {26, "CCDIS"},
    {-26, "CCDIS"},

    // Neutral-current
    {31, "NCRes1#pi^{0}"},
    {-31, "NCRes1#pi^{0}"},
    {32, "NCRes1#pi^{0}"},
    {-32, "NCRes1#pi^{0}"},
    {33, "NCRes1#pi^{-}"},
    {-33, "NCRes1#pi^{-}"},
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

std::map<TString, TString> variable_labels = {
    {"merprmslg", "Most energetic ring PID"}};

std::map<TString, TString> sample_labels = {
    {"CCnue", "#nu_{e} CC"},
    {"CCnuebar", "#bar{#nu}_{e} CC"},
    {"CCnue_nuebar", "#nu_{e} & #bar{#nu}_{e} CC"},
    {"CCnumu", "#nu_{#mu} CC"},
    {"CCnumubar", "#bar{#nu}_{#mu} CC"},
    {"CCnumu_numubar", "#nu_{#mu} & #bar{#nu}_{#mu} CC"},
    {"CCnutau", "#nu_{#tau} CC"},
    {"NC", "NC"}};

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
    stack->GetXaxis()->SetTitleOffset(0.92);
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

    return n_bins;
}

void plot_abs_scale_mod(TString sample = "CCnue")
{
    // Open the ROOT file and get the tree
    TFile *file = new TFile("outfileAppMulti_sk4_ip.root", "READ");
    TTree *tree;
    file->GetObject(sample, tree);

    // Declare branches
    Float_t merprmslg;
    Int_t mode;
    tree->SetBranchAddress("merprmslg", &merprmslg);
    tree->SetBranchAddress("mode", &mode);

    // Determine the range and binning for the histograms
    auto [x_min, x_max] = find_min_max_X(tree, "merprmslg");
    if (sample.BeginsWith("CCnue"))
    {
        x_min = -5500;
        x_max = 500;
    }
    else if (sample.BeginsWith("CCnumu"))
    {
        x_min = -3000;
        x_max = 2350;
    }
    else if (sample.BeginsWith("CCnutau"))
    {
        x_min = -13000;
        x_max = 750;
    }
    else if (sample.BeginsWith("NC"))
    {
        x_min = -5000;
        x_max = 600;
    }
    Int_t n_bins = find_n_bins(tree, "merprmslg");

    // Step 1: Count occurrences of modes with the same absolute value
    std::map<Int_t, Float_t> mode_frequencies;
    tree->Draw("mode");
    TH1F *mode_hist = (TH1F *)gPad->GetPrimitive("htemp");

    for (int i = 1; i <= mode_hist->GetNbinsX(); ++i)
    {
        Int_t mode_value = mode_hist->GetBinCenter(i); // Get mode value
        if (mode_value <= 0)
        {
            mode_value -= 1; // Adjust mode value for negative modes
        }
        Float_t frequency = mode_hist->GetBinContent(i); // Get mode frequency
        if (frequency > 0)
        {
            Int_t abs_mode_value = std::abs(mode_value); // Use absolute value for merging
            mode_frequencies[abs_mode_value] += frequency;
        }
    }

    delete mode_hist;

    // Step 2: Sort merged modes by frequency in descending order
    std::vector<std::pair<Int_t, Float_t>> sorted_modes(mode_frequencies.begin(), mode_frequencies.end());
    std::sort(sorted_modes.begin(), sorted_modes.end(),
              [](const auto &a, const auto &b)
              { return a.second > b.second; });

    // Step 3: Select the top 5 merged modes with labels
    std::vector<Int_t> top_modes;
    for (const auto &[mode_value, frequency] : sorted_modes)
    {
        if (mode_labels.find(mode_value) != mode_labels.end())
        {
            top_modes.push_back(mode_value);
        }
        if (top_modes.size() >= 5)
            break; // Stop after 5 modes
    }

    // Step 4: Create a THStack for merged modes
    THStack *stack = new THStack("stack", "Stacked Interaction Modes");

    // Step 5: Create a legend
    TLegend *legend = new TLegend(0.15, 0.6, 0.3, 0.85);
    legend->SetTextSize(0.03);
    legend->SetBorderSize(0);

    // Step 6: Add histograms for the top 5 merged modes
    std::map<Int_t, TH1F *> mode_histograms;
    for (const auto &mode_value : top_modes)
    {
        TString mode_label = mode_labels[mode_value];

        // Create a histogram for this merged mode
        TString hist_name = Form("hist_mode_%d", mode_value);
        TH1F *hist = new TH1F(hist_name, "", n_bins, x_min, x_max);

        // Fill the histogram with events passing the condition
        TString draw_condition = Form("abs(mode) == %d", mode_value);
        tree->Draw(Form("merprmslg>>%s", hist_name.Data()), draw_condition);
        delete gROOT->FindObject("c1");

        // Skip empty histograms
        if (hist->GetEntries() == 0)
        {
            delete hist;
            continue;
        }

        // hist->Scale(1.0 / hist->Integral());

        // hist->Scale(hist->Integral());

        // Style the histogram
        int color_index;
        if (sample.BeginsWith("NC"))
        {
            color_index = mode_value * 2 - 20; // Subtract 20 for NC samples
        }
        else
        {
            color_index = mode_value * 2 + 20; // Add 20 for other samples
        }

        set_plotting_options(hist, "merprmslg", color_index);
        hist->SetFillColor(color_index);
        hist->SetFillStyle(1001); // Solid fill

        // Add the histogram to the stack
        stack->Add(hist);
        mode_histograms[mode_value] = hist;

        // Add an entry to the legend with frequency percentage
        Float_t percentage = (mode_frequencies[mode_value] / tree->GetEntries()) * 100;
        legend->AddEntry(hist, Form("%s (%.1f%%)", mode_label.Data(), percentage), "f");
    }

    // Step 7: Add a histogram for "Others"
    TH1F *others_hist = new TH1F("others_hist", "", n_bins, x_min, x_max);
    for (Int_t i = 0; i < tree->GetEntries(); ++i)
    {
        tree->GetEntry(i);

        // Check if the absolute mode is NOT in the top_modes list
        if (std::find(top_modes.begin(), top_modes.end(), std::abs(mode)) == top_modes.end())
        {
            others_hist->Fill(merprmslg);
        }
    }

    // Style the "Others" histogram
    int others_color = 1; // Black
    // others_hist->Scale(1.0 / others_hist->Integral());
    // others_hist->Scale(others_hist->Integral());

    set_plotting_options(others_hist, "merprmslg", others_color);
    others_hist->SetFillColor(others_color);
    others_hist->SetFillStyle(3004); // Dotted fill

    // Add "Others" to the stack and legend
    stack->Add(others_hist);

    Float_t others_percentage = (others_hist->GetEntries() / tree->GetEntries()) * 100;
    legend->AddEntry(others_hist, Form("Other (%.1f%%)", others_percentage), "f");

    // Step 8: Draw the stack on a canvas
    TCanvas *canvas = new TCanvas("canvas", "PID", 800, 600);
    stack->Draw("HIST");
    stack->GetXaxis()->SetTitle(variable_labels["merprmslg"]);
    stack->GetYaxis()->SetTitle("events");


     if (sample.BeginsWith("NC"))
        {
                stack->SetTitle("                   " + sample_labels[sample]);
        }
        else
        {
               stack->SetTitle("    " + sample_labels[sample]);
        }


    set_plotting_options(stack);

    // Draw the legend
    legend->Draw();
    canvas->SaveAs(Form("%s_PID_abs_scale_modes.pdf", sample.Data()));

    // Clean up
    delete file;
}
