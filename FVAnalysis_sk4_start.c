#define FVAnalysis_sk4_start_cxx
#include "FVAnalysis_sk4_start.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
using namespace std;

//#include "Event.h"

   // 0 data files
   // >0 MC files
   // 1 all MC
   // 2 vue CC
   //3  bar nue CC
   // 4 numu CC
   //5 bar numu CC
   // 6 NC all
// may be added
   //7 nutau CC
   // 8 bar nu tau CC
   //

FVAnalysis_sk4_start *E = new FVAnalysis_sk4_start();
//E->SetVersion(4); // for SK4, 3- sk3 etc
//E->SetDatatype(2); //0- fcData, 2 fcMC

//const int analysed_files=1; //0 data, 1 mc change in the FVAnalysis_sk4.h

//numbers taken from ~/osc/Card/sk.16b/standard/q13-fixed.full.nh.card

const double mc_livetime=182625.0; //eq 500years sk4
const double data_livetime=3244.4;// sk4 fcdata


void FVAnalysis_sk4_start::Loop()
{
    //   In a ROOT session, you can do:
    //      root> .L FVAnalysis_sk4_start.C
    //      root> FVAnalysis_sk4_start t
    //      root> t.GetEntry(12); // Fill t data members with entry number 12
    //      root> t.Show();       // Show values of entry 12
    //      root> t.Show(16);     // Read and show values of entry 16
    //      root> t.Loop();       // Loop on all entries
    //
    
    //     This is the loop skeleton where:
    //    jentry is the global entry number in the chain
    //    ientry is the entry number in the current Tree
    //  Note that the argument to GetEntry must be:
    //    jentry for TChain::GetEntry
    //    ientry for TTree::GetEntry and TBranch::GetEntry
    //
    //       To read only selected branches, Insert statements like:
    // METHOD1:
    //    fChain->SetBranchStatus("*",0);  // disable all branches
    //    fChain->SetBranchStatus("branchname",1);  // activate branchname
    // METHOD2: replace line
    //    fChain->GetEntry(jentry);       //read all branches
    //by  b_branchname->GetEntry(ientry); //read only this branch
    
    gROOT->Time();
    gROOT->Reset();
    
    //some extras
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1, 0);
    gStyle->SetOptStat("e");	//prints only number of events, underflows and overflows
    
    //gStyle->SetOptStat(0);
    gStyle->SetOptFit(111);
    gStyle->SetTitleSize(0.07, "X");
    gStyle->SetTitleSize(0.07, "Y");
    gStyle->SetLabelSize(0.05, "X");
    gStyle->SetLabelSize(0.05, "Y");
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadBottomMargin(0.2);
    gStyle->SetTitleFontSize(0.08);
    //gStyle->SetTitleX(0.4);
    gStyle->SetTitleFontSize(0.08);
    gStyle->SetMarkerSize(1.0);
    //gStyle->SetTitle0ffset(1.5, "Y");
    gStyle->SetPaintTextFormat("2.2f");
    gStyle->SetFitFormat("2.2f");
    gStyle->SetStatFormat("2.2f");
    
    
    Char_t run[300];
    sprintf(run,"%s_fc_data_mc_all",sk_per);
    
    Char_t input[300];
    sprintf(input,"%s_fc_data_mc_all", sk_per);
    //OUTPUT ROOT
    
    
    
    // OUTPUT PDF
    Char_t outfile_data[300];
    sprintf(outfile_data, "/home/mk/Documents/SK/analysis/outputData/APFIT_momentum_data_%s.pdf",run);
    
    Char_t outfile_mc[300];
    sprintf(outfile_mc, "/home/mk/Documents/SK/analysis/outputMC/APFIT_momentum_mc_%s.pdf",run);
    
    
    
    
    Int_t nbins_mom=200;
    Float_t min_p=0;
    Float_t max_p=2000;
    
    
    TH1F *mom_data = new TH1F("mom_data","data", nbins_mom, min_p,max_p);
    TH1F *mom_mc = new TH1F("mom_mc","MC", nbins_mom, min_p,max_p);
    
    mom_mc->GetXaxis()->SetNdivisions(505, 4);
	//mom_mc->GetXaxis()->SetTitle("z (mm)");
	//mom_mc->GetYaxis()->SetTitle("counts");
	mom_mc->GetYaxis()->SetNdivisions(409);
	mom_mc->GetXaxis()->SetLabelSize(0.04);
	mom_mc->GetYaxis()->SetLabelSize(0.04);
	mom_mc->GetYaxis()->SetLabelOffset(0.005);
	mom_mc->GetXaxis()->SetLabelOffset(0.005);
	mom_mc->GetXaxis()->CenterTitle(true);
	mom_mc->GetYaxis()->CenterTitle(true);
	mom_mc->GetXaxis()->SetTitleOffset(1.2);
        
    if (fChain == 0) return;
    
    
    Long64_t nentries = fChain->GetEntries();
    cout<<" number of nentries "<< nentries<<endl;
    
    
    double ratio= data_livetime/mc_livetime;
    cout<<"ratio data/mc  = "<<ratio<<endl;
    
    // TCut eSingleRingQuasiElasticCut = "(((osc_tuple.nring == 1) && (osc_tuple.ip == 2)) && ((osc_tuple.mode == 1) || (osc_tuple.mode == -1)))";
    // TTreeFormula selectedTreeFormula("eSingleRingQuasiElasticSelected",eSingleRingQuasiElasticCut,fChain);


    for (Long64_t i=0; i<nentries;i++)
    {
        
        E->GetEntry(i);
        
        if(analysed_files==0) //data files
        {
            
                mom_data->Fill(E->amom);
            
        }
           
            if(analysed_files==1)//MC files
            {
                
                    mom_mc->Fill(E->amom);
            }

    }
    
 TCanvas *c1 = new TCanvas("c1", "Momentum - MC", 15,15,900,900);
    
    c1->SetLeftMargin(0.05);
    c1->SetRightMargin(0.05);
    c1->SetTopMargin(0.05);
    c1->SetBottomMargin(0.05);
    c1->SetFillColor(kWhite);
    c1->Draw();
    
    c1->cd();
    mom_mc->Draw();
    c1->Print(outfile_mc);
}
