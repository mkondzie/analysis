#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

void setPlottingOptions(TH1F *hist){

	hist->GetYaxis()->SetNdivisions(505, 4);
	hist->GetXaxis()->SetLabelSize(0.04);
	hist->GetYaxis()->SetLabelSize(0.04);
	hist->GetYaxis()->SetLabelOffset(0.005);
	hist->GetXaxis()->SetLabelOffset(0.005);
	hist->GetXaxis()->CenterTitle(true);
	hist->GetYaxis()->CenterTitle(true);
	hist->GetXaxis()->SetTitleOffset(1.2);
	//gStyle->SetOptStat(0);

}

void loadFile(TString filePath){
  TFile *rootFile = new TFile(filePath, "READ");	
	TTree *oscillationTuple;
	rootFile->GetObject("osc_tuple", oscillationTuple);

  Int_t           nn;
  Int_t           nn_mctruth;
  Int_t           nring;
  Int_t           ipnu;
  Float_t         dirnu[3];
  Float_t         pnu;
  Int_t           mode;
  Int_t           ip;
  Float_t         dprob;
  Float_t         dir[3];
  //Float_t         pos[3]; //no pos leaf in MC files
  Float_t         amom;
  Float_t         path;
  Float_t         wall;
  Int_t           itype;
  Int_t           muedk;
  Float_t         flxg[3];
  Float_t         flxgo[3];
  Float_t         flxh[3];
  Float_t         flxho[3];
  Float_t         weightx;
  Float_t         ring;
  Float_t         probms[5];   //[nring] single ring
  Float_t         prmslg[10][5];   //[nring] multiring
  Float_t         ra;
  Float_t         dec;
  Float_t         time[3];
  Float_t         date[3];
  Float_t         oscweight3f;
  Float_t         cp1p2oscweight3f;
  Float_t         cp3p2oscweight3f;
  Float_t         invoscweight3f;
  Float_t         true_lepmom;
  Float_t         true_lepdir[3];

	oscillationTuple->SetBranchAddress("nn", &nn);
  oscillationTuple->SetBranchAddress("nn_mctruth", &nn_mctruth);
  oscillationTuple->SetBranchAddress("nring", &nring);
  oscillationTuple->SetBranchAddress("ipnu", &ipnu);
  oscillationTuple->SetBranchAddress("dirnu", dirnu);
  oscillationTuple->SetBranchAddress("pnu", &pnu);
  oscillationTuple->SetBranchAddress("mode", &mode);
  oscillationTuple->SetBranchAddress("ip", &ip);
  oscillationTuple->SetBranchAddress("dprob", &dprob);
  oscillationTuple->SetBranchAddress("dir", dir);
  // oscillationTuple->SetBranchAddress("pos", pos);
  oscillationTuple->SetBranchAddress("amom", &amom);
  oscillationTuple->SetBranchAddress("path", &path);
  oscillationTuple->SetBranchAddress("wall", &wall);
  oscillationTuple->SetBranchAddress("itype", &itype);
  oscillationTuple->SetBranchAddress("muedk", &muedk);
  oscillationTuple->SetBranchAddress("flxg", flxg);
  oscillationTuple->SetBranchAddress("flxgo", flxgo);
  oscillationTuple->SetBranchAddress("flxh", flxh);
  oscillationTuple->SetBranchAddress("flxho", flxho);
  oscillationTuple->SetBranchAddress("weightx", &weightx);
  oscillationTuple->SetBranchAddress("ringc", &ring);
  oscillationTuple->SetBranchAddress("probms", &probms);
  oscillationTuple->SetBranchAddress("prmslg", &prmslg);
  oscillationTuple->SetBranchAddress("ra", &ra);
  oscillationTuple->SetBranchAddress("dec", &dec);
  oscillationTuple->SetBranchAddress("time", time);
  oscillationTuple->SetBranchAddress("date", date);
  oscillationTuple->SetBranchAddress("oscweight3f", &oscweight3f);
  oscillationTuple->SetBranchAddress("cp1p2oscweight3f", &cp1p2oscweight3f);
  oscillationTuple->SetBranchAddress("cp3p2oscweight3f", &cp3p2oscweight3f);
  oscillationTuple->SetBranchAddress("invoscweight3f", &invoscweight3f);
  oscillationTuple->SetBranchAddress("true_lepmom", &true_lepmom);
  oscillationTuple->SetBranchAddress("true_lepdir", true_lepdir);
   
}

void makePlots(TString filePath = "fcmc.sk3.19b.mrbdt-2020.0.root"){

	TFile *rootFile = new TFile(filePath, "READ");	
	TTree *oscillationTuple;
	rootFile->GetObject("osc_tuple", oscillationTuple);
	
    TObjArray *branches = oscillationTuple->GetListOfBranches();
    TBranch *branch;

      for(int i = 0; i < branches->GetSize(); i++){
        branch = static_cast<TBranch *>(branches->At(i));
        TString branchName = branch->GetName();
        
				auto leaf = branch->GetLeaf(branchName);
        auto nDim = leaf->GetNdata();
				std::cout << branchName << " " << nDim << std::endl;
		}

}
