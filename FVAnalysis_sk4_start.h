 //////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Mar 15 22:55:21 2020 by ROOT version 6.12/04
// from TChain osc_tuple/
//////////////////////////////////////////////////////////

#ifndef FVAnalysis_sk4_start_h
#define FVAnalysis_sk4_start_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
//#include "global.h"
#include <iostream>
using namespace std;

const int nHistos = 10;
const int analysed_files=1;

//for FVanalysis_reaction_types_sk4.c option 0 for data and 1 for MC and NHistos = 10
//for FVanalysis_sk4.c add option 2 tau neutrinos and Nhsitos=9

string sk_period("sk4");

Char_t sk_per[200]={"sk4"};

// Header file for the classes stored in the TTree if any.

class FVAnalysis_sk4_start {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
    
    
    
    
    
     
    
    
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
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
    Float_t        pos[3];
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
    Float_t        ringc;
    Float_t        probms[5];   //[nring] single ring
    Float_t        prmslg[10][5];   //[nring] multiring
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

   // List of branches
   TBranch        *b_nn;   //!
   TBranch        *b_nn_mctruth;   //!
   TBranch        *b_nring;   //!
   TBranch        *b_ipnu;   //!
   TBranch        *b_dirnu;   //!
   TBranch        *b_pnu;   //!
   TBranch        *b_mode;   //!
   TBranch        *b_ip;   //!
   TBranch        *b_dprob;   //!
   TBranch        *b_dir;   //!
   TBranch        *b_pos;   //!
   TBranch        *b_amom;   //!
   TBranch        *b_path;   //!
   TBranch        *b_wall;   //!
   TBranch        *b_itype;   //!
   TBranch        *b_muedk;   //!
   TBranch        *b_flxg;   //!
   TBranch        *b_flxgo;   //!
   TBranch        *b_flxh;   //!
   TBranch        *b_flxho;   //!
   TBranch        *b_weightx;   //!
    TBranch        *b_ringc;   //!
    TBranch        *b_probms;
    TBranch        *b_prmslg;
   TBranch        *b_ra;   //!
   TBranch        *b_dec;   //!
   TBranch        *b_time;   //!
   TBranch        *b_date;   //!
   TBranch        *b_oscweight3f;   //!
   TBranch        *b_cp1p2oscweight3f;   //!
   TBranch        *b_cp3p2oscweight3f;   //!
   TBranch        *b_invoscweight3f;   //!
   TBranch        *b_true_lepmom;   //!
   TBranch        *b_true_lepdir;   //!

   FVAnalysis_sk4_start(TTree *tree=0);
   // FVAnalysis_sk4 *E = new FVAnalysis_sk4();
    
   virtual ~FVAnalysis_sk4_start();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
    
    enum DataTypes{
        fcData,
        pcData,
        fcMC,
        pcMC,
        upmuMC,
        
        SK1fcData,
        SK1pcData,
        SK1fcMC,
        SK1pcMC,
        SK1upmuMC,
        
        SK2fcData,
        SK2pcData,
        SK2fcMC,
        SK2pcMC,
        SK2upmuMC,
        
        SK3fcData,
        SK3pcData,
        SK3fcMC,
        SK3pcMC,
        SK3upmuMC,
        
        SK4fcData,
        SK4pcData,
        SK4fcMC,
        SK4pcMC,
        SK4upmuMC,
    };
    
    enum EventTypes{
        SubGeVE0dcy,
        SubGeVE1dcy,
        SubGeVPi01R,
        SubGeVMu0dcy,
        SubGeVMu1dcy,
        SubGeVMu2dcy,
        SubGeVPi0MR,
        MulGeVE,
        MulGeVMu,
        MultiRingE,
        MultiRingMu,
        PCstop,
        PCthrough,
        UPstopMu,
        NonShoweringMu,
        ShoweringMu
    };

    
    int  _version;  // SK detector version
    int _datatype;
    
    void SetVersion( int x)  { _version = x; }
    int  GetVersion( void)   { return _version;}
    
    int  GetDatatype(void)   { return _datatype;}
    void SetDatatype(int x)  { _datatype = x;}
    
    double weight();
    
    int _datalivetime;
    
    void  set_datalivetime ( double x ) { _datalivetime = x; cout << _datalivetime << endl;}
    double get_datalivetime( void ) { return _datalivetime;}
    
   
  
    
   
    
    
   };

#endif

#ifdef FVAnalysis_sk4_start_cxx
FVAnalysis_sk4_start::FVAnalysis_sk4_start(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0 && analysed_files==0)
   {

    #ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/mk/Documents/SK/analysis/expanded_fv/data/fcdt.sk4.19b.mrbdt-2020.9.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/home/mk/Documents/SK/analysis/expanded_fv/data/fcdt.sk4.19b.mrbdt-2020.9.root");
          cout<<"reading one SK4 data file"<<f->GetName()<<endl;
      }
      f->GetObject("osc_tuple",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("osc_tuple","");
      chain->Add("/home/mk/Documents/SK/analysis/expanded_fv/data/fcdt.sk4.19b.mrbdt-2020.*.root/osc_tuple");

       cout<<"reading chain SK4 data file "<<chain->GetName()<<endl;
      tree = chain;
#endif // SINGLE_TREE

   }
    
    if (tree == 0 && analysed_files==1)
    {
        
#ifdef SINGLE_TREE
        // The following code should be used if you want this class to access
        // a single tree instead of a chain
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/mk/Documents/SK/analysis/expanded_fv/MC/fcmc.sk4.19b.mrbdt-2020.9.root");
        if (!f || !f->IsOpen()) {
            f = new TFile("/home/mk/Documents/SK/analysis/expanded_fv/MC/fcmc.sk4.19b.mrbdt-2020.9.root");
            cout<<"reading one MC file"<<endl;
        }
        f->GetObject("osc_tuple",tree);
        
#else // SINGLE_TREE
        
        // The following code should be used if you want this class to access a chain
        // of trees.
        TChain * chain = new TChain("osc_tuple","");
        chain->Add("/home/mk/Documents/SK/analysis/expanded_fv/MC/fcmc.sk4.19b.mrbdt-2020.*.root/osc_tuple");
        
        
        cout<<"reading chain file for SK4 MC"<<endl;
        tree = chain;
#endif // SINGLE_TREE
        
    }
    if (tree == 0 && analysed_files==2)
        {
            
    #ifdef SINGLE_TREE
            // The following code should be used if you want this class to access
            // a single tree instead of a chain
            TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/Users/magdapz/Documents/analiza_sk/ntuple_files/1_2_m_FV/taumc.sk4.19b.skmeeting-202002.9.root");
            if (!f || !f->IsOpen()) {
                f = new TFile("/Users/magdapz/Documents/analiza_sk/ntuple_files/1_2_m_FV/taumc.sk4.19b.mrbdt-2020.9.root");
                cout<<"reading one MC SK4 tau file"<<endl;
            }
            f->GetObject("osc_tuple",tree);
            
    #else // SINGLE_TREE



    
            
            // The following code should be used if you want this class to access a chain
            // of trees.
            TChain * chain = new TChain("osc_tuple","");
            chain->Add("~magdapz/Documents/analiza_sk/ntuple_files/1_2_m_FV/taumc.sk4.19b.mrbdt-2020.*.root/osc_tuple");
            
            
            cout<<"reading chain file for SK4 tau MC"<<endl;
            tree = chain;
    #endif // SINGLE_TREE
            
        }
 
   Init(tree);
}

FVAnalysis_sk4_start::~FVAnalysis_sk4_start()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t FVAnalysis_sk4_start::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t FVAnalysis_sk4_start::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void FVAnalysis_sk4_start::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nn", &nn, &b_nn);
   fChain->SetBranchAddress("nn_mctruth", &nn_mctruth, &b_nn_mctruth);
   fChain->SetBranchAddress("nring", &nring, &b_nring);
   fChain->SetBranchAddress("ipnu", &ipnu, &b_ipnu);
   fChain->SetBranchAddress("dirnu", dirnu, &b_dirnu);
   fChain->SetBranchAddress("pnu", &pnu, &b_pnu);
   fChain->SetBranchAddress("mode", &mode, &b_mode);
   fChain->SetBranchAddress("ip", &ip, &b_ip);
   fChain->SetBranchAddress("dprob", &dprob, &b_dprob);
   fChain->SetBranchAddress("dir", dir, &b_dir);
    //fChain->SetBranchAddress("pos", pos, &b_pos);
   fChain->SetBranchAddress("amom", &amom, &b_amom);
   fChain->SetBranchAddress("path", &path, &b_path);
   fChain->SetBranchAddress("wall", &wall, &b_wall);
   fChain->SetBranchAddress("itype", &itype, &b_itype);
   fChain->SetBranchAddress("muedk", &muedk, &b_muedk);
   fChain->SetBranchAddress("flxg", flxg, &b_flxg);
   fChain->SetBranchAddress("flxgo", flxgo, &b_flxgo);
   fChain->SetBranchAddress("flxh", flxh, &b_flxh);
   fChain->SetBranchAddress("flxho", flxho, &b_flxho);
   fChain->SetBranchAddress("weightx", &weightx, &b_weightx);
    fChain->SetBranchAddress("ringc", &ringc, &b_ringc);
    fChain->SetBranchAddress("probms", &probms, &b_probms);
    fChain->SetBranchAddress("prmslg", &prmslg, &b_prmslg);
   fChain->SetBranchAddress("ra", &ra, &b_ra);
   fChain->SetBranchAddress("dec", &dec, &b_dec);
   fChain->SetBranchAddress("time", time, &b_time);
   fChain->SetBranchAddress("date", date, &b_date);
   fChain->SetBranchAddress("oscweight3f", &oscweight3f, &b_oscweight3f);
   fChain->SetBranchAddress("cp1p2oscweight3f", &cp1p2oscweight3f, &b_cp1p2oscweight3f);
   fChain->SetBranchAddress("cp3p2oscweight3f", &cp3p2oscweight3f, &b_cp3p2oscweight3f);
   fChain->SetBranchAddress("invoscweight3f", &invoscweight3f, &b_invoscweight3f);
   fChain->SetBranchAddress("true_lepmom", &true_lepmom, &b_true_lepmom);
   fChain->SetBranchAddress("true_lepdir", true_lepdir, &b_true_lepdir);
   Notify();
}

Bool_t FVAnalysis_sk4_start::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

    cout<<"New file is opened"<<endl;
   return kTRUE;
}

void FVAnalysis_sk4_start::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t FVAnalysis_sk4_start::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef FVAnalysis_sk4_cxx
