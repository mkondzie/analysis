#include <dirent.h>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <cmath>

#include "APFit_variables.h"
#include "fiTQun_variables.h"

#include "APFit_functions.h"
#include "fiTQun_functions.h"

float dirZaxis[3] = {0.,0.,-1.};

void setoutputBranch(TTree *output_tree){
  output_tree->Branch("mc_years", &mc_years);

  output_tree->Branch("mode", &mode);
  output_tree->Branch("wall", &wall);
  output_tree->Branch("nhitac", &nhitac);

  output_tree->Branch("fracmom", &fracmom);
  output_tree->Branch("transmom", &transmom);
  output_tree->Branch("dpose", &dpose);
  output_tree->Branch("mercos", &mercos);
  output_tree->Branch("merprmslg1", &merprmslg1);
  output_tree->Branch("merprmslg2", &merprmslg2);
  output_tree->Branch("merprmslg", &merprmslg);

  //output_tree->Branch("ntag_nn", &ntag_nn);
  //output_tree->Branch("ntag_mctruth_nn", &ntag_mctruth_nn);

  //output_tree->Branch("itype", &itype);
  output_tree->Branch("ip", ip);
  output_tree->Branch("nmue", &nmue);
  output_tree->Branch("evis", &evis);
  output_tree->Branch("nring", &nring);
  
  output_tree->Branch("oscwgt", &oscwgt);
  output_tree->Branch("flxh11", flxh11, "flxh11[3]/F");
  output_tree->Branch("flxh06", flxh06, "flxh06[3]/F");
  output_tree->Branch("live", &live);

  output_tree->Branch("fracmom_2_fiTQun", &fracmom_2_fiTQun);
  output_tree->Branch("likelihood_fiTQun", &likelihood_fiTQun);

  
}

void FillTree(TFile *outfile, TTree **out_trees,  char *dir_path, std::string run_mode, std::string period){

  struct dirent **namelist;
  int r = scandir(dir_path, &namelist, NULL, alphasort);
  //r is # of files or directries (include ./ and ../).
  //If there are errors, -1 is returned.

  mc_years = -1;
  for (int i = 0; i < r; ++i) {
    char *filename = namelist[i]->d_name;
    if( ! strstr(filename, ".root") ) continue;

    mc_years += 1;
    // rvw hax 
    //if ( mc_years > 1 ) continue ; 

    if( run_mode=="test" && !( 100 <= mc_years && mc_years < 200 ) ) continue;

    char allpath[1024];
    sprintf(allpath, "%s%s", dir_path, filename);
    if(mc_years%10==0) std::cout<<"progress: "<<allpath<<std::endl;
    
    TFile *fin = new TFile(allpath,"read");
    TTree *tree = (TTree*)fin->Get("h1");

    Int_t nEntry = tree->GetEntries();
    
    tree->SetBranchAddress("nmue",&nmue);
    tree->SetBranchAddress("evis",&evis);
    tree->SetBranchAddress("nring",&nring);
    tree->SetBranchAddress("amome",amome);
    tree->SetBranchAddress("amomm",amomm);
    
    tree->SetBranchAddress("etime",etime);
    tree->SetBranchAddress("etype",etype);
    tree->SetBranchAddress("epos",epos);
    tree->SetBranchAddress("ehit",ehit);
    tree->SetBranchAddress("egood",egood);
    tree->SetBranchAddress("edir",edir);
    tree->SetBranchAddress("pos",pos);
    
    tree->SetBranchAddress("ip",ip);
    tree->SetBranchAddress("prmslg",prmslg);
    
    tree->SetBranchAddress("dir",dir);
    
    tree->SetBranchAddress("mode",&mode);
    tree->SetBranchAddress("ipnu",ipnu);
    tree->SetBranchAddress("wall",&wall);
    tree->SetBranchAddress("nhitac",&nhitac);
    
    //tree->SetBranchAddress("ntag_nn",&ntag_nn);
    //tree->SetBranchAddress("ntag_mctruth_nn",&ntag_mctruth_nn);
        
    tree->SetBranchAddress("oscwgt", &oscwgt);
    tree->SetBranchAddress("flxh11", flxh11);
    tree->SetBranchAddress("flxh06", flxh06);
    tree->SetBranchAddress("live", &live);
    

    if(period=="sk4_fiTQun"){
      tree->SetBranchAddress("fqnmrfit", &fqnmrfit);
      tree->SetBranchAddress("fqmrifit", fqmrifit);
      tree->SetBranchAddress("fqmrnring", fqmrnring);
      tree->SetBranchAddress("fqmrpcflg", fqmrpcflg);
      tree->SetBranchAddress("fqmrnll", fqmrnll);
      tree->SetBranchAddress("fqmrtotmu", fqmrnll);
      tree->SetBranchAddress("fqmrpid", fqmrpid);
      tree->SetBranchAddress("fqmrmom", fqmrmom);
      tree->SetBranchAddress("fqmrdconv", fqmrdconv);
      tree->SetBranchAddress("fqmreloss", fqmreloss);
      tree->SetBranchAddress("fqmrt0", fqmrt0);
      tree->SetBranchAddress("fqmrpos", fqmrpos);
      tree->SetBranchAddress("fqmrdir", fqmrdir);

      tree->SetBranchAddress("fqnse", &fqnse);
      tree->SetBranchAddress("fq1rpcflg", fq1rpcflg);
      tree->SetBranchAddress("fq1rmom", fq1rmom);
      tree->SetBranchAddress("fq1rt0", fq1rt0);
      tree->SetBranchAddress("fq1rtotmu", fq1rtotmu);
      tree->SetBranchAddress("fq1rnll", fq1rnll);
      tree->SetBranchAddress("fq1rpos", fq1rpos);
      tree->SetBranchAddress("fq1rdir", fq1rdir);
      tree->SetBranchAddress("fq1rdconv", fq1rdconv);
      tree->SetBranchAddress("fq1reloss", fq1reloss);
    }

    int skgen;
    if(period=="sk1") skgen=0;
    if(period=="sk2") skgen=1;
    if(period=="sk3") skgen=2;
    if(period=="sk4" || period=="sk4_fiTQun") skgen=3;
    if(period=="sk5") skgen=4;
        
    for(int i=0;i<nEntry;i++)
      {	
	fin->cd();
	tree->GetEntry(i);

	int mer = 0;
	int mer_id = 0;
	float merp = 0;

	if(period=="sk4_fiTQun"){
	  //--- fiTQun variable calculation ---------------------------
	  const int ifQtoUse = 0;       //0 means best MR fit result
	  int i_fQring[20];     //index of the ring in the order of mom
	  SortfQRing(ifQtoUse,i_fQring);
	  int iPID_fQMER = fqmrpid[ifQtoUse][i_fQring[0]];
	  
	  int ifQnDcy= 0 ;
	  int nMueDcy_fQ = 0;
	  double dMcl_Val = 0.0;
	  GetDcye(ifQtoUse, i_fQring, ifQnDcy,nMueDcy_fQ,dMcl_Val);
	  if (ifQnDcy==0) dMcl_Val = -1.;
	  //double dMcl_norm_Val = TMath::Min(dMcl_Val/fQEvistot,0.999);
	  int nDcy_Val = TMath::Min(10,ifQnDcy);
	  double evisfQ = GetEvistot(0);
	  double fQEvistot = GetEvistot(ifQtoUse);
	  double mom_fQMER = fqmrmom[ifQtoUse][i_fQring[0]];
	  
	  //0 = GAMMA, 1 = ELECTRON, 2 = MUON, 3 = PION, 4 = KAON, 5 = PROTON,  6 = CONE GENERATOR
	  int iMEpipRng = GetfQMER(ifQtoUse,ipip, i_fQring);
	  int iMEmuRng = GetfQMER(ifQtoUse,imu,i_fQring);
	  
	  double EfrcMER_Val = GetRingEvis(ifQtoUse,i_fQring[0])/fQEvistot;
	  double EfrcpipR_Val = TMath::Max(GetRingEvis(ifQtoUse,iMEpipRng),GetRingEvis(ifQtoUse,iMEmuRng))/fQEvistot;
	  int nRing_Val = fqmrnring[ifQtoUse];
	  double pTrans_Val = GetTransMom(ifQtoUse, i_fQring);
	  
	  double wallfq = 0, towallfq = 0;
	  GetfQVtx(ifQtoUse, i_fQring, wallfq, towallfq);

	  //--- substitute fiTQun variables into usual variable (used in APFit) -----
	  evis = fQEvistot;
	  wall = wallfq;
	  nring = nRing_Val;
	  fracmom = EfrcMER_Val;
	  fracmom_2_fiTQun = EfrcpipR_Val;
	  transmom = pTrans_Val;
	  nmue = nDcy_Val;
	  dpose = dMcl_Val;
	  mercos = CalcCos(fqmrdir[ifQtoUse][ i_fQring[0] ], dirZaxis);
	  mer_id = iPID_fQMER;
	  merp = mom_fQMER;
	  
	  likelihood_fiTQun = sigle_ring_likelihood_search(fqnmrfit, fqmrifit, fqmrnll);
	    
	}

	else{ 
	  SearchMER(nring, amome, amomm, prmslg, &mer, &merp);
	  mer_id = ( prmslg[mer][1] < prmslg[mer][2] ? 1 : 2 ); // 1:e 2:mu
	  mercos = CalcCos(dir[mer], dirZaxis);
	}
	
	//-----data reduction and select analysis sample ---------------
	if(period=="sk1"){
	  if(!(evis > 30.0 && nhitac<10))continue;//fc cut
	}
	else {
	  if(!(evis > 30.0 && nhitac<16))continue;//fc cut
	}
	//if(!(200<wall))continue;//fv cut
	//if(!(150<wall))continue;//fv cut (fiTQun?)
		
	//-----SubGeV SingleRing selection---------------
	/*if ( !(evis<1330.0 && nring==1) ) continue;
	  if ( !(  (ip[0]==3 && amomm[0]>200) || (ip[0]==2 &&amome[0]>100)  ) ) continue;
	*/
	//-----MultiGeV SingleRing selection-------------
	//if ( !(evis>1330.0 && nring==1) ) continue;
	 
	//-----MultiGeV MultiRing selection--------------
	if ( !(nring>1) ) continue;
	if ( !(evis>1330.0 || (mer_id==2 && evis>600. && merp>600) ) )continue;
	
	//-----MultiGeV----------------------------------
	//if ( !(evis>1330.0 || (mer_id==2 && evis>600. && merp>600 && nring>1) )  )continue;

	//-----blind analysis------------------------------
	if (  run_mode=="data" && !(mercos>0.4) ) continue;
	//-------------------------------------------------

	
	if (period != "sk4_fiTQun"){
	  fracmom = Fracmom(evis, amome, mer);
	  transmom = Transmom(dir, nring, evis, amome, mer );
         
	  dpose = Dpose(nmue, evis, etime, etype, epos, ehit, pos, egood, amome[mer], skgen );

	  merprmslg1 = prmslg[mer][1];
	  merprmslg2 = prmslg[mer][2];
	  merprmslg = merprmslg1 - merprmslg2;
	  
	}
		
	//------write trees on output files----------------------
	outfile->cd();
	
	if(run_mode == "data"){
	  out_trees[8]->Fill();//DST->Fill();
	  continue;
	}
	
	bool CC_event =  std::abs(mode)<30;
	bool nue_event =  ipnu[0]==12;
	bool nuebar_event =  ipnu[0]==-12;
	bool numu_event =  ipnu[0]==14;
	bool numubar_event =  ipnu[0]==-14;
	bool nutau_nutaubar_event =  std::abs(ipnu[0])==16;

	if (CC_event){
	  if(nue_event){
	    out_trees[0]->Fill();
	    out_trees[2]->Fill();
	  }
	  else if(nuebar_event){
	    out_trees[1]->Fill();
	    out_trees[2]->Fill();
	  }
	  else if(numu_event){
	    out_trees[3]->Fill();
	    out_trees[5]->Fill();
	  }
	  else if(numubar_event){
	    out_trees[4]->Fill();
	    out_trees[5]->Fill();
	  }
	  else if(nutau_nutaubar_event)out_trees[6]->Fill();
	}
	//else out_trees[7]->Fill();
	else if( !CC_event && !nutau_nutaubar_event) out_trees[7]->Fill();
	//-----------------------------------------------------
	
      }
    tree->Delete();
    fin->Close();
  }
  
}



int main(int argc,char *argv[]){
  if(argc != 4){
    std::cout<<"select [run_mode (test, train or data)] [period (sk4, sk4_fiTQun, sk3, sk2 or sk1)] [path to output directory]"<<std::endl;
    return -1;
  }
  std::string run_mode = std::string(argv[1]); //// test train data
  std::string period = std::string(argv[2]); //// sk4 sk4_fiTQun sk3 sk2 sk1
  std::string output_path = std::string(argv[3]); //// path to output directory


  //-----define names of output files and trees----------------
  std::string outfile_name;
  if (run_mode == "data" ) outfile_name = "outfileDST_";
  else if (run_mode == "test" ) outfile_name = "outfileAppMulti_";
  else if (run_mode == "train" ) outfile_name = "outfileMulti_";
  else{
    std::cout<<" select run mode from (data, test, train)"<<std::endl;
    return -1;
  }

  TFile *outfile = new TFile((output_path + "/" + outfile_name + period +".root").c_str(),"recreate");

  const int class_size = 9;
  const char *tree_names[class_size] = {"CCnue","CCnuebar","CCnue_nuebar","CCnumu","CCnumubar","CCnumu_numubar","CCnutau","NC","DST"};
  TTree *out_trees[class_size];
  for(int i=0; i<class_size; i++){
    out_trees[i] = new TTree(tree_names[i], tree_names[i]);
    setoutputBranch(out_trees[i]);
  }
  //-----------------------------------------------------------

  char *input_fcmc;
  char *input_tau;
  char *input_fcmc_add;
  char *input_tau_add;
  char *input_data;
    
  if(period == "sk4"){
    input_fcmc = "/disk01/atmpd5/sk4_dst/apr18/fc_mc/APFit_ntuple/root/";
    input_tau = "/disk01/atmpd5/sk4_dst/apr18/tau_mc/fc_mc/APFit_ntuple/root/";
    input_fcmc_add = "/disk01/atmpd5/sk4_dst/apr18/AdditionalMC/ResonanceMC/fc_mc/root/";
    input_tau_add = "/disk01/atmpd5/sk4_dst/apr18/AdditionalMC/ResonanceMC/tau_mc/fc_mc/root/";
    input_data = "/disk01/atmpd5/sk4_dst/apr18/fc_dst/ntuple/";
  }
  
  else if(period=="sk4_fiTQun"){
    input_fcmc = "/disk01/atmpd5/sk4_dst/apr18/fc_mc/fiTQun_ntuple/root/";
    input_tau = "/disk01/atmpd5/sk4_dst/apr18/tau_mc/fc_mc/fiTQun_ntuple/root/";
    input_data = "/disk01/atmpd5/sk4_dst/apr18/fc_dst/fiTQunv6r0/root/";
    // there are no additional MC used fiTQun so far.
  }

  else if(period=="sk3"){
    input_fcmc = "/disk01/atmpd5/sk3_dst/apr18/fc_mc/ntuple/";
    input_tau = "/disk01/atmpd5/sk3_dst/apr18/tau_mc/fc_mc/root/";
    input_fcmc_add = "/disk01/atmpd5/sk3_dst/apr18/AdditionalMC/ResonanceMC/fc_mc/root/";
    input_tau_add = "/disk01/atmpd5/sk3_dst/apr18/AdditionalMC/ResonanceMC/tau_mc/fc_mc/root/";
    input_data = "/disk01/atmpd5/sk3_dst/apr18/fc_dst/ntuple/";
  }
  
  
  else if(period=="sk2"){
    input_fcmc = "/disk01/atmpd5/sk2_dst/apr18/fc_mc/ntuple/";
    input_tau = "/disk01/atmpd5/sk2_dst/apr18/tau_mc/fc_mc/root/";
    input_fcmc_add = "/disk01/atmpd5/sk2_dst/apr18/AdditionalMC/ResonanceMC/fc_mc/root/";
    input_tau_add = "/disk01/atmpd5/sk2_dst/apr18/AdditionalMC/ResonanceMC/tau_mc/fc_mc/root/";
    input_data = "/disk01/atmpd5/sk2_dst/apr18/fc_dst/ntuple/";
  }

  else if(period=="sk1"){
    input_fcmc = "/disk01/atmpd5/sk1_dst/apr18/fc_mc/ntuple/";
    input_tau = "/disk01/atmpd5/sk1_dst/apr18/tau_mc/fc_mc/root/";
    input_fcmc_add = "/disk01/atmpd5/sk1_dst/apr18/AdditionalMC/ResonanceMC/fc_mc/root/";
    input_tau_add = "/disk01/atmpd5/sk1_dst/apr18/AdditionalMC/ResonanceMC/tau_mc/fc_mc/root/";
    input_data = "/disk01/atmpd5/sk1_dst/apr18/fc_dst/ntuple/";
  }
  else{
    std::cout<<" select run period (sk1, sk2, sk3, sk4, sk4_fiTQun)"<<std::endl;
    return -1;
  }
  
  //-----Fill event--------------------------------------
  if(run_mode == "train"){
    FillTree(outfile, out_trees, input_fcmc_add, run_mode, period);
    FillTree(outfile, out_trees, input_tau_add, run_mode, period);
  }
  if(run_mode == "test"){
    FillTree(outfile, out_trees, input_fcmc, run_mode, period);
    FillTree(outfile, out_trees, input_tau, run_mode, period);
  }
  else if(run_mode == "data"){
    FillTree(outfile, out_trees, input_data, run_mode, period);
  }
  //-----------------------------------------------------
  
  
  outfile->cd();
  if(run_mode == "data") out_trees[class_size-1]->Write();
  else for(int i=0; i<class_size-1; i++) out_trees[i]->Write();
  
  std::cout<<"written"<<std::endl;
  outfile->Close();
  return 0;
}





