//--- my fitqun function ---//
float sigle_ring_likelihood_search(float fqnmrfit, int *fqmrifit, float *fqmrnll){
  //float fqnmrfit; // num of fit
  //int fqmrifit[50]; //fit ID
  //float fqmrnll[50]; //-lnL

  int ID[2];
  for(int i = 0; i<fqnmrfit; i++){
    //select 8bit and 1ring
    if (TMath::Abs(fqmrifit[i]) == 10000003) ID[0] = i;// pi-like
    else if (TMath::Abs(fqmrifit[i]) == 10000001) ID[1] = i;// e-like
  }
  return fqmrnll[ID[0]] - fqmrnll[ID[1]];

}
//--- fitqun function ---//
int ie = 1; int ipip=3; int imu=2;
#define nPID 7
const double m0[nPID]={0.,0.5110,105.6584,139.5702,493.677,938.2720,0.};//Particle rest mass

double pCherenkovThr(double m0tmp) {
  const double nphase=1.334;
  return m0tmp/sqrt(nphase*nphase-1.);
}

double GetElossMax(int iPID, double mom) {
  return TMath::Max(0.,sqrt(mom*mom+m0[iPID]*m0[iPID])-sqrt(pow(pCherenkovThr(m0[iPID]),2.)+m0[iPID]*m0[iPID]));
}

double GetRingEvis(int ifit, int iring) {

  if (iring<0) {
    return 0.01;
  }

  int iPID = fqmrpid[ifit][iring];
  if (iPID==ie) {//electron ring
    return fqmrmom[ifit][iring];
  }
  else if (iPID==imu) {//muon ring
    return GetElossMax(iPID,fqmrmom[ifit][iring]);
  }
  else {//upstream track
    return fqmreloss[ifit][iring];
  }
}


double GetEvistot(int ifit) {
  double Evistmp=0.;
  for (int iring=0; iring<fqmrnring[ifit]; iring++) {
    Evistmp+=GetRingEvis(ifit,iring);
  }
  return Evistmp;
}

int SortfQRing(int ifit, int *idxSrt) {
  /*
    if(nev(0) == 1089678)
    {
        cout << "HERE" << endl;
        }
  */
  double RingEvis[20];
  for (int iring=0; iring<fqmrnring[ifit]; iring++) {// Sort ring energy in decending order
    
    double momtmp=0.;
    int iPID = fqmrpid[ifit][iring];
    if (iPID==ie) {//electron ring
      momtmp = fqmrmom[ifit][iring];
    }
    else if (iPID==imu) {//muon ring
      momtmp = fqmrmom[ifit][iring];
    }
    else {//upstream track - use muon-equiValent momentum!
      //        double tmpEmu = fQMR_Eloss[ifit][iring]+sqrt(pow(pCherenkovThr(m0[imu]),2.)+m0[imu]*m0[imu]);
      //        momtmp = sqrt(tmpEmu*tmpEmu-m0[imu]*m0[imu]);
      momtmp = fqmreloss[ifit][iring];// use visible energy for pi+ rings!
    }
    RingEvis[iring] = momtmp;

    idxSrt[iring]=iring;
    for (int jring=iring-1; jring>=0; jring--) {
      if (RingEvis[idxSrt[jring]]<RingEvis[idxSrt[jring+1]]) {
	int idxBuf=idxSrt[jring];
	idxSrt[jring]=idxSrt[jring+1];
	idxSrt[jring+1]=idxBuf;
      }
      else break;
    }
  }

  return fqmrpid[ifit][idxSrt[0]];// MER PID
}

void GetDcye(int ifQtoUse, int *i_fQring, int& ifQnDcy, int & nMueDcy_fQ, double &dMcl_Val)
{
  double fQdMclMR[10];//distance to the farthest decay electron from the interaction vertex
  for (int iring=0; iring<fqmrnring[ifQtoUse]; iring++) {
    fQdMclMR[iring]=-1.;
    for (int ise=1; ise<fqnse; ise++) {

      double tDcy_fQ = fq1rt0[ise][ie]-fqmrt0[ifQtoUse][iring];
      if (tDcy_fQ<600.) continue;

      double dmcltmp=0.;
      for (int i=0; i<3; i++) {
	dmcltmp += (fq1rpos[ise][ie][i]-fqmrpos[ifQtoUse][iring][i])*(fq1rpos[ise][ie][i]-fqmrpos[ifQtoUse][iring][i]);
      }
      dmcltmp = sqrt(dmcltmp);
      if (fQdMclMR[iring]<dmcltmp) fQdMclMR[iring]=dmcltmp;
    }
  }
  
  
  for (int ise=1; ise<fqnse; ise++) {
    //      cout << Form(" iDcy=%d, t=%f, p=%f",ise,fq1rt0[ise][ie],fq1rmom[ise][ie]) << endl;
    double tDcy_fQ = fq1rt0[ise][ie]-fqmrt0[ifQtoUse][i_fQring[0]];
    if (tDcy_fQ>10000.) {
      //        cout << "  Skipping!" << endl;
      continue;
    }
    ifQnDcy++;
    if (tDcy_fQ>100.) {
      nMueDcy_fQ++;
    }
  }

  dMcl_Val = fQdMclMR[i_fQring[0]];

}

int GetfQMER(int ifit, int iPID, int *i_fQring) {
  
  int iMER=-1;
  
  //      int i_fQring[20];
  //      SortfQRing(ifit,i_fQring,byMom);
  
  for (int i=0; i<fqmrnring[ifit]; i++) {
    if (iPID>=0) {
      if (iPID!=fqmrpid[ifit][i_fQring[i]]) continue;
    }
    
    iMER = i_fQring[i];
    break;
  }

  return iMER;
}


double GetDwall(double *tmppos){
  double dwall;
  double dwalls=1690.-sqrt(tmppos[0]*tmppos[0]+tmppos[1]*tmppos[1]);
  double dwallt=1810.-fabs(tmppos[2]);

  if (dwalls<dwallt) {
    dwall=dwalls;
  }
  else {
    dwall=dwallt;
  }

  return dwall;
}

double GetToWall(Double_t *postmp, Double_t *dirtmp){

  int flgToWall=0;
  double ToWall = 0;
  if(postmp[0]*postmp[0]+postmp[1]*postmp[1] > 1690*1690) return -sqrt(postmp[0]*postmp[0]+postmp[1]*postmp[1] - 1690*1690);      //Outside
  if(fabs(dirtmp[0]*dirtmp[0] + dirtmp[1]*dirtmp[1]) < 1e-6) ToWall = dirtmp[2]>0? 10000: -10000;  //Vertical
  else
    ToWall=(-(postmp[0]*dirtmp[0]+postmp[1]*dirtmp[1])+sqrt((postmp[0]*dirtmp[0]+postmp[1]*dirtmp[1])
							    *(postmp[0]*dirtmp[0]+postmp[1]*dirtmp[1])+(dirtmp[0]*dirtmp[0]+dirtmp[1]*dirtmp[1])
							    *(1690.*1690.-postmp[0]*postmp[0]-postmp[1]*postmp[1])))/(dirtmp[0]*dirtmp[0]+dirtmp[1]*dirtmp[1]);

  if ((postmp[2]+ToWall*dirtmp[2])>1810.) {//penetrates top wall!
    ToWall=(1810.-postmp[2])/dirtmp[2];
    flgToWall=1;
  }
  else if ((postmp[2]+ToWall*dirtmp[2])<-1810.) {//penetrates bottom wall!
    ToWall=(-1810.-postmp[2])/dirtmp[2];
    flgToWall=-1;
  }
  
  return ToWall;
}

void GetfQVtx(int ifQtoUse, int* i_fQring, double &dwall, double &towall)
{
  double fqpos[10][3];
  double fqdir[10][3];

  /*if(nev(0) ==  33583611)
  {
  cout << "HERE" << endl;
  }
  */
  for (int iring=0; iring<fqmrnring[ifQtoUse]; iring++) {
    for (int i=0; i<3; i++) {
      fqpos[iring][i]= fqmrpos[ifQtoUse][iring][i];
      fqdir[iring][i] = fqmrdir[ifQtoUse][iring][i];
    }
  }

  double fqdwall_MR=GetDwall(fqpos[i_fQring[0]]);
  double twallfq_MR = GetToWall(fqpos[i_fQring[0]], fqdir[i_fQring[0]]);
  double twallpipm = -1;

  int iMEpipRng = GetfQMER(ifQtoUse,ipip, i_fQring);
  int iMEmuRng = GetfQMER(ifQtoUse,imu,i_fQring);

  if(iMEpipRng >= 0)
    {
      twallpipm = GetToWall(fqpos[iMEpipRng],fqdir[iMEpipRng]);
    }
  if(iMEmuRng >= 0)
    {
      if(iMEpipRng < 0 ||( fqmrmom[ifQtoUse][iMEmuRng] > fqmrmom[ifQtoUse][iMEpipRng]))
	{
	  twallpipm =GetToWall(fqpos[iMEmuRng],fqdir[iMEmuRng]);
	}
    }
  if(twallpipm > 0)
    {
      twallfq_MR = TMath::Min(twallfq_MR,twallpipm);
    }
  dwall = fqdwall_MR;
  towall = twallfq_MR;

  /*
  os->MERmom = fqmrmom[ifQtoUse][i_fQring[0]];
  os->SERmom=0;
  os->SERDwall = 0;
  if (fqmrnring[ifQtoUse]>1)
    {
      os->SERmom = fqmrmom[ifQtoUse][i_fQring[1]];
      os->SERDwall = GetDwall(fqpos[i_fQring[1]]);
    }
  */
}

double GetTransMom(int ifit, int * i_fQring) {// Calculate the transverse momentum / total energy factor

  //      int i_fQring[20];
  //      SortfQRing(ifit,i_fQring,byMom);

  int iMER = i_fQring[0];

  double ltransmom = 0.;
  for (int iring=0; iring<fqmrnring[ifit]; iring++) {

    double costmp = fqmrdir[ifit][iMER][0]*fqmrdir[ifit][iring][0] + fqmrdir[ifit][iMER][1]*fqmrdir[ifit][iring][1] + fqmrdir[ifit][iMER][2]*fqmrdir[ifit][iring][2];


    if (costmp>1.) costmp = 1.;
    if (costmp<-1.) costmp = -1.;

    double sintmp = sqrt(1.-costmp*costmp);

    ltransmom = ltransmom + GetRingEvis(ifit,iring)*sintmp;
  }
  
  return ltransmom/GetEvistot(ifit);
}


