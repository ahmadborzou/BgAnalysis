#include <cstdio>
using namespace std;

////Define a class to calculate the significance
double zbi(double n_on, double mu_b_hat, double sigma_b){

  //double n_on     = 140.;                         // total events in signal region (S+B)
  //double mu_b_hat = 83.33;                        // mean num of BG events expected in sig. region
  //double sigma_b  = 8.333;                        // uncertainty of mu_b_hat

  double tau      = mu_b_hat / (sigma_b*sigma_b); // scale factor to corresp. Noff/Non              
  double n_off    = tau*mu_b_hat;
  double P_Bi     = TMath::BetaIncomplete(1./(1.+tau), n_on, n_off+1);
  double Z_Bi     = sqrt(2.)*TMath::ErfInverse(1 - 2.*P_Bi);
/*
  cout  <<"  total events in signal region (S+B)               - n_on     " <<n_on      <<endl
        <<"  mean num of BG events expected in sig. region     - mu_b_hat " <<mu_b_hat  <<endl
        <<"  uncertainty of mu_b_hat                           - sigma_b  " <<sigma_b   <<endl
        <<"  scale factor to corresp. Noff/Non                 - tau      " <<tau       <<endl
        <<"  tau*mu_b_hat                                      - n_off    " <<n_off     <<endl
        <<"  TMath::BetaIncomplete(1./(1.+tau), n_on, n_off+1) - P_Bi     " <<P_Bi      <<endl
        <<"  sqrt(2.)*TMath::ErfInverse(1 - 2.*P_Bi)           - Z_Bi     " <<Z_Bi      <<endl;
*/
  return Z_Bi;
}

CheckSignif(string cutname="RA2nocut", string histname="NJet", double xmin=0.0, double lumi=3000.){

  //
  TFile * BJfile  =new TFile("../PhaseII4_BJ_14TEV_140PileUp.root","R");
  TFile * TTfile  =new TFile("../PhaseII4_TT_14TEV_140PileUp.root","R");
  TFile * Bfile   =new TFile("../PhaseII4_B_14TEV_140PileUp.root","R");  
  TFile * BBfile  =new TFile("../PhaseII4_BB_14TEV_140PileUp.root","R");
  TFile * BBBfile =new TFile("../PhaseII4_BBB_14TEV_140PileUp.root","R");
  TFile * BJJfile =new TFile("../PhaseII4_BJJ_14TEV_140PileUp.root","R");
  TFile * Hfile   =new TFile("../PhaseII4_H_14TEV_140PileUp.root","R");
  TFile * LLfile  =new TFile("../PhaseII4_LL_14TEV_140PileUp.root","R");
  TFile * LLBfile =new TFile("../PhaseII4_LLB_14TEV_140PileUp.root","R");
  TFile * TBfile  =new TFile("../PhaseII4_TB_14TEV_140PileUp.root","R");
  TFile * TJfile  =new TFile("../PhaseII4_TJ_14TEV_140PileUp.root","R");
  //  TFile * TTBfile =new TFile("../PhaseII4_TTB_14TEV_140PileUp.root","R");
  //
  THStack * tempstack;
  THStack * finalstack = new THStack("finalstack","Final Plot");
  TH1D * temphist;
  char tempname[200];

  ///add Wlv to the finalstack
  sprintf(tempname,"Wlv/%s/%s_%s_Wlv",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BJfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
if(temphist->Integral()!=temphist->GetSumOfWeights())printf("+++++\n+++++\nWARNING\nWARNING\n++++\n+++++\nIntegral is not the same as sum of weights");
printf("Wlv: +++++ Integral From Xmin: %f, Integral: %f, SumofWeights: %f \n",temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast()),temphist->Integral(), temphist->GetSumOfWeights());  
double Nwlv= temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast());


  ///add TTbar to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)TTfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
if(temphist->Integral()!=temphist->GetSumOfWeights())printf("+++++\n+++++\nWARNING\nWARNING\n++++\n+++++\nIntegral is not the same as sum of weights");
printf("TTbar: +++++ Integral From Xmin: %f, Integral: %f, SumofWeights: %f \n",temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast()),temphist->Integral(), temphist->GetSumOfWeights());
double Ntt= temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast());

  ///add Zvv to the finalstack
  sprintf(tempname,"Zvv/%s/%s_%s_Zvv",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BJfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
if(temphist->Integral()!=temphist->GetSumOfWeights())printf("+++++\n+++++\nWARNING\nWARNING\n++++\n+++++\nIntegral is not the same as sum of weights");
printf("Zvv: +++++ Integral From Xmin: %f, Integral: %f, SumofWeights: %f \n",temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast()),temphist->Integral(), temphist->GetSumOfWeights());
double Nzvv= temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast());

  ///do NOT add BJ allEvents, but check the number of events
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BJfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
if(temphist->Integral()!=temphist->GetSumOfWeights())printf("+++++\n+++++\nWARNING\nWARNING\n++++\n+++++\nIntegral is not the same as sum of weights");
printf("BJ: +++++ Integral From Xmin: %f, Integral: %f, SumofWeights: %f \n",temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast()),temphist->Integral(), temphist->GetSumOfWeights());
double Nbj= temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast());

  ///add B to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)Bfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
if(temphist->Integral()!=temphist->GetSumOfWeights())printf("+++++\n+++++\nWARNING\nWARNING\n++++\n+++++\nIntegral is not the same as sum of weights");
printf("B: +++++ Integral From Xmin: %f, Integral: %f, SumofWeights: %f \n",temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast()),temphist->Integral(), temphist->GetSumOfWeights());
double Nb= temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast());


  ///add BB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BBfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
if(temphist->Integral()!=temphist->GetSumOfWeights())printf("+++++\n+++++\nWARNING\nWARNING\n++++\n+++++\nIntegral is not the same as sum of weights");
printf("BB: +++++ Integral From Xmin: %f, Integral: %f, SumofWeights: %f \n",temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast()),temphist->Integral(), temphist->GetSumOfWeights());
double Nbb= temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast());


  ///add BBB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BBBfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
if(temphist->Integral()!=temphist->GetSumOfWeights())printf("+++++\n+++++\nWARNING\nWARNING\n++++\n+++++\nIntegral is not the same as sum of weights");
printf("BBB: +++++ Integral From Xmin: %f, Integral: %f, SumofWeights: %f \n",temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast()),temphist->Integral(), temphist->GetSumOfWeights());
double Nbbb= temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast());

  ///add BJJ to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BJJfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
if(temphist->Integral()!=temphist->GetSumOfWeights())printf("+++++\n+++++\nWARNING\nWARNING\n++++\n+++++\nIntegral is not the same as sum of weights");
printf("BJJ: +++++ Integral From Xmin: %f, Integral: %f, SumofWeights: %f \n",temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast()),temphist->Integral(), temphist->GetSumOfWeights());
double Nbjj= temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast());

  ///add H to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)Hfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
if(temphist->Integral()!=temphist->GetSumOfWeights())printf("+++++\n+++++\nWARNING\nWARNING\n++++\n+++++\nIntegral is not the same as sum of weights");
printf("H: +++++ Integral From Xmin: %f, Integral: %f, SumofWeights: %f \n",temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast()),temphist->Integral(), temphist->GetSumOfWeights());
double Nh= temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast());


  ///add LL to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)LLfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
if(temphist->Integral()!=temphist->GetSumOfWeights())printf("+++++\n+++++\nWARNING\nWARNING\n++++\n+++++\nIntegral is not the same as sum of weights");
printf("LL: +++++ Integral From Xmin: %f, Integral: %f, SumofWeights: %f \n",temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast()),temphist->Integral(), temphist->GetSumOfWeights());
double Nll= temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast());


  ///add LLB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)LLBfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
if(temphist->Integral()!=temphist->GetSumOfWeights())printf("+++++\n+++++\nWARNING\nWARNING\n++++\n+++++\nIntegral is not the same as sum of weights");
printf("LLB: +++++ Integral From Xmin: %f, Integral: %f, SumofWeights: %f \n",temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast()),temphist->Integral(), temphist->GetSumOfWeights());
double Nllb= temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast());


  ///add TB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)TBfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
if(temphist->Integral()!=temphist->GetSumOfWeights())printf("+++++\n+++++\nWARNING\nWARNING\n++++\n+++++\nIntegral is not the same as sum of weights");
printf("TB: +++++ Integral From Xmin: %f, Integral: %f, SumofWeights: %f \n",temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast()),temphist->Integral(), temphist->GetSumOfWeights());
double Ntb= temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast());


  ///add TJ to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)TJfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
if(temphist->Integral()!=temphist->GetSumOfWeights())printf("+++++\n+++++\nWARNING\nWARNING\n++++\n+++++\nIntegral is not the same as sum of weights");
printf("TJ: +++++ Integral From Xmin: %f, Integral: %f, SumofWeights: %f \n ",temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast()),temphist->Integral(), temphist->GetSumOfWeights());
double Ntj= temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast());


  // Signal
  //--------
  TFile * Sigfile =new TFile("../PhaseII4_FullExceptStopv4_14TEV_140PileUp.root","R");

  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());

  tempstack = (THStack *) Sigfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  if(temphist->Integral()!=temphist->GetSumOfWeights())printf("+++++\n+++++\nWARNING\nWARNING\n++++\n+++++\nIntegral is not the same as sum of weights");
  printf("Signal: +++++ Integral From Xmin: %f, Integral: %f, SumofWeights: %f \n ",temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast()),temphist->Integral(), temphist->GetSumOfWeights());
  double Nsig= temphist->Integral(temphist->GetXaxis()->FindBin(xmin),temphist->GetXaxis()->GetLast());
 
  //KH
  //double totalBG=Nwlv+Nzvv+Ntt+Nb+Nbb+Nbbb+Nh+Nll+Nllb+Nbjj+Ntb+Ntj;
  double totalBG=Nbj+Ntt+Nb+Nbb+Nbbb+Nh+Nll+Nllb+Nbjj+Ntb+Ntj;
  double delBG=0.3*totalBG;
  double delBGsquare=pow(delBG,2);

  Nsig    = Nsig * lumi/3000.;
  totalBG = totalBG * lumi/3000.;
  delBG   = delBG * lumi/3000.;

  //KH 
  //double significance=Nsig/sqrt(Nsig+totalBG+delBGsquare);
  double significance=Nsig/sqrt(totalBG+delBGsquare);
  
  double zbii = zbi((totalBG+Nsig),totalBG ,sqrt(delBGsquare));

  printf("\n\nSignificance: %f, ZBI: %f \n============================================================\n", significance, zbii);
  printf("S: %8.3f, B: %8.3f, deltaB: %8.3f\n",Nsig,totalBG,delBG);

  
  
}
