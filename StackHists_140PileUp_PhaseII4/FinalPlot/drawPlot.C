#include <cstdio>
using namespace std;

drawPlot(string cutname="RA2nocut", string histname="NJet"){


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////Some cosmetic work for official documents. 
 gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  gROOT->LoadMacro("CMS_lumi_v2.C");

  writeExtraText = false;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"

  int iPeriod = 14;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 14= PU=140,14TeV 
  int iPos = 33;       // left-aligned
  
  int W = 600;
  int H = 600;

  // 
  // Simple example of macro: plot with CMS name and lumi text
  //  (this script does not pretend to work in all configurations)
  // iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV) 
  // For instance: 
  //               iPeriod = 3 means: 7 TeV + 8 TeV
  //               iPeriod = 7 means: 7 TeV + 8 TeV + 13 TeV 
  // Initiated by: Gautier Hamel de Monchenault (Saclay)
  //
  int H_ref = 600; 
  int W_ref = 800; 

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




  //
  TFile * BJfile =new TFile("../PhaseII4_BJ_14TEV_140PileUp.root","R");
  TFile * TTfile =new TFile("../PhaseII4_TT_14TEV_140PileUp.root","R");
  TFile * Bfile =new TFile("../PhaseII4_B_14TEV_140PileUp.root","R");  
  TFile * BBfile =new TFile("../PhaseII4_BB_14TEV_140PileUp.root","R");
  TFile * BBBfile =new TFile("../PhaseII4_BBB_14TEV_140PileUp.root","R");
  TFile * BJJfile =new TFile("../PhaseII4_BJJ_14TEV_140PileUp.root","R");
  TFile * Hfile =new TFile("../PhaseII4_H_14TEV_140PileUp.root","R");
  TFile * LLfile =new TFile("../PhaseII4_LL_14TEV_140PileUp.root","R");
  TFile * LLBfile =new TFile("../PhaseII4_LLB_14TEV_140PileUp.root","R");
  TFile * TBfile =new TFile("../PhaseII4_TB_14TEV_140PileUp.root","R");
  TFile * TJfile =new TFile("../PhaseII4_TJ_14TEV_140PileUp.root","R");
//  TFile * TTBfile =new TFile("../PhaseII4_TTB_14TEV_140PileUp.root","R");
//
  THStack * tempstack;
  THStack * finalstack = new THStack("finalstack","");
  TH1D * temphist, * VJhist, * Thist, * Otherhist;
  char tempname[200];
  TCanvas * c1 = new TCanvas("c1", "c1", 800, 800);
 
  Float_t legendX1 = .60;  //.50;
  Float_t legendX2 = .90;  //.70;
  Float_t legendY1 = .45;  //.65;
  Float_t legendY2 = .85;
  TLegend* catLeg1 = new TLegend(legendX1,legendY1,legendX2,legendY2);
  catLeg1->SetTextSize(0.032);
  catLeg1->SetTextFont(42);


/*
  ///add TTB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)TTBfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetFillColor(7);
  catLeg1->AddEntry(temphist,"t#bar{t} + V","f");
  finalstack->Add(temphist);
*/



/////////All Other Backgrounds

  ///add BBB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BBBfile->Get(tempname)->Clone();
  Otherhist = (TH1D *) tempstack->GetStack()->Last();



  ///add H to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)Hfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  Otherhist->Add(temphist);


  ///add LL to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)LLfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  Otherhist->Add(temphist);


  ///add LLB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)LLBfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  Otherhist->Add(temphist);
  Otherhist->SetFillColor(4);
  catLeg1->AddEntry(Otherhist,"Other SM","f");
  finalstack->Add(Otherhist);

/////endl of other backgrounds


  ///add BB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BBfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetFillColor(6);
  catLeg1->AddEntry(temphist,"VV","f");
  finalstack->Add(temphist);


  ///add TB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)TBfile->Get(tempname)->Clone();
  Thist = (TH1D *) tempstack->GetStack()->Last();


  ///add TJ to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)TJfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  Thist->Add(temphist);
  Thist->SetFillColor(9);
  catLeg1->AddEntry(Thist,"Single Top","f");
  finalstack->Add(Thist);


  ///add B to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)Bfile->Get(tempname)->Clone();
  VJhist = (TH1D *) tempstack->GetStack()->Last();

  ///add BJ to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BJfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  VJhist->Add(temphist);

  ///add BJJ to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BJJfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  VJhist->Add(temphist);
  VJhist->SetFillColor(46);
  catLeg1->AddEntry(VJhist,"V + Jets","f");
  finalstack->Add(VJhist);


  ///add TTbar to the finalstack
  sprintf(tempname,"TTbar/%s/%s_%s_TTbar",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)TTfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetFillColor(8);
  catLeg1->AddEntry(temphist,"t#bar{t}","f");
  finalstack->Add(temphist);



  ///Draw the background stack
finalstack->Draw();
//  finalstack->SetMinimum(200.);
   //  finalstack->GetHistogram()->GetXaxis()->SetTitle("p_{T}(jet1) [GeV]");
if(histname=="NBtagLoose"){
gPad->SetLogy();
sprintf(tempname,"Number of b-tags");
finalstack->GetXaxis()->SetLimits(0., 10.);
}
else if(histname=="NJet"){
gPad->SetLogy();
sprintf(tempname,"Number of jets");
finalstack->GetXaxis()->SetLimits(0., 20.);finalstack->SetMaximum(100.);
}
else if(histname=="MHT"){sprintf(tempname,"#slash{H}_{T} (GeV)");finalstack->GetXaxis()->SetLimits(500., 3000.);}
else if(histname=="HT"){sprintf(tempname,"H_{T} (GeV)");finalstack->GetXaxis()->SetLimits(1000., 5000.);}
else{sprintf(tempname,"%s",histname.c_str());}
  finalstack->GetHistogram()->GetXaxis()->SetTitle(tempname);
if(histname=="MHT"){sprintf(tempname,"Number of Events / 100 GeV");}
else if(histname=="HT"){sprintf(tempname,"Number of Events / 200 GeV");}
else{sprintf(tempname,"Number of Events");}
  finalstack->GetHistogram()->GetYaxis()->SetTitle(tempname);
//  c1->Modified();

 
  //------------------------------
  // Signal
TFile * Sigfile =new TFile("../PhaseII4_FullExceptStopv4_14TEV_140PileUp.root","R");

  //Draw the signal on the same canvas
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());

  tempstack = (THStack *) Sigfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetLineColor(4);
  temphist->SetLineWidth(2);
  //temphist->SetLineStyle(2);
  temphist->SetFillStyle(0);
  temphist->Draw("SAME");
  std::cout << temphist->GetSumOfWeights() << std::endl;
  catLeg1->AddEntry(temphist,"STOC","l");

  gPad->RedrawAxis();

  //
  catLeg1->SetFillColor(kWhite);
  catLeg1->SetBorderSize(0);
  catLeg1->Draw();


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //For official documents
  // writing the lumi information and the CMS "logo"
  //
  {
    CMS_lumi_v2( c1, iPeriod, iPos );
  }

  c1->Update();
  c1->RedrawAxis();
  c1->GetFrame()->Draw();
sprintf(tempname,"%s_%s.pdf",cutname.c_str(), histname.c_str());
  c1->Print(tempname,".pdf");
sprintf(tempname,"%s_%s.png",cutname.c_str(), histname.c_str());
//  c1->Print(tempname,".png");

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}
