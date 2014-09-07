/*
  
  Usage:
  .L plot_significance_vs_lumi.C
  plot_significance_vs_lumi()

*/

double Significance_Gauss(double *x, double *p) {

  Float_t lumi =x[0]; // luminosity
  Double_t S = p[0]*lumi/3000.;                   // Signal yield
  Double_t B = p[1]*lumi/3000.;                   // Background yield
  Double_t deltaB_const = p[2]*lumi/3000.;        // Background uncertainty - fractionally constrant  (linearly increase with lumi)
  Double_t deltaB_scale = p[3]*sqrt(lumi/3000.);  // Background uncertainty - scale down w/ 1/sqrt(L) (scale with sqrt(lumi) )
  Double_t Significance_Gauss = S / sqrt(B + deltaB_const*deltaB_const + deltaB_scale*deltaB_scale) ;
  return Significance_Gauss;

}

double Significance_Zbi(double *x, double *p) {

  Float_t lumi =x[0]; // luminosity
  Double_t S = p[0]*lumi/3000.;                   // Signal yield
  Double_t B = p[1]*lumi/3000.;                   // Background yield
  Double_t deltaB_const = p[2]*lumi/3000.;        // Background uncertainty - fractionally constrant  (linearly increase with lumi)
  Double_t deltaB_scale = p[3]*sqrt(lumi/3000.);  // Background uncertainty - scale down w/ 1/sqrt(L) (scale with sqrt(lumi) )
  Double_t deltaB_tot = pow( pow(deltaB_const,2) + pow(deltaB_scale,2), 0.5);
  Double_t Significance_Zbi = RooStats::NumberCountingUtils::BinomialObsZ(S+B,B,deltaB_tot/B) ;
  return Significance_Zbi;

}

void plot_significance_vs_lumi(){

  //
  // General setting
  //-----------------
  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  gROOT->LoadMacro("CMS_lumi_v2_forLumiPlot.C");

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
  int W_ref = 600; 

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  // Margins:
  //KH tdrStyle->SetPadTopMargin(0.05);
  //tdrStyle->SetPadTopMargin(0.08);
  //tdrStyle->SetPadBottomMargin(0.13);
  //tdrStyle->SetPadLeftMargin(0.16);
  //KH tdrStyle->SetPadRightMargin(0.05);
  //tdrStyle->SetPadRightMargin(0.07);

  // Updated on August 27, 2014
  //double Sig,   Sig_3000   = 4381.63; //2587.0;
  //double Sig,   Sig_3000   = 2127.12; // pt1>1100 GeV
  //double Sig,   Sig_3000   = 3722.46; // pt1>1000 GeV
  //double Sig,   Sig_3000   = 6623.08; // pt1>900 GeV
  //double Sig,   Sig_3000   = 11239.9; // pt1>800 GeV
  double Sig_3000 = 23.732;
  
  // Bkg
  // 1000 GeV
  //double Wlv,   Wlv_3000   = 2656.42; //2602.3;
  //double Ttbar, Ttbar_3000 = 56.6897; //54.2;
  //double Zvv,   Zvv_3000   = 15038.9; //18137.2;
  // 1100 GeV
  //double Wlv,   Wlv_3000   = 1140.62; //
  //double Ttbar, Ttbar_3000 = 27.4762; //
  //double Zvv,   Zvv_3000   = 8209.28; //
  // 900 GeV
  //double Wlv,   Wlv_3000   = 4990.07; //
  //double Ttbar, Ttbar_3000 = 137.919; //
  //double Zvv,   Zvv_3000   = 28385.6; //
  // 800 GeV
  //double Wlv,   Wlv_3000   = 10427.6; //
  //double Ttbar, Ttbar_3000 = 359.675; //
  //double Zvv,   Zvv_3000   = 55917.8; //

  double Tot_Bkg_3000 = 4.631;

  // Uncertainty
  double BkgError_const_3000 = 1.389;
  double BkgError_scale_3000 = 0.0;

   // Significance vs Lumi

  TF1 *sig_1000GeV_gauss = new TF1("significance",Significance_Gauss,100.,3000.,4);
  sig_1000GeV_gauss->SetParameters(Sig_3000, Tot_Bkg_3000, BkgError_const_3000, 0.);

  TF1 *sig_1000GeV_zbi = new TF1("significance",Significance_Zbi,100.,3000.,4);
  sig_1000GeV_zbi->SetParameters(Sig_3000, Tot_Bkg_3000, BkgError_const_3000, 0.);

  //
  std::cout << "3000/fb" << std::endl;
  std::cout << "Signal:        " << Sig_3000            << std::endl
	    << "Bkg:           " << Tot_Bkg_3000        << std::endl
	    << "BkgErr(const): " << BkgError_const_3000 << " " << BkgError_const_3000/Tot_Bkg_3000 << "(/totalBg)" << std::endl
	    << "BkgErr(scale): " << BkgError_scale_3000 << " " << BkgError_scale_3000/Tot_Bkg_3000 << "(/totalBg)" << std::endl << std::endl;

  //
  std::cout << "300/fb" << std::endl;
  std::cout << "Signal:        " << Sig_3000/10.        << std::endl
	    << "Bkg:           " << Tot_Bkg_3000/10.    << std::endl
	    << "BkgErr(const): " << BkgError_const_3000/10.       << " " << BkgError_const_3000/Tot_Bkg_3000          << "(/totalBg)" << std::endl
	    << "BkgErr(scale): " << BkgError_scale_3000/sqrt(10.) << " " << BkgError_scale_3000/Tot_Bkg_3000*sqrt(10) << "(/totalBg)" << std::endl << std::endl;

  double lumis[10] = {100,200,300,500,700,1000,1500,2000,2500,3000};
  double combine[10] = {1.3424, 1.8764, 2.2714, 2.8646, 3.3115, 3.8238, 4.4283, 4.8478, 5.154, 5.3858};
  TGraph *sig_1000GeV_combine = new TGraph(10,lumis,combine);

  double x[1]; x[0]=3000.;
  double p[3]; p[0]=Sig_3000;  p[1]=Tot_Bkg_3000;  p[2]=BkgError_const_3000;  p[3]=BkgError_scale_3000; 
  std::cout << "GaussianEstimator: " << Significance_Gauss(x,p) << std::endl;
  std::cout << "GaussianZbi:       " << Significance_Zbi(x,p)   << std::endl;

  // 
  TH1F* tbox = new TH1F("tbox","tbox",100.,100.,3100.);
  tbox->SetMinimum(1.);
  tbox->SetMaximum(7.);
  tbox->GetXaxis()->SetTitle("Luminosity (fb^{-1})");
  tbox->GetYaxis()->SetTitleOffset(1);
  tbox->GetYaxis()->SetTitle("Significance");


  //
  // Significance
  //-------------------
  TString canvName = "plot_significance_vs_lumi_comparison";
  TCanvas* canv1 = new TCanvas(canvName,canvName,10,10,W,H);
  canv1->SetFillColor(0);
  canv1->SetBorderMode(0);
  canv1->SetFrameFillStyle(0);
  canv1->SetFrameBorderMode(0);
  canv1->SetLeftMargin( L/W );
  canv1->SetRightMargin( R/W );
  canv1->SetTopMargin( T/H );
  canv1->SetBottomMargin( B/H );
  canv1->SetTickx(0);
  canv1->SetTicky(0);

  tbox->Draw();
  
  sig_1000GeV_gauss->SetLineColor(2);  
  sig_1000GeV_gauss->SetLineWidth(3);  
  sig_1000GeV_gauss->Draw("same");  

  sig_1000GeV_zbi->SetLineColor(4);  
  sig_1000GeV_zbi->SetLineWidth(3);  
  sig_1000GeV_zbi->Draw("same");  

  sig_1000GeV_combine->SetLineColor(6);  
  sig_1000GeV_combine->SetLineWidth(3);  
  sig_1000GeV_combine->Draw("same");  

  /*
  std::string text = "CMS preliminary  #sqrt{s}=14TeV";
  TLatex *cmsprel = new TLatex(0.18,0.93,text.c_str());
  cmsprel->SetNDC();
  cmsprel->SetTextSize (0.045);
  cmsprel->SetTextFont(42);
  cmsprel->Draw();
  */

  std::string text2 = "Search region: p_{T}^{jet1} > 1 TeV";
  TLatex *label2 = new TLatex(0.55,0.35,text2.c_str());
  label2->SetNDC();
  label2->SetTextSize (0.045);
  label2->SetTextFont(42);
  //label2->Draw();

  Double_t xl1=.15, yl1=0.75, xl2=xl1+.3, yl2=yl1+.125;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->AddEntry(sig_1000GeV_gauss,"S/sqrt(B+#DeltaB^{2})","l");   // h1 and h2 are histogram pointers
  leg->AddEntry(sig_1000GeV_zbi,"Zbi","l");
  leg->AddEntry(sig_1000GeV_combine,"Combine","l");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->Draw();

  //
  // writing the lumi information and the CMS "logo"
  //
  {
    CMS_lumi_v2( canv1, iPeriod, iPos );
  }

  canv1->Update();
  canv1->RedrawAxis();
  canv1->GetFrame()->Draw();

  canv1->Print(canvName+".pdf",".pdf");
  canv1->Print(canvName+".png",".png");

  //
  // Significance
  //-------------------
  TString canvName = "plot_significance_vs_lumi";
  TCanvas* canv2 = new TCanvas(canvName,canvName,10,10,W,H);
  canv2->SetFillColor(0);
  canv2->SetBorderMode(0);
  canv2->SetFrameFillStyle(0);
  canv2->SetFrameBorderMode(0);
  canv2->SetLeftMargin( L/W );
  canv2->SetRightMargin( R/W );
  canv2->SetTopMargin( T/H );
  canv2->SetBottomMargin( B/H );
  canv2->SetTickx(0);
  canv2->SetTicky(0);

  tbox->Draw();
  
  // sig_1000GeV_gauss->Draw("same");  
  sig_1000GeV_zbi->SetRange(300.,3000.);
  sig_1000GeV_zbi->Draw("same,R");  
  //sig_1000GeV_combine->Draw("same");  

  /*
  std::string text = "CMS preliminary  #sqrt{s}=14TeV";
  TLatex *cmsprel = new TLatex(0.18,0.93,text.c_str());
  cmsprel->SetNDC();
  cmsprel->SetTextSize (0.045);
  cmsprel->SetTextFont(42);
  cmsprel->Draw();
  */

  std::string text2 = "Search region: p_{T}^{jet1} > 1 TeV";
  TLatex *label2 = new TLatex(0.45,0.35,text2.c_str());
  label2->SetNDC();
  label2->SetTextSize (0.042);
  label2->SetTextFont(42);
  //label2->Draw();

  //
  // writing the lumi information and the CMS "logo"
  //
  {
    CMS_lumi_v2( canv2, iPeriod, iPos );
  }

  canv2->Update();
  canv2->RedrawAxis();
  canv2->GetFrame()->Draw();

  canv2->Print(canvName+".pdf",".pdf");
  canv2->Print(canvName+".png",".png");



}

