//There some points here that should be cared about
//First, cutname must be initialized the same as the one in main.cpp. 
//The same thing is true about type. Also,  
//there is a loop in the part that loads histograms into vec_map_map. 
//Number of loops should be the same as size of the vector of histograms
//in main.cpp, i.e. number of histograms in each branch. Currently,
//we have 4 histograms Weight,HT,MHT, NJets 

#include <cassert>
#include "TH1.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include "TSystem.h"
#include "THStack.h"
#include "TFile.h"
#include <vector>
#include <map>
#include "TDirectory.h"
#include "TCanvas.h"
#include "TAxis.h"

using namespace std;

class mainClass{
double tempvalue;
vector<double> Sig_scalevec, BJ_scalevec, TT_scalevec;
char tempname[200];
char newhistname[200];
vector<double> Sig_xs_vec,BJ_xs_vec,TT_xs_vec;
map<int, string> cutname; 
map<int, string> sigtype, BJtype, TTtype;
map<int, string> histname;
map<int, vector<double> > yieldmap; //the first component must correspond to the components of the cutname map. For instance yieldmap[0]=0 would correspond to nocut while yieldmap[0]=1 to Asys. The other components represnt the # of events of Wlv, Zvv, TTbar, and Signal.
map<string, map<string , vector<TH1D> > > map_map; //here the string shows the type of events and map is between cutname and a vector of histograms like HT, MHT,...
vector< map<string, map<string , vector<TH1D> > >  > vec_map_map; //this is a vector of previous map. Each component of the vector corresponds to one HTbin. 
                                                                  //In BJ case there are 7 of them.
TFile *file;
vector<TFile *> Sig_inputfilevec, BJ_inputfilevec, TT_inputfilevec;
TH1D *temphist;
//KKH TH1D temphist;
vector<TH1D > vtemphist;
vector<TH1D > temphistvec;
THStack * tempstack;
TDirectory *cdtoitt;
TDirectory *cdtoit;


public:
mainClass(int luminosity){//constructor
//Importnat
//make sure this initialization of the 
//maps is the same as that in main.cpp
     cutname[0]="RA2nocut";
    cutname[1]="RA2Asys";
    cutname[2]="RA2Inc3Jetcut";
    cutname[3]="RA2HT500cut";
    cutname[4]="RA2MHT200cut";
    cutname[5]="RA2delphicut";
    cutname[6]="RA2noleptoncut";
    cutname[7]="noPhotoncut";
    cutname[8]="RA2Inc4Jetcut";
    cutname[9]="RA2Inc5Jetcut";
    cutname[10]="RA2Inc6Jetcut";
    cutname[11]="RA2allbutHT2500cut";
    cutname[12]="RA2allbutMHT1000cut";
    cutname[13]="RA2allcut";
    cutname[14]="RA2noleptoncutMHT1000";
    cutname[15]="RA2noleptoncutBtag2";
    cutname[16]="RA2noleptoncutBtag2MHT1000";
    cutname[17]="RA2Inc4JetcutMHT1000";
    cutname[18]="RA2Inc4JetcutBtag2";
    cutname[19]="RA2Inc4JetcutBtag2MHT1000";
    cutname[20]="RA2Inc5JetcutMHT1000";
    cutname[21]="RA2Inc5JetcutBtag2";
    cutname[22]="RA2Inc5JetcutBtag2MHT1000";
    cutname[23]="RA2Inc6JetcutMHT1000";
    cutname[24]="RA2Inc6JetcutBtag2";
    cutname[25]="RA2Inc6JetcutBtag2MHT1000";


    sigtype[0]="allEvents";
    sigtype[1]="glgl";

    BJtype[0]="allEvents";
    BJtype[1]="W";
    BJtype[2]="Wlv";
    BJtype[3]="Wjj";
    BJtype[4]="Z";
    BJtype[5]="Zll";
    BJtype[6]="Zvv";
    BJtype[7]="Zjj";
    BJtype[8]="photon";
    BJtype[9]="H";

   TTtype[0]="allEvents";
   TTtype[1]="TTbar";
   TTtype[2]="TTSingLep";
   TTtype[3]="TTdiLep";
   TTtype[4]="TThadronic";



  //KH
  histname[0]="weight";
  histname[1]="HT";
  histname[2]="MHT";
  histname[3]="NJet";
  histname[4]="NBtagLoose";
  histname[5]="NBtagTight";
  histname[6]="BtagLoose1Pt";
histname[7]="BtagLoose1Eta";
histname[8]="BtagLoose1Phi";
histname[9]="BtagLoose2Pt";
histname[10]="BtagLoose2Eta";
histname[11]="BtagLoose2Phi";
histname[12]="BtagTight1Pt";
histname[13]="BtagTight1Eta";
histname[14]="BtagTight1Phi";
histname[15]="BtagTight2Pt";
histname[16]="BtagTight2Eta";
histname[17]="BtagTight2Phi";

  ///end of initialization of the maps

yieldmap.clear();



//Signal Section//Signal Section//Signal Section//Signal Section//Signal Section//Signal Section//Signal Section//Signal Section

  //build a vector of scale factors
  //first load the cross sections into a vector
//Sig_xs_vec.push_back(0.757); /// v1
//Sig_xs_vec.push_back(1.12); // v2
//Sig_xs_vec.push_back(1.15); // v3
//Sig_xs_vec.push_back(1.14); // M(Stop,LSP)=(450,410) and also M(Stop,LSP)=(450,440)
//Sig_xs_vec.push_back(2.18); // M(Stop,LSP)=(400,390) and also M(Stop,LSP)=(400,360)
//Sig_xs_vec.push_back(4.41); // M(Stop,LSP)=(350,340) and also M(Stop,LSP)=(350,310)
//Sig_xs_vec.push_back(0.009635); //STOCv4
Sig_xs_vec.push_back(1.58); //StauC

  double Sig_numberofevents =0;//this will use GetSumOfWeights() 
  const int Sig_nHT = 1;   // Total number of HT bin samples
  const int nHist = 18; // Number of histograms in each TDirectory


  for(int i=1; i<=Sig_nHT ; i++){
//sprintf(tempname,"../Results/results_PhaseII4_Stop_CharmLSP_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_Stop_CharmLSPv2_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_Stop_CharmLSPv3_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc450410_14TEV_140PileUp_00.root");  
//sprintf(tempname,"../Results/results_PhaseII4_t2cc450440_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc400390_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc400360_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc350340_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc350310_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_Stop_CharmLSPv4_14TEV_140PileUp_00.root");
sprintf(tempname,"../Results/results_PhaseII4_StauC_14TEV_140PileUp.root");
  file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/RA2nocut/MHT_RA2nocut_allEvents");
    tempvalue = (luminosity*Sig_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    Sig_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=Sig_nHT; i++){
//sprintf(tempname,"../Results/results_PhaseII4_Stop_CharmLSP_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_Stop_CharmLSPv2_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_Stop_CharmLSPv3_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc450410_14TEV_140PileUp_00.root");  
//sprintf(tempname,"../Results/results_PhaseII4_t2cc450440_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc400390_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc400360_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc350340_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc350310_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_Stop_CharmLSPv4_14TEV_140PileUp_00.root");
sprintf(tempname,"../Results/results_PhaseII4_StauC_14TEV_140PileUp.root");
Sig_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

tempstack = new THStack("stack","Binned Sample Stack");
//sprintf(tempname,"PhaseII4_Stop_CharmLSP_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_Stop_CharmLSPv2_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_Stop_CharmLSPv3_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_t2cc450410_14TEV_140PileUp_00.root");  
//sprintf(tempname,"PhaseII4_t2cc450440_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_t2cc400390_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_t2cc400360_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_t2cc350340_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_t2cc350310_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_Stop_CharmLSPv4_14TEV_140PileUp_00.root");
sprintf(tempname,"PhaseII4_StauC_14TEV_140PileUp.root");
file = new TFile(tempname,"RECREATE");
 for(map<int , string >::iterator itt=sigtype.begin(); itt!=sigtype.end();itt++){        // loop over different event types
    cdtoitt = file->mkdir((itt->second).c_str());
    cdtoitt->cd();
int c=0;
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
      cdtoit =  cdtoitt->mkdir((it->second).c_str());
      cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<Sig_nHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) Sig_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(Sig_scalevec[i]);
if(histname[j]=="MHT"){
Sig_numberofevents+=(double)temphist->GetSumOfWeights();
} 
temphist->SetFillColor(i+2);
tempstack->Add(temphist);


               }//end of loop over HTbins 1..7
if(histname[j]=="MHT"){
if(itt->second=="allEvents"){
yieldmap[c].push_back(Sig_numberofevents);
}

}
Sig_numberofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
        tempstack->Write(tempname);
        delete tempstack;
        tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
 c+=1;   }//end of loop over cutnames
  }//end of loop over event types
file->Close();

//BJ Section//BJ Section//BJ Section//BJ Section//BJ Section//BJ Section//BJ Section//BJ Section//BJ Section//BJ Section//BJ Section

 //build a vector of scale factors
  //first load the cross sections into a vector
  BJ_xs_vec.push_back(34409.92339);
  BJ_xs_vec.push_back(2642.85309);
  BJ_xs_vec.push_back(294.12311);
  BJ_xs_vec.push_back(25.95000);
  BJ_xs_vec.push_back(2.42111);
  BJ_xs_vec.push_back(0.22690);
  BJ_xs_vec.push_back(0.02767);

  double BJ_numberofevents =0;
  const int bjnHT = 7;   // Total number of HT bin samples


  for(int i=1; i<=bjnHT ; i++){
    sprintf(tempname,"../Results/results_PhaseII4_BJ_14TEV_HT%d_140PileUp.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/RA2nocut/MHT_RA2nocut_allEvents");
    tempvalue = (luminosity*BJ_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    BJ_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=bjnHT; i++){
sprintf(tempname,"../Results/results_PhaseII4_BJ_14TEV_HT%d_140PileUp.root",i);
BJ_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

//tempstack = new THStack("stack","Binned Sample Stack");
sprintf(tempname,"PhaseII4_BJ_14TEV_140PileUp.root");
//file = new TFile(tempname,"RECREATE");
 for(map<int , string >::iterator itt=BJtype.begin(); itt!=BJtype.end();itt++){        // loop over different event types
//    cdtoitt = file->mkdir((itt->second).c_str());
//    cdtoitt->cd();
int c=0;
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
//      cdtoit =  cdtoitt->mkdir((it->second).c_str());
//      cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<bjnHT ; i++){                                                  // loop over different HT bins
//cout << "================================" << endl;
//cout << "HT#: " <<i << ", BJtype: " << itt->second << ", cutname: " << it->second << ", hist#: " << j << endl;  
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) BJ_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(BJ_scalevec[i]);
if(histname[j]=="MHT"){
BJ_numberofevents+=(double)temphist->GetSumOfWeights();
}
temphist->SetFillColor(i+2);
//tempstack->Add(temphist);
               }//end of loop over HTbins 1..7
if(histname[j]=="MHT"){
if(itt->second=="Wlv" || itt->second=="Zvv"){
yieldmap[c].push_back(BJ_numberofevents);
}

}
BJ_numberofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
  //      tempstack->Write(tempname);
  //      delete tempstack;
  //      tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
c+=1;    }//end of loop over cutnames
  }//end of loop over event types
//file->Close();


//TTbar Section//TTbar Section//TTbar Section//TTbar Section//TTbar Section//TTbar Section//TTbar Section//TTbar Section

  //build a vector of scale factors
  //first load the cross sections into a vector
  TT_xs_vec.push_back(530.89358);
  TT_xs_vec.push_back(42.55351);
  TT_xs_vec.push_back(4.48209);
  TT_xs_vec.push_back(0.52795);
  TT_xs_vec.push_back(0.05449);

  double TT_numberofevents =0;
  const int ttnHT = 5;   // Total number of HT bin samples


  for(int i=1; i<=ttnHT ; i++){
    sprintf(tempname,"../Results/results_PhaseII4_TT_14TEV_HT%d_140PileUp.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/RA2nocut/MHT_RA2nocut_allEvents");
    tempvalue = (luminosity*TT_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    TT_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=ttnHT; i++){
sprintf(tempname,"../Results/results_PhaseII4_TT_14TEV_HT%d_140PileUp.root",i);
TT_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

//tempstack = new THStack("stack","Binned Sample Stack");
sprintf(tempname,"PhaseII4_TT_14TEV_140PileUp.root");
//file = new TFile(tempname,"RECREATE");
 for(map<int , string >::iterator itt=TTtype.begin(); itt!=TTtype.end();itt++){        // loop over different event types
//    cdtoitt = file->mkdir((itt->second).c_str());
//    cdtoitt->cd();
int c=0;
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
//      cdtoit =  cdtoitt->mkdir((it->second).c_str());
//      cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<ttnHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) TT_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(TT_scalevec[i]);
if(histname[j]=="MHT"){
TT_numberofevents+=(double)temphist->GetSumOfWeights();
}
temphist->SetFillColor(i+2);
//tempstack->Add(temphist);
               }//end of loop over HTbins 1..5
if(histname[j]=="MHT"){
if(itt->second=="allEvents"){
yieldmap[c].push_back(TT_numberofevents);
}

}
TT_numberofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
//        tempstack->Write(tempname);
//        delete tempstack;
//        tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
c+=1;   }//end of loop over cutnames
  }//end of loop over event types
//file->Close();



///write the output in a file
fstream ff;
ff.open("CutFlow.txt", std::fstream::out);
ff << " Cut Name,    " << "  Signal,      " << "  Wlv,      " << "  Zvv,     " << "  TTbar,      "<< "    Total BG,   " << " % Signal/Background,   "  <<  "    Significance " << endl; 
double totalBG=0, delWlv=0, delZvv=0, delTT=0, delB=0, delBsquare=0;
for(int i=0; i<yieldmap.size(); i++){
totalBG=(double) (yieldmap[i].at(1)+yieldmap[i].at(2)+yieldmap[i].at(3));
delWlv= 0.08*yieldmap[i].at(1);///uncrtainty for Wlv is 8%
delZvv= 0.05*yieldmap[i].at(2);
delTT= 0.5*yieldmap[i].at(3);///uncrtainty for TTbar is 50%
delBsquare=pow(delWlv,2)+pow(delZvv,2)+pow(delTT,2);///delta_background = sqrt(delWlv^2+delZvv^2+delTT^2)
ff << "  " <<cutname[i]<<",     " << yieldmap[i].at(0) << ",     " << yieldmap[i].at(1) <<",     " << yieldmap[i].at(2) <<",     " <<yieldmap[i].at(3) << ",      "<< totalBG << ",      " << yieldmap[i].at(0)/totalBG*100  <<  ",       " << yieldmap[i].at(0)/sqrt(delBsquare+totalBG+yieldmap[i].at(0))  <<endl;  
}
ff.close();

/*
///write the output in a file
fstream ff;
ff.open("CutFlow.txt", std::fstream::out);
ff << " Cut Name,    " << "  Signal,      "  << endl; 
for(int i=0; i<yieldmap.size(); i++){
ff << "  " <<cutname[i]<<",     " << yieldmap[i].at(0) << endl;  
}
ff.close();
*/
}//end of the constructor
};
int main(){
mainClass mainObj(3000000);
//mainClass mainObj(19700);
cout << " done :) " << endl;

}

