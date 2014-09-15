/*
The purpose of this code is to add the k-factors which is the ratio of NLO/LO cross section. Our sample contains different interactions like glgl, slsl, sqsq, etc. They have different kfactors. However, there are some interactions of "unknown" type. Their contribution can be calculated by subtracting all these known types from allEvents. This will be stored in the "other_hist". We then add all the types with kfactor to the "final_hist". At the end we add the "other_hist" to the "final_hist". This final_hist will be stored in the allEvents folder. 


*/
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
#include "TMath.h"

using namespace std;

class mainClass{
double tempvalue;
char tempname[200];
map<int, string> cutname; 
map<int, string> sigtype;
map<int, string> histname;
map<string, double> kfactor;
TFile *file;
vector<TFile *> Sig_inputfilevec;
TH1D *temphist, *other_hist, *final_hist;
//KKH TH1D temphist;
TDirectory *cdtoitt;
TDirectory *cdtoit;


public:
mainClass(int luminosity){//constructor
//Importnat
//make sure this initialization of the 
//maps is the same as that in main.cpp
 
cout << "hello " << endl;

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
    cutname[12]="RA2allbutMHT1300cut";
    cutname[13]="RA2allcut";
    cutname[14]="RA2noleptoncutMHT1300";
    cutname[15]="RA2noleptoncutBtag2";
    cutname[16]="RA2noleptoncutBtag2MHT1300";
    cutname[17]="RA2Inc4JetcutMHT1300";
    cutname[18]="RA2Inc4JetcutBtag2";
    cutname[19]="RA2Inc4JetcutBtag2MHT1300";
    cutname[20]="RA2Inc5JetcutMHT1300";
    cutname[21]="RA2Inc5JetcutBtag2";
    cutname[22]="RA2Inc5JetcutBtag2MHT1300";
    cutname[23]="RA2Inc6JetcutMHT1300";
    cutname[24]="RA2Inc6JetcutBtag2";
    cutname[25]="RA2Inc6JetcutBtag2MHT1300";
cutname[26]="RA2noleptoncutht1500";
cutname[27]="RA2noleptoncutht2000";
cutname[28]="RA2noleptoncutht2500";
cutname[29]="RA2noleptoncutht1500mht700";
cutname[30]="RA2noleptoncutht2000mht800";
cutname[31]="RA2noleptoncutht2000mht1300";
cutname[32]="RA2noleptoncutht2500mht1300";
cutname[33]="RA2noleptoncutbtag2ht1500";
cutname[34]="RA2noleptoncutbtag2ht2000";
cutname[35]="RA2noleptoncutbtag2ht2500";
cutname[36]="RA2noleptoncutbtag2ht1500mht700";
cutname[37]="RA2noleptoncutbtag2ht2000mht800";
cutname[38]="RA2noleptoncutbtag2ht2000mht1300";
cutname[39]="RA2noleptoncutbtag2ht2500mht1300";
/*
cutname[40]="RA2Inc4Jetcutht1500";
cutname[41]="RA2Inc4Jetcutht2000";
cutname[42]="RA2Inc4Jetcutht2500";
cutname[43]="RA2Inc4Jetcutht1500mht700";
cutname[44]="RA2Inc4Jetcutht2000mht800";
cutname[45]="RA2Inc4Jetcutht2000mht1300";
cutname[46]="RA2Inc4Jetcutht2500mht1300";
cutname[47]="RA2Inc4Jetcutbtag2ht1500";
cutname[48]="RA2Inc4Jetcutbtag2ht2000";
cutname[49]="RA2Inc4Jetcutbtag2ht2500";
cutname[50]="RA2Inc4Jetcutbtag2ht1500mht700";
cutname[51]="RA2Inc4Jetcutbtag2ht2000mht800";
cutname[52]="RA2Inc4Jetcutbtag2ht2000mht1300";
cutname[53]="RA2Inc4Jetcutbtag2ht2500mht1300";

cutname[54]="RA2Inc5Jetcutht1500";
cutname[55]="RA2Inc5Jetcutht2000";
cutname[56]="RA2Inc5Jetcutht2500";
cutname[57]="RA2Inc5Jetcutht1500mht700";
cutname[58]="RA2Inc5Jetcutht2000mht800";
cutname[59]="RA2Inc5Jetcutht2000mht1300";
cutname[60]="RA2Inc5Jetcutht2500mht1300";
cutname[61]="RA2Inc5Jetcutbtag2ht1500";
cutname[62]="RA2Inc5Jetcutbtag2ht2000";
cutname[63]="RA2Inc5Jetcutbtag2ht2500";
cutname[64]="RA2Inc5Jetcutbtag2ht1500mht700";
cutname[65]="RA2Inc5Jetcutbtag2ht2000mht800";
cutname[66]="RA2Inc5Jetcutbtag2ht2000mht1300";
cutname[67]="RA2Inc5Jetcutbtag2ht2500mht1300";

cutname[68]="RA2Inc6Jetcutht1500";
cutname[69]="RA2Inc6Jetcutht2000";
cutname[70]="RA2Inc6Jetcutht2500";
cutname[71]="RA2Inc6Jetcutht1500mht700";
cutname[72]="RA2Inc6Jetcutht2000mht800";
cutname[73]="RA2Inc6Jetcutht2000mht1300";
cutname[74]="RA2Inc6Jetcutht2500mht1300";
cutname[75]="RA2Inc6Jetcutbtag2ht1500";
cutname[76]="RA2Inc6Jetcutbtag2ht2000";
cutname[77]="RA2Inc6Jetcutbtag2ht2500";
cutname[78]="RA2Inc6Jetcutbtag2ht1500mht700";
cutname[79]="RA2Inc6Jetcutbtag2ht2000mht800";
cutname[80]="RA2Inc6Jetcutbtag2ht2000mht1300";
cutname[81]="RA2Inc6Jetcutbtag2ht2500mht1300";
*/


sigtype[0]="glgl";
sigtype[1]="sqsq";
sigtype[2]="glsq";
sigtype[3]="slsl";
sigtype[4]="t1t1";
sigtype[5]="b1b1";
sigtype[6]="ewew";
sigtype[7]="Rest";
/*
///FullExceptStop
kfactor["glgl"]=1.5;
kfactor["sqsq"]=1.28;
kfactor["glsq"]=2.64;
kfactor["slsl"]=1;
kfactor["t1t1"]=1;
kfactor["b1b1"]=2.43;
kfactor["ewew"]=1;
assert(kfactor.size()==sigtype.size());
*/


///Stau
kfactor["glgl"]=1.3;
kfactor["sqsq"]=1.21;
kfactor["glsq"]=2.24;
kfactor["slsl"]=1.2;
kfactor["t1t1"]=1.7;
kfactor["b1b1"]=1.75;
kfactor["ewew"]=1;
kfactor["Rest"]=1;
assert(kfactor.size()==sigtype.size());




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
/*histname[12]="BtagTight1Pt";
histname[13]="BtagTight1Eta";
histname[14]="BtagTight1Phi";
histname[15]="BtagTight2Pt";
histname[16]="BtagTight2Eta";
histname[17]="BtagTight2Phi";
*/
  ///end of initialization of the maps


//Signal Section//Signal Section//Signal Section//Signal Section//Signal Section//Signal Section//Signal Section//Signal Section

  const int Sig_nHT = 1;   // Total number of HT bin samples
  const int nHist = (int) histname.size(); // Number of histograms in each TDirectory


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
//sprintf(tempname,"../Results/results_PhaseII4_Stop_CharmLSPv4_14TEV_140PileUp.root");
//sprintf(tempname,"../Results/results_PhaseII4_FullExceptStopv4_14TEV_140PileUp.root");
sprintf(tempname,"../Results/results_PhaseII4_StauC_14TEV_140PileUp.root");
Sig_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

//sprintf(tempname,"PhaseII4_Stop_CharmLSP_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_Stop_CharmLSPv2_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_Stop_CharmLSPv3_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_t2cc450410_14TEV_140PileUp_00.root");  
//sprintf(tempname,"PhaseII4_t2cc450440_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_t2cc400390_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_t2cc400360_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_t2cc350340_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_t2cc350310_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_Stop_CharmLSPv4_14TEV_140PileUp.root");
//sprintf(tempname,"results_PhaseII4_FullExceptStopv4_14TEV_140PileUp.root");
sprintf(tempname,"results_PhaseII4_StauC_14TEV_140PileUp.root");
file = new TFile(tempname,"RECREATE"); 
cdtoitt = file->mkdir("allEvents");
cdtoitt->cd();

    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<Sig_nHT ; i++){                                                  // loop over different HT bins
if(i==0 && j==0){cdtoit =  cdtoitt->mkdir((it->second).c_str());cdtoit->cd();}
else{sprintf(tempname,"allEvents/%s",(it->second).c_str());file->cd(tempname);}

for(map<int , string >::iterator itt=sigtype.begin(); itt!=sigtype.end();itt++){        // loop over different event types

sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) Sig_inputfilevec.at(i)->Get(tempname)->Clone();
if(itt->second=="glgl"){final_hist=temphist;final_hist->Add(temphist,(1.0-kfactor["glgl"]));}
else{final_hist->Add(temphist,kfactor[itt->second]);}
               }//end of loop over event types
        sprintf(tempname,"%s_%s_allEvents",histname[j].c_str(),(it->second).c_str());
        final_hist->Write(tempname);
      }//end of loop over HT Bins. Here, no loop.
    }//end of loop over histograms 
  }//end of loop over cutnames
file->Close();

}//end of the constructor
};


int main(){
mainClass mainObj(3000000);
cout << " done :) " << endl;
}

