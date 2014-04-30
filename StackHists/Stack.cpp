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

using namespace std;

class mainClass{
double tempvalue;
vector<double> scalevec;
char tempname[200];
char newhistname[200];
vector<double> xs_vec;
map<int, string> cutname; 
map<int, string> type;
map<int, string> histname;
map<string, map<string , vector<TH1D> > > map_map; //here the string shows the type of events and map is between cutname and a vector of histograms like HT, MHT,...
vector< map<string, map<string , vector<TH1D> > >  > vec_map_map; //this is a vector of previous map. Each component of the vector corresponds to one HTbin. 
                                                                  //In BJ case there are 7 of them.
TFile *file;
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
  cutname[5]="RA2delphicut" ;
  cutname[6]="RA2noleptoncut" ;
  cutname[7]="RA2Inc4Jetcut";
  cutname[8]="RA2Inc5Jetcut";
  cutname[9]="RA2Inc6Jetcut" ;
  cutname[10]="RA2allbutHT2500cut";
  cutname[11]="RA2allbutMHT1000cut";
  cutname[12]="RA2allcut";
  cutname[13]="RA2_Asys_allbutHT2500";
  cutname[14]="RA2_Asys_allbutMHT1000";

  type[0]="allEvents";
  type[1]="W";
  type[2]="Wlv";
  type[3]="Wjj";
  type[4]="Z";
  type[5]="Zll";
  type[6]="Zvv";
  type[7]="Zjj";
  type[8]="photon";
  type[9]="H";

  //KH
  histname[0]="weight";
  histname[1]="HT";
  histname[2]="MHT";
  histname[3]="NJet";
  ///end of initialization of the maps

  //build a vector of scale factors
  //first load the cross sections into a vector
  xs_vec.push_back(34409.92339);
  xs_vec.push_back(2642.85309);
  xs_vec.push_back(294.12311);
  xs_vec.push_back(25.95000);
  xs_vec.push_back(2.42111);
  xs_vec.push_back(0.22690);
  xs_vec.push_back(0.02767);
  for(int i=1; i<=7 ; i++){
    sprintf(tempname,"/home/hatake/ana/Delphes/CMSSW_5_2_6/src/DelphesAnalysis/Results/results_PhaseI_BJ_14TEV_HT%d_NoPileUp.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/RA2nocut/HT_RA2nocut_allEvents");
    tempvalue = (luminosity*xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    scalevec.push_back(tempvalue);
  }//end of loop over HTbins 1..7
  std::cout << "normalization scale factor determination done" << std::endl;

  const int nHT = 7;   // Total number of HT bin samples
  const int nHist = 4; // Number of histograms in each TDirectory

  //load the histograms into the vec_map_map  
  for(int i=1; i<=nHT ; i++){    // loop over different HT bins
    map_map.clear();
    //open the file
    sprintf(tempname,"/home/hatake/ana/Delphes/CMSSW_5_2_6/src/DelphesAnalysis/Results/results_PhaseI_BJ_14TEV_HT%d_NoPileUp.root",i);
    file = new TFile(tempname, "R");
    for(map<int , string >::iterator itt=type.begin(); itt!=type.end();itt++){       // loop over different event types
      for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){  // loop over diffrent cut names
	temphistvec.clear();
	vtemphist.clear();
	//for(int j=0; j<nHist; j++){
	for(map<int , string >::iterator ithist=histname.begin(); ithist!=histname.end();ithist++){ // loop over different histograms
	  //KH
	  //KH std::cout << i << " " << (itt->second).c_str() << " " << (it->second).c_str() << " " << (ithist->second).c_str() << std::endl;
	  sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),
		  (ithist->second).c_str(),(it->second).c_str(),(itt->second).c_str());
	  sprintf(newhistname,"%s/%s/%s_%s_%s_HT%d",(itt->second).c_str(),(it->second).c_str(),
		  (ithist->second).c_str(),(it->second).c_str(),(itt->second).c_str(),i);
	  //sprintf(tempname,"allEvents/RA2nocut/allEvents_RA2nocut_hist%d",j);
	  vtemphist.push_back(*(TH1D*) file->Get(tempname));
	  vtemphist.back().SetName(newhistname);//scale the hists before loading them 
	  vtemphist.back().SetTitle(newhistname);//scale the hists before loading them 
	  vtemphist.back().Scale(scalevec[i-1]);//scale the hists before loading them 
	  //KH
	  //KH Print out the mean HT values for different HT samples
	  if ((ithist->second)==histname[1] && (itt->second)==type[0] && (it->second)==cutname[0]){
	    std::cout << i << " " << (itt->second).c_str() << " " << (it->second).c_str() << " " << (ithist->second).c_str() << std::endl;
	    std::cout << vtemphist.back().GetMean() << std::endl;
	    //std::cout << vtemphist[1].GetEntries() << std::endl(); 	
	  }
	  //KH
	  temphistvec.push_back(vtemphist.back());
	  //KH
	  //KH std::cout << i << " " << (itt->second).c_str() << " " << (it->second).c_str() << " " << (ithist->second).c_str() << std::endl;
	}//end of loop over j=0..3
	map_map[itt->second][it->second]=temphistvec;   // vector of histograms - 
      }//end of loop over map
    }//end of loop over map
    vec_map_map.push_back(map_map);                     // vector of histograms - 
  }//end of loop over HTbins 1..7
  std::cout << "loading histograms done" << std::endl;

  tempstack = new THStack("stack","Binned Sample Stack");
 
  file = new TFile("stack.root","RECREATE");
  //stack all 7 HTbins 
  //open a root file to write the stacked hists
  for(map<int , string >::iterator itt=type.begin(); itt!=type.end();itt++){        // loop over different event types
    cdtoitt = file->mkdir((itt->second).c_str());
    cdtoitt->cd();
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
      cdtoit =  cdtoitt->mkdir((it->second).c_str());
      cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
	for(int i=0; i<nHT ; i++){                                                  // loop over different HT bins
	  //std::cout << "i " << i << " j " << j << " itt->second " << itt->second << " it->second " << it->second << std::endl;
	  //KKH temphist = vec_map_map[i][itt->second][it->second][j];
	  temphist = (TH1D*)vec_map_map[i][itt->second][it->second][j].Clone();
	  temphist->SetFillColor(i+2);
	  //KKH tempstack->Add(&temphist);
	  tempstack->Add(temphist);
	  //KH
	  //KH Print out the mean HT values for different HT samples
	  if (j==1 && (itt->second)==type[0] && (it->second)==cutname[0]){
	    std::cout << i << " " << (itt->second).c_str() << " " << (it->second).c_str() << " " << histname[j] << std::endl;
	    std::cout << temphist->GetMean() << std::endl;
	    //std::cout << vtemphist[1].GetEntries() << std::endl(); 	
	    tempstack->GetHists()->Print("all");
	  }
	  //KH
	}//end of loop over HTbins 1..7
	sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
	tempstack->Write(tempname);
	if (j==1 && (itt->second)==type[0] && (it->second)==cutname[0]){
	}
	delete tempstack;
	tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over j=0..3
    }//end of loop over map
  }//end of loop over map

}//end of the constructor
};

int main(){
mainClass mainObj(3000000);

return 0;
}
