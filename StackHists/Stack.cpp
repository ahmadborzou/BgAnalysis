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
vector<double> xs_vec;
map<int, string> cutname; 
map<int, string> type;
map<string, map<string , vector<TH1D> > > map_map; //here the string shows the type of events and map is between cutname and a vector of histograms like HT, MHT,...
vector< map<string, map<string , vector<TH1D> > >  > vec_map_map; //this is a vector of previous map. Each component of the vector corresponds to one HTbin. 
                                                                  //In BJ case there are 7 of them.
TFile * file;
TH1D temphist;
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
cutname[0]="RA2nocut";cutname[1]="RA23Jetcut";cutname[2]="RA2HT500cut" ;cutname[3]="RA2MHT200cut" ;cutname[4]="RA2delphicut" ;cutname[5]="RA2noleptoncut" ;cutname[6]="RA24Jetcut" ;cutname[7]="RA25Jetcut" ;cutname[8]="RA26Jetcut" ;cutname[9]="RA2allbutHT2500cut" ;cutname[10]="RA2allbutMHT1000cut";cutname[11]= "RA2allcut";cutname[12]="RA2_Asys_allbutHT2500";
cutname[13]="RA2_Asys_allbutMHT1000";

type[0]="allEvents";type[1]="W";type[2]="Wlv";type[3]="Wjj";type[4]="Z";type[5]="Zll";type[6]="Zvv";type[7]="Zjj";type[8]="photon";type[9]="H";
///end of initialization of the maps

//build a vector of scale factors
//first load the cross sections into a vector
xs_vec.push_back(34409.92339);xs_vec.push_back(2642.85309);xs_vec.push_back(294.12311);xs_vec.push_back(25.95000);xs_vec.push_back(2.42111);xs_vec.push_back(0.22690);xs_vec.push_back(0.02767);
for(int i=1; i<=7 ; i++){
sprintf(tempname,"../../BgAnalysisV2/Results/results_PhaseI_BJ_14TEV_HT%d_NoPileUp_00.root",i);
file = new TFile(tempname, "R");
sprintf(tempname,"allEvents/RA2nocut/allEvents_RA2nocut_hist1");
tempvalue = (luminosity*xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
scalevec.push_back(tempvalue);
}//end of loop over HTbins 1..7



//load the histograms into the vec_map_map
for(int i=1; i<=7 ; i++){
map_map.clear();
//open the file
sprintf(tempname,"../../BgAnalysisV2/Results/results_PhaseI_BJ_14TEV_HT%d_NoPileUp_00.root",i);
file = new TFile(tempname, "R");
for(map<int , string >::iterator itt=type.begin(); itt!=type.end();itt++){
for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){
temphistvec.clear();
vtemphist.clear();
for(int j=0; j<=3; j++){
sprintf(tempname,"%s/%s/%s_%s_hist%d",(itt->second).c_str(),(it->second).c_str(),(itt->second).c_str(),(it->second).c_str(),j);
//sprintf(tempname,"allEvents/RA2nocut/allEvents_RA2nocut_hist%d",j);
vtemphist.push_back(*(TH1D* ) file->Get(tempname));
vtemphist.back().Scale(scalevec[i-1]);//scale the hists before loading them 
temphistvec.push_back(vtemphist.back());
}//end of loop over j=0..3
map_map[itt->second][it->second]=temphistvec;
}//end of loop over map
}//end of loop over map
vec_map_map.push_back(map_map);
}//end of loop over HTbins 1..7

tempstack = new THStack("stack","Binned Sample Stack");


file = new TFile("stack.root","RECREATE");
//stack all 7 HTbins 
//open a root file to write the stacked hists
for(map<int , string >::iterator itt=type.begin(); itt!=type.end();itt++){
cdtoitt = file->mkdir((itt->second).c_str());
cdtoitt->cd();
for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){
cdtoit =  cdtoitt->mkdir((it->second).c_str());
cdtoit->cd();
for(int j=0; j<=3; j++){
for(int i=0; i<7 ; i++){

//cout << "i " << i << " j " << j << " itt->second " << itt->second << " it->second " << it->second << endl;
temphist = vec_map_map[i][itt->second][it->second][j];
tempstack->Add( & temphist);
}//end of loop over HTbins 1..7
sprintf(tempname,"%s_%s_hist%d",(itt->second).c_str(),(it->second).c_str(),j);
tempstack->Write(tempname);
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
