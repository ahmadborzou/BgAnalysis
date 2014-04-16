#include "TChain.h"
#include "TH1.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include "TSystem.h"
#include "TClonesArray.h"
#include <vector>
#include <map>
#include "Delphes-3.0.10/external/ExRootAnalysis/ExRootTreeReader.h"
#include "Delphes-3.0.10/classes/DelphesClasses.h"
using namespace std;


class histClass{
double * a;
TH1D * b_hist;
public:
histClass(double * eveinfarr_, TH1D * hist_){
a = eveinfarr_;
b_hist=hist_;
(*b_hist).Fill(*a);
for(int i=1; i<=3 ; i++){
(*(b_hist+i)).Fill(*(a+i),*a);
}
}
};

//define a function to evaluate delta phi
double delphi(vector<double> a, double tPx, double tPy,double mht){
//-totpx is the ptx comp of MHT
double jetpt = sqrt((a[1]*cos(a[2]))*(a[1]*cos(a[2]))+(a[1]*sin(a[2]))*(a[1]*sin(a[2])));
double MHT_Jet_Dot = (-tPx*(a[1]*cos(a[2]))-tPy*(a[1]*sin(a[2])));
double deltaphi = acos(MHT_Jet_Dot/(mht*jetpt));
return deltaphi;
///end of function deltaphi
}

//this function is exclusively written for BJ processes with emphesis on one B.
bool bg_type(string bg_ ,vector<GenParticle*> pvec){

if(bg_=="allEvents"){return 1;}

if(bg_=="Zvv"){
vector<int> vvvec;
for(int i = 0; i < pvec.size(); ++i){
GenParticle * p = pvec.at(i);
if(fabs(p->PID)==23){//23 is the PID code of Z boson.
vvvec.clear();
for(int j = 0; j < pvec.size(); ++j){
GenParticle * pp = pvec.at(j);
//if (pp->Status != 3 || pp->M1 > pvec.size() || pp->M2 > pvec.size() ){continue;}// Assuming the status 3 electron directly from Z. Delphes has broken mother links, can't track back to mother Z
if (pp->Status == 3 && pp->M1 < pvec.size() && pp->M2 < pvec.size() && fabs(pp->PID) == 12 || fabs(pp->PID) == 14 || fabs(pp->PID) == 16 ){
vvvec.push_back(pp->PID);
}//end of if pp->PID == 12, 14, 16 = nutrinos
}//end of second loop
if(vvvec.size()==2){
return true;}//end of if 
}//end of if PID==23=Z boson
}//end of loop
return false;
}//end of if Zvv


if(bg_=="Wlv"){
vector<int> llvec;
for(int i = 0; i < pvec.size(); ++i){
GenParticle * pa = pvec.at(i);
if(fabs(pa->PID)==24){//+-24 are the PID codes of W bosons.
llvec.clear();
for(int j = 0; j < pvec.size(); ++j){
GenParticle * ppa = pvec.at(j);
//if (ppa->Status != 3 || ppa->M1 > pvec.size() || ppa->M2 > pvec.size() ){continue;}// Assuming the status 3 electron directly from W. Delphes has broken mother links, can't track back to mother W
if (ppa->Status == 3 && ppa->M1 < pvec.size() && ppa->M2 < pvec.size() && fabs(ppa->PID) == 11 || fabs(ppa->PID) == 13 || fabs(ppa->PID) == 15 ){
llvec.push_back(ppa->PID);
}//end of if ppa->PID == 11, 13, 15 = electron , muon, tau
}//end of second loop
if(llvec.size()==1){//llvec.size() ==1 since W decays to one lepton and one nutrino
return true;}//end of if 
}//end of if PID==24=W boson
}//end of loop
return false;
}//end of if Wlv


} //end of function bg_type


///////////////////////////////////////////
//Begining of the main()//Begining of the main()//Begining of the main()//Begining of the main()//Begining of the main()//Begining of the main()//Begining of the main()
////////////////////////////////////////////
class mainClass{
//List of variables
int terminator, threejetev, fourjetev, fivejetev, sixjetev, HTev, MHTev, delphicutn, finalHTev, finalMHTev, nolePev;
float xs, xserr;
double weight, CrossSection, CrossSectionError, totPx, totPy, HT, MHT, cutHT, cutMHT, pt, coss, sinn;
vector<vector<double> > vecjetvec, vecelecvec, vecmuvec;
vector<double> jetvec, elecvec, muvec;
vector<GenParticle*> GenParticlevec;
vector<TH1D > vec;
vector<Photon> photonvec;
char TreeList[200], tempname[200];
string pro, line;
fstream file, input;
map<string , vector<TH1D> > cut_histvec_map;  
map<string, map<string , vector<TH1D> > > map_map;
//constructor
public:
mainClass(string Pileup, string Process, string Detector, string Outdir, string inputnumber){
terminator=0;threejetev=0;fourjetev=0;fivejetev=0;sixjetev=0;HTev=0;MHTev=0;delphicutn=0;finalHTev=0;finalMHTev=0;nolePev=0;
weight=0; CrossSection=-999.0; CrossSectionError=0.0; totPx=0; totPy=0; HT=0; MHT=0; cutHT=0; cutMHT=0; pt=0; coss=0; sinn=0;
TChain chain("Delphes");
string cutname[]={"RA2nocut", "RA23Jetcut", "RA2HT500cut", "RA2MHT200cut", "RA2delphicut", "RA2noleptoncut", "RA24Jetcut", "RA25Jetcut", "RA26Jetcut", "RA2allbutHT2500cut", "RA2allbutMHT1000cut", "RA2allcut"};
//build a vector of histograms
TH1D  weight_hist = TH1D("weight", "Weight Distribution", 5,0,5);
vec.push_back(weight_hist);
TH1D  RA2HT_hist =  TH1D("HT","HT Distribution",50,0,5000);
vec.push_back(RA2HT_hist);
TH1D  RA2MHT_hist =  TH1D("MHT","MHT Distribution",100,0,5000);
vec.push_back(RA2MHT_hist);
TH1D  RA2NJet_hist = TH1D("NJet","Number of Jets Distribution",10,0,20);
vec.push_back(RA2NJet_hist);
//initialize a map between string=cutnames and histvecs. copy one histvec into all of them. The histograms, though, will be filled differently.
cut_histvec_map["RA2nocut"]=vec;cut_histvec_map["RA23Jetcut"]=vec;cut_histvec_map["RA2HT500cut"]=vec;cut_histvec_map["RA2MHT200cut"]=vec;cut_histvec_map["RA2delphicut"]=vec;
cut_histvec_map["RA2noleptoncut"]=vec;cut_histvec_map["RA24Jetcut"]=vec;cut_histvec_map["RA25Jetcut"]=vec;cut_histvec_map["RA26Jetcut"]=vec;
cut_histvec_map["RA2allbutHT2500cut"]=vec;cut_histvec_map["RA2allbutMHT1000cut"]=vec;cut_histvec_map["RA2allcut"]=vec;
//initialize a map between string and maps. copy the map of histvecs into each
map_map["allEvents"]=cut_histvec_map; map_map["Zvv"]=cut_histvec_map; map_map["Wlv"]=cut_histvec_map;

///Add the root files to a chain called Delphes
sprintf(TreeList,"./FileList/%s/%s_%s_%s",Detector.c_str(),Process.c_str(),Pileup.c_str(),inputnumber.c_str());
input.open(TreeList,std::fstream::in);
if(!input.is_open()){sprintf(TreeList,"./FileList/%s/%s_%s.list",Detector.c_str(),Process.c_str(),Pileup.c_str());input.open(TreeList,std::fstream::in);}
cout << "file name " << TreeList << endl; 
for(std::string linee; getline(input, linee);)
    {
      if (linee[0] == '#') continue;
      std::cout << "Add File: " << linee << std::endl;
      chain.Add(linee.c_str());
    }
if (chain.GetListOfFiles()->GetEntries() == 0)
  {
    std::cout << "No files attached! Exiting ...."  << std::endl;
    
  }
//end of adding file

//Get the cross-section from the file ./FileList/CrossSection.list
//open the file
file.open("FileList/CrossSection.list", std::fstream::in);
 if (!file.is_open())
  {
    std::cout << " Error to open the Cross Section file!" << std::endl;
    
  }
while (getline(file, line)){
if (line.empty()) continue;
    if (line[0] == '#') continue;

stringstream ss; ss << line;
ss >>  pro >> xs >> xserr;
if (pro==Process){CrossSection = xs; CrossSectionError = xserr ; break;}
}
   
  if (CrossSection == -999.)
  {
    std::cerr << "Unable to find the process and its cross section!" << std::endl;
    
  }
cout<<"\nCrossSection : "<<CrossSection<<" +- "<<CrossSectionError<<endl<<endl;
//end of acquiring XS

// Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
// Get pointers to branches used in this analysis
TClonesArray * branchEvent  = treeReader->UseBranch("Event");
TClonesArray * branchJet = treeReader->UseBranch("Jet");
TClonesArray * branchElectron = treeReader->UseBranch("Electron");
TClonesArray * branchMuon = treeReader->UseBranch("Muon");
TClonesArray * branchPhoton = treeReader->UseBranch("Photon");
TClonesArray * branchMet = treeReader->UseBranch("MissingET");
TClonesArray * branchHT = treeReader->UseBranch("ScalarHT");
TClonesArray * branchParticle = treeReader->UseBranch("Particle"); 
//report the total number of events
cout << "the total number of events: " << treeReader->GetEntries() << endl; 

//Loop Over all Events//Loop Over all Events//Loop Over all Events//Loop Over all Events//Loop Over all Events//Loop Over all Events//Loop Over all Events/
for(int entry = 0; entry < treeReader->GetEntries() ; entry++ ){
//for(int entry = 0; entry < 1000; entry++){
treeReader->ReadEntry(entry);

//Set Weight
LHEFEvent* event = (LHEFEvent*) branchEvent->At(0);
weight= event->Weight;

//a counter
if (entry % 5000 == 0){cout << "--------------------" << entry << endl;}

//MET and HT of the event
MissingET* met =(MissingET*) branchHT->At(0); 
//cout << "MET : " << met->MET << endl;
ScalarHT* ht= (ScalarHT*) branchHT->At(0);
//cout << "HT : " << ht->HT << endl;

GenParticlevec.clear();
///loop over all the particles in the history of an event. load them to a vector
for (int i = 0; i < branchParticle->GetEntries(); ++i)
  {
///define a pointer of class GenParticle and point it to the "entry"th event in the branch particle. Hence, we have access to PID, status and other properties of the particles in the event. 
GenParticle * particle = (GenParticle*)branchParticle->At(i);
GenParticlevec.push_back(particle);
}//end of loop over "particles in history" 

//////////////////loop over photons and load them to a vector
for (int i = 0; i < branchPhoton->GetEntries(); ++i)
  {
    Photon* pho = (Photon*)branchPhoton->At(i);
    if (fabs(pho->Eta) > 5)
      continue;
    photonvec.push_back(*pho);
  }
/////////////////end of loop over photons

////loop over electrons (load them to a vector)
elecvec.clear();
vecelecvec.clear();
for(int elecn=0; elecn <branchElectron->GetEntries();elecn++)
{

Electron* elec = (Electron*) branchElectron->At(elecn);
///for HT we want events with all elecs pt > 10 and |eta|< 2.5
if(elec->PT > 10 && elec->Eta < 2.5 && elec->Eta > (-2.5))
{
//the zeroth component is the tag of the elec/first:pt /second:phi/third:eta
elecvec.push_back((double) elecn);
elecvec.push_back(elec->PT);
elecvec.push_back((double)elec->Phi);
elecvec.push_back((double)elec->Eta);
vecelecvec.push_back(elecvec);

/// end of if over pt and eta for HT
} 
////end of loop over electrons
}

////loop over muons    (load them to a vector)
muvec.clear();
vecmuvec.clear();
for(int mun=0; mun <branchMuon->GetEntries();mun++)
{

Muon* mu = (Muon*) branchMuon->At(mun);
///for HT we want events with all muons pt > 10 and |eta|< 2.4
if(mu->PT > 10 && mu->Eta < 2.4 && mu->Eta > (-2.4))
{
//the zeroth component is the tag of the mu/first:pt /second:phi/third:eta
muvec.push_back((double) mun);
muvec.push_back(mu->PT);
muvec.push_back((double)mu->Phi);
muvec.push_back((double)mu->Eta);
vecmuvec.push_back(muvec);

/// end of if over pt and eta for HT
}
////end of loop over muons
}
///making the values zero for each event
totPx=0;
totPy=0;
HT=0;
///////////////////////////////////////////////////////////////////loop over jets    (load them to a vector)
vecjetvec.clear();
jetvec.clear();
for(int jetn=0; jetn <branchJet->GetEntries();jetn++){

Jet* jet = (Jet*) branchJet->At(jetn);
sinn = (double) sin(jet->Phi);
coss = (double) cos(jet->Phi);
pt = jet->PT;

///for HT we want events with all jets pt > 50 and |eta|< 2.5
if(pt>50 && jet->Eta < 2.5 && jet->Eta > (-2.5))
{
//the zeroth component is the tag of the jet/first:pt /second:phi/third:eta
jetvec.push_back((double) jetn);
jetvec.push_back(pt);
jetvec.push_back((double)jet->Phi);
jetvec.push_back((double)jet->Eta);
vecjetvec.push_back(jetvec);

///calculate HT
HT+=pt;
/// end of if over pt and eta for HT
}

//// for MHT we want events with all jets pt > 30 and |eta|< 5
if(pt>30 && jet->Eta < 5 && jet->Eta > (-5))
{
///calculate total jet-px and jet-py
totPx += pt * coss;  
totPy += pt * sinn;
///end of if over pt and eta for MHT
}
}///////////////////////////////////////////////////////////////////end of loop over jets

///find the three most energetic jets
while(terminator!=0){
terminator=0;
for(int iv=0; iv<vecjetvec.size()-1;iv++){

if(vecjetvec[iv][1]<vecjetvec[iv+1][1]){
swap(vecjetvec[iv],vecjetvec[iv+1]);
terminator+=1;
}
//end of the for
}
//end of the while
}
///end of find the three most energetic jets
///calculate MHT
MHT = sqrt( totPx*totPx + totPy*totPy );



//build and array that contains the quantities we need a histogram for. Here order is important and must be the same as RA2nocutvec
double eveinfvec[] = {weight, HT, MHT , vecjetvec.size()}; //the last one gives the RA2 defined number of jets.


//loop over all the different backgrounds: "allEvents", "Wlv", "Zvv"
for(map<string, map<string , vector<TH1D> > >::iterator itt=map_map.begin(); itt!=map_map.end();itt++){//this will be terminated after the cuts

//cout << "bg_type:  " << itt->first << ", bool:  " << bg_type(itt->first , GenParticlevec) << endl;
//determine what type of background should pass
if(bg_type(itt->first , GenParticlevec)==true){//all the cuts are inside this







//call the constructor of the histClass to fill the histograms before cuts are applied.
histClass nocut_Object( &eveinfvec[0] , &itt->second["RA2nocut"][0] );

//Cuts//Cuts//Cuts//Cuts//Cuts//Cuts//Cuts//Cuts//Cuts//Cuts//Cuts//Cuts//Cuts//Cuts//Cuts//Cuts//Cuts//Cuts//Cuts//Cuts

//there are jets missing in the event, which cause much larger MHT than expected. In order to supress this problem,
// we calculate the METMHTAsys by |MHT-MET|/(MHT+MET)
//if(fabs(MHT-met->MET)/(MHT+met->MET)<0.2){

/// jetcut
if(vecjetvec.size() >= 3 && vecjetvec[0][1]> 50 ){
threejetev+=1;
histClass Jet3cut_Object( &eveinfvec[0] , &itt->second["RA23Jetcut"][0] );
///HT cut
if(HT>=500){
HTev+=1;
histClass HT500cut_Object( &eveinfvec[0] , &itt->second["RA2HT500cut"][0] );
///MHT cut
if(MHT>=200){
MHTev+=1;
histClass MHT200cut_Object( &eveinfvec[0] , &itt->second["RA2MHT200cut"][0] );
///delta phi cut
if(delphi(vecjetvec[0],totPx,totPy,MHT)>0.5 && delphi(vecjetvec[1],totPx,totPy,MHT)>0.3 && delphi(vecjetvec[2],totPx,totPy,MHT)>0.3)
{
delphicutn+=1;
histClass delphicut_Object( &eveinfvec[0] , &itt->second["RA2delphicut"][0] );
/// electron veto
if(vecelecvec.size()==0)
{
///muon veto
if(vecmuvec.size()==0)
{
nolePev+=1;
histClass noleptoncut_Object( &eveinfvec[0] , &itt->second["RA2noleptoncut"][0] );

if(vecjetvec.size() >= 4){
fourjetev+=1;
histClass Jet4cut_Object( &eveinfvec[0] , &itt->second["RA24Jetcut"][0] );

if(vecjetvec.size() >= 5){
fivejetev+=1;
histClass Jet5cut_Object( &eveinfvec[0] , &itt->second["RA25Jetcut"][0] );

if(vecjetvec.size() >= 6){
sixjetev+=1;
histClass Jet6cut_Object( &eveinfvec[0] , &itt->second["RA26Jetcut"][0] );

if(MHT>=1000){ histClass allbutHT2500cut_Object( &eveinfvec[0] , &itt->second["RA2allbutHT2500cut"][0] );
}


if(HT>=2500){
finalHTev+=1;
histClass allbutMHT1000cut_Object(&eveinfvec[0] , &itt->second["RA2allbutMHT1000cut"][0]);
///Search region one MHT cut
if(MHT>=1000){
finalMHTev+=1;

histClass allcut_Object(&eveinfvec[0] , &itt->second["RA2allcut"][0]);

///Search region one MHT cut
}
///end of Search region one HT cut
}
///end of 6jetcut if
}
///end of 5jetcut if
}
///end of 4jetcut if
}
//end of muon if
}
//end of lepton if
}
///end of deltaphi if
}
//end of MHT if
}
//end of HT if
}
}///end of 3jetcut if
//}// end of |MHT-MET|/(MHT+MET)<0.2
//EndOfCuts//EndOfCuts//EndOfCuts//EndOfCuts//EndOfCuts//EndOfCuts//EndOfCuts//EndOfCuts//EndOfCuts//EndOfCuts

}//end of bg_type determination
}//end of loop over all the different backgrounds: "allEvents", "Wlv", "Zvv"

}///End of loop over all Events//end of loop over events//end of loop over events//end of loop over events//end of loop over events//end of loop over events//

//open a file to write the histograms
sprintf(tempname,"%s/results_%s_%s_%s_%s.root",Outdir.c_str(),Detector.c_str(),Process.c_str(),Pileup.c_str(),inputnumber.c_str());
TFile *resFile = new TFile(tempname, "RECREATE");
TDirectory *cdtoitt;
TDirectory *cdtoit;

for(map<string, map<string , vector<TH1D> > >::iterator itt=map_map.begin(); itt!=map_map.end();itt++){
cdtoitt = resFile->mkdir((itt->first).c_str());
cdtoitt->cd();
for(map<string , vector<TH1D> >::iterator it=itt->second.begin(); it!=itt->second.end();it++){ 
cdtoit =  cdtoitt->mkdir((it->first).c_str());
cdtoit->cd();
for(int i=0; i<=3; i++){//since we only have 4 type of histograms 
sprintf(tempname,"%s_%s_hist%d",(itt->first).c_str(),(it->first).c_str(),i);
it->second[i].Write(tempname);
}

cdtoitt->cd();
}
}
file.close();
input.close(); 




}//end of the constructor
};//end of class mainClass








int main()
{

mainClass mainObj("NoPileUp","T1qqqq_14TEV","PhaseI", "Results","00");
//mainClass mainObj_BJ("NoPileUp","BJ_14TEV_HT1","PhaseI", "Results","00");


return 0;
}

