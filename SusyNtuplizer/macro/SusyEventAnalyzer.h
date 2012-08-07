// -*- C++ -*-
//
// Package:    SusyNtuplizer
// Class:      SusyEventAnalyzer.h
// 
/*

 Description: an analyzer for susy::Event

 Implementation:

*/
//
// Original Author:  Dongwook Jang
// $Id: SusyEventAnalyzer.h,v 1.3 2011/10/27 13:16:01 dmorse Exp $
//

#ifndef SusyEventAnalyzer_h
#define SusyEventAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <map>

#include "../src/SusyEvent.h"

class SusyEventAnalyzer {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  susy::Event     *event;

  // List of branches
  TBranch        *b_Event;

  SusyEventAnalyzer(TTree *tree=0);
  virtual ~SusyEventAnalyzer();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();                          // event loop for main analysis
  virtual void     DR03();
  virtual void     Pileup();
  virtual void     Filter();
  virtual void     PhotonId();


  // utility functions
  bool isSameObject(TVector3& p1, TLorentzVector& p2, float dR_Cut);
  bool isSameObject(TLorentzVector& p1, TVector3& p2, float dR_Cut);
  bool isSameObject(TLorentzVector& p1, TLorentzVector& p2, float dR_Cut);
  bool isSameObject(TVector3& p1, TVector3& p2, float dR_Cut);
  float getDR(TVector3& p1, TLorentzVector& p2);
  float getDR(TLorentzVector& p1, TLorentzVector& p2);
  float getDphi(float p1, float p2);
  float d0correction(TVector3& beamSpot, susy::Track& track) const;
  void IncludeAJson(std::string jsonfile);  // Call to pull in a json file 
  bool isInJson(Int_t run,Int_t lumi);      // JSON based good run list cut...
  bool PassTrigger(TString v); // return true if path v is fired
  bool PassTriggers(); // return true if any of names in hltNames are fired
  float GetDiEmPt(susy::Photon* Pho1, susy::Photon* Pho2);//calculates and returns diEmPt
  float GetDiJetPt(susy::PFJet* Jet1, susy::PFJet* Jet2);//calculates and returns diJetPt
  float GetTriEmPt(susy::Photon* Pho1, susy::Photon* Pho2, susy::Photon* Pho3);//calculates and returns TriEmPt
  float InvariantMass(TLorentzVector P1, TLorentzVector P2);//calculates and returns Invariant Mass
  bool tooClosePhi(TVector3& p1, TVector3& p2, float phi_Cut);//checks whether gg or ff are too close to each other in phi
  bool tooClosePhi(TVector3& p1, TVector2& p2, float phi_Cut);//checks whether gg or ff are too close to MET in phi
  std::pair<float,float> GetMetReweight(float diEMPT,std::string type,std::vector< std::pair<float,float> > binEE[4],std::vector< std::pair<float,float> > binFF[4],std::vector< std::pair<float,float> > binEEsidebandLowJet[4],std::vector< std::pair<float,float> > binEEsidebandHighJet[4],int numJets);//gets reweighting for met plots from diempt ratio
  float GetRazrMr(susy::Photon* Pho1, susy::Photon* Pho2);
  float GetRazrR2(susy::Photon* Pho1, susy::Photon* Pho2, susy::MET* Met);
  void MatchPhosToJets(susy::Photon* pOne, susy::Photon* pTwo, std::vector<susy::PFJet*> jets, susy::PFJet* &jet1, susy::PFJet* &jet2, bool &hasdijetpt, float dR);
  float GetAlphaT(TLorentzVector pOne,TLorentzVector pTwo);
  float GetPhotonLessHt(float Ht, TLorentzVector pOne, TLorentzVector pTwo);
  // parameter configuration functions
  void Initialize();         // global variables needed to be initialized just once
  void InitializePerEvent(); // global variables needed to be initialized per event
  void SetDataset(TString& v) {          ds = v; }
  void SetPrintInterval(int v) {         printInterval = v; }
  void SetPrintLevel(int v) {            printLevel = v; }
  void SetProcessNEvents(int v) {        processNEvents = v; }
  void SetUseTrigger(bool v) {           useTrigger = v; }
  void SetUseJSON(bool v) {              useJSON = v; }
  void AddHltName(TString v) {           hltNames.push_back(v); }
  void SetFilter(bool v) {               enableFilter = v; }
  void SetFilteredFileName(TString v) {  filtered_file_name = v; }
  void SetOutputEventNumbers(bool v) {   outputEventNumbers = v; }
  void DoRhoCorrection(bool v) {         doRhoCorrection = v; }
  void SetDR03Rho25Corr(float ecal, float hcal, float track){ PUCorr_ECAL=ecal;PUCorr_HCAL=hcal;PUCorr_TRACK = track; }
  void SetPFisoRho25Corr(float ch, float nh, float ph){ PUCorr_chargedHadron=ch;PUCorr_neutralHadron=nh;PUCorr_photon=ph; }

 private:

  TString ds;               // dataset name to be used for output histfile name

  // printLevel
  // 0 : default - no printout
  // 1 : print functional step in every event
  // 2 : print values in collections
  int printLevel;           // print frequency

  int printInterval;        // print level for event content: defined in Event.h
  int processNEvents;       // number of events to be processed
  bool useTrigger;          // flag for using trigger bit selection.
  bool useJSON;             // flag for using JSON selection
  std::vector<TString> hltNames;          // HLT trigger path names
  bool enableFilter;        // filter events of interest
  TString filtered_file_name; // filtered output file name
  bool outputEventNumbers;  // print run and event numbers for gg, eg, ee, ff to txt 
  bool doRhoCorrection;
  bool doNVertexCorrection;
  float PUCorr_ECAL;
  float PUCorr_HCAL;
  float PUCorr_TRACK;
  float PUCorr_chargedHadron;
  float PUCorr_neutralHadron;
  float PUCorr_photon;

  typedef std::map<int,std::map<int,bool> > RunLumiFlagHolder;  //define map that holds json list
  RunLumiFlagHolder goodrunlumilist;  // instantiate it

};

#endif

#ifdef SusyEventAnalyzer_cxx
SusyEventAnalyzer::SusyEventAnalyzer(TTree *tree)
{
  if (tree == 0) {
    std::cout << "Error!!! There is no file containing a tree." << std::endl;
  }
  Init(tree);
  Initialize();
}

SusyEventAnalyzer::~SusyEventAnalyzer()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t SusyEventAnalyzer::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry,0);
}
Long64_t SusyEventAnalyzer::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
  }
  return centry;
}

void SusyEventAnalyzer::Init(TTree *tree)
{
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  //fChain->SetMakeClass(1);

  event = new susy::Event;

  fChain->SetBranchAddress("susyEvent", &event, &b_Event);
  fChain->SetBranchStatus("l1Map*",0);
  fChain->SetBranchStatus("superClusters*",0);
  fChain->SetBranchStatus("clusters*",0);
  //fChain->SetBranchStatus("electrons*",0);
  //fChain->SetBranchStatus("muons*",0);
  fChain->SetBranchStatus("caloJets*",0);
  fChain->SetBranchStatus("jptJets*",0);
  fChain->SetBranchStatus("generalTracks*",0);
}

void SusyEventAnalyzer::Initialize() {

  ds = "test";
  printLevel = 0;
  printInterval = 1000;
  processNEvents = -1;
  useTrigger = false;
  useJSON = false;
  enableFilter = false;
  filtered_file_name = "filtered.root";

}

void SusyEventAnalyzer::IncludeAJson(std::string jsonfile) {


// Fairly primitive brute force json parser -- opens the json file named in the argument
// and adds that to the goodrunlumilist map.  Overlapping jsons are merged inclusively.

 char thing;

 ifstream jsonInput;

 std::cout << "Sucking in Json file: " << jsonfile << " which includes: " << std::endl;

 jsonInput.open(jsonfile.c_str());
 
 if (!jsonInput.good()) {
   std::cout << "Problem reading Json file...  Didn't suck anything in... " << std::endl;
   return;
 }
 
 jsonInput.get(thing);
 
 while (jsonInput.good()) {
   if (thing=='{') {  // start of list
     while (thing != '}') {
       int runnum;
       if (thing == '"') {
         std::string srunnum;
         jsonInput.get(thing); // get stuff inside ""

         while (thing != '"') {
           srunnum+=thing; // get stuff inside ""
           jsonInput.get(thing);

	   }
         sscanf(srunnum.c_str(),"%i",&runnum);
         std::cout << " runnum: " << runnum << std::endl;
         bool newrun=true;
         
       } // inside ""
       if (thing == '[') {
          jsonInput.get(thing); // get stuff inside []
	 while (thing != ']') {
           if (thing == '[') {
             jsonInput.get(thing); // get stuff inside series []

             std::string lumiseries;
             int firstlumi,lastlumi;
             while (thing !=']') {
               lumiseries+=thing;
                jsonInput.get(thing); // get stuff inside series []
             }
             sscanf(lumiseries.c_str(),"%i,%i",&firstlumi,&lastlumi);
             std::cout << "  lumis  " << firstlumi << " to " << lastlumi << std::endl;

	     // At this point have runnum, first lumi, last lumi -- so can fill map here...
	     for (int l=firstlumi;l<=lastlumi;l++) {
               goodrunlumilist[runnum][l]=true;
	     }

           } // inside actual series []
             jsonInput.get(thing); // get stuff inside []
         }
       } // inside []
         jsonInput.get(thing); // get another char looking for "

     } 
   } // inside {}
    jsonInput.get(thing); // get another char looking for {

 } // EOF 

 jsonInput.close();

}


bool SusyEventAnalyzer::isInJson(Int_t run,Int_t lumi) {

//#ifdef MC
//  return 1;
//#endif

if (goodrunlumilist[run][lumi]) return true;

return false;

}


bool SusyEventAnalyzer::isSameObject(TVector3& p1, TVector3& p2, float dR_Cut) {

  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
  if(dR < dR_Cut) return true;
  //if(dR < 0.5) return true;
  return false;
}

bool SusyEventAnalyzer::isSameObject(TVector3& p1, TLorentzVector& p2, float dR_Cut) {
  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
  if(dR < dR_Cut) return true;
  //if(dR < 0.5) return true;
  return false;
}
bool SusyEventAnalyzer::isSameObject(TLorentzVector& p1, TVector3& p2, float dR_Cut) {
  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
  if(dR < dR_Cut) return true;
  //if(dR < 0.5) return true;
  return false;
}
bool SusyEventAnalyzer::isSameObject(TLorentzVector& p1, TLorentzVector& p2, float dR_Cut) {
  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
  if(dR < dR_Cut) return true;
  //if(dR < 0.5) return true;
  return false;
}


float SusyEventAnalyzer::getDR(TVector3& p1, TLorentzVector& p2) {

  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
  return dR;
}
float SusyEventAnalyzer::getDR(TLorentzVector& p1, TLorentzVector& p2) {

  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
  return dR;
}

float SusyEventAnalyzer::getDphi(float p1, float p2)
{
  float dPhi = std::fabs(TVector2::Phi_mpi_pi(p1 - p2));
  return dPhi;
}


float SusyEventAnalyzer::d0correction(TVector3& beamSpot, susy::Track& track) const {
  float d0 = track.d0() - beamSpot.X()*std::sin(track.phi()) + beamSpot.Y()*std::cos(track.phi());
  return d0;
}

bool SusyEventAnalyzer::PassTrigger(TString path) {
  bool pass = false;
  for(susy::TriggerMap::iterator it = event->hltMap.begin(); it != event->hltMap.end(); it++) {
    if(it->first.Contains(path) && (int(it->second.second)) && it->second.first==1 ) {
      //std::cout<<"Path "<<it->first<<" passed!  Prescale: "<<it->second.first<<endl;
      pass = true;
      break;
    }
  }
  return pass;
}


bool SusyEventAnalyzer::PassTriggers() {
  bool pass = false;
  for(std::vector<TString>::iterator it = hltNames.begin(); it != hltNames.end(); it++) {
    if(PassTrigger(*it)) {
      pass = true;
      break;
    }
  }
  return pass;
}

bool SusyEventAnalyzer::tooClosePhi(TVector3& p1, TVector3& p2, float phi_Cut) {
  float dPhi = std::fabs(TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi()));
  if(dPhi<phi_Cut){
    //std::cout<<"p1Phi= "<<p1.Phi()<<"  p2Phi= "<<p2.Phi()<<"  dPhi= "<<dPhi<<endl; 
    return true;
  }
  return false;
}

bool SusyEventAnalyzer::tooClosePhi(TVector3& p1, TVector2& p2, float phi_Cut) {
  float dPhi = std::fabs(TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi()));
  if(dPhi<phi_Cut){
    //std::cout<<"p1Phi= "<<p1.Phi()<<"  p2Phi= "<<p2.Phi()<<"  dPhi= "<<dPhi<<endl; 
    return true;
  }
  return false;
}

float SusyEventAnalyzer::InvariantMass(TLorentzVector P1, TLorentzVector P2){
  //TLorentzVector P4=Pho1->momentum + Pho2->momentum;
  float InvarMass = (P1 + P2).M();
  return InvarMass;
}

float SusyEventAnalyzer::GetDiEmPt(susy::Photon* Pho1, susy::Photon* Pho2){
  /*float PtX = Pho1->momentum.Pt()*cos(Pho1->caloPosition.Phi()) + Pho2->momentum.Pt()*cos(Pho2->caloPosition.Phi());
    float PtY = Pho1->momentum.Pt()*sin(Pho1->caloPosition.Phi()) + Pho2->momentum.Pt()*sin(Pho2->caloPosition.Phi());
    float diEmPt = sqrt(PtX*PtX + PtY*PtY);*/
  float diEmPt = (Pho1->momentum + Pho2->momentum).Pt();
  /*cout<<"PhoOnePT: "<<Pho1->momentum.Pt()<<endl;
    cout<<"PhoTwoPT: "<<Pho2->momentum.Pt()<<endl;
    cout<<"PtX: "<<PtX<<endl;
    cout<<"PtY: "<<PtY<<endl;
    cout<<"diEmPt: "<<diEmPt<<endl;*/
  return diEmPt;
}

float SusyEventAnalyzer::GetDiJetPt(susy::PFJet* Jet1, susy::PFJet* Jet2){
  /*float PtX = Jet1->momentum.Pt()*cos(Jet1->momentum.Phi()) + Jet2->momentum.Pt()*cos(Jet2->momentum.Phi());
    float PtY = Jet1->momentum.Pt()*sin(Jet1->momentum.Phi()) + Jet2->momentum.Pt()*sin(Jet2->momentum.Phi());
    float diJetPt = sqrt(PtX*PtX + PtY*PtY);*/
  float diJetPt = (Jet1->momentum + Jet2->momentum).Pt();
  /*cout<<"PhoOnePT: "<<Pho1->momentum.Pt()<<endl;
    cout<<"PhoTwoPT: "<<Pho2->momentum.Pt()<<endl;
    cout<<"PtX: "<<PtX<<endl;
    cout<<"PtY: "<<PtY<<endl;
    cout<<"diEmPt: "<<diEmPt<<endl;*/
  return diJetPt;
}

float SusyEventAnalyzer::GetTriEmPt(susy::Photon* Pho1, susy::Photon* Pho2,susy::Photon* Pho3){
  /*float PtX = Pho1->momentum.Pt()*cos(Pho1->caloPosition.Phi()) + Pho2->momentum.Pt()*cos(Pho2->caloPosition.Phi()) + Pho3->momentum.Pt()*cos(Pho3->caloPosition.Phi());
    float PtY = Pho1->momentum.Pt()*sin(Pho1->caloPosition.Phi()) + Pho2->momentum.Pt()*sin(Pho2->caloPosition.Phi()) + Pho3->momentum.Pt()*sin(Pho3->caloPosition.Phi());*/
  float TriEmPt = (Pho1->momentum + Pho2->momentum + Pho3->momentum).Pt();
  /*cout<<"PhoOnePT: "<<Pho1->momentum.Pt()<<endl;
    cout<<"PhoTwoPT: "<<Pho2->momentum.Pt()<<endl;
    cout<<"PtX: "<<PtX<<endl;
    cout<<"PtY: "<<PtY<<endl;
    cout<<"diEmPt: "<<diEmPt<<endl;*/
  return TriEmPt;
}

void SusyEventAnalyzer::InitializePerEvent() {

}


float SusyEventAnalyzer::GetRazrMr(susy::Photon* Pho1, susy::Photon* Pho2){
  float E1 = Pho1->momentum.E();
  float E2 = Pho2->momentum.E();
  float pz1 = Pho1->momentum.Pz();
  float pz2 = Pho2->momentum.Pz();
  float MR = sqrt((E1+E2)*(E1+E2) - (pz1+pz2)*(pz1+pz2));
  return MR;
}
float SusyEventAnalyzer::GetRazrR2(susy::Photon* Pho1, susy::Photon* Pho2, susy::MET* Met){
  float Mr = GetRazrMr(Pho1,Pho2);
  float met = Met->met();
  float pt1 = Pho1->momentum.Pt();
  float pt2 = Pho2->momentum.Pt();
  float Pt1plusPt2 = pt1+pt2;
  TLorentzVector p = Pho1->momentum + Pho2->momentum;
  float MxdotPx = Met->metX()*p.Px();
  float MydotPy = Met->metY()*p.Py();
  float MdotP = MxdotPx + MydotPy;
  float MT = sqrt((met*Pt1plusPt2 - MdotP)/2);
  float R = MT/Mr;
  float R2 = R*R;
  return R2;
}


float SusyEventAnalyzer::GetAlphaT(TLorentzVector pOne, TLorentzVector pTwo)
{
  float minEt = (( pOne.Et()<pTwo.Et())?pOne.Et():pTwo.Et() );
  float Et2 = (pOne.Et()+pTwo.Et())*(pOne.Et()+pTwo.Et());
  float Px2 = (pOne.Px()+pTwo.Px())*(pOne.Px()+pTwo.Px());
  float Py2 = (pOne.Py()+pTwo.Py())*(pOne.Py()+pTwo.Py());
  float Mt = sqrt(Et2 - Px2 - Py2);
  float alphaT = minEt/Mt;
  return alphaT;
}

float SusyEventAnalyzer::GetPhotonLessHt(float Ht, TLorentzVector pOne, TLorentzVector pTwo)
{
  float pLessHt = Ht - pOne.Et() - pTwo.Et();
  return pLessHt;
}

template<typename T> bool EtGreater(const T* p1, const T* p2) {
  return (p1->momentum.Et() > p2->momentum.Et());
}


#endif // #ifdef SusyEventAnalyzer_cxx
