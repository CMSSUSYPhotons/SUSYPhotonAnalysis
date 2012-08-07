//
// Package:    SusyNtuplizer
// Class:      SusyEventAnalyzer.cc
// 
//Description: an analyzer for susy::Event

#define SusyEventAnalyzer_cxx

#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TGraphAsymmErrors.h>
#include <TRandom.h>
#include <TLegend.h>

#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <utility>
#include <fstream>
#include <vector>

#include "SusyEventAnalyzer.h"
#include "SusyEventPrinter.h"

#include "../jec/JetMETObjects/interface/JetCorrectorParameters.h"
#include "../jec/JetMETObjects/interface/FactorizedJetCorrector.h"

using namespace std;


int numBins=15;
Double_t DiEMPtBins[16]={0.,5.,10.,15.,20.,30.,40.,50.,60.,70.,80.,100.,120.,150.,200.,600.,};

/* //old way
   pair<float,float> SusyEventAnalyzer::GetMetReweight(float diEMPT,string type,vector< pair<float,float> > binEE[4],vector< pair<float,float> > binFF[4],int numJets){
   if(printLevel>5)cout<<"inside GetMetReweight"<<endl;
   float reweight=1., error=0.;//,lowff=0.,highff=0.,lowee=0.,highee=0.;
   int bin=999;
   if(printLevel>5)cout<<"Line 167"<<endl;
   if(numJets>2)numJets=2;
   int numJetsBin=numJets+1;
   if(diEMPT<92) bin=int(diEMPT/4.);//----4 GeV bins so divide by that to get bin number (0-100)
   else if(diEMPT>=92 && diEMPT<600) bin=int((diEMPT+92)/8.);//----8 GeV bins so divide by that to get bin number (100-600)
   //else if(diEMPT>=200) bin=int((diEMPT+1200)/20.);//----20 GeV bins so divide by 20 to get bin number (200-600)
   if( type=="ff" ){
   reweight = binFF[numJetsBin][bin].first;
   error = binFF[numJetsBin][bin].second;  
   }
   if( type=="ee" ){
   if(printLevel>5)cout<<"Line 177"<<endl;
   reweight = binEE[numJetsBin][bin].first;
   if(printLevel>5)cout<<"Line 179"<<endl;
   error = binEE[numJetsBin][bin].second;
   if(printLevel>5)cout<<"Line 181"<<endl;
   }
   if(bin==999 || reweight==0) reweight=1.;
   //cout<<"reweight: "<<reweight<<endl;
   if(printLevel>5)cout<<"End of GetMetReweight"<<endl;
   return make_pair(reweight,error);
   }

   pair<float,float> SusyEventAnalyzer::GetMetReweightJet(float diJETPT,string typeJet,vector< pair<float,float> > binEEjet[4],vector< pair<float,float> > binFFjet[4],int numJets){
   float reweight=1., error=0.;//,lowff=0.,highff=0.,lowee=0.,highee=0.;
   int bin=999;
   if(numJets>2)numJets=2;
   int numJetsBin=numJets+1;
   if(diJETPT<92) bin=int(diJETPT/4.);//----4 GeV bins so divide by that to get bin number (0-100)
   else if(diJETPT>=92 && diJETPT<600) bin=int((diJETPT+92)/8.);//----8 GeV bins so divide by that to get bin number (100-600)
   //else if(diJETPT>=200) bin=int((diJETPT+1200)/20.);//----20 GeV bins so divide by 20 to get bin number (200-600)
   if( typeJet=="ff" ){
   reweight = binFFjet[numJetsBin][bin].first;
   error = binFFjet[numJetsBin][bin].second;  
   }
   if( typeJet=="ee" ){
   reweight = binEEjet[numJetsBin][bin].first;
   error = binEEjet[numJetsBin][bin].second;
   }
   if(bin==999 || reweight==0) reweight=1.;
   //cout<<"reweight: "<<reweight<<endl;
   return make_pair(reweight,error);
   }
*/
pair<float,float> SusyEventAnalyzer::GetMetReweight(float diEMPT,string type,vector< pair<float,float> > binEE[4],vector< pair<float,float> > binFF[4],vector< pair<float,float> > binEEsidebandLowJet[4],vector< pair<float,float> > binEEsidebandHighJet[4],int numJets){
  //if(printLevel>5)cout<<"inside GetMetReweight"<<endl;
  float reweight=1., error=0.;//,lowff=0.,highff=0.,lowee=0.,highee=0.;
  int bin=999;
  //if(printLevel>5)cout<<"Line 167"<<endl;
  if(numJets>2)numJets=2;
  int numJetsBin=numJets+1;
  //this gives the bins from 0-(numBins-1) even though they really go from 1-numBins, so that it lines up with the vector numbering
  if(diEMPT<20) bin=int(diEMPT/5.);//----5 GeV bins so divide by that to get bin number (0-20)
  else if(diEMPT>=20 && diEMPT<80) bin=int((diEMPT+20)/10.);//----10 GeV bins so divide by that to get bin number (20-80)
  else if(diEMPT>=80 && diEMPT<120) bin=int((diEMPT+120)/20.);//----20 GeV bins so divide by 20 to get bin number (80-120)
  else if(diEMPT>=120 && diEMPT<150) bin=int((diEMPT+240)/30.);
  else if(diEMPT>=150 && diEMPT<200) bin=int((diEMPT+500)/50.);
  else if(diEMPT>=200) bin=14;
  //This was an attempt to be smart - it is incredibly slow!!! don't try to be smart...
  /*TH1F* binHist = new TH1F("binHist","",numBins,DiEMPtBins);
    bin=binHist->FindBin(diEMPT)-1;*/
  if( type=="ff" ){
    reweight = binFF[numJetsBin][bin].first;
    error = binFF[numJetsBin][bin].second;  
  }
  else if( type=="ee" ){
    reweight = binEE[numJetsBin][bin].first;
    error = binEE[numJetsBin][bin].second;
  }
  else if( type=="eeSbLow" ){
    reweight = binEEsidebandLowJet[numJetsBin][bin].first;
    error = binEEsidebandLowJet[numJetsBin][bin].second;
  }
  else if( type=="eeSbHigh" ){
    reweight = binEEsidebandHighJet[numJetsBin][bin].first;
    error = binEEsidebandHighJet[numJetsBin][bin].second;
  }
  else cout<<"Trying to reweight but have no type ee,ff,eeSbLow,eeSbHigh"<<endl;
  if(bin==-1 || bin==999 || reweight==0) {reweight=1.;error=0.;}
  //cout<<"reweight: "<<reweight<<endl;
  //if(printLevel>5)cout<<"End of GetMetReweight"<<endl;
  //delete binHist;
  return make_pair(reweight,error);
}

void SusyEventAnalyzer::MatchPhosToJets(susy::Photon* pOne, susy::Photon* pTwo, std::vector<susy::PFJet*> jets, susy::PFJet* &jet1, susy::PFJet* &jet2, bool &hasdijetpt, float dR)
{
  hasdijetpt=false;
  //cout<<"hasdijetpt="<<hasdijetpt<<endl;
  for(std::vector<susy::PFJet*>::iterator jet_it1 = jets.begin(); jet_it1 != jets.end(); jet_it1++){	
    if(isSameObject(pOne->caloPosition,(*jet_it1)->momentum,dR)){
      for(std::vector<susy::PFJet*>::iterator jet_it2 = jets.begin(); jet_it2 != jets.end(); jet_it2++){
	if(/*jet_it2!=jet_it1*/ !isSameObject((*jet_it1)->momentum,(*jet_it2)->momentum,0.1) && isSameObject(pTwo->caloPosition,(*jet_it2)->momentum,dR)){
	  jet1=*jet_it1;jet2=*jet_it2;
	  //cout<<"Inside function Jet pT: "<<jet1->momentum.Pt()<<" , "<<jet2->momentum.Pt()<<endl;
	  hasdijetpt=true;break;
	}//end jet2 match
      }//end jet2 iterator
      if(hasdijetpt)break;
    }//end jet1 match
  }//end jet1 iterator
  //cout<<"hasdijetpt="<<hasdijetpt<<endl;
  //return;
}

void SusyEventAnalyzer::Loop() {

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  cout << "total events in files  : " << nentries << endl;

  if(processNEvents <= 0 || processNEvents > nentries) processNEvents = nentries;

  cout << "events to be processed : " << processNEvents << endl; 
 
  if(printLevel > 1) cout << "Initialize event counters." << endl;
  const int NCNT = 20;
  int nCnt[NCNT];
  for(int i=0; i<NCNT; i++) nCnt[i] = 0;

  int nFiltered = 0;
  TTree* filterTree = 0;
  
  if(printLevel > 1) cout << "Define Filter File" << endl;
  if(enableFilter) {
    TFile* filterFile = new TFile("/data/ndpc3/b/dmorse/RA3/SusyNtuples/FilteredNtuples/"+ds+"-SelectedEvents.root","RECREATE");
    //TFile* filterFile = new TFile("/Users/dmorse/RA3/SusyNtuples/FilteredNtuples/"+ds+"-SelectedEvents.root","RECREATE");
    filterTree = (TTree*) fChain->GetTree()->CloneTree(0);
    filterTree->SetAutoSave();
  }
  
  if(printLevel > 1) cout<<"Open hist file" << endl;
  //TFile* fout = new TFile("/tmp/dmorse/hist_"+ds+".root","RECREATE");
  TFile* fout = new TFile("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_"+ds+".root","RECREATE");
  //TFile* fout = new TFile("/Users/dmorse/RA3/AnalysisOutput/hist_"+ds+".root","RECREATE");
  if(printLevel > 1) cout<<"Hist file opened" << endl;

  fout->cd();

  if(printLevel > 1) cout<<"Define histograms" << endl;
  TH1F* h_vtxZ = new TH1F("vtxZ","Z position of the primary vertex;Z (cm);Events",100,-50.,50.);
  TH1F* h_bsZ = new TH1F("bsZ","Z position of the beam spot;Z (cm);Events",20,-10.,10.);
  TH1F* h_met = new TH1F("met","missing transverse energy;#slash{E}_{T} (GeV);Events",220,0.0,1100.);
  TH1F* h_sumEt = new TH1F("sumEt","Scalar sum of all calorimeter energy;#sigmaE_{T} (GeV);Events",500,0.0,5000.);
  TH1F* h_sumEt_gg = new TH1F("sumEt_gg","Scalar sum of all calorimeter energy in gg events;#sigmaE_{T} (GeV);Events",500,0.0,5000.);
  TH1F* h_sumEt_eg = new TH1F("sumEt_eg","Scalar sum of all calorimeter energy in eg events;#sigmaE_{T} (GeV);Events",500,0.0,5000.);
  TH1F* h_sumEt_ee = new TH1F("sumEt_ee","Scalar sum of all calorimeter energy in ee events with 81<InvariantMass<101GeV;#sigmaE_{T} (GeV);Events",500,0.0,5000.);
  TH1F* h_sumEt_eeFullRange = new TH1F("sumEt_eeFullRange","Scalar sum of all calorimeter energy in ee events with full InvarianMass range;#sigmaE_{T} (GeV);Events",500,0.0,5000.);
  TH1F* h_sumEt_ff = new TH1F("sumEt_ff","Scalar sum of all calorimeter energy in ff events;#sigmaE_{T} (GeV);Events",500,0.0,5000.);

  TH1F* h_PhoPt = new TH1F("PhoPt","All Photon Transverse Momentum (GeV)", 200,0.,1000.);
  TH1F* h_ggPt = new TH1F("ggPt","gg Transverse Momentum;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_egPt = new TH1F("egPt","eg Transverse Momentum;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_eePt = new TH1F("eePt","ee Transverse Momentum;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_ffPt = new TH1F("ffPt","ff Transverse Momentum;p_{T} (GeV);Events", 200,0.,1000.);

  TH1F* h_ggdPhi = new TH1F("ggdPhi","#Delta#phi between gg objects;#Delta#phi;Events",32,0.,3.2);
  TH1F* h_egdPhi = new TH1F("egdPhi","#Delta#phi between eg objects;#Delta#phi;Events",32,0.,3.2);
  TH1F* h_eedPhi = new TH1F("eedPhi","#Delta#phi between ee objects;#Delta#phi;Events",32,0.,3.2);
  TH1F* h_ffdPhi = new TH1F("ffdPhi","#Delta#phi between ff objects;#Delta#phi;Events",32,0.,3.2);

  TH1F* h_ggMETdPhiLead = new TH1F("ggMETdPhiLead","#Delta#phi between Lead gg object and M_{E_{T}};#Delta#phi(pho,MET);Events",32,0.,3.2);
  TH1F* h_egMETdPhiLead = new TH1F("egMETdPhiLead","#Delta#phi between Lead eg object and M_{E_{T}};#Delta#phi(pho,MET);Events",32,0.,3.2);
  TH1F* h_eeMETdPhiLead = new TH1F("eeMETdPhiLead","#Delta#phi between Lead ee object and M_{E_{T}};#Delta#phi(pho,MET);Events",32,0.,3.2);
  TH1F* h_ffMETdPhiLead = new TH1F("ffMETdPhiLead","#Delta#phi between Lead ff object and M_{E_{T}};#Delta#phi(pho,MET);Events",32,0.,3.2);

  TH1F* h_ggMETdPhiTrail = new TH1F("ggMETdPhiTrail","#Delta#phi between Trail gg object and M_{E_{T}};#Delta#phi(pho,MET);Events",32,0.,3.2);
  TH1F* h_egMETdPhiTrail = new TH1F("egMETdPhiTrail","#Delta#phi between Trail eg object and M_{E_{T}};#Delta#phi(pho,MET);Events",32,0.,3.2);
  TH1F* h_eeMETdPhiTrail = new TH1F("eeMETdPhiTrail","#Delta#phi between Trail ee object and M_{E_{T}};#Delta#phi(pho,MET);Events",32,0.,3.2);
  TH1F* h_ffMETdPhiTrail = new TH1F("ffMETdPhiTrail","#Delta#phi between Trail ff object and M_{E_{T}};#Delta#phi(pho,MET);Events",32,0.,3.2);

  TH1F* h_ggAlphaT = new TH1F("ggAlphaT","#alpha_{T} of gg sample;#alpha_{T};Events",200,0.,2.); 
  TH1F* h_egAlphaT = new TH1F("egAlphaT","#alpha_{T} of eg sample;#alpha_{T};Events",200,0.,2.);
  TH1F* h_eeAlphaT = new TH1F("eeAlphaT","#alpha_{T} of ee sample;#alpha_{T};Events",200,0.,2.);
  TH1F* h_ffAlphaT = new TH1F("ffAlphaT","#alpha_{T} of ff sample;#alpha_{T};Events",200,0.,2.);

  TH1F* h_ggAlphaT_0Jet = new TH1F("ggAlphaT_0Jet","#alpha_{T} of gg sample in 0Jet case;#alpha_{T};Events",200,0.,2.); 
  TH1F* h_egAlphaT_0Jet = new TH1F("egAlphaT_0Jet","#alpha_{T} of eg sample in 0Jet case;#alpha_{T};Events",200,0.,2.);
  TH1F* h_eeAlphaT_0Jet = new TH1F("eeAlphaT_0Jet","#alpha_{T} of ee sample in 0Jet case;#alpha_{T};Events",200,0.,2.);
  TH1F* h_ffAlphaT_0Jet = new TH1F("ffAlphaT_0Jet","#alpha_{T} of ff sample in 0Jet case;#alpha_{T};Events",200,0.,2.);

  TH1F* h_ggPhotonLessHt = new TH1F("ggPhotonLessHt","Scalar sum of all calorimeter energy minus two EM objects;#SigmaE_{T} (GeV);Events",500,0.0,5000.);
  TH1F* h_egPhotonLessHt = new TH1F("egPhotonLessHt","Scalar sum of all calorimeter energy minus two EM objects;#SigmaE_{T} (GeV);Events",500,0.0,5000.);
  TH1F* h_eePhotonLessHt = new TH1F("eePhotonLessHt","Scalar sum of all calorimeter energy minus two EM objects;#SigmaE_{T} (GeV);Events",500,0.0,5000.);
  TH1F* h_ffPhotonLessHt = new TH1F("ffPhotonLessHt","Scalar sum of all calorimeter energy minus two EM objects;#SigmaE_{T} (GeV);Events",500,0.0,5000.);

  TH2F* h_ggPhotonLessHtVsMET = new TH2F("ggPhotonLessHtVsMET","Scalar sum of all calorimeter energy minus two EM objects VS MET;MET;#SigmaE_{T} (GeV)",200,0.,1000.,500,0.0,5000.);
  TH2F* h_egPhotonLessHtVsMET = new TH2F("egPhotonLessHtVsMET","Scalar sum of all calorimeter energy minus two EM objects VS MET;MET;#SigmaE_{T} (GeV)",200,0.,1000.,500,0.0,5000.);
  TH2F* h_eePhotonLessHtVsMET = new TH2F("eePhotonLessHtVsMET","Scalar sum of all calorimeter energy minus two EM objects VS MET;MET;#SigmaE_{T} (GeV)",200,0.,1000.,500,0.0,5000.);
  TH2F* h_ffPhotonLessHtVsMET = new TH2F("ffPhotonLessHtVsMET","Scalar sum of all calorimeter energy minus two EM objects VS MET;MET;#SigmaE_{T} (GeV)",200,0.,1000.,500,0.0,5000.);

  TH1F* h_ggePt = new TH1F("ggePt","gge Transverse Momentum;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_ggeElePt = new TH1F("ggeElePt","gge Electron Transverse Momentum;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_ggeDiEMPt = new TH1F("ggeDiEMPt","(gg)e DiEMPt;DiEMPt (GeV);Events", 120,0.,600.);h_ggeDiEMPt->Sumw2();
  TH1F* h_ggeTriEMPt = new TH1F("ggeTriEMPt","gge TriEMPt;TriEMPt (GeV);Events", 120,0.,600.);h_ggeTriEMPt->Sumw2();
  TH1F* h_ggeInvarMass = new TH1F("ggeInvarMass","(gg)e Invariant Mass;(GeV);Events", 1000,0.,1000.);
  TH1F* h_ggeInvarMassMET30 = new TH1F("ggeInvarMassMET30","(gg)e Invariant Mass in Events with MET>30;(GeV);Events", 200,0.,200.);
  TH1F* h_ggeSigIetaIeta = new TH1F("ggeSigIetaIeta","gge SigmaIetaIeta;#sigma_{i#etai#eta};Events", 20,0.,.02);
  TH1F* h_ggeMet = new TH1F("ggeMet","Missing transverse energy in events with two photons with Et>40,25 GeV and one electron;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_rho_gge = new TH1F("rho_gge","Fastjet correction rho for all photonWevents in gge sample",200,0,40);
  TH1F* h_NVertex_gge = new TH1F("NVertex_gge","Number of good vertices for all events in gge sample",40,0,40);
  
  TH1F* h_ggeePt = new TH1F("ggeePt","ggee Transverse Momentum;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_ggeeElePt = new TH1F("ggeeElePt","ggee Electron Transverse Momentum;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_ggeeDiEMPt = new TH1F("ggeeDiEMPt","(gg)e DiEMPt;DiEMPt (GeV);Events", 120,0.,600.);h_ggeeDiEMPt->Sumw2();
  TH1F* h_ggeeInvarMass = new TH1F("ggeeInvarMass","gg Invariant Mass in ggee events;(GeV);Events", 1000,0.,1000.);
  TH1F* h_ggeeSigIetaIeta = new TH1F("ggeeSigIetaIeta","ggee SigmaIetaIeta;#sigma_{i#etai#eta};Events", 20,0.,.02);
  TH1F* h_ggeeMet = new TH1F("ggeeMet","Missing transverse energy in events with two photons with Et>40,25 GeV and one electron;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_rho_ggee = new TH1F("rho_ggee","Fastjet correction rho for all photonWevents in ggee sample",200,0,40);
  TH1F* h_NVertex_ggee = new TH1F("NVertex_ggee","Number of good vertices for all events in ggee sample",40,0,40);

  TH1F* h_ggmMet = new TH1F("ggmMet","Missing transverse energy in events with two photons with Et>40,25 GeV and one pfMuon;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_ggmInvarMass = new TH1F("ggmInvarMass","(gg)m Invariant Mass;(GeV);Events", 1000,0.,1000.);
  TH1F* h_ggmInvarMassMET30 = new TH1F("ggmInvarMassMET30","(gg)m Invariant Mass in Events with MET>30;(GeV);Events", 200,0.,200.);

  TH1F* h_eegPt = new TH1F("eegPt","eeg Transverse Momentum;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_eegElePt = new TH1F("eegElePt","eeg Electron Transverse Momentum;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_eegDiEMPt = new TH1F("eegDiEMPt","(gg)e DiEMPt;DiEMPt (GeV);Events", 120,0.,600.);h_eegDiEMPt->Sumw2();
  TH1F* h_eegTriEMPt = new TH1F("eegTriEMPt","eeg TriEMPt;TriEMPt (GeV);Events", 120,0.,600.);h_eegTriEMPt->Sumw2();
  TH1F* h_eegInvarMass = new TH1F("eegInvarMass","(gg)e Invariant Mass;(GeV);Events", 1000,0.,1000.);
  TH1F* h_eegSigIetaIeta = new TH1F("eegSigIetaIeta","eeg SigmaIetaIeta;#sigma_{i#etai#eta};Events", 20,0.,.02);
  TH1F* h_eegMet = new TH1F("eegMet","Missing transverse energy in events with two photons with Et>40,25 GeV and one electron;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_rho_eeg = new TH1F("rho_eeg","Fastjet correction rho for all photonWevents in eeg sample",200,0,40);
  TH1F* h_NVertex_eeg = new TH1F("NVertex_eeg","Number of good vertices for all events in eeg sample",40,0,40);

  TH1F* h_eeggPt = new TH1F("eeggPt","eegg Transverse Momentum;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_eeggElePt = new TH1F("eeggElePt","eegg Electron Transverse Momentum;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_eeggDiEMPt = new TH1F("eeggDiEMPt","(gg)e DiEMPt;DiEMPt (GeV);Events", 120,0.,600.);h_eeggDiEMPt->Sumw2();
  TH1F* h_eeggInvarMass = new TH1F("eeggInvarMass","(gg)e Invariant Mass;(GeV);Events", 1000,0.,1000.);
  TH1F* h_eeggSigIetaIeta = new TH1F("eeggSigIetaIeta","eegg SigmaIetaIeta;#sigma_{i#etai#eta};Events", 20,0.,.02);
  TH1F* h_eeggMet = new TH1F("eeggMet","Missing transverse energy in events with two photons with Et>40,25 GeV and one electron;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_rho_eegg = new TH1F("rho_eegg","Fastjet correction rho for all photonWevents in eegg sample",200,0,40);
  TH1F* h_NVertex_eegg = new TH1F("NVertex_eegg","Number of good vertices for all events in eegg sample",40,0,40);

  TH1F* h_ggDiEMPt = new TH1F("ggDiEMPt","gg DiEMPt;DiEMPt (GeV);Events", 120,0.,600.);h_ggDiEMPt->Sumw2();
  TH1F* h_egDiEMPt = new TH1F("egDiEMPt","eg DiEMPt;DiEMPt (GeV);Events", 120,0.,600.);h_egDiEMPt->Sumw2();
  TH1F* h_eeDiEMPt = new TH1F("eeDiEMPt","ee DiEMPt;DiEMPt (GeV);Events", 120,0.,600.);h_eeDiEMPt->Sumw2();
  TH1F* h_ffDiEMPt = new TH1F("ffDiEMPt","ff DiEMPt;DiEMPt (GeV);Events", 120,0.,600.);h_ffDiEMPt->Sumw2();

  TH1F* h_eeDiEMPt_reweight_binned = new TH1F("eeDiEMPt_reweight_binned","ee DiEMPt;DiEMPt (GeV);Events", numBins,DiEMPtBins);h_eeDiEMPt_reweight_binned->Sumw2();
  TH1F* h_eeSideBandDiEMPt = new TH1F("eeSideBandDiEMPt","eeSideBand DiEMPt;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandDiEMPt->Sumw2();
  TH1F* h_eeSideBandDiEMPt_0Jet = new TH1F("eeSideBandDiEMPt_0Jet","eeSideBand_0Jet DiEMPt;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandDiEMPt_0Jet->Sumw2();
  TH1F* h_eeSideBandDiEMPt_1Jet = new TH1F("eeSideBandDiEMPt_1Jet","eeSideBand_1Jet DiEMPt;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandDiEMPt_1Jet->Sumw2();
  TH1F* h_eeSideBandDiEMPt_2Jet = new TH1F("eeSideBandDiEMPt_2Jet","eeSideBand_2Jet DiEMPt;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandDiEMPt_2Jet->Sumw2();
  TH1F* h_eeSideBandLowDiEMPt = new TH1F("eeSideBandLowDiEMPt","eeSideBandLow DiEMPt - 71<InvMass<81;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandLowDiEMPt->Sumw2();
  TH1F* h_eeSideBandLowDiEMPt_0Jet = new TH1F("eeSideBandLowDiEMPt_0Jet","eeSideBandLow_0Jet DiEMPt - 71<InvMass<81;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandLowDiEMPt_0Jet->Sumw2();
  TH1F* h_eeSideBandLowDiEMPt_1Jet = new TH1F("eeSideBandLowDiEMPt_1Jet","eeSideBandLow_1Jet DiEMPt - 71<InvMass<81;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandLowDiEMPt_1Jet->Sumw2();
  TH1F* h_eeSideBandLowDiEMPt_2Jet = new TH1F("eeSideBandLowDiEMPt_2Jet","eeSideBandLow_2Jet DiEMPt - 71<InvMass<81;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandLowDiEMPt_2Jet->Sumw2();
  TH1F* h_eeSideBandHighDiEMPt = new TH1F("eeSideBandHighDiEMPt","eeSideBandHigh DiEMPt - 101<InvMass<111;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandHighDiEMPt->Sumw2();
  TH1F* h_eeSideBandHighDiEMPt_0Jet = new TH1F("eeSideBandHighDiEMPt_0Jet","eeSideBandHigh_0Jet DiEMPt - 101<InvMass<111;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandHighDiEMPt_0Jet->Sumw2();
  TH1F* h_eeSideBandHighDiEMPt_1Jet = new TH1F("eeSideBandHighDiEMPt_1Jet","eeSideBandHigh_1Jet DiEMPt - 101<InvMass<111;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandHighDiEMPt_1Jet->Sumw2();
  TH1F* h_eeSideBandHighDiEMPt_2Jet = new TH1F("eeSideBandHighDiEMPt_2Jet","eeSideBandHigh_2Jet DiEMPt - 101<InvMass<111;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandHighDiEMPt_2Jet->Sumw2();
  TH1F* h_eeSideBandDiJetPt = new TH1F("eeSideBandDiJetPt","eeSideBand DiJetPt;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandDiJetPt->Sumw2();
  TH1F* h_eeSideBandDiJetPt_0Jet = new TH1F("eeSideBandDiJetPt_0Jet","eeSideBand_0Jet DiJetPt;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandDiJetPt_0Jet->Sumw2();
  TH1F* h_eeSideBandDiJetPt_1Jet = new TH1F("eeSideBandDiJetPt_1Jet","eeSideBand_1Jet DiJetPt;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandDiJetPt_1Jet->Sumw2();
  TH1F* h_eeSideBandDiJetPt_2Jet = new TH1F("eeSideBandDiJetPt_2Jet","eeSideBand_2Jet DiJetPt;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandDiJetPt_2Jet->Sumw2();
  TH1F* h_eeSideBandLowDiJetPt = new TH1F("eeSideBandLowDiJetPt","eeSideBandLow DiJetPt - 71<InvMass<81;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandLowDiJetPt->Sumw2();
  TH1F* h_eeSideBandLowDiJetPt_0Jet = new TH1F("eeSideBandLowDiJetPt_0Jet","eeSideBandLow_0Jet DiJetPt - 71<InvMass<81;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandLowDiJetPt_0Jet->Sumw2();
  TH1F* h_eeSideBandLowDiJetPt_1Jet = new TH1F("eeSideBandLowDiJetPt_1Jet","eeSideBandLow_1Jet DiJetPt - 71<InvMass<81;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandLowDiJetPt_1Jet->Sumw2();
  TH1F* h_eeSideBandLowDiJetPt_2Jet = new TH1F("eeSideBandLowDiJetPt_2Jet","eeSideBandLow_2Jet DiJetPt - 71<InvMass<81;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandLowDiJetPt_2Jet->Sumw2();
  TH1F* h_eeSideBandHighDiJetPt = new TH1F("eeSideBandHighDiJetPt","eeSideBandHigh DiJetPt - 101<InvMass<111;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandHighDiJetPt->Sumw2();
  TH1F* h_eeSideBandHighDiJetPt_0Jet = new TH1F("eeSideBandHighDiJetPt_0Jet","eeSideBandHigh_0Jet DiJetPt - 101<InvMass<111;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandHighDiJetPt_0Jet->Sumw2();
  TH1F* h_eeSideBandHighDiJetPt_1Jet = new TH1F("eeSideBandHighDiJetPt_1Jet","eeSideBandHigh_1Jet DiJetPt - 101<InvMass<111;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandHighDiJetPt_1Jet->Sumw2();
  TH1F* h_eeSideBandHighDiJetPt_2Jet = new TH1F("eeSideBandHighDiJetPt_2Jet","eeSideBandHigh_2Jet DiJetPt - 101<InvMass<111;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandHighDiJetPt_2Jet->Sumw2();

  TH1F* h_eeDiEMPt_reweight_binned_JetReq = new TH1F("eeDiEMPt_reweight_binned_JetReq","ee_JetReq DiEMPt;DiEMPt (GeV);Events", numBins,DiEMPtBins);h_eeDiEMPt_reweight_binned_JetReq->Sumw2();
  TH1F* h_eeSideBandLowDiEMPt_JetReq = new TH1F("eeSideBandLowDiEMPt_JetReq","eeSideBandLow_JetReq DiEMPt - 71<InvMass<81;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandLowDiEMPt_JetReq->Sumw2();
  TH1F* h_eeSideBandLowDiEMPt_JetReq_0Jet = new TH1F("eeSideBandLowDiEMPt_JetReq_0Jet","eeSideBandLow_JetReq_0Jet DiEMPt - 71<InvMass<81;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandLowDiEMPt_JetReq_0Jet->Sumw2();
  TH1F* h_eeSideBandLowDiEMPt_JetReq_1Jet = new TH1F("eeSideBandLowDiEMPt_JetReq_1Jet","eeSideBandLow_JetReq_1Jet DiEMPt - 71<InvMass<81;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandLowDiEMPt_JetReq_1Jet->Sumw2();
  TH1F* h_eeSideBandLowDiEMPt_JetReq_2Jet = new TH1F("eeSideBandLowDiEMPt_JetReq_2Jet","eeSideBandLow_JetReq_2Jet DiEMPt - 71<InvMass<81;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandLowDiEMPt_JetReq_2Jet->Sumw2();
  TH1F* h_eeSideBandHighDiEMPt_JetReq = new TH1F("eeSideBandHighDiEMPt_JetReq","eeSideBandHigh_JetReq DiEMPt - 101<InvMass<111;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandHighDiEMPt_JetReq->Sumw2();
  TH1F* h_eeSideBandHighDiEMPt_JetReq_0Jet = new TH1F("eeSideBandHighDiEMPt_JetReq_0Jet","eeSideBandHigh_JetReq_0Jet DiEMPt - 101<InvMass<111;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandHighDiEMPt_JetReq_0Jet->Sumw2();
  TH1F* h_eeSideBandHighDiEMPt_JetReq_1Jet = new TH1F("eeSideBandHighDiEMPt_JetReq_1Jet","eeSideBandHigh_JetReq_1Jet DiEMPt - 101<InvMass<111;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandHighDiEMPt_JetReq_1Jet->Sumw2();
  TH1F* h_eeSideBandHighDiEMPt_JetReq_2Jet = new TH1F("eeSideBandHighDiEMPt_JetReq_2Jet","eeSideBandHigh_JetReq_2Jet DiEMPt - 101<InvMass<111;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandHighDiEMPt_JetReq_2Jet->Sumw2();
  TH1F* h_eeSideBandDiJetPt_JetReq = new TH1F("eeSideBandDiJetPt_JetReq","eeSideBand_JetReq DiJetPt;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandDiJetPt_JetReq->Sumw2();
  TH1F* h_eeSideBandDiJetPt_JetReq_0Jet = new TH1F("eeSideBandDiJetPt_JetReq_0Jet","eeSideBand_JetReq_0Jet DiJetPt;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandDiJetPt_JetReq_0Jet->Sumw2();
  TH1F* h_eeSideBandDiJetPt_JetReq_1Jet = new TH1F("eeSideBandDiJetPt_JetReq_1Jet","eeSideBand_JetReq_1Jet DiJetPt;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandDiJetPt_JetReq_1Jet->Sumw2();
  TH1F* h_eeSideBandDiJetPt_JetReq_2Jet = new TH1F("eeSideBandDiJetPt_JetReq_2Jet","eeSideBand_JetReq_2Jet DiJetPt;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandDiJetPt_JetReq_2Jet->Sumw2();
  TH1F* h_eeSideBandLowDiJetPt_JetReq = new TH1F("eeSideBandLowDiJetPt_JetReq","eeSideBandLow_JetReq DiJetPt - 71<InvMass<81;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandLowDiJetPt_JetReq->Sumw2();
  TH1F* h_eeSideBandLowDiJetPt_JetReq_0Jet = new TH1F("eeSideBandLowDiJetPt_JetReq_0Jet","eeSideBandLow_JetReq_0Jet DiJetPt - 71<InvMass<81;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandLowDiJetPt_JetReq_0Jet->Sumw2();
  TH1F* h_eeSideBandLowDiJetPt_JetReq_1Jet = new TH1F("eeSideBandLowDiJetPt_JetReq_1Jet","eeSideBandLow_JetReq_1Jet DiJetPt - 71<InvMass<81;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandLowDiJetPt_JetReq_1Jet->Sumw2();
  TH1F* h_eeSideBandLowDiJetPt_JetReq_2Jet = new TH1F("eeSideBandLowDiJetPt_JetReq_2Jet","eeSideBandLow_JetReq_2Jet DiJetPt - 71<InvMass<81;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandLowDiJetPt_JetReq_2Jet->Sumw2();
  TH1F* h_eeSideBandHighDiJetPt_JetReq = new TH1F("eeSideBandHighDiJetPt_JetReq","eeSideBandHigh_JetReq DiJetPt - 101<InvMass<111;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandHighDiJetPt_JetReq->Sumw2();
  TH1F* h_eeSideBandHighDiJetPt_JetReq_0Jet = new TH1F("eeSideBandHighDiJetPt_JetReq_0Jet","eeSideBandHigh_JetReq_0Jet DiJetPt - 101<InvMass<111;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandHighDiJetPt_JetReq_0Jet->Sumw2();
  TH1F* h_eeSideBandHighDiJetPt_JetReq_1Jet = new TH1F("eeSideBandHighDiJetPt_JetReq_1Jet","eeSideBandHigh_JetReq_1Jet DiJetPt - 101<InvMass<111;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandHighDiJetPt_JetReq_1Jet->Sumw2();
  TH1F* h_eeSideBandHighDiJetPt_JetReq_2Jet = new TH1F("eeSideBandHighDiJetPt_JetReq_2Jet","eeSideBandHigh_JetReq_2Jet DiJetPt - 101<InvMass<111;DiJetPt (GeV);Events", 120,0.,600.);h_eeSideBandHighDiJetPt_JetReq_2Jet->Sumw2();

  TH1F* h_ggDiJetPt = new TH1F("ggDiJetPt","gg DiJetPt;DiJetPt (GeV);Events", 120,0.,600.);h_ggDiJetPt->Sumw2();
  TH1F* h_eeDiJetPt = new TH1F("eeDiJetPt","ee DiJetPt;DiJetPt (GeV);Events", 120,0.,600.);h_eeDiJetPt->Sumw2();
  TH1F* h_ggDiJetPt_JetReq = new TH1F("ggDiJetPt_JetReq","gg_JetReq DiJetPt;DiJetPt (GeV);Events", 120,0.,600.);h_ggDiJetPt_JetReq->Sumw2();
  TH1F* h_eeDiJetPt_JetReq = new TH1F("eeDiJetPt_JetReq","ee_JetReq DiJetPt;DiJetPt (GeV);Events", 120,0.,600.);h_eeDiJetPt_JetReq->Sumw2();
  TH1F* h_ggDiJetPt_JetFail = new TH1F("ggDiJetPt_JetFail","gg DiJetPt in events that have a 50GeV pfJet failing loose id;DiJetPt (GeV);Events", 120,0.,600.);h_ggDiJetPt_JetFail->Sumw2();
  TH1F* h_egDiJetPt_JetFail = new TH1F("egDiJetPt_JetFail","eg DiJetPt in events that have a 50GeV pfJet failing loose id;DiJetPt (GeV);Events", 120,0.,600.);h_egDiJetPt_JetFail->Sumw2();
  TH1F* h_eeDiJetPt_JetFail = new TH1F("eeDiJetPt_JetFail","ee DiJetPt in events that have a 50GeV pfJet failing loose id;DiJetPt (GeV);Events", 120,0.,600.);h_eeDiJetPt_JetFail->Sumw2();
  TH1F* h_ffDiJetPt_JetFail = new TH1F("ffDiJetPt_JetFail","ff DiJetPt in events that have a 50GeV pfJet failing loose id;DiJetPt (GeV);Events", 120,0.,600.);h_ffDiJetPt_JetFail->Sumw2();
  TH1F* h_ffDiJetPt = new TH1F("ffDiJetPt","ff DiJetPt;DiJetPt (GeV);Events", 120,0.,600.);h_ffDiJetPt->Sumw2();
  TH1F* h_ffDiJetPt_JetReq = new TH1F("ffDiJetPt_JetReq","ff_JetReq DiJetPt;DiJetPt (GeV);Events", 120,0.,600.);h_ffDiJetPt_JetReq->Sumw2();
  TH1F* h_eeDiJetPt_reweight_binned = new TH1F("eeDiJetPt_reweight_binned","ee DiJetPt;DiJetPt (GeV);Events", 120,0.,600.);h_eeDiJetPt_reweight_binned->Sumw2();
  TH1F* h_ffDiJetPt_reweight_binned= new TH1F("ffDiJetPt_reweight_binned","ff DiJetPt;DiJetPt (GeV);Events", 120,0.,600.);h_ffDiJetPt_reweight_binned->Sumw2();
  TH1F* h_eeDiJetPt_reweight_binned_JetReq = new TH1F("eeDiJetPt_reweight_binned_JetReq","ee_JetReq DiJetPt;DiJetPt (GeV);Events", 120,0.,600.);h_eeDiJetPt_reweight_binned_JetReq->Sumw2();
  TH1F* h_ffDiJetPt_reweight_binned_JetReq = new TH1F("ffDiJetPt_reweight_binned_JetReq","ff_JetReq DiJetPt;DiJetPt (GeV);Events", 120,0.,600.);h_ffDiJetPt_reweight_binned_JetReq->Sumw2();

  TH2F* h_ggDiJetPtOverDiEMPtVsDiEMPt = new TH2F("h_ggDiJetPtOverDiEMPtVsDiEMPt","DiJetPt/DiEMPt vs DiEMPt for gg sample",120,0.,600.,100,0,100);h_ggDiJetPtOverDiEMPtVsDiEMPt->Sumw2();
  TH2F* h_eeDiJetPtOverDiEMPtVsDiEMPt = new TH2F("h_eeDiJetPtOverDiEMPtVsDiEMPt","DiJetPt/DiEMPt vs DiEMPt for ee sample",120,0.,600.,100,0,100);h_eeDiJetPtOverDiEMPtVsDiEMPt->Sumw2();
  TH2F* h_ffDiJetPtOverDiEMPtVsDiEMPt = new TH2F("h_ffDiJetPtOverDiEMPtVsDiEMPt","DiJetPt/DiEMPt vs DiEMPt for ff sample",120,0.,600.,100,0,100);h_ffDiJetPtOverDiEMPtVsDiEMPt->Sumw2();
  TH2F* h_ggDiJetPtOverDiEMPtVsDiEMPt_JetReq = new TH2F("h_ggDiJetPtOverDiEMPtVsDiEMPt_JetReq","DiJetPt/DiEMPt vs DiEMPt for gg_JetReq sample",120,0.,600.,100,0,100);h_ggDiJetPtOverDiEMPtVsDiEMPt_JetReq->Sumw2();
  TH2F* h_eeDiJetPtOverDiEMPtVsDiEMPt_JetReq = new TH2F("h_eeDiJetPtOverDiEMPtVsDiEMPt_JetReq","DiJetPt/DiEMPt vs DiEMPt for ee_JetReq sample",120,0.,600.,100,0,100);h_eeDiJetPtOverDiEMPtVsDiEMPt_JetReq->Sumw2();
  TH2F* h_ffDiJetPtOverDiEMPtVsDiEMPt_JetReq = new TH2F("h_ffDiJetPtOverDiEMPtVsDiEMPt_JetReq","DiJetPt/DiEMPt vs DiEMPt for ff_JetReq sample",120,0.,600.,100,0,100);h_ffDiJetPtOverDiEMPtVsDiEMPt_JetReq->Sumw2();

  TH1F* h_ggDiEMPtInvMassBelow40 = new TH1F("ggDiEMPtInvMassBelow40","gg DiEMPt in events with InvariantMass<40GeV;DiEMPt (GeV);Events", 120,0.,600.);h_ggDiEMPtInvMassBelow40->Sumw2();
  TH1F* h_egDiEMPtInvMassBelow40 = new TH1F("egDiEMPtInvMassBelow40","eg DiEMPt in events with InvariantMass<40GeV;DiEMPt (GeV);Events", 120,0.,600.);h_egDiEMPtInvMassBelow40->Sumw2();
  TH1F* h_eeDiEMPtInvMassBelow40 = new TH1F("eeDiEMPtInvMassBelow40","ee DiEMPt in events with InvariantMass<40GeV;DiEMPt (GeV);Events", 120,0.,600.);h_eeDiEMPtInvMassBelow40->Sumw2();
  TH1F* h_ffDiEMPtInvMassBelow40 = new TH1F("ffDiEMPtInvMassBelow40","ff DiEMPt in events with InvariantMass<40GeV;DiEMPt (GeV);Events", 120,0.,600.);h_ffDiEMPtInvMassBelow40->Sumw2();
  TH1F* h_ggDiEMPtInvMassAbove40 = new TH1F("ggDiEMPtInvMassAbove40","gg DiEMPt in events with InvariantMass>40GeV;DiEMPt (GeV);Events", 120,0.,600.);h_ggDiEMPtInvMassAbove40->Sumw2();
  TH1F* h_egDiEMPtInvMassAbove40 = new TH1F("egDiEMPtInvMassAbove40","eg DiEMPt in events with InvariantMass>40GeV;DiEMPt (GeV);Events", 120,0.,600.);h_egDiEMPtInvMassAbove40->Sumw2();
  TH1F* h_eeDiEMPtInvMassAbove40 = new TH1F("eeDiEMPtInvMassAbove40","ee DiEMPt in events with InvariantMass>40GeV;DiEMPt (GeV);Events", 120,0.,600.);h_eeDiEMPtInvMassAbove40->Sumw2();
  TH1F* h_ffDiEMPtInvMassAbove40 = new TH1F("ffDiEMPtInvMassAbove40","ff DiEMPt in events with InvariantMass>40GeV;DiEMPt (GeV);Events", 120,0.,600.);h_ffDiEMPtInvMassAbove40->Sumw2();

  TH1F* h_ffDiEMPt_reweight_binned = new TH1F("ffDiEMPt_reweight_binned","ff DiEMPt;DiEMPt (GeV);Events", numBins,DiEMPtBins);h_ffDiEMPt_reweight_binned->Sumw2();
  TH1F* h_ffDiEMPt_reweight_binned_JetReq = new TH1F("ffDiEMPt_reweight_binned_JetReq","ff_JetReq DiEMPt;DiEMPt (GeV);Events", numBins,DiEMPtBins);h_ffDiEMPt_reweight_binned_JetReq->Sumw2();

  TH1F* h_ggDiEMPt_JetReq = new TH1F("ggDiEMPt_JetReq","gg DiEMPt with One Jet Requirement;DiEMPt (GeV);Events", 120,0.,600.);h_ggDiEMPt_JetReq->Sumw2();
  TH1F* h_egDiEMPt_JetReq = new TH1F("egDiEMPt_JetReq","eg DiEMPt with One Jet Requirement;DiEMPt (GeV);Events", 120,0.,600.);h_egDiEMPt_JetReq->Sumw2();
  TH1F* h_eeDiEMPt_JetReq = new TH1F("eeDiEMPt_JetReq","ee DiEMPt with One Jet Requirement;DiEMPt (GeV);Events", 120,0.,600.);h_eeDiEMPt_JetReq->Sumw2();
  TH1F* h_ffDiEMPt_JetReq = new TH1F("ffDiEMPt_JetReq","ff DiEMPt with One Jet Requirement;DiEMPt (GeV);Events", 120,0.,600.);h_ffDiEMPt_JetReq->Sumw2();
  
  TH1F* h_ggDiEMPt_0Jet = new TH1F("ggDiEMPt_0Jet","gg DiEMPt with Exactly Zero Jets;DiEMPt (GeV);Events", 120,0.,600.);h_ggDiEMPt_0Jet->Sumw2();
  TH1F* h_egDiEMPt_0Jet = new TH1F("egDiEMPt_0Jet","eg DiEMPt with Exactly Zero Jets;DiEMPt (GeV);Events", 120,0.,600.);h_egDiEMPt_0Jet->Sumw2();
  TH1F* h_eeDiEMPt_0Jet = new TH1F("eeDiEMPt_0Jet","ee DiEMPt with Exactly Zero Jets;DiEMPt (GeV);Events", 120,0.,600.);h_eeDiEMPt_0Jet->Sumw2();
  TH1F* h_ffDiEMPt_0Jet = new TH1F("ffDiEMPt_0Jet","ff DiEMPt with Exactly Zero Jets;DiEMPt (GeV);Events", 120,0.,600.);h_ffDiEMPt_0Jet->Sumw2();

  TH1F* h_ggDiEMPt_1Jet = new TH1F("ggDiEMPt_1Jet","gg DiEMPt with Exactly One Jet;DiEMPt (GeV);Events", 120,0.,600.);h_ggDiEMPt_1Jet->Sumw2();
  TH1F* h_egDiEMPt_1Jet = new TH1F("egDiEMPt_1Jet","eg DiEMPt with Exactly One Jet;DiEMPt (GeV);Events", 120,0.,600.);h_egDiEMPt_1Jet->Sumw2();
  TH1F* h_eeDiEMPt_1Jet = new TH1F("eeDiEMPt_1Jet","ee DiEMPt with Exactly One Jet;DiEMPt (GeV);Events", 120,0.,600.);h_eeDiEMPt_1Jet->Sumw2();
  TH1F* h_ffDiEMPt_1Jet = new TH1F("ffDiEMPt_1Jet","ff DiEMPt with Exactly One Jet;DiEMPt (GeV);Events", 120,0.,600.);h_ffDiEMPt_1Jet->Sumw2();

  TH1F* h_ggDiEMPt_2Jet = new TH1F("ggDiEMPt_2Jet","gg DiEMPt with at least Two Jets;DiEMPt (GeV);Events", 120,0.,600.);h_ggDiEMPt_2Jet->Sumw2();
  TH1F* h_egDiEMPt_2Jet = new TH1F("egDiEMPt_2Jet","eg DiEMPt with at least Two Jets;DiEMPt (GeV);Events", 120,0.,600.);h_egDiEMPt_2Jet->Sumw2();
  TH1F* h_eeDiEMPt_2Jet = new TH1F("eeDiEMPt_2Jet","ee DiEMPt with at least Two Jets;DiEMPt (GeV);Events", 120,0.,600.);h_eeDiEMPt_2Jet->Sumw2();
  TH1F* h_ffDiEMPt_2Jet = new TH1F("ffDiEMPt_2Jet","ff DiEMPt with at least Two Jets;DiEMPt (GeV);Events", 120,0.,600.);h_ffDiEMPt_2Jet->Sumw2();

  TH1F* h_ggDiEMPt_MoreThan2Jets = new TH1F("ggDiEMPt_MoreThan2Jets","gg DiEMPt with Greater Than 2 Jets;DiEMPt (GeV);Events", 120,0.,600.);h_ggDiEMPt_MoreThan2Jets->Sumw2();
  TH1F* h_egDiEMPt_MoreThan2Jets = new TH1F("egDiEMPt_MoreThan2Jets","eg DiEMPt with Greater Than 2 Jets;DiEMPt (GeV);Events", 120,0.,600.);h_egDiEMPt_MoreThan2Jets->Sumw2();
  TH1F* h_eeDiEMPt_MoreThan2Jets = new TH1F("eeDiEMPt_MoreThan2Jets","ee DiEMPt with Greater Than 2 Jets;DiEMPt (GeV);Events", 120,0.,600.);h_eeDiEMPt_MoreThan2Jets->Sumw2();
  TH1F* h_ffDiEMPt_MoreThan2Jets = new TH1F("ffDiEMPt_MoreThan2Jets","ff DiEMPt with Greater Than 2 Jets;DiEMPt (GeV);Events", 120,0.,600.);h_ffDiEMPt_MoreThan2Jets->Sumw2();

  TH1F* h_ggDiJetPt_0Jet = new TH1F("ggDiJetPt_0Jet","gg DiJetPt with Exactly Zero Jets;DiJetPt (GeV);Events", 120,0.,600.);h_ggDiJetPt_0Jet->Sumw2();
  TH1F* h_egDiJetPt_0Jet = new TH1F("egDiJetPt_0Jet","eg DiJetPt with Exactly Zero Jets;DiJetPt (GeV);Events", 120,0.,600.);h_egDiJetPt_0Jet->Sumw2();
  TH1F* h_eeDiJetPt_0Jet = new TH1F("eeDiJetPt_0Jet","ee DiJetPt with Exactly Zero Jets;DiJetPt (GeV);Events", 120,0.,600.);h_eeDiJetPt_0Jet->Sumw2();
  TH1F* h_ffDiJetPt_0Jet = new TH1F("ffDiJetPt_0Jet","ff DiJetPt with Exactly Zero Jets;DiJetPt (GeV);Events", 120,0.,600.);h_ffDiJetPt_0Jet->Sumw2();

  TH1F* h_ggDiJetPt_1Jet = new TH1F("ggDiJetPt_1Jet","gg DiJetPt with Exactly One Jet;DiJetPt (GeV);Events", 120,0.,600.);h_ggDiJetPt_1Jet->Sumw2();
  TH1F* h_egDiJetPt_1Jet = new TH1F("egDiJetPt_1Jet","eg DiJetPt with Exactly One Jet;DiJetPt (GeV);Events", 120,0.,600.);h_egDiJetPt_1Jet->Sumw2();
  TH1F* h_eeDiJetPt_1Jet = new TH1F("eeDiJetPt_1Jet","ee DiJetPt with Exactly One Jet;DiJetPt (GeV);Events", 120,0.,600.);h_eeDiJetPt_1Jet->Sumw2();
  TH1F* h_ffDiJetPt_1Jet = new TH1F("ffDiJetPt_1Jet","ff DiJetPt with Exactly One Jet;DiJetPt (GeV);Events", 120,0.,600.);h_ffDiJetPt_1Jet->Sumw2();

  TH1F* h_ggDiJetPt_2Jet = new TH1F("ggDiJetPt_2Jet","gg DiJetPt with at least Two Jets;DiJetPt (GeV);Events", 120,0.,600.);h_ggDiJetPt_2Jet->Sumw2();
  TH1F* h_egDiJetPt_2Jet = new TH1F("egDiJetPt_2Jet","eg DiJetPt with at least Two Jets;DiJetPt (GeV);Events", 120,0.,600.);h_egDiJetPt_2Jet->Sumw2();
  TH1F* h_eeDiJetPt_2Jet = new TH1F("eeDiJetPt_2Jet","ee DiJetPt with at least Two Jets;DiJetPt (GeV);Events", 120,0.,600.);h_eeDiJetPt_2Jet->Sumw2();
  TH1F* h_ffDiJetPt_2Jet = new TH1F("ffDiJetPt_2Jet","ff DiJetPt with at least Two Jets;DiJetPt (GeV);Events", 120,0.,600.);h_ffDiJetPt_2Jet->Sumw2();

  TH1F* h_ggDiJetPt_MoreThan2Jets = new TH1F("ggDiJetPt_MoreThan2Jets","gg DiJetPt with Greater Than 2 Jets;DiJetPt (GeV);Events", 120,0.,600.);h_ggDiJetPt_MoreThan2Jets->Sumw2();
  TH1F* h_egDiJetPt_MoreThan2Jets = new TH1F("egDiJetPt_MoreThan2Jets","eg DiJetPt with Greater Than 2 Jets;DiJetPt (GeV);Events", 120,0.,600.);h_egDiJetPt_MoreThan2Jets->Sumw2();
  TH1F* h_eeDiJetPt_MoreThan2Jets = new TH1F("eeDiJetPt_MoreThan2Jets","ee DiJetPt with Greater Than 2 Jets;DiJetPt (GeV);Events", 120,0.,600.);h_eeDiJetPt_MoreThan2Jets->Sumw2();
  TH1F* h_ffDiJetPt_MoreThan2Jets = new TH1F("ffDiJetPt_MoreThan2Jets","ff DiJetPt with Greater Than 2 Jets;DiJetPt (GeV);Events", 120,0.,600.);h_ffDiJetPt_MoreThan2Jets->Sumw2();

  TH1F* h_ggDiEMPt_JetReq_0Jet = new TH1F("ggDiEMPt_JetReq_0Jet","gg_JetReq DiEMPt with Exactly Zero Jets;DiEMPt (GeV);Events", 120,0.,600.);h_ggDiEMPt_JetReq_0Jet->Sumw2();
  TH1F* h_egDiEMPt_JetReq_0Jet = new TH1F("egDiEMPt_JetReq_0Jet","eg_JetReq DiEMPt with Exactly Zero Jets;DiEMPt (GeV);Events", 120,0.,600.);h_egDiEMPt_JetReq_0Jet->Sumw2();
  TH1F* h_eeDiEMPt_JetReq_0Jet = new TH1F("eeDiEMPt_JetReq_0Jet","ee_JetReq DiEMPt with Exactly Zero Jets;DiEMPt (GeV);Events", 120,0.,600.);h_eeDiEMPt_JetReq_0Jet->Sumw2();
  TH1F* h_ffDiEMPt_JetReq_0Jet = new TH1F("ffDiEMPt_JetReq_0Jet","ff_JetReq DiEMPt with Exactly Zero Jets;DiEMPt (GeV);Events", 120,0.,600.);h_ffDiEMPt_JetReq_0Jet->Sumw2();

  TH1F* h_ggDiEMPt_JetReq_1Jet = new TH1F("ggDiEMPt_JetReq_1Jet","gg_JetReq DiEMPt with Exactly One Jet;DiEMPt (GeV);Events", 120,0.,600.);h_ggDiEMPt_JetReq_1Jet->Sumw2();
  TH1F* h_egDiEMPt_JetReq_1Jet = new TH1F("egDiEMPt_JetReq_1Jet","eg_JetReq DiEMPt with Exactly One Jet;DiEMPt (GeV);Events", 120,0.,600.);h_egDiEMPt_JetReq_1Jet->Sumw2();
  TH1F* h_eeDiEMPt_JetReq_1Jet = new TH1F("eeDiEMPt_JetReq_1Jet","ee_JetReq DiEMPt with Exactly One Jet;DiEMPt (GeV);Events", 120,0.,600.);h_eeDiEMPt_JetReq_1Jet->Sumw2();
  TH1F* h_ffDiEMPt_JetReq_1Jet = new TH1F("ffDiEMPt_JetReq_1Jet","ff_JetReq DiEMPt with Exactly One Jet;DiEMPt (GeV);Events", 120,0.,600.);h_ffDiEMPt_JetReq_1Jet->Sumw2();

  TH1F* h_ggDiEMPt_JetReq_2Jet = new TH1F("ggDiEMPt_JetReq_2Jet","gg_JetReq DiEMPt with at least Two Jets;DiEMPt (GeV);Events", 120,0.,600.);h_ggDiEMPt_JetReq_2Jet->Sumw2();
  TH1F* h_egDiEMPt_JetReq_2Jet = new TH1F("egDiEMPt_JetReq_2Jet","eg_JetReq DiEMPt with at least Two Jets;DiEMPt (GeV);Events", 120,0.,600.);h_egDiEMPt_JetReq_2Jet->Sumw2();
  TH1F* h_eeDiEMPt_JetReq_2Jet = new TH1F("eeDiEMPt_JetReq_2Jet","ee_JetReq DiEMPt with at least Two Jets;DiEMPt (GeV);Events", 120,0.,600.);h_eeDiEMPt_JetReq_2Jet->Sumw2();
  TH1F* h_ffDiEMPt_JetReq_2Jet = new TH1F("ffDiEMPt_JetReq_2Jet","ff_JetReq DiEMPt with at least Two Jets;DiEMPt (GeV);Events", 120,0.,600.);h_ffDiEMPt_JetReq_2Jet->Sumw2();

  TH1F* h_ggDiEMPt_JetReq_MoreThan2Jets = new TH1F("ggDiEMPt_JetReq_MoreThan2Jets","gg_JetReq DiEMPt with Greater Than 2 Jets;DiEMPt (GeV);Events", 120,0.,600.);h_ggDiEMPt_JetReq_MoreThan2Jets->Sumw2();
  TH1F* h_egDiEMPt_JetReq_MoreThan2Jets = new TH1F("egDiEMPt_JetReq_MoreThan2Jets","eg_JetReq DiEMPt with Greater Than 2 Jets;DiEMPt (GeV);Events", 120,0.,600.);h_egDiEMPt_JetReq_MoreThan2Jets->Sumw2();
  TH1F* h_eeDiEMPt_JetReq_MoreThan2Jets = new TH1F("eeDiEMPt_JetReq_MoreThan2Jets","ee_JetReq DiEMPt with Greater Than 2 Jets;DiEMPt (GeV);Events", 120,0.,600.);h_eeDiEMPt_JetReq_MoreThan2Jets->Sumw2();
  TH1F* h_ffDiEMPt_JetReq_MoreThan2Jets = new TH1F("ffDiEMPt_JetReq_MoreThan2Jets","ff_JetReq DiEMPt with Greater Than 2 Jets;DiEMPt (GeV);Events", 120,0.,600.);h_ffDiEMPt_JetReq_MoreThan2Jets->Sumw2();

  TH1F* h_ggDiJetPt_JetReq_0Jet = new TH1F("ggDiJetPt_JetReq_0Jet","gg_JetReq DiJetPt with Exactly Zero Jets;DiJetPt (GeV);Events", 120,0.,600.);h_ggDiJetPt_JetReq_0Jet->Sumw2();
  TH1F* h_egDiJetPt_JetReq_0Jet = new TH1F("egDiJetPt_JetReq_0Jet","eg_JetReq DiJetPt with Exactly Zero Jets;DiJetPt (GeV);Events", 120,0.,600.);h_egDiJetPt_JetReq_0Jet->Sumw2();
  TH1F* h_eeDiJetPt_JetReq_0Jet = new TH1F("eeDiJetPt_JetReq_0Jet","ee_JetReq DiJetPt with Exactly Zero Jets;DiJetPt (GeV);Events", 120,0.,600.);h_eeDiJetPt_JetReq_0Jet->Sumw2();
  TH1F* h_ffDiJetPt_JetReq_0Jet = new TH1F("ffDiJetPt_JetReq_0Jet","ff_JetReq DiJetPt with Exactly Zero Jets;DiJetPt (GeV);Events", 120,0.,600.);h_ffDiJetPt_JetReq_0Jet->Sumw2();

  TH1F* h_ggDiJetPt_JetReq_1Jet = new TH1F("ggDiJetPt_JetReq_1Jet","gg_JetReq DiJetPt with Exactly One Jet;DiJetPt (GeV);Events", 120,0.,600.);h_ggDiJetPt_JetReq_1Jet->Sumw2();
  TH1F* h_egDiJetPt_JetReq_1Jet = new TH1F("egDiJetPt_JetReq_1Jet","eg_JetReq DiJetPt with Exactly One Jet;DiJetPt (GeV);Events", 120,0.,600.);h_egDiJetPt_JetReq_1Jet->Sumw2();
  TH1F* h_eeDiJetPt_JetReq_1Jet = new TH1F("eeDiJetPt_JetReq_1Jet","ee_JetReq DiJetPt with Exactly One Jet;DiJetPt (GeV);Events", 120,0.,600.);h_eeDiJetPt_JetReq_1Jet->Sumw2();
  TH1F* h_ffDiJetPt_JetReq_1Jet = new TH1F("ffDiJetPt_JetReq_1Jet","ff_JetReq DiJetPt with Exactly One Jet;DiJetPt (GeV);Events", 120,0.,600.);h_ffDiJetPt_JetReq_1Jet->Sumw2();

  TH1F* h_ggDiJetPt_JetReq_2Jet = new TH1F("ggDiJetPt_JetReq_2Jet","gg_JetReq DiJetPt with at least Two Jets;DiJetPt (GeV);Events", 120,0.,600.);h_ggDiJetPt_JetReq_2Jet->Sumw2();
  TH1F* h_egDiJetPt_JetReq_2Jet = new TH1F("egDiJetPt_JetReq_2Jet","eg_JetReq DiJetPt with at least Two Jets;DiJetPt (GeV);Events", 120,0.,600.);h_egDiJetPt_JetReq_2Jet->Sumw2();
  TH1F* h_eeDiJetPt_JetReq_2Jet = new TH1F("eeDiJetPt_JetReq_2Jet","ee_JetReq DiJetPt with at least Two Jets;DiJetPt (GeV);Events", 120,0.,600.);h_eeDiJetPt_JetReq_2Jet->Sumw2();
  TH1F* h_ffDiJetPt_JetReq_2Jet = new TH1F("ffDiJetPt_JetReq_2Jet","ff DiJetPt with at least Two Jets;DiJetPt (GeV);Events", 120,0.,600.);h_ffDiJetPt_JetReq_2Jet->Sumw2();

  TH1F* h_ggDiJetPt_JetReq_MoreThan2Jets = new TH1F("ggDiJetPt_JetReq_MoreThan2Jets","gg_JetReq DiJet with Greater Than 2 Jets;DiJet (GeV);Events", 120,0.,600.);h_ggDiJetPt_JetReq_MoreThan2Jets->Sumw2();
  TH1F* h_egDiJetPt_JetReq_MoreThan2Jets = new TH1F("egDiJetPt_JetReq_MoreThan2Jets","eg_JetReq DiJet with Greater Than 2 Jets;DiJet (GeV);Events", 120,0.,600.);h_egDiJetPt_JetReq_MoreThan2Jets->Sumw2();
  TH1F* h_eeDiJetPt_JetReq_MoreThan2Jets = new TH1F("eeDiJetPt_JetReq_MoreThan2Jets","ee_JetReq DiJet with Greater Than 2 Jets;DiJet (GeV);Events", 120,0.,600.);h_eeDiJetPt_JetReq_MoreThan2Jets->Sumw2();
  TH1F* h_ffDiJetPt_JetReq_MoreThan2Jets = new TH1F("ffDiJetPt_JetReq_MoreThan2Jets","ff_JetReq DiJet with Greater Than 2 Jets;DiJet (GeV);Events", 120,0.,600.);h_ffDiJetPt_JetReq_MoreThan2Jets->Sumw2();

  TH1F* h_eeDiEMPt_2JetReq = new TH1F("eeDiEMPt_2JetReq","ee DiEMPt with Two Jet Requirement;DiEMPt (GeV);Events", 120,0.,600.);h_eeDiEMPt_2JetReq->Sumw2();
  TH1F* h_eeSideBandDiEMPt_JetReq = new TH1F("eeSideBandDiEMPt_JetReq","eeSideBand DiEMPt with One Jet Requirement;DiEMPt (GeV);Events", 120,0.,600.);h_eeSideBandDiEMPt_JetReq->Sumw2();
  TH1F* h_ggInvarMass = new TH1F("ggInvarMass","gg Invariant Mass;(GeV);Events", 1000,0.,1000.);
  TH1F* h_ggMETInvarMass110_135 = new TH1F("ggMETInvarMass110_135","MET in gg events with 110<InvarMass<135;(GeV);Events", 300,0.,300.);
  TH1F* h_ggMETInvarMass85_110 = new TH1F("ggMETInvarMass85_110","MET in gg events with 85<InvarMass<110;(GeV);Events", 300,0.,300.);
  TH1F* h_ggMETInvarMass135_160 = new TH1F("ggMETInvarMass135_160","MET in gg events with 135<InvarMass<160;(GeV);Events", 300,0.,300.);
  TH1F* h_ggInvarMassMET30 = new TH1F("ggInvarMassMET30","gg Invariant Mass in events with met>30;(GeV);Events", 200,0.,200.);
  TH1F* h_egInvarMassMET30 = new TH1F("egInvarMassMET30","eg Invariant Mass in events with met>30;(GeV);Events", 200,0.,200.);
  TH1F* h_ffInvarMassMET30 = new TH1F("ffInvarMassMET30","ff Invariant Mass in events with met>30;(GeV);Events", 200,0.,200.);
  TH1F* h_ggInvarMassMET40 = new TH1F("ggInvarMassMET40","gg Invariant Mass in events with met>40;(GeV);Events", 200,0.,200.);
  TH1F* h_ggInvarMassMET50 = new TH1F("ggInvarMassMET50","gg Invariant Mass in events with met>50;(GeV);Events", 200,0.,200.);
  TH1F* h_ggInvarMassMET60 = new TH1F("ggInvarMassMET60","gg Invariant Mass in events with met>60;(GeV);Events", 200,0.,200.);
  TH1F* h_ggInvarMassMET70 = new TH1F("ggInvarMassMET70","gg Invariant Mass in events with met>70;(GeV);Events", 200,0.,200.);
  TH1F* h_ggInvarMassMET80 = new TH1F("ggInvarMassMET80","gg Invariant Mass in events with met>80;(GeV);Events", 200,0.,200.);
  TH1F* h_ggInvarMassMET100 = new TH1F("ggInvarMassMET100","gg Invariant Mass in events with met>100;(GeV);Events", 200,0.,200.);
  TH1F* h_ggInvarMassMET30_JetReq = new TH1F("ggInvarMassMET30_JetReq","gg_JetReq Invariant Mass in events with met>30;(GeV);Events", 200,0.,200.);
  TH1F* h_ggInvarMassMET40_JetReq = new TH1F("ggInvarMassMET40_JetReq","gg_JetReq Invariant Mass in events with met>40;(GeV);Events", 200,0.,200.);
  TH1F* h_ggInvarMassMET50_JetReq = new TH1F("ggInvarMassMET50_JetReq","gg_JetReq Invariant Mass in events with met>50;(GeV);Events", 200,0.,200.);
  TH1F* h_ggInvarMassMET60_JetReq = new TH1F("ggInvarMassMET60_JetReq","gg_JetReq Invariant Mass in events with met>60;(GeV);Events", 200,0.,200.);
  TH1F* h_ggInvarMassMET70_JetReq = new TH1F("ggInvarMassMET70_JetReq","gg_JetReq Invariant Mass in events with met>70;(GeV);Events", 200,0.,200.);
  TH1F* h_ggInvarMassMET80_JetReq = new TH1F("ggInvarMassMET80_JetReq","gg_JetReq Invariant Mass in events with met>80;(GeV);Events", 200,0.,200.);
  TH1F* h_ggInvarMassMET100_JetReq = new TH1F("ggInvarMassMET100_JetReq","gg_JetReq Invariant Mass in events with met>100;(GeV);Events", 200,0.,200.);
  TH1F* h_egInvarMass = new TH1F("egInvarMass","eg Invariant Mass;(GeV);Events", 1000,0.,1000.);
  TH1F* h_egInvarMass_Pt25to40 = new TH1F("egInvarMass_Pt25to40","eg Invariant Mass for trailing object pt<40;(GeV);Events", 1000,0.,1000.);
  TH1F* h_egInvarMass_Pt40to45 = new TH1F("egInvarMass_Pt40to45","eg Invariant Mass for leading object 40<pt<45;(GeV);Events", 1000,0.,1000.);
  TH1F* h_egInvarMass_Pt45to50 = new TH1F("egInvarMass_Pt45to50","eg Invariant Mass for leading object 45<pt<50;(GeV);Events", 1000,0.,1000.);
  TH1F* h_egInvarMass_Pt50to60 = new TH1F("egInvarMass_Pt50to60","eg Invariant Mass for leading object 50<pt<60;(GeV);Events", 1000,0.,1000.);
  TH1F* h_egInvarMass_Pt60to80 = new TH1F("egInvarMass_Pt60to80","eg Invariant Mass for leading object 60<pt<80;(GeV);Events", 1000,0.,1000.);
  TH1F* h_egInvarMass_Pt80 = new TH1F("egInvarMass_Pt80","eg Invariant Mass for leading object pt>80;(GeV);Events", 1000,0.,1000.);
  TH1F* h_eeInvarMass = new TH1F("eeInvarMass","ee Invariant Mass for 81<InvarMass<101;(GeV);Events", 1000,0.,1000.);
  TH1F* h_eeInvarMassFullRange = new TH1F("eeInvarMassFullRange","ee Invariant Mass for all InvarMass;(GeV);Events", 1000,0.,1000.);
  TH1F* h_eeInvarMassFullRange_Pt25to40 = new TH1F("eeInvarMassFullRange_Pt25to40","ee Invariant Mass for trailing electron pt<40;(GeV);Events", 1000,0.,1000.);
  TH1F* h_eeInvarMassFullRange_Pt40to45 = new TH1F("eeInvarMassFullRange_Pt40to45","ee Invariant Mass for leading electron 40<pt<45;(GeV);Events", 1000,0.,1000.);
  TH1F* h_eeInvarMassFullRange_Pt45to50 = new TH1F("eeInvarMassFullRange_Pt45to50","ee Invariant Mass for leading electron 45<pt<50;(GeV);Events", 1000,0.,1000.);
  TH1F* h_eeInvarMassFullRange_Pt50to60 = new TH1F("eeInvarMassFullRange_Pt50to60","ee Invariant Mass for leading electron 50<pt<60;(GeV);Events", 1000,0.,1000.);
  TH1F* h_eeInvarMassFullRange_Pt60to80 = new TH1F("eeInvarMassFullRange_Pt60to80","ee Invariant Mass for leading electron 60<pt<80;(GeV);Events", 1000,0.,1000.);
  TH1F* h_eeInvarMassFullRange_Pt80 = new TH1F("eeInvarMassFullRange_Pt80","ee Invariant Mass for leading electron pt>80;(GeV);Events", 1000,0.,1000.);
  TH1F* h_ffInvarMass = new TH1F("ffInvarMass","ff Invariant Mass;(GeV);Events", 1000,0.,1000.);
  TH1F* h_ggInvarMass_JetReq = new TH1F("ggInvarMass_JetReq","gg Invariant Mass with One Jet Requirement;(GeV);Events", 1000,0.,1000.);
  TH1F* h_egInvarMass_JetReq = new TH1F("egInvarMass_JetReq","eg Invariant Mass with One Jet Requirement;(GeV);Events", 1000,0.,1000.);
  TH1F* h_eeInvarMass_JetReq = new TH1F("eeInvarMass_JetReq","ee Invariant Mass for 81<InvarMass<101 with One Jet Requirement;(GeV);Events", 1000,0.,1000.);
  TH1F* h_eeInvarMassFullRange_JetReq = new TH1F("eeInvarMassFullRange_JetReq","ee Invariant Mass for all InvarMass with One Jet Requirement;(GeV);Events", 1000,0.,1000.);
  TH1F* h_eeInvarMassFullRange_2JetReq = new TH1F("eeInvarMassFullRange_2JetReq","ee Invariant Mass for all InvarMass with Two Jet Requirement;(GeV);Events", 1000,0.,1000.);
  TH1F* h_eeInvarMass_2JetReq = new TH1F("eeInvarMass_2JetReq","ee Invariant Mass with Two Jet Requirement;(GeV);Events", 1000,0.,1000.);
  TH1F* h_ffInvarMass_JetReq = new TH1F("ffInvarMass_JetReq","ff Invariant Mass with One Jet Requirement;(GeV);Events", 1000,0.,1000.);
  TH1F* h_SigIetaIeta_allPho = new TH1F("SigIetaIeta_allPho","SigmaIetaIeta for all photons;#sigma_{i#etai#eta};Events", 100,0.,.1);

  TH1F* h_EcalIsoDR03_allPho = new TH1F("EcalIsoDR03_allPho","Ecal Isolation in DR03 cone for all photons",150,-5,25);
  TH1F* h_HcalIsoDR03_allPho = new TH1F("HcalIsoDR03_allPho","Hcal Isolation in DR03 cone for all photons",150,-5,25);
  TH1F* h_TrackIsoDR03_allPho = new TH1F("TrackIsoDR03_allPho","Track Isolation in DR03 cone for all photons",150,-5,25);
  TH1F* h_CombIsoDR03_allPho = new TH1F("CombIsoDR03_allPho","Combined Isolation in DR03 cone for all photons",600,-5,100);
  TH1F* h_EcalIsoDR03RhoCorr_allPho = new TH1F("EcalIsoDR03RhoCorr_allPho","RhoCorrected Ecal Isolation in DR03 cone for all photons",150,-5,25);
  TH1F* h_HcalIsoDR03RhoCorr_allPho = new TH1F("HcalIsoDR03RhoCorr_allPho","RhoCorrected Hcal Isolation in DR03 cone for all photons",150,-5,25);
  TH1F* h_CombIsoDR03RhoCorr_allPho = new TH1F("CombIsoDR03RhoCorr_allPho","RhoCorrected Combined Isolation in DR03 cone for all photons",600,-5,100);

  //TH1F* h_SigIetaIeta_twoPho = new TH1F("SigIetaIeta_twoPho","SigmaIetaIeta for two highest p_{T} photons;#sigma_{i#etai#eta};Events", 100,0.,.1);
  TH1I* h_PixelSeeds = new TH1I("PixelSeeds","Number of pixel seeds for photon candidate after ecal iso, hcal iso, h/e",150,0,150);
  TH1I* h_eePixelSeeds = new TH1I("eePixelSeeds","Number of pixel seeds for ee sample",150,0,150);
  TH1I* h_egPixelSeeds = new TH1I("egPixelSeeds","Number of pixel seeds for eg sample",150,0,150);
  //TH1I* h_ffPixelSeeds = new TH1I("ffPixelSeeds","Number of pixel seeds for ff sample",150,0,150);
  TH1I* h_PixelSeeds_JetReq = new TH1I("PixelSeeds_JetReq","Number of pixel seeds for photon candidate after ecal iso, hcal iso, h/e with One Jet Requirement",150,0,150);
  TH1I* h_eePixelSeeds_JetReq = new TH1I("eePixelSeeds_JetReq","Number of pixel seeds for ee sample with One Jet Requirement",150,0,150);
  TH1I* h_egPixelSeeds_JetReq = new TH1I("egPixelSeeds_JetReq","Number of pixel seeds for eg sample with One Jet Requirement",150,0,150);
  //TH1I* h_ffPixelSeeds_JetReq = new TH1I("ffPixelSeeds_JetReq","Number of pixel seeds for ff sample with One Jet Requirement",150,0,150);
  TH1F* h_ggSigIetaIeta = new TH1F("ggSigIetaIeta","gg SigmaIetaIeta;#sigma_{i#etai#eta};Events", 20,0.,.02);
  TH1F* h_egSigIetaIeta = new TH1F("egSigIetaIeta","eg SigmaIetaIeta;#sigma_{i#etai#eta};Events", 20,0.,.02);
  TH1F* h_eeSigIetaIeta = new TH1F("eeSigIetaIeta","ee SigmaIetaIeta;#sigma_{i#etai#eta};Events", 20,0.,.02);
  TH1F* h_ffSigIetaIeta = new TH1F("ffSigIetaIeta","ff SigmaIetaIeta;#sigma_{i#etai#eta};Events", 20,0.,.02);

  TH1F* h_ggMet = new TH1F("ggMet","Missing transverse energy in events with two photons with Et>40,25 GeV;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_egMet = new TH1F("egMet","Missing transverse energy in events with one photon and one electron with Et>40,25 GeV;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeMet = new TH1F("eeMet","Missing transverse energy in events with two electrons with Et>40,25 GeV;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_ggMet_JetFail = new TH1F("ggMet_JetFail","Missing transverse energy in gg events that have a 50GeV pfJet failing loose jet id;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_egMet_JetFail = new TH1F("egMet_JetFail","Missing transverse energy in eg events that have a 50GeV pfJet failing loose jet id;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeMet_JetFail = new TH1F("eeMet_JetFail","Missing transverse energy in ee events that have a 50GeV pfJet failing loose jet id;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_ffMet_JetFail = new TH1F("ffMet_JetFail","Missing transverse energy in ff events that have a 50GeV pfJet failing loose jet id;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeSideBandMet = new TH1F("eeSideBandMet","Missing transverse energy in events with two electrons in sideband with Et>40,25 GeV;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_ffMet = new TH1F("ffMet","Missing transverse energy in events with two fake photons with Et>40,25 GeV;#slash{E}_{T} (GeV);Events",200,0.,1000.);

  TH1F* h_eeSidebandMet_reweight_binned = new TH1F("eeSidebandMet_reweight_binned","Missing transverse energy in events with two electrons in sideband with Et>40,25 GeV reweighted from gg/ee diEMPt bins ;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeSidebandLowMet_reweight_binned = new TH1F("eeSidebandLowMet_reweight_binned","Missing transverse energy in events with two electrons in Low sideband with Et>40,25 GeV reweighted from gg/ee diEMPt bins ;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeSidebandHighMet_reweight_binned = new TH1F("eeSidebandHighMet_reweight_binned","Missing transverse energy in events with two electrons in Highsideband with Et>40,25 GeV reweighted from gg/ee diEMPt bins ;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeSidebandLowMet_reweightJet_binned = new TH1F("eeSidebandLowMet_reweightJet_binned","Missing transverse energy in events with two electrons in Low sideband with Et>40,25 GeV reweighted from gg/ee diJetPt bins ;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeSidebandHighMet_reweightJet_binned = new TH1F("eeSidebandHighMet_reweightJet_binned","Missing transverse energy in events with two electrons in Highsideband with Et>40,25 GeV reweighted from gg/ee diJetPt bins ;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeSidebandLowMet_reweightJet_binned_JetReq = new TH1F("eeSidebandLowMet_reweightJet_binned_JetReq","Missing transverse energy in events with two electrons in Low sideband with Et>40,25 GeV reweighted from gg/ee JetReq diJetPt bins ;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeSidebandHighMet_reweightJet_binned_JetReq = new TH1F("eeSidebandHighMet_reweightJet_binned_JetReq","Missing transverse energy in events with two electrons in Highsideband with Et>40,25 GeV reweighted from gg/ee JetReq diJetPt bins ;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeSidebandLowMet_reweight_binned_JetReq = new TH1F("eeSidebandLowMet_reweight_binned_JetReq","Missing transverse energy in events with two electrons in Low sideband with Et>40,25 GeV reweighted from gg/ee JetReq diEMPt bins ;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeSidebandHighMet_reweight_binned_JetReq = new TH1F("eeSidebandHighMet_reweight_binned_JetReq","Missing transverse energy in events with two electrons in Highsideband with Et>40,25 GeV reweighted from gg/ee JetReq diEMPt bins ;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeMet_reweight_binned = new TH1F("eeMet_reweight_binned","Missing transverse energy in events with two electrons with Et>40,25 GeV reweighted from gg/ee diEMPt bins;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeMet_reweightJet_binned = new TH1F("eeMet_reweightJet_binned","Missing transverse energy in events with two electrons with Et>40,25 GeV reweighted from gg/ee diJetPt bins;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_ffMet_reweight_binned = new TH1F("ffMet_reweight_binned","Missing transverse energy in events with two fake photons with Et>40,25 GeV reweighted from gg/ff diEMPt bins;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_ffMet_reweightJet_binned = new TH1F("ffMet_reweightJet_binned","Missing transverse energy in events with two electrons with Et>40,25 GeV reweighted from gg/ff diJetPt bins;#slash{E}_{T} (GeV);Events",200,0.,1000.);

  TH1F* h_eeMet_reweight_binned_JetReq = new TH1F("eeMet_reweight_binned_JetReq","Missing transverse energy in events with two electrons with Et>40,25 GeV reweighted from gg/ee JetReq diEMPt bins;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeMet_reweightJet_binned_JetReq = new TH1F("eeMet_reweightJet_binned_JetReq","Missing transverse energy in events with two electrons with Et>40,25 GeV reweighted from gg/ee JetReq diJetPt bins;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_ffMet_reweight_binned_JetReq = new TH1F("ffMet_reweight_binned_JetReq","Missing transverse energy in events with two fake photons with Et>40,25 GeV reweighted from gg/ff JetReq diEMPt bins;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_ffMet_reweightJet_binned_JetReq = new TH1F("ffMet_reweightJet_binned_JetReq","Missing transverse energy in events with two electrons with Et>40,25 GeV reweighted from gg/ff JetReq diJetPt bins;#slash{E}_{T} (GeV);Events",200,0.,1000.);

  TH1F* h_ggMet_JetReq = new TH1F("ggMet_JetReq","Missing transverse energy in events with two photons with Et>40,25 GeV with One Jet Requirement;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_egMet_JetReq = new TH1F("egMet_JetReq","Missing transverse energy in events with one photon and one electron with Et>40,25 GeV with One Jet Requirement;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeMet_JetReq = new TH1F("eeMet_JetReq","Missing transverse energy in events with two electrons with Et>40,25 GeV with One Jet Requirement;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeMet_2JetReq = new TH1F("eeMet_2JetReq","Missing transverse energy in events with two electrons with Et>40,25 GeV with Two Jet Requirement;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeSideBandMet_JetReq = new TH1F("eeSideBandMet_JetReq","Missing transverse energy in events with two electrons in sideband with Et>40,25 GeV with One Jet Requirement;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_ffMet_JetReq = new TH1F("ffMet_JetReq","Missing transverse energy in events with two fake photons with Et>40,25 GeV with One Jet Requirement;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeSidebandMet_reweight_binned_JetReq = new TH1F("eeSidebandMet_reweight_binned_JetReq","Missing transverse energy in events with two electrons in sideband with Et>40,25 GeV with One Jet Requirement reweighted from gg/ee diEMPt bins;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeMet_reweight_binned_2JetReq = new TH1F("eeMet_reweight_binned_2JetReq","Missing transverse energy in events with two electrons with Et>40,25 GeV with Two Jet Requirement reweighted from gg/ee diEMPt bins;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1I* h_NPhoCands = new TH1I("NPhoCands","Number of photon candidates after Et, ecalIso, hcalIso, h/e, trackIso, SigIetaIeta cuts; ;Events",5,0,5); 
  TH1I* h_NFakeCands = new TH1I("NFakeCands","Number of fake candidates after Et, ecalIso, hcalIso, h/e, !trackIso||!SigIetaIeta cuts; ;Events",5,0,5);
  TH2F* h_EcalIsoVsNVertex = new TH2F("EcalIsoVsNVertex","Ecal Isolation of barrel photons VS Number of Vertices after Et, h/e cuts;nVertex;ecalRecHitSumEtConeDR04",20,0,20,100,-5.,20.);
  TH2F* h_HcalIsoVsNVertex = new TH2F("HcalIsoVsNVertex","Hcal Isolation of barrel photons VS Number of Vertices after Et, h/e cuts;nVertex;hcalTowerSumEtConeDR04",20,0,20,100,-5.,20.);
  TH2F* h_TrackIsoVsNVertex = new TH2F("TrackIsoVsNVertex","Tracker Isolation of barrel photons VS Number of Vertices after Et, h/e cuts;nVertex;trkSumPtHollowConeDR04",20,0,20,100,-5.,20.);
  TH1F* h_EcalIso = new TH1F("EcalIso","Ecal Isolation for barrel photons after Et, hcal iso, h/e cuts;ecalRecHitSumEtConeDR04;Events",50,-5.,20.);
  TH1F* h_HcalIso = new TH1F("HcalIso","Hcal Isolation for barrel photons after Et, ecal iso, h/e cuts;hcalTowerSumEtConeDR04;Events",50,-5.,20.);
  TH1F* h_TrackIso = new TH1F("TrackIso","Tracker Isolation for barrel photons after Et, ecal iso, hcal iso, h/e cuts;trkSumPtHollowConeDR04;Events",50,-5.,20.);
  TH1F* h_EcalIso_Preselection = new TH1F("EcalIso_Preselection","Ecal Isolation for barrel photons after Et, ecal iso, hcal iso, h/e cuts;ecalRecHitSumEtConeDR04;Events",50,-5.,20.);
  TH1F* h_HcalIso_Preselection = new TH1F("HcalIso_Preselection","Hcal Isolation for barrel photons after Et, ecal iso, hcal iso, h/e cuts;hcalTowerSumEtConeDR04;Events",50,-5.,20.);
  TH1F* h_TrackIso_Preselection = new TH1F("TrackIso_Preselection","Tracker Isolation for barrel photons after Et, ecal iso, hcal iso, h/e cuts;trkSumPtHollowConeDR04;Events",50,-5.,20.);
  TH1F* h_EcalIso_Preselection_NVertex_Corrected = new TH1F("EcalIso_Preselection_NVertex_Corrected","Ecal NVertex Corrected Isolation for barrel photons after Et, ecal iso, hcal iso, h/e cuts;ecalRecHitSumEtConeDR04;Events",50,-5.,20.);
  TH1F* h_HcalIso_Preselection_NVertex_Corrected = new TH1F("HcalIso_Preselection_NVertex_Corrected","Hcal NVertex Corrected Isolation for barrel photons after Et, ecal iso, hcal iso, h/e cuts;hcalTowerSumEtConeDR04;Events",50,-5.,20.);
  TH1F* h_TrackIso_Preselection_NVertex_Corrected = new TH1F("TrackIso_Preselection_NVertex_Corrected","Tracker NVertex Corrected Isolation for barrel photons after Et, ecal iso, hcal iso, h/e cuts;trkSumPtHollowConeDR04;Events",50,-5.,20.);
  TH1F* h_EcalIso_Preselection_Rho_Corrected = new TH1F("EcalIso_Preselection_Rho_Corrected","Ecal Rho Corrected Isolation for barrel photons after Et, ecal iso, hcal iso, h/e cuts;ecalRecHitSumEtConeDR04;Events",50,-5.,20.);
  TH1F* h_HcalIso_Preselection_Rho_Corrected = new TH1F("HcalIso_Preselection_Rho_Corrected","Hcal Rho Corrected Isolation for barrel photons after Et, ecal iso, hcal iso, h/e cuts;hcalTowerSumEtConeDR04;Events",50,-5.,20.);
  TH1F* h_TrackIso_Preselection_Rho_Corrected = new TH1F("TrackIso_Preselection_Rho_Corrected","Tracker Rho Corrected Isolation for barrel photons after Et, ecal iso, hcal iso, h/e cuts;trkSumPtHollowConeDR04;Events",50,-5.,20.);
  /* TH2F* h_FailEB_PhiVsEta = new TH2F("FailEB_PhiVsEta","Phi VS Eta for Photons that fail isEB;Eta;Phi",1000,-5.,5.,1000,-5.,5.);
     TH2F* h_FailEBEtaGap_PhiVsEta = new TH2F("FailEBEtaGap_PhiVsEta","Phi VS Eta for Photons that fail isEBEtaGap;Eta;Phi",1000,-5.,5.,1000,-5.,5.);
     TH2F* h_FailEBPhiGap_PhiVsEta = new TH2F("FailEBPhiGap_PhiVsEta","Phi VS Eta for Photons that fail isEBPhiGap;Eta;Phi",1000,-5.,5.,1000,-5.,5.);
     TH2F* h_FailEBEEGap_PhiVsEta = new TH2F("FailEBEEGap_PhiVsEta","Phi VS Eta for Photons that fail isEBEEGap;Eta;Phi",1000,-5.,5.,1000,-5.,5.);*/
  TH1F* h_SeedTime = new TH1F("SeedTime","SeedTime",4200,-200.,4000.);
  TH1F* h_SeedTime_afterR9 = new TH1F("SeedTime_afterR9","SeedTime after R9 cut",4200,-200.,4000.);
  TH2F* h_SeedTimeVsEta = new TH2F("SeedTimeVsEta","SeedTime VS Photon Eta",350,-3.5,3.5,4200,-200.,4000.);
  TH2F* h_SeedTimeVSE = new TH2F("SeedTimeVSE","Photon Energy VS SeedTime ",4200,-200.,4000.,200,0.,1000.);
  TH2F* h_SeedTimeVSE_afterR9 = new TH2F("SeedTimeVSE_afterR9","Photon Energy VS SeedTime after R9 cut",4200,-200.,4000.,200,0.,1000.);
  TH2F* h_MetVsSeedTime = new TH2F("MetVsSeedTime","MET VS Photon SeedTime ",4200,-200.,4000.,200,0.,1000.);
  TH2F* h_MetVsSeedTime_gg = new TH2F("MetVsSeedTime_gg","MET VS Photon SeedTime - gg",4200,-200.,4000.,200,0.,1000.);
  TH2F* h_MetVsSeedTime_ee = new TH2F("MetVsSeedTime_ee","MET VS Photon SeedTime - ee",4200,-200.,4000.,200,0.,1000.);
  TH2F* h_MetVsSeedTime_ff = new TH2F("MetVsSeedTime_ff","MET VS Photon SeedTime - ff",4200,-200.,4000.,200,0.,1000.);
  TH2F* h_SeedTimeVsEta_gg = new TH2F("SeedTimeVsEta_gg","SeedTime VS Photon Eta - gg",170,-1.7,1.7,4200,-200.,4000.);
  TH2F* h_SeedTimeVsEta_eg = new TH2F("SeedTimeVsEta_eg","SeedTime VS Photon Eta - eg",170,-1.7,1.7,4200,-200.,4000.);
  TH2F* h_SeedTimeVsEta_ee = new TH2F("SeedTimeVsEta_ee","SeedTime VS Photon Eta - ee",170,-1.7,1.7,4200,-200.,4000.);
  TH2F* h_SeedTimeVsEta_ff = new TH2F("SeedTimeVsEta_ff","SeedTime VS Photon Eta - ff",170,-1.7,1.7,4200,-200.,4000.);
  TH2F* h_SeedTimeVsEta_gg_JetReq = new TH2F("SeedTimeVsEta_gg_JetReq","SeedTime VS Photon Eta - gg_JetReq",170,-1.7,1.7,4200,-200.,4000.);
  TH2F* h_SeedTimeVsEta_eg_JetReq = new TH2F("SeedTimeVsEta_eg_JetReq","SeedTime VS Photon Eta - eg_JetReq",170,-1.7,1.7,4200,-200.,4000.);
  TH2F* h_SeedTimeVsEta_ee_JetReq = new TH2F("SeedTimeVsEta_ee_JetReq","SeedTime VS Photon Eta - ee_JetReq",170,-1.7,1.7,4200,-200.,4000.);
  TH2F* h_SeedTimeVsEta_ff_JetReq = new TH2F("SeedTimeVsEta_ff_JetReq","SeedTime VS Photon Eta - ff_JetReq",170,-1.7,1.7,4200,-200.,4000.);

  TH1F* h_TrigPhosEta = new TH1F("TrigPhosEta","Eta of all photons passing triggers",400,-4,4);
  TH1F* h_TrigPhosEta_TopTwo = new TH1F("TrigPhosEta_TopTwo","Eta of top 2 p_{T} photons passing triggers",400,-4,4);

  TH1F* h_Pho_ChargedHadronIso = new TH1F("Pho_ChargedHadronIso","",500,0.,100.);
  TH1F* h_Pho_NeutralHadronIso = new TH1F("Pho_NeutralHadronIso","",500,0.,100.);
  TH1F* h_Pho_PhotonIso = new TH1F("Pho_PhotonIso","",500,0.,100.);
  TH1F* h_Pho_PfCombinedIso = new TH1F("Pho_PfCombinedIso","",1500,0.,300.);
  TH1F* h_Pho_ChargedHadronIsoDeposit = new TH1F("Pho_ChargedHadronIsoDeposit","",500,0.,100.);
  TH1F* h_Pho_NeutralHadronIsoDeposit = new TH1F("Pho_NeutralHadronIsoDeposit","",500,0.,100.);
  TH1F* h_Pho_PhotonIsoDeposit = new TH1F("Pho_PhotonIsoDeposit","",500,0.,100.);
  TH1F* h_Pho_PfCombinedIsoDeposit = new TH1F("Pho_PfCombinedIsoDeposit","",1500,0.,300.);
  TH1F* h_Pho_CombinedIsoDR03 = new TH1F("Pho_CombinedIsoDR03","",1500,0.,300.);

  TH1F* h_r9 = new TH1F("r9","R9 of all photons",250,0.,5.);
  TH1F* h_r9_gg = new TH1F("r9_gg","R9 of all photons in gg sample",250,0.,5.);
  TH1F* h_r9_eg = new TH1F("r9_eg","R9 of all photons in eg sample",250,0.,5.);
  TH1F* h_r9_ee = new TH1F("r9_ee","R9 of all photons in ee sample",250,0.,5.);
  TH1F* h_r9_ff = new TH1F("r9_ff","R9 of all photons in ff sample",250,0.,5.);
  TH1F* h_R9fromR9trig = new TH1F("R9fromR9trig","R9 of all photons from R9 trigger",250,0.,5.);
  TH1F* h_rho = new TH1F("rho","Fastjet correction rho for all events that pass JSON, HLT",150,0,150);
  TH1F* h_rhoBarrel = new TH1F("rhoBarrel","Fastjet correction rhoBarrel for all events that pass JSON, HLT",150,0,150);
  TH1F* h_EcalIsoDR03 = new TH1F("EcalIsoDR03","Ecal Isolation in DR03 cone for all photons after JSON, HLT, R9, Et>30, h/e",150,-5,25);
  TH1F* h_HcalIsoDR03 = new TH1F("HcalIsoDR03","Hcal Isolation in DR03 cone for all photons after JSON, HLT, R9, Et>30, h/e",150,-5,25);
  TH1F* h_TrackIsoDR03 = new TH1F("TrackIsoDR03","Track Isolation in DR03 cone for all photons after JSON, HLT, R9, Et>30, h/e",150,-5,25);
  TH1F* h_EcalIsoDR04 = new TH1F("EcalIsoDR04","Ecal Isolation in DR04 cone for all photons after JSON, HLT, R9, Et>30, h/e",150,-5,25);
  TH1F* h_HcalIsoDR04 = new TH1F("HcalIsoDR04","Hcal Isolation in DR04 cone for all photons after JSON, HLT, R9, Et>30, h/e",150,-5,25);
  TH1F* h_TrackIsoDR04 = new TH1F("TrackIsoDR04","Track Isolation in DR04 cone for all photons after JSON, HLT, R9, Et>30, h/e",150,-5,25);
  TH1F* h_EcalIsoDR03_gg = new TH1F("EcalIsoDR03_gg","Ecal Isolation in DR03 cone for all photons in gg sample",150,-5,25);
  TH1F* h_EcalIsoDR03RhoCorr_eg = new TH1F("EcalIsoDR03RhoCorr_eg","Ecal Isolation in DR03 cone for all photons in eg sample",150,-5,25);
  TH1F* h_EcalIsoDR03RhoCorr_ee = new TH1F("EcalIsoDR03RhoCorr_ee","Ecal Isolation in DR03 cone for all photons in ee sample",150,-5,25);
  TH1F* h_EcalIsoDR03RhoCorr_ff = new TH1F("EcalIsoDR03RhoCorr_ff","Ecal Isolation in DR03 cone for all photons in ff sample",150,-5,25);
  TH1F* h_HcalIsoDR03_gg = new TH1F("HcalIsoDR03_gg","Hcal Isolation in DR03 cone for all photons in gg sample",150,-5,25);
  TH1F* h_HcalIsoDR03RhoCorr_eg = new TH1F("HcalIsoDR03RhoCorr_eg","Hcal Isolation in DR03 cone for all photons in eg sample",150,-5,25);
  TH1F* h_HcalIsoDR03RhoCorr_ee = new TH1F("HcalIsoDR03RhoCorr_ee","Hcal Isolation in DR03 cone for all photons in ee sample",150,-5,25);
  TH1F* h_HcalIsoDR03RhoCorr_ff = new TH1F("HcalIsoDR03RhoCorr_ff","Hcal Isolation in DR03 cone for all photons in ff sample",150,-5,25);
  TH1F* h_TrackIsoDR03_gg = new TH1F("TrackIsoDR03_gg","Track Isolation in DR03 cone for all photons in gg sample",150,-5,25);
  TH1F* h_TrackIsoDR03RhoCorr_eg = new TH1F("TrackIsoDR03RhoCorr_eg","Track Isolation in DR03 cone for all photons in eg sample",150,-5,25);
  TH1F* h_TrackIsoDR03RhoCorr_ee = new TH1F("TrackIsoDR03RhoCorr_ee","Track Isolation in DR03 cone for all photons in ee sample",150,-5,25);
  TH1F* h_TrackIsoDR03RhoCorr_ff = new TH1F("TrackIsoDR03RhoCorr_ff","Track Isolation in DR03 cone for all photons in ff sample",150,-5,25);

  TH1F* h_EcalIsoDR03RhoCorr_gg = new TH1F("EcalIsoDR03RhoCorr_gg","RhoCorrected Ecal Isolation in DR03 cone for all photons in gg sample",150,-5,25);
  TH1F* h_HcalIsoDR03RhoCorr_gg = new TH1F("HcalIsoDR03RhoCorr_gg","RhoCorrected Hcal Isolation in DR03 cone for all photons in gg sample",150,-5,25);
  TH1F* h_TrackIsoDR03RhoCorr_gg = new TH1F("TrackIsoDR03RhoCorr_gg","RhoCorrected Track Isolation in DR03 cone for all photons in gg sample",150,-5,25);
  TH1F* h_CombIsoDR03RhoCorr_gg = new TH1F("CombIsoDR03RhoCorr_gg","RhoCorrected Combined Isolation in DR03 cone for all photons in gg sample",600,-5,100);

  TH1F* h_CombIsoDR03_gg = new TH1F("CombIsoDR03_gg","RhoCorrected Combined Isolation in DR03 cone for all photons in gg sample",600,-5,100);
  TH1F* h_CombIsoDR03RhoCorr_eg = new TH1F("CombIsoDR03RhoCorr_eg","RhoCorrected Combined Isolation in DR03 cone for all photons in eg sample",600,-5,100);
  TH1F* h_CombIsoDR03RhoCorr_ee = new TH1F("CombIsoDR03RhoCorr_ee","RhoCorrected Combined Isolation in DR03 cone for all photons in ee sample",600,-5,100);
  TH1F* h_CombIsoDR03RhoCorr_ff = new TH1F("CombIsoDR03RhoCorr_ff","RhoCorrected Combined Isolation in DR03 cone for all photons in ff sample",600,-5,100);
  TH1F* h_rho_gg = new TH1F("rho_gg","Fastjet correction rho for all events in gg sample",200,0,40);
  TH1F* h_rho_eg = new TH1F("rho_eg","Fastjet correction rho for all events in eg sample",200,0,40);
  TH1F* h_rho_ee = new TH1F("rho_ee","Fastjet correction rho for all events in ee sample",200,0,40);
  TH1F* h_rho_ff = new TH1F("rho_ff","Fastjet correction rho for all events in ff sample",200,0,40);
  TH1F* h_NVertex = new TH1F("NVertex","Number of good vertices for all events",40,0,40);
  TH1F* h_NVertex_gg = new TH1F("NVertex_gg","Number of good vertices for all events in gg sample",40,0,40);
  TH1F* h_NVertex_eg = new TH1F("NVertex_eg","Number of good vertices for all events in eg sample",40,0,40);
  TH1F* h_NVertex_ee = new TH1F("NVertex_ee","Number of good vertices for all events in ee sample",40,0,40);
  TH1F* h_NVertex_ff = new TH1F("NVertex_ff","Number of good vertices for all events in ff sample",40,0,40);

  Double_t MetBins[13]={0,5,10,15,20,25,30,35,40,50,70,100,275,};
  TH2F* h_NVertexVsMET_gg = new TH2F("NVertexVsMET_gg","Number of good vertices VS MET for all events in gg sample",12,MetBins,40,0,40);
  TH2F* h_NVertexVsMET_eg = new TH2F("NVertexVsMET_eg","Number of good vertices VS MET for all events in eg sample",12,MetBins,40,0,40);
  TH2F* h_NVertexVsMET_ee = new TH2F("NVertexVsMET_ee","Number of good vertices VS MET for all events in ee sample",12,MetBins,40,0,40);
  TH2F* h_NVertexVsMET_ff = new TH2F("NVertexVsMET_ff","Number of good vertices VS MET for all events in ff sample",12,MetBins,40,0,40);

  //-------------------
  TH1F *h_DenomVsNVertex = new TH1F("DenomVsNVertex","",40,0,40);
  TH1F *h_DenomVsNVertex_afterLooseHE = new TH1F("DenomVsNVertex_afterHE","",40,0,40);
  TH1F *h_NumerVsNVertex = new TH1F("NumerVsNVertex","",40,0,40);
  TH1F *h_RhoCorrectedNumerVsNVertex = new TH1F("RhoCorrectedNumerVsNVertex","",40,0,40);
  TH1F *h_NVertexCorrectedNumerVsNVertex = new TH1F("NVertexCorrectedNumerVsNVertex","",40,0,40);
  TH1F *h_DenomVsEt = new TH1F("DenomVsEt","",100,0,500);
  TH1F *h_DenomVsEt_afterLooseHE = new TH1F("DenomVsEt_afterHE","",100,0,500);
  TH1F *h_NumerVsEt = new TH1F("NumerVsEt","",100,0,500);
  TH1F *h_RhoCorrectedNumerVsEt = new TH1F("RhoCorrectedNumerVsEt","",100,0,500);
  TH1F *h_NVertexCorrectedNumerVsEt = new TH1F("NVertexCorrectedNumerVsEt","",100,0,500);
  TH1F *h_DenomVsEta = new TH1F("DenomVsEta","",15,-1.5,1.5);
  TH1F *h_DenomVsEta_afterLooseHE = new TH1F("DenomVsEta_afterHE","",15,-1.5,1.5);
  TH1F *h_NumerVsEta = new TH1F("NumerVsEta","",15,-1.5,1.5);
  TH1F *h_RhoCorrectedNumerVsEta = new TH1F("RhoCorrectedNumerVsEta","",15,-1.5,1.5);
  TH1F *h_NVertexCorrectedNumerVsEta = new TH1F("NVertexCorrectedNumerVsEta","",15,-1.5,1.5);
  TH1F *h_DenomVsMET = new TH1F("DenomVsMET","",100,0,500);
  TH1F *h_DenomVsMET_afterLooseHE = new TH1F("DenomVsMET_afterHE","",100,0,500);
  TH1F *h_NumerVsMET = new TH1F("NumerVsMET","",100,0,500);
  TH1F *h_RhoCorrectedNumerVsMET = new TH1F("RhoCorrectedNumerVsMET","",100,0,500);
  TH1F *h_NVertexCorrectedNumerVsMET = new TH1F("NVertexCorrectedNumerVsMET","",100,0,500);
  //--------------------
  TH2F *h_DiEMPtVsMet_gg = new TH2F("DiEMPtVsMet_gg","",100,0.,500.,120,0.,600.);
  TH2F *h_DiEMPtVsMet_eg = new TH2F("DiEMPtVsMet_eg","",100,0.,500.,120,0.,600.);
  TH2F *h_DiEMPtVsMet_ee = new TH2F("DiEMPtVsMet_ee","",100,0.,500.,120,0.,600.);
  TH2F *h_DiEMPtVsMet_eeSideBand = new TH2F("DiEMPtVsMet_eeSideBand","",100,0.,500.,120,0.,600.);
  TH2F *h_DiEMPtVsMet_ff = new TH2F("DiEMPtVsMet_ff","",100,0.,500.,120,0.,600.);
  //------------------
  TH1F *h_BgroundCombIsoDR03Nminus3 = new TH1F("BgroundCombIsoDR03Nminus3","Combined Iso after Et>30,R9,h/e,sIetaIeta,Pixel=0,MET<30;combined Isolation (DR03 cone);Events",210,-5.,100.);
  TH1F *h_BgroundCombIsoDR04Nminus3 = new TH1F("BgroundCombIsoDR04Nminus3","Combined Iso after Et>30,R9,h/e,sIetaIeta,Pixel=0,MET<30;combined Isolation (DR04 cone);Events",210,-5.,100.);
  TH1F *h_BgroundCombIsoDR03Nminus3Trail = new TH1F("BgroundCombIsoDR03Nminus3Trail","Combined Iso of trail leg after Et>30,R9,h/e,sIetaIeta,Pixel=0,MET<30;combined Isolation (DR03 cone);Events",210,-5.,100.);
  TH1F *h_BgroundCombIsoDR04Nminus3Trail = new TH1F("BgroundCombIsoDR04Nminus3Trail","Combined Iso of trail leg after Et>30,R9,h/e,sIetaIeta,Pixel=0,MET<30;combined Isolation (DR04 cone);Events",210,-5.,100.);
  TH1F *h_BgroundEcalIsoDR03Nminus1 = new TH1F("BgroundEcalIsoDR03Nminus1","EcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,hcalIso,trackIso,MET<30;ecalRecHitSumEtConeDR03;Events",150,-5,25);
  TH1F *h_BgroundHcalIsoDR03Nminus1 = new TH1F("BgroundHcalIsoDR03Nminus1","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,trackIso,MET<30;hcalTowerSumEtConeDR03();Events",150,-5,25);
  TH1F *h_BgroundTrackIsoDR03Nminus1 = new TH1F("BgroundTrackIsoDR03Nminus1","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,hcalIso,MET<30;trkSumPtHollowConeDR03;Events",150,-5,25);
  //for trailing leg
  TH1F *h_BgroundCombIsoDR03Nminus3trail = new TH1F("BgroundCombIsoDR03Nminus3trail","Combined Iso after Et>30,R9,h/e,sIetaIeta,Pixel=0,MET<30;combined Isolation (DR03 cone);Events",210,-5.,100.);
  TH1F *h_BgroundCombIsoDR04Nminus3trail = new TH1F("BgroundCombIsoDR04Nminus3trail","Combined Iso after Et>30,R9,h/e,sIetaIeta,Pixel=0,MET<30;combined Isolation (DR04 cone);Events",210,-5.,100.);
  TH1F *h_BgroundEcalIsoDR03Nminus1trail = new TH1F("BgroundEcalIsoDR03Nminus1trail","EcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,hcalIso,trackIso,MET<30;ecalRecHitSumEtConeDR03;Events",150,-5,25);
  TH1F *h_BgroundHcalIsoDR03Nminus1trail = new TH1F("BgroundHcalIsoDR03Nminus1trail","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,trackIso,MET<30;hcalTowerSumEtConeDR03();Events",150,-5,25);
  TH1F *h_BgroundTrackIsoDR03Nminus1trail = new TH1F("BgroundTrackIsoDR03Nminus1trail","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,hcalIso,MET<30;trkSumPtHollowConeDR03;Events",150,-5,25);
  
  TH1F *h_EcalIsoDR03Nminus1 = new TH1F("EcalIsoDR03Nminus1","EcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,hcalIso,trackIso;ecalRecHitSumEtConeDR03;Events",150,-5,25);
  TH1F *h_EcalIsoDR04Nminus1 = new TH1F("EcalIsoDR04Nminus1","EcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,hcalIso,trackIso;ecalRecHitSumEtConeDR04;Events",150,-5,25);
  TH1F *h_HcalIsoDR03Nminus3 = new TH1F("HcalIsoDR03Nminus3","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0;hcalTowerSumEtConeDR03();Events",150,-5,25);
  TH1F *h_HcalIsoDR03Nminus1 = new TH1F("HcalIsoDR03Nminus1","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,trackIso;hcalTowerSumEtConeDR03();Events",150,-5,25);
  TH1F *h_HcalIsoDR04Nminus1 = new TH1F("HcalIsoDR04Nminus1","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,trackIso;hcalTowerSumEtConeDR04();Events",150,-5,25);
  TH1F *h_TrackIsoDR03Nminus3 = new TH1F("TrackIsoDR03Nminus3","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0;trkSumPtHollowConeDR03;Events",150,-5,25);
  TH1F *h_TrackIsoDR03Nminus1 = new TH1F("TrackIsoDR03Nminus1","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,hcalIso;trkSumPtHollowConeDR03;Events",150,-5,25);
  TH1F *h_TrackIsoDR04Nminus1 = new TH1F("TrackIsoDR04Nminus1","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,hcalIso;trkSumPtHollowConeDR04;Events",150,-5,25);
  TH1F *h_CombIsoNminus4 = new TH1F("CombIsoNminus4","RhoCorrected Combined Iso after Et>30,R9,h/e;Combined Isolation;Events",210,-5.,100.);
  TH1F *h_SigIetaIetaNminus4 = new TH1F("SigIetaIetaNminus4","SigmaIetaIeta after Et>30,R9,h/e;#sigma_{i#etai#eta};Events",20,0.,.02);
  TH2F *h_EcalIsoDR03VsCombIsoDR03Nminus4 = new TH2F("EcalIsoDR03VsCombIsoDR03Nminus4","EcalIsoDR03 Vs CombIsoDR03 - both RhoCorrected - after eta,Et>30,R9,h/e",150,-5.,25.,210,-5.,100.); 
  TH2F *h_HcalIsoDR03VsCombIsoDR03Nminus4 = new TH2F("HcalIsoDR03VsCombIsoDR03Nminus4","HcalIsoDR03 Vs CombIsoDR03 - both RhoCorrected - after eta,Et>30,R9,h/e",150,-5.,25.,210,-5.,100.); 
  TH2F *h_TrackIsoDR03VsCombIsoDR03Nminus4 = new TH2F("TrackIsoDR03VsCombIsoDR03Nminus4","TrackIsoDR03 Vs CombIsoDR03 - both RhoCorrected - after eta,Et>30,R9,h/e",150,-5.,25.,210,-5.,100.); 
  //Razr variables
  TH1F* h_ggMR = new TH1F("ggMR","Razor MR Variable in events with two photons with Et>40,25 GeV;M_{R} (GeV);Events",300,0.,1500.);
  TH1F* h_egMR = new TH1F("egMR","Razor MR Variable in events with one photon and one electron with Et>40,25 GeV;M_{R} (GeV);Events",300,0.,1500.);
  TH1F* h_eeMR = new TH1F("eeMR","Razor MR Variable in events with two electrons with Et>40,25 GeV;M_{R} (GeV);Events",300,0.,1500.);
  TH1F* h_eeSideBandMR = new TH1F("eeSideBandMR","Razor MR Variable in events with two electrons in sideband with Et>40,25 GeV;M_{R} (GeV);Events",300,0.,1500.);
  TH1F* h_eeMR_reweight_binned = new TH1F("eeMR_reweight_binned","Razor MR Variable in events with two electrons with Et>40,25 GeV;M_{R} (GeV);Events",300,0.,1500.);
  TH1F* h_eeSideBandMR_reweight_binned = new TH1F("eeSideBandMR_reweight_binned","Razor MR Variable in events with two electrons in sideband with Et>40,25 GeV;M_{R} (GeV);Events",300,0.,1500.);
  TH1F* h_ffMR = new TH1F("ffMR","Razor MR Variable in events with two fake photons with Et>40,25 GeV;M_{R} (GeV);Events",300,0.,1500.);
  TH1F* h_ffMR_reweight_binned = new TH1F("ffMR_reweight_binned","Razor MR Variable in events with two fake photons with Et>40,25 GeV;M_{R} (GeV);Events",300,0.,1500.);

  TH1F* h_ggR2 = new TH1F("ggR2","Razor R Variable in events with two photons with Et>40,25 GeV;R^{2};Events",50,0.,1.);
  TH1F* h_egR2 = new TH1F("egR2","Razor R Variable in events with one photon and one electron with Et>40,25 GeV;R^{2};Events",50,0.,1.);
  TH1F* h_eeR2 = new TH1F("eeR2","Razor R Variable in events with two electrons with Et>40,25 GeV;R^{2};Events",50,0.,1.);
  TH1F* h_eeSideBandR2 = new TH1F("eeSideBandR2","Razor R Variable in events with two electrons in sideband with Et>40,25 GeV;R^{2};Events",50,0.,1.);
  TH1F* h_eeR2_reweight_binned = new TH1F("eeR2_reweight_binned","Razor R Variable in events with two electrons with Et>40,25 GeV;R^{2};Events",50,0.,1.);
  TH1F* h_eeSideBandR2_reweight_binned = new TH1F("eeSideBandR2_reweight_binned","Razor R Variable in events with two electrons in sideband with Et>40,25 GeV;R^{2};Events",50,0.,1.);
  TH1F* h_ffR2 = new TH1F("ffR2","Razor R Variable in events with two fake photons with Et>40,25 GeV;R^{2};Events",50,0.,1.);
  TH1F* h_ffR2_reweight_binned = new TH1F("ffR2_reweight_binned","Razor R Variable in events with two fake photons with Et>40,25 GeV;R^{2};Events",50,0.,1.);

  TH2F* h_ggR2vsMR = new TH2F("ggR2vsMR","gg R^{2} vs M_{R}",300,0.,1500.,50,0.,1.);
  TH2F* h_egR2vsMR = new TH2F("egR2vsMR","eg R^{2} vs M_{R}",300,0.,1500.,50,0.,1.);
  TH2F* h_eeR2vsMR = new TH2F("eeR2vsMR","ee R^{2} vs M_{R}",300,0.,1500.,50,0.,1.);
  TH2F* h_eeR2vsMR_reweight_binned = new TH2F("eeR2vsMR_reweight_binned","ee_reweight_binned R^{2} vs M_{R}",300,0.,1500.,50,0.,1.);
  TH2F* h_eeSideBandR2vsMR = new TH2F("eeSideBandR2vsMR","eeSideBand R^{2} vs M_{R}",300,0.,1500.,50,0.,1.);
  TH2F* h_eeSideBandR2vsMR_reweight_binned = new TH2F("eeSideBandR2vsMR_reweight_binned","eeSideBand_reweight_binned R^{2} vs M_{R}",300,0.,1500.,50,0.,1.);
  TH2F* h_ffR2vsMR = new TH2F("ffR2vsMR","ff R^{2} vs M_{R}",300,0.,1500.,50,0.,1.);
  TH2F* h_ffR2vsMR_reweight_binned = new TH2F("ffR2vsMR_reweight_binned","ff_reweight_binned R^{2} vs M_{R}",300,0.,1500.,50,0.,1.);

  TH2F* h_NumE_NumG = new TH2F("NumE_NumG","Number of e(x), and  Number of g(y) from Pho_Cands sample;Number of #gamma;Number of e",7,0,7,7,0,7);
  TH3F* h_NumE_NumG_NumF = new TH3F("NumE_NumG_NumF","Number of e(x), Number of g(y), Number of f(z);Number of #gamma;Number of e; Number of f",7,0,7,7,0,7,7,0,7);
  TH2I* h_NPho_NFake = new TH2I("NPho_NFake","Number of photon (e/g) candidates(x) and Number of fake candidates(y);Number of Pho(e/g) candidates;Number of Fake candidates",5,0,5,5,0,5); 

  TH1I* h_NumGsfEles = new TH1I("NumGsfEles","Number of GsfElectrons that pass quality cuts;;Events",20,0,20);
  TH1I* h_NumPfEles = new TH1I("NumPfEles","Number of pfElectrons that pass quality cuts;;Events",20,0,20);
  TH1I* h_NumGsfEles_NoQcut = new TH1I("NumGsfEles_NoQcut","Number of GsfElectrons;;Events",20,0,20);
  TH1I* h_NumPfEles_NoQcut = new TH1I("NumPfEles_NoQcut","Number of pfElectrons;;Events",20,0,20);;
  TH1I* h_NumMvaEles_NoQcut = new TH1I("NumMvaEles_NoQcut","Number of MvaElectrons;;Events",20,0,20);
  TH1I* h_NumMuons = new TH1I("NumMuons","Number of muons that pass quality cuts;;Events",20,0,20);

  TH1F* h_MET_JetFail = new TH1F("MET_JetFail","Missing transverse energy in events with a L1FastL2L3 corrected PFJet with Pt>50 and #eta<2.6 that fails loose jet id;#slash{E}_{T} (GeV);Events",200,0.,1000.);

  TH1F * h_eMu_InvarMass_FullRange = new TH1F("h_eMu_InvarMass_FullRange","eMu 71<Invariant Mass<111;InvariantMass (GeV);Events", 1000,0.,1000.);
  TH1F * h_eMu_InvarMass = new TH1F("h_eMu_InvarMass","eMu 81<Invariant Mass<101;InvariantMass (GeV);Events", 1000,0.,1000.);
  TH1F * h_MuMu_InvarMass_FullRange = new TH1F("h_MuMu_InvarMass_FullRange","MuMu 71<Invariant Mass<111;InvariantMass (GeV);Events", 1000,0.,1000.);
  TH1F * h_MuMu_InvarMass = new TH1F("h_MuMu_InvarMass","MuMu 81<Invariant Mass<101;InvariantMass (GeV);Events", 1000,0.,1000.);
  TH1F * h_eMu_MET_FullRange = new TH1F("h_eMu_MET_FullRange","eMu 71<Invariant Mass<111;InvariantMass (GeV);Events", 1000,0.,1000.);
  TH1F * h_eMu_MET = new TH1F("h_eMu_MET","eMu 81<Invariant Mass<101;InvariantMass (GeV);Events", 1000,0.,1000.);
  TH1F * h_MuMu_MET_FullRange = new TH1F("h_MuMu_MET_FullRange","MuMu 71<Invariant Mass<111;InvariantMass (GeV);Events", 1000,0.,1000.);
  TH1F * h_MuMu_MET = new TH1F("h_MuMu_MET","MuMu 81<Invariant Mass<101;InvariantMass (GeV);Events", 1000,0.,1000.);

  TH1F* h_MuPt = new TH1F("MuPt","Pt of Muons;P_{T};Events",200,0.,1000.);
  TH1F* h_Mud0 = new TH1F("Mud0","d0 of Muons;d0;Events",100,-20.,20.);
  TH1F* h_MudZ = new TH1F("MudZ","dZ of Muons;dZ;Events",100,-20.,20.);
  TH1F* h_MuCombIso = new TH1F("MuCombIso","CombIso of Muons;CombIso;Events",550,-10.,100.);
  TH1F* h_MuRelIso = new TH1F("MuRelIso","RelIso of Muons;RelIso;Events",100,0.,2.);
  TH1F* h_MuEcalIso = new TH1F("MuEcalIso","EcalIso of Muons;EcalIso;Events",550,-10.,100.);
  TH1F* h_MuHcalIso = new TH1F("MuHcalIso","HcalIso of Muons;HcalIso;Events",550,-10.,100.);
  TH1F* h_MuTrackIso = new TH1F("MuTrackIso","TrackIso of Muons;TrackIso;Events",550,-10.,100.);

  TH1I* MuCount = new TH1I("MuCount","Num of muons passing different cuts;Cut;Events",15,0,15);

  TH1F* h_InstLumi = new TH1F("InstLumi","Number of events for a given instantaneous luminosity;Avg Inst Lumi (Hz/#mub);Events",200,0,2000);  
  TH2F* h_NVertexVsInstLumi = new TH2F("NVertexVsInstLumi","Number of vertices for a given average instantaneous luminosity;Avg Inst Lumi (Hz/#mub);Events",200,0,2000,50,0,50);
  TH1F* h_ggPerInstLumi = new TH1F("ggPerInstLumi","Number of gg events for a given instantaneous luminosity;Avg Inst Lumi (Hz/#mub);Events",200,0,2000);
  TH1F* h_egPerInstLumi = new TH1F("egPerInstLumi","Number of gg events for a given instantaneous luminosity;Avg Inst Lumi (Hz/#mub);Events",200,0,2000);
  TH1F* h_eePerInstLumi = new TH1F("eePerInstLumi","Number of gg events for a given instantaneous luminosity;Avg Inst Lumi (Hz/#mub);Events",200,0,2000);
  TH1F* h_ffPerInstLumi = new TH1F("ffPerInstLumi","Number of gg events for a given instantaneous luminosity;Avg Inst Lumi (Hz/#mub);Events",200,0,2000);

  TH1F* h_ggMet_NV0_10 = new TH1F("ggMet_NV0_10","gg MET, 0<NVertex<10;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_egMet_NV0_10 = new TH1F("egMet_NV0_10","eg MET, 0<NVertex<10;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeMet_NV0_10 = new TH1F("eeMet_NV0_10","ee MET, 0<NVertex<10;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeSBMet_NV0_10 = new TH1F("eeSBMet_NV0_10","ee SideBand MET, 0<NVertex<10;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_ffMet_NV0_10 = new TH1F("ffMet_NV0_10","ff MET, 0<NVertex<10;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  
  TH1F* h_ggMet_NV10_15 = new TH1F("ggMet_NV10_15","gg MET, 10<=NVertex<15;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_egMet_NV10_15 = new TH1F("egMet_NV10_15","eg MET, 10<=NVertex<15;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeMet_NV10_15 = new TH1F("eeMet_NV10_15","ee MET, 10<=NVertex<15;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeSBMet_NV10_15 = new TH1F("eeSBMet_NV10_15","ee SideBand MET, 10<=NVertex<15;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_ffMet_NV10_15 = new TH1F("ffMet_NV10_15","ff MET, 10<=NVertex<15;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  
  TH1F* h_ggMet_NV15up = new TH1F("ggMet_NV15up","gg MET, NVertex>=15;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_egMet_NV15up = new TH1F("egMet_NV15up","eg MET, NVertex>=15;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeMet_NV15up = new TH1F("eeMet_NV15up","ee MET, NVertex>=15;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeSBMet_NV15up = new TH1F("eeSBMet_NV15up","ee SideBand MET, NVertex>=15;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_ffMet_NV15up = new TH1F("ffMet_NV15up","ff MET, NVertex>=15;#slash{E}_{T} (GeV);Events",200,0.,1000.);

  TH1F* h_ggMet_NoJetMatch = new TH1F("ggMet_NoJetMatch","gg MET, at least one has no pf jet matched;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_egMet_NoJetMatch = new TH1F("egMet_NoJetMatch","eg MET, at least one has no pf jet matched;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeMet_NoJetMatch = new TH1F("eeMet_NoJetMatch","ee MET, at least one has no pf jet matched;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_ffMet_NoJetMatch = new TH1F("ffMet_NoJetMatch","ff MET, at least one has no pf jet matched;#slash{E}_{T} (GeV);Events",200,0.,1000.);

  TH1F* h_ggMet_NoJetMatch_JetReq = new TH1F("ggMet_NoJetMatch_JetReq","gg_JetReq MET, at least one has no pf jet matched;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_egMet_NoJetMatch_JetReq = new TH1F("egMet_NoJetMatch_JetReq","eg_JetReq MET, at least one has no pf jet matched;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_eeMet_NoJetMatch_JetReq = new TH1F("eeMet_NoJetMatch_JetReq","ee_JetReq MET, at least one has no pf jet matched;#slash{E}_{T} (GeV);Events",200,0.,1000.);
  TH1F* h_ffMet_NoJetMatch_JetReq = new TH1F("ffMet_NoJetMatch_JetReq","ff_JetReq MET, at least one has no pf jet matched;#slash{E}_{T} (GeV);Events",200,0.,1000.);

  TH1F* h_numTrueInt = new TH1F("numTrueInt","",60,0,60);h_numTrueInt->Sumw2();

  TH1F* h_hadOverEm = new TH1F("hadOverEm","2011 h/e for all photons",50,0,.5);
  TH1F* h_hadTowOverEm = new TH1F("hadTowOverEm","2012 h/e for all photons",50,0,.5);

  if(printLevel > 1) cout<<"Histograms defined" << endl;

  //vectors for Run&Event numbers - print out at end
  vector< pair<int, int> > ggRunEvent,eeRunEvent,ffRunEvent,egRunEvent;
  vector< pair<int, int> > ggJetReqRunEvent,eeJetReqRunEvent,ffJetReqRunEvent,egJetReqRunEvent;
  vector< pair<int, int> > failMetFilterRunEvent;

  //toy diempt plots to get errors 
  TString ggtitle="",eetitle="",fftitle="",eeSBLowtitle="",eeSBHightitle="";
  TString ggtitle_JetReq="",eetitle_JetReq="",fftitle_JetReq="",eeSBLowtitle_JetReq="",eeSBHightitle_JetReq="";
  float Val=0., Err=0.;

  if(printLevel > 1) cout<<"Create Toys directory" << endl;
  TDirectory *toys = fout->mkdir("Toys");
  toys->cd();
  if(printLevel > 1) cout<<"Toys directory created" << endl;

  vector<TH1F*> eeMetValues,ffMetValues;
  vector<TH1F*> eeMet_reweight_binned_toy;
  vector<TH1F*> eeSidebandLowMet_reweight_binned_toy;
  vector<TH1F*> eeSidebandHighMet_reweight_binned_toy;
  vector<TH1F*> ffMet_reweight_binned_toy,ffMet_reweight_binned_toy_JetReq;
  vector<TH1F*> eeMetValues_JetReq,ffMetValues_JetReq;
  vector<TH1F*> eeMet_reweight_binned_toy_JetReq;
  vector<TH1F*> eeSidebandLowMet_reweight_binned_toy_JetReq;
  vector<TH1F*> eeSidebandHighMet_reweight_binned_toy_JetReq;

  vector<TH1F*> h_ffMetChi2,h_ffMetChi2_JetReq;
  vector<TH1F*> h_ffMetChi2PF,h_ffMetChi2PF_JetReq;

  for(int i=0;i<1000;i++){
    eetitle="eeMet_toy_";eetitle+=i+1;
    eeSBLowtitle="eeSidebandLowMet_toy_";eeSBLowtitle+=i+1;
    eeSBHightitle="eeSidebandHighMet_toy_";eeSBHightitle+=i+1;
    fftitle="ffMet_toy_";fftitle+=i+1;

    eeMet_reweight_binned_toy.push_back( new TH1F(eetitle,eetitle,200,0.,1000.));
    eeSidebandLowMet_reweight_binned_toy.push_back( new TH1F(eeSBLowtitle,eeSBLowtitle,200,0.,1000.));
    eeSidebandHighMet_reweight_binned_toy.push_back( new TH1F(eeSBHightitle,eeSBHightitle,200,0.,1000.));
    ffMet_reweight_binned_toy.push_back( new TH1F(fftitle,fftitle,200,0.,1000.));

    eetitle_JetReq="ee_JetReqMet_toy_";eetitle_JetReq+=i+1;
    eeSBLowtitle_JetReq="eeSideband_JetReqLowMet_toy_";eeSBLowtitle_JetReq+=i+1;
    eeSBHightitle_JetReq="eeSideband_JetReqHighMet_toy_";eeSBHightitle_JetReq+=i+1;
    fftitle_JetReq="ff_JetReqMet_toy_";fftitle_JetReq+=i+1;

    eeMet_reweight_binned_toy_JetReq.push_back( new TH1F(eetitle_JetReq,eetitle_JetReq,200,0.,1000.));
    eeSidebandLowMet_reweight_binned_toy_JetReq.push_back( new TH1F(eeSBLowtitle_JetReq,eeSBLowtitle_JetReq,200,0.,1000.));
    eeSidebandHighMet_reweight_binned_toy_JetReq.push_back( new TH1F(eeSBHightitle_JetReq,eeSBHightitle_JetReq,200,0.,1000.));
    ffMet_reweight_binned_toy_JetReq.push_back( new TH1F(fftitle_JetReq,fftitle_JetReq,200,0.,1000.));


    if(i<51){
      eetitle="ffMet_Cut_";eetitle+=i+1;eetitle+="GeV";
      h_ffMetChi2.push_back( new TH1F(eetitle,eetitle,200,0.,1000.));
      eetitle="ffMet_Cut_";eetitle+=i+1;eetitle+="GeV_JetReq";
      h_ffMetChi2_JetReq.push_back( new TH1F(eetitle,eetitle,200,0.,1000.));

      eetitle="ffMet_Cut_";eetitle+=i+1;eetitle+="GeV_PF";
      h_ffMetChi2PF.push_back( new TH1F(eetitle,eetitle,200,0.,1000.));
      eetitle="ffMet_Cut_";eetitle+=i+1;eetitle+="GeV_PF_JetReq";
      h_ffMetChi2PF_JetReq.push_back( new TH1F(eetitle,eetitle,200,0.,1000.));

    }

  }
  fout->cd();

  if(printLevel > 1) cout<<"Define job-wide variables " << endl;
  //two susy::photon pointers, later they will be used for the top two photon candidates
  susy::Photon* PhoOne = new susy::Photon; PhoOne->Init();
  susy::Photon* PhoTwo = new susy::Photon; PhoTwo->Init();
  susy::PFJet* DiJetPtJet1 = new susy::PFJet; DiJetPtJet1->Init();
  susy::PFJet* DiJetPtJet2 = new susy::PFJet; DiJetPtJet2->Init();
  std::pair<float,float> ggPTs=make_pair(0.,0.),eePTs=make_pair(0.,0.),egPTs=make_pair(0.,0.);
  std::pair<float,float> ggPTs_JetReq=make_pair(0.,0.),eePTs_JetReq=make_pair(0.,0.),egPTs_JetReq=make_pair(0.,0.);

  //susy::Photon* PhoThree = new susy::Photon; PhoThree->Init();
  //counters for the number of events
  int ngg=0, ngge=0, nggee=0, neg=0, nee=0, neef=0, neeg=0, neegg=0, neeSideBand=0, nff=0, ndiPho=0, ndiPhoCands=0, neMu=0, nMuMu=0;
  int ngg_JetReq=0, neg_JetReq=0, nee_JetReq=0, neeSideBand_JetReq=0, nee_2JetReq=0, nff_JetReq=0, ndiPho_JetReq=0, ndiPhoCands_JetReq=0;
  int pfPho=0,notPfPho=0;
  //int nfgg=0, nfeg=0, nfee=0;
  //int nfg=0, nfe=0;
  //counters for random things I want to keep track of
  int nTwoPhosAndTwoFakes=0;
  int nThreePhoEvents=0, nThreeFakeEvents=0;
  int nFourPhoEvents=0, nFourFakeEvents=0;
  int nPhosFailDR=0, nFakesFailDR=0;
  int nPhosFailDphi=0, nFakesFailDphi=0;
  int eleCount=0,pfEleCount=0,MvaEleCount=0;
  int oldP=0,newP=0;
  TRandom rando;
  if(printLevel > 1) cout<<"Job-wide variables defined" << endl;

  //define bin-wise MET corrections from DiEMPt ratios
  //TFile f_weights("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_Data2011_Filter_JsonHLTnVertexTwo40-25GeVBarrelPhosWithR9-HoverE-CaloIdL-R9idOrIsoVL_RhoPileupCorr_Photon_Analysis_OldRatios_NewJetDefWithMuons.root","READ");
  if(printLevel>1)cout<<"Read in weights histogram" << endl;
  TFile f_weights("reweights.root","READ");
  if(printLevel>1)cout<<"Weights histogram read in" << endl;
  f_weights.cd();
  TString eeWeightTitle="",ffWeightTitle="",eeWeightTitleJet="",ffWeightTitleJet="",eeSideBandLowWeightTitleJet="",eeSideBandHighWeightTitleJet="",eeSideBandLowWeightTitle="",eeSideBandHighWeightTitle="";  
  TString eeWeightTitle_JetReq="",ffWeightTitle_JetReq="",eeWeightTitleJet_JetReq="",ffWeightTitleJet_JetReq="",eeSideBandLowWeightTitleJet_JetReq="",eeSideBandHighWeightTitleJet_JetReq="",eeSideBandLowWeightTitle_JetReq="",eeSideBandHighWeightTitle_JetReq="";
  //TH1F *h_weightsEE = (TH1F*)f_weights.Get("ggeeDiEMPtRatio");
  //TH1F *h_weightsFF = (TH1F*)f_weights.Get("ggffDiEMPtRatio");
  vector< pair<float,float> > EEreweights[4];//full,0jet,1jet,>=2jet  filled by bins of diempt
  vector< pair<float,float> > FFreweights[4];
  vector< pair<float,float> >  EESideBandLowReweights[4];
  vector< pair<float,float> >  EESideBandHighReweights[4];
  vector< pair<float,float> >  EESideBandLowReweightsJet[4];
  vector< pair<float,float> >  EESideBandHighReweightsJet[4];
  vector< pair<float,float> >  EEreweightsJet[4];
  vector< pair<float,float> >  FFreweightsJet[4];
  vector<TH1F*> h_weightsEEVec;
  vector<TH1F*> h_weightsFFVec;
  vector<TH1F*> h_weightsEEJetVec;
  vector<TH1F*> h_weightsFFJetVec;
  vector<TH1F*> h_weightsEEsidebandLowVec;
  vector<TH1F*> h_weightsEEsidebandHighVec;
  vector<TH1F*> h_weightsEEsidebandLowJetVec;
  vector<TH1F*> h_weightsEEsidebandHighJetVec;
  vector< pair<float,float> > EEreweights_JetReq[4];//full,0jet,1jet,>=2jet  filled by bins of diempt
  vector< pair<float,float> > FFreweights_JetReq[4];
  vector< pair<float,float> >  EEreweightsJet_JetReq[4];
  vector< pair<float,float> >  FFreweightsJet_JetReq[4];
  vector< pair<float,float> >  EESideBandLowReweights_JetReq[4];
  vector< pair<float,float> >  EESideBandHighReweights_JetReq[4];
  vector< pair<float,float> >  EESideBandLowReweightsJet_JetReq[4];
  vector< pair<float,float> >  EESideBandHighReweightsJet_JetReq[4];
  vector<TH1F*> h_weightsEEVec_JetReq;
  vector<TH1F*> h_weightsEEJetVec_JetReq;
  vector<TH1F*> h_weightsFFVec_JetReq;
  vector<TH1F*> h_weightsFFJetVec_JetReq;
  vector<TH1F*> h_weightsEEsidebandLowVec_JetReq;
  vector<TH1F*> h_weightsEEsidebandHighVec_JetReq;
  vector<TH1F*> h_weightsEEsidebandLowJetVec_JetReq;
  vector<TH1F*> h_weightsEEsidebandHighJetVec_JetReq;
  if(printLevel>5)cout<<"before grabbing weights histograms"<<endl;
  for(int i=0;i<4;i++){
    eeWeightTitle="ggeeDiEMPtRatio";ffWeightTitle="ggffDiEMPtRatio";
    eeWeightTitleJet="ggeeDiJetPtRatio";ffWeightTitleJet="ggffDiJetPtRatio";
    eeSideBandLowWeightTitle="ggeeSideBandLowDiEMPtRatio";eeSideBandHighWeightTitle="ggeeSideBandHighDiEMPtRatio";
    eeSideBandLowWeightTitleJet="ggeeSideBandLowDiJetPtRatio";eeSideBandHighWeightTitleJet="ggeeSideBandHighDiJetPtRatio";
    eeWeightTitle_JetReq="ggeeDiEMPtRatio_JetReq";ffWeightTitle_JetReq="ggffDiEMPtRatio_JetReq";
    eeWeightTitleJet_JetReq="ggeeDiJetPtRatio_JetReq";ffWeightTitleJet_JetReq="ggffDiJetPtRatio_JetReq";
    eeSideBandLowWeightTitle_JetReq="ggeeSideBandLowDiEMPtRatio_JetReq";eeSideBandHighWeightTitle_JetReq="ggeeSideBandHighDiEMPtRatio_JetReq";
    eeSideBandLowWeightTitleJet_JetReq="ggeeSideBandLowDiJetPtRatio_JetReq";eeSideBandHighWeightTitleJet_JetReq="ggeeSideBandHighDiJetPtRatio_JetReq";
    if(i==0){
      h_weightsEEVec.push_back((TH1F*)f_weights.Get(eeWeightTitle));
      h_weightsFFVec.push_back((TH1F*)f_weights.Get(ffWeightTitle));
      h_weightsEEJetVec.push_back((TH1F*)f_weights.Get(eeWeightTitleJet));
      h_weightsFFJetVec.push_back((TH1F*)f_weights.Get(ffWeightTitleJet));
      h_weightsEEsidebandLowVec.push_back((TH1F*)f_weights.Get(eeSideBandLowWeightTitle));
      h_weightsEEsidebandHighVec.push_back((TH1F*)f_weights.Get(eeSideBandHighWeightTitle));
      h_weightsEEsidebandLowJetVec.push_back((TH1F*)f_weights.Get(eeSideBandLowWeightTitleJet));
      h_weightsEEsidebandHighJetVec.push_back((TH1F*)f_weights.Get(eeSideBandHighWeightTitleJet));
      h_weightsEEVec_JetReq.push_back((TH1F*)f_weights.Get(eeWeightTitle_JetReq));
      h_weightsFFVec_JetReq.push_back((TH1F*)f_weights.Get(ffWeightTitle_JetReq));
      h_weightsEEJetVec_JetReq.push_back((TH1F*)f_weights.Get(eeWeightTitleJet_JetReq));
      h_weightsFFJetVec_JetReq.push_back((TH1F*)f_weights.Get(ffWeightTitleJet_JetReq));
      h_weightsEEsidebandLowVec_JetReq.push_back((TH1F*)f_weights.Get(eeSideBandLowWeightTitle_JetReq));
      h_weightsEEsidebandHighVec_JetReq.push_back((TH1F*)f_weights.Get(eeSideBandHighWeightTitle_JetReq));
      h_weightsEEsidebandLowJetVec_JetReq.push_back((TH1F*)f_weights.Get(eeSideBandLowWeightTitleJet_JetReq));
      h_weightsEEsidebandHighJetVec_JetReq.push_back((TH1F*)f_weights.Get(eeSideBandHighWeightTitleJet_JetReq));
    }
    else{
      //first diEMPt
      eeWeightTitle+="_";eeWeightTitle+=i-1;eeWeightTitle+="Jet";//i-1 to get 0Jet,1Jet,2Jet instead of 1Jet, etc.  (0th element is full dist.)
      ffWeightTitle+="_";ffWeightTitle+=i-1;ffWeightTitle+="Jet";
      if(printLevel>5)cout<<"eeWeightTitle="<<eeWeightTitle<<endl;
      h_weightsEEVec.push_back((TH1F*)f_weights.Get(eeWeightTitle));
      h_weightsFFVec.push_back((TH1F*)f_weights.Get(ffWeightTitle));
      //now diJetPt
      eeWeightTitleJet+="_";eeWeightTitleJet+=i-1;eeWeightTitleJet+="Jet";
      ffWeightTitleJet+="_";ffWeightTitleJet+=i-1;ffWeightTitleJet+="Jet";
      eeSideBandLowWeightTitleJet+="_";eeSideBandLowWeightTitleJet+=i-1;eeSideBandLowWeightTitleJet+="Jet";
      eeSideBandHighWeightTitleJet+="_";eeSideBandHighWeightTitleJet+=i-1;eeSideBandHighWeightTitleJet+="Jet";
      eeSideBandLowWeightTitle+="_";eeSideBandLowWeightTitle+=i-1;eeSideBandLowWeightTitle+="Jet";
      eeSideBandHighWeightTitle+="_";eeSideBandHighWeightTitle+=i-1;eeSideBandHighWeightTitle+="Jet";
      if(printLevel>5)cout<<"eeWeightTitleJet="<<eeWeightTitleJet<<endl;
      h_weightsEEJetVec.push_back((TH1F*)f_weights.Get(eeWeightTitleJet));
      h_weightsFFJetVec.push_back((TH1F*)f_weights.Get(ffWeightTitleJet));
      h_weightsEEsidebandLowVec.push_back((TH1F*)f_weights.Get(eeSideBandLowWeightTitle));
      h_weightsEEsidebandHighVec.push_back((TH1F*)f_weights.Get(eeSideBandHighWeightTitle));
      h_weightsEEsidebandLowJetVec.push_back((TH1F*)f_weights.Get(eeSideBandLowWeightTitleJet));
      h_weightsEEsidebandHighJetVec.push_back((TH1F*)f_weights.Get(eeSideBandHighWeightTitleJet));
      //first diEMPt
      eeWeightTitle_JetReq+="_";eeWeightTitle_JetReq+=i-1;eeWeightTitle_JetReq+="Jet";//i-1 to get 0Jet,1Jet,2Jet instead of 1Jet, etc.  (0th element is full dist.)
      ffWeightTitle_JetReq+="_";ffWeightTitle_JetReq+=i-1;ffWeightTitle_JetReq+="Jet";
      if(printLevel>5)cout<<"eeWeightTitle_JetReq="<<eeWeightTitle_JetReq<<endl;
      h_weightsEEVec_JetReq.push_back((TH1F*)f_weights.Get(eeWeightTitle_JetReq));
      h_weightsFFVec_JetReq.push_back((TH1F*)f_weights.Get(ffWeightTitle_JetReq));
      //now diJetPt
      eeWeightTitleJet_JetReq+="_";eeWeightTitleJet_JetReq+=i-1;eeWeightTitleJet_JetReq+="Jet";
      ffWeightTitleJet_JetReq+="_";ffWeightTitleJet_JetReq+=i-1;ffWeightTitleJet_JetReq+="Jet";
      eeSideBandLowWeightTitle_JetReq+="_";eeSideBandLowWeightTitle_JetReq+=i-1;eeSideBandLowWeightTitle_JetReq+="Jet";
      eeSideBandHighWeightTitle_JetReq+="_";eeSideBandHighWeightTitle_JetReq+=i-1;eeSideBandHighWeightTitle_JetReq+="Jet";
      eeSideBandLowWeightTitleJet_JetReq+="_";eeSideBandLowWeightTitleJet_JetReq+=i-1;eeSideBandLowWeightTitleJet_JetReq+="Jet";
      eeSideBandHighWeightTitleJet_JetReq+="_";eeSideBandHighWeightTitleJet_JetReq+=i-1;eeSideBandHighWeightTitleJet_JetReq+="Jet";
      if(printLevel>5)cout<<"eeWeightTitleJet_JetReq="<<eeWeightTitleJet_JetReq<<endl;
      h_weightsEEJetVec_JetReq.push_back((TH1F*)f_weights.Get(eeWeightTitleJet_JetReq));
      h_weightsFFJetVec_JetReq.push_back((TH1F*)f_weights.Get(ffWeightTitleJet_JetReq));
      h_weightsEEsidebandLowVec_JetReq.push_back((TH1F*)f_weights.Get(eeSideBandLowWeightTitle_JetReq));
      h_weightsEEsidebandHighVec_JetReq.push_back((TH1F*)f_weights.Get(eeSideBandHighWeightTitle_JetReq));
      h_weightsEEsidebandLowJetVec_JetReq.push_back((TH1F*)f_weights.Get(eeSideBandLowWeightTitleJet_JetReq));
      h_weightsEEsidebandHighJetVec_JetReq.push_back((TH1F*)f_weights.Get(eeSideBandHighWeightTitleJet_JetReq));
    }
  }
  if(printLevel>5)cout<<"after grabbing weights histograms"<<endl;
 
  for(int j=0;j<4;j++)
    {
      EEreweights[j].clear();FFreweights[j].clear();EEreweightsJet[j].clear();FFreweightsJet[j].clear();
      EESideBandLowReweightsJet[j].clear();EESideBandHighReweightsJet[j].clear();EESideBandLowReweights[j].clear();EESideBandHighReweights[j].clear();
      EEreweights_JetReq[j].clear();FFreweights_JetReq[j].clear();EEreweightsJet_JetReq[j].clear();FFreweightsJet_JetReq[j].clear();
      EESideBandLowReweightsJet_JetReq[j].clear();EESideBandHighReweightsJet_JetReq[j].clear();EESideBandLowReweights_JetReq[j].clear();EESideBandHighReweights_JetReq[j].clear();
      if(printLevel>5)cout<<"Inside filling weight vector, j= "<<j<<endl;
      for(int i=0;i<16;i++){
	//0-th index should be bin 1
	EEreweights[j].push_back(make_pair(h_weightsEEVec[j]->GetBinContent(i+1),h_weightsEEVec[j]->GetBinError(i+1)));
	FFreweights[j].push_back(make_pair(h_weightsFFVec[j]->GetBinContent(i+1),h_weightsFFVec[j]->GetBinError(i+1)));
	
	EEreweightsJet[j].push_back(make_pair(h_weightsEEJetVec[j]->GetBinContent(i+1),h_weightsEEJetVec[j]->GetBinError(i+1)));
	FFreweightsJet[j].push_back(make_pair(h_weightsFFJetVec[j]->GetBinContent(i+1),h_weightsFFJetVec[j]->GetBinError(i+1)));
	EESideBandLowReweights[j].push_back(make_pair(h_weightsEEsidebandLowVec[j]->GetBinContent(i+1),h_weightsEEsidebandLowVec[j]->GetBinError(i+1)));
	EESideBandHighReweights[j].push_back(make_pair(h_weightsEEsidebandHighVec[j]->GetBinContent(i+1),h_weightsEEsidebandHighVec[j]->GetBinError(i+1)));
	EESideBandLowReweightsJet[j].push_back(make_pair(h_weightsEEsidebandLowJetVec[j]->GetBinContent(i+1),h_weightsEEsidebandLowJetVec[j]->GetBinError(i+1)));
	EESideBandHighReweightsJet[j].push_back(make_pair(h_weightsEEsidebandHighJetVec[j]->GetBinContent(i+1),h_weightsEEsidebandHighJetVec[j]->GetBinError(i+1)));

	EEreweights_JetReq[j].push_back(make_pair(h_weightsEEVec_JetReq[j]->GetBinContent(i+1),h_weightsEEVec_JetReq[j]->GetBinError(i+1)));
	FFreweights_JetReq[j].push_back(make_pair(h_weightsFFVec_JetReq[j]->GetBinContent(i+1),h_weightsFFVec_JetReq[j]->GetBinError(i+1)));
	
	EEreweightsJet_JetReq[j].push_back(make_pair(h_weightsEEJetVec_JetReq[j]->GetBinContent(i+1),h_weightsEEJetVec_JetReq[j]->GetBinError(i+1)));
	FFreweightsJet_JetReq[j].push_back(make_pair(h_weightsFFJetVec_JetReq[j]->GetBinContent(i+1),h_weightsFFJetVec_JetReq[j]->GetBinError(i+1)));
	EESideBandLowReweights_JetReq[j].push_back(make_pair(h_weightsEEsidebandLowVec_JetReq[j]->GetBinContent(i+1),h_weightsEEsidebandLowVec_JetReq[j]->GetBinError(i+1)));
	EESideBandHighReweights_JetReq[j].push_back(make_pair(h_weightsEEsidebandHighVec_JetReq[j]->GetBinContent(i+1),h_weightsEEsidebandHighVec_JetReq[j]->GetBinError(i+1)));
	EESideBandLowReweightsJet_JetReq[j].push_back(make_pair(h_weightsEEsidebandLowJetVec_JetReq[j]->GetBinContent(i+1),h_weightsEEsidebandLowJetVec_JetReq[j]->GetBinError(i+1)));
	EESideBandHighReweightsJet_JetReq[j].push_back(make_pair(h_weightsEEsidebandHighJetVec_JetReq[j]->GetBinContent(i+1),h_weightsEEsidebandHighJetVec_JetReq[j]->GetBinError(i+1)));
      }
    }
  
  f_weights.Close();
  fout->cd();

  TH1F* puweights;
  //Get MC pileup weight
  //TFile f_PUweights("PUweightsForBino.root","READ");
  TFile f_PUweights("PUweightsVBF.root","READ");
  f_PUweights.cd();
  //puweights = (TH1F*)f_PUweights.Get("PUweightsBino");
  puweights = (TH1F*)f_PUweights.Get("PUweights");
  //f_PUweights.Close();
  fout->cd();
  float PUweight=1.;


  //define NVertex and Rho Correction factors
  if(doRhoCorrection)cout<<"Applying Rho Pileup corrections!"<<endl;
  else cout<<"Applying NO Pileup corrections!"<<endl;
  
  // to check duplicated events
  std::map<int, std::set<int> > allEvents;

  // start event looping
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < processNEvents; jentry++) {

    //make sure we are not trying to do multiple corrections at once
    if(doNVertexCorrection && doRhoCorrection){
      cout<<"Trying to correct for both NVertex and Rho!!!!!"<<endl;
      break;
    }

    //bool filterThis = false;//want to keep subset of events

    if(printLevel > 1) cout << endl << "Get the tree contents." << endl;

    // if( (int)(((float)jentry/processNEvents)*100)%10==0 )cout << (int)(((float)jentry/processNEvents)*100) << " percent done..." << endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    if(printLevel > 1 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0)) ) {
      cout <</* endl <<*/ int(jentry) << " events processed with Run:"
	   << event->runNumber << ", Event:" << event->eventNumber << endl;
    }

    //cout<<event->runNumber<<"  "<<event->luminosityBlockNumber<<"  "<<event->avgInsRecLumi<<"  "<<event->intgRecLumi<<endl;
    //if(printLevel > 1) cout << "Initialize any global variables to be reset per event." << endl;

    //InitializePerEvent();

    if(!event->isRealData &&jentry==0){
      cout<<"Neutralino Mass: " <<event->gridParams["mChi0"]<<endl;
      cout<<"Gluino Mass: " <<event->gridParams["mGluino"]<<endl;
      cout<<"Squark Mass: " <<event->gridParams["mSquark"]<<endl;
      cout<<"Cross Section: " <<event->gridParams["xsec"]<<endl;
      //cout<<"ptHat: " <<event->gridParams["ptHat"]<<endl;
    }

    float Rho = event->rho, RhoB = event->rhoBarrel, Rho25 = event->rho25;

    if(event->runNumber>=0/* && event->runNumber<=175831*/){
      
      //check from gg events that I don't have 
      /*
	if(event->isRealData){
	if(
	(event->runNumber==176797 && event->eventNumber==282232912)   ||
	(event->runNumber==172949 && event->eventNumber==367225623) 
	){
	  
	cout<<"Run: "<<event->runNumber<<"  Event: "<<event->eventNumber<<endl;
	if(!isInJson(event->runNumber,event->luminosityBlockNumber))cout<<"Not In JSON!"<<endl;
	else{ cout<<"Is In Json"<<endl;
	    
	bool passHLT_Test =  (useTrigger ? PassTriggers() : true);
	if(!passHLT_Test)cout<<"Does NOT Pass HLT!"<<endl;
	else{ cout<<"Passes HLT!"<<endl;
	      
	int NVertex=0;
	for(std::vector<susy::Vertex>::iterator Vtx_it = event->vertices.begin(); Vtx_it<event->vertices.end(); Vtx_it++){
	if(    !Vtx_it->isFake() 
	&& Vtx_it->ndof>4 
	&& fabs(Vtx_it->position.z()<24.0) 
	&& sqrt(Vtx_it->position.x()*Vtx_it->position.x()+Vtx_it->position.y()*Vtx_it->position.y())<2.0 ) NVertex++;
	}
	if(NVertex>=1)cout<<"At least one good vertex!"<<endl;
	else cout<<"-----No good vertex!"<<endl;
	std::map<TString, std::vector<susy::Photon> >::iterator phoMap = event->photons.find("photons");
	if(phoMap == event->photons.end()) {
	if(event->photons.size() > 0) cout << "photon collection is not available!" << endl;
	}
	else {
		
	//cout<<"pho_Cands size: "<<pho_Cands.size()<<endl;
	//cout<<"fake_Cands size: "<<fake_Cands.size()<<endl;
		
	for(std::vector<susy::Photon>::iterator it_Pho_Test = phoMap->second.begin(); it_Pho_Test!= phoMap->second.end(); it_Pho_Test++) {
	if(it_Pho_Test == phoMap->second.begin()){
	cout << "Number of photons: "<<phoMap->second.size()<<endl;
	}
	cout<<"--------Photon Et: "<<it_Pho_Test->momentum.Et()<<"  Eta: "<<it_Pho_Test->caloPosition.Eta()<<"  Phi: "<<it_Pho_Test->caloPosition.Phi()<<endl;
	if(it_Pho_Test->isEB()&& !it_Pho_Test->isEBEtaGap() 
	&& !it_Pho_Test->isEBPhiGap() && !it_Pho_Test->isEERingGap() 
	&& !it_Pho_Test->isEEDeeGap() && !it_Pho_Test->isEBEEGap() )cout<<"Is EB Fed!"<<endl;
	else{ cout<<"Is NOT EB Fed!"<<endl;
	if(!it_Pho_Test->isEB()) cout<<"Is NOT EB!"<<endl;
	if(it_Pho_Test->isEBEtaGap()) cout<<"Is isEBEtaGap!"<<endl;
	if(it_Pho_Test->isEBPhiGap()) cout<<"Is isEBPhiGap!"<<endl;
	if(it_Pho_Test->isEERingGap()) cout<<"Is isEERingGap!"<<endl;
	if(it_Pho_Test->isEEDeeGap()) cout<<"Is isEEDeeGap!"<<endl;
	if(it_Pho_Test->isEBEEGap()) cout<<"Is isEBEEGap!"<<endl;
	}
	cout<<"Ecal Iso: "<<it_Pho_Test->ecalRecHitSumEtConeDR04<<" Cut: "<<4.2+0.006*it_Pho_Test->momentum.Et()<<endl;
	cout<<"Hcal Iso: "<<it_Pho_Test->hcalTowerSumEtConeDR04()<<" Cut: "<<2.2+0.0025*it_Pho_Test->momentum.Et()<<endl;
	cout<<"Track Iso: "<<it_Pho_Test->trkSumPtHollowConeDR04<<" Cut: "<<2.0+0.001*it_Pho_Test->momentum.Et()<<endl;
	cout<<"Comb Iso: "<<it_Pho_Test->ecalRecHitSumEtConeDR03 - PUCorr_ECAL*Rho25 + it_Pho_Test->hcalTowerSumEtConeDR03() - PUCorr_HCAL*Rho25 + it_Pho_Test->trkSumPtHollowConeDR03 - PUCorr_TRACK*Rho25<< " Cut: "<<6.0<<endl;
	cout<<"SigmaIetaIeta: "<<it_Pho_Test->sigmaIetaIeta<<endl;
	cout<<"H/E: "<<it_Pho_Test->hadronicOverEm<<endl;
	cout<<"Number of Pixel Seeds: "<<it_Pho_Test->nPixelSeeds<<endl;
		  
	}
	}
	}
	}
	Print(*event);
	}
	} //End of test 
      */ 
      
      
      nCnt[0]++; // total number of events   
      
      if(printLevel > 1) cout << "Apply good run list." << endl;
      if(printLevel > 1) cout<<"runNumber="<<event->runNumber<<"  lumiNumber="<<event->luminosityBlockNumber<<endl;
      if(useJSON){
	if(!isInJson(event->runNumber,event->luminosityBlockNumber)) continue;
      }

      nCnt[1]++; // total number of events that pass Json

      // uncomment this to print all ntuple variables
      //Print(*event);

      if(printLevel > 1) cout << "Check duplicated events for data only." << endl;
     
      if(event->isRealData){
	bool duplicateEvent = ! (allEvents[event->runNumber].insert(event->eventNumber)).second;
	if(duplicateEvent){
	  cout<<"Duplicate Event! Run "<<event->runNumber<<" Event "<<event->eventNumber<<endl;
	  continue;
	}
      }
      nCnt[2]++;//number of events that pass duplicate check
    
      if(printLevel > 1) cout << "Apply trigger selection in the event." << endl;
      bool passHLT = (useTrigger ? PassTriggers() : true);
      if(!passHLT )continue;//only accept events which pass our trigger(s)
      nCnt[3]++;// number of events that pass trigger

      if(printLevel > 1) cout << "Find primary vertex in the event." << endl;
      TVector3* primVtx = 0;
      if(event->vertices.size() > 0) primVtx = &(event->vertices[0].position);
      if(primVtx) h_vtxZ->Fill(primVtx->Z());
      h_bsZ->Fill(event->beamSpot.Z());

      //  Get NVertex and Rho for event
      int NVertex=0;
      for(std::vector<susy::Vertex>::iterator Vtx_it = event->vertices.begin(); Vtx_it<event->vertices.end(); Vtx_it++){
	if(    !Vtx_it->isFake() 
	       && Vtx_it->ndof>4 
	       && fabs(Vtx_it->position.z()<24.0) 
	       && sqrt(Vtx_it->position.x()*Vtx_it->position.x()+Vtx_it->position.y()*Vtx_it->position.y())<2.0 ) NVertex++;
      }
      h_rho->Fill(Rho);h_NVertex->Fill(NVertex);
      h_rhoBarrel->Fill(RhoB);

      //Require at least 1 good vertex
      if(NVertex<1){
	if(printLevel > 1){cout<<"No Good Vertex!!!!  Run: "<<event->runNumber<<"  Event: "<<event->eventNumber<<endl;}
	continue;
      }
      nCnt[4]++;// number of events that pass nVertex>=1

      int numTrueInt = -1.;
      if(!event->isRealData){
	susy::PUSummaryInfoCollection::const_iterator iBX = event->pu.begin();
	bool foundInTimeBX = false;
	while((iBX != event->pu.end()) && !foundInTimeBX) {
	  if(iBX->BX == 0) {
	    numTrueInt = iBX->trueNumInteractions;
	    foundInTimeBX = true;
	  }
	  ++iBX;
	}
      }
      h_numTrueInt->Fill(numTrueInt);

      if(printLevel>1)cout<<"Line 1387, numTrueInt="<<numTrueInt<<endl;

      int PUbin=0;
      if(!event->isRealData)PUbin = puweights->FindBin(double(numTrueInt));
      if(printLevel>1)cout<<"Line 1390, PUbin="<<PUbin<<endl;

      if(!event->isRealData)PUweight = puweights->GetBinContent(PUbin);

      if(printLevel>1)cout<<"Line 1394, PUweight="<<PUweight<<endl;



      if(event->isRealData && !event->passMetFilters()){failMetFilterRunEvent.push_back(make_pair(event->runNumber,event->eventNumber));continue;}
      nCnt[5]++;// number of events that pass all Met filters //HcalNoiseFilter

      if(printLevel > 1) cout << "Select which met will be used in the event." << endl;
      std::map<TString, susy::MET>::iterator met_it = event->metMap.find("pfMet");
      if(met_it == event->metMap.end()) {
	cout << "MET map is not available!!!" << endl;
	continue;
      }
      susy::MET* met = &(met_it->second);
      h_met->Fill(met->met());
      h_sumEt->Fill(met->sumEt);


      h_InstLumi->Fill(event->avgInsRecLumi);
      h_NVertexVsInstLumi->Fill(event->avgInsRecLumi,NVertex);
      if(printLevel > 1) cout << "Find PFelectrons in the event." << endl;
      //----------
      //find PFelectrons, apply quality cuts
      //----------
      std::vector<susy::Electron*>   pfEles;
      int nGsf=0,nPF=0,nGsf_NoQcut=0,nPF_NoQcut=0,nMva_NoQcut=0;
      std::map<TString, std::vector<susy::Electron> >::iterator eleMap = event->electrons.find("gsfElectrons");
      if(eleMap != event->electrons.end()) {
	nGsf=0;nPF=0;nGsf_NoQcut=0;nPF_NoQcut=0;nMva_NoQcut=0;
	//loop over electron collection 
	for(std::vector<susy::Electron>::iterator it_Ele = eleMap->second.begin(); it_Ele != eleMap->second.end(); it_Ele++) {
	  if(printLevel > 1) cout<<"looping over electron collection"<<endl;
	  eleCount++;nGsf_NoQcut++;
	  if(it_Ele->isPF()){ pfEleCount++;nPF_NoQcut++;}//continue;
	  if(it_Ele->passingMvaPreselection()){ MvaEleCount++;nMva_NoQcut++;}
	  else continue;
	  float Iso=0.;
	  /*if(it_Ele->isPF())*/Iso=it_Ele->chargedHadronIso + it_Ele->neutralHadronIso + it_Ele->photonIso;
	  //else Iso=it_Ele->dr03EcalRecHitSumEt + it_Ele->dr03HcalDepth1TowerSumEt + it_Ele->dr03HcalDepth2TowerSumEt + it_Ele->dr03TkSumPt;
	  float pt=it_Ele->momentum.Pt();
	  
	  if(it_Ele->momentum.Eta()>2.6) continue;
	  if(pt<15) continue;
	  if(Iso/pt>0.2) continue;
	  
	  //nGsf++;
	  
	  pfEles.push_back(&*it_Ele);
	}//end it_Ele electron loop
	//sort pfEles by Pt
	std::sort(pfEles.begin(), pfEles.end(), EtGreater<susy::Electron>);
	if(printLevel>1)cout<<"pfEles size= "<<pfEles.size()<<endl;
	h_NumGsfEles_NoQcut->Fill(nGsf_NoQcut);
	h_NumPfEles_NoQcut->Fill(nPF_NoQcut);
	h_NumMvaEles_NoQcut->Fill(nMva_NoQcut);
	h_NumGsfEles->Fill(nGsf);
	h_NumPfEles->Fill(nPF);
      }//end eleMap

      if(printLevel > 1) cout << "Find Muons in the event." << endl;
      //----------
      //find Muons, apply quality cuts
      //----------
      std::vector<susy::Muon*>   Muons;
      //loop over muon collection 
      //for(std::vector<susy::Muon>::iterator it_Mu = muMap->second.begin(); it_Mu != muMap->second.end(); it_Mu++) {
      for(std::vector<susy::Muon>::iterator it_Mu = event->muons.begin();
	  it_Mu != event->muons.end(); it_Mu++) {
	if(printLevel > 1) cout<<"looping over muon collection"<<endl;
	
	MuCount->Fill(0);
	//std::map<TString,UChar_t>::iterator idPairMu = it_Mu->idPairs.find("muidGlobalMuonPromptTight");
	//if(idPairMu == it_Mu->idPairs.end()) continue;MuCount->Fill(1);
	//if(int(idPairMu->second) == 0) continue;MuCount->Fill(2);
	
	if(!it_Mu->isGlobalMuon()) continue;MuCount->Fill(1);
	if(!it_Mu->isPFMuon()) continue;MuCount->Fill(2);

	float pt = it_Mu->momentum.Pt();
	float combIso=(it_Mu->ecalIsoR03 + it_Mu->hcalIsoR03 + it_Mu->trackIsoR03);
	float combIsoPF=( it_Mu->sumChargedHadronPt04 + std::max(0.,it_Mu->sumNeutralHadronEt04+it_Mu->sumPhotonEt04-0.5*it_Mu->sumPUPt04) );
	//float combIsoPF=( it_Mu->sumChargedHadronPt04 + it_Mu->sumNeutralHadronEt04 + it_Mu->sumPhotonEt04 );
	float relIso = combIso/pt;
	float relIsoPF = combIsoPF/pt;
	float eta = fabs(it_Mu->momentum.Eta());
	susy::Track& innerTrack = event->tracks[it_Mu->trackIndex];
	//float d0 = d0correction(event->beamSpot,combinedTrack);
	float numberOfValidPixelHits = innerTrack.numberOfValidPixelHits;
	float d0 = innerTrack.d0();
	float dZ = innerTrack.dz();
	float numberOfValidTrackerHits = it_Mu->nValidTrackerHits;
	float numberOfValidMuonHits = it_Mu->nValidMuonHits;
	float chi2OverNdof = innerTrack.normChi2();
	//cout<<"Muon d0: "<<d0<<endl;
	//cout<<"Muon dZ: "<<dZ<<endl;

	h_MuPt->Fill(pt);
	h_Mud0->Fill(d0);
	h_MudZ->Fill(dZ);
	h_MuCombIso->Fill(combIso);
	h_MuRelIso->Fill(relIso);
	h_MuEcalIso->Fill(it_Mu->ecalIsoR03);
	h_MuHcalIso->Fill(it_Mu->hcalIsoR03);
	h_MuTrackIso->Fill(it_Mu->trackIsoR03);
	if(pt < 15) continue;MuCount->Fill(3);
	if(eta > 2.6) continue;MuCount->Fill(4);
	if(chi2OverNdof > 10)continue;MuCount->Fill(5);
	if(numberOfValidMuonHits<=0)continue;MuCount->Fill(6);
	if(std::fabs(d0) > 0.2) continue;MuCount->Fill(7);
	if(std::fabs(dZ) > 0.5) continue;MuCount->Fill(8);
	if(numberOfValidPixelHits<=0)continue;MuCount->Fill(9);
	if(numberOfValidTrackerHits<11)continue;MuCount->Fill(10);
	//don't have trackerlayerswithmeasurement or numberofmatchedstations
	//id from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId

	if(relIsoPF > 0.2) continue;MuCount->Fill(11);

	Muons.push_back(&*it_Mu);
      }//end it_Mu muon loop
      //sort Muons by Pt
      std::sort(Muons.begin(), Muons.end(), EtGreater<susy::Muon>);
      if(printLevel>1)cout<<"Muons size= "<<Muons.size()<<endl;
      h_NumMuons->Fill((int)Muons.size());
      
      //----------JETS------------
      std::pair< std::vector<susy::PFJet*>,std::vector<susy::PFJet*> >   pfJets;
      std::vector<susy::PFJet*> pfJetsFail;
      bool isJetFail=false;
      /*
      // setup on-the-fly jet corrections for PFJets
      std::vector<JetCorrectorParameters> pfJECs;
    
      pfJECs.push_back(JetCorrectorParameters("../jec/Jec11_V1_AK5PF_L1FastJet.txt"));
      pfJECs.push_back(JetCorrectorParameters("../jec/Jec11_V1_AK5PF_L2Relative.txt"));
      pfJECs.push_back(JetCorrectorParameters("../jec/Jec11_V1_AK5PF_L3Absolute.txt"));
    
      FactorizedJetCorrector pfJetCorrector(pfJECs);
      */
      if(printLevel > 1) cout << "Find pfJets in the event." << endl;
  
      std::map<TString,susy::PFJetCollection>::iterator pfJets_it = event->pfJets.find("ak5");

      if(pfJets_it != event->pfJets.end()){

	susy::PFJetCollection& jetColl = pfJets_it->second;

	for(std::vector<susy::PFJet>::iterator it = jetColl.begin();
	    it != jetColl.end(); it++) {
	  /*
	    pfJetCorrector.setJetEta(it->momentum.Eta());
	    pfJetCorrector.setJetPt(it->momentum.Pt());
	    pfJetCorrector.setJetA(it->jetArea);
	    pfJetCorrector.setRho(Rho);
	    double corr = pfJetCorrector.getCorrection();
	  */
	  // set up corrections for PFJets
	  std::map<TString,Float_t>::iterator s_it_L1FastL2L3 = it->jecScaleFactors.find("L1FastL2L3");
	  std::map<TString,Float_t>::iterator s_it_L2L3       = it->jecScaleFactors.find("L2L3");
	  if (s_it_L1FastL2L3 == it->jecScaleFactors.end() || s_it_L2L3 == it->jecScaleFactors.end()) {
	    cout << "JEC is not available for this jet!!!" << endl;
	    continue;
	  }
	  float scale = s_it_L1FastL2L3->second;
	  float scaleMatch = s_it_L1FastL2L3->second/s_it_L2L3->second;
	  //if(printLevel > 2) cout << "PFJet stored (" << scale << "), onTheFly (" << corr << ")" << endl;

	  TLorentzVector corrP4 = scale * it->momentum;
	  TLorentzVector corrP4Match = scaleMatch * it->momentum;
	  bool go=true;
	  //This is for diJetPt Jet Matching
	  if( corrP4Match.Pt()>=20. && std::abs(corrP4Match.Eta())<=2.6 ){
	    it->momentum=corrP4Match;
	    pfJets.second.push_back(&*it); 
	  }
	  //This is for pfJets
	  if(   corrP4.Pt()>=30.
		&& std::abs(corrP4.Eta()) <= 2.6
		&& it->nConstituents>1
		//use uncorrected jet E for fractions
		&& it->neutralHadronEnergy/it->momentum.E()<0.99
		&& it->neutralEmEnergy/it->momentum.E()<0.99 ) {
	    if(std::abs(corrP4.Eta()) >= 2.4){
	      //cout<<"Jet Eta:            "<<it->momentum.Eta()<<endl;
	      //set jet momentum to L1FastL2L3 corrected one
	      it->momentum=corrP4;
	      //cout<<"Jet corrP4 Eta:     "<<it->momentum.Eta()<<endl;
	      pfJets.first.push_back(&*it);
	    }
	    else if( it->chargedMultiplicity>0 
		     && it->chargedHadronEnergy/it->momentum.E()>0
		     && it->chargedEmEnergy/it->momentum.E()<0.99){
	      it->momentum=corrP4;
	      pfJets.first.push_back(&*it);
	    }
	  }
	  else if(corrP4.Pt()>=30. && std::abs(corrP4.Eta()) <= 2.6){
	    go=true;
	    for(std::vector<susy::Electron*>::iterator ele_it = pfEles.begin();ele_it!=pfEles.end();ele_it++){
	      if(isSameObject((*ele_it)->momentum,it->momentum,0.5)){go=false;}
	    }
	    for(std::vector<susy::Muon*>::iterator mu_it = Muons.begin();mu_it!=Muons.end();mu_it++){
	      if(isSameObject((*mu_it)->momentum,it->momentum,0.5)){go=false;}
	    }
	    //cout<<"Jet with Pt>50 fails id criteria!  Event:"<<event->eventNumber<<"  Run:"<<event->runNumber<<endl;
	    if(go){
	      isJetFail=true;
	      /*cout<<"Jet         Pt:     "<<it->momentum.Pt()<<endl;
		cout<<"Jet CorrP4  Pt:     "<<corrP4.Pt()<<endl<<endl;*/
	      pfJetsFail.push_back(&*it);
	    }
	  }//end JetFail
	}// for jet
      }// if not end
      
      //std::sort(pfJets.begin(),pfJets.end(),EtGreater<susy::PFJet>);
  
      bool AtLeastOnePFJet = false;
      bool doJetReq=false,doJetReqFake=false;
      if(pfJets.first.size()>0) AtLeastOnePFJet=true;
      //------end JETS--------------

      if(printLevel > 1) cout << "Find loose and tight photons in the event." << endl;
      //----------
      //find photons, sort by energy
      //----------
      std::map<TString, std::vector<susy::Photon> >::iterator phoMap = event->photons.find("photons");
      //std::map<TString, std::vector<susy::Photon> >::iterator phoMap = event->photons.find("pfPhotonTranslator:pfphot");
      if(phoMap != event->photons.end()) {
      
	//   Need to sort Photons and select two with largest Pt.  
      
	std::vector<susy::Photon*>   pho_Cands;
	std::vector<susy::Photon*>   pho_CandsNminus3;
	std::vector<susy::Photon*>   pho_Cands_NoEtCut;
	std::vector<susy::Photon*>   pho_Cands_EBEE;
	std::vector<susy::Photon*>   fake_Cands;

	//Criteria to Filter events
	int filterCheck=0, GreaterThan43=0;
    
	//loop over photon collection 
	for(std::vector<susy::Photon>::iterator it_Pho = phoMap->second.begin(); it_Pho != phoMap->second.end(); it_Pho++) {
	  if(printLevel > 1) cout<<"looping over photon collection"<<endl;

	
	  //do things that want all photons with no cuts here:
	  if(it_Pho->isPF())pfPho++;
	  else notPfPho++;

	  h_SeedTime->Fill(it_Pho->seedTime);
	  h_MetVsSeedTime->Fill(it_Pho->seedTime, met->met());
	  h_SeedTimeVsEta->Fill(it_Pho->caloPosition.Eta(),it_Pho->seedTime);
	  h_SeedTimeVSE->Fill(it_Pho->seedTime,it_Pho->momentum.E());
	  h_SigIetaIeta_allPho->Fill(it_Pho->sigmaIetaIeta);
	  h_PhoPt->Fill(it_Pho->momentum.Pt());
	 
	  h_Pho_CombinedIsoDR03->Fill(it_Pho->ecalRecHitSumEtConeDR03+it_Pho->hcalTowerSumEtConeDR03()+it_Pho->trkSumPtHollowConeDR03);
	  h_Pho_ChargedHadronIso->Fill(it_Pho->chargedHadronIso);
	  h_Pho_NeutralHadronIso->Fill(it_Pho->neutralHadronIso); 
	  h_Pho_PhotonIso->Fill(it_Pho->photonIso);
	  h_Pho_PfCombinedIso->Fill(it_Pho->chargedHadronIso+it_Pho->neutralHadronIso+it_Pho->photonIso);
	  h_Pho_ChargedHadronIsoDeposit->Fill(it_Pho->chargedHadronIsoDeposit);
	  h_Pho_NeutralHadronIsoDeposit->Fill(it_Pho->neutralHadronIsoDeposit); 
	  h_Pho_PhotonIsoDeposit->Fill(it_Pho->photonIsoDeposit);
	  h_Pho_PfCombinedIsoDeposit->Fill(it_Pho->chargedHadronIsoDeposit+it_Pho->neutralHadronIsoDeposit+it_Pho->photonIsoDeposit);

	  bool passR9  = (useTrigger ? PassTrigger("HLT_Photon36_R9Id85_Photon22_R9Id85_v"  ) : true);
	  if(passR9)h_R9fromR9trig->Fill(it_Pho->r9);

	  h_EcalIsoDR03_allPho->Fill(it_Pho->ecalRecHitSumEtConeDR03);
	  h_EcalIsoDR03RhoCorr_allPho->Fill(it_Pho->ecalRecHitSumEtConeDR03-PUCorr_ECAL*Rho25);
	  h_HcalIsoDR03_allPho->Fill(it_Pho->hcalTowerSumEtConeDR03());
	  h_HcalIsoDR03RhoCorr_allPho->Fill(it_Pho->hcalTowerSumEtConeDR03()-PUCorr_HCAL*Rho25);
	  h_TrackIsoDR03_allPho->Fill(it_Pho->trkSumPtHollowConeDR03);
	  h_CombIsoDR03_allPho->Fill(it_Pho->ecalRecHitSumEtConeDR03+it_Pho->hcalTowerSumEtConeDR03()+it_Pho->trkSumPtHollowConeDR03);
	  h_CombIsoDR03RhoCorr_allPho->Fill(it_Pho->ecalRecHitSumEtConeDR03-PUCorr_ECAL*Rho25+it_Pho->hcalTowerSumEtConeDR03()-PUCorr_HCAL*Rho25+it_Pho->trkSumPtHollowConeDR03-PUCorr_TRACK*Rho25);


	  //ignore event if there are not at least 2 photon candidates 
	  if(event->isRealData){
	    if(phoMap->second.size()<2) break;//break not continue since we want to exit phoMap iterator if there are not at least two photons in event
	  }

	  if(std::abs(it_Pho->caloPosition.Eta())>1.0){
	    PUCorr_chargedHadron=0.010;
	    PUCorr_neutralHadron=0.057;
	    PUCorr_photon=0.130;
	  }

	  bool phoCand=false, fakeCand=false;
	  //----------------set up cuts-------------------
	  float ecalIsoDR03=it_Pho->ecalRecHitSumEtConeDR03;
	  float hcalIsoDR03=it_Pho->hcalTowerSumEtConeDR03();
	  float trackIsoDR03=it_Pho->trkSumPtHollowConeDR03;
	  float ecalIsoDR04=it_Pho->ecalRecHitSumEtConeDR04;
	  float hcalIsoDR04=it_Pho->hcalTowerSumEtConeDR04();
	  float trackIsoDR04=it_Pho->trkSumPtHollowConeDR04;
	  float chargedHadronIso=it_Pho->chargedHadronIso;
	  float neutralHadronIso=it_Pho->neutralHadronIso;
	  float photonIso=it_Pho->photonIso;
	  float PhoEt=it_Pho->momentum.Et();
	  // fiducial cuts. Look for only barrel now
	  bool etaCut = (std::abs(it_Pho->caloPosition.Eta()) < susy::etaGapBegin);

	  // Spike cleaning
	  h_r9->Fill(it_Pho->r9);
	  h_hadOverEm->Fill(it_Pho->hadronicOverEm);
	  h_hadTowOverEm->Fill(it_Pho->hadTowOverEm);
	  bool isSpike = (it_Pho->r9 > 1.0 || it_Pho->sigmaIetaIeta<=0.001 || it_Pho->sigmaIphiIphi<=0.001);
	  //if(isSpike) continue;
	  h_SeedTime_afterR9->Fill(it_Pho->seedTime);
	  h_SeedTimeVSE_afterR9->Fill(it_Pho->seedTime,it_Pho->momentum.E());

	  // Et cuts, 25 GeV for trailing photons. Will apply tighter for the leading one.
	  bool EtCut = (PhoEt > 25.0);
	  // cuts containing EE cases for later use, but EE is not used for the time being.

	  // H/E (in trigger, 0.15 for EB, 0.10 for EE)
	  bool heCut = (it_Pho->hadronicOverEm < 0.05);
	  bool heTowCut = (it_Pho->hadTowOverEm < 0.05);
	  bool heCutLoose = (it_Pho->hadronicOverEm < 0.1);
	  // sigma_ietaieta (in trigger 0.014 for EB, 0.034 for EE)
	  bool sIetaCut = (it_Pho->sigmaIetaIeta < 0.011);

	  // Ecal Isolation
	  bool ecalIsoCut = (ecalIsoDR04 < 4.2 + 0.006*PhoEt + 0.183*RhoB);
	  bool ecalIsoCutdR03 = (ecalIsoDR03 < 4.2 + 0.006*PhoEt);
	  // Hcal Isolation
	  bool hcalIsoCut = (hcalIsoDR04 < 2.2 + 0.0025*PhoEt + 0.062*RhoB);
	  bool hcalIsoCutdR03 = (hcalIsoDR03 < 2.2 + 0.0025*PhoEt);
	  // Track Isolation
	  bool trackIsoCut = (trackIsoDR04 < 2.0 + 0.001*PhoEt + 0.0167*RhoB);
	  bool trackIsoCutdR03 = (trackIsoDR03 < 2.0 + 0.001*PhoEt);

	  //PF Isolation
	  bool PfCHisoLoose = (chargedHadronIso - PUCorr_chargedHadron*Rho25) < 2.6;
	  bool PfCHisoMed = (chargedHadronIso - PUCorr_chargedHadron*Rho25) < 1.5;
	  bool PfCHisoTight = (chargedHadronIso - PUCorr_chargedHadron*Rho25) < 0.7;
	  bool PfNHisoLoose = (neutralHadronIso - 0.04*PhoEt - PUCorr_neutralHadron*Rho25) < 3.5;
	  bool PfNHisoMed = (neutralHadronIso - 0.04*PhoEt - PUCorr_neutralHadron*Rho25) < 1.0;
	  bool PfNHisoTight = (neutralHadronIso - 0.04*PhoEt - PUCorr_neutralHadron*Rho25) < 0.4;
	  bool PfPisoLoose = (photonIso - 0.005*PhoEt - PUCorr_photon*Rho25) < 1.3;
	  bool PfPisoMed = (photonIso - 0.005*PhoEt - PUCorr_photon*Rho25) < 0.7;
	  bool PfPisoTight = (photonIso - 0.005*PhoEt - PUCorr_photon*Rho25) < 0.5;


	  bool PfIsoCutLoose = PfCHisoLoose && PfNHisoLoose && PfPisoLoose;
	  bool PfIsoCutMed = PfCHisoMed && PfNHisoMed && PfPisoMed;
	  bool PfIsoCutTight = PfCHisoTight && PfNHisoTight && PfPisoTight;

	  //combined cuts
	  //  bool combIsoCut =( ( ecalIsoDR03 + hcalIsoDR03 + trackIsoDR03 ) < 6.0 );
	  bool combIsoCut =( ( ecalIsoDR03+hcalIsoDR03+trackIsoDR03 - (PUCorr_ECAL+PUCorr_HCAL+PUCorr_TRACK)*Rho25 ) < 6.0 );
	  bool combIsoFakeHighCut =( ( ecalIsoDR03+hcalIsoDR03+trackIsoDR03 - (PUCorr_ECAL+PUCorr_HCAL+PUCorr_TRACK)*Rho25 ) < 20.0 );
	  bool PfIsoFakeHighCut = (chargedHadronIso+neutralHadronIso+photonIso - (PUCorr_chargedHadron+PUCorr_neutralHadron+PUCorr_photon)*Rho25)<27.0;
	  bool stdIsoCut = (ecalIsoCut+hcalIsoCut+trackIsoCut);
	  //hlt cuts
	  bool ecalIsoVLcut = ( ecalIsoDR03  < 6.0 + 0.012*PhoEt );
	  bool hcalIsoVLcut = ( hcalIsoDR03 < 4.0 + 0.005*PhoEt );
	  bool trackIsoVLcut = ( trackIsoDR03  < 4.0 + 0.002*PhoEt );
	  bool IsoVLCut = ( ecalIsoVLcut && hcalIsoVLcut && trackIsoVLcut);
	  bool CaloIdLCut = ( (it_Pho->hadronicOverEm < 0.15) && (it_Pho->sigmaIetaIeta < 0.014) );
	  bool R9TriggerCut = (it_Pho->r9 > 0.8);
	  //Pixel cut
	  bool pixelCut = (it_Pho->nPixelSeeds == 0);
	  //bool pixelCut=(it_Pho->passelectronveto);
	  bool passEleVeto = it_Pho->passelectronveto;

	  bool R9Id85=it_Pho->r9>0.85;
	  bool CaloId10=( it_Pho->hadronicOverEm<0.1 && it_Pho->sigmaIetaIeta < 0.014 );
	  bool Iso50=(ecalIsoDR03<5.0 + 0.012*PhoEt
		      && hcalIsoDR03<5.0 + 0.005*PhoEt
		      && trackIsoDR03<5.0+ 0.002*PhoEt
		      );
	  bool R9Id85orCaloId10andIso50Cut = R9Id85 || ( CaloId10 && Iso50 );	 
	  bool PhoCutEBEE = (/*etaCut && */EtCut && !isSpike && heCut && combIsoCut && sIetaCut /*&& CaloId10 && Iso50*/);
	  bool PhoCut = (etaCut && EtCut && !isSpike && heCut && combIsoCut && sIetaCut /*&& CaloId10 && Iso50*/);
	  bool PhoCutPF = (etaCut && EtCut && !isSpike && heTowCut && PfIsoCutLoose && sIetaCut /*&& CaloId10 && Iso50*/);
	  //bool PhoCut = (etaCut && EtCut && !isSpike && heCut && combIsoCut && sIetaCut && R9Id85orCaloId10andIso50Cut);
	  bool StdPhoCut = (etaCut && EtCut && !isSpike && heCut && stdIsoCut && sIetaCut && CaloId10 && Iso50);
	  bool PhoCutNoEtCut = (etaCut && !isSpike && heCut && combIsoCut && sIetaCut && CaloIdLCut && !pixelCut && IsoVLCut);
	  bool PhoCutNminus3 = (etaCut && EtCut && !isSpike && heCut && sIetaCut && pixelCut && CaloId10 && Iso50);
	  bool FakeCut = ( etaCut && EtCut && !isSpike && heCut && pixelCut && combIsoFakeHighCut /*&& R9Id85orCaloId10andIso50Cut*/ ) && (!combIsoCut || !sIetaCut);
	  bool FakeCutPF = ( etaCut && EtCut && !isSpike && heTowCut && passEleVeto && PfIsoFakeHighCut/* && R9Id85orCaloId10andIso50Cut*/ ) && (!PfIsoCutLoose || !sIetaCut);	  
	  bool StdFakeCut = ( etaCut && EtCut && !isSpike && heCut && pixelCut && combIsoFakeHighCut && R9Id85orCaloId10andIso50Cut ) && (!hcalIsoCut || !trackIsoCut || !sIetaCut);

	  //bool FakeCut = ( etaCut && EtCut && !isSpike && heCut && pixelCut && combIsoFakeHighCut && R9Id85orCaloId10andIso50Cut && (!combIsoCut || !sIetaCut) );


	  bool PhoCutNew = (etaCut && EtCut && !isSpike && heCut && combIsoCut && sIetaCut && CaloId10 && Iso50);
	  bool PhoCutOld = (etaCut && EtCut && !isSpike && heCut && combIsoCut && sIetaCut && CaloIdLCut && IsoVLCut);

	  if(PhoCutOld)oldP++;
	  if(PhoCutNew)newP++;

	  bool passHLT_Pho  = (useTrigger ? PassTrigger("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v"  ) : true);
	  
	  //bool passHLT_Pho  = (useTrigger ? ( PassTrigger("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v") ) : true );

	  bool passHLT_Fake = (useTrigger ? ( PassTrigger("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v")
					      || PassTrigger("HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v") 
					      || PassTrigger("HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v")
					      || PassTrigger("HLT_Photon36_R9Id85_Photon22_R9Id85_v") ) : true);

	  //----set up filtered events-----------
	  if(PhoEt > 30.0 && heCut) filterCheck++;
	  if(PhoEt > 43.0) GreaterThan43++;
	  //-------------------------------------
	
	  //Now fill pho_Cands and fake_Cands
	  //pho_Cands is for e and g
	  //fake_Cands is for ff only
	  if( PhoCut/*PF/* && passHLT_Pho*/){
	    phoCand=true;  pho_Cands.push_back(&*it_Pho);
	  }
	  if( FakeCut/*PF/* && passHLT_Fake*/){
	    fakeCand=true; fake_Cands.push_back(&*it_Pho);
	  }
	  if( PhoCutNoEtCut && passHLT_Pho){
	    pho_Cands_NoEtCut.push_back(&*it_Pho);
	  }

	  if( PhoCutNminus3 && passHLT_Pho){
	    /*phoCandNminus3=true;  */pho_CandsNminus3.push_back(&*it_Pho);
	  }
	  if(PhoCutEBEE){pho_Cands_EBEE.push_back(&*it_Pho);}

	  if(etaCut && EtCut && !isSpike && heCut && CaloIdLCut && (IsoVLCut || R9TriggerCut)){
	    h_CombIsoNminus4->Fill(ecalIsoDR03-PUCorr_ECAL*Rho25+hcalIsoDR03-PUCorr_HCAL*Rho25+trackIsoDR03-PUCorr_TRACK*Rho25);
	    h_SigIetaIetaNminus4->Fill(it_Pho->sigmaIetaIeta);
	    h_EcalIsoDR03VsCombIsoDR03Nminus4->Fill(ecalIsoDR03-PUCorr_ECAL*Rho25+hcalIsoDR03-PUCorr_HCAL*Rho25+trackIsoDR03-PUCorr_TRACK*Rho25,ecalIsoDR03-PUCorr_ECAL*Rho25);
	    h_HcalIsoDR03VsCombIsoDR03Nminus4->Fill(ecalIsoDR03-PUCorr_ECAL*Rho25+hcalIsoDR03-PUCorr_HCAL*Rho25+trackIsoDR03-PUCorr_TRACK*Rho25,hcalIsoDR03-PUCorr_HCAL*Rho25);
	    h_TrackIsoDR03VsCombIsoDR03Nminus4->Fill(ecalIsoDR03-PUCorr_ECAL*Rho25+hcalIsoDR03-PUCorr_HCAL*Rho25+trackIsoDR03-PUCorr_TRACK*Rho25,trackIsoDR03-PUCorr_TRACK*Rho25);
	  }

	  if( etaCut && EtCut && !isSpike && heCut && sIetaCut && pixelCut && CaloIdLCut && (IsoVLCut || R9TriggerCut)){
	    if(met->met()<30){
	      h_BgroundCombIsoDR03Nminus3->Fill(ecalIsoDR03+hcalIsoDR03+trackIsoDR03);
	      h_BgroundCombIsoDR04Nminus3->Fill(ecalIsoDR04+hcalIsoDR04+trackIsoDR04);
	    }
	    if(ecalIsoCut && trackIsoCut){
	      h_HcalIsoDR03Nminus1->Fill(hcalIsoDR03);
	      h_HcalIsoDR04Nminus1->Fill(hcalIsoDR04);
	      if(met->met()<30){
		h_BgroundHcalIsoDR03Nminus1->Fill(hcalIsoDR03);
	      }
	    }
	    if(hcalIsoCut && trackIsoCut){
	      h_EcalIsoDR03Nminus1->Fill(ecalIsoDR03);
	      h_EcalIsoDR04Nminus1->Fill(ecalIsoDR04);
	      if(met->met()<30){
		h_BgroundEcalIsoDR03Nminus1->Fill(ecalIsoDR03);
	      }
	    }
	    if(ecalIsoCut && hcalIsoCut){
	      h_TrackIsoDR03Nminus1->Fill(trackIsoDR03);
	      h_TrackIsoDR04Nminus1->Fill(trackIsoDR04);
	      if(met->met()<30){
		h_BgroundTrackIsoDR03Nminus1->Fill(trackIsoDR03);
	      }
	    }
	  }
	

	  // Fill lots of plots
	  if(etaCut && EtCut && !isSpike){
	    h_EcalIsoVsNVertex->Fill(NVertex+1,ecalIsoDR04);
	    h_HcalIsoVsNVertex->Fill(NVertex+1,hcalIsoDR04);
	    h_TrackIsoVsNVertex->Fill(NVertex+1,trackIsoDR04);
	    if(heCut){
	      /*h_Rho_VS_NVertex_afterHE->Fill(NVertex,Rho);
		h_Rho_VS_NVertex_Prof_afterHE->Fill(NVertex,Rho);*/
	      //DR03 cone study
	      h_EcalIsoDR03->Fill(ecalIsoDR03);
	      h_HcalIsoDR03->Fill(hcalIsoDR03);
	      h_TrackIsoDR03->Fill(trackIsoDR03);
	      h_EcalIsoDR04->Fill(ecalIsoDR04);
	      h_HcalIsoDR04->Fill(hcalIsoDR04);
	      if(hcalIsoCut){
		h_EcalIso->Fill(ecalIsoDR04);
	      }
	      if(ecalIsoCut){
		h_HcalIso->Fill(hcalIsoDR04);
		if(hcalIsoCut){
		  h_EcalIso_Preselection->Fill(ecalIsoDR04);
		  h_HcalIso_Preselection->Fill(hcalIsoDR04);
		  h_TrackIso_Preselection->Fill(trackIsoDR04);
		  h_EcalIso_Preselection_NVertex_Corrected->Fill(ecalIsoDR04 - 0.127*NVertex);
		  h_HcalIso_Preselection_NVertex_Corrected->Fill(hcalIsoDR04 - 0.0428*NVertex);
		  h_TrackIso_Preselection_NVertex_Corrected->Fill(trackIsoDR04 - (-0.03527)*NVertex);
		  h_EcalIso_Preselection_Rho_Corrected->Fill(ecalIsoDR04 - 0.2888*Rho);
		  h_HcalIso_Preselection_Rho_Corrected->Fill(hcalIsoDR04 - 0.106*Rho);
		  h_TrackIso_Preselection_Rho_Corrected->Fill(trackIsoDR04 - 0.09845*Rho);
		  h_TrackIso->Fill(trackIsoDR04);
		  h_PixelSeeds->Fill(it_Pho->nPixelSeeds);
		  if(AtLeastOnePFJet){h_PixelSeeds_JetReq->Fill(it_Pho->nPixelSeeds);}		  
		}
	      }
	    }//if(heCut)
	  }//if(etaCut && EtCut && !isSpike)
     
	  if(printLevel > 1) cout<<"End of Photon Loop"<<endl;
	}//for(it_Pho)

	//sort pho_Cands and fake_Cands by Pt
	std::sort(pho_Cands.begin(), pho_Cands.end(), EtGreater<susy::Photon>);
	if(printLevel>1)cout<<"phoCands size= "<<pho_Cands.size()<<endl;
	std::sort(pho_CandsNminus3.begin(), pho_CandsNminus3.end(), EtGreater<susy::Photon>);
	std::sort(pho_Cands_EBEE.begin(), pho_Cands_EBEE.end(), EtGreater<susy::Photon>);
	std::sort(fake_Cands.begin(), fake_Cands.end(), EtGreater<susy::Photon>);
	if(printLevel>1)cout<<"fakeCands size= "<<fake_Cands.size()<<endl;

	//fill number of candidate plots
	int g=0,e=0;
	for(int k=0;k<int(pho_Cands.size());k++){
	  if(pho_Cands[k]->nPixelSeeds==0)g++;
	  else e++;
	}
	h_NumE_NumG->Fill(e,g);
	h_NumE_NumG_NumF->Fill(e,g,fake_Cands.size());
	h_NPhoCands->Fill(pho_Cands.size());
	h_NFakeCands->Fill(fake_Cands.size());
	h_NPho_NFake->Fill(pho_Cands.size(),fake_Cands.size());
      
	for(std::vector<susy::Photon*>::iterator pho_it_EBEE = pho_Cands_EBEE.begin();pho_it_EBEE!=pho_Cands_EBEE.end();pho_it_EBEE++){
	  h_TrigPhosEta->Fill((*pho_it_EBEE)->caloPosition.Eta());
	  if(pho_it_EBEE==pho_Cands_EBEE.begin() || pho_it_EBEE==++(pho_Cands_EBEE.begin()))h_TrigPhosEta_TopTwo->Fill((*pho_it_EBEE)->caloPosition.Eta());
	}

	//Make e-mu and mu-mu samples
	if(printLevel>0)cout<<"Fill eMu sample"<<endl;
	if( pho_Cands.size()>0 && Muons.size()>0 && ( (*pho_Cands.begin())->momentum.Et()>40. || (*Muons.begin())->momentum.Et()>40. ) ){
	  //First eMu
	  for(std::vector<susy::Photon*>::iterator pho_it_eMu = pho_Cands.begin();pho_it_eMu!=pho_Cands.end();pho_it_eMu++){
	    if((*pho_it_eMu)->nPixelSeeds>0){//want electrons only
	      for(std::vector<susy::Muon*>::iterator mu_it_eMu = Muons.begin();mu_it_eMu!=Muons.end();mu_it_eMu++){
		if((*pho_it_eMu)->momentum.Et()>40. || (*mu_it_eMu)->momentum.Et()>40.){
		  if( InvariantMass((*pho_it_eMu)->momentum,(*mu_it_eMu)->momentum)>71 && InvariantMass((*pho_it_eMu)->momentum,(*mu_it_eMu)->momentum)<111){
		    h_eMu_InvarMass_FullRange->Fill(InvariantMass((*pho_it_eMu)->momentum,(*mu_it_eMu)->momentum));
		    h_eMu_MET_FullRange->Fill(met->met());
		    if( InvariantMass((*pho_it_eMu)->momentum,(*mu_it_eMu)->momentum)>81 && InvariantMass((*pho_it_eMu)->momentum,(*mu_it_eMu)->momentum)<101){
		      h_eMu_InvarMass->Fill(InvariantMass((*pho_it_eMu)->momentum,(*mu_it_eMu)->momentum));
		      h_eMu_MET->Fill(met->met());
		      neMu++;
		    }//end 81-101 InvarMass
		  }//end 71-111 InvarMass
		}//end 40 gev requirements
	      }//end mu_it iterator
	    }//end pixel seed requirement
	  }//end pho_it_eMu iterator -- end eMu
	}//end e-mu
	if(printLevel>0)cout<<"Fill MuMu sample"<<endl;
	//Now mu-mu	
	if(Muons.size()>1 && (*Muons.begin())->momentum.Et()>40.){
	  for(std::vector<susy::Muon*>::iterator mu_it1_mumu = Muons.begin();mu_it1_mumu<(Muons.end())-1;mu_it1_mumu++){
	    for(std::vector<susy::Muon*>::iterator mu_it2_mumu = ++mu_it1_mumu;mu_it2_mumu!=Muons.end();mu_it2_mumu++){
	      if((*mu_it1_mumu)->momentum.Et()>40. || (*mu_it2_mumu)->momentum.Et()>40.){
		float muInvMass=InvariantMass((*mu_it1_mumu)->momentum,(*mu_it2_mumu)->momentum);
		if( muInvMass>71 && muInvMass<111){
		  h_MuMu_InvarMass_FullRange->Fill(muInvMass);
		  h_MuMu_MET_FullRange->Fill(met->met());
		  if( muInvMass>81 && muInvMass<101){
		    h_MuMu_InvarMass->Fill(muInvMass);
		    h_MuMu_MET->Fill(met->met());
		    nMuMu++;
		  }//end 81-101 InvarMass
		}//end 71-111 InvarMass
	      }//end 40 gev requirements
	    }//end second iterator ( mu_it2_mumu )
	  }//end first iterator ( mu_it1_mumu ) -- end MuMu
	}//end mu-mu
	  
	//Now back to photons

	if(printLevel>1)cout<<"Just sorted pho and fake Cands, about to create PhoOne and PhoTwo"<<endl;
	//Now we have sorted into two candidate samples, pho_Cands and fake_Cands, which are mutually exclusive on photon to photon basis
	//
	//need to make sure I can't get e.g. a gg and ff from one event.
  
	bool isgg=false, isee=false, iseg=false, isff=false, isgg_JetReq=false, isee_JetReq=false, iseg_JetReq=false, isff_JetReq=false;
	bool PhosTooCloseDR=true, FakesTooCloseDR=true, PhosTooCloseDphi=true, FakesTooCloseDphi=true;
	bool JetIsolatedFromBothPhos=false, JetIsolatedFromBothFakes=false;
	bool TwoPhosAndTwoFakes=false;
	bool ThreePhoEvent=false, ThreeFakeEvent=false, FourPhoEvent=false, FourFakeEvent=false;
	bool UsingThirdPho=false,UsingThirdFake=false;
	bool breakPho=false, breakPho_JetReq=false, breakFake=false, breakFake_JetReq=false;
	bool HasDiJetPt=false;
	int JetCounter=0;
	float RazrMr=0.,RazrR2=0.;
	float PhotonLessHt=0.,alphaT=0.;
	bool IsGoodJet=true;
	//------------check of how many events have >=3 pho_Cands or fake_Cands------------	
	if(pho_Cands.size()>=3){nThreePhoEvents++;ThreePhoEvent=true;if(pho_Cands.size()>=4){nFourPhoEvents++;FourPhoEvent=true;}}
	if(fake_Cands.size()>=3){nThreeFakeEvents++;ThreeFakeEvent=true;if(fake_Cands.size()>=4){nFourFakeEvents++;FourFakeEvent=true;}}
	//---------------------------------------------------------------------------------
        
	//If there are at least 2 photon candidates, check pixel seeds and sort into gg, eg, ee
	if( pho_Cands.size()>=2 && (*pho_Cands.begin())->momentum.Et()>40.){ //make sure top pT object has at least 40GeV
	  ndiPhoCands++;if(printLevel>1)cout<<"pho_Cands.size()>=2 && (*pho_Cands.begin())->momentum.Et()>40."<<endl;
	  if(fake_Cands.size()>=2){TwoPhosAndTwoFakes=true;nTwoPhosAndTwoFakes++;}
	  //Now do no jet req case
	  if(printLevel>5)cout<<"Now do phoCands no jet req case"<<endl;
	  //This mades sure the two pho objects are separated by dPhi>0.05 for no jet req case
	  PhosTooCloseDR=true;PhosTooCloseDphi=true;breakPho=false;doJetReq=false;JetIsolatedFromBothPhos=false;
	  for(std::vector<susy::Photon*>::iterator pho_it = pho_Cands.begin();pho_it<(pho_Cands.end())-1;pho_it++){
	    if(breakPho)break;
	    PhoOne=*pho_it;if(PhoOne->momentum.Pt()<40)break;
	    for(std::vector<susy::Photon*>::iterator pho_it2 = ++pho_it;pho_it2!=pho_Cands.end();pho_it2++){
	      //require dPhi>0.05,dR>0.6
	      if(   !isSameObject(PhoOne->caloPosition,(*pho_it2)->caloPosition,0.6) 
		    && !tooClosePhi(PhoOne->caloPosition,(*pho_it2)->caloPosition,0.05)    ){
		PhoTwo=*pho_it2;PhosTooCloseDR=false;PhosTooCloseDphi=false;breakPho=true;
		//This checks if there is at least one loose jet separated by 0.5 from both phos
		if(!PhosTooCloseDR && !PhosTooCloseDphi && AtLeastOnePFJet){
		  //use pfJets.first because this is for all jets
		  for(std::vector<susy::PFJet*>::iterator jet_it = pfJets.first.begin(); jet_it != pfJets.first.end(); jet_it++){
		    IsGoodJet=true;
		    //first clean from eles and muons
		    for(std::vector<susy::Electron*>::iterator ele_it = pfEles.begin();ele_it!=pfEles.end();ele_it++){
		      if(isSameObject((*ele_it)->momentum,(*jet_it)->momentum,0.5)){IsGoodJet=false;break;}
		    }
		    for(std::vector<susy::Muon*>::iterator mu_it = Muons.begin();mu_it!=Muons.end();mu_it++){
		      if(isSameObject((*mu_it)->momentum,(*jet_it)->momentum,0.5)){IsGoodJet=false;break;}
		    }
		    //now clean from the 2 pho objects
		    if(IsGoodJet){
		      if(!isSameObject(PhoOne->caloPosition,(*jet_it)->momentum,0.5) && !isSameObject(PhoTwo->caloPosition,(*jet_it)->momentum,0.5)){
			JetIsolatedFromBothPhos=true;ndiPhoCands_JetReq++;break;
		      }
		    }
		  }
		}
		if(JetIsolatedFromBothPhos==true)doJetReq=true;
		break;
	      }
	      else{
		if(isSameObject(PhoOne->caloPosition,(*pho_it2)->caloPosition,0.6)){
		  PhosTooCloseDR=true;nPhosFailDR++;
		  if(printLevel>0)cout <<"PhosFailDR!"<< "  Run: "<<event->runNumber<<"  Event: "<<event->eventNumber<<"\n";
		}
		if(tooClosePhi(PhoOne->caloPosition,(*pho_it2)->caloPosition,0.05)){
		  PhosTooCloseDphi=true;nPhosFailDphi++;
		  if(printLevel>0)cout <<"PhosFailDphi!"<< "  Run: "<<event->runNumber<<"  Event: "<<event->eventNumber<<"\n";
		}	     
	      }
	    }
	  }
	  if(printLevel>1)cout<<"doJetReq: "<<doJetReq<<endl;
	
	  if(PhoOne->momentum.Et()>40.){
	    if(!PhosTooCloseDR && !PhosTooCloseDphi){
	      if(printLevel>5)cout<<"Now inside phoCands no jet req case"<<endl;
	      RazrMr = GetRazrMr(PhoOne,PhoTwo);
	      RazrR2 = GetRazrR2(PhoOne,PhoTwo,met);
	      //try invariant mass cut to get rid of crap in tail
	      float InvMass=InvariantMass(PhoOne->momentum,PhoTwo->momentum);
	      float diEMpt=GetDiEmPt(PhoOne,PhoTwo);
	      float diJETpt=0.;
	      if(InvMass>-9999.){
		//find out how many isolated jets
		JetCounter=0;
		//use pfJets.first here because we want all jets
		for(std::vector<susy::PFJet*>::iterator jet_it = pfJets.first.begin(); jet_it != pfJets.first.end(); jet_it++){	
		  //to be counted as a jet, it must be dr>0.5 separated from the 2 objects, plus pfEles and pfMuons	
		  IsGoodJet=true;
		  bool check=true;
		  if(isSameObject(PhoOne->caloPosition,(*jet_it)->momentum,0.5)){IsGoodJet=false;}		  
		  if(isSameObject(PhoTwo->caloPosition,(*jet_it)->momentum,0.5)){IsGoodJet=false;}
		  check=IsGoodJet;
		  for(std::vector<susy::Electron*>::iterator ele_it = pfEles.begin();ele_it!=pfEles.end();ele_it++){
		    if(isSameObject((*ele_it)->momentum,(*jet_it)->momentum,0.5)){IsGoodJet=false;}
		  }
		  if(check && !IsGoodJet) nCnt[6]++; //number of jets that fail ele cleaning from e/g
		  check=true;
		  for(std::vector<susy::Muon*>::iterator mu_it = Muons.begin();mu_it!=Muons.end();mu_it++){
		    if(isSameObject((*mu_it)->momentum,(*jet_it)->momentum,0.5)){IsGoodJet=false;}
		  }
		  if(check && !IsGoodJet) nCnt[7]++; //number of jets that fail muon cleaning from e/g
		  if(IsGoodJet)JetCounter++;
		}
		//Try to match pho objects to jets
		//here use pfJets.second because we want to use l1FastL2L2/L2L3 for jet matching
		HasDiJetPt=false;
		MatchPhosToJets(PhoOne,PhoTwo,pfJets.second,DiJetPtJet1,DiJetPtJet2,HasDiJetPt,0.3);
		if(HasDiJetPt) diJETpt=GetDiJetPt(DiJetPtJet1,DiJetPtJet2);
		//look for pixel seeds - if neither have, -> gg
		float dPhi = std::fabs(TVector2::Phi_mpi_pi(PhoOne->caloPosition.Phi() - PhoTwo->caloPosition.Phi()));
		if(HasDiJetPt) alphaT = GetAlphaT(DiJetPtJet1->momentum, DiJetPtJet2->momentum);
		else alphaT = GetAlphaT(PhoOne->momentum, PhoTwo->momentum);
		if(HasDiJetPt) PhotonLessHt=GetPhotonLessHt(met->sumEt,DiJetPtJet1->momentum,DiJetPtJet2->momentum);
		else PhotonLessHt=GetPhotonLessHt(met->sumEt,PhoOne->momentum,PhoTwo->momentum);
		if( PhoOne->nPixelSeeds==0 && PhoTwo->nPixelSeeds==0 ){
		  if(printLevel>1)cout<<"Inside gg"<<endl;
		  //require seed time within +-3ns for no jet requirement case
		  //if(fabs(PhoOne->seedTime)<3. && fabs(PhoTwo->seedTime)<3.){
		  ngg++;isgg=true;
		  ggPTs=make_pair(PhoOne->momentum.Pt(),PhoTwo->momentum.Pt());
		  h_NVertexVsMET_gg->Fill(met->met(),NVertex);
		  h_ggPerInstLumi->Fill(event->avgInsRecLumi);
		  if(NVertex<10)h_ggMet_NV0_10->Fill(met->met());
		  else if(NVertex>=10 && NVertex<15)h_ggMet_NV10_15->Fill(met->met());
		  else if(NVertex>=15)h_ggMet_NV15up->Fill(met->met());
		  if(event->isRealData && met->met()>100)cout<<"gg Event with MET="<<met->met()<<"  Run: "<<event->runNumber<<"  Lumi: "<<event->luminosityBlockNumber<<"  Event: "<<event->eventNumber<<"  PhoOne pT:"<<PhoOne->momentum.Pt()<<"  PhoTwo pT:"<<PhoTwo->momentum.Pt()<<endl;
		  ggRunEvent.push_back(make_pair(event->runNumber,event->eventNumber));
		  h_SeedTimeVsEta_gg->Fill(PhoOne->caloPosition.Eta(),PhoOne->seedTime);
		  h_SeedTimeVsEta_gg->Fill(PhoTwo->caloPosition.Eta(),PhoTwo->seedTime);
		  h_MetVsSeedTime_gg->Fill(PhoOne->seedTime, met->met());
		  h_MetVsSeedTime_gg->Fill(PhoTwo->seedTime, met->met());
		  if(isJetFail){bool go2=true;
		    for(std::vector<susy::PFJet*>::iterator jet_it = pfJetsFail.begin(); jet_it != pfJetsFail.end(); jet_it++){	
		      if( isSameObject(PhoOne->caloPosition,(*jet_it)->momentum,0.5) || isSameObject(PhoTwo->caloPosition,(*jet_it)->momentum,0.5) ) go2=false;
		    }
		    if(go2){nCnt[8]++;h_ggMet_JetFail->Fill(met->met());h_MET_JetFail->Fill(met->met());}
		  }
		  // Di-Jet Pt
		  if(HasDiJetPt){
		    h_ggDiJetPt->Fill(diJETpt);
		    h_ggDiJetPtOverDiEMPtVsDiEMPt->Fill(diEMpt,diJETpt/diEMpt);
		  }
		  else{h_ggDiJetPt->Fill(diEMpt);}
		  //trailing leg photon for background
		  if(met->met()<30.){
		    h_BgroundCombIsoDR03Nminus3trail->Fill(PhoTwo->ecalRecHitSumEtConeDR03+PhoTwo->hcalTowerSumEtConeDR03()+PhoTwo->trkSumPtHollowConeDR03);
		    h_BgroundCombIsoDR04Nminus3trail->Fill(PhoTwo->ecalRecHitSumEtConeDR04+PhoTwo->hcalTowerSumEtConeDR04()+PhoTwo->trkSumPtHollowConeDR04);
		    h_BgroundEcalIsoDR03Nminus1trail->Fill(PhoTwo->ecalRecHitSumEtConeDR03);
		    h_BgroundHcalIsoDR03Nminus1trail->Fill(PhoTwo->hcalTowerSumEtConeDR03());
		    h_BgroundTrackIsoDR03Nminus1trail->Fill(PhoTwo->trkSumPtHollowConeDR03);
		  }
		  h_r9_gg->Fill(PhoOne->r9);h_r9_gg->Fill(PhoTwo->r9);
		  h_ggPt->Fill(PhoOne->momentum.Pt());h_ggPt->Fill(PhoTwo->momentum.Pt());
		  h_sumEt_gg->Fill(met->sumEt);
		  h_ggMet->Fill(met->met(),PUweight);
		  if(InvMass>110 && InvMass<135)h_ggMETInvarMass110_135->Fill(met->met(),PUweight);
		  if(InvMass>85 && InvMass<110)h_ggMETInvarMass85_110->Fill(met->met(),PUweight);
		  if(InvMass>135 && InvMass<160)h_ggMETInvarMass135_160->Fill(met->met(),PUweight);
		  h_ggSigIetaIeta->Fill(PhoOne->sigmaIetaIeta);h_ggSigIetaIeta->Fill(PhoTwo->sigmaIetaIeta);
		  h_ggDiEMPt->Fill(diEMpt);
		  if(InvMass<40.)h_ggDiEMPtInvMassBelow40->Fill(diEMpt);
		  else h_ggDiEMPtInvMassAbove40->Fill(diEMpt);
		  h_ggInvarMass->Fill(InvMass,PUweight);
		  if(met->met()>30){
		    h_ggInvarMassMET30->Fill(InvMass,PUweight);
		    if(met->met()>40){
		      h_ggInvarMassMET40->Fill(InvMass,PUweight);
		      if(met->met()>50){
			h_ggInvarMassMET50->Fill(InvMass,PUweight);
			if(met->met()>60){
			  h_ggInvarMassMET60->Fill(InvMass,PUweight);
			  if(met->met()>70){
			    h_ggInvarMassMET70->Fill(InvMass,PUweight);
			    if(met->met()>80){
			      h_ggInvarMassMET80->Fill(InvMass,PUweight);
			      if(met->met()>100){
				h_ggInvarMassMET100->Fill(InvMass,PUweight);
			      }
			    }
			  }
			}
		      }
		    }
		  }
		  h_rho_gg->Fill(Rho);h_NVertex_gg->Fill(NVertex);
		  h_EcalIsoDR03_gg->Fill(PhoOne->ecalRecHitSumEtConeDR03);h_EcalIsoDR03_gg->Fill(PhoTwo->ecalRecHitSumEtConeDR03);
		  h_EcalIsoDR03RhoCorr_gg->Fill(PhoOne->ecalRecHitSumEtConeDR03-PUCorr_ECAL*Rho25);h_EcalIsoDR03RhoCorr_gg->Fill(PhoTwo->ecalRecHitSumEtConeDR03-PUCorr_ECAL*Rho25);
		  h_HcalIsoDR03_gg->Fill(PhoOne->hcalTowerSumEtConeDR03());h_HcalIsoDR03_gg->Fill(PhoTwo->hcalTowerSumEtConeDR03());
		  h_HcalIsoDR03RhoCorr_gg->Fill(PhoOne->hcalTowerSumEtConeDR03()-PUCorr_HCAL*Rho25);h_HcalIsoDR03RhoCorr_gg->Fill(PhoTwo->hcalTowerSumEtConeDR03()-PUCorr_HCAL*Rho25);
		  h_TrackIsoDR03_gg->Fill(PhoOne->trkSumPtHollowConeDR03);h_TrackIsoDR03_gg->Fill(PhoTwo->trkSumPtHollowConeDR03);
		  h_TrackIsoDR03RhoCorr_gg->Fill(PhoOne->trkSumPtHollowConeDR03-PUCorr_TRACK*Rho25);h_TrackIsoDR03RhoCorr_gg->Fill(PhoTwo->trkSumPtHollowConeDR03-PUCorr_TRACK*Rho25);
		  h_CombIsoDR03_gg->Fill(PhoOne->ecalRecHitSumEtConeDR03+PhoOne->hcalTowerSumEtConeDR03()+PhoOne->trkSumPtHollowConeDR03);
		  h_CombIsoDR03_gg->Fill(PhoTwo->ecalRecHitSumEtConeDR03+PhoTwo->hcalTowerSumEtConeDR03()+PhoTwo->trkSumPtHollowConeDR03);
		  h_CombIsoDR03RhoCorr_gg->Fill(PhoOne->ecalRecHitSumEtConeDR03-PUCorr_ECAL*Rho25+PhoOne->hcalTowerSumEtConeDR03()-PUCorr_HCAL*Rho25+PhoOne->trkSumPtHollowConeDR03-PUCorr_TRACK*Rho25);
		  h_CombIsoDR03RhoCorr_gg->Fill(PhoTwo->ecalRecHitSumEtConeDR03-PUCorr_ECAL*Rho25+PhoTwo->hcalTowerSumEtConeDR03()-PUCorr_HCAL*Rho25+PhoTwo->trkSumPtHollowConeDR03-PUCorr_TRACK*Rho25);
		  h_DiEMPtVsMet_gg->Fill(met->met(),diEMpt);
		  if(HasDiJetPt){
		    if(JetCounter==0){h_ggDiJetPt_0Jet->Fill(diJETpt);h_ggDiEMPt_0Jet->Fill(diEMpt);}
		    else if(JetCounter==1){h_ggDiJetPt_1Jet->Fill(diJETpt);h_ggDiEMPt_1Jet->Fill(diEMpt);}
		    else if(JetCounter>=2){h_ggDiJetPt_2Jet->Fill(diJETpt);h_ggDiEMPt_2Jet->Fill(diEMpt);}
		    if(JetCounter>2){h_ggDiJetPt_MoreThan2Jets->Fill(diJETpt);}
		  }
		  else{
		    if(JetCounter==0){h_ggDiJetPt_0Jet->Fill(diEMpt);h_ggDiEMPt_0Jet->Fill(diEMpt);}
		    else if(JetCounter==1){h_ggDiJetPt_1Jet->Fill(diEMpt);h_ggDiEMPt_1Jet->Fill(diEMpt);}
		    else if(JetCounter>=2){h_ggDiJetPt_2Jet->Fill(diEMpt);h_ggDiEMPt_2Jet->Fill(diEMpt);}
		    if(JetCounter>2){h_ggDiJetPt_MoreThan2Jets->Fill(diEMpt);}
		    h_ggMet_NoJetMatch->Fill(met->met());
		  }
		  h_ggMR->Fill(RazrMr);
		  h_ggR2->Fill(RazrR2);
		  h_ggR2vsMR->Fill(RazrMr,RazrR2);
		  int n_e=0;
		  for(std::vector<susy::Photon*>::iterator pho_it = pho_Cands.begin();pho_it<pho_Cands.end();pho_it++){
		    if(!isSameObject(PhoOne->caloPosition,(*pho_it)->caloPosition,0.1) 
		       && !isSameObject(PhoTwo->caloPosition,(*pho_it)->caloPosition,0.1) ){
		      if((*pho_it)->nPixelSeeds>0){
			if(n_e==0){
			  ngge++;
			  h_ggePt->Fill(PhoOne->momentum.Pt());h_ggePt->Fill(PhoTwo->momentum.Pt());
			  h_ggePt->Fill((*pho_it)->momentum.Pt());
			  h_ggeElePt->Fill((*pho_it)->momentum.Pt());
			  h_ggeMet->Fill(met->met());
			  h_ggeSigIetaIeta->Fill(PhoOne->sigmaIetaIeta);h_ggeSigIetaIeta->Fill(PhoTwo->sigmaIetaIeta);
			  h_ggeSigIetaIeta->Fill((*pho_it)->sigmaIetaIeta);
			  h_ggeDiEMPt->Fill(diEMpt);
			  h_ggeTriEMPt->Fill(GetTriEmPt(PhoOne,PhoTwo,(*pho_it)));
			  h_ggeInvarMass->Fill(InvMass);
			  if(met->met()>30)h_ggeInvarMassMET30->Fill(InvMass);
			  h_rho_gge->Fill(Rho);h_NVertex_gge->Fill(NVertex);
			  n_e++;
			}
			else if(n_e==1){
			  nggee++;
			  h_ggeePt->Fill(PhoOne->momentum.Pt());h_ggeePt->Fill(PhoTwo->momentum.Pt());
			  h_ggeePt->Fill((*pho_it)->momentum.Pt());
			  h_ggeeElePt->Fill((*pho_it)->momentum.Pt());
			  h_ggeeMet->Fill(met->met());
			  h_ggeeSigIetaIeta->Fill(PhoOne->sigmaIetaIeta);h_ggeeSigIetaIeta->Fill(PhoTwo->sigmaIetaIeta);
			  h_ggeeSigIetaIeta->Fill((*pho_it)->sigmaIetaIeta);
			  h_ggeeDiEMPt->Fill(diEMpt);
			  h_ggeeInvarMass->Fill(InvMass);
			  h_rho_ggee->Fill(Rho);h_NVertex_ggee->Fill(NVertex);
			}
		      }
		    }
		  }//end gge
		  //now ggMuon
		  for(std::vector<susy::Muon*>::iterator mu_it = Muons.begin();mu_it!=Muons.end();mu_it++){
		    h_ggmMet->Fill(met->met());
		    h_ggmInvarMass->Fill(InvMass);
		    if(met->met()>30)h_ggmInvarMassMET30->Fill(InvMass);
		  }//end ggm
		  h_ggdPhi->Fill(dPhi);
		  h_ggAlphaT->Fill(alphaT);
		  if(JetCounter==0)h_ggAlphaT_0Jet->Fill(alphaT);
		  h_ggPhotonLessHt->Fill(PhotonLessHt);
		  h_ggPhotonLessHtVsMET->Fill(met->met(),PhotonLessHt);
		  h_ggMETdPhiLead->Fill(getDphi(PhoOne->caloPosition.Phi(),met->mEt.Phi()));
		  h_ggMETdPhiTrail->Fill(getDphi(PhoTwo->caloPosition.Phi(),met->mEt.Phi()));
		  if(doJetReq){
		    if(printLevel>1)cout<<"Inside gg_JetReq"<<endl;
		    ngg_JetReq++;isgg_JetReq=true;
		    ggPTs_JetReq=make_pair(PhoOne->momentum.Pt(),PhoTwo->momentum.Pt());
		    ggJetReqRunEvent.push_back(make_pair(event->runNumber,event->eventNumber));
		    h_SeedTimeVsEta_gg_JetReq->Fill(PhoOne->caloPosition.Eta(),PhoOne->seedTime);
		    h_SeedTimeVsEta_gg_JetReq->Fill(PhoTwo->caloPosition.Eta(),PhoTwo->seedTime);
		    h_ggMet_JetReq->Fill(met->met(),PUweight);
		    h_ggDiEMPt_JetReq->Fill(diEMpt);		  
		    if(event->isRealData && met->met()>100)cout<<"gg_JetReq Event with MET="<<met->met()<<"  Run: "<<event->runNumber<<"  Lumi: "<<event->luminosityBlockNumber<<"  Event: "<<event->eventNumber<<"  PhoOne pT:"<<PhoOne->momentum.Pt()<<"  PhoTwo pT:"<<PhoTwo->momentum.Pt()<<endl;
		    
		    // Di-Jet Pt
		    if(HasDiJetPt){
		      h_ggDiJetPt_JetReq->Fill(diJETpt);
		      h_ggDiJetPtOverDiEMPtVsDiEMPt_JetReq->Fill(diEMpt,diJETpt/diEMpt);
		    }
		    else{h_ggDiJetPt_JetReq->Fill(diEMpt);}
		    h_ggInvarMass_JetReq->Fill(InvMass,PUweight);
		    if(met->met()>30){
		      h_ggInvarMassMET30_JetReq->Fill(InvMass,PUweight);
		      if(met->met()>40){
			h_ggInvarMassMET40_JetReq->Fill(InvMass,PUweight);
			if(met->met()>50){
			  h_ggInvarMassMET50_JetReq->Fill(InvMass,PUweight);
			  if(met->met()>60){
			    h_ggInvarMassMET60_JetReq->Fill(InvMass,PUweight);
			    if(met->met()>70){
			      h_ggInvarMassMET70_JetReq->Fill(InvMass,PUweight);
			      if(met->met()>80){
				h_ggInvarMassMET80_JetReq->Fill(InvMass,PUweight);
				if(met->met()>100){
				  h_ggInvarMassMET100_JetReq->Fill(InvMass,PUweight);
				}
			      }
			    }
			  }
			}
		      }
		    }
		    if(HasDiJetPt){
		      if(JetCounter==0){h_ggDiJetPt_JetReq_0Jet->Fill(diJETpt);h_ggDiEMPt_JetReq_0Jet->Fill(diEMpt);}
		      else if(JetCounter==1){h_ggDiJetPt_JetReq_1Jet->Fill(diJETpt);h_ggDiEMPt_JetReq_1Jet->Fill(diEMpt);}
		      else if(JetCounter>=2){h_ggDiJetPt_JetReq_2Jet->Fill(diJETpt);h_ggDiEMPt_JetReq_2Jet->Fill(diEMpt);}
		      if(JetCounter>2){h_ggDiJetPt_JetReq_MoreThan2Jets->Fill(diJETpt);}
		    }
		    else{
		      if(JetCounter==0){h_ggDiJetPt_JetReq_0Jet->Fill(diEMpt);h_ggDiEMPt_JetReq_0Jet->Fill(diEMpt);}
		      else if(JetCounter==1){h_ggDiJetPt_JetReq_1Jet->Fill(diEMpt);h_ggDiEMPt_JetReq_1Jet->Fill(diEMpt);}
		      else if(JetCounter>=2){h_ggDiJetPt_JetReq_2Jet->Fill(diEMpt);h_ggDiEMPt_JetReq_2Jet->Fill(diEMpt);}
		      if(JetCounter>2){h_ggDiJetPt_JetReq_MoreThan2Jets->Fill(diEMpt);}
		      h_ggMet_NoJetMatch_JetReq->Fill(met->met());
		    }
		  }
		  //}//seedTime
		}//end gg (nPixelSeeds==0)
		// look for pixel seeds - if both have at least one, -> ee
		else if( PhoOne->nPixelSeeds>0 && PhoTwo->nPixelSeeds>0 ){
		  if(printLevel>1)cout<<"Inside ee"<<endl;
		  isee=true;//need this here for filtering purposes, so sideband ends up in dataset
		  eePTs=make_pair(PhoOne->momentum.Pt(),PhoTwo->momentum.Pt());
		  //Require 81<InvariantMass<101
		  h_eeInvarMassFullRange->Fill(InvMass);
		  if(PhoTwo->momentum.Pt()>25 && PhoTwo->momentum.Pt()<40)h_eeInvarMassFullRange_Pt25to40->Fill(InvMass);
		  if(PhoOne->momentum.Pt()>40 && PhoOne->momentum.Pt()<45)h_eeInvarMassFullRange_Pt40to45->Fill(InvMass);
		  else if(PhoOne->momentum.Pt()>45 && PhoOne->momentum.Pt()<50)h_eeInvarMassFullRange_Pt45to50->Fill(InvMass);
		  else if(PhoOne->momentum.Pt()>50 && PhoOne->momentum.Pt()<69)h_eeInvarMassFullRange_Pt50to60->Fill(InvMass);
		  else if(PhoOne->momentum.Pt()>60 && PhoOne->momentum.Pt()<80)h_eeInvarMassFullRange_Pt60to80->Fill(InvMass);
		  else if(PhoOne->momentum.Pt()>80)h_eeInvarMassFullRange_Pt80->Fill(InvMass);
		  h_sumEt_eeFullRange->Fill(met->sumEt);
		  if(InvMass<40.)h_eeDiEMPtInvMassBelow40->Fill(diEMpt);
		  else h_eeDiEMPtInvMassAbove40->Fill(diEMpt);
		  if(JetCounter>=2){
		    h_eeInvarMassFullRange_2JetReq->Fill(InvMass);
		    if( InvMass > 81. && InvMass < 101.){
		      nee_2JetReq++;
		      h_eeMet_2JetReq->Fill(met->met());
		      h_eeDiEMPt_2JetReq->Fill(diEMpt);
		      h_eeInvarMass_2JetReq->Fill(InvMass);
		      if(printLevel>1)cout<<"Line1675"<<endl;
		      h_eeMet_reweight_binned_2JetReq->Fill(met->met(),GetMetReweight(diEMpt,"ee",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,-1).first);
		      if(printLevel>1)cout<<"Line1677"<<endl;
		    }
		  }
		  if( InvMass > 81. && InvMass < 101.){
		  if(printLevel>1)cout<<"Inside ee  Z peak"<<endl;
		    //require seed time within +-3ns for no jet requirement case
		    //if(fabs(PhoOne->seedTime)<3. && fabs(PhoTwo->seedTime)<3.){
		    nee++;
		    if(fake_Cands.size()>0)neef++;
		    eeRunEvent.push_back(make_pair(event->runNumber,event->eventNumber));
		    h_SeedTimeVsEta_ee->Fill(PhoOne->caloPosition.Eta(),PhoOne->seedTime);
		    h_SeedTimeVsEta_ee->Fill(PhoTwo->caloPosition.Eta(),PhoTwo->seedTime);
		    h_MetVsSeedTime_ee->Fill(PhoOne->seedTime, met->met());
		    h_MetVsSeedTime_ee->Fill(PhoTwo->seedTime, met->met());
		    h_NVertexVsMET_ee->Fill(met->met(),NVertex);
		    h_eePerInstLumi->Fill(event->avgInsRecLumi);
		    if(isJetFail){bool go2=true;
		      for(std::vector<susy::PFJet*>::iterator jet_it = pfJetsFail.begin(); jet_it != pfJetsFail.end(); jet_it++){	
			if( isSameObject(PhoOne->caloPosition,(*jet_it)->momentum,0.5) || isSameObject(PhoTwo->caloPosition,(*jet_it)->momentum,0.5) ) go2=false;
		      }
		      if(go2){
			nCnt[9]++;h_MET_JetFail->Fill(met->met());
			h_eeMet_JetFail->Fill(met->met());
			if(HasDiJetPt)h_eeDiJetPt_JetFail->Fill(diJETpt);
			else h_eeDiJetPt_JetFail->Fill(diEMpt);
		      }
		    }
		    // Di-Jet Pt
		    if(HasDiJetPt){
		      h_eeDiJetPt->Fill(diJETpt);
		      h_eeDiJetPtOverDiEMPtVsDiEMPt->Fill(diEMpt,diJETpt/diEMpt);
		    }
		    else{h_eeDiJetPt->Fill(diEMpt);}
		    h_r9_ee->Fill(PhoOne->r9);h_r9_ee->Fill(PhoTwo->r9);
		    h_eePt->Fill(PhoOne->momentum.Pt());h_eePt->Fill(PhoTwo->momentum.Pt());
		    h_sumEt_ee->Fill(met->sumEt);
		    h_eeMet->Fill(met->met());
		    h_eeMet_reweight_binned->Fill(met->met(),GetMetReweight(diEMpt,"ee",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		    if(HasDiJetPt){
		      h_eeMet_reweightJet_binned->Fill(met->met(),GetMetReweight(diJETpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		      h_eeDiJetPt_reweight_binned->Fill(diJETpt,GetMetReweight(diJETpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		      if(JetCounter==0){h_eeDiJetPt_0Jet->Fill(diJETpt);h_eeDiEMPt_0Jet->Fill(diEMpt);}
		      else if(JetCounter==1){h_eeDiJetPt_1Jet->Fill(diJETpt);h_eeDiEMPt_1Jet->Fill(diEMpt);}
		      else if(JetCounter>=2){h_eeDiJetPt_2Jet->Fill(diJETpt);h_eeDiEMPt_2Jet->Fill(diEMpt);}
		      if(JetCounter>2){h_eeDiJetPt_MoreThan2Jets->Fill(diJETpt);}
		      if(NVertex<10)h_eeMet_NV0_10->Fill(met->met(),GetMetReweight(diJETpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		      else if(NVertex>=10 && NVertex<15)h_eeMet_NV10_15->Fill(met->met(),GetMetReweight(diJETpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		      else if(NVertex>=15)h_eeMet_NV15up->Fill(met->met(),GetMetReweight(diJETpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		    }
		    else{
		      h_eeMet_reweightJet_binned->Fill(met->met(),GetMetReweight(diEMpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		      h_eeDiJetPt_reweight_binned->Fill(diEMpt,GetMetReweight(diEMpt,"ee",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		      if(JetCounter==0){h_eeDiJetPt_0Jet->Fill(diEMpt);h_eeDiEMPt_0Jet->Fill(diEMpt);}
		      else if(JetCounter==1){h_eeDiJetPt_1Jet->Fill(diEMpt);h_eeDiEMPt_1Jet->Fill(diEMpt);}
		      else if(JetCounter>=2){h_eeDiJetPt_2Jet->Fill(diEMpt);h_eeDiEMPt_2Jet->Fill(diEMpt);}
		      if(JetCounter>2){h_eeDiJetPt_MoreThan2Jets->Fill(diEMpt);}
		      if(NVertex<10)h_eeMet_NV0_10->Fill(met->met(),GetMetReweight(diEMpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		      else if(NVertex>=10 && NVertex<15)h_eeMet_NV10_15->Fill(met->met(),GetMetReweight(diEMpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		      else if(NVertex>=15)h_eeMet_NV15up->Fill(met->met(),GetMetReweight(diEMpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		      h_eeMet_NoJetMatch->Fill(met->met());
		    }
		    h_eeMR->Fill(RazrMr);
		    h_eeR2->Fill(RazrR2);
		    h_eeR2vsMR->Fill(RazrMr,RazrR2);
		    h_eeMR_reweight_binned->Fill(RazrMr,GetMetReweight(diEMpt,"ee",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,-1).first);
		    h_eeR2_reweight_binned->Fill(RazrR2,GetMetReweight(diEMpt,"ee",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,-1).first);
		    h_eeR2vsMR_reweight_binned->Fill(RazrMr,RazrR2,GetMetReweight(diEMpt,"ee",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,-1).first);
		    toys->cd();
		    if(printLevel>1)cout<<"Inside ee  Z peak Toys"<<endl;
		    if(HasDiJetPt){
		      for(int i=0;i<1000;i++){
			//title="eeMet_reweight_binned_toy";title+=i;
			eeMet_reweight_binned_toy[i]->Fill(met->met(),rando.Gaus(GetMetReweight(diJETpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first,GetMetReweight(diJETpt,"ee",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).second));
			//eeMet_reweight_binned_toy[i]->Fill(met->met(),rando.Poisson(GetMetReweight(diEMpt,"ee",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,-1).first));
		      }
		    }
		    else{
		      for(int i=0;i<1000;i++){
			//title="eeMet_reweight_binned_toy";title+=i;
			eeMet_reweight_binned_toy[i]->Fill(met->met(),rando.Gaus(GetMetReweight(diEMpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first,GetMetReweight(diEMpt,"ee",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).second));
			//eeMet_reweight_binned_toy[i]->Fill(met->met(),rando.Poisson(GetMetReweight(diEMpt,"ee",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,-1).first));
		      }
		    }
		    if(printLevel>1)cout<<"Outside ee  Z peak Toys"<<endl;
		    fout->cd();
		    h_eeSigIetaIeta->Fill(PhoOne->sigmaIetaIeta);h_eeSigIetaIeta->Fill(PhoTwo->sigmaIetaIeta);
		    h_eeDiEMPt->Fill(diEMpt);
		    h_eeDiEMPt_reweight_binned->Fill(diEMpt,GetMetReweight(diEMpt,"ee",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		    h_eePixelSeeds->Fill(PhoOne->nPixelSeeds);h_eePixelSeeds->Fill(PhoTwo->nPixelSeeds);
		    h_eeInvarMass->Fill(InvMass);
		    h_rho_ee->Fill(Rho);h_NVertex_ee->Fill(NVertex);
		    h_EcalIsoDR03RhoCorr_ee->Fill(PhoOne->ecalRecHitSumEtConeDR03-PUCorr_ECAL*Rho25);h_EcalIsoDR03RhoCorr_ee->Fill(PhoTwo->ecalRecHitSumEtConeDR03-PUCorr_ECAL*Rho25);
		    h_HcalIsoDR03RhoCorr_ee->Fill(PhoOne->hcalTowerSumEtConeDR03()-PUCorr_HCAL*Rho25);h_HcalIsoDR03RhoCorr_ee->Fill(PhoTwo->hcalTowerSumEtConeDR03()-PUCorr_HCAL*Rho25);
		    h_TrackIsoDR03RhoCorr_ee->Fill(PhoOne->trkSumPtHollowConeDR03-PUCorr_TRACK*Rho25);h_TrackIsoDR03RhoCorr_ee->Fill(PhoTwo->trkSumPtHollowConeDR03-PUCorr_TRACK*Rho25);
		    h_CombIsoDR03RhoCorr_ee->Fill(PhoOne->ecalRecHitSumEtConeDR03-PUCorr_ECAL*Rho25+PhoOne->hcalTowerSumEtConeDR03()-PUCorr_HCAL*Rho25+PhoOne->trkSumPtHollowConeDR03-PUCorr_TRACK*Rho25);
		    h_CombIsoDR03RhoCorr_ee->Fill(PhoTwo->ecalRecHitSumEtConeDR03-PUCorr_ECAL*Rho25+PhoTwo->hcalTowerSumEtConeDR03()-PUCorr_HCAL*Rho25+PhoTwo->trkSumPtHollowConeDR03-PUCorr_TRACK*Rho25);
		    h_DiEMPtVsMet_ee->Fill(met->met(),diEMpt);
		    int n_g=0;
		    for(std::vector<susy::Photon*>::iterator pho_it = pho_Cands.begin();pho_it<pho_Cands.end();pho_it++){
		      if(!isSameObject(PhoOne->caloPosition,(*pho_it)->caloPosition,0.1) 
			 && !isSameObject(PhoTwo->caloPosition,(*pho_it)->caloPosition,0.1) ){
			if((*pho_it)->nPixelSeeds==0){
			  if(n_g==0){
			    neeg++;
			    h_eegPt->Fill(PhoOne->momentum.Pt());h_eegPt->Fill(PhoTwo->momentum.Pt());
			    h_eegPt->Fill((*pho_it)->momentum.Pt());
			    h_eegElePt->Fill((*pho_it)->momentum.Pt());
			    h_eegMet->Fill(met->met());
			    h_eegSigIetaIeta->Fill(PhoOne->sigmaIetaIeta);h_eegSigIetaIeta->Fill(PhoTwo->sigmaIetaIeta);
			    h_eegSigIetaIeta->Fill((*pho_it)->sigmaIetaIeta);
			    h_eegDiEMPt->Fill(diEMpt);
			    h_eegTriEMPt->Fill(GetTriEmPt(PhoOne,PhoTwo,(*pho_it)));
			    h_eegInvarMass->Fill(InvMass);
			    h_rho_eeg->Fill(Rho);h_NVertex_eeg->Fill(NVertex);
			    n_g++;
			  }
			  else if(n_g==1){
			    neegg++;
			    h_eeggPt->Fill(PhoOne->momentum.Pt());h_eeggPt->Fill(PhoTwo->momentum.Pt());
			    h_eeggPt->Fill((*pho_it)->momentum.Pt());
			    h_eeggElePt->Fill((*pho_it)->momentum.Pt());
			    h_eeggMet->Fill(met->met());
			    h_eeggSigIetaIeta->Fill(PhoOne->sigmaIetaIeta);h_eeggSigIetaIeta->Fill(PhoTwo->sigmaIetaIeta);
			    h_eeggSigIetaIeta->Fill((*pho_it)->sigmaIetaIeta);
			    h_eeggDiEMPt->Fill(diEMpt);
			    h_eeggInvarMass->Fill(InvMass);
			    h_rho_eegg->Fill(Rho);h_NVertex_eegg->Fill(NVertex);
			  }
			}
		      }
		    }//end eeg
		    h_eedPhi->Fill(dPhi);
		    h_eeAlphaT->Fill(alphaT);
		    if(JetCounter==0)h_eeAlphaT_0Jet->Fill(alphaT);
		    h_eePhotonLessHt->Fill(PhotonLessHt);
		    h_eePhotonLessHtVsMET->Fill(met->met(),PhotonLessHt);
		    h_eeMETdPhiLead->Fill(getDphi(PhoOne->caloPosition.Phi(),met->mEt.Phi()));
		    h_eeMETdPhiTrail->Fill(getDphi(PhoTwo->caloPosition.Phi(),met->mEt.Phi()));
		    //}//seedTime
		  }//81<InvMass<101
		  //now do sideband
		  if(   (InvMass > 71.  && InvMass < 81.)
			|| (InvMass > 101. && InvMass < 111.) ){
		    neeSideBand++;
		    h_eeSideBandDiEMPt->Fill(diEMpt);
		    h_eeSideBandMet->Fill(met->met());
		    h_eeSidebandMet_reweight_binned->Fill(met->met(),GetMetReweight(diEMpt,"ee",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,-1).first);
		    h_DiEMPtVsMet_eeSideBand->Fill(met->met(),diEMpt);
		    h_eeSideBandMR->Fill(RazrMr);
		    h_eeSideBandR2->Fill(RazrR2);
		    h_eeSideBandR2vsMR->Fill(RazrMr,RazrR2);
		    h_eeSideBandMR_reweight_binned->Fill(RazrMr,GetMetReweight(diEMpt,"ee",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,-1).first);
		    h_eeSideBandR2_reweight_binned->Fill(RazrR2,GetMetReweight(diEMpt,"ee",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,-1).first);
		    h_eeSideBandR2vsMR_reweight_binned->Fill(RazrMr,RazrR2,GetMetReweight(diEMpt,"ee",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,-1).first);
		    if(HasDiJetPt){
		      h_eeSideBandDiJetPt->Fill(diJETpt);
		      if(NVertex<10)h_eeSBMet_NV0_10->Fill(met->met(),GetMetReweight(diJETpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		      else if(NVertex>=10 && NVertex<15)h_eeSBMet_NV10_15->Fill(met->met(),GetMetReweight(diJETpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		      else if(NVertex>=15)h_eeSBMet_NV15up->Fill(met->met(),GetMetReweight(diJETpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		      if(InvMass > 71.  && InvMass < 81.){
			h_eeSidebandLowMet_reweightJet_binned->Fill(met->met(),GetMetReweight(diJETpt,"eeSbLow",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
			h_eeSideBandLowDiEMPt->Fill(diEMpt);h_eeSideBandLowDiJetPt->Fill(diJETpt);
			if(JetCounter==0){h_eeSideBandLowDiEMPt_0Jet->Fill(diEMpt);h_eeSideBandLowDiJetPt_0Jet->Fill(diJETpt);}
			else if(JetCounter==1){h_eeSideBandLowDiEMPt_1Jet->Fill(diEMpt);h_eeSideBandLowDiJetPt_1Jet->Fill(diJETpt);}
			else if(JetCounter>=2){h_eeSideBandLowDiEMPt_2Jet->Fill(diEMpt);h_eeSideBandLowDiJetPt_2Jet->Fill(diJETpt);}
			toys->cd();
			for(int i=0;i<1000;i++){
			  eeSidebandLowMet_reweight_binned_toy[i]->Fill(met->met(),rando.Gaus(GetMetReweight(diJETpt,"eeSbLow",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first,GetMetReweight(diJETpt,"eeSbLow",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).second));
			}
			fout->cd();
		      }
		      if(InvMass > 101.  && InvMass < 111.){
			if(printLevel>1)cout<<"Inside ee SBhigh"<<endl;
			h_eeSidebandHighMet_reweightJet_binned->Fill(met->met(),GetMetReweight(diJETpt,"eeSbHigh",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
			h_eeSidebandHighMet_reweight_binned->Fill(met->met(),GetMetReweight(diEMpt,"eeSbHigh",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
			h_eeSideBandHighDiEMPt->Fill(diEMpt);h_eeSideBandHighDiJetPt->Fill(diJETpt);
			if(JetCounter==0){h_eeSideBandHighDiEMPt_0Jet->Fill(diEMpt);h_eeSideBandHighDiJetPt_0Jet->Fill(diJETpt);}
			else if(JetCounter==1){h_eeSideBandHighDiEMPt_1Jet->Fill(diEMpt);h_eeSideBandHighDiJetPt_1Jet->Fill(diJETpt);}
			else if(JetCounter>=2){h_eeSideBandHighDiEMPt_2Jet->Fill(diEMpt);h_eeSideBandHighDiJetPt_2Jet->Fill(diJETpt);}
			toys->cd();
			for(int i=0;i<1000;i++){
			  eeSidebandHighMet_reweight_binned_toy[i]->Fill(met->met(),rando.Gaus(GetMetReweight(diJETpt,"eeSbHigh",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first,GetMetReweight(diJETpt,"eeSbHigh",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).second));
			}
			fout->cd();
		      }
		    }//end if(HasDiJet)
		    else{ 
		      h_eeSideBandDiJetPt->Fill(diEMpt);
		      if(NVertex<10)h_eeSBMet_NV0_10->Fill(met->met(),GetMetReweight(diEMpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		      else if(NVertex>=10 && NVertex<15)h_eeSBMet_NV10_15->Fill(met->met(),GetMetReweight(diEMpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		      else if(NVertex>=15)h_eeSBMet_NV15up->Fill(met->met(),GetMetReweight(diEMpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		      if(InvMass > 71.  && InvMass < 81.){
			h_eeSidebandLowMet_reweightJet_binned->Fill(met->met(),GetMetReweight(diEMpt,"eeSbLow",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
			h_eeSideBandLowDiEMPt->Fill(diEMpt);h_eeSideBandLowDiJetPt->Fill(diEMpt);
			if(JetCounter==0){h_eeSideBandLowDiEMPt_0Jet->Fill(diEMpt);h_eeSideBandLowDiJetPt_0Jet->Fill(diEMpt);}
			else if(JetCounter==1){h_eeSideBandLowDiEMPt_1Jet->Fill(diEMpt);h_eeSideBandLowDiJetPt_1Jet->Fill(diEMpt);}
			else if(JetCounter>=2){h_eeSideBandLowDiEMPt_2Jet->Fill(diEMpt);h_eeSideBandLowDiJetPt_2Jet->Fill(diEMpt);}
			toys->cd();
		  if(printLevel>1)cout<<"Inside ee SBhigh Toys"<<endl;
			for(int i=0;i<1000;i++){
			  eeSidebandLowMet_reweight_binned_toy[i]->Fill(met->met(),rando.Gaus(GetMetReweight(diEMpt,"eeSbLow",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first,GetMetReweight(diEMpt,"eeSbLow",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).second));
			}
			fout->cd();
		      }
		      if(InvMass > 101.  && InvMass < 111.){
		  if(printLevel>1)cout<<"Inside ee SBlow"<<endl;
			h_eeSidebandLowMet_reweightJet_binned->Fill(met->met(),GetMetReweight(diEMpt,"eeSbHigh",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
			h_eeSideBandHighDiEMPt->Fill(diEMpt);h_eeSideBandHighDiJetPt->Fill(diEMpt);
			if(JetCounter==0){h_eeSideBandHighDiEMPt_0Jet->Fill(diEMpt);h_eeSideBandHighDiJetPt_0Jet->Fill(diEMpt);}
			else if(JetCounter==1){h_eeSideBandHighDiEMPt_1Jet->Fill(diEMpt);h_eeSideBandHighDiJetPt_1Jet->Fill(diEMpt);}
			else if(JetCounter>=2){h_eeSideBandHighDiEMPt_2Jet->Fill(diEMpt);h_eeSideBandHighDiJetPt_2Jet->Fill(diEMpt);}
			toys->cd();	
			if(printLevel>1)cout<<"Inside ee SBlow Toys"<<endl;
			for(int i=0;i<1000;i++){
			  eeSidebandHighMet_reweight_binned_toy[i]->Fill(met->met(),rando.Gaus(GetMetReweight(diEMpt,"eeSbHigh",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first,GetMetReweight(diEMpt,"eeSbHigh",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).second));
			}
			fout->cd();
		      }

		    }//end else
		  }//end sideband
		  if(printLevel>1)cout<<"Finished ee, right before ee_JetReq"<<endl;
		  //now jetreq
		  if(doJetReq){
		    isee_JetReq=true;if(printLevel>1)cout<<"Inside ee_JetReq"<<endl;
		    eePTs_JetReq=make_pair(PhoOne->momentum.Pt(),PhoTwo->momentum.Pt());
		    //Require 81<InvariantMass<101
		    h_eeInvarMassFullRange_JetReq->Fill(InvMass);
		    if( InvMass > 81. && InvMass < 101.){
		      nee_JetReq++;
		      eeJetReqRunEvent.push_back(make_pair(event->runNumber,event->eventNumber));
		      //if(pfJets.size()>=2){nee_2JetReq++;}
		      h_eeMet_JetReq->Fill(met->met());
		      h_SeedTimeVsEta_ee_JetReq->Fill(PhoOne->caloPosition.Eta(),PhoOne->seedTime);
		      h_SeedTimeVsEta_ee_JetReq->Fill(PhoTwo->caloPosition.Eta(),PhoTwo->seedTime);
		      h_eeInvarMass_JetReq->Fill(InvMass);
		      h_eeMet_reweight_binned_JetReq->Fill(met->met(),GetMetReweight(diEMpt,"ee",EEreweights_JetReq,FFreweights_JetReq,EESideBandLowReweightsJet_JetReq,EESideBandHighReweightsJet_JetReq,-1).first);
		      if(HasDiJetPt){
			h_eeMet_reweightJet_binned_JetReq->Fill(met->met(),GetMetReweight(diJETpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
			h_eeDiJetPt_reweight_binned_JetReq->Fill(diJETpt,GetMetReweight(diJETpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
			if(JetCounter==0){h_eeDiJetPt_JetReq_0Jet->Fill(diJETpt);h_eeDiEMPt_JetReq_0Jet->Fill(diEMpt);}
			else if(JetCounter==1){h_eeDiJetPt_JetReq_1Jet->Fill(diJETpt);h_eeDiEMPt_JetReq_1Jet->Fill(diEMpt);}
			else if(JetCounter>=2){h_eeDiJetPt_JetReq_2Jet->Fill(diJETpt);h_eeDiEMPt_JetReq_2Jet->Fill(diEMpt);}
			if(JetCounter>2){h_eeDiJetPt_JetReq_MoreThan2Jets->Fill(diJETpt);}
		      }
		      else{
			h_eeMet_reweightJet_binned_JetReq->Fill(met->met(),GetMetReweight(diEMpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
			h_eeDiJetPt_reweight_binned_JetReq->Fill(diEMpt,GetMetReweight(diEMpt,"ee",EEreweights_JetReq,FFreweights_JetReq,EESideBandLowReweightsJet_JetReq,EESideBandHighReweightsJet_JetReq,JetCounter).first);
			if(JetCounter==0){h_eeDiJetPt_JetReq_0Jet->Fill(diEMpt);h_eeDiEMPt_JetReq_0Jet->Fill(diEMpt);}
			else if(JetCounter==1){h_eeDiJetPt_JetReq_1Jet->Fill(diEMpt);h_eeDiEMPt_JetReq_1Jet->Fill(diEMpt);}
			else if(JetCounter>=2){h_eeDiJetPt_JetReq_2Jet->Fill(diEMpt);h_eeDiEMPt_JetReq_2Jet->Fill(diEMpt);}
			else if(JetCounter>2){h_eeDiJetPt_JetReq_MoreThan2Jets->Fill(diEMpt);}
			h_eeMet_NoJetMatch_JetReq->Fill(met->met());
		      }
		      h_eeDiEMPt_JetReq->Fill(diEMpt);
		      // Di-Jet Pt
		      if(HasDiJetPt){
			h_eeDiJetPt_JetReq->Fill(diJETpt);
			h_eeDiJetPtOverDiEMPtVsDiEMPt_JetReq->Fill(diEMpt,diJETpt/diEMpt);
		      }
		      else{h_eeDiJetPt_JetReq->Fill(diEMpt);}
		      h_eePixelSeeds_JetReq->Fill(PhoOne->nPixelSeeds);h_eePixelSeeds_JetReq->Fill(PhoTwo->nPixelSeeds);
		      //fill toys
		      toys->cd();
		      if(HasDiJetPt){
			for(int i=0;i<1000;i++){
			  //title="eeMet_reweight_binned_toy";title+=i;
			  eeMet_reweight_binned_toy_JetReq[i]->Fill(met->met(),rando.Gaus(GetMetReweight(diJETpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first,GetMetReweight(diJETpt,"ee",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).second));
			  //eeMet_reweight_binned_toy[i]->Fill(met->met(),rando.Poisson(GetMetReweight(diEMpt,"ee",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,-1).first));
			}
		      }
		      else{
			for(int i=0;i<1000;i++){
			  //title="eeMet_reweight_binned_toy";title+=i;
			  eeMet_reweight_binned_toy_JetReq[i]->Fill(met->met(),rando.Gaus(GetMetReweight(diEMpt,"ee",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first,GetMetReweight(diEMpt,"ee",EEreweights_JetReq,FFreweights_JetReq,EESideBandLowReweightsJet_JetReq,EESideBandHighReweightsJet_JetReq,JetCounter).second));
			  //eeMet_reweight_binned_toy[i]->Fill(met->met(),rando.Poisson(GetMetReweight(diEMpt,"ee",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,-1).first));
			}
		      }
		      fout->cd();
		    }//end if(81<InvMass<101)
		    //now do sideband
		    if(   (InvMass > 71.  && InvMass < 81.)
			  || (InvMass > 101. && InvMass < 111.) ){
		      neeSideBand_JetReq++;
		      h_eeSideBandDiEMPt_JetReq->Fill(diEMpt);
		      h_eeSideBandMet_JetReq->Fill(met->met());
		      h_eeSidebandMet_reweight_binned_JetReq->Fill(met->met(),GetMetReweight(diEMpt,"ee",EEreweights_JetReq,FFreweights_JetReq,EESideBandLowReweightsJet_JetReq,EESideBandHighReweightsJet_JetReq,-1).first);
		      if(HasDiJetPt){
			h_eeSideBandDiJetPt_JetReq->Fill(diJETpt);
			if(InvMass > 71.  && InvMass < 81.){
			  h_eeSidebandLowMet_reweightJet_binned_JetReq->Fill(met->met(),GetMetReweight(diJETpt,"eeSbLow",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
			  h_eeSideBandLowDiEMPt_JetReq->Fill(diEMpt);h_eeSideBandLowDiJetPt_JetReq->Fill(diJETpt);
			  if(JetCounter==0){h_eeSideBandLowDiEMPt_JetReq_0Jet->Fill(diEMpt);h_eeSideBandLowDiJetPt_JetReq_0Jet->Fill(diJETpt);}
			  else if(JetCounter==1){h_eeSideBandLowDiEMPt_JetReq_1Jet->Fill(diEMpt);h_eeSideBandLowDiJetPt_JetReq_1Jet->Fill(diJETpt);}
			  else if(JetCounter>=2){h_eeSideBandLowDiEMPt_JetReq_2Jet->Fill(diEMpt);h_eeSideBandLowDiJetPt_JetReq_2Jet->Fill(diJETpt);}
			  toys->cd();
			  for(int i=0;i<1000;i++){
			    eeSidebandLowMet_reweight_binned_toy_JetReq[i]->Fill(met->met(),rando.Gaus(GetMetReweight(diJETpt,"eeSbLow",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first,GetMetReweight(diJETpt,"eeSbLow",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).second));
			  }
			  fout->cd();
			}
			else if(InvMass > 101.  && InvMass < 111.){
			  h_eeSidebandHighMet_reweightJet_binned_JetReq->Fill(met->met(),GetMetReweight(diJETpt,"eeSbHigh",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
			  h_eeSideBandHighDiEMPt_JetReq->Fill(diEMpt);h_eeSideBandHighDiJetPt_JetReq->Fill(diJETpt);
			  if(JetCounter==0){h_eeSideBandHighDiEMPt_JetReq_0Jet->Fill(diEMpt);h_eeSideBandHighDiJetPt_JetReq_0Jet->Fill(diJETpt);}
			  else if(JetCounter==1){h_eeSideBandHighDiEMPt_JetReq_1Jet->Fill(diEMpt);h_eeSideBandHighDiJetPt_JetReq_1Jet->Fill(diJETpt);}
			  else if(JetCounter>=2){h_eeSideBandHighDiEMPt_JetReq_2Jet->Fill(diEMpt);h_eeSideBandHighDiJetPt_JetReq_2Jet->Fill(diJETpt);}
			  toys->cd();
			  for(int i=0;i<1000;i++){
			    eeSidebandHighMet_reweight_binned_toy_JetReq[i]->Fill(met->met(),rando.Gaus(GetMetReweight(diJETpt,"eeSbHigh",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first,GetMetReweight(diJETpt,"eeSbHigh",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).second));
			  }
			  fout->cd();
			}
		      }//end if(HasDiJet)
		      else{ 
			h_eeSideBandDiJetPt_JetReq->Fill(diEMpt);
			if(InvMass > 71.  && InvMass < 81.){
			  h_eeSidebandLowMet_reweightJet_binned_JetReq->Fill(met->met(),GetMetReweight(diEMpt,"eeSbLow",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
			  h_eeSideBandLowDiEMPt_JetReq->Fill(diEMpt);h_eeSideBandLowDiJetPt_JetReq->Fill(diEMpt);
			  if(JetCounter==0){h_eeSideBandLowDiEMPt_JetReq_0Jet->Fill(diEMpt);h_eeSideBandLowDiJetPt_JetReq_0Jet->Fill(diEMpt);}
			  else if(JetCounter==1){h_eeSideBandLowDiEMPt_JetReq_1Jet->Fill(diEMpt);h_eeSideBandLowDiJetPt_JetReq_1Jet->Fill(diEMpt);}
			  else if(JetCounter>=2){h_eeSideBandLowDiEMPt_JetReq_2Jet->Fill(diEMpt);h_eeSideBandLowDiJetPt_JetReq_2Jet->Fill(diEMpt);}
			  toys->cd();
			  for(int i=0;i<1000;i++){
			    eeSidebandLowMet_reweight_binned_toy_JetReq[i]->Fill(met->met(),rando.Gaus(GetMetReweight(diEMpt,"eeSbLow",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first,GetMetReweight(diEMpt,"eeSbLow",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).second));
			  }
			  fout->cd();
			}
			if(InvMass > 101.  && InvMass < 111.){
			  h_eeSidebandLowMet_reweightJet_binned->Fill(met->met(),GetMetReweight(diEMpt,"eeSbHigh",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
			  h_eeSideBandHighDiEMPt->Fill(diEMpt);h_eeSideBandHighDiJetPt->Fill(diEMpt);
			  if(JetCounter==0){h_eeSideBandHighDiEMPt_0Jet->Fill(diEMpt);h_eeSideBandHighDiJetPt_0Jet->Fill(diEMpt);}
			  else if(JetCounter==1){h_eeSideBandHighDiEMPt_1Jet->Fill(diEMpt);h_eeSideBandHighDiJetPt_1Jet->Fill(diEMpt);}
			  else if(JetCounter>=2){h_eeSideBandHighDiEMPt_2Jet->Fill(diEMpt);h_eeSideBandHighDiJetPt_2Jet->Fill(diEMpt);}
			  toys->cd();
			  for(int i=0;i<1000;i++){
			    eeSidebandHighMet_reweight_binned_toy[i]->Fill(met->met(),rando.Gaus(GetMetReweight(diEMpt,"eeSbHigh",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first,GetMetReweight(diEMpt,"eeSbHigh",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).second));
			  }
			  fout->cd();
			}
		      }//end else
		    }//end sideband
		  }//end ee jet req case		  
		}//end ee (nPixelSeeds>0)
		// look for pixel seeds - if one has, -> eg
		else if( (PhoOne->nPixelSeeds==0 && PhoTwo->nPixelSeeds>0)
			 || (PhoOne->nPixelSeeds>0 && PhoTwo->nPixelSeeds==0) ){
		  if(printLevel>1)cout<<"Inside eg"<<endl;
		  //require seed time within +-3ns for no jet requirement case
		  //if(fabs(PhoOne->seedTime)<3. && fabs(PhoTwo->seedTime)<3.){
		  neg++;iseg=true;
		  egPTs=make_pair(PhoOne->momentum.Pt(),PhoTwo->momentum.Pt());
		  egRunEvent.push_back(make_pair(event->runNumber,event->eventNumber));
		  h_SeedTimeVsEta_eg->Fill(PhoOne->caloPosition.Eta(),PhoOne->seedTime);
		  h_SeedTimeVsEta_eg->Fill(PhoTwo->caloPosition.Eta(),PhoTwo->seedTime);
		  h_NVertexVsMET_eg->Fill(met->met(),NVertex);
		  h_egPerInstLumi->Fill(event->avgInsRecLumi);
		  if(NVertex<10)h_egMet_NV0_10->Fill(met->met());
		  else if(NVertex>=10 && NVertex<15)h_egMet_NV10_15->Fill(met->met());
		  else if(NVertex>=15)h_egMet_NV15up->Fill(met->met());
		  if(isJetFail){bool go2=true;
		    for(std::vector<susy::PFJet*>::iterator jet_it = pfJetsFail.begin(); jet_it != pfJetsFail.end(); jet_it++){	
		      if( isSameObject(PhoOne->caloPosition,(*jet_it)->momentum,0.5) || isSameObject(PhoTwo->caloPosition,(*jet_it)->momentum,0.5) ) go2=false;
		    }
		    if(go2){nCnt[10]++;h_egMet_JetFail->Fill(met->met());h_MET_JetFail->Fill(met->met());}
		  }
		  h_r9_eg->Fill(PhoOne->r9);h_r9_eg->Fill(PhoTwo->r9);
		  h_egPt->Fill(PhoOne->momentum.Pt());h_egPt->Fill(PhoTwo->momentum.Pt());
		  h_sumEt_eg->Fill(met->sumEt);
		  h_egMet->Fill(met->met());
		  h_egSigIetaIeta->Fill(PhoOne->sigmaIetaIeta);h_egSigIetaIeta->Fill(PhoTwo->sigmaIetaIeta);
		  h_egDiEMPt->Fill(diEMpt);
		  if(InvMass<40.)h_egDiEMPtInvMassBelow40->Fill(diEMpt);
		  else h_egDiEMPtInvMassAbove40->Fill(diEMpt);
		  h_egPixelSeeds->Fill(PhoOne->nPixelSeeds);h_egPixelSeeds->Fill(PhoTwo->nPixelSeeds);
		  h_egInvarMass->Fill(InvMass);
		  if(met->met()>30)h_egInvarMassMET30->Fill(InvMass);
		  if(PhoTwo->nPixelSeeds==0 && PhoTwo->momentum.Pt()>25 && PhoTwo->momentum.Pt()<40)h_egInvarMass_Pt25to40->Fill(InvMass);
		  if(PhoOne->nPixelSeeds==0 && PhoOne->momentum.Pt()>40 && PhoOne->momentum.Pt()<45)h_egInvarMass_Pt40to45->Fill(InvMass);
		  else if(PhoOne->nPixelSeeds==0 && PhoOne->momentum.Pt()>45 && PhoOne->momentum.Pt()<50)h_egInvarMass_Pt45to50->Fill(InvMass);
		  else if(PhoOne->nPixelSeeds==0 && PhoOne->momentum.Pt()>50 && PhoOne->momentum.Pt()<60)h_egInvarMass_Pt50to60->Fill(InvMass);
		  else if(PhoOne->nPixelSeeds==0 && PhoOne->momentum.Pt()>60 && PhoOne->momentum.Pt()<80)h_egInvarMass_Pt60to80->Fill(InvMass);
		  else if(PhoOne->nPixelSeeds==0 && PhoOne->momentum.Pt()>80)h_egInvarMass_Pt80->Fill(InvMass);
		  h_rho_eg->Fill(Rho);h_NVertex_eg->Fill(NVertex);
		  h_EcalIsoDR03RhoCorr_eg->Fill(PhoOne->ecalRecHitSumEtConeDR03-PUCorr_ECAL*Rho25);h_EcalIsoDR03RhoCorr_eg->Fill(PhoTwo->ecalRecHitSumEtConeDR03-PUCorr_ECAL*Rho25);
		  h_HcalIsoDR03RhoCorr_eg->Fill(PhoOne->hcalTowerSumEtConeDR03()-PUCorr_HCAL*Rho25);h_HcalIsoDR03RhoCorr_eg->Fill(PhoTwo->hcalTowerSumEtConeDR03()-PUCorr_HCAL*Rho25);
		  h_TrackIsoDR03RhoCorr_eg->Fill(PhoOne->trkSumPtHollowConeDR03-PUCorr_TRACK*Rho25);h_TrackIsoDR03RhoCorr_eg->Fill(PhoTwo->trkSumPtHollowConeDR03-PUCorr_TRACK*Rho25);
		  h_CombIsoDR03RhoCorr_eg->Fill(PhoOne->ecalRecHitSumEtConeDR03-PUCorr_ECAL*Rho25+PhoOne->hcalTowerSumEtConeDR03()-PUCorr_HCAL*Rho25+PhoOne->trkSumPtHollowConeDR03-PUCorr_TRACK*Rho25);
		  h_CombIsoDR03RhoCorr_eg->Fill(PhoTwo->ecalRecHitSumEtConeDR03-PUCorr_ECAL*Rho25+PhoTwo->hcalTowerSumEtConeDR03()-PUCorr_HCAL*Rho25+PhoTwo->trkSumPtHollowConeDR03-PUCorr_TRACK*Rho25);
		  h_DiEMPtVsMet_eg->Fill(met->met(),diEMpt);
		  if(HasDiJetPt){
		    if(JetCounter==0){h_egDiJetPt_0Jet->Fill(diJETpt);h_egDiEMPt_0Jet->Fill(diEMpt);}
		    else if(JetCounter==1){h_egDiJetPt_1Jet->Fill(diJETpt);h_egDiEMPt_1Jet->Fill(diEMpt);}
		    else if(JetCounter>=2){h_egDiJetPt_2Jet->Fill(diJETpt);h_egDiEMPt_2Jet->Fill(diEMpt);}
		    if(JetCounter>2){h_egDiJetPt_MoreThan2Jets->Fill(diJETpt);}
		  }
		  else{
		    if(JetCounter==0){h_egDiJetPt_0Jet->Fill(diEMpt);h_egDiEMPt_0Jet->Fill(diEMpt);}
		    else if(JetCounter==1){h_egDiJetPt_1Jet->Fill(diEMpt);h_egDiEMPt_1Jet->Fill(diEMpt);}
		    else if(JetCounter>=2){h_egDiJetPt_2Jet->Fill(diEMpt);h_egDiEMPt_2Jet->Fill(diEMpt);}
		    if(JetCounter>2){h_egDiJetPt_MoreThan2Jets->Fill(diEMpt);}
		    h_egMet_NoJetMatch->Fill(met->met());
		  }		  
		  h_egMR->Fill(RazrMr);
		  h_egR2->Fill(RazrR2);
		  h_egR2vsMR->Fill(RazrMr,RazrR2);
		  h_egdPhi->Fill(dPhi);
		  h_egAlphaT->Fill(alphaT);
		  if(JetCounter==0)h_egAlphaT_0Jet->Fill(alphaT);
		  h_egPhotonLessHt->Fill(PhotonLessHt);
		  h_egPhotonLessHtVsMET->Fill(met->met(),PhotonLessHt);
		  h_egMETdPhiLead->Fill(getDphi(PhoOne->caloPosition.Phi(),met->mEt.Phi()));
		  h_egMETdPhiTrail->Fill(getDphi(PhoTwo->caloPosition.Phi(),met->mEt.Phi()));
		  //}//seedTime
		  if(doJetReq){
		    if(printLevel>1)cout<<"Inside eg_JetReq"<<endl;
		    neg_JetReq++;iseg_JetReq=true;
		    egPTs_JetReq=make_pair(PhoOne->momentum.Pt(),PhoTwo->momentum.Pt());
		    egJetReqRunEvent.push_back(make_pair(event->runNumber,event->eventNumber));
		    h_egMet_JetReq->Fill(met->met());
		    h_SeedTimeVsEta_eg_JetReq->Fill(PhoOne->caloPosition.Eta(),PhoOne->seedTime);
		    h_SeedTimeVsEta_eg_JetReq->Fill(PhoTwo->caloPosition.Eta(),PhoTwo->seedTime);
		    h_egDiEMPt_JetReq->Fill(diEMpt);
		    h_egPixelSeeds_JetReq->Fill(PhoOne->nPixelSeeds);h_egPixelSeeds_JetReq->Fill(PhoTwo->nPixelSeeds);
		    h_egInvarMass_JetReq->Fill(InvMass);
		    if(HasDiJetPt){
		      if(JetCounter==0){h_egDiJetPt_JetReq_0Jet->Fill(diJETpt);h_egDiEMPt_JetReq_0Jet->Fill(diEMpt);}
		      else if(JetCounter==1){h_egDiJetPt_JetReq_1Jet->Fill(diJETpt);h_egDiEMPt_JetReq_1Jet->Fill(diEMpt);}
		      else if(JetCounter>=2){h_egDiJetPt_JetReq_2Jet->Fill(diJETpt);h_egDiEMPt_JetReq_2Jet->Fill(diEMpt);}
		      if(JetCounter>2){h_egDiJetPt_JetReq_MoreThan2Jets->Fill(diJETpt);}
		    }
		    else{
		      if(JetCounter==0){h_egDiJetPt_JetReq_0Jet->Fill(diEMpt);h_egDiEMPt_JetReq_0Jet->Fill(diEMpt);}
		      else if(JetCounter==1){h_egDiJetPt_JetReq_1Jet->Fill(diEMpt);h_egDiEMPt_JetReq_1Jet->Fill(diEMpt);}
		      else if(JetCounter>=2){h_egDiJetPt_JetReq_2Jet->Fill(diEMpt);h_egDiEMPt_JetReq_2Jet->Fill(diEMpt);}
		      if(JetCounter>2){h_egDiJetPt_MoreThan2Jets->Fill(diEMpt);}
		      h_egMet_NoJetMatch_JetReq->Fill(met->met());
		    }
		
		  }//end eg jet req case
		}//end eg (nPixelSeeds)
	      }//end if(InvMass>40)
	    }//end if(!PhosTooCloseDR && !PhosTooCloseDphi)
	  }//end if(PhoOne->momentum.Et()>40.) - end gg,ee,eg no jet req case
	}//if pho_Cands>=2 && ...
 
	//now do fakes
	if(fake_Cands.size()>=2 && (*fake_Cands.begin())->momentum.Et()>40.){  //can't make this an else if, what to do if there are gg and ff?
	  if(pho_Cands.size()>=2){TwoPhosAndTwoFakes=true;nTwoPhosAndTwoFakes++;}
	  //if(checkDouble){cout<<"PhoOne and PhoTwo getting reassigned!!!!"<<endl;}
	  ndiPhoCands++;
	  //Loop through pho_Cands to assign PhoOne and PhoTwo
	  //This mades sure the two pho objects are separated by dR>0.8 && dPhi>0.05 for no jet req case
	  FakesTooCloseDR=true;FakesTooCloseDphi=true;breakFake=false;
	  for(std::vector<susy::Photon*>::iterator fake_it = fake_Cands.begin();fake_it<(fake_Cands.end())-1;fake_it++){
	    if(breakFake)break;
	    PhoOne=*fake_it;
	    for(std::vector<susy::Photon*>::iterator fake_it2 = ++fake_it;fake_it2!=fake_Cands.end();fake_it2++){
	      if(   !isSameObject(PhoOne->caloPosition,(*fake_it2)->caloPosition,0.6) 
		    && !tooClosePhi(PhoOne->caloPosition,(*fake_it2)->caloPosition,0.05)    ){
		PhoTwo=*fake_it2;FakesTooCloseDR=false;FakesTooCloseDphi=false;breakFake=true;
		//This checks if there is at least one loose jet separated by 0.5 from both phos
		doJetReqFake=false;JetIsolatedFromBothPhos=false;
		if(!FakesTooCloseDR && !FakesTooCloseDphi && AtLeastOnePFJet){
		  //use pfJets.first here because we want all jets
		  for(std::vector<susy::PFJet*>::iterator jet_it = pfJets.first.begin(); jet_it != pfJets.first.end(); jet_it++){	
		    IsGoodJet=true;
		    //first clean from eles and muons
		    for(std::vector<susy::Electron*>::iterator ele_it = pfEles.begin();ele_it!=pfEles.end();ele_it++){
		      if(isSameObject((*ele_it)->momentum,(*jet_it)->momentum,0.5)){IsGoodJet=false;break;}
		    }
		    for(std::vector<susy::Muon*>::iterator mu_it = Muons.begin();mu_it!=Muons.end();mu_it++){
		      if(isSameObject((*mu_it)->momentum,(*jet_it)->momentum,0.5)){IsGoodJet=false;break;}
		    }
		    //now clean from the 2 pho objects
		    if(IsGoodJet){
		      if(!isSameObject(PhoOne->caloPosition,(*jet_it)->momentum,0.5) && !isSameObject(PhoTwo->caloPosition,(*jet_it)->momentum,0.5)){
			JetIsolatedFromBothPhos=true;ndiPhoCands_JetReq++;break;
		      }
		      else{JetIsolatedFromBothPhos=false;}
		    }
		  }
		}
		if(JetIsolatedFromBothPhos==true)doJetReqFake=true;
		break;
	      }
	      else{
		if(isSameObject(PhoOne->caloPosition,(*fake_it2)->caloPosition,0.6)){
		  FakesTooCloseDR=true;nFakesFailDR++;
		  if(printLevel>0)cout <<"FakesFailDR!"<< "  Run: "<<event->runNumber<<"  Event: "<<event->eventNumber<<"\n";
		}
		if(tooClosePhi(PhoOne->caloPosition,(*fake_it2)->caloPosition,0.05)){
		  FakesTooCloseDphi=true;nFakesFailDphi++;
		  if(printLevel>0)cout <<"FakesFailDphi!"<< "  Run: "<<event->runNumber<<"  Event: "<<event->eventNumber<<"\n";
		}	     
	      }
	    }
	  }
	  if(printLevel>1)cout<<"doJetReqFake: "<<doJetReqFake<<endl;

	  if(PhoOne->momentum.Et()>40. || PhoTwo->momentum.Et()>40.){
	    if(!FakesTooCloseDR && !FakesTooCloseDphi){
	      RazrMr = GetRazrMr(PhoOne,PhoTwo);
	      RazrR2 = GetRazrR2(PhoOne,PhoTwo,met);
	      float InvMass=InvariantMass(PhoOne->momentum,PhoTwo->momentum);
	      float diEMpt=GetDiEmPt(PhoOne,PhoTwo);
	      float diJETpt=0.;
	      //try invariant mass cut to get rid of crap in tail
	      if(InvMass>-9999.){
		//find out how many isolated jets
		JetCounter=0;
		//use pfJets.first here because we want all jets
		for(std::vector<susy::PFJet*>::iterator jet_it = pfJets.first.begin(); jet_it != pfJets.first.end(); jet_it++){	
		  //to be counted as a jet, it must be dr>0.5 separated from the 2 objects, plus pfEles and pfMuons	
		  IsGoodJet=true;
		  /*
		    for(std::vector<susy::Photon*>::iterator eg_it = pho_Cands.begin();eg_it!=pho_Cands.end();eg_it++){
		    if(isSameObject((*eg_it)->caloPosition,(*jet_it)->momentum,0.5)){IsGoodJet=false;}
		    }
		    for(std::vector<susy::Photon*>::iterator f_it = fake_Cands.begin();f_it!=fake_Cands.end();f_it++){
		    if(isSameObject((*f_it)->caloPosition,(*jet_it)->momentum,0.5)){IsGoodJet=false;}
		    }
		  */
		  bool check=true;
		  if(isSameObject(PhoOne->caloPosition,(*jet_it)->momentum,0.5)){IsGoodJet=false;}		  
		  if(isSameObject(PhoTwo->caloPosition,(*jet_it)->momentum,0.5)){IsGoodJet=false;}
		  check=IsGoodJet;
		  for(std::vector<susy::Electron*>::iterator ele_it = pfEles.begin();ele_it!=pfEles.end();ele_it++){
		    if(isSameObject((*ele_it)->momentum,(*jet_it)->momentum,0.5)){IsGoodJet=false;}
		  }
		  if(check && !IsGoodJet) nCnt[11]++; //number of jets that fail ele cleaning from e/g
		  check=true;
		  for(std::vector<susy::Muon*>::iterator mu_it = Muons.begin();mu_it!=Muons.end();mu_it++){
		    if(isSameObject((*mu_it)->momentum,(*jet_it)->momentum,0.5)){IsGoodJet=false;}
		  }
		  if(check && !IsGoodJet) nCnt[12]++; //number of jets that fail muon cleaning from e/g
		  if(IsGoodJet)JetCounter++;
		}
		/*  //old way of counting jets - only require jet to be isolated from PhoOne and PhoTwo
		    for(std::vector<susy::PFJet*>::iterator jet_it = pfJets.begin(); jet_it != pfJets.end(); jet_it++){	
		    if(!isSameObject(PhoOne->caloPosition,(*jet_it)->momentum,0.5) && !isSameObject(PhoTwo->caloPosition,(*jet_it)->momentum,0.5)){
		    JetCounter++;
		    }
		    }*/
		//require seed time within +-3ns for no jet requirement case
		//if(fabs(PhoOne->seedTime)<3. && fabs(PhoTwo->seedTime)<3.){
		if(!(event->runNumber==194115 && event->eventNumber==385483032)){//overlap with gg/ee/eg
		  if(isgg || isee || iseg){
		    if(isgg)cout<<"gg and ff event! ggPTs:"<<ggPTs.first<<" , "<<ggPTs.second<<endl;
		    if(isee)cout<<"ee and ff event! eePTs:"<<eePTs.first<<" , "<<eePTs.second<<endl;
		    if(iseg)cout<<"eg and ff event! egPTs:"<<egPTs.first<<" , "<<egPTs.second<<endl;
		    cout<<"                 ffPTs:"<<PhoOne->momentum.Pt()<<" , "<<PhoTwo->momentum.Pt()<<"   Run:"<<event->runNumber<<"  Event:"<<event->eventNumber<<endl;
		  }
		  if(printLevel>1)cout<<"Inside ff"<<endl;
		  nff++;isff=true;
		  ffRunEvent.push_back(make_pair(event->runNumber,event->eventNumber));
		  h_SeedTimeVsEta_ff->Fill(PhoOne->caloPosition.Eta(),PhoOne->seedTime);
		  h_SeedTimeVsEta_ff->Fill(PhoTwo->caloPosition.Eta(),PhoTwo->seedTime);
		  h_MetVsSeedTime_ff->Fill(PhoOne->seedTime, met->met());
		  h_MetVsSeedTime_ff->Fill(PhoTwo->seedTime, met->met());
		  h_NVertexVsMET_ff->Fill(met->met(),NVertex);
		  h_ffPerInstLumi->Fill(event->avgInsRecLumi);
		  if(isJetFail){bool go2=true;
		    for(std::vector<susy::PFJet*>::iterator jet_it = pfJetsFail.begin(); jet_it != pfJetsFail.end(); jet_it++){	
		      if( isSameObject(PhoOne->caloPosition,(*jet_it)->momentum,0.5) || isSameObject(PhoTwo->caloPosition,(*jet_it)->momentum,0.5) ) go2=false;
		    }
		    if(go2){nCnt[13]++;h_ffMet_JetFail->Fill(met->met());h_MET_JetFail->Fill(met->met());}
		  }
		  // Di-Jet Pt
		  //Try to match pho objects to jets
		  //use pfJets.second here because we want L1FastL2L3/L2L3 for e/g/f matching
		  MatchPhosToJets(PhoOne,PhoTwo,pfJets.second,DiJetPtJet1,DiJetPtJet2,HasDiJetPt,0.3);
		  if(HasDiJetPt){
		    diJETpt = GetDiJetPt(DiJetPtJet1,DiJetPtJet2);
		    h_ffDiJetPt->Fill(diJETpt);
		    h_ffDiJetPtOverDiEMPtVsDiEMPt->Fill(diEMpt,diJETpt/diEMpt);
		  }
		  else{h_ffDiJetPt->Fill(diEMpt);}
		  //cout<<"HasDiJetPt: " << HasDiJetPt<<"  DiEMPt: "<<diEMpt<<"  DiJETpt: "<<diJETpt<<endl;
		  h_r9_ff->Fill(PhoOne->r9);h_r9_ff->Fill(PhoTwo->r9);
		  h_ffPt->Fill(PhoOne->momentum.Pt());h_ffPt->Fill(PhoTwo->momentum.Pt());
		  h_sumEt_ff->Fill(met->sumEt);
		  h_ffMet->Fill(met->met(),PUweight);
		  for(int i=0;i<50;i++){
		    if((PhoOne->ecalRecHitSumEtConeDR03+PhoOne->hcalTowerSumEtConeDR03()+PhoOne->trkSumPtHollowConeDR03-PUCorr_ECAL*Rho25-PUCorr_HCAL*Rho25-PUCorr_TRACK*Rho25)<(float)(i+1) 
		       && (PhoTwo->ecalRecHitSumEtConeDR03+PhoTwo->hcalTowerSumEtConeDR03()+PhoTwo->trkSumPtHollowConeDR03-PUCorr_ECAL*Rho25-PUCorr_HCAL*Rho25-PUCorr_TRACK*Rho25)<(float)(i+1)){
		      h_ffMetChi2[i]->Fill(met->met());
		    }
		    if((PhoOne->chargedHadronIso+PhoOne->neutralHadronIso+PhoOne->photonIso-(PUCorr_chargedHadron+PUCorr_neutralHadron+PUCorr_photon)*Rho25)<(float)(i+1)
		       && (PhoTwo->chargedHadronIso+PhoTwo->neutralHadronIso+PhoTwo->photonIso-(PUCorr_chargedHadron+PUCorr_neutralHadron+PUCorr_photon)*Rho25)<(float)(i+1) ){
		      h_ffMetChi2PF[i]->Fill(met->met());
		    }
		  }
		  h_ffMet_reweight_binned->Fill(met->met(),GetMetReweight(diEMpt,"ff",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		  if(HasDiJetPt){
		    h_ffMet_reweightJet_binned->Fill(met->met(),GetMetReweight(diJETpt,"ff",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		    h_ffDiJetPt_reweight_binned->Fill(diJETpt,GetMetReweight(diJETpt,"ff",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		    if(NVertex<10)h_ffMet_NV0_10->Fill(met->met(),GetMetReweight(diJETpt,"ff",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		    else if(NVertex>=10 && NVertex<15)h_ffMet_NV10_15->Fill(met->met(),GetMetReweight(diJETpt,"ff",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		    else if(NVertex>=15)h_ffMet_NV15up->Fill(met->met(),GetMetReweight(diJETpt,"ff",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		  }
		  else{
		    h_ffMet_reweightJet_binned->Fill(met->met(),GetMetReweight(diEMpt,"ff",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		    h_ffDiJetPt_reweight_binned->Fill(diEMpt,GetMetReweight(diEMpt,"ff",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		    if(NVertex<10)h_ffMet_NV0_10->Fill(met->met(),GetMetReweight(diEMpt,"ff",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		    else if(NVertex>=10 && NVertex<15)h_ffMet_NV10_15->Fill(met->met(),GetMetReweight(diEMpt,"ff",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		    else if(NVertex>=15)h_ffMet_NV15up->Fill(met->met(),GetMetReweight(diEMpt,"ff",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		    h_ffMet_NoJetMatch->Fill(met->met());
		  }	
		  toys->cd();
		  if(HasDiJetPt){
		    for(int i=0;i<1000;i++){
		      //title="eeMet_reweight_binned_toy";title+=i;
		      ffMet_reweight_binned_toy[i]->Fill(met->met(),rando.Gaus(GetMetReweight(diJETpt,"ff",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first,GetMetReweight(diJETpt,"ff",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).second));
		      //ffMet_reweight_binned_toy[i]->Fill(met->met(),rando.Poisson(GetMetReweight(diEMpt,"ff",EEreweights,FFreweights).first));
		    }
		  }
		  else{
		    for(int i=0;i<1000;i++){
		      //title="eeMet_reweight_binned_toy";title+=i;
		      ffMet_reweight_binned_toy[i]->Fill(met->met(),rando.Gaus(GetMetReweight(diEMpt,"ff",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first,GetMetReweight(diEMpt,"ff",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).second));
		      //ffMet_reweight_binned_toy[i]->Fill(met->met(),rando.Poisson(GetMetReweight(diEMpt,"ff",EEreweights,FFreweights).first));
		    }
		  }
		  fout->cd();
		  h_ffSigIetaIeta->Fill(PhoOne->sigmaIetaIeta);h_ffSigIetaIeta->Fill(PhoTwo->sigmaIetaIeta);
		  h_ffDiEMPt->Fill(diEMpt);
		  if(InvMass<40.)h_ffDiEMPtInvMassBelow40->Fill(diEMpt);
		  else h_ffDiEMPtInvMassAbove40->Fill(diEMpt);
		  h_ffDiEMPt_reweight_binned->Fill(diEMpt,GetMetReweight(diEMpt,"ff",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		  h_ffInvarMass->Fill(InvMass);	
		  if(met->met()>30)h_ffInvarMassMET30->Fill(InvMass);
		  h_rho_ff->Fill(Rho);h_NVertex_ff->Fill(NVertex);
		  h_EcalIsoDR03RhoCorr_ff->Fill(PhoOne->ecalRecHitSumEtConeDR03-PUCorr_ECAL*Rho25);h_EcalIsoDR03RhoCorr_ff->Fill(PhoTwo->ecalRecHitSumEtConeDR03-PUCorr_ECAL*Rho25);
		  h_HcalIsoDR03RhoCorr_ff->Fill(PhoOne->hcalTowerSumEtConeDR03()-PUCorr_HCAL*Rho25);h_HcalIsoDR03RhoCorr_ff->Fill(PhoTwo->hcalTowerSumEtConeDR03()-PUCorr_HCAL*Rho25);
		  h_TrackIsoDR03RhoCorr_ff->Fill(PhoOne->trkSumPtHollowConeDR03-PUCorr_TRACK*Rho25);h_TrackIsoDR03RhoCorr_ff->Fill(PhoTwo->trkSumPtHollowConeDR03-PUCorr_TRACK*Rho25);
		  h_CombIsoDR03RhoCorr_ff->Fill(PhoOne->ecalRecHitSumEtConeDR03-PUCorr_ECAL*Rho25+PhoOne->hcalTowerSumEtConeDR03()-PUCorr_HCAL*Rho25+PhoOne->trkSumPtHollowConeDR03-PUCorr_TRACK*Rho25);
		  h_CombIsoDR03RhoCorr_ff->Fill(PhoTwo->ecalRecHitSumEtConeDR03-PUCorr_ECAL*Rho25+PhoTwo->hcalTowerSumEtConeDR03()-PUCorr_HCAL*Rho25+PhoTwo->trkSumPtHollowConeDR03-PUCorr_TRACK*Rho25);
		  h_DiEMPtVsMet_ff->Fill(met->met(),diEMpt);
		  if(HasDiJetPt){
		    if(JetCounter==0){h_ffDiJetPt_0Jet->Fill(diJETpt);h_ffDiEMPt_0Jet->Fill(diEMpt);}
		    else if(JetCounter==1){h_ffDiJetPt_1Jet->Fill(diJETpt);h_ffDiEMPt_1Jet->Fill(diEMpt);}
		    else if(JetCounter>=2){h_ffDiJetPt_2Jet->Fill(diJETpt);h_ffDiEMPt_2Jet->Fill(diEMpt);}
		    if(JetCounter>2){h_ffDiJetPt_MoreThan2Jets->Fill(diJETpt);}
		  }
		  else{
		    if(JetCounter==0){h_ffDiJetPt_0Jet->Fill(diEMpt);h_ffDiEMPt_0Jet->Fill(diEMpt);}
		    else if(JetCounter==1){h_ffDiJetPt_1Jet->Fill(diEMpt);h_ffDiEMPt_1Jet->Fill(diEMpt);}
		    else if(JetCounter>=2){h_ffDiJetPt_2Jet->Fill(diEMpt);h_ffDiEMPt_2Jet->Fill(diEMpt);}
		    if(JetCounter>2){h_ffDiJetPt_MoreThan2Jets->Fill(diEMpt);}
		  }
		  h_ffMR->Fill(RazrMr);
		  h_ffR2->Fill(RazrR2);
		  h_ffR2vsMR->Fill(RazrMr,RazrR2);
		  h_ffMR_reweight_binned->Fill(RazrMr,GetMetReweight(diEMpt,"ff",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,-1).first);
		  h_ffR2_reweight_binned->Fill(RazrR2,GetMetReweight(diEMpt,"ff",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,-1).first);
		  h_ffR2vsMR_reweight_binned->Fill(RazrMr,RazrR2,GetMetReweight(diEMpt,"ff",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,-1).first);
		  float dPhi = std::fabs(TVector2::Phi_mpi_pi(PhoOne->caloPosition.Phi() - PhoTwo->caloPosition.Phi()));
		  if(HasDiJetPt) alphaT = GetAlphaT(DiJetPtJet1->momentum, DiJetPtJet2->momentum);
		  else alphaT = GetAlphaT(PhoOne->momentum, PhoTwo->momentum);
		  if(HasDiJetPt) PhotonLessHt=GetPhotonLessHt(met->sumEt,DiJetPtJet1->momentum,DiJetPtJet2->momentum);
		  else PhotonLessHt=GetPhotonLessHt(met->sumEt,PhoOne->momentum,PhoTwo->momentum);
		  h_ffdPhi->Fill(dPhi);
		  h_ffAlphaT->Fill(alphaT);
		  if(JetCounter==0)h_ffAlphaT_0Jet->Fill(alphaT);
		  h_ffPhotonLessHt->Fill(PhotonLessHt);
		  h_ffPhotonLessHtVsMET->Fill(met->met(),PhotonLessHt);
		  h_ffMETdPhiLead->Fill(getDphi(PhoOne->caloPosition.Phi(),met->mEt.Phi()));
		  h_ffMETdPhiTrail->Fill(getDphi(PhoTwo->caloPosition.Phi(),met->mEt.Phi()));
		  //}//seedTime	
		}//end gg/ee/eg and ff events
		//}//end if(PhoOne->momentum.Et()>40.) - end ff no jet req case
		if(doJetReqFake){
		  if(printLevel>1)cout<<"Inside ff_JetReq"<<endl;
		  nff_JetReq++;isff_JetReq=true;
		  ffJetReqRunEvent.push_back(make_pair(event->runNumber,event->eventNumber));
		  h_ffMet_JetReq->Fill(met->met(),PUweight);
		  for(int i=0;i<50;i++){
		    if((PhoOne->ecalRecHitSumEtConeDR03+PhoOne->hcalTowerSumEtConeDR03()+PhoOne->trkSumPtHollowConeDR03-PUCorr_ECAL*Rho25-PUCorr_HCAL*Rho25-PUCorr_TRACK*Rho25)<(float)(i+1) 
		       && (PhoTwo->ecalRecHitSumEtConeDR03+PhoTwo->hcalTowerSumEtConeDR03()+PhoTwo->trkSumPtHollowConeDR03-PUCorr_ECAL*Rho25-PUCorr_HCAL*Rho25-PUCorr_TRACK*Rho25)<(float)(i+1) ){
		      h_ffMetChi2_JetReq[i]->Fill(met->met());
		    }
		    if((PhoOne->chargedHadronIso+PhoOne->neutralHadronIso+PhoOne->photonIso-(PUCorr_chargedHadron+PUCorr_neutralHadron+PUCorr_photon)*Rho25)<(float)(i+1)
		       && (PhoTwo->chargedHadronIso+PhoTwo->neutralHadronIso+PhoTwo->photonIso-(PUCorr_chargedHadron+PUCorr_neutralHadron+PUCorr_photon)*Rho25)<(float)(i+1) ){
		      h_ffMetChi2PF_JetReq[i]->Fill(met->met());
		    }
		  }
		  h_SeedTimeVsEta_ff_JetReq->Fill(PhoOne->caloPosition.Eta(),PhoOne->seedTime);
		  h_SeedTimeVsEta_ff_JetReq->Fill(PhoTwo->caloPosition.Eta(),PhoTwo->seedTime);
		  h_ffMet_reweight_binned_JetReq->Fill(met->met(),GetMetReweight(diEMpt,"ff",EEreweights_JetReq,FFreweights_JetReq,EESideBandLowReweightsJet_JetReq,EESideBandHighReweightsJet_JetReq,-1).first);
		  h_ffDiEMPt_reweight_binned_JetReq->Fill(diEMpt,GetMetReweight(diEMpt,"ff",EEreweights_JetReq,FFreweights_JetReq,EESideBandLowReweightsJet_JetReq,EESideBandHighReweightsJet_JetReq,JetCounter).first);
		  h_ffDiEMPt_JetReq->Fill(diEMpt);
		  //h_ffPixelSeeds_JetReq->Fill(PhoOne->nPixelSeeds);h_ffPixelSeeds_JetReq->Fill(PhoTwo->nPixelSeeds);
		  h_ffInvarMass_JetReq->Fill(InvMass);
		  // Di-Jet Pt
		  //Try to match pho objects to jets
		  //use pfJets.second here because we want L1FastL2L3/L2L3 for e/g/f matching
		  MatchPhosToJets(PhoOne,PhoTwo,pfJets.second,DiJetPtJet1,DiJetPtJet2,HasDiJetPt,0.3);
		  if(HasDiJetPt){
		    diJETpt = GetDiJetPt(DiJetPtJet1,DiJetPtJet2);
		    h_ffMet_reweightJet_binned_JetReq->Fill(met->met(),GetMetReweight(diJETpt,"ff",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		    h_ffDiJetPt_reweight_binned_JetReq->Fill(diJETpt,GetMetReweight(diJETpt,"ff",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		    h_ffDiJetPt_JetReq->Fill(diJETpt);
		    h_ffDiJetPtOverDiEMPtVsDiEMPt_JetReq->Fill(diEMpt,diJETpt/diEMpt);
		    if(JetCounter==0){h_ffDiJetPt_JetReq_0Jet->Fill(diJETpt);h_ffDiEMPt_JetReq_0Jet->Fill(diEMpt);}
		    else if(JetCounter==1){h_ffDiJetPt_JetReq_1Jet->Fill(diJETpt);h_ffDiEMPt_JetReq_1Jet->Fill(diEMpt);}
		    else if(JetCounter>=2){h_ffDiJetPt_JetReq_2Jet->Fill(diJETpt);h_ffDiEMPt_JetReq_2Jet->Fill(diEMpt);}
		    if(JetCounter>2){h_ffDiJetPt_JetReq_MoreThan2Jets->Fill(diJETpt);}
		  }
		  else{
		    h_ffDiJetPt_JetReq->Fill(diEMpt);
		    h_ffMet_reweightJet_binned_JetReq->Fill(met->met(),GetMetReweight(diEMpt,"ff",EEreweightsJet,FFreweightsJet,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first);
		    h_ffDiJetPt_reweight_binned_JetReq->Fill(diEMpt,GetMetReweight(diEMpt,"ff",EEreweights_JetReq,FFreweights_JetReq,EESideBandLowReweightsJet_JetReq,EESideBandHighReweightsJet_JetReq,JetCounter).first);
		    if(JetCounter==0){h_ffDiJetPt_JetReq_0Jet->Fill(diEMpt);h_ffDiEMPt_JetReq_0Jet->Fill(diEMpt);}
		    else if(JetCounter==1){h_ffDiJetPt_JetReq_1Jet->Fill(diEMpt);h_ffDiEMPt_JetReq_1Jet->Fill(diEMpt);}
		    else if(JetCounter>=2){h_ffDiJetPt_JetReq_2Jet->Fill(diEMpt);h_ffDiEMPt_JetReq_2Jet->Fill(diEMpt);}
		    if(JetCounter>2){h_ffDiJetPt_JetReq_MoreThan2Jets->Fill(diEMpt);}
		    h_ffMet_NoJetMatch_JetReq->Fill(met->met());
		  }	
		  toys->cd();
		  if(HasDiJetPt){
		    for(int i=0;i<1000;i++){
		      //title="eeMet_reweight_binned_toy";title+=i;
		      ffMet_reweight_binned_toy_JetReq[i]->Fill(met->met(),rando.Gaus(GetMetReweight(diJETpt,"ff",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).first,GetMetReweight(diJETpt,"ff",EEreweights,FFreweights,EESideBandLowReweightsJet,EESideBandHighReweightsJet,JetCounter).second));
		      //ffMet_reweight_binned_toy[i]->Fill(met->met(),rando.Poisson(GetMetReweight(diEMpt,"ff",EEreweights,FFreweights).first));
		    }
		  }
		  else{
		    for(int i=0;i<1000;i++){
		      //title="eeMet_reweight_binned_toy";title+=i;
		      ffMet_reweight_binned_toy_JetReq[i]->Fill(met->met(),rando.Gaus(GetMetReweight(diEMpt,"ff",EEreweights_JetReq,FFreweights_JetReq,EESideBandLowReweightsJet_JetReq,EESideBandHighReweightsJet_JetReq,JetCounter).first,GetMetReweight(diEMpt,"ff",EEreweights_JetReq,FFreweights_JetReq,EESideBandLowReweightsJet_JetReq,EESideBandHighReweightsJet_JetReq,JetCounter).second));
		      //ffMet_reweight_binned_toy[i]->Fill(met->met(),rando.Poisson(GetMetReweight(diEMpt,"ff",EEreweights,FFreweights).first));
		    }
		  }
		  fout->cd();
		}
	      }//end if(InvMass>40)
	    }//end if(!FakesTooCloseDR && !FakesTooCloseDphi)
	  }//if fake_Cands>=2 && ... - end ff
	}	
	if(isgg || isee || iseg || isff) ndiPho++;
	if(isgg_JetReq || isee_JetReq || iseg_JetReq || isff_JetReq) ndiPho_JetReq++;
      
      	
	//Filter events that pass json and hlt, have >=1 good vertex
	//----and have at least two photons with Et>30 (at least one with Et>43) that pass R9, h/e
	if(enableFilter) {
	  bool filterThis = (isgg || isee || iseg || isff || isgg_JetReq || isee_JetReq || iseg_JetReq || isff_JetReq);
	  //bool filterThis = ( (isgg || isgg_JetReq) && met->met()>100 );
	  if(filterThis) {
	    nFiltered++;
	    filterTree->Fill();
	  }
	}// if(enableFilter)
    
	//Need to delete everything we new'ed:
	/*   cout<<" pho_Cands size after creation before final init: "<<pho_Cands.size()<<endl;
	     cout<<"fake_Cands size after creation before final init: "<<fake_Cands.size()<<endl;*/
    
	PhoOne=NULL;
	//delete PhoOne;
	PhoTwo=NULL;
	//delete PhoTwo;
     
	for(size_t i=0;i<pho_Cands.size();++i){
	  pho_Cands[i]=NULL;
	}
	for(size_t i=0;i<pho_Cands_NoEtCut.size();++i){
	  pho_Cands_NoEtCut[i]=NULL;
	}
	for(size_t i=0;i<fake_Cands.size();++i){
	  fake_Cands[i]=NULL;
	}
	pho_Cands.clear();pho_Cands_NoEtCut.clear();fake_Cands.clear();
  	
      }//if phoMap!=end (after phoMap def)

      if(printLevel > 1) {
	cout << "------------------------------------------" << endl;
	cout << "              event summary" << endl;
	cout << "------------------------------------------" << endl;
	cout << "pfJets            : " << pfJets.first.size() << endl;
	cout << "------------------------------------------" << endl;
	cout << "------------------------------------------" << endl;
	cout << "met               : " << met->met() << endl;
	if(useTrigger){
	  for(size_t trig_it=0;trig_it<hltNames.size();trig_it++){
	    cout<<hltNames[trig_it]<<endl;
	  }
	}
      }


      if(printLevel > 1) cout << "Apply event level cuts from now on..." << endl;
      // filter conditions

      if(met->met() < 20.0) continue;
      nCnt[14]++;

    }//end run number selection
  }// for jentry
 
  //make DiEMPt Ratio
  TH1F *h_ggeeDiEMPtRatio = new TH1F("ggeeDiEMPtRatio","ggDiEMPt/eeDiEMPt - No Jet Requirement;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeDiEMPtRatio->GetSumw2N()==0)h_ggeeDiEMPtRatio->Sumw2();
  TH1F *h_ggeeDiEMPtRatio_0Jet = new TH1F("ggeeDiEMPtRatio_0Jet","ggDiEMPt/eeDiEMPt - Require exactly 0 Jets;DiEMPtP_{T};DiEMPtP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeDiEMPtRatio_0Jet->GetSumw2N()==0)h_ggeeDiEMPtRatio_0Jet->Sumw2();
  TH1F *h_ggeeDiEMPtRatio_1Jet = new TH1F("ggeeDiEMPtRatio_1Jet","ggDiEMPt/eeDiEMPt - Require exactly 1 Jet;DiEMPtP_{T};DiEMPtP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeDiEMPtRatio_1Jet->GetSumw2N()==0)h_ggeeDiEMPtRatio_1Jet->Sumw2();
  TH1F *h_ggeeDiEMPtRatio_2Jet = new TH1F("ggeeDiEMPtRatio_2Jet","ggDiEMPt/eeDiEMPt - Require at least 2 Jets;DiEMPtP_{T};DiEMPtP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeDiEMPtRatio_2Jet->GetSumw2N()==0)h_ggeeDiEMPtRatio_2Jet->Sumw2();
  TH1F *h_ggeeDiJetPtRatio = new TH1F("ggeeDiJetPtRatio","ggDiJetPt/eeDiJetPt - No Jet Requirement;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeDiJetPtRatio->GetSumw2N()==0)h_ggeeDiJetPtRatio->Sumw2();
  TH1F *h_ggeeDiJetPtRatio_0Jet = new TH1F("ggeeDiJetPtRatio_0Jet","ggDiJetPt/eeDiJetPt - Require exactly 0 Jets;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeDiJetPtRatio_0Jet->GetSumw2N()==0)h_ggeeDiJetPtRatio_0Jet->Sumw2();
  TH1F *h_ggeeDiJetPtRatio_1Jet = new TH1F("ggeeDiJetPtRatio_1Jet","ggDiJetPt/eeDiJetPt - Require exactly 1 Jet;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeDiJetPtRatio_1Jet->GetSumw2N()==0)h_ggeeDiJetPtRatio_1Jet->Sumw2();
  TH1F *h_ggeeDiJetPtRatio_2Jet = new TH1F("ggeeDiJetPtRatio_2Jet","ggDiJetPt/eeDiJetPt - Require at least 2 Jets;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeDiJetPtRatio_2Jet->GetSumw2N()==0)h_ggeeDiJetPtRatio_2Jet->Sumw2();
  TH1F *h_ggeeSideBandDiEMPtRatio = new TH1F("ggeeSideBandDiEMPtRatio","ggDiEMPt/eeSideBandDiEMPt - No Jet Requirement;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandDiEMPtRatio->GetSumw2N()==0)h_ggeeSideBandDiEMPtRatio->Sumw2();
  TH1F *h_ggeeSideBandDiJetPtRatio = new TH1F("ggeeSideBandDiJetPtRatio","ggDiJetPt/eeSideBandDiJetPt - No Jet Requirement;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandDiJetPtRatio->GetSumw2N()==0)h_ggeeSideBandDiJetPtRatio->Sumw2();
  TH1F *h_ggeeSideBandLowDiJetPtRatio = new TH1F("ggeeSideBandLowDiJetPtRatio","ggDiJetPt/eeSideBandLowDiJetPt - No Jet Requirement;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandLowDiJetPtRatio->GetSumw2N()==0)h_ggeeSideBandLowDiJetPtRatio->Sumw2();
  TH1F *h_ggeeSideBandLowDiJetPtRatio_0Jet = new TH1F("ggeeSideBandLowDiJetPtRatio_0Jet","ggDiJetPt/eeSideBandLowDiJetPt - Require exactly 0 Jets;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandLowDiJetPtRatio_0Jet->GetSumw2N()==0)h_ggeeSideBandLowDiJetPtRatio_0Jet->Sumw2();
  TH1F *h_ggeeSideBandLowDiJetPtRatio_1Jet = new TH1F("ggeeSideBandLowDiJetPtRatio_1Jet","ggDiJetPt/eeSideBandLowDiJetPt - Require exactly 1 Jet;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandLowDiJetPtRatio_1Jet->GetSumw2N()==0)h_ggeeSideBandLowDiJetPtRatio_1Jet->Sumw2();
  TH1F *h_ggeeSideBandLowDiJetPtRatio_2Jet = new TH1F("ggeeSideBandLowDiJetPtRatio_2Jet","ggDiJetPt/eeSideBandLowDiJetPt - Require at least 2 Jets;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandLowDiJetPtRatio_2Jet->GetSumw2N()==0)h_ggeeSideBandLowDiJetPtRatio_2Jet->Sumw2();
  TH1F *h_ggeeSideBandHighDiEMPtRatio = new TH1F("ggeeSideBandHighDiEMPtRatio","ggDiEMPt/eeSideBandHighDiEMPt - No Jet Requirement;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandHighDiEMPtRatio->GetSumw2N()==0)h_ggeeSideBandHighDiEMPtRatio->Sumw2();
  TH1F *h_ggeeSideBandHighDiEMPtRatio_0Jet = new TH1F("ggeeSideBandHighDiEMPtRatio_0Jet","ggDiEMPt/eeSideBandHighDiEMPt - Require exactly 0 Jets;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandHighDiEMPtRatio_0Jet->GetSumw2N()==0)h_ggeeSideBandHighDiEMPtRatio_0Jet->Sumw2();
  TH1F *h_ggeeSideBandHighDiEMPtRatio_1Jet = new TH1F("ggeeSideBandHighDiEMPtRatio_1Jet","ggDiEMPt/eeSideBandHighDiEMPt - Require exactly 1 Jet;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandHighDiEMPtRatio_1Jet->GetSumw2N()==0)h_ggeeSideBandHighDiEMPtRatio_1Jet->Sumw2();
  TH1F *h_ggeeSideBandHighDiEMPtRatio_2Jet = new TH1F("ggeeSideBandHighDiEMPtRatio_2Jet","ggDiEMPt/eeSideBandHighDiEMPt - Require at least 2 Jets;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandHighDiEMPtRatio_2Jet->GetSumw2N()==0)h_ggeeSideBandHighDiEMPtRatio_2Jet->Sumw2();
  TH1F *h_ggeeSideBandLowDiEMPtRatio = new TH1F("ggeeSideBandLowDiEMPtRatio","ggDiEMPt/eeSideBandLowDiEMPt - No Jet Requirement;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandLowDiEMPtRatio->GetSumw2N()==0)h_ggeeSideBandLowDiEMPtRatio->Sumw2();
  TH1F *h_ggeeSideBandLowDiEMPtRatio_0Jet = new TH1F("ggeeSideBandLowDiEMPtRatio_0Jet","ggDiEMPt/eeSideBandLowDiEMPt - Require exactly 0 Jets;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandLowDiEMPtRatio_0Jet->GetSumw2N()==0)h_ggeeSideBandLowDiEMPtRatio_0Jet->Sumw2();
  TH1F *h_ggeeSideBandLowDiEMPtRatio_1Jet = new TH1F("ggeeSideBandLowDiEMPtRatio_1Jet","ggDiEMPt/eeSideBandLowDiEMPt - Require exactly 1 Jet;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandLowDiEMPtRatio_1Jet->GetSumw2N()==0)h_ggeeSideBandLowDiEMPtRatio_1Jet->Sumw2();
  TH1F *h_ggeeSideBandLowDiEMPtRatio_2Jet = new TH1F("ggeeSideBandLowDiEMPtRatio_2Jet","ggDiEMPt/eeSideBandLowDiEMPt - Require at least 2 Jets;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandLowDiEMPtRatio_2Jet->GetSumw2N()==0)h_ggeeSideBandLowDiEMPtRatio_2Jet->Sumw2();
  TH1F *h_ggeeSideBandHighDiEMPtRatio_JetReq = new TH1F("ggeeSideBandHighDiEMPtRatio_JetReq","ggDiEMPt/eeSideBandHighDiEMPt - No Jet Requirement;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandHighDiEMPtRatio->GetSumw2N()==0)h_ggeeSideBandHighDiEMPtRatio_JetReq->Sumw2();
  TH1F *h_ggeeSideBandHighDiEMPtRatio_JetReq_0Jet = new TH1F("ggeeSideBandHighDiEMPtRatio_JetReq_0Jet","ggDiEMPt/eeSideBandHighDiEMPt - Require exactly 0 Jets;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandHighDiEMPtRatio_JetReq_0Jet->GetSumw2N()==0)h_ggeeSideBandHighDiEMPtRatio_JetReq_0Jet->Sumw2();
  TH1F *h_ggeeSideBandHighDiEMPtRatio_JetReq_1Jet = new TH1F("ggeeSideBandHighDiEMPtRatio_JetReq_1Jet","ggDiEMPt/eeSideBandHighDiEMPt - Require exactly 1 Jet;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandHighDiEMPtRatio_JetReq_1Jet->GetSumw2N()==0)h_ggeeSideBandHighDiEMPtRatio_JetReq_1Jet->Sumw2();
  TH1F *h_ggeeSideBandHighDiEMPtRatio_JetReq_2Jet = new TH1F("ggeeSideBandHighDiEMPtRatio_JetReq_2Jet","ggDiEMPt/eeSideBandHighDiEMPt - Require at least 2 Jets;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandHighDiEMPtRatio_JetReq_2Jet->GetSumw2N()==0)h_ggeeSideBandHighDiEMPtRatio_JetReq_2Jet->Sumw2();
  TH1F *h_ggeeSideBandLowDiEMPtRatio_JetReq = new TH1F("ggeeSideBandLowDiEMPtRatio_JetReq","ggDiEMPt/eeSideBandLowDiEMPt - No Jet Requirement;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandLowDiEMPtRatio->GetSumw2N()==0)h_ggeeSideBandLowDiEMPtRatio_JetReq->Sumw2();
  TH1F *h_ggeeSideBandLowDiEMPtRatio_JetReq_0Jet = new TH1F("ggeeSideBandLowDiEMPtRatio_JetReq_0Jet","ggDiEMPt/eeSideBandLowDiEMPt - Require exactly 0 Jets;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandLowDiEMPtRatio_JetReq_0Jet->GetSumw2N()==0)h_ggeeSideBandLowDiEMPtRatio_JetReq_0Jet->Sumw2();
  TH1F *h_ggeeSideBandLowDiEMPtRatio_JetReq_1Jet = new TH1F("ggeeSideBandLowDiEMPtRatio_JetReq_1Jet","ggDiEMPt/eeSideBandLowDiEMPt - Require exactly 1 Jet;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandLowDiEMPtRatio_JetReq_1Jet->GetSumw2N()==0)h_ggeeSideBandLowDiEMPtRatio_JetReq_1Jet->Sumw2();
  TH1F *h_ggeeSideBandLowDiEMPtRatio_JetReq_2Jet = new TH1F("ggeeSideBandLowDiEMPtRatio_JetReq_2Jet","ggDiEMPt/eeSideBandLowDiEMPt - Require at least 2 Jets;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandLowDiEMPtRatio_JetReq_2Jet->GetSumw2N()==0)h_ggeeSideBandLowDiEMPtRatio_JetReq_2Jet->Sumw2();
  TH1F *h_ggeeSideBandHighDiJetPtRatio = new TH1F("ggeeSideBandHighDiJetPtRatio","ggDiJetPt/eeSideBandHighDiJetPt - No Jet Requirement;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandHighDiJetPtRatio->GetSumw2N()==0)h_ggeeSideBandHighDiJetPtRatio->Sumw2();
  TH1F *h_ggeeSideBandHighDiJetPtRatio_0Jet = new TH1F("ggeeSideBandHighDiJetPtRatio_0Jet","ggDiJetPt/eeSideBandHighDiJetPt - Require exactly 0 Jets;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandHighDiJetPtRatio_0Jet->GetSumw2N()==0)h_ggeeSideBandHighDiJetPtRatio_0Jet->Sumw2();
  TH1F *h_ggeeSideBandHighDiJetPtRatio_1Jet = new TH1F("ggeeSideBandHighDiJetPtRatio_1Jet","ggDiJetPt/eeSideBandHighDiJetPt - Require exactly 1 Jet;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandHighDiJetPtRatio_1Jet->GetSumw2N()==0)h_ggeeSideBandHighDiJetPtRatio_1Jet->Sumw2();
  TH1F *h_ggeeSideBandHighDiJetPtRatio_2Jet = new TH1F("ggeeSideBandHighDiJetPtRatio_2Jet","ggDiJetPt/eeSideBandHighDiJetPt - Require at least 2 Jets;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandHighDiJetPtRatio_2Jet->GetSumw2N()==0)h_ggeeSideBandHighDiJetPtRatio_2Jet->Sumw2();
  TH1F *h_ggffDiEMPtRatio = new TH1F("ggffDiEMPtRatio","ggDiEMPt/ffDiEMPt - No Jet Requirement;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggffDiEMPtRatio->GetSumw2N()==0)h_ggffDiEMPtRatio->Sumw2();
  TH1F *h_ggffDiEMPtRatio_0Jet = new TH1F("ggffDiEMPtRatio_0Jet","ggDiEMPt/ffDiEMPt - Require exactly 0 Jets;DiEMPtP_{T};DiEMPtP_{T} Ratio",numBins,DiEMPtBins);if(h_ggffDiEMPtRatio_0Jet->GetSumw2N()==0)h_ggffDiEMPtRatio_0Jet->Sumw2();
  TH1F *h_ggffDiEMPtRatio_1Jet = new TH1F("ggffDiEMPtRatio_1Jet","ggDiEMPt/ffDiEMPt - Require exactly 1 Jet;DiEMPtP_{T};DiEMPtP_{T} Ratio",numBins,DiEMPtBins);if(h_ggffDiEMPtRatio_1Jet->GetSumw2N()==0)h_ggffDiEMPtRatio_1Jet->Sumw2();
  TH1F *h_ggffDiEMPtRatio_2Jet = new TH1F("ggffDiEMPtRatio_2Jet","ggDiEMPt/ffDiEMPt - Require at least 2 Jets;DiEMPtP_{T};DiEMPtP_{T} Ratio",numBins,DiEMPtBins);if(h_ggffDiEMPtRatio_2Jet->GetSumw2N()==0)h_ggffDiEMPtRatio_2Jet->Sumw2();
  TH1F *h_ggffDiJetPtRatio = new TH1F("ggffDiJetPtRatio","ggDiJetPt/ffDiJetPt - No Jet Requirement;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggffDiJetPtRatio->GetSumw2N()==0)h_ggffDiJetPtRatio->Sumw2();
  TH1F *h_ggffDiJetPtRatio_0Jet = new TH1F("ggffDiJetPtRatio_0Jet","ggDiJetPt/ffDiJetPt - Require exactly 0 Jets;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggffDiJetPtRatio_0Jet->GetSumw2N()==0)h_ggffDiJetPtRatio_0Jet->Sumw2();
  TH1F *h_ggffDiJetPtRatio_1Jet = new TH1F("ggffDiJetPtRatio_1Jet","ggDiJetPt/ffDiJetPt - Require exactly 1 Jet;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggffDiJetPtRatio_1Jet->GetSumw2N()==0)h_ggffDiJetPtRatio_1Jet->Sumw2();
  TH1F *h_ggffDiJetPtRatio_2Jet = new TH1F("ggffDiJetPtRatio_2Jet","ggDiJetPt/ffDiJetPt - Require at least 2 Jets;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggffDiJetPtRatio_2Jet->GetSumw2N()==0)h_ggffDiJetPtRatio_2Jet->Sumw2();

  TH1F *h_ggffDiEMPtRatio_JetReq = new TH1F("ggffDiEMPtRatio_JetReq","ggDiEMPt/ffDiEMPt - 1+ Jet Requirement;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggffDiEMPtRatio_JetReq->GetSumw2N()==0)h_ggffDiEMPtRatio_JetReq->Sumw2();
  TH1F *h_ggeeDiEMPtRatio_JetReq = new TH1F("ggeeDiEMPtRatio_JetReq","ggDiEMPt/eeDiEMPt - 1+ Jet Requirement;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeDiEMPtRatio_JetReq->GetSumw2N()==0)h_ggeeDiEMPtRatio_JetReq->Sumw2();
  TH1F *h_ggeeDiEMPtRatio_JetReq_0Jet = new TH1F("ggeeDiEMPtRatio_JetReq_0Jet","ggDiEMPt/eeDiEMPt - 1+ Jet Requirement Require exactly 0 Jets;DiEMPtP_{T};DiEMPtP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeDiEMPtRatio_JetReq_0Jet->GetSumw2N()==0)h_ggeeDiEMPtRatio_JetReq_0Jet->Sumw2();
  TH1F *h_ggeeDiEMPtRatio_JetReq_1Jet = new TH1F("ggeeDiEMPtRatio_JetReq_1Jet","ggDiEMPt/eeDiEMPt - 1+ Jet Requirement Require exactly 1 Jet;DiEMPtP_{T};DiEMPtP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeDiEMPtRatio_JetReq_1Jet->GetSumw2N()==0)h_ggeeDiEMPtRatio_JetReq_1Jet->Sumw2();
  TH1F *h_ggeeDiEMPtRatio_JetReq_2Jet = new TH1F("ggeeDiEMPtRatio_JetReq_2Jet","ggDiEMPt/eeDiEMPt - 1+ Jet Requirement Require at least 2 Jets;DiEMPtP_{T};DiEMPtP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeDiEMPtRatio_JetReq_2Jet->GetSumw2N()==0)h_ggeeDiEMPtRatio_JetReq_2Jet->Sumw2();
  TH1F *h_ggeeDiJetPtRatio_JetReq = new TH1F("ggeeDiJetPtRatio_JetReq","ggDiJetPt/eeDiJetPt - 1+ Jet Requirement;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeDiJetPtRatio_JetReq->GetSumw2N()==0)h_ggeeDiJetPtRatio_JetReq->Sumw2();
  TH1F *h_ggeeDiJetPtRatio_JetReq_0Jet = new TH1F("ggeeDiJetPtRatio_JetReq_0Jet","ggDiJetPt/eeDiJetPt - 1+ Jet Requirement Require exactly 0 Jets;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeDiJetPtRatio_JetReq_0Jet->GetSumw2N()==0)h_ggeeDiJetPtRatio_JetReq_0Jet->Sumw2();
  TH1F *h_ggeeDiJetPtRatio_JetReq_1Jet = new TH1F("ggeeDiJetPtRatio_JetReq_1Jet","ggDiJetPt/eeDiJetPt - 1+ Jet Requirement Require exactly 1 Jet;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeDiJetPtRatio_JetReq_1Jet->GetSumw2N()==0)h_ggeeDiJetPtRatio_JetReq_1Jet->Sumw2();
  TH1F *h_ggeeDiJetPtRatio_JetReq_2Jet = new TH1F("ggeeDiJetPtRatio_JetReq_2Jet","ggDiJetPt/eeDiJetPt - 1+ Jet Requirement Require at least 2 Jets;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeDiJetPtRatio_JetReq_2Jet->GetSumw2N()==0)h_ggeeDiJetPtRatio_JetReq_2Jet->Sumw2();
  TH1F *h_ggeeSideBandDiEMPtRatio_JetReq = new TH1F("ggeeSideBandDiEMPtRatio_JetReq","ggDiEMPt/eeSideBandDiEMPt - 1+ Jet Requirement;DiEMP_{T};DiEMP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandDiEMPtRatio_JetReq->GetSumw2N()==0)h_ggeeSideBandDiEMPtRatio_JetReq->Sumw2();
  TH1F *h_ggeeSideBandDiJetPtRatio_JetReq = new TH1F("ggeeSideBandDiJetPtRatio_JetReq","ggDiJetPt/eeSideBandDiJetPt - 1+ Jet Requirement;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandDiJetPtRatio_JetReq->GetSumw2N()==0)h_ggeeSideBandDiJetPtRatio_JetReq->Sumw2();
  TH1F *h_ggeeSideBandLowDiJetPtRatio_JetReq = new TH1F("ggeeSideBandLowDiJetPtRatio_JetReq","ggDiJetPt/eeSideBandLowDiJetPt - 1+ Jet Requirement;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandLowDiJetPtRatio_JetReq->GetSumw2N()==0)h_ggeeSideBandLowDiJetPtRatio_JetReq->Sumw2();
  TH1F *h_ggeeSideBandLowDiJetPtRatio_JetReq_0Jet = new TH1F("ggeeSideBandLowDiJetPtRatio_JetReq_0Jet","ggDiJetPt/eeSideBandLowDiJetPt - 1+ Jet Requirement Require exactly 0 Jets;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandLowDiJetPtRatio_JetReq_0Jet->GetSumw2N()==0)h_ggeeSideBandLowDiJetPtRatio_JetReq_0Jet->Sumw2();
  TH1F *h_ggeeSideBandLowDiJetPtRatio_JetReq_1Jet = new TH1F("ggeeSideBandLowDiJetPtRatio_JetReq_1Jet","ggDiJetPt/eeSideBandLowDiJetPt - 1+ Jet Requirement Require exactly 1 Jet;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandLowDiJetPtRatio_JetReq_1Jet->GetSumw2N()==0)h_ggeeSideBandLowDiJetPtRatio_JetReq_1Jet->Sumw2();
  TH1F *h_ggeeSideBandLowDiJetPtRatio_JetReq_2Jet = new TH1F("ggeeSideBandLowDiJetPtRatio_JetReq_2Jet","ggDiJetPt/eeSideBandLowDiJetPt - 1+ Jet Requirement Require at least 2 Jets;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandLowDiJetPtRatio_JetReq_2Jet->GetSumw2N()==0)h_ggeeSideBandLowDiJetPtRatio_JetReq_2Jet->Sumw2();
  TH1F *h_ggeeSideBandHighDiJetPtRatio_JetReq = new TH1F("ggeeSideBandHighDiJetPtRatio_JetReq","ggDiJetPt/eeSideBandHighDiJetPt - 1+ Jet Requirement;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandHighDiJetPtRatio_JetReq->GetSumw2N()==0)h_ggeeSideBandHighDiJetPtRatio_JetReq->Sumw2();
  TH1F *h_ggeeSideBandHighDiJetPtRatio_JetReq_0Jet = new TH1F("ggeeSideBandHighDiJetPtRatio_JetReq_0Jet","ggDiJetPt/eeSideBandHighDiJetPt - 1+ Jet Requirement Require exactly 0 Jets;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandHighDiJetPtRatio_JetReq_0Jet->GetSumw2N()==0)h_ggeeSideBandHighDiJetPtRatio_JetReq_0Jet->Sumw2();
  TH1F *h_ggeeSideBandHighDiJetPtRatio_JetReq_1Jet = new TH1F("ggeeSideBandHighDiJetPtRatio_JetReq_1Jet","ggDiJetPt/eeSideBandHighDiJetPt - 1+ Jet Requirement Require exactly 1 Jet;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandHighDiJetPtRatio_JetReq_1Jet->GetSumw2N()==0)h_ggeeSideBandHighDiJetPtRatio_JetReq_1Jet->Sumw2();
  TH1F *h_ggeeSideBandHighDiJetPtRatio_JetReq_2Jet = new TH1F("ggeeSideBandHighDiJetPtRatio_JetReq_2Jet","ggDiJetPt/eeSideBandHighDiJetPt - 1+ Jet Requirement Require at least 2 Jets;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggeeSideBandHighDiJetPtRatio_JetReq_2Jet->GetSumw2N()==0)h_ggeeSideBandHighDiJetPtRatio_JetReq_2Jet->Sumw2();
  TH1F *h_ggffDiEMPtRatio_JetReq_0Jet = new TH1F("ggffDiEMPtRatio_JetReq_0Jet","ggDiEMPt/ffDiEMPt - 1+ Jet Requirement Require exactly 0 Jets;DiEMPtP_{T};DiEMPtP_{T} Ratio",numBins,DiEMPtBins);if(h_ggffDiEMPtRatio_JetReq_0Jet->GetSumw2N()==0)h_ggffDiEMPtRatio_JetReq_0Jet->Sumw2();
  TH1F *h_ggffDiEMPtRatio_JetReq_1Jet = new TH1F("ggffDiEMPtRatio_JetReq_1Jet","ggDiEMPt/ffDiEMPt - 1+ Jet Requirement Require exactly 1 Jet;DiEMPtP_{T};DiEMPtP_{T} Ratio",numBins,DiEMPtBins);if(h_ggffDiEMPtRatio_JetReq_1Jet->GetSumw2N()==0)h_ggffDiEMPtRatio_JetReq_1Jet->Sumw2();
  TH1F *h_ggffDiEMPtRatio_JetReq_2Jet = new TH1F("ggffDiEMPtRatio_JetReq_2Jet","ggDiEMPt/ffDiEMPt - 1+ Jet Requirement Require at least 2 Jets;DiEMPtP_{T};DiEMPtP_{T} Ratio",numBins,DiEMPtBins);if(h_ggffDiEMPtRatio_JetReq_2Jet->GetSumw2N()==0)h_ggffDiEMPtRatio_JetReq_2Jet->Sumw2();
  TH1F *h_ggffDiJetPtRatio_JetReq = new TH1F("ggffDiJetPtRatio_JetReq","ggDiJetPt/ffDiJetPt - 1+ Jet Requirement;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggffDiJetPtRatio_JetReq->GetSumw2N()==0)h_ggffDiJetPtRatio_JetReq->Sumw2();
  TH1F *h_ggffDiJetPtRatio_JetReq_0Jet = new TH1F("ggffDiJetPtRatio_JetReq_0Jet","ggDiJetPt/ffDiJetPt - 1+ Jet Requirement Require exactly 0 Jets;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggffDiJetPtRatio_JetReq_0Jet->GetSumw2N()==0)h_ggffDiJetPtRatio_JetReq_0Jet->Sumw2();
  TH1F *h_ggffDiJetPtRatio_JetReq_1Jet = new TH1F("ggffDiJetPtRatio_JetReq_1Jet","ggDiJetPt/ffDiJetPt - 1+ Jet Requirement Require exactly 1 Jet;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggffDiJetPtRatio_JetReq_1Jet->GetSumw2N()==0)h_ggffDiJetPtRatio_JetReq_1Jet->Sumw2();
  TH1F *h_ggffDiJetPtRatio_JetReq_2Jet = new TH1F("ggffDiJetPtRatio_JetReq_2Jet","ggDiJetPt/ffDiJetPt - 1+ Jet Requirement Require at least 2 Jets;DiJetP_{T};DiJetP_{T} Ratio",numBins,DiEMPtBins);if(h_ggffDiJetPtRatio_JetReq_2Jet->GetSumw2N()==0)h_ggffDiJetPtRatio_JetReq_2Jet->Sumw2();


  //rebin diempt for ratio plots
  TH1F* h_ggDiEMPt_new=(TH1F*)h_ggDiEMPt->Rebin(numBins,"h_ggDiEMPt_new",DiEMPtBins);
  TH1F* h_ggDiJetPt_new=(TH1F*)h_ggDiJetPt->Rebin(numBins,"h_ggDiJetPt_new",DiEMPtBins);
  TH1F* h_ggDiEMPt_0Jet_new=(TH1F*)h_ggDiEMPt_0Jet->Rebin(numBins,"h_ggDiEMPt_0Jet_new",DiEMPtBins);
  TH1F* h_ggDiEMPt_1Jet_new=(TH1F*)h_ggDiEMPt_1Jet->Rebin(numBins,"h_ggDiEMPt_1Jet_new",DiEMPtBins);
  TH1F* h_ggDiEMPt_2Jet_new=(TH1F*)h_ggDiEMPt_2Jet->Rebin(numBins,"h_ggDiEMPt_2Jet_new",DiEMPtBins);
  TH1F* h_ggDiJetPt_0Jet_new=(TH1F*)h_ggDiJetPt_0Jet->Rebin(numBins,"h_ggDiJetPt_0Jet_new",DiEMPtBins);
  TH1F* h_ggDiJetPt_1Jet_new=(TH1F*)h_ggDiJetPt_1Jet->Rebin(numBins,"h_ggDiJetPt_1Jet_new",DiEMPtBins);
  TH1F* h_ggDiJetPt_2Jet_new=(TH1F*)h_ggDiJetPt_2Jet->Rebin(numBins,"h_ggDiJetPt_2Jet_new",DiEMPtBins);
  TH1F* h_eeDiEMPt_new=(TH1F*)h_eeDiEMPt->Rebin(numBins,"h_eeDiEMPt_new",DiEMPtBins);
  TH1F* h_eeDiJetPt_new=(TH1F*)h_eeDiJetPt->Rebin(numBins,"h_eeDiJetPt_new",DiEMPtBins);
  TH1F* h_eeDiEMPt_0Jet_new=(TH1F*)h_eeDiEMPt_0Jet->Rebin(numBins,"h_eeDiEMPt_0Jet_new",DiEMPtBins);
  TH1F* h_eeDiEMPt_1Jet_new=(TH1F*)h_eeDiEMPt_1Jet->Rebin(numBins,"h_eeDiEMPt_1Jet_new",DiEMPtBins);
  TH1F* h_eeDiEMPt_2Jet_new=(TH1F*)h_eeDiEMPt_2Jet->Rebin(numBins,"h_eeDiEMPt_2Jet_new",DiEMPtBins);
  TH1F* h_eeDiJetPt_0Jet_new=(TH1F*)h_eeDiJetPt_0Jet->Rebin(numBins,"h_eeDiJetPt_0Jet_new",DiEMPtBins);
  TH1F* h_eeDiJetPt_1Jet_new=(TH1F*)h_eeDiJetPt_1Jet->Rebin(numBins,"h_eeDiJetPt_1Jet_new",DiEMPtBins);
  TH1F* h_eeDiJetPt_2Jet_new=(TH1F*)h_eeDiJetPt_2Jet->Rebin(numBins,"h_eeDiJetPt_2Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandDiEMPt_new=(TH1F*)h_eeSideBandDiEMPt->Rebin(numBins,"h_eeSideBandDiEMPt_new",DiEMPtBins);
  TH1F* h_eeSideBandDiJetPt_new=(TH1F*)h_eeSideBandDiJetPt->Rebin(numBins,"h_eeSideBandDiJetPt_new",DiEMPtBins);
  TH1F* h_eeSideBandLowDiJetPt_new=(TH1F*)h_eeSideBandLowDiJetPt->Rebin(numBins,"h_eeSideBandLowDiJetPt_new",DiEMPtBins);
  TH1F* h_eeSideBandHighDiJetPt_new=(TH1F*)h_eeSideBandHighDiJetPt->Rebin(numBins,"h_eeSideBandHighDiJetPt_new",DiEMPtBins);
  TH1F* h_eeSideBandLowDiEMPt_new=(TH1F*)h_eeSideBandLowDiEMPt->Rebin(numBins,"h_eeSideBandLowDiEMPt_new",DiEMPtBins);
  TH1F* h_eeSideBandHighDiEMPt_new=(TH1F*)h_eeSideBandHighDiEMPt->Rebin(numBins,"h_eeSideBandHighDiEMPt_new",DiEMPtBins);
  TH1F* h_eeSideBandLowDiEMPt_0Jet_new=(TH1F*)h_eeSideBandLowDiEMPt_0Jet->Rebin(numBins,"h_eeSideBandLowDiEMPt_0Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandLowDiEMPt_1Jet_new=(TH1F*)h_eeSideBandLowDiEMPt_1Jet->Rebin(numBins,"h_eeSideBandLowDiEMPt_1Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandLowDiEMPt_2Jet_new=(TH1F*)h_eeSideBandLowDiEMPt_2Jet->Rebin(numBins,"h_eeSideBandLowDiEMPt_2Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandHighDiEMPt_0Jet_new=(TH1F*)h_eeSideBandHighDiEMPt_0Jet->Rebin(numBins,"h_eeSideBandHighDiEMPt_0Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandHighDiEMPt_1Jet_new=(TH1F*)h_eeSideBandHighDiEMPt_1Jet->Rebin(numBins,"h_eeSideBandHighDiEMPt_1Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandHighDiEMPt_2Jet_new=(TH1F*)h_eeSideBandHighDiEMPt_2Jet->Rebin(numBins,"h_eeSideBandHighDiEMPt_2Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandLowDiJetPt_0Jet_new=(TH1F*)h_eeSideBandLowDiJetPt_0Jet->Rebin(numBins,"h_eeSideBandLowDiJetPt_0Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandLowDiJetPt_1Jet_new=(TH1F*)h_eeSideBandLowDiJetPt_1Jet->Rebin(numBins,"h_eeSideBandLowDiJetPt_1Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandLowDiJetPt_2Jet_new=(TH1F*)h_eeSideBandLowDiJetPt_2Jet->Rebin(numBins,"h_eeSideBandLowDiJetPt_2Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandHighDiJetPt_0Jet_new=(TH1F*)h_eeSideBandHighDiJetPt_0Jet->Rebin(numBins,"h_eeSideBandHighDiJetPt_0Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandHighDiJetPt_1Jet_new=(TH1F*)h_eeSideBandHighDiJetPt_1Jet->Rebin(numBins,"h_eeSideBandHighDiJetPt_1Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandHighDiJetPt_2Jet_new=(TH1F*)h_eeSideBandHighDiJetPt_2Jet->Rebin(numBins,"h_eeSideBandHighDiJetPt_2Jet_new",DiEMPtBins);
  TH1F* h_ffDiEMPt_new=(TH1F*)h_ffDiEMPt->Rebin(numBins,"h_ffDiEMPt_new",DiEMPtBins);
  TH1F* h_ffDiJetPt_new=(TH1F*)h_ffDiJetPt->Rebin(numBins,"h_ffDiJetPt_new",DiEMPtBins);
  TH1F* h_ffDiEMPt_0Jet_new=(TH1F*)h_ffDiEMPt_0Jet->Rebin(numBins,"h_ffDiEMPt_0Jet_new",DiEMPtBins);
  TH1F* h_ffDiEMPt_1Jet_new=(TH1F*)h_ffDiEMPt_1Jet->Rebin(numBins,"h_ffDiEMPt_1Jet_new",DiEMPtBins);
  TH1F* h_ffDiEMPt_2Jet_new=(TH1F*)h_ffDiEMPt_2Jet->Rebin(numBins,"h_ffDiEMPt_2Jet_new",DiEMPtBins);
  TH1F* h_ffDiJetPt_0Jet_new=(TH1F*)h_ffDiJetPt_0Jet->Rebin(numBins,"h_ffDiJetPt_0Jet_new",DiEMPtBins);
  TH1F* h_ffDiJetPt_1Jet_new=(TH1F*)h_ffDiJetPt_1Jet->Rebin(numBins,"h_ffDiJetPt_1Jet_new",DiEMPtBins);
  TH1F* h_ffDiJetPt_2Jet_new=(TH1F*)h_ffDiJetPt_2Jet->Rebin(numBins,"h_ffDiJetPt_2Jet_new",DiEMPtBins);

  TH1F* h_ggDiEMPt_JetReq_new=(TH1F*)h_ggDiEMPt_JetReq->Rebin(numBins,"h_ggDiEMPt_JetReq_new",DiEMPtBins);
  TH1F* h_ggDiJetPt_JetReq_new=(TH1F*)h_ggDiJetPt_JetReq->Rebin(numBins,"h_ggDiJetPt_JetReq_new",DiEMPtBins);
  TH1F* h_ggDiEMPt_JetReq_0Jet_new=(TH1F*)h_ggDiEMPt_JetReq_0Jet->Rebin(numBins,"h_ggDiEMPt_JetReq_0Jet_new",DiEMPtBins);
  TH1F* h_ggDiEMPt_JetReq_1Jet_new=(TH1F*)h_ggDiEMPt_JetReq_1Jet->Rebin(numBins,"h_ggDiEMPt_JetReq_1Jet_new",DiEMPtBins);
  TH1F* h_ggDiEMPt_JetReq_2Jet_new=(TH1F*)h_ggDiEMPt_JetReq_2Jet->Rebin(numBins,"h_ggDiEMPt_JetReq_2Jet_new",DiEMPtBins);
  TH1F* h_ggDiJetPt_JetReq_0Jet_new=(TH1F*)h_ggDiJetPt_JetReq_0Jet->Rebin(numBins,"h_ggDiJetPt_JetReq_0Jet_new",DiEMPtBins);
  TH1F* h_ggDiJetPt_JetReq_1Jet_new=(TH1F*)h_ggDiJetPt_JetReq_1Jet->Rebin(numBins,"h_ggDiJetPt_JetReq_1Jet_new",DiEMPtBins);
  TH1F* h_ggDiJetPt_JetReq_2Jet_new=(TH1F*)h_ggDiJetPt_JetReq_2Jet->Rebin(numBins,"h_ggDiJetPt_JetReq_2Jet_new",DiEMPtBins);
  TH1F* h_eeDiEMPt_JetReq_new=(TH1F*)h_eeDiEMPt_JetReq->Rebin(numBins,"h_eeDiEMPt_JetReq_new",DiEMPtBins);
  TH1F* h_eeDiJetPt_JetReq_new=(TH1F*)h_eeDiJetPt_JetReq->Rebin(numBins,"h_eeDiJetPt_JetReq_new",DiEMPtBins);
  TH1F* h_eeDiEMPt_JetReq_0Jet_new=(TH1F*)h_eeDiEMPt_JetReq_0Jet->Rebin(numBins,"h_eeDiEMPt_JetReq_0Jet_new",DiEMPtBins);
  TH1F* h_eeDiEMPt_JetReq_1Jet_new=(TH1F*)h_eeDiEMPt_JetReq_1Jet->Rebin(numBins,"h_eeDiEMPt_JetReq_1Jet_new",DiEMPtBins);
  TH1F* h_eeDiEMPt_JetReq_2Jet_new=(TH1F*)h_eeDiEMPt_JetReq_2Jet->Rebin(numBins,"h_eeDiEMPt_JetReq_2Jet_new",DiEMPtBins);
  TH1F* h_eeDiJetPt_JetReq_0Jet_new=(TH1F*)h_eeDiJetPt_JetReq_0Jet->Rebin(numBins,"h_eeDiJetPt_JetReq_0Jet_new",DiEMPtBins);
  TH1F* h_eeDiJetPt_JetReq_1Jet_new=(TH1F*)h_eeDiJetPt_JetReq_1Jet->Rebin(numBins,"h_eeDiJetPt_JetReq_1Jet_new",DiEMPtBins);
  TH1F* h_eeDiJetPt_JetReq_2Jet_new=(TH1F*)h_eeDiJetPt_JetReq_2Jet->Rebin(numBins,"h_eeDiJetPt_JetReq_2Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandDiEMPt_JetReq_new=(TH1F*)h_eeSideBandDiEMPt_JetReq->Rebin(numBins,"h_eeSideBandDiEMPt_JetReq_new",DiEMPtBins);
  TH1F* h_eeSideBandDiJetPt_JetReq_new=(TH1F*)h_eeSideBandDiJetPt_JetReq->Rebin(numBins,"h_eeSideBandDiJetPt_JetReq_new",DiEMPtBins);
  TH1F* h_eeSideBandLowDiEMPt_JetReq_new=(TH1F*)h_eeSideBandLowDiEMPt_JetReq->Rebin(numBins,"h_eeSideBandLowDiEMPt_JetReq_new",DiEMPtBins);
  TH1F* h_eeSideBandHighDiEMPt_JetReq_new=(TH1F*)h_eeSideBandHighDiEMPt_JetReq->Rebin(numBins,"h_eeSideBandHighDiEMPt_JetReq_new",DiEMPtBins);
  TH1F* h_eeSideBandLowDiEMPt_JetReq_0Jet_new=(TH1F*)h_eeSideBandLowDiEMPt_JetReq_0Jet->Rebin(numBins,"h_eeSideBandLowDiEMPt_JetReq_0Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandLowDiEMPt_JetReq_1Jet_new=(TH1F*)h_eeSideBandLowDiEMPt_JetReq_1Jet->Rebin(numBins,"h_eeSideBandLowDiEMPt_JetReq_1Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandLowDiEMPt_JetReq_2Jet_new=(TH1F*)h_eeSideBandLowDiEMPt_JetReq_2Jet->Rebin(numBins,"h_eeSideBandLowDiEMPt_JetReq_2Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandHighDiEMPt_JetReq_0Jet_new=(TH1F*)h_eeSideBandHighDiEMPt_JetReq_0Jet->Rebin(numBins,"h_eeSideBandHighDiEMPt_JetReq_0Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandHighDiEMPt_JetReq_1Jet_new=(TH1F*)h_eeSideBandHighDiEMPt_JetReq_1Jet->Rebin(numBins,"h_eeSideBandHighDiEMPt_JetReq_1Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandHighDiEMPt_JetReq_2Jet_new=(TH1F*)h_eeSideBandHighDiEMPt_JetReq_2Jet->Rebin(numBins,"h_eeSideBandHighDiEMPt_JetReq_2Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandLowDiJetPt_JetReq_new=(TH1F*)h_eeSideBandLowDiJetPt_JetReq->Rebin(numBins,"h_eeSideBandLowDiJetPt_JetReq_new",DiEMPtBins);
  TH1F* h_eeSideBandHighDiJetPt_JetReq_new=(TH1F*)h_eeSideBandHighDiJetPt_JetReq->Rebin(numBins,"h_eeSideBandHighDiJetPt_JetReq_new",DiEMPtBins);
  TH1F* h_eeSideBandLowDiJetPt_JetReq_0Jet_new=(TH1F*)h_eeSideBandLowDiJetPt_JetReq_0Jet->Rebin(numBins,"h_eeSideBandLowDiJetPt_JetReq_0Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandLowDiJetPt_JetReq_1Jet_new=(TH1F*)h_eeSideBandLowDiJetPt_JetReq_1Jet->Rebin(numBins,"h_eeSideBandLowDiJetPt_JetReq_1Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandLowDiJetPt_JetReq_2Jet_new=(TH1F*)h_eeSideBandLowDiJetPt_JetReq_2Jet->Rebin(numBins,"h_eeSideBandLowDiJetPt_JetReq_2Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandHighDiJetPt_JetReq_0Jet_new=(TH1F*)h_eeSideBandHighDiJetPt_JetReq_0Jet->Rebin(numBins,"h_eeSideBandHighDiJetPt_JetReq_0Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandHighDiJetPt_JetReq_1Jet_new=(TH1F*)h_eeSideBandHighDiJetPt_JetReq_1Jet->Rebin(numBins,"h_eeSideBandHighDiJetPt_JetReq_1Jet_new",DiEMPtBins);
  TH1F* h_eeSideBandHighDiJetPt_JetReq_2Jet_new=(TH1F*)h_eeSideBandHighDiJetPt_JetReq_2Jet->Rebin(numBins,"h_eeSideBandHighDiJetPt_JetReq_2Jet_new",DiEMPtBins);
  TH1F* h_ffDiEMPt_JetReq_new=(TH1F*)h_ffDiEMPt_JetReq->Rebin(numBins,"h_ffDiEMPt_JetReq_new",DiEMPtBins);
  TH1F* h_ffDiJetPt_JetReq_new=(TH1F*)h_ffDiJetPt_JetReq->Rebin(numBins,"h_ffDiJetPt_JetReq_new",DiEMPtBins);
  TH1F* h_ffDiEMPt_JetReq_0Jet_new=(TH1F*)h_ffDiEMPt_JetReq_0Jet->Rebin(numBins,"h_ffDiEMPt_JetReq_0Jet_new",DiEMPtBins);
  TH1F* h_ffDiEMPt_JetReq_1Jet_new=(TH1F*)h_ffDiEMPt_JetReq_1Jet->Rebin(numBins,"h_ffDiEMPt_JetReq_1Jet_new",DiEMPtBins);
  TH1F* h_ffDiEMPt_JetReq_2Jet_new=(TH1F*)h_ffDiEMPt_JetReq_2Jet->Rebin(numBins,"h_ffDiEMPt_JetReq_2Jet_new",DiEMPtBins);
  TH1F* h_ffDiJetPt_JetReq_0Jet_new=(TH1F*)h_ffDiJetPt_JetReq_0Jet->Rebin(numBins,"h_ffDiJetPt_JetReq_0Jet_new",DiEMPtBins);
  TH1F* h_ffDiJetPt_JetReq_1Jet_new=(TH1F*)h_ffDiJetPt_JetReq_1Jet->Rebin(numBins,"h_ffDiJetPt_JetReq_1Jet_new",DiEMPtBins);
  TH1F* h_ffDiJetPt_JetReq_2Jet_new=(TH1F*)h_ffDiJetPt_JetReq_2Jet->Rebin(numBins,"h_ffDiJetPt_JetReq_2Jet_new",DiEMPtBins);



  float ggInt=h_ggDiEMPt_new->Integral(),eeInt=h_eeDiEMPt_new->Integral(),eeSideBandInt=h_eeSideBandDiEMPt_new->Integral(),ffInt=h_ffDiEMPt_new->Integral();
  float eeSideBandLowInt=h_eeSideBandLowDiEMPt_new->Integral(),eeSideBandHighInt=h_eeSideBandHighDiEMPt_new->Integral();
  float eeSideBandLowInt_0Jet=h_eeSideBandLowDiEMPt_0Jet_new->Integral(),eeSideBandLowInt_1Jet=h_eeSideBandLowDiEMPt_1Jet_new->Integral(),eeSideBandLowInt_2Jet=h_eeSideBandLowDiEMPt_2Jet_new->Integral();
  float eeSideBandHighInt_0Jet=h_eeSideBandHighDiEMPt_0Jet_new->Integral(),eeSideBandHighInt_1Jet=h_eeSideBandHighDiEMPt_1Jet_new->Integral(),eeSideBandHighInt_2Jet=h_eeSideBandHighDiEMPt_2Jet_new->Integral();
  float eeSideBandLowJetInt=h_eeSideBandLowDiJetPt_new->Integral(),eeSideBandHighJetInt=h_eeSideBandHighDiJetPt_new->Integral();
  float eeSideBandLowJetInt_0Jet=h_eeSideBandLowDiJetPt_0Jet_new->Integral(),eeSideBandLowJetInt_1Jet=h_eeSideBandLowDiJetPt_1Jet_new->Integral(),eeSideBandLowJetInt_2Jet=h_eeSideBandLowDiJetPt_2Jet_new->Integral();
  float eeSideBandHighJetInt_0Jet=h_eeSideBandHighDiJetPt_0Jet_new->Integral(),eeSideBandHighJetInt_1Jet=h_eeSideBandHighDiJetPt_1Jet_new->Integral(),eeSideBandHighJetInt_2Jet=h_eeSideBandHighDiJetPt_2Jet_new->Integral();
  float ggJetInt=h_ggDiJetPt_new->Integral(),eeJetInt=h_eeDiJetPt_new->Integral(),ffJetInt=h_ffDiJetPt_new->Integral();
  float ggInt_0Jet=h_ggDiEMPt_0Jet_new->Integral(),eeInt_0Jet=h_eeDiEMPt_0Jet_new->Integral(),ffInt_0Jet=h_ffDiEMPt_0Jet_new->Integral();
  float ggInt_1Jet=h_ggDiEMPt_1Jet_new->Integral(),eeInt_1Jet=h_eeDiEMPt_1Jet_new->Integral(),ffInt_1Jet=h_ffDiEMPt_1Jet_new->Integral();
  float ggInt_2Jet=h_ggDiEMPt_2Jet_new->Integral(),eeInt_2Jet=h_eeDiEMPt_2Jet_new->Integral(),ffInt_2Jet=h_ffDiEMPt_2Jet_new->Integral();
  float ggJetInt_0Jet=h_ggDiJetPt_0Jet_new->Integral(),eeJetInt_0Jet=h_eeDiJetPt_0Jet_new->Integral(),ffJetInt_0Jet=h_ffDiJetPt_0Jet_new->Integral();
  float ggJetInt_1Jet=h_ggDiJetPt_1Jet_new->Integral(),eeJetInt_1Jet=h_eeDiJetPt_1Jet_new->Integral(),ffJetInt_1Jet=h_ffDiJetPt_1Jet_new->Integral();
  float ggJetInt_2Jet=h_ggDiJetPt_2Jet_new->Integral(),eeJetInt_2Jet=h_eeDiJetPt_2Jet_new->Integral(),ffJetInt_2Jet=h_ffDiJetPt_2Jet_new->Integral();
  float scaleggee=ggInt/eeInt;
  float scaleggeeSideBand=ggInt/eeSideBandInt,scaleggeeSideBandLow=ggInt/eeSideBandLowInt,scaleggeeSideBandHigh=ggInt/eeSideBandHighInt;
  float scaleggeeSideBandLowJet=ggJetInt/eeSideBandLowJetInt,scaleggeeSideBandHighJet=ggJetInt/eeSideBandHighJetInt;
  float scaleggff=ggInt/ffInt;
  float scaleggeeJet=ggJetInt/eeJetInt,scaleggffJet=ggJetInt/ffJetInt;

  float ggInt_JetReq=h_ggDiEMPt_JetReq_new->Integral(),eeInt_JetReq=h_eeDiEMPt_JetReq_new->Integral(),eeSideBandInt_JetReq=h_eeSideBandDiEMPt_JetReq_new->Integral(),ffInt_JetReq=h_ffDiEMPt_JetReq_new->Integral();
  float eeSideBandLowInt_JetReq=h_eeSideBandLowDiEMPt_JetReq_new->Integral(),eeSideBandHighInt_JetReq=h_eeSideBandHighDiEMPt_JetReq_new->Integral();
  float eeSideBandLowInt_JetReq_0Jet=h_eeSideBandLowDiEMPt_JetReq_0Jet_new->Integral(),eeSideBandLowInt_JetReq_1Jet=h_eeSideBandLowDiEMPt_JetReq_1Jet_new->Integral(),eeSideBandLowInt_JetReq_2Jet=h_eeSideBandLowDiEMPt_JetReq_2Jet_new->Integral();
  float eeSideBandHighInt_JetReq_0Jet=h_eeSideBandHighDiEMPt_JetReq_0Jet_new->Integral(),eeSideBandHighInt_JetReq_1Jet=h_eeSideBandHighDiEMPt_JetReq_1Jet_new->Integral(),eeSideBandHighInt_JetReq_2Jet=h_eeSideBandHighDiEMPt_JetReq_2Jet_new->Integral();
  float eeSideBandLowJetInt_JetReq=h_eeSideBandLowDiJetPt_JetReq_new->Integral(),eeSideBandHighJetInt_JetReq=h_eeSideBandHighDiJetPt_JetReq_new->Integral();
  float eeSideBandLowJetInt_JetReq_0Jet=h_eeSideBandLowDiJetPt_JetReq_0Jet_new->Integral(),eeSideBandLowJetInt_JetReq_1Jet=h_eeSideBandLowDiJetPt_JetReq_1Jet_new->Integral(),eeSideBandLowJetInt_JetReq_2Jet=h_eeSideBandLowDiJetPt_JetReq_2Jet_new->Integral();
  float eeSideBandHighJetInt_JetReq_0Jet=h_eeSideBandHighDiJetPt_JetReq_0Jet_new->Integral(),eeSideBandHighJetInt_JetReq_1Jet=h_eeSideBandHighDiJetPt_JetReq_1Jet_new->Integral(),eeSideBandHighJetInt_JetReq_2Jet=h_eeSideBandHighDiJetPt_JetReq_2Jet_new->Integral();
  float ggJetInt_JetReq=h_ggDiJetPt_JetReq_new->Integral(),eeJetInt_JetReq=h_eeDiJetPt_JetReq_new->Integral(),ffJetInt_JetReq=h_ffDiJetPt_JetReq_new->Integral();
  float ggInt_JetReq_0Jet=h_ggDiEMPt_JetReq_0Jet_new->Integral(),eeInt_JetReq_0Jet=h_eeDiEMPt_JetReq_0Jet_new->Integral(),ffInt_JetReq_0Jet=h_ffDiEMPt_JetReq_0Jet_new->Integral();
  float ggInt_JetReq_1Jet=h_ggDiEMPt_JetReq_1Jet_new->Integral(),eeInt_JetReq_1Jet=h_eeDiEMPt_JetReq_1Jet_new->Integral(),ffInt_JetReq_1Jet=h_ffDiEMPt_JetReq_1Jet_new->Integral();
  float ggInt_JetReq_2Jet=h_ggDiEMPt_JetReq_2Jet_new->Integral(),eeInt_JetReq_2Jet=h_eeDiEMPt_JetReq_2Jet_new->Integral(),ffInt_JetReq_2Jet=h_ffDiEMPt_JetReq_2Jet_new->Integral();
  float ggJetInt_JetReq_0Jet=h_ggDiJetPt_JetReq_0Jet_new->Integral(),eeJetInt_JetReq_0Jet=h_eeDiJetPt_JetReq_0Jet_new->Integral(),ffJetInt_JetReq_0Jet=h_ffDiJetPt_JetReq_0Jet_new->Integral();
  float ggJetInt_JetReq_1Jet=h_ggDiJetPt_JetReq_1Jet_new->Integral(),eeJetInt_JetReq_1Jet=h_eeDiJetPt_JetReq_1Jet_new->Integral(),ffJetInt_JetReq_1Jet=h_ffDiJetPt_JetReq_1Jet_new->Integral();
  float ggJetInt_JetReq_2Jet=h_ggDiJetPt_JetReq_2Jet_new->Integral(),eeJetInt_JetReq_2Jet=h_eeDiJetPt_JetReq_2Jet_new->Integral(),ffJetInt_JetReq_2Jet=h_ffDiJetPt_JetReq_2Jet_new->Integral();
  float scaleggee_JetReq=ggInt_JetReq/eeInt;
  float scaleggeeSideBand_JetReq=ggInt_JetReq/eeSideBandInt_JetReq,scaleggeeSideBandLow_JetReq=ggInt_JetReq/eeSideBandLowInt_JetReq,scaleggeeSideBandHigh_JetReq=ggInt_JetReq/eeSideBandHighInt_JetReq;
  float scaleggeeSideBandLowJet_JetReq=ggJetInt_JetReq/eeSideBandLowJetInt_JetReq,scaleggeeSideBandHighJet_JetReq=ggJetInt_JetReq/eeSideBandHighJetInt_JetReq;
  float scaleggff_JetReq=ggInt_JetReq/ffInt_JetReq;
  float scaleggeeJet_JetReq=ggJetInt_JetReq/eeJetInt_JetReq,scaleggffJet_JetReq=ggJetInt_JetReq/ffJetInt_JetReq;

  h_eeDiEMPt_new->Scale(scaleggee);
  h_eeSideBandDiEMPt_new->Scale(scaleggeeSideBand);
  h_eeSideBandDiJetPt_new->Scale(scaleggeeSideBand);
  h_eeSideBandLowDiEMPt_0Jet_new->Scale(scaleggeeSideBandLow);h_eeSideBandLowDiEMPt_1Jet_new->Scale(scaleggeeSideBandLow);h_eeSideBandLowDiEMPt_2Jet_new->Scale(scaleggeeSideBandLow);
  h_eeSideBandHighDiEMPt_0Jet_new->Scale(scaleggeeSideBandHigh);h_eeSideBandHighDiEMPt_1Jet_new->Scale(scaleggeeSideBandHigh);h_eeSideBandHighDiEMPt_2Jet_new->Scale(scaleggeeSideBandHigh);
  h_eeSideBandLowDiJetPt_0Jet_new->Scale(scaleggeeSideBandLowJet);h_eeSideBandLowDiJetPt_1Jet_new->Scale(scaleggeeSideBandLowJet);h_eeSideBandLowDiJetPt_2Jet_new->Scale(scaleggeeSideBandLowJet);
  h_eeSideBandHighDiJetPt_0Jet_new->Scale(scaleggeeSideBandHighJet);h_eeSideBandHighDiJetPt_1Jet_new->Scale(scaleggeeSideBandHighJet);h_eeSideBandHighDiJetPt_2Jet_new->Scale(scaleggeeSideBandHighJet);
  h_ffDiEMPt_new->Scale(scaleggff);
  h_eeDiJetPt_new->Scale(scaleggeeJet);
  h_eeDiEMPt_0Jet_new->Scale(scaleggee);h_eeDiEMPt_1Jet_new->Scale(scaleggee);h_eeDiEMPt_2Jet_new->Scale(scaleggee);//need to multiply by full ratios
  h_eeDiJetPt_0Jet_new->Scale(scaleggeeJet);h_eeDiJetPt_1Jet_new->Scale(scaleggeeJet);h_eeDiJetPt_2Jet_new->Scale(scaleggeeJet);//need to multiply by full ratios
  h_ffDiJetPt_new->Scale(scaleggffJet);
  h_ffDiEMPt_0Jet_new->Scale(scaleggff);h_ffDiEMPt_1Jet_new->Scale(scaleggff);h_ffDiEMPt_2Jet_new->Scale(scaleggff);
  h_ffDiJetPt_0Jet_new->Scale(scaleggffJet);h_ffDiJetPt_1Jet_new->Scale(scaleggffJet);h_ffDiJetPt_2Jet_new->Scale(scaleggffJet);
  //JetReq
  h_eeDiEMPt_JetReq_new->Scale(scaleggee_JetReq);
  h_eeSideBandDiEMPt_JetReq_new->Scale(scaleggeeSideBand_JetReq);
  h_eeSideBandDiJetPt_JetReq_new->Scale(scaleggeeSideBand_JetReq);
  h_eeSideBandLowDiEMPt_JetReq_0Jet_new->Scale(scaleggeeSideBandLow_JetReq);h_eeSideBandLowDiEMPt_JetReq_1Jet_new->Scale(scaleggeeSideBandLow_JetReq);h_eeSideBandLowDiEMPt_JetReq_2Jet_new->Scale(scaleggeeSideBandLow_JetReq);
  h_eeSideBandHighDiEMPt_JetReq_0Jet_new->Scale(scaleggeeSideBandHigh_JetReq);h_eeSideBandHighDiEMPt_JetReq_1Jet_new->Scale(scaleggeeSideBandHigh_JetReq);h_eeSideBandHighDiEMPt_JetReq_2Jet_new->Scale(scaleggeeSideBandHigh_JetReq);
  h_eeSideBandLowDiJetPt_JetReq_0Jet_new->Scale(scaleggeeSideBandLowJet_JetReq);h_eeSideBandLowDiJetPt_JetReq_1Jet_new->Scale(scaleggeeSideBandLowJet_JetReq);h_eeSideBandLowDiJetPt_JetReq_2Jet_new->Scale(scaleggeeSideBandLowJet_JetReq);
  h_eeSideBandHighDiJetPt_JetReq_0Jet_new->Scale(scaleggeeSideBandHighJet_JetReq);h_eeSideBandHighDiJetPt_JetReq_1Jet_new->Scale(scaleggeeSideBandHighJet_JetReq);h_eeSideBandHighDiJetPt_JetReq_2Jet_new->Scale(scaleggeeSideBandHighJet_JetReq);
  h_ffDiEMPt_JetReq_new->Scale(scaleggff_JetReq);
  h_eeDiJetPt_JetReq_new->Scale(scaleggeeJet_JetReq);
  h_eeDiEMPt_JetReq_0Jet_new->Scale(scaleggee_JetReq);h_eeDiEMPt_JetReq_1Jet_new->Scale(scaleggee_JetReq);h_eeDiEMPt_JetReq_2Jet_new->Scale(scaleggee_JetReq);//need to multiply by full ratios
  h_eeDiJetPt_JetReq_0Jet_new->Scale(scaleggeeJet_JetReq);h_eeDiJetPt_JetReq_1Jet_new->Scale(scaleggeeJet_JetReq);h_eeDiJetPt_JetReq_2Jet_new->Scale(scaleggeeJet_JetReq);//need to multiply by full ratios
  h_ffDiJetPt_JetReq_new->Scale(scaleggffJet_JetReq);
  h_ffDiEMPt_JetReq_0Jet_new->Scale(scaleggff_JetReq);h_ffDiEMPt_JetReq_1Jet_new->Scale(scaleggff_JetReq);h_ffDiEMPt_JetReq_2Jet_new->Scale(scaleggff_JetReq);
  h_ffDiJetPt_JetReq_0Jet_new->Scale(scaleggffJet_JetReq);h_ffDiJetPt_JetReq_1Jet_new->Scale(scaleggffJet_JetReq);h_ffDiJetPt_JetReq_2Jet_new->Scale(scaleggffJet_JetReq);
  


  h_ggeeDiEMPtRatio->Divide(h_ggDiEMPt_new,h_eeDiEMPt_new,1,1,"");
  h_ggeeSideBandDiEMPtRatio->Divide(h_ggDiEMPt_new,h_eeSideBandDiEMPt_new,1,1,"");
  h_ggeeSideBandDiJetPtRatio->Divide(h_ggDiJetPt_new,h_eeSideBandDiJetPt_new,1,1,"");
  h_ggeeSideBandLowDiEMPtRatio->Divide(h_ggDiEMPt_new,h_eeSideBandLowDiEMPt_new,1,1,"");
  h_ggeeSideBandHighDiEMPtRatio->Divide(h_ggDiEMPt_new,h_eeSideBandHighDiEMPt_new,1,1,"");
  h_ggeeSideBandLowDiEMPtRatio_0Jet->Divide(h_ggDiEMPt_0Jet_new,h_eeSideBandLowDiEMPt_0Jet_new,1,1,"");
  h_ggeeSideBandLowDiEMPtRatio_1Jet->Divide(h_ggDiEMPt_1Jet_new,h_eeSideBandLowDiEMPt_1Jet_new,1,1,"");
  h_ggeeSideBandLowDiEMPtRatio_2Jet->Divide(h_ggDiEMPt_2Jet_new,h_eeSideBandLowDiEMPt_2Jet_new,1,1,"");
  h_ggeeSideBandHighDiEMPtRatio_0Jet->Divide(h_ggDiEMPt_0Jet_new,h_eeSideBandHighDiEMPt_0Jet_new,1,1,"");
  h_ggeeSideBandHighDiEMPtRatio_1Jet->Divide(h_ggDiEMPt_1Jet_new,h_eeSideBandHighDiEMPt_1Jet_new,1,1,"");
  h_ggeeSideBandHighDiEMPtRatio_2Jet->Divide(h_ggDiEMPt_2Jet_new,h_eeSideBandHighDiEMPt_2Jet_new,1,1,"");
  h_ggeeSideBandLowDiJetPtRatio->Divide(h_ggDiJetPt_new,h_eeSideBandLowDiJetPt_new,1,1,"");
  h_ggeeSideBandHighDiJetPtRatio->Divide(h_ggDiJetPt_new,h_eeSideBandHighDiJetPt_new,1,1,"");
  h_ggeeSideBandLowDiJetPtRatio_0Jet->Divide(h_ggDiJetPt_0Jet_new,h_eeSideBandLowDiJetPt_0Jet_new,1,1,"");
  h_ggeeSideBandLowDiJetPtRatio_1Jet->Divide(h_ggDiJetPt_1Jet_new,h_eeSideBandLowDiJetPt_1Jet_new,1,1,"");
  h_ggeeSideBandLowDiJetPtRatio_2Jet->Divide(h_ggDiJetPt_2Jet_new,h_eeSideBandLowDiJetPt_2Jet_new,1,1,"");
  h_ggeeSideBandHighDiJetPtRatio_0Jet->Divide(h_ggDiJetPt_0Jet_new,h_eeSideBandHighDiJetPt_0Jet_new,1,1,"");
  h_ggeeSideBandHighDiJetPtRatio_1Jet->Divide(h_ggDiJetPt_1Jet_new,h_eeSideBandHighDiJetPt_1Jet_new,1,1,"");
  h_ggeeSideBandHighDiJetPtRatio_2Jet->Divide(h_ggDiJetPt_2Jet_new,h_eeSideBandHighDiJetPt_2Jet_new,1,1,"");
  h_ggffDiEMPtRatio->Divide(h_ggDiEMPt_new,h_ffDiEMPt_new,1,1,"");
  h_ggeeDiJetPtRatio->Divide(h_ggDiJetPt_new,h_eeDiJetPt_new,1,1,"");
  h_ggeeDiEMPtRatio_0Jet->Divide(h_ggDiEMPt_0Jet_new,h_eeDiEMPt_0Jet_new,1,1,"");
  h_ggeeDiEMPtRatio_1Jet->Divide(h_ggDiEMPt_1Jet_new,h_eeDiEMPt_1Jet_new,1,1,"");
  h_ggeeDiEMPtRatio_2Jet->Divide(h_ggDiEMPt_2Jet_new,h_eeDiEMPt_2Jet_new,1,1,"");
  h_ggeeDiJetPtRatio_0Jet->Divide(h_ggDiJetPt_0Jet_new,h_eeDiJetPt_0Jet_new,1,1,"");
  h_ggeeDiJetPtRatio_1Jet->Divide(h_ggDiJetPt_1Jet_new,h_eeDiJetPt_1Jet_new,1,1,"");
  h_ggeeDiJetPtRatio_2Jet->Divide(h_ggDiJetPt_2Jet_new,h_eeDiJetPt_2Jet_new,1,1,"");
  h_ggffDiJetPtRatio->Divide(h_ggDiJetPt_new,h_ffDiJetPt_new,1,1,"");
  h_ggffDiEMPtRatio_0Jet->Divide(h_ggDiEMPt_0Jet_new,h_ffDiEMPt_0Jet_new,1,1,"");
  h_ggffDiEMPtRatio_1Jet->Divide(h_ggDiEMPt_1Jet_new,h_ffDiEMPt_1Jet_new,1,1,"");
  h_ggffDiEMPtRatio_2Jet->Divide(h_ggDiEMPt_2Jet_new,h_ffDiEMPt_2Jet_new,1,1,"");
  h_ggffDiJetPtRatio_0Jet->Divide(h_ggDiJetPt_0Jet_new,h_ffDiJetPt_0Jet_new,1,1,"");
  h_ggffDiJetPtRatio_1Jet->Divide(h_ggDiJetPt_1Jet_new,h_ffDiJetPt_1Jet_new,1,1,"");
  h_ggffDiJetPtRatio_2Jet->Divide(h_ggDiJetPt_2Jet_new,h_ffDiJetPt_2Jet_new,1,1,"");

  h_ggeeDiEMPtRatio_JetReq->Divide(h_ggDiEMPt_JetReq_new,h_eeDiEMPt_JetReq_new,1,1,"");
  h_ggeeSideBandDiEMPtRatio_JetReq->Divide(h_ggDiEMPt_JetReq_new,h_eeSideBandDiEMPt_JetReq_new,1,1,"");
  h_ggeeSideBandDiJetPtRatio_JetReq->Divide(h_ggDiJetPt_JetReq_new,h_eeSideBandDiJetPt_JetReq_new,1,1,"");
  h_ggeeSideBandLowDiEMPtRatio_JetReq->Divide(h_ggDiEMPt_JetReq_new,h_eeSideBandLowDiEMPt_JetReq_new,1,1,"");
  h_ggeeSideBandHighDiEMPtRatio_JetReq->Divide(h_ggDiEMPt_JetReq_new,h_eeSideBandHighDiEMPt_JetReq_new,1,1,"");
  h_ggeeSideBandLowDiEMPtRatio_JetReq_0Jet->Divide(h_ggDiEMPt_JetReq_0Jet_new,h_eeSideBandLowDiEMPt_JetReq_0Jet_new,1,1,"");
  h_ggeeSideBandLowDiEMPtRatio_JetReq_1Jet->Divide(h_ggDiEMPt_JetReq_1Jet_new,h_eeSideBandLowDiEMPt_JetReq_1Jet_new,1,1,"");
  h_ggeeSideBandLowDiEMPtRatio_JetReq_2Jet->Divide(h_ggDiEMPt_JetReq_2Jet_new,h_eeSideBandLowDiEMPt_JetReq_2Jet_new,1,1,"");
  h_ggeeSideBandHighDiEMPtRatio_JetReq_0Jet->Divide(h_ggDiEMPt_JetReq_0Jet_new,h_eeSideBandHighDiEMPt_JetReq_0Jet_new,1,1,"");
  h_ggeeSideBandHighDiEMPtRatio_JetReq_1Jet->Divide(h_ggDiEMPt_JetReq_1Jet_new,h_eeSideBandHighDiEMPt_JetReq_1Jet_new,1,1,"");
  h_ggeeSideBandHighDiEMPtRatio_JetReq_2Jet->Divide(h_ggDiEMPt_JetReq_2Jet_new,h_eeSideBandHighDiEMPt_JetReq_2Jet_new,1,1,"");
  h_ggeeSideBandLowDiJetPtRatio_JetReq->Divide(h_ggDiJetPt_JetReq_new,h_eeSideBandLowDiJetPt_JetReq_new,1,1,"");
  h_ggeeSideBandHighDiJetPtRatio_JetReq->Divide(h_ggDiJetPt_JetReq_new,h_eeSideBandHighDiJetPt_JetReq_new,1,1,"");
  h_ggeeSideBandLowDiJetPtRatio_JetReq_0Jet->Divide(h_ggDiJetPt_JetReq_0Jet_new,h_eeSideBandLowDiJetPt_JetReq_0Jet_new,1,1,"");
  h_ggeeSideBandLowDiJetPtRatio_JetReq_1Jet->Divide(h_ggDiJetPt_JetReq_1Jet_new,h_eeSideBandLowDiJetPt_JetReq_1Jet_new,1,1,"");
  h_ggeeSideBandLowDiJetPtRatio_JetReq_2Jet->Divide(h_ggDiJetPt_JetReq_2Jet_new,h_eeSideBandLowDiJetPt_JetReq_2Jet_new,1,1,"");
  h_ggeeSideBandHighDiJetPtRatio_JetReq_0Jet->Divide(h_ggDiJetPt_JetReq_0Jet_new,h_eeSideBandHighDiJetPt_JetReq_0Jet_new,1,1,"");
  h_ggeeSideBandHighDiJetPtRatio_JetReq_1Jet->Divide(h_ggDiJetPt_JetReq_1Jet_new,h_eeSideBandHighDiJetPt_JetReq_1Jet_new,1,1,"");
  h_ggeeSideBandHighDiJetPtRatio_JetReq_2Jet->Divide(h_ggDiJetPt_JetReq_2Jet_new,h_eeSideBandHighDiJetPt_JetReq_2Jet_new,1,1,"");
  h_ggffDiEMPtRatio_JetReq->Divide(h_ggDiEMPt_JetReq_new,h_ffDiEMPt_JetReq_new,1,1,"");
  h_ggeeDiJetPtRatio_JetReq->Divide(h_ggDiJetPt_JetReq_new,h_eeDiJetPt_JetReq_new,1,1,"");
  h_ggeeDiEMPtRatio_JetReq_0Jet->Divide(h_ggDiEMPt_JetReq_0Jet_new,h_eeDiEMPt_JetReq_0Jet_new,1,1,"");
  h_ggeeDiEMPtRatio_JetReq_1Jet->Divide(h_ggDiEMPt_JetReq_1Jet_new,h_eeDiEMPt_JetReq_1Jet_new,1,1,"");
  h_ggeeDiEMPtRatio_JetReq_2Jet->Divide(h_ggDiEMPt_JetReq_2Jet_new,h_eeDiEMPt_JetReq_2Jet_new,1,1,"");
  h_ggeeDiJetPtRatio_JetReq_0Jet->Divide(h_ggDiJetPt_JetReq_0Jet_new,h_eeDiJetPt_JetReq_0Jet_new,1,1,"");
  h_ggeeDiJetPtRatio_JetReq_1Jet->Divide(h_ggDiJetPt_JetReq_1Jet_new,h_eeDiJetPt_JetReq_1Jet_new,1,1,"");
  h_ggeeDiJetPtRatio_JetReq_2Jet->Divide(h_ggDiJetPt_JetReq_2Jet_new,h_eeDiJetPt_JetReq_2Jet_new,1,1,"");
  h_ggffDiJetPtRatio_JetReq->Divide(h_ggDiJetPt_JetReq_new,h_ffDiJetPt_JetReq_new,1,1,"");
  h_ggffDiEMPtRatio_JetReq_0Jet->Divide(h_ggDiEMPt_JetReq_0Jet_new,h_ffDiEMPt_JetReq_0Jet_new,1,1,"");
  h_ggffDiEMPtRatio_JetReq_1Jet->Divide(h_ggDiEMPt_JetReq_1Jet_new,h_ffDiEMPt_JetReq_1Jet_new,1,1,"");
  h_ggffDiEMPtRatio_JetReq_2Jet->Divide(h_ggDiEMPt_JetReq_2Jet_new,h_ffDiEMPt_JetReq_2Jet_new,1,1,"");
  h_ggffDiJetPtRatio_JetReq_0Jet->Divide(h_ggDiJetPt_JetReq_0Jet_new,h_ffDiJetPt_JetReq_0Jet_new,1,1,"");
  h_ggffDiJetPtRatio_JetReq_1Jet->Divide(h_ggDiJetPt_JetReq_1Jet_new,h_ffDiJetPt_JetReq_1Jet_new,1,1,"");
  h_ggffDiJetPtRatio_JetReq_2Jet->Divide(h_ggDiJetPt_JetReq_2Jet_new,h_ffDiJetPt_JetReq_2Jet_new,1,1,"");




  //-----make background efficiency plots------------
  if(!event->isRealData){
    TGraphAsymmErrors *BgroundEffVsNVertex = new TGraphAsymmErrors(h_NumerVsNVertex,h_DenomVsNVertex,"");
    BgroundEffVsNVertex->SetMarkerSize(0.6);
    BgroundEffVsNVertex->GetXaxis()->SetTitle("NVertex");
    BgroundEffVsNVertex->GetYaxis()->SetTitle("Efficiency");
    BgroundEffVsNVertex->Write("BgroundEffVsNVertex");
    TGraphAsymmErrors *RhoCorrectedBgroundEffVsNVertex = new TGraphAsymmErrors(h_RhoCorrectedNumerVsNVertex,h_DenomVsNVertex,"");
    RhoCorrectedBgroundEffVsNVertex->SetMarkerSize(0.6);
    RhoCorrectedBgroundEffVsNVertex->GetXaxis()->SetTitle("NVertex");
    RhoCorrectedBgroundEffVsNVertex->GetYaxis()->SetTitle("Efficiency");
    RhoCorrectedBgroundEffVsNVertex->SetLineColor(kRed);
    RhoCorrectedBgroundEffVsNVertex->Write("RhoCorrectedBgroundEffVsNVertex");
    TGraphAsymmErrors *NVertexCorrectedBgroundEffVsNVertex = new TGraphAsymmErrors(h_NVertexCorrectedNumerVsNVertex,h_DenomVsNVertex,"");
    NVertexCorrectedBgroundEffVsNVertex->SetMarkerSize(0.6);
    NVertexCorrectedBgroundEffVsNVertex->GetXaxis()->SetTitle("NVertex");
    NVertexCorrectedBgroundEffVsNVertex->GetYaxis()->SetTitle("Efficiency");
    NVertexCorrectedBgroundEffVsNVertex->SetLineColor(kBlue);
    NVertexCorrectedBgroundEffVsNVertex->Write("NVertexCorrectedBgroundEffVsNVertex");
    TGraphAsymmErrors *BgroundEffVsEt = new TGraphAsymmErrors(h_NumerVsEt,h_DenomVsEt,"");
    BgroundEffVsEt->SetMarkerSize(0.6);
    BgroundEffVsEt->GetXaxis()->SetTitle("Et");
    BgroundEffVsEt->GetYaxis()->SetTitle("Efficiency");
    BgroundEffVsEt->Write("BgroundEffVsEt");
    TGraphAsymmErrors *RhoCorrectedBgroundEffVsEt = new TGraphAsymmErrors(h_RhoCorrectedNumerVsEt,h_DenomVsEt,"");
    RhoCorrectedBgroundEffVsEt->SetMarkerSize(0.6);
    RhoCorrectedBgroundEffVsEt->GetXaxis()->SetTitle("Et");
    //RhoCorrectedBgroundEffVsEt->GetYaxis()->SetTitleOffset(.8);
    RhoCorrectedBgroundEffVsEt->GetYaxis()->SetTitle("Efficiency");
    RhoCorrectedBgroundEffVsEt->SetLineColor(kRed);
    RhoCorrectedBgroundEffVsEt->Write("RhoCorrectedBgroundEffVsEt");
    TGraphAsymmErrors *NVertexCorrectedBgroundEffVsEt = new TGraphAsymmErrors(h_NVertexCorrectedNumerVsEt,h_DenomVsEt,"");
    NVertexCorrectedBgroundEffVsEt->SetMarkerSize(0.6);
    NVertexCorrectedBgroundEffVsEt->GetXaxis()->SetTitle("Et");
    NVertexCorrectedBgroundEffVsEt->GetYaxis()->SetTitle("Efficiency");
    NVertexCorrectedBgroundEffVsEt->SetLineColor(kBlue);
    NVertexCorrectedBgroundEffVsEt->Write("NVertexCorrectedBgroundEffVsEt");
    TGraphAsymmErrors *BgroundEffVsEta = new TGraphAsymmErrors(h_NumerVsEta,h_DenomVsEta,"");
    BgroundEffVsEta->SetMarkerSize(0.6);
    BgroundEffVsEta->GetXaxis()->SetTitle("Eta");
    BgroundEffVsEta->GetYaxis()->SetTitle("Efficiency");
    BgroundEffVsEta->Write("BgroundEffVsEta");
    TGraphAsymmErrors *RhoCorrectedBgroundEffVsEta = new TGraphAsymmErrors(h_RhoCorrectedNumerVsEta,h_DenomVsEta,"");
    RhoCorrectedBgroundEffVsEta->SetMarkerSize(0.6);
    RhoCorrectedBgroundEffVsEta->GetXaxis()->SetTitle("Eta");
    RhoCorrectedBgroundEffVsEta->GetYaxis()->SetTitle("Efficiency");
    RhoCorrectedBgroundEffVsEta->SetLineColor(kRed);
    RhoCorrectedBgroundEffVsEta->Write("RhoCorrectedBgroundEffVsEta");
    TGraphAsymmErrors *NVertexCorrectedBgroundEffVsEta = new TGraphAsymmErrors(h_NVertexCorrectedNumerVsEta,h_DenomVsEta,"");
    NVertexCorrectedBgroundEffVsEta->SetMarkerSize(0.6);
    NVertexCorrectedBgroundEffVsEta->GetXaxis()->SetTitle("Eta");
    NVertexCorrectedBgroundEffVsEta->GetYaxis()->SetTitle("Efficiency");
    NVertexCorrectedBgroundEffVsEta->SetLineColor(kBlue);
    NVertexCorrectedBgroundEffVsEta->Write("NVertexCorrectedBgroundEffVsEta");
    TGraphAsymmErrors *BgroundEffVsMET = new TGraphAsymmErrors(h_NumerVsMET,h_DenomVsMET,"");
    BgroundEffVsMET->SetMarkerSize(0.6);
    BgroundEffVsMET->GetXaxis()->SetTitle("MET");
    BgroundEffVsMET->GetYaxis()->SetTitle("Efficiency");
    BgroundEffVsMET->Write("BgroundEffVsMET");
    TGraphAsymmErrors *RhoCorrectedBgroundEffVsMET = new TGraphAsymmErrors(h_RhoCorrectedNumerVsMET,h_DenomVsMET,"");
    RhoCorrectedBgroundEffVsMET->SetMarkerSize(0.6);
    RhoCorrectedBgroundEffVsMET->GetXaxis()->SetTitle("MET");
    RhoCorrectedBgroundEffVsMET->GetYaxis()->SetTitle("Efficiency");
    RhoCorrectedBgroundEffVsMET->SetLineColor(kRed);
    RhoCorrectedBgroundEffVsMET->Write("RhoCorrectedBgroundEffVsMET");
    TGraphAsymmErrors *NVertexCorrectedBgroundEffVsMET = new TGraphAsymmErrors(h_NVertexCorrectedNumerVsMET,h_DenomVsMET,"");
    NVertexCorrectedBgroundEffVsMET->SetMarkerSize(0.6);
    NVertexCorrectedBgroundEffVsMET->GetXaxis()->SetTitle("MET");
    NVertexCorrectedBgroundEffVsMET->GetYaxis()->SetTitle("Efficiency");
    NVertexCorrectedBgroundEffVsMET->SetLineColor(kBlue);
    NVertexCorrectedBgroundEffVsMET->Write("NVertexCorrectedBgroundEffVsMET");
  }
  //output event number to txt file
  ofstream EventWriteOut;
  if(outputEventNumbers){
    EventWriteOut.open("TextOutput/RhoCorrection/Events_"+ds+".txt",ios::app);
    EventWriteOut << ggRunEvent.size() << " gg Events\n";
    EventWriteOut << ffRunEvent.size() << " ff Events\n";
    EventWriteOut << egRunEvent.size() << " eg Events\n";
    EventWriteOut << eeRunEvent.size() << " ee Events\n\n";

    EventWriteOut << ggJetReqRunEvent.size() << " gg Events with >=1 Jet Requirement\n";
    EventWriteOut << ffJetReqRunEvent.size() << " ff Events with >=1 Jet Requirement\n";
    EventWriteOut << egJetReqRunEvent.size() << " eg Events with >=1 Jet Requirement\n";
    EventWriteOut << eeJetReqRunEvent.size() << " ee Events with >=1 Jet Requirement\n\n\n";

    for(unsigned int i=0;i<ggRunEvent.size();i++){
      if(i==0)EventWriteOut << ggRunEvent.size() << " gg Events (Run   Event): \n\n";
      EventWriteOut << ggRunEvent[i].first << "  " << ggRunEvent[i].second<<"\n";
    }
    EventWriteOut << "\n\n";
    for(unsigned int i=0;i<ffRunEvent.size();i++){
      if(i==0)EventWriteOut << ffRunEvent.size() << " ff Events (Run   Event): \n\n";
      EventWriteOut << ffRunEvent[i].first << "  " << ffRunEvent[i].second<<"\n";
    }
    EventWriteOut << "\n\n";
    for(unsigned int i=0;i<egRunEvent.size();i++){
      if(i==0)EventWriteOut << egRunEvent.size() << " eg Events (Run   Event): \n\n";
      EventWriteOut << egRunEvent[i].first << "  " << egRunEvent[i].second<<"\n";
    }
    EventWriteOut << "\n\n";
    for(unsigned int i=0;i<eeRunEvent.size();i++){
      if(i==0)EventWriteOut << eeRunEvent.size() << " ee Events (Run   Event): \n\n";
      EventWriteOut << eeRunEvent[i].first << "  " << eeRunEvent[i].second<<"\n";
    }
    EventWriteOut << "\n\n";
    for(unsigned int i=0;i<ggJetReqRunEvent.size();i++){
      if(i==0)EventWriteOut << ggJetReqRunEvent.size() << " gg Events with >=1 Jet Requirement (Run   Event): \n\n";
      EventWriteOut << ggJetReqRunEvent[i].first << "  " << ggJetReqRunEvent[i].second<<"\n";
    }
    EventWriteOut << "\n\n";
    for(unsigned int i=0;i<ffJetReqRunEvent.size();i++){
      if(i==0)EventWriteOut << ffJetReqRunEvent.size() << " ff Events with >=1 Jet Requirement (Run   Event): \n\n";
      EventWriteOut << ffJetReqRunEvent[i].first << "  " << ffJetReqRunEvent[i].second<<"\n";
    }
    EventWriteOut << "\n\n";
    for(unsigned int i=0;i<egJetReqRunEvent.size();i++){
      if(i==0)EventWriteOut << egJetReqRunEvent.size() << " eg Events with >=1 Jet Requirement (Run   Event): \n\n";
      EventWriteOut << egJetReqRunEvent[i].first << "  " << egJetReqRunEvent[i].second<<"\n";
    }
    EventWriteOut << "\n\n";
    for(unsigned int i=0;i<eeJetReqRunEvent.size();i++){
      if(i==0)EventWriteOut << eeJetReqRunEvent.size() << " ee Events with >=1 Jet Requirement (Run   Event): \n\n";
      EventWriteOut << eeJetReqRunEvent[i].first << "  " << eeJetReqRunEvent[i].second<<"\n";
    }
    EventWriteOut << "\n\n";
    for(unsigned int i=0;i<failMetFilterRunEvent.size();i++){
      if(i==0)EventWriteOut << failMetFilterRunEvent.size() << " Fail HcalNoiseFilter (Run   Event): \n\n";
      EventWriteOut << failMetFilterRunEvent[i].first << "  " << failMetFilterRunEvent[i].second<<"\n";
    }
    EventWriteOut << "\n\n";
    EventWriteOut.close();
  }

  // end of event loop and print summary

  cout << " ----------------- Job Summary ----------------- " << endl;
  cout << " Total events            : " << nCnt[0] << " (" << 100*nCnt[0]/float(nCnt[0]) << "%)" << endl;
  cout << " JSON passed             : " << nCnt[1] << " (" << 100*nCnt[1]/float(nCnt[0]) << "%)" << endl;
  cout << " Duplicate check passed  : " << nCnt[2] << " (" << 100*nCnt[2]/float(nCnt[0]) << "%)" << endl;
  cout << " HLT passed              : " << nCnt[3] << " (" << 100*nCnt[3]/float(nCnt[0]) << "%)" << endl;
  cout << " nVertex passed              : " << nCnt[4] << " (" << 100*nCnt[4]/float(nCnt[0]) << "%)" << endl;
  cout << " All MET Filters passed     : " << nCnt[5] << " (" << 100*nCnt[5]/float(nCnt[0]) << "%)" << endl;
  cout << " Number Jets that fail ele cleaning in gg,eg,ee events  : " << nCnt[6] << " (" << 100*nCnt[6]/float(nCnt[0]) << "%)" << endl;
  cout << " Number Jets that fail muon cleaning in gg,eg,ee events : " << nCnt[7] << " (" << 100*nCnt[7]/float(nCnt[0]) << "%)" << endl;
  cout << " Number Jets that fail ele cleaning in ff events        : " << nCnt[11] << " (" << 100*nCnt[11]/float(nCnt[0]) << "%)" << endl;
  cout << " Number Jets that fail muon cleaning in ff events       : " << nCnt[12] << " (" << 100*nCnt[12]/float(nCnt[0]) << "%)" << endl;
  cout << " met > 20 GeV            : " << nCnt[14] << " (" << 100*nCnt[14]/float(nCnt[0]) << "%)" << endl;
  cout << " Number diPho Candidate Events     : " << ndiPhoCands  << " (" <<  100*ndiPhoCands/float(nCnt[0]) << "%)" << endl;
  cout << " Number diPho Events     : " << ndiPho  << " (" <<  100*ndiPho/float(nCnt[0]) << "%)" << endl;
  cout << " Number gg               : " << ngg     << " (" <<     100*ngg/float(nCnt[0]) << "%)" << endl;
  cout << " Number eg               : " << neg     << " (" <<     100*neg/float(nCnt[0]) << "%)" << endl;
  cout << " Number ee               : " << nee     << " (" <<     100*nee/float(nCnt[0]) << "%)" << endl;
  cout << " Number ff               : " << nff     << " (" <<     100*nff/float(nCnt[0]) << "%)" << endl;
  cout << " Number gg that have a L1FastL2L3 corrected pfJet with Pt>50 and eta<2.6 and fails loose jet id  :" << nCnt[8] << " (" << 100*nCnt[8]/float(ngg) << "% of gg)" << endl;
  cout << " Number eg that have a L1FastL2L3 corrected pfJet with Pt>50 and eta<2.6 and fails loose jet id  :" << nCnt[10] << " (" << 100*nCnt[10]/float(neg) << "% of eg)" << endl;
  cout << " Number ee that have a L1FastL2L3 corrected pfJet with Pt>50 and eta<2.6 and fails loose jet id  :" << nCnt[9] << " (" << 100*nCnt[9]/float(nee) << "% of ee)" << endl;
  cout << " Number ff that have a L1FastL2L3 corrected pfJet with Pt>50 and eta<2.6 and fails loose jet id  :" << nCnt[13] << " (" << 100*nCnt[13]/float(nff) << "% of ff)" << endl;
  cout << " Number pho_Cands that fail dR  : " << nPhosFailDR    << " (" << 100*nPhosFailDR/float(nCnt[0]) << "%)" << endl;
  cout << " Number fake_Cands that fail dR  : " << nFakesFailDR    << " (" << 100*nFakesFailDR/float(nCnt[0]) << "%)" << endl;
  cout << " Number pho_Cands that fail dPhi : " << nPhosFailDphi << " (" << 100*nPhosFailDphi/float(nCnt[0]) << "%)" << endl;
  cout << " Number fake_Cands that fail dPhi : " << nFakesFailDphi << " (" << 100*nFakesFailDphi/float(nCnt[0]) << "%)" << endl;
  cout << " Number Events with >=2 phos and >=2 fakes: "<<nTwoPhosAndTwoFakes<< " (" << 100*nTwoPhosAndTwoFakes/float(nCnt[0])<<"%)"<< endl;   
  cout << " Number diPho Candidate Events With One Jet Requirement : " << ndiPhoCands_JetReq  << " (" <<  100*ndiPhoCands_JetReq/float(nCnt[0]) << "%)" << endl;
  cout << " Number diPho Events With One Jet Requirement : " << ndiPho_JetReq  << " (" <<  100*ndiPho_JetReq/float(nCnt[0]) << "%)" << endl;
  cout << " Number gg With One Jet Requirement : " << ngg_JetReq     << " (" <<     100*ngg_JetReq/float(nCnt[0]) << "%)" << endl;
  cout << " Number eg With One Jet Requirement : " << neg_JetReq     << " (" <<     100*neg_JetReq/float(nCnt[0]) << "%)" << endl;
  cout << " Number ee With One Jet Requirement : " << nee_JetReq     << " (" <<     100*nee_JetReq/float(nCnt[0]) << "%)" << endl;
  cout << " Number ff With One Jet Requirement : " << nff_JetReq     << " (" <<     100*nff_JetReq/float(nCnt[0]) << "%)" << endl;
  cout << " Number ee With Two Jet Requirement : " << nee_2JetReq    << " (" <<     100*nee_2JetReq/float(nCnt[0]) << "%)" << endl;
  cout << " Number ee in Sideband : "	              << neeSideBand << " (" <<     100*neeSideBand/float(nCnt[0]) << "%)" << endl;
  cout << " Number ee in Sideband With One Jet Requirement : " << neeSideBand_JetReq<< " (" <<     100*neeSideBand_JetReq/float(nCnt[0]) << "%)" << endl;
  cout << " Number gge               : " << ngge     << " (" <<     100*ngge/float(nCnt[0]) << "%)" << endl;
  cout << " Number ggee              : " << nggee     << " (" <<     100*nggee/float(nCnt[0]) << "%)" << endl;
  cout << " Number eef               : " << neef     << " (" <<     100*neef/float(nCnt[0]) << "%)" << endl;
  cout << " Number eeg               : " << neeg     << " (" <<     100*neeg/float(nCnt[0]) << "%)" << endl;
  cout << " Number eegg              : " << neegg    << " (" <<     100*neegg/float(nCnt[0]) << "%)" << endl;
  cout << " Number 3Pho Events  : " << nThreePhoEvents << " (" <<     100*nThreePhoEvents/float(nCnt[0]) << "%)" << endl;
  cout << " Number 4Pho Events  : " << nFourPhoEvents << " (" <<     100*nFourPhoEvents/float(nCnt[0]) << "%)" << endl;
  cout << " Number 3Fake Events : " << nThreeFakeEvents << " (" <<     100*nThreeFakeEvents/float(nCnt[0]) << "%)" << endl;
  cout << " Number 4Fake Events : " << nFourFakeEvents << " (" <<     100*nFourFakeEvents/float(nCnt[0]) << "%)" << endl;
  cout << " Number pf electrons, no quality cuts : " << pfEleCount << endl;
  cout << " Number Mva electrons, no quality cuts : " << MvaEleCount << endl;
  cout << " Number pf electrons over gsf electrons : " << pfEleCount << " / " << eleCount << " = " << 100*pfEleCount/float(eleCount) << "%" << endl;
  cout << " Number phos pass old PhoCut    : " << oldP << endl;
  cout << " Number phos pass new PhoCut    : " << newP << endl;
  cout << " Number ifPF() photons pass : " << pfPho << endl;
  cout << " Number ifPF() photons fail : " << notPfPho << endl;


  if(enableFilter) {
    cout << " --------------- Filtered events --------------- " << endl;
    cout << " filtered events         : " << nFiltered << " (" << nFiltered/float(nCnt[0]) << ")" << endl;
  }
  cout << " ----------------------------------------------- " << endl;

  //cout<<"Writing root output to: /tmp/dmorse/hist_"<<ds<<".root"<<endl;
  cout<<"Writing root output to: /data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_"<<ds<<".root"<<endl;
  // cout<<"Writing root output to: /Users/dmorse/RA3/AnalysisOutput/hist_"<<ds<<".root"<<endl;

  // close the output file

  fout->cd();
  fout->Write();
  fout->Close();
  delete fout;

  if(enableFilter) {
    filterTree->GetCurrentFile()->cd();
    filterTree->GetCurrentFile()->Write();
    filterTree->GetCurrentFile()->Close();
    //delete filterTree;
  }

  
}


//---------------------------------------------
//       Filter!!!!!!
//---------------------------------------------


void SusyEventAnalyzer::Filter() {
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  std::cout << "total events in files  : " << nentries << std::endl;

  if(processNEvents <= 0 || processNEvents > nentries) processNEvents = nentries;

  std::cout << "events to be processed : " << processNEvents << std::endl; 


  int nFiltered = 0;
  TTree* filterTree = 0;

  if(enableFilter) {    
    TFile* filterFile = new TFile("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms525v3_jec2012/"+ds+".root","RECREATE");
    //TFile* filterFile = new TFile("/tmp/dmorse/"+ds+".root","RECREATE");
    //TFile* filterFile = new TFile("/Users/dmorse/RA3/SusyNtuples/FilteredNtuples/"+ds+".root","RECREATE");
    filterTree = (TTree*) fChain->GetTree()->CloneTree(0);
    filterTree->SetAutoSave();
  }

  int nCnt[20]={0};

  // open hist file and define histograms

  // to check duplicated events
  std::map<int, std::set<int> > allEvents;

  // start event looping

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < processNEvents; jentry++) {

    nCnt[0]++;//num total events
    if(printLevel > 0) std::cout << "Get the tree contents." << std::endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    if(printLevel > 0 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0)) ) {
      std::cout << int(jentry) << " events processed with run="
		<< event->runNumber << ", event=" << event->eventNumber << std::endl;
    }

    //if( !(event->runNumber==195112 || event->runNumber==195378) )continue;

    if(printLevel > 0) std::cout << "Initialize any global variables to be reset per event." << std::endl;

    InitializePerEvent();
    
    if(printLevel > 0) std::cout << "Check duplicated events for data only." << std::endl;
    bool duplicateEvent = ! (allEvents[event->runNumber].insert(event->eventNumber)).second;
    if(event->isRealData && duplicateEvent){
      cout<<"Duplicate Event! Run "<<event->runNumber<<" Event "<<event->eventNumber<<endl;
      continue;
    }
    nCnt[1]++;//num events pass not duplicate

    if(printLevel > 0) std::cout << "Apply good run list." << std::endl;
    // uncomment this to use the Json file to flag good data (or bad depending on your outlook)  
    if(useJSON){
      if(!isInJson(event->runNumber,event->luminosityBlockNumber)) continue;
    }
    nCnt[2]++;//num events pass json

    if(printLevel > 0) std::cout << "Apply trigger selection in the event." << std::endl;
    bool passHLT = (useTrigger ? PassTriggers() : true);
    if(!passHLT) continue;    
    nCnt[3]++;//num events pass hlt

    //  Get NVertex for event
    int NVertex=0;
    for(std::vector<susy::Vertex>::iterator Vtx_it = event->vertices.begin(); Vtx_it<event->vertices.end(); Vtx_it++){
      if(    !Vtx_it->isFake() 
	     && Vtx_it->ndof>4 
	     && fabs(Vtx_it->position.z()<24.0) 
	     && sqrt(Vtx_it->position.x()*Vtx_it->position.x()+Vtx_it->position.y()*Vtx_it->position.y())<2.0 ) NVertex++;
    }
    if(NVertex<1)continue;
    nCnt[4]++;//num events pass nvertex

    int filterCheck=0;
    int GreaterThan40=0;
    bool StudyCut=false,AnalysisCut=false;
    int AnalysisPhos=0;
    std::map<TString, std::vector<susy::Photon> >::iterator phoMap = event->photons.find("photons");
    if(phoMap != event->photons.end()) {
      for(std::vector<susy::Photon>::iterator it_Pho = phoMap->second.begin(); it_Pho != phoMap->second.end(); it_Pho++) {
	if(printLevel > 1) cout<<"looping over photon collection"<<endl;

	float ecalIsoDR03=it_Pho->ecalRecHitSumEtConeDR03;
	float hcalIsoDR03=it_Pho->hcalTowerSumEtConeDR03();
	float trackIsoDR03=it_Pho->trkSumPtHollowConeDR03;
	float PhoEt=it_Pho->momentum.Et();
	bool Et25Cut=PhoEt>25.0;	
	bool Et40Cut=PhoEt>40.0;
	bool r9Cut=it_Pho->r9<1.0;
	bool HoverECut=it_Pho->hadronicOverEm<0.05;
	bool EtaCut=std::fabs(it_Pho->caloPosition.Eta())<1.4442;
	bool CaloIdLCut=(it_Pho->hadronicOverEm<0.15 && it_Pho->sigmaIetaIeta<0.014);
	bool NVertexCut=NVertex>0;
	bool r9idCut=it_Pho->r9>0.8;
	bool IsoVLCut=( ecalIsoDR03  < 6.0 + 0.012*PhoEt ) 
	  && ( hcalIsoDR03 < 4.0 + 0.005*PhoEt ) 
	  && ( trackIsoDR03  < 4.0 + 0.002*PhoEt );
	bool R9Id85=it_Pho->r9>0.85;
	bool CaloId10=( it_Pho->hadronicOverEm<0.1 && it_Pho->sigmaIetaIeta < 0.014 );
	bool Iso50=(ecalIsoDR03<5.0 + 0.012*PhoEt
		    && hcalIsoDR03<5.0 + 0.005*PhoEt
		    && trackIsoDR03<5.0+ 0.002*PhoEt
		    );
	bool R9Id85orCaloId10andIso50Cut = R9Id85 || ( CaloId10 && Iso50 );
	
	//bool AnalysisCut = (Et25Cut && r9Cut && EtaCut && HoverECut && NVertexCut && CaloIdLCut && (r9idCut || IsoVLCut));

	//if(StudyCut)filterCheck++;
	//if(AnalysisCut)filterCheck++;
	//if(Et25Cut && r9Cut && EtaCut && CaloIdL)filterCheck++;
	//if(Et25Cut && r9Cut && HoverECut && EtaCut && CaloIdL)filterCheck++;
	if(printLevel>2)cout<<"#phos: "<<phoMap->second.size()<<" r9: "<<it_Pho->r9<<" Et: "<<PhoEt<<" h/e: "<<it_Pho->hadronicOverEm<<" filterCheck: "<<filterCheck<<std::endl;

	//StudyCut = ((Et25Cut && r9Cut && EtaCut && NVertexCut && CaloIdLCut && (r9idCut || IsoVLCut) ) || (Et25Cut && r9Cut && EtaCut && NVertexCut && R9Id85orCaloId10andIso50Cut));
	StudyCut = Et25Cut && EtaCut && NVertexCut;
	bool AnalysisPho = (Et25Cut && EtaCut && HoverECut && NVertexCut);
	if(AnalysisPho)AnalysisPhos++;
	if(Et40Cut)GreaterThan40++;
	AnalysisCut = (AnalysisPhos>1 && GreaterThan40>0);

	//if(StudyCut)break; //study
	if(AnalysisCut)break;  //analysis
      }
    }
  


    // filter conditions

    if(enableFilter) {
      //bool filterThis = StudyCut;  //study
      bool filterThis = AnalysisCut;  //analysis
      if(filterThis) {
	nFiltered++;
	filterTree->Fill();
      }
    }// if(enableFilter)

  } // for jentry


    // end of event loop and print summary

  std::cout << " ----------------- ALL DONE! ----------------- " << std::endl;

  //if(enableFilter) {
  //std::cout << " all events      : " << processNEvents << std::endl;
  cout << " Total events            : " << nCnt[0] << " (" << 100*nCnt[0]/float(nCnt[0]) << "%)" << endl;
  cout << " Duplicate check passed  : " << nCnt[1] << " (" << 100*nCnt[1]/float(nCnt[0]) << "%)" << endl;
  cout << " JSON passed             : " << nCnt[2] << " (" << 100*nCnt[2]/float(nCnt[0]) << "%)" << endl;
  cout << " HLT passed              : " << nCnt[3] << " (" << 100*nCnt[3]/float(nCnt[0]) << "%)" << endl;
  cout << " NVertex passed          : " << nCnt[4] << " (" << 100*nCnt[4]/float(nCnt[0]) << "%)" << endl;
  std::cout << " --------------- Filtered events --------------- " << std::endl;
  std::cout << " filtered events : " << nFiltered << "   ("<<nFiltered/float(nCnt[0])<<"%)"<<std::endl;
  //}
  std::cout << " ----------------------------------------------- " << std::endl;
  std::cout << " Filtered file written to: /data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms525v3_jec2012/"<<ds<<".root"<<std::endl;

  if(enableFilter) {
    filterTree->GetCurrentFile()->cd();
    filterTree->GetCurrentFile()->Write();
    filterTree->GetCurrentFile()->Close();
  }

}

//---------------------------------------------
//       DR03 Study!!!!!!
//---------------------------------------------


void SusyEventAnalyzer::DR03() {

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  cout << "total events in files  : " << nentries << endl;

  if(processNEvents <= 0 || processNEvents > nentries) processNEvents = nentries;

  cout << "events to be processed : " << processNEvents << endl; 
  


  if(printLevel > 1) cout << "Initialize event counters." << endl;
  const int NCNT = 20;
  int nCnt[NCNT];
  for(int i=0; i<NCNT; i++) nCnt[i] = 0;

  // open hist file and define histograms

  TFile* fout = new TFile("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_"+ds+".root","RECREATE");
  //TFile* fout = new TFile("/Users/dmorse/RA3/AnalysisOutput/hist_"+ds+".root","RECREATE");

  fout->cd();

  TH1F *h_SignalCombIsoDR03Nminus3 = new TH1F("SignalCombIsoDR03Nminus3","Combined Iso after Et>30,R9,h/e,sIetaIeta,Pixel=0,MET<30;combined Isolation (DR03 cone);Events",210,-5.,100.);
  TH1F *h_SignalCombRelPFIsoDR03Nminus3 = new TH1F("SignalCombRelPFIsoDR03Nminus3","Combined relative PfIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,MET<30;PF combined Isolation (DR03 cone) / E_{T};Events",210,-5.,100.);
  TH1F *h_BgroundCombRelPFIsoDR03Nminus3 = new TH1F("BgroundCombRelPFIsoDR03Nminus3","Combined relative PfIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,MET<30;PF combined Isolation (DR03 cone) / E_{T};Events",210,-5.,100.);

  TH1F *h_SignalCombIsoDR04Nminus3 = new TH1F("SignalCombIsoDR04Nminus3","Combined Iso after Et>30,R9,h/e,sIetaIeta,Pixel=0,MET<30;combined Isolation (DR04 cone);Events",210,-5.,100.);

  TH1F *h_BgroundCombIsoDR03Nminus3 = new TH1F("BgroundCombIsoDR03Nminus3","Combined Iso after Et>30,R9,h/e,sIetaIeta,Pixel=0,MET<30;combined Isolation (DR03 cone);Events",210,-5.,100.);
  TH1F *h_BgroundCombIsoDR04Nminus3 = new TH1F("BgroundCombIsoDR04Nminus3","Combined Iso after Et>30,R9,h/e,sIetaIeta,Pixel=0,MET<30;combined Isolation (DR04 cone);Events",210,-5.,100.);
  TH1F *h_BgroundEcalIsoDR03Nminus1 = new TH1F("BgroundEcalIsoDR03Nminus1","EcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,hcalIso,trackIso,MET<30;ecalRecHitSumEtConeDR03;Events",150,-5,25);
  TH1F *h_BgroundHcalIsoDR03Nminus1 = new TH1F("BgroundHcalIsoDR03Nminus1","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,trackIso,MET<30;hcalTowerSumEtConeDR03();Events",150,-5,25);
  TH1F *h_BgroundTrackIsoDR03Nminus1 = new TH1F("BgroundTrackIsoDR03Nminus1","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,hcalIso,MET<30;trkSumPtHollowConeDR03;Events",150,-5,25);
  
  TH1F* h_EcalIsoDR03 = new TH1F("EcalIsoDR03","Ecal Isolation in DR03 cone for all photons after JSON, HLT, R9, Et>30, h/e",150,-5,25);
  TH1F* h_HcalIsoDR03 = new TH1F("HcalIsoDR03","Hcal Isolation in DR03 cone for all photons after JSON, HLT, R9, Et>30, h/e",150,-5,25);
  TH1F* h_TrackIsoDR03 = new TH1F("TrackIsoDR03","Track Isolation in DR03 cone for all photons after JSON, HLT, R9, Et>30, h/e",150,-5,25);
  TH1F* h_EcalIsoDR04 = new TH1F("EcalIsoDR04","Ecal Isolation in DR04 cone for all photons after JSON, HLT, R9, Et>30, h/e",150,-5,25);
  TH1F* h_HcalIsoDR04 = new TH1F("HcalIsoDR04","Hcal Isolation in DR04 cone for all photons after JSON, HLT, R9, Et>30, h/e",150,-5,25);
  TH1F* h_TrackIsoDR04 = new TH1F("TrackIsoDR04","Track Isolation in DR04 cone for all photons after JSON, HLT, R9, Et>30, h/e",150,-5,25);
  TH1F *h_EcalIsoDR03Nminus3 = new TH1F("EcalIsoDR03Nminus3","EcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0;ecalRecHitSumEtConeDR03;Events",150,-5,25);
  TH1F *h_EcalIsoDR03Nminus1 = new TH1F("EcalIsoDR03Nminus1","EcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,hcalIso,trackIso;ecalRecHitSumEtConeDR03;Events",150,-5,25);
  TH1F *h_EcalIsoDR04Nminus1 = new TH1F("EcalIsoDR04Nminus1","EcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,hcalIso,trackIso;ecalRecHitSumEtConeDR04;Events",150,-5,25);
  TH1F *h_HcalIsoDR03Nminus3 = new TH1F("HcalIsoDR03Nminus3","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0;hcalTowerSumEtConeDR03();Events",150,-5,25);
  TH1F *h_HcalIsoDR03Nminus1 = new TH1F("HcalIsoDR03Nminus1","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,trackIso;hcalTowerSumEtConeDR03();Events",150,-5,25);
  TH1F *h_HcalIsoDR04Nminus1 = new TH1F("HcalIsoDR04Nminus1","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,trackIso;hcalTowerSumEtConeDR04();Events",150,-5,25);
  TH1F *h_TrackIsoDR03Nminus3 = new TH1F("TrackIsoDR03Nminus3","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0;trkSumPtHollowConeDR03;Events",150,-5,25);
  TH1F *h_TrackIsoDR03Nminus1 = new TH1F("TrackIsoDR03Nminus1","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,hcalIso;trkSumPtHollowConeDR03;Events",150,-5,25);
  TH1F *h_TrackIsoDR04Nminus1 = new TH1F("TrackIsoDR04Nminus1","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,hcalIso;trkSumPtHollowConeDR04;Events",150,-5,25);
  TH1F *h_EcalIsoDR03Nminus3Signal = new TH1F("EcalIsoDR03Nminus3Signal","EcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0;ecalRecHitSumEtConeDR03;Events",150,-5,25);
  TH1F *h_EcalIsoDR03Nminus1Signal = new TH1F("EcalIsoDR03Nminus1Signal","EcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,hcalIso,trackIso;ecalRecHitSumEtConeDR03;Events",150,-5,25);
  TH1F *h_EcalIsoDR04Nminus1Signal = new TH1F("EcalIsoDR04Nminus1Signal","EcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,hcalIso,trackIso;ecalRecHitSumEtConeDR04;Events",150,-5,25);
  TH1F *h_EcalIsoDR03Nminus1SignalRhoCorr = new TH1F("EcalIsoDR03Nminus1SignalRhoCorr","EcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,hcalIsoRhoCorr,trackIsoRhoCorr;ecalRecHitSumEtConeDR03;Events",150,-5,25);
  TH1F *h_EcalIsoDR04Nminus1SignalRhoCorr = new TH1F("EcalIsoDR04Nminus1SignalRhoCorr","EcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,hcalIsoRhoCorr,trackIsoRhoCorr;ecalRecHitSumEtConeDR04;Events",150,-5,25);
  TH1F *h_HcalIsoDR03Nminus3Signal = new TH1F("HcalIsoDR03Nminus3Signal","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0;hcalTowerSumEtConeDR03();Events",150,-5,25);
  TH1F *h_HcalIsoDR03Nminus1Signal = new TH1F("HcalIsoDR03Nminus1Signal","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,trackIso;hcalTowerSumEtConeDR03();Events",150,-5,25);
  TH1F *h_HcalIsoDR04Nminus1Signal = new TH1F("HcalIsoDR04Nminus1Signal","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,trackIso;hcalTowerSumEtConeDR04();Events",150,-5,25);
  TH1F *h_HcalIsoDR03Nminus1SignalRhoCorr = new TH1F("HcalIsoDR03Nminus1SignalRhoCorr","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIsoRhoCorr,trackIsoRhoCorr;hcalTowerSumEtConeDR03();Events",150,-5,25);
  TH1F *h_HcalIsoDR04Nminus1SignalRhoCorr = new TH1F("HcalIsoDR04Nminus1SignalRhoCorr","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIsoRhoCorr,trackIsoRhoCorr;hcalTowerSumEtConeDR04();Events",150,-5,25);
  TH1F *h_TrackIsoDR03Nminus3Signal = new TH1F("TrackIsoDR03Nminus3Signal","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0;trkSumPtHollowConeDR03;Events",150,-5,25);
  TH1F *h_TrackIsoDR03Nminus1Signal = new TH1F("TrackIsoDR03Nminus1Signal","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,hcalIso;trkSumPtHollowConeDR03;Events",150,-5,25);
  TH1F *h_TrackIsoDR04Nminus1Signal = new TH1F("TrackIsoDR04Nminus1Signal","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,hcalIso;trkSumPtHollowConeDR04;Events",150,-5,25);
  TH1F *h_EcalIsoDR03Nminus1Bground = new TH1F("EcalIsoDR03Nminus1Bground","EcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,hcalIso,trackIso;ecalRecHitSumEtConeDR03;Events",150,-5,25);
  TH1F *h_HcalIsoDR03Nminus1Bground = new TH1F("HcalIsoDR03Nminus1Bground","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,trackIso;hcalTowerSumEtConeDR03();Events",150,-5,25);
  TH1F *h_TrackIsoDR03Nminus1Bground = new TH1F("TrackIsoDR03Nminus1Bground","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,hcalIso;trkSumPtHollowConeDR03;Events",150,-5,25);

  TH1F *h_chargedHadronIsoNminus3Bground = new TH1F("chargedHadronIsoNminus3Bground","chargedHadronIso after Et>25,R9,h/e,sIetaIeta,Pixel=0;chargedHadronIso;Events",150,-5,30);
  TH1F *h_neutralHadronIsoNminus3Bground = new TH1F("neutralHadronIsoNminus3Bground","neutralHadronIso after Et>25,R9,h/e,sIetaIeta,Pixel=0;neutralHadronIso;Events",150,-5,30);
  TH1F *h_photonIsoNminus3Bground = new TH1F("photonIsoNminus3Bground","photonIso after Et>25,R9,h/e,sIetaIeta,Pixel=0;photonIso;Events",150,-5,30);
  TH1F *h_chargedHadronIsoDepositNminus3Bground = new TH1F("chargedHadronIsoDepositNminus3Bground","chargedHadronIsoDeposit after Et>25,R9,h/e,sIetaIeta,Pixel=0;chargedHadronIsoDeposit;Events",150,-5,30);
  TH1F *h_neutralHadronIsoDepositNminus3Bground = new TH1F("neutralHadronIsoDepositNminus3Bground","neutralHadronIsoDeposit after Et>25,R9,h/e,sIetaIeta,Pixel=0;neutralHadronIsoDeposit;Events",150,-5,30);
  TH1F *h_photonIsoDepositNminus3Bground = new TH1F("photonIsoDepositNminus3Bground","PfTrackIsoDeposit after Et>25,R9,h/e,sIetaIeta,Pixel=0;photonIsoDeposit;Events",150,-5,30);
  TH1F* h_PfCombinedIsoNminus3 = new TH1F("PfCombinedIsoNminus3","Combined Pf Isolation after Et>25,R9,h/e,sIetaIeta,Pixel=0;chargedHadronIso+neutralHadronIso+photonIso;Events",210,-5.,100.);
  TH1F* h_PfCombinedIsoNminus3Bground = new TH1F("PfCombinedIsoNminus3Bground","Combined Pf Isolation after Et>25,R9,h/e,sIetaIeta,Pixel=0 in events with MET<30GeV;chargedHadronIso+neutralHadronIso+photonIso;Events",210,-5.,100.);
  TH1F* h_PfCombinedIsoDepositNminus3 = new TH1F("PfCombinedIsoDepositNminus3","Combined Pf IsolationDeposit after Et>25,R9,h/e,sIetaIeta,Pixel=0;chargedHadronIsoDeposit+neutralHadronIsoDeposit+photonIsoDeposit;Events",210,-5.,100.);
  TH1F* h_PfCombinedIsoDepositNminus3Bground = new TH1F("PfCombinedIsoDepositNminus3Bground","Combined Pf IsolationDeposit after Et>25,R9,h/e,sIetaIeta,Pixel=0 in events with MET<30GeV;chargedHadronIsoDeposit+neutralHadronIsoDeposit+photonIsoDeposit;Events",210,-5.,100.);

  TH1F *h_chargedHadronIsoNminus3Signal = new TH1F("chargedHadronIsoNminus3Signal","chargedHadronIso after Et>25,R9,h/e,sIetaIeta,Pixel=0;chargedHadronIso;Events",150,-5,30);
  TH1F *h_neutralHadronIsoNminus3Signal = new TH1F("neutralHadronIsoNminus3Signal","neutralHadronIso after Et>25,R9,h/e,sIetaIeta,Pixel=0;neutralHadronIso;Events",150,-5,30);
  TH1F *h_photonIsoNminus3Signal = new TH1F("photonIsoNminus3Signal","photonIso after Et>25,R9,h/e,sIetaIeta,Pixel=0;photonIso;Events",150,-5,30);
  TH1F *h_chargedHadronIsoDepositNminus3Signal = new TH1F("chargedHadronIsoDepositNminus3Signal","chargedHadronIsoDeposit after Et>25,R9,h/e,sIetaIeta,Pixel=0;chargedHadronIsoDeposit;Events",150,-5,30);
  TH1F *h_neutralHadronIsoDepositNminus3Signal = new TH1F("neutralHadronIsoDepositNminus3Signal","neutralHadronIsoDeposit after Et>25,R9,h/e,sIetaIeta,Pixel=0;neutralHadronIsoDeposit;Events",150,-5,30);
  TH1F *h_photonIsoDepositNminus3Signal = new TH1F("photonIsoDepositNminus3Signal","PfTrackIsoDeposit after Et>25,R9,h/e,sIetaIeta,Pixel=0;photonIsoDeposit;Events",150,-5,30);
  TH1F* h_PfCombinedIsoNminus3Signal = new TH1F("PfCombinedIsoNminus3Signal","Combined Pf Isolation after Et>25,R9,h/e,sIetaIeta,Pixel=0;chargedHadronIso+neutralHadronIso+photonIso;Events",210,-5.,100.);
  TH1F* h_PfCombinedIsoDepositNminus3Signal = new TH1F("PfCombinedIsoDepositNminus3Signal","Combined Pf IsolationDeposit after Et>25,R9,h/e,sIetaIeta,Pixel=0;chargedHadronIsoDeposit+neutralHadronIsoDeposit+photonIsoDeposit;Events",210,-5.,100.);

  TH1F *h_EcalIsoDR04Nminus1Bground = new TH1F("EcalIsoDR04Nminus1Bground","EcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,hcalIso,trackIso;ecalRecHitSumEtConeDR04;Events",150,-5,25);
  TH1F *h_HcalIsoDR04Nminus1Bground = new TH1F("HcalIsoDR04Nminus1Bground","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,trackIso;hcalTowerSumEtConeDR04();Events",150,-5,25);
  TH1F *h_TrackIsoDR04Nminus1Bground = new TH1F("TrackIsoDR04Nminus1Bground","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIso,hcalIso;trkSumPtHollowConeDR04;Events",150,-5,25);
  TH1F *h_TrackIsoDR03Nminus1SignalRhoCorr = new TH1F("TrackIsoDR03Nminus1SignalRhoCorr","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIsoRhoCorr,hcalIsoRhoCorr;trkSumPtHollowConeDR03;Events",150,-5,25);
  TH1F *h_TrackIsoDR04Nminus1SignalRhoCorr = new TH1F("TrackIsoDR04Nminus1SignalRhoCorr","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIsoRhoCorr,hcalIsoRhoCorr;trkSumPtHollowConeDR04;Events",150,-5,25);
  TH1F *h_EcalIsoDR03Nminus1BgroundRhoCorr = new TH1F("EcalIsoDR03Nminus1BgroundRhoCorr","EcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,hcalIsoRhoCorr,trackIsoRhoCorr;ecalRecHitSumEtConeDR03;Events",150,-5,25);
  TH1F *h_HcalIsoDR03Nminus1BgroundRhoCorr = new TH1F("HcalIsoDR03Nminus1BgroundRhoCorr","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIsoRhoCorr,trackIsoRhoCorr;hcalTowerSumEtConeDR03();Events",150,-5,25);
  TH1F *h_TrackIsoDR03Nminus1BgroundRhoCorr = new TH1F("TrackIsoDR03Nminus1BgroundRhoCorr","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIsoRhoCorr,hcalIsoRhoCorr;trkSumPtHollowConeDR03;Events",150,-5,25);
  TH1F *h_EcalIsoDR04Nminus1BgroundRhoCorr = new TH1F("EcalIsoDR04Nminus1BgroundRhoCorr","EcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,hcalIsoRhoCorr,trackIsoRhoCorr;ecalRecHitSumEtConeDR04;Events",150,-5,25);
  TH1F *h_HcalIsoDR04Nminus1BgroundRhoCorr = new TH1F("HcalIsoDR04Nminus1BgroundRhoCorr","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIsoRhoCorr,trackIsoRhoCorr;hcalTowerSumEtConeDR04();Events",150,-5,25);
  TH1F *h_TrackIsoDR04Nminus1BgroundRhoCorr = new TH1F("TrackIsoDR04Nminus1BgroundRhoCorr","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0,ecalIsoRhoCorr,hcalIsoRhoCorr;trkSumPtHollowConeDR04;Events",150,-5,25);
  TH1F *h_EcalIsoDR03Nminus3Bground = new TH1F("EcalIsoDR03Nminus3Bground","EcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0;ecalRecHitSumEtConeDR03;Events",150,-5,25);
  TH1F *h_HcalIsoDR03Nminus3Bground = new TH1F("HcalIsoDR03Nminus3Bground","HcalIso after Et>30,R9,h/e,sIetaIeta,Pixel=0;hcalTowerSumEtConeDR03();Events",150,-5,25);
  TH1F *h_TrackIsoDR03Nminus3Bground = new TH1F("TrackIsoDR03Nminus3Bground","TrackIso after Et>30,R9,h/e,sIetaIeta,Pixel=0;trkSumPtHollowConeDR03;Events",150,-5,25);
  TH1F *h_DenomAllButIsosVsNVertex = new TH1F("DenomAllButIsosVsNVertex","",40,0,40);
  TH1F *h_NumerEcalNormalHcalNormalVsNVertex = new TH1F("NumerEcalNormalHcalNormalVsNVertex","",40,0,40);
  TH1F *h_NumerEcalNormalHcalRhoVsNVertex = new TH1F("NumerEcalNormalHcalRhoVsNVertex","",40,0,40);
  TH1F *h_NumerEcalRhoHcalNormalVsNVertex = new TH1F("NumerEcalRhoHcalNormalVsNVertex","",40,0,40);
  TH1F *h_NumerEcalRhoHcalRhoVsNVertex = new TH1F("NumerEcalRhoHcalRhoVsNVertex","",40,0,40);
  TH1F *h_DrDenomSignal = new TH1F("DrDenomSignal","",64,0,4);h_DrDenomSignal->Sumw2();
  TH1F *h_DrDenomBground = new TH1F("DrDenomBground","",64,0,4);h_DrDenomBground->Sumw2();
  TH1F *h_DrGenDenomSignal = new TH1F("DrGenDenomSignal","",64,0,4);h_DrGenDenomSignal->Sumw2();
  TH1F *h_DrGenDenomBground = new TH1F("DrGenDenomBground","",64,0,4);h_DrGenDenomBground->Sumw2();
  TH1F *h_DrNumerSignalDR03 = new TH1F("DrNumerSignalDR03","",64,0,4);h_DrNumerSignalDR03->Sumw2();
  TH1F *h_DrNumerSignalDR04 = new TH1F("DrNumerSignalDR04","",64,0,4);h_DrNumerSignalDR04->Sumw2();
  TH1F *h_DrNumerBgroundDR03 = new TH1F("DrNumerBgroundDR03","",64,0,4);h_DrNumerBgroundDR03->Sumw2();
  TH1F *h_DrNumerBgroundDR04 = new TH1F("DrNumerBgroundDR04","",64,0,4);h_DrNumerBgroundDR04->Sumw2();
  TH1F *h_DrNumerSignalDR03RhoCorr = new TH1F("DrNumerSignalDR03RhoCorr","",64,0,4);h_DrNumerSignalDR03RhoCorr->Sumw2();
  TH1F *h_DrNumerSignalDR04RhoCorr = new TH1F("DrNumerSignalDR04RhoCorr","",64,0,4);h_DrNumerSignalDR04RhoCorr->Sumw2();
  TH1F *h_DrNumerBgroundDR03RhoCorr = new TH1F("DrNumerBgroundDR03RhoCorr","",64,0,4);h_DrNumerBgroundDR03RhoCorr->Sumw2();
  TH1F *h_DrNumerBgroundDR04RhoCorr = new TH1F("DrNumerBgroundDR04RhoCorr","",64,0,4);h_DrNumerBgroundDR04RhoCorr->Sumw2();
  TH1F *h_DrGenNumerSignalDR03 = new TH1F("DrGenNumerSignalDR03","",64,0,4);h_DrGenNumerSignalDR03->Sumw2();
  TH1F *h_DrGenNumerSignalDR04 = new TH1F("DrGenNumerSignalDR04","",64,0,4);h_DrGenNumerSignalDR04->Sumw2();
  TH1F *h_DrGenNumerBgroundDR03 = new TH1F("DrGenNumerBgroundDR03","",64,0,4);h_DrGenNumerBgroundDR03->Sumw2();
  TH1F *h_DrGenNumerBgroundDR04 = new TH1F("DrGenNumerBgroundDR04","",64,0,4);h_DrGenNumerBgroundDR04->Sumw2();
  TH1F *h_DrGenNumerSignalDR03RhoCorr = new TH1F("DrGenNumerSignalDR03RhoCorr","",64,0,4);h_DrGenNumerSignalDR03RhoCorr->Sumw2();
  TH1F *h_DrGenNumerSignalDR04RhoCorr = new TH1F("DrGenNumerSignalDR04RhoCorr","",64,0,4);h_DrGenNumerSignalDR04RhoCorr->Sumw2();
  TH1F *h_DrGenNumerBgroundDR03RhoCorr = new TH1F("DrGenNumerBgroundDR03RhoCorr","",64,0,4);h_DrGenNumerBgroundDR03RhoCorr->Sumw2();
  TH1F *h_DrGenNumerBgroundDR04RhoCorr = new TH1F("DrGenNumerBgroundDR04RhoCorr","",64,0,4);h_DrGenNumerBgroundDR04RhoCorr->Sumw2();
  TH1F *h_EtDenomSignal = new TH1F("EtDenomSignal","",100,0,500);h_EtDenomSignal->Sumw2();
  TH1F *h_EtDenomBground = new TH1F("EtDenomBground","",100,0,500);h_EtDenomBground->Sumw2();
  TH1F *h_EtGenDenomSignal = new TH1F("EtGenDenomSignal","",100,0,500);h_EtGenDenomSignal->Sumw2();
  TH1F *h_EtGenDenomBground = new TH1F("EtGenDenomBground","",100,0,500);h_EtGenDenomBground->Sumw2();
  TH1F *h_EtNumerSignalDR03 = new TH1F("EtNumerSignalDR03","",100,0,500);h_EtNumerSignalDR03->Sumw2();
  TH1F *h_EtNumerSignalDR04 = new TH1F("EtNumerSignalDR04","",100,0,500);h_EtNumerSignalDR04->Sumw2();
  TH1F *h_EtNumerBgroundDR03 = new TH1F("EtNumerBgroundDR03","",100,0,500);h_EtNumerBgroundDR03->Sumw2();
  TH1F *h_EtNumerBgroundDR04 = new TH1F("EtNumerBgroundDR04","",100,0,500);h_EtNumerBgroundDR04->Sumw2();
  TH1F *h_EtNumerSignalDR03RhoCorr = new TH1F("EtNumerSignalDR03RhoCorr","",100,0,500);h_EtNumerSignalDR03RhoCorr->Sumw2();
  TH1F *h_EtNumerSignalDR04RhoCorr = new TH1F("EtNumerSignalDR04RhoCorr","",100,0,500);h_EtNumerSignalDR04RhoCorr->Sumw2();
  TH1F *h_EtNumerBgroundDR03RhoCorr = new TH1F("EtNumerBgroundDR03RhoCorr","",100,0,500);h_EtNumerBgroundDR03RhoCorr->Sumw2();
  TH1F *h_EtNumerBgroundDR04RhoCorr = new TH1F("EtNumerBgroundDR04RhoCorr","",100,0,500);h_EtNumerBgroundDR04RhoCorr->Sumw2();
  TH1F *h_EtGenNumerSignalDR03 = new TH1F("EtGenNumerSignalDR03","",100,0,500);h_EtGenNumerSignalDR03->Sumw2();
  TH1F *h_EtGenNumerSignalDR04 = new TH1F("EtGenNumerSignalDR04","",100,0,500);h_EtGenNumerSignalDR04->Sumw2();
  TH1F *h_EtGenNumerBgroundDR03 = new TH1F("EtGenNumerBgroundDR03","",100,0,500);h_EtGenNumerBgroundDR03->Sumw2();
  TH1F *h_EtGenNumerBgroundDR04 = new TH1F("EtGenNumerBgroundDR04","",100,0,500);h_EtGenNumerBgroundDR04->Sumw2();
  TH1F *h_EtGenNumerSignalDR03RhoCorr = new TH1F("EtGenNumerSignalDR03RhoCorr","",100,0,500);h_EtGenNumerSignalDR03RhoCorr->Sumw2();
  TH1F *h_EtGenNumerSignalDR04RhoCorr = new TH1F("EtGenNumerSignalDR04RhoCorr","",100,0,500);h_EtGenNumerSignalDR04RhoCorr->Sumw2();
  TH1F *h_EtGenNumerBgroundDR03RhoCorr = new TH1F("EtGenNumerBgroundDR03RhoCorr","",100,0,500);h_EtGenNumerBgroundDR03RhoCorr->Sumw2();
  TH1F *h_EtGenNumerBgroundDR04RhoCorr = new TH1F("EtGenNumerBgroundDR04RhoCorr","",100,0,500);h_EtGenNumerBgroundDR04RhoCorr->Sumw2();
  TH2I *h_genPdgId = new TH2I("genPdgId","",300,0,300,4,0,4);
  TProfile *h_SigCombIsoDR03VsGenDR = new TProfile("SigCombIsoDR03VsGenDR","",32,0,4,"");h_SigCombIsoDR03VsGenDR->Sumw2();
  TProfile *h_SigCombIsoDR04VsGenDR = new TProfile("SigCombIsoDR04VsGenDR","",32,0,4,"");h_SigCombIsoDR04VsGenDR->Sumw2();

  TH1F *h_Pt = new TH1F("Pt","",400,0,2000);

  susy::Particle* genPart = new susy::Particle; genPart->Init();
  susy::Particle* genJet = new susy::Particle; genJet->Init();

 
  if(doRhoCorrection)cout<<"Applying Rho Pileup corrections!"<<endl;
  else cout<<"Applying NO Pileup corrections!"<<endl;
  
  // to check duplicated events
  //std::map<int, std::set<int> > allEvents;

  // start event looping
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < processNEvents; jentry++) {

    if(printLevel > 1) cout << "Get the tree contents." << endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    if(printLevel > 1 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0)) ) {
      cout << endl << int(jentry) << " events processed with Run:"
	   << event->runNumber << ", Event:" << event->eventNumber << endl;
    }


    if(printLevel > 1) cout << "Initialize any global variables to be reset per event." << endl;

    InitializePerEvent();

    if(!event->isRealData &&jentry==0){
      cout<<"Neutralino Mass: " <<event->gridParams["mChi0"]<<endl;
      cout<<"Gluino Mass: " <<event->gridParams["mGluino"]<<endl;
      cout<<"Squark Mass: " <<event->gridParams["mSquark"]<<endl;
      cout<<"Cross Section: " <<event->gridParams["xsec"]<<endl;
      //cout<<"ptHat: " <<event->gridParams["ptHat"]<<endl;
    }
    //LumiWeight = (x-sec)(IntegratedLumi)/(# events)
    float LumiWeight = 1;
    
    //LumiWeight=(4.18E+07*5000./6094670.)*2.60E-04;        //30-40 double EM Enriched QCD
    //LumiWeight=(1.87E+07*5000./40287002.)*0.00233; //40 double EM Enriched QCD
    //LumiWeight=8.16E+08*5000./11000000; //15to30 
    //LumiWeight=5.31E+07*5000./6583068; //30to50 
    //LumiWeight=6360000.0*5000./6600000; //50to80 
    //LumiWeight=784000.0*5000./6589956; //80to120 
    //LumiWeight=115000.0*5000./6127528; //120to170 
    //LumiWeight=24300.0*5000./6220160; //170to300 
    //LumiWeight=1170.0*5000./6432669; //300to470 
    //LumiWeight=70.2*5000./3990085; //470to600 
    //LumiWeight=15.6*5000./4245695; //600to800 
    //LumiWeight=1.84*5000./4053888; //800to1000 
    //LumiWeight=0.332*5000./2093222; //1000to1400 
    //LumiWeight=0.0109*5000./2196200; //1400to1800
    //LumiWeight=3.58E-04*5000./293139; //1800 
    if(printLevel>1 && jentry%10000==0)cout<<"\nLumiWeight: "<<LumiWeight<<endl;
    //cout<<"\nisRealData: "<<event->isRealData<<endl;
  
    bool isRealData=event->isRealData;

    if(printLevel > 1) cout << "Apply good run list." << endl;
    // uncomment this to use the Json file to flag good data (or bad depending on your outlook)  
    if(printLevel > 1)cout<<"runNumber="<<event->runNumber<<"  lumiNumber="<<event->luminosityBlockNumber<<endl;
    if(useJSON){
      if(!isInJson(event->runNumber,event->luminosityBlockNumber)) continue;
    }
    // uncomment this to print all ntuple variables
    //Print(*event);
    /*
      if(printLevel > 1) cout << "Check duplicated events for data only." << endl;
      if(isRealData){
      bool duplicateEvent = ! (allEvents[event->runNumber].insert(event->eventNumber)).second;
      if(duplicateEvent){
      cout<<"Duplicate Event! Run "<<event->runNumber<<" Event "<<event->eventNumber<<endl;
      continue;
      }
      }
    */
    nCnt[0]++; // total number of events    

    if(printLevel > 1) cout << "Apply trigger selection in the event." << endl;
    bool passHLT = (useTrigger ? PassTriggers() : true);
    if(!passHLT )continue;//only accept events which pass our trigger(s)
    nCnt[1]++;// number of events that pass trigger

    if(printLevel > 1) cout << "Find primary vertex in the event." << endl;
    TVector3* primVtx = 0;
    if(event->vertices.size() > 0) primVtx = &(event->vertices[0].position);

    //  Get NVertex and Rho for event
    int NVertex=0;
    for(std::vector<susy::Vertex>::iterator Vtx_it = event->vertices.begin(); Vtx_it<event->vertices.end(); Vtx_it++){
      if(    !Vtx_it->isFake() 
	     && Vtx_it->ndof>4 
	     && fabs(Vtx_it->position.z()<24.0) 
	     && sqrt(Vtx_it->position.x()*Vtx_it->position.x()+Vtx_it->position.y()*Vtx_it->position.y())<2.0 ) NVertex++;
    }
    float Rho = event->rho;

    //Require at least 1 good vertex
    if(NVertex<1){
      cout<<"No Good Vertex!!!!  Run: "<<event->runNumber<<"  Event: "<<event->eventNumber<<endl;
      continue;
    }

    if(printLevel > 1) cout << "Select which met will be used in the event." << endl;
    
    std::map<TString, susy::MET>::iterator met_it = event->metMap.find("pfMet");
    if(met_it == event->metMap.end()) {
      cout << "MET map is not available!!!" << endl;
      continue;
    }
    susy::MET* met = &(met_it->second);
    if(printLevel > 1) cout << "Find loose and tight photons in the event." << endl;
    //----------
    //find photons, sort by energy so that later we can just grab top two
    //----------
    std::map<TString, std::vector<susy::Photon> >::iterator phoMap = event->photons.find("photons");
    if(phoMap != event->photons.end()) {
      
      //   Need to sort Photons and select two with largest Pt.  
    
      //loop over photon collection 
      for(std::vector<susy::Photon>::iterator it_Pho = phoMap->second.begin(); it_Pho != phoMap->second.end(); it_Pho++) {
	if(printLevel > 1) cout<<"looping over photon collection"<<endl;
	
	
	if(printLevel > 1)cout<<"chargedHadronIso: "<<it_Pho->chargedHadronIso<<"  neutralHadronIso:"<<it_Pho->neutralHadronIso<<"  photonIso: "<<it_Pho->photonIso<<endl<<endl;
	//----------------set up cuts-------------------
	float ecalIsoDR03=it_Pho->ecalRecHitSumEtConeDR03;
	float hcalIsoDR03=it_Pho->hcalTowerSumEtConeDR03();
	float trackIsoDR03=it_Pho->trkSumPtHollowConeDR03;
	float ecalIsoDR04=it_Pho->ecalRecHitSumEtConeDR04;
	float hcalIsoDR04=it_Pho->hcalTowerSumEtConeDR04();
	float trackIsoDR04=it_Pho->trkSumPtHollowConeDR04;
	float PhoEt=it_Pho->momentum.Et();
	float chargedHadronIso=it_Pho->chargedHadronIso;
	float neutralHadronIso=it_Pho->neutralHadronIso;
	float photonIso=it_Pho->photonIso;
	float chargedHadronIsoDeposit=it_Pho->chargedHadronIsoDeposit;
	float neutralHadronIsoDeposit=it_Pho->neutralHadronIsoDeposit;
	float photonIsoDeposit=it_Pho->photonIsoDeposit;
	// fiducial cuts. Look for only barrel now
	bool etaCut = (std::abs(it_Pho->caloPosition.Eta()) < susy::etaGapBegin);
	
	// Spike cleaning
	bool isSpike = (it_Pho->r9 > 1.0);
	//if(isSpike) continue;
	
	// Et cuts, 30 GeV for trailing photons. Will apply tighter for the leading one.
	bool EtCut = (PhoEt > 30.0);
	// cuts containing EE cases for later use, but EE is not used for the time being.
	
	// H/E (in trigger, 0.15 for EB, 0.10 for EE)
	bool heCut = (it_Pho->hadronicOverEm < 0.05);
	bool heCutLoose = (it_Pho->hadronicOverEm < 0.1);
	// sigma_ietaieta (in trigger 0.014 for EB, 0.034 for EE)
	bool sIetaCut = (it_Pho->sigmaIetaIeta < 0.011);
	
	// Ecal Isolation
	bool ecalIsoCut = (ecalIsoDR04 < 4.2 + 0.006*PhoEt);
	bool ecalIsoCutdR03 = (ecalIsoDR03 < 4.2 + 0.006*PhoEt);
	bool ecalIsoCutRhoCorr = (ecalIsoDR04 - PUCorr_ECAL*Rho < 4.2 + 0.006*PhoEt);
	bool ecalIsoCutdR03RhoCorr = (ecalIsoDR03 - PUCorr_ECAL*Rho < 4.2 + 0.006*PhoEt);
	
	// Hcal Isolation
	bool hcalIsoCut = (hcalIsoDR04 < 2.2 + 0.0025*PhoEt);
	bool hcalIsoCutdR03 = (hcalIsoDR03 < 2.2 + 0.0025*PhoEt);
	bool hcalIsoCutRhoCorr  = (hcalIsoDR04 - PUCorr_HCAL*Rho < 2.2 + 0.0025*PhoEt);
	bool hcalIsoCutdR03RhoCorr  = (hcalIsoDR03 - PUCorr_HCAL*Rho < 2.2 + 0.0025*PhoEt);

	// Track Isolation
	bool trackIsoCut = (trackIsoDR04 < 2.0 + 0.001*PhoEt);
	bool trackIsoCutdR03 = (trackIsoDR03 < 2.0 + 0.001*PhoEt);
	bool trackIsoCutRhoCorr = (trackIsoDR04 - PUCorr_TRACK*Rho < 2.0 + 0.001*PhoEt);
	bool trackIsoCutdR03RhoCorr = (trackIsoDR03 - PUCorr_TRACK*Rho < 2.0 + 0.001*PhoEt);

	//Combined Isolation
	bool combIsoCutDR03 = ecalIsoDR03 + hcalIsoDR03 + trackIsoDR03 < 6.0;
	bool combIsoCutDR03RhoCorr = ecalIsoDR03 - PUCorr_ECAL*Rho + hcalIsoDR03 - PUCorr_HCAL*Rho + trackIsoDR03 < 6.0;

	//PF Isolation
	bool PfIsoCut = (chargedHadronIso+neutralHadronIso+photonIso)/PhoEt < 0.2;

	//hlt cuts	
	bool ecalIsoVLcut = ( ecalIsoDR03  < 6.0 + 0.012*PhoEt );
	bool hcalIsoVLcut = ( hcalIsoDR03 < 4.0 + 0.005*PhoEt );
	bool trackIsoVLcut = ( trackIsoDR03  < 4.0 + 0.002*PhoEt );
	bool IsoVLCut = ( ecalIsoVLcut && hcalIsoVLcut && trackIsoVLcut );
	bool CaloIdLCut = ( (it_Pho->hadronicOverEm < 0.15) && (it_Pho->sigmaIetaIeta < 0.014) );
	bool R9TriggerCut = (it_Pho->r9 > 0.8);

	//Pixel cut
	bool pixelCut = (it_Pho->nPixelSeeds == 0);
	h_Pt->Fill(it_Pho->momentum.Pt(),LumiWeight);
	//dR03 cone study
	if(!/*event->*/isRealData){
	  if( etaCut && EtCut && !isSpike && heCut && sIetaCut && pixelCut && CaloIdLCut && (IsoVLCut || R9TriggerCut) ){
	    //h_Pt->Fill(it_Pho->momentum.Pt(),LumiWeight);
	    //loop over all genParticles and grab the nearest one to the photon
	    float nearestDR=999.;
	    float SmallestDR=999.;
	    float genJetDR=999.;
	    for(std::vector<susy::Particle>::iterator Part_it = event->genParticles.begin(); 
		Part_it != event->genParticles.end(); Part_it++){
	      if(Part_it->pdgId==1000022){
		h_genPdgId->Fill(299,Part_it->status);
	      }
	      else{
		h_genPdgId->Fill(fabs(Part_it->pdgId),Part_it->status);
	      }
	      if(getDR(it_Pho->caloPosition,Part_it->momentum)<nearestDR){
		nearestDR=getDR(it_Pho->caloPosition,Part_it->momentum);
		genPart=&*Part_it;
	      }
	    }	  
	    for(std::vector<susy::Particle>::iterator genJet_it = event->genParticles.begin(); 
		genJet_it != event->genParticles.end(); genJet_it++){
	      if(fabs(genJet_it->pdgId)<7 || fabs(genJet_it->pdgId)==21){
		if(getDR(it_Pho->caloPosition,genJet_it->momentum)<genJetDR){
		  genJetDR=getDR(it_Pho->caloPosition,genJet_it->momentum);
		  genJet=&*genJet_it;
		}
	      }
	    }
	    if(nearestDR<900){
	      if(genPart->pdgId==22/* && genPart->motherId==1000022*/ ){
		if(isSameObject(it_Pho->caloPosition,genPart->momentum,0.1)){
		  if(SmallestDR<90){
		    h_DrDenomSignal->Fill(SmallestDR);
		    h_EtDenomSignal->Fill(PhoEt);
		  }
		  h_HcalIsoDR03Nminus3Signal->Fill(hcalIsoDR03,/*1*/LumiWeight);
		  h_EcalIsoDR03Nminus3Signal->Fill(ecalIsoDR03,/*1*/LumiWeight);
		  h_TrackIsoDR03Nminus3Signal->Fill(trackIsoDR03,/*1*/LumiWeight);
		  h_SignalCombIsoDR03Nminus3->Fill(ecalIsoDR03+hcalIsoDR03+trackIsoDR03,/*1*/LumiWeight);
		  h_SignalCombIsoDR04Nminus3->Fill(ecalIsoDR04+hcalIsoDR04+trackIsoDR04,/*1*/LumiWeight);
		  h_SignalCombRelPFIsoDR03Nminus3->Fill((it_Pho->chargedHadronIso+it_Pho->neutralHadronIso+it_Pho->photonIso),LumiWeight);
		  h_chargedHadronIsoNminus3Signal->Fill(chargedHadronIso);
		  h_neutralHadronIsoNminus3Signal->Fill(neutralHadronIso);
		  h_photonIsoNminus3Signal->Fill(photonIso);
		  h_chargedHadronIsoDepositNminus3Signal->Fill(chargedHadronIsoDeposit);
		  h_neutralHadronIsoDepositNminus3Signal->Fill(neutralHadronIsoDeposit);
		  h_photonIsoDepositNminus3Signal->Fill(photonIsoDeposit);
		  h_PfCombinedIsoNminus3Signal->Fill(chargedHadronIso+neutralHadronIso+photonIso);
		  h_PfCombinedIsoDepositNminus3Signal->Fill(chargedHadronIsoDeposit+neutralHadronIsoDeposit+photonIsoDeposit);
		  if(genJetDR<900){
		    h_DrGenDenomSignal->Fill(genJetDR);
		    h_EtGenDenomSignal->Fill(PhoEt);  
		  }
		  if(ecalIsoCut && trackIsoCut){
		    h_HcalIsoDR03Nminus1Signal->Fill(hcalIsoDR03);
		    h_HcalIsoDR04Nminus1Signal->Fill(hcalIsoDR04);
		  }
		  if(hcalIsoCut && trackIsoCut){
		    h_EcalIsoDR03Nminus1Signal->Fill(ecalIsoDR03);
		    h_EcalIsoDR04Nminus1Signal->Fill(ecalIsoDR04);
		  }
		  if(ecalIsoCut && hcalIsoCut){
		    h_TrackIsoDR03Nminus1Signal->Fill(trackIsoDR03);
		    h_TrackIsoDR04Nminus1Signal->Fill(trackIsoDR04);
		  } 
		  if(SmallestDR<900 && ecalIsoCut && hcalIsoCut && trackIsoCut){
		    h_DrNumerSignalDR04->Fill(SmallestDR);
		    h_EtNumerSignalDR04->Fill(PhoEt);
		  }
		  if(SmallestDR<900 && ecalIsoCutdR03 && hcalIsoCutdR03 && trackIsoCutdR03){
		    h_DrNumerSignalDR03->Fill(SmallestDR);
		    h_EtNumerSignalDR03->Fill(PhoEt);
		  }
		  if(genJetDR<900 && ecalIsoCut && hcalIsoCut && trackIsoCut){
		    h_DrGenNumerSignalDR04->Fill(genJetDR);
		    h_EtGenNumerSignalDR04->Fill(PhoEt);
		    h_SigCombIsoDR03VsGenDR->Fill(genJetDR,trackIsoDR03+ecalIsoDR03+hcalIsoDR03);
		    h_SigCombIsoDR04VsGenDR->Fill(genJetDR,trackIsoDR04+ecalIsoDR04+hcalIsoDR04);
		  }
		  if(genJetDR<900 && ecalIsoCutdR03 && hcalIsoCutdR03 && trackIsoCutdR03){
		    h_DrGenNumerSignalDR03->Fill(genJetDR);
		    h_EtGenNumerSignalDR03->Fill(PhoEt);
		  }
		  //now Rho corrected
		  if(ecalIsoCutRhoCorr && trackIsoCutRhoCorr){
		    h_HcalIsoDR03Nminus1SignalRhoCorr->Fill(hcalIsoDR03);
		    h_HcalIsoDR04Nminus1SignalRhoCorr->Fill(hcalIsoDR04);
		  }
		  if(hcalIsoCutRhoCorr && trackIsoCutRhoCorr){
		    h_EcalIsoDR03Nminus1SignalRhoCorr->Fill(ecalIsoDR03);
		    h_EcalIsoDR04Nminus1SignalRhoCorr->Fill(ecalIsoDR04);
		  }
		  if(ecalIsoCutRhoCorr && hcalIsoCutRhoCorr){
		    h_TrackIsoDR03Nminus1SignalRhoCorr->Fill(trackIsoDR03);
		    h_TrackIsoDR04Nminus1SignalRhoCorr->Fill(trackIsoDR04);
		  } 
		  if(SmallestDR<900 && ecalIsoCutRhoCorr && hcalIsoCutRhoCorr && trackIsoCutRhoCorr){
		    h_DrNumerSignalDR04RhoCorr->Fill(SmallestDR);
		    h_EtNumerSignalDR04RhoCorr->Fill(PhoEt);
		  }
		  if(SmallestDR<900 && ecalIsoCutdR03RhoCorr && hcalIsoCutdR03RhoCorr && trackIsoCutdR03RhoCorr){
		    h_DrNumerSignalDR03RhoCorr->Fill(SmallestDR);
		    h_EtNumerSignalDR03RhoCorr->Fill(PhoEt);
		  }
		  if(genJetDR<900 && ecalIsoCutRhoCorr && hcalIsoCutRhoCorr && trackIsoCutRhoCorr){
		    h_DrGenNumerSignalDR04RhoCorr->Fill(genJetDR);
		    h_EtGenNumerSignalDR04RhoCorr->Fill(PhoEt);
		  }
		  if(genJetDR<900 && ecalIsoCutdR03RhoCorr && hcalIsoCutdR03RhoCorr && trackIsoCutdR03RhoCorr){
		    h_DrGenNumerSignalDR03RhoCorr->Fill(genJetDR);
		    h_EtGenNumerSignalDR03RhoCorr->Fill(PhoEt);
		  }
		}
	      }
	      
	      //now background
	      //cout<<"pdgId: "<<genPart->pdgId<<"\nstatus: "<<genPart->status<<"\nmotherId: "<<genPart->motherId<<endl<<endl;
	      if(genPart->pdgId!=22 && fabs(genPart->pdgId)!=11// && genPart->motherId>22/*genPart->motherId!=22 && fabs(genPart->motherId)!=21 && fabs(genPart->motherId)>=7*/ && genPart->motherId!=1000022/* || (genPart->pdgId==22 && (genPart->motherId<=6 || genPart->motherId==21))*/
		 ){
		if(isSameObject(it_Pho->caloPosition,genPart->momentum,0.1)){
		  float ptHatWeight=1./pow((event->gridParams["ptHat"]/15.),4.5);
		  if(SmallestDR<900){
		    h_DrDenomBground->Fill(SmallestDR,/*1*/LumiWeight);
		    h_EtDenomBground->Fill(PhoEt,/*1*/LumiWeight);
		  }
		  h_BgroundCombIsoDR03Nminus3->Fill(ecalIsoDR03+hcalIsoDR03+trackIsoDR03,/*1*/LumiWeight);
		  h_BgroundCombIsoDR04Nminus3->Fill(ecalIsoDR04+hcalIsoDR04+trackIsoDR04,/*1*/LumiWeight);
		  h_BgroundCombRelPFIsoDR03Nminus3->Fill((chargedHadronIso+neutralHadronIso+photonIso)/PhoEt,LumiWeight);
		  h_HcalIsoDR03Nminus3Bground->Fill(hcalIsoDR03,/*1*/LumiWeight);
		  h_EcalIsoDR03Nminus3Bground->Fill(ecalIsoDR03,/*1*/LumiWeight);
		  h_TrackIsoDR03Nminus3Bground->Fill(trackIsoDR03,/*1*/LumiWeight);
		  if(genJetDR<900){
		    h_DrGenDenomBground->Fill(genJetDR,/*1*/LumiWeight);
		    h_EtGenDenomBground->Fill(PhoEt,/*1*/LumiWeight);
		  }
		  if(ecalIsoCut && trackIsoCut){
		    h_HcalIsoDR03Nminus1Bground->Fill(hcalIsoDR03,/*1*/LumiWeight);
		    h_HcalIsoDR04Nminus1Bground->Fill(hcalIsoDR04,/*1*/LumiWeight);
		  }
		  if(hcalIsoCut && trackIsoCut){
		    h_EcalIsoDR03Nminus1Bground->Fill(ecalIsoDR03,/*1*/LumiWeight);
		    h_EcalIsoDR04Nminus1Bground->Fill(ecalIsoDR04,/*1*/LumiWeight);
		  }
		  if(ecalIsoCut && hcalIsoCut){
		    h_TrackIsoDR03Nminus1Bground->Fill(trackIsoDR03,/*1*/LumiWeight);
		    h_TrackIsoDR04Nminus1Bground->Fill(trackIsoDR04,/*1*/LumiWeight);
		  } 
		  
		  if(SmallestDR<900 && ecalIsoCut && hcalIsoCut && trackIsoCut){
		    h_DrNumerBgroundDR04->Fill(SmallestDR,/*1*/LumiWeight);
		    h_EtNumerBgroundDR04->Fill(PhoEt,/*1*/LumiWeight);
		  }
		  if(SmallestDR<900 && ecalIsoCutdR03 && hcalIsoCutdR03 && trackIsoCutdR03){
		    h_DrNumerBgroundDR03->Fill(SmallestDR,/*1*/LumiWeight);
		    h_EtNumerBgroundDR03->Fill(PhoEt,/*1*/LumiWeight);
		  }
		  if(genJetDR<900 && ecalIsoCut && hcalIsoCut && trackIsoCut){
		    h_DrGenNumerBgroundDR04->Fill(genJetDR,/*1*/LumiWeight);
		    h_EtGenNumerBgroundDR04->Fill(PhoEt,/*1*/LumiWeight);
		  }
		  if(genJetDR<900 && ecalIsoCutdR03 && hcalIsoCutdR03 && trackIsoCutdR03){
		    h_DrGenNumerBgroundDR03->Fill(genJetDR,/*1*/LumiWeight);
		    h_EtGenNumerBgroundDR03->Fill(PhoEt,/*1*/LumiWeight);
		  }
		  //now Rho corrected 
		  if(ecalIsoCutRhoCorr && trackIsoCutRhoCorr){
		    h_HcalIsoDR03Nminus1BgroundRhoCorr->Fill(hcalIsoDR03,/*1*/LumiWeight);
		    h_HcalIsoDR04Nminus1BgroundRhoCorr->Fill(hcalIsoDR04,/*1*/LumiWeight);
		  }
		  if(hcalIsoCutRhoCorr && trackIsoCutRhoCorr){
		    h_EcalIsoDR03Nminus1BgroundRhoCorr->Fill(ecalIsoDR03,/*1*/LumiWeight);
		    h_EcalIsoDR04Nminus1BgroundRhoCorr->Fill(ecalIsoDR04,/*1*/LumiWeight);
		  }
		  if(ecalIsoCutRhoCorr && hcalIsoCutRhoCorr){
		    h_TrackIsoDR03Nminus1BgroundRhoCorr->Fill(trackIsoDR03,/*1*/LumiWeight);
		    h_TrackIsoDR04Nminus1BgroundRhoCorr->Fill(trackIsoDR04,/*1*/LumiWeight);
		  } 
		  
		  if(SmallestDR<900 && ecalIsoCutRhoCorr && hcalIsoCutRhoCorr && trackIsoCutRhoCorr){
		    h_DrNumerBgroundDR04RhoCorr->Fill(SmallestDR,/*1*/LumiWeight);
		    h_EtNumerBgroundDR04RhoCorr->Fill(PhoEt,/*1*/LumiWeight);
		  }
		  if(SmallestDR<900 && ecalIsoCutdR03RhoCorr && hcalIsoCutdR03RhoCorr && trackIsoCutdR03RhoCorr){
		    h_DrNumerBgroundDR03RhoCorr->Fill(SmallestDR,/*1*/LumiWeight);
		    h_EtNumerBgroundDR03RhoCorr->Fill(PhoEt,/*1*/LumiWeight);
		  }
		  if(genJetDR<900 && ecalIsoCutRhoCorr && hcalIsoCutRhoCorr && trackIsoCutRhoCorr){
		    h_DrGenNumerBgroundDR04RhoCorr->Fill(genJetDR,/*1*/LumiWeight);
		    h_EtGenNumerBgroundDR04RhoCorr->Fill(PhoEt,/*1*/LumiWeight);
		  }
		  if(genJetDR<900 && ecalIsoCutdR03RhoCorr && hcalIsoCutdR03RhoCorr && trackIsoCutdR03RhoCorr){
		    h_DrGenNumerBgroundDR03RhoCorr->Fill(genJetDR,/*1*/LumiWeight);
		    h_EtGenNumerBgroundDR03RhoCorr->Fill(PhoEt,/*1*/LumiWeight);
		  }
		}
	      }//end bground
	    }//end if nearestDR<90
	  }
	}//end !isRealData
	
	if(/*event->*/isRealData){
	  if( etaCut && EtCut && !isSpike && heCut && sIetaCut && pixelCut/* && CaloIdLCut && (IsoVLCut || R9TriggerCut)*/ ){
	    if(met->met()<30){
	      h_BgroundCombIsoDR03Nminus3->Fill(ecalIsoDR03+hcalIsoDR03+trackIsoDR03);
	      h_BgroundCombIsoDR04Nminus3->Fill(ecalIsoDR04+hcalIsoDR04+trackIsoDR04);
	      h_chargedHadronIsoNminus3Bground->Fill(chargedHadronIso);
	      h_neutralHadronIsoNminus3Bground->Fill(neutralHadronIso);
	      h_photonIsoNminus3Bground->Fill(photonIso);
	      h_chargedHadronIsoDepositNminus3Bground->Fill(chargedHadronIsoDeposit);
	      h_neutralHadronIsoDepositNminus3Bground->Fill(neutralHadronIsoDeposit);
	      h_photonIsoDepositNminus3Bground->Fill(photonIsoDeposit);
	      h_PfCombinedIsoNminus3Bground->Fill(chargedHadronIso+neutralHadronIso+photonIso);
	      h_PfCombinedIsoDepositNminus3Bground->Fill(chargedHadronIsoDeposit+neutralHadronIsoDeposit+photonIsoDeposit);
	    }
	    h_EcalIsoDR03Nminus3->Fill(ecalIsoDR03);
	    h_HcalIsoDR03Nminus3->Fill(hcalIsoDR03);
	    h_TrackIsoDR03Nminus3->Fill(trackIsoDR03);
	    h_PfCombinedIsoNminus3->Fill(chargedHadronIso+neutralHadronIso+photonIso);
	    h_PfCombinedIsoDepositNminus3->Fill(chargedHadronIsoDeposit+neutralHadronIsoDeposit+photonIsoDeposit);
	    //	  h_dRphoJetNminus3->Fill(getDR());
	    if(ecalIsoCut && trackIsoCut){
	      h_HcalIsoDR03Nminus1->Fill(hcalIsoDR03);
	      h_HcalIsoDR04Nminus1->Fill(hcalIsoDR04);
	      if(met->met()<30){
		h_BgroundHcalIsoDR03Nminus1->Fill(hcalIsoDR03);
	      }
	    }
	    if(hcalIsoCut && trackIsoCut){
	      h_EcalIsoDR03Nminus1->Fill(ecalIsoDR03);
	      h_EcalIsoDR04Nminus1->Fill(ecalIsoDR04);
	      if(met->met()<30){
		h_BgroundEcalIsoDR03Nminus1->Fill(ecalIsoDR03);
	      }
	    }
	    if(ecalIsoCut && hcalIsoCut){
	      h_TrackIsoDR03Nminus1->Fill(trackIsoDR03);
	      h_TrackIsoDR04Nminus1->Fill(trackIsoDR04);
	      if(met->met()<30){
		h_BgroundTrackIsoDR03Nminus1->Fill(trackIsoDR03);
	      }
	    }
	  }
	}
	
	if(printLevel > 1) cout<<"End of Photon Loop"<<endl;
	
      }//for(it_Pho)
      
    }//if phoMap!=end (after phoMap def)

    //}//end diff check
  }// for jentry
 
  TH1F *h_DrSigEffDr04 = new TH1F("DrSigEffDr04","",64,0,4);
  h_DrSigEffDr04->Divide(h_DrNumerSignalDR04,h_DrDenomSignal,1,1,"B");
  TH1F *h_DrSigEffDr03 = new TH1F("DrSigEffDr03","",64,0,4);
  h_DrSigEffDr03->Divide(h_DrNumerSignalDR03,h_DrDenomSignal,1,1,"B");
  TH1F *h_DrGenSigEffDr04 = new TH1F("DrGenSigEffDr04","",64,0,4);
  h_DrGenSigEffDr04->Divide(h_DrGenNumerSignalDR04,h_DrGenDenomSignal,1,1,"B");
  TH1F *h_DrGenSigEffDr03 = new TH1F("DrGenSigEffDr03","",64,0,4);
  h_DrGenSigEffDr03->Divide(h_DrGenNumerSignalDR03,h_DrGenDenomSignal,1,1,"B");
  TH1F *h_DrBackEffDr04 = new TH1F("DrBackEffDr04","",64,0,4);
  h_DrBackEffDr04->Divide(h_DrNumerBgroundDR04,h_DrDenomBground,1,1,"B");
  TH1F *h_DrBackEffDr03 = new TH1F("DrBackEffDr03","",64,0,4);
  h_DrBackEffDr03->Divide(h_DrNumerBgroundDR03,h_DrGenDenomBground,1,1,"B");
  TH1F *h_DrGenBackEffDr04 = new TH1F("DrGenBackEffDr04","",64,0,4);
  h_DrGenBackEffDr04->Divide(h_DrGenNumerBgroundDR04,h_DrGenDenomBground,1,1,"B");
  TH1F *h_DrGenBackEffDr03 = new TH1F("DrGenBackEffDr03","",64,0,4);
  h_DrGenBackEffDr03->Divide(h_DrGenNumerBgroundDR03,h_DrGenDenomBground,1,1,"B");
  TH1F *h_EtSigEffDr04 = new TH1F("EtSigEffDr04","",100,0,500);
  h_EtSigEffDr04->Divide(h_EtNumerSignalDR04,h_EtDenomSignal,1,1,"B");
  TH1F *h_EtSigEffDr03 = new TH1F("EtSigEffDr03","",100,0,500);
  h_EtSigEffDr03->Divide(h_EtNumerSignalDR03,h_EtDenomSignal,1,1,"B");
  TH1F *h_EtGenSigEffDr04 = new TH1F("EtGenSigEffDr04","",100,0,500);
  h_EtGenSigEffDr04->Divide(h_EtGenNumerSignalDR04,h_EtGenDenomSignal,1,1,"B");
  TH1F *h_EtGenSigEffDr03 = new TH1F("EtGenSigEffDr03","",100,0,500);
  h_EtGenSigEffDr03->Divide(h_EtGenNumerSignalDR03,h_EtGenDenomSignal,1,1,"B");
  TH1F *h_EtBackEffDr04 = new TH1F("EtBackEffDr04","",100,0,500);
  h_EtBackEffDr04->Divide(h_EtNumerBgroundDR04,h_EtDenomBground,1,1,"B");
  TH1F *h_EtBackEffDr03 = new TH1F("EtBackEffDr03","",100,0,500);
  h_EtBackEffDr03->Divide(h_EtNumerBgroundDR03,h_EtDenomBground,1,1,"B");
  TH1F *h_EtGenBackEffDr04 = new TH1F("EtGenBackEffDr04","",100,0,500);
  h_EtGenBackEffDr04->Divide(h_EtGenNumerBgroundDR04,h_EtGenDenomBground,1,1,"B");
  TH1F *h_EtGenBackEffDr03 = new TH1F("EtGenBackEffDr03","",100,0,500);
  h_EtGenBackEffDr03->Divide(h_EtGenNumerBgroundDR03,h_EtGenDenomBground,1,1,"B");
  //now Rho corrected
  TH1F *h_DrSigEffDr04RhoCorr = new TH1F("DrSigEffDr04RhoCorr","",64,0,4);
  h_DrSigEffDr04RhoCorr->Divide(h_DrNumerSignalDR04RhoCorr,h_DrDenomSignal,1,1,"B");
  TH1F *h_DrSigEffDr03RhoCorr = new TH1F("DrSigEffDr03RhoCorr","",64,0,4);
  h_DrSigEffDr03RhoCorr->Divide(h_DrNumerSignalDR03RhoCorr,h_DrDenomSignal,1,1,"B");
  TH1F *h_DrGenSigEffDr04RhoCorr = new TH1F("DrGenSigEffDr04RhoCorr","",64,0,4);
  h_DrGenSigEffDr04RhoCorr->Divide(h_DrGenNumerSignalDR04RhoCorr,h_DrGenDenomSignal,1,1,"B");
  TH1F *h_DrGenSigEffDr03RhoCorr = new TH1F("DrGenSigEffDr03RhoCorr","",64,0,4);
  h_DrGenSigEffDr03RhoCorr->Divide(h_DrGenNumerSignalDR03RhoCorr,h_DrGenDenomSignal,1,1,"B");
  TH1F *h_DrBackEffDr04RhoCorr = new TH1F("DrBackEffDr04RhoCorr","",64,0,4);
  h_DrBackEffDr04RhoCorr->Divide(h_DrNumerBgroundDR04RhoCorr,h_DrDenomBground,1,1,"B");
  TH1F *h_DrBackEffDr03RhoCorr = new TH1F("DrBackEffDr03RhoCorr","",64,0,4);
  h_DrBackEffDr03RhoCorr->Divide(h_DrNumerBgroundDR03RhoCorr,h_DrGenDenomBground,1,1,"B");
  TH1F *h_DrGenBackEffDr04RhoCorr = new TH1F("DrGenBackEffDr04RhoCorr","",64,0,4);
  h_DrGenBackEffDr04RhoCorr->Divide(h_DrGenNumerBgroundDR04RhoCorr,h_DrGenDenomBground,1,1,"B");
  TH1F *h_DrGenBackEffDr03RhoCorr = new TH1F("DrGenBackEffDr03RhoCorr","",64,0,4);
  h_DrGenBackEffDr03RhoCorr->Divide(h_DrGenNumerBgroundDR03RhoCorr,h_DrGenDenomBground,1,1,"B");
  TH1F *h_EtSigEffDr04RhoCorr = new TH1F("EtSigEffDr04RhoCorr","",100,0,500);
  h_EtSigEffDr04RhoCorr->Divide(h_EtNumerSignalDR04RhoCorr,h_EtDenomSignal,1,1,"B");
  TH1F *h_EtSigEffDr03RhoCorr = new TH1F("EtSigEffDr03RhoCorr","",100,0,500);
  h_EtSigEffDr03RhoCorr->Divide(h_EtNumerSignalDR03RhoCorr,h_EtDenomSignal,1,1,"B");
  TH1F *h_EtGenSigEffDr04RhoCorr = new TH1F("EtGenSigEffDr04RhoCorr","",100,0,500);
  h_EtGenSigEffDr04RhoCorr->Divide(h_EtGenNumerSignalDR04RhoCorr,h_EtGenDenomSignal,1,1,"B");
  TH1F *h_EtGenSigEffDr03RhoCorr = new TH1F("EtGenSigEffDr03RhoCorr","",100,0,500);
  h_EtGenSigEffDr03RhoCorr->Divide(h_EtGenNumerSignalDR03RhoCorr,h_EtGenDenomSignal,1,1,"B");
  TH1F *h_EtBackEffDr04RhoCorr = new TH1F("EtBackEffDr04RhoCorr","",100,0,500);
  h_EtBackEffDr04RhoCorr->Divide(h_EtNumerBgroundDR04RhoCorr,h_EtDenomBground,1,1,"B");
  TH1F *h_EtBackEffDr03RhoCorr = new TH1F("EtBackEffDr03RhoCorr","",100,0,500);
  h_EtBackEffDr03RhoCorr->Divide(h_EtNumerBgroundDR03RhoCorr,h_EtDenomBground,1,1,"B");
  TH1F *h_EtGenBackEffDr04RhoCorr = new TH1F("EtGenBackEffDr04RhoCorr","",100,0,500);
  h_EtGenBackEffDr04RhoCorr->Divide(h_EtGenNumerBgroundDR04RhoCorr,h_EtGenDenomBground,1,1,"B");
  TH1F *h_EtGenBackEffDr03RhoCorr = new TH1F("EtGenBackEffDr03RhoCorr","",100,0,500);
  h_EtGenBackEffDr03RhoCorr->Divide(h_EtGenNumerBgroundDR03RhoCorr,h_EtGenDenomBground,1,1,"B");

  cout << " ----------------- Job Summary ----------------- " << endl;
  cout << " Total events            : " << nCnt[0] << " (" << 100*nCnt[0]/float(nCnt[0]) << "%)" << endl;
  cout << " HLT passed              : " << nCnt[1] << " (" << 100*nCnt[1]/float(nCnt[0]) << "%)" << endl;
 
  cout<<"Writing root output to: /data/ndpc3/b/dmorse/RA3/AnalysisOutput/hist_"<<ds<<".root"<<endl;
  //cout<<"Writing root output to: /Users/dmorse/RA3/AnalysisOutput/hist_"<<ds<<".root"<<endl;

  // close the output file

  fout->cd();
  fout->Write();
  fout->Close();
  delete fout;

}

//---------------------------------------------
//       PILEUP!!!!!!
//---------------------------------------------

void SusyEventAnalyzer::Pileup() {
  
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  cout << "total events in files  : " << nentries << endl;

  if(processNEvents <= 0 || processNEvents > nentries) processNEvents = nentries;

  cout << "events to be processed : " << processNEvents << endl; 
  


  if(printLevel > 1) cout << "Initialize event counters." << endl;
  const int NCNT = 20;
  int nCnt[NCNT];
  for(int i=0; i<NCNT; i++) nCnt[i] = 0;
  int ne=0,n2e=0,n3e=0;

  // open hist file and define histograms

  TFile* fout = new TFile("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/Pileup/cms525v3_jec2012/hist_"+ds+".root","RECREATE");
  //TFile* fout = new TFile("hist_"+ds+".root","RECREATE");
  //TFile* fout = new TFile("/Users/dmorse/RA3/AnalysisOutput/hist_"+ds+".root","RECREATE");

  fout->cd();

  ofstream EventWriteOut;
 
  //--------------------Pileup studies--------
  TH1F* h_rho = new TH1F("rho","Fastjet correction rho for all events that pass JSON, HLT",150,0,150);h_rho->Sumw2(); 
  TH1F* h_rho25 = new TH1F("rho25","Fastjet correction rho25 for all events that pass JSON, HLT",150,0,150);h_rho25->Sumw2();
  TH1F* h_rhoB = new TH1F("rhoB","Fastjet correction rho with Rho_EtaMax=1.4442 for all events that pass JSON, HLT",150,0,150);h_rhoB->Sumw2();
  TH2F* h_NvertexVsRho25 = new TH2F("NvertexVsRho25","Number of good Vertices vs. FastJet rho for all events that pass JSON, HLT;rho;nVertex",50,0,50,50,0,50);h_NvertexVsRho25->Sumw2();
  TH2F* h_NvertexVsRhoB = new TH2F("NvertexVsRhoB","Number of good Vertices vs. FastJet rho with Rho_EtaMax=1.4442 for all events that pass JSON, HLT;nVertex;rho",50,0,50,50,0,50);h_NvertexVsRhoB->Sumw2();
  //------------
  TH2F* h_HOverEVsRho25 = new TH2F("HOverEVsRho25","H/E of barrel photons VS Rho after Et, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;Rho (Et/ unit area);h/E",50,0,50,41,-.05,.2);h_HOverEVsRho25->Sumw2(); 
  TH2F* h_HOverEVsRhoB = new TH2F("HOverEVsRhoB","H/E of barrel photons VS RhoB after Et, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;RhoB (Et/ unit area);h/E",50,0,50,41,-.05,.2);h_HOverEVsRhoB->Sumw2(); 
  TH2F* h_HOverEVsNVertex = new TH2F("HOverEVsNVertex","H/E of barrel photons VS NVertex after Et, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;NVertex ;h/E",50,0,50,41,-.05,.2);h_HOverEVsNVertex->Sumw2();

  TH2F* h_EcalIsoVsRho25DR03 = new TH2F("EcalIsoVsRho25DR03","Ecal Isolation of barrel photons VS Rho after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;Rho (Et/ unit area);ecalRecHitSumEtConeDR03",50,0,50,2700,-40,500);h_EcalIsoVsRho25DR03->Sumw2();
  TH2F* h_HcalIsoVsRho25DR03 = new TH2F("HcalIsoVsRho25DR03"," Hcal Isolation of barrel photons VS Rho after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;Rho (Et/ unit area);hcalTowerSumEtConeDR03",50,0,50,1525,-5,300);h_HcalIsoVsRho25DR03->Sumw2();
  TH2F* h_TrackIsoVsRho25DR03 = new TH2F("TrackIsoVsRho25DR03"," Tracker Isolation of barrel photons VS Rho after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;Rho (Et/ unit area);trkSumPtHollowConeDR03",50,0,50,2025,-5,400); h_TrackIsoVsRho25DR03->Sumw2();
  TH2F* h_CombIsoVsRho25DR03 = new TH2F("CombIsoVsRho25DR03"," Combined Isolation of barrel photons VS Rho after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;Rho (Et/ unit area);ecalRecHitSumEtConeDR03+hcalTowerSumEtConeDR03+trkSumPtHollowConeDR03",50,0,50,5200,-40,1000);h_CombIsoVsRho25DR03->Sumw2();


  TH2F* h_chargedHadronIsoVsRho25_ee = new TH2F("chargedHadronIsoVsRho25_ee","chargedHadron Isolation of barrel photons VS Rho25 after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;Rho (E_{T}/ unit area);chargedHadronIso",50,0,50,2700,-40,500);h_chargedHadronIsoVsRho25_ee->Sumw2();
  TH2F* h_neutralHadronIsoVsRho25_ee = new TH2F("neutralHadronIsoVsRho25_ee","neutralHadron Isolation of barrel photons VS Rho25 after E_{T}, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;Rho (E_{T}/ unit area);neutralHadronIso",50,0,50,2700,-40,500);h_neutralHadronIsoVsRho25_ee->Sumw2();
  TH2F* h_photonIsoVsRho25_ee = new TH2F("photonIsoVsRho25_ee","photon Isolation of barrel photons VS Rho25 after E_{T}, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;Rho25 (E_{T}/ unit area);photonIso",50,0,50,2700,-40,500);h_photonIsoVsRho25_ee->Sumw2();
  TH2F* h_PfCombIsoVsRho25_ee = new TH2F("PfCombIsoVsRho25_ee","PfComb Isolation of barrel photons VS Rho25 after E_{T}, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;Rho25 (E_{T}/ unit area);PfCombIso",50,0,50,2700,-40,500);h_PfCombIsoVsRho25_ee->Sumw2();

  //single pho
  TH2F* h_EcalIsoVsRho25DR03_SinglePho = new TH2F("EcalIsoVsRho25DR03_SinglePho","Ecal Isolation of barrel photons VS Rho after Et, h/e, eta, r9<1.0, sIetaIeta<0.011, pixel=0 cuts;Rho (Et/ unit area);ecalRecHitSumEtConeDR03",50,0,50,2700,-40,500);h_EcalIsoVsRho25DR03_SinglePho->Sumw2();
  TH2F* h_HcalIsoVsRho25DR03_SinglePho = new TH2F("HcalIsoVsRho25DR03_SinglePho"," Hcal Isolation of barrel photons VS Rho after Et, h/e, eta, r9<1.0, sIetaIeta<0.011, pixel=0 cuts;Rho (Et/ unit area);hcalTowerSumEtConeDR03",50,0,50,1525,-5,300);h_HcalIsoVsRho25DR03_SinglePho->Sumw2();
  TH2F* h_TrackIsoVsRho25DR03_SinglePho = new TH2F("TrackIsoVsRho25DR03_SinglePho"," Tracker Isolation of barrel photons VS Rho after Et, h/e, eta, r9<1.0, sIetaIeta<0.011, pixel=0 cuts;Rho (Et/ unit area);trkSumPtHollowConeDR03_SinglePho",50,0,50,2025,-5,400); h_TrackIsoVsRho25DR03_SinglePho->Sumw2();
  TH2F* h_CombIsoVsRho25DR03_SinglePho = new TH2F("CombIsoVsRho25DR03_SinglePho"," Combined Isolation of barrel photons VS Rho after Et, h/e, eta, r9<1.0, sIetaIeta<0.011, pixel=0 cuts;Rho (Et/ unit area);ecalRecHitSumEtConeDR03+hcalTowerSumEtConeDR03+trkSumPtHollowConeDR03",50,0,50,5200,-40,1000);h_CombIsoVsRho25DR03_SinglePho->Sumw2();
  TH2F* h_EcalIsoVsRhoBDR03_SinglePho = new TH2F("EcalIsoVsRhoBDR03_SinglePho","Ecal Isolation of barrel photons VS RhoB after Et, h/e, eta, r9<1.0, sIetaIeta<0.011, pixel=0 cuts;RhoB (Et/ unit area);ecalRecHitSumEtConeDR03",50,0,50,2700,-40,500);h_EcalIsoVsRhoBDR03_SinglePho->Sumw2();
  TH2F* h_HcalIsoVsRhoBDR03_SinglePho = new TH2F("HcalIsoVsRhoBDR03_SinglePho"," Hcal Isolation of barrel photons VS RhoB after Et, h/e, eta, r9<1.0, sIetaIeta<0.011, pixel=0 cuts;RhoB (Et/ unit area);hcalTowerSumEtConeDR03",50,0,50,1525,-5,300);h_HcalIsoVsRhoBDR03_SinglePho->Sumw2();
  TH2F* h_TrackIsoVsRhoBDR03_SinglePho = new TH2F("TrackIsoVsRhoBDR03_SinglePho"," Tracker Isolation of barrel photons VS RhoB after Et, h/e, eta, r9<1.0, sIetaIeta<0.011, pixel=0 cuts;RhoB (Et/ unit area);trkSumPtHollowConeDR03_SinglePho",50,0,50,2025,-5,400); h_TrackIsoVsRhoBDR03_SinglePho->Sumw2();
  TH2F* h_CombIsoVsRhoBDR03_SinglePho = new TH2F("CombIsoVsRhoBDR03_SinglePho"," Combined Isolation of barrel photons VS RhoB after Et, h/e, eta, r9<1.0, sIetaIeta<0.011, pixel=0 cuts;RhoB (Et/ unit area);ecalRecHitSumEtConeDR03+hcalTowerSumEtConeDR03+trkSumPtHollowConeDR03",50,0,50,5200,-40,1000);h_CombIsoVsRhoBDR03_SinglePho->Sumw2();

  TH2F* h_EcalIsoVsRho25DR04 = new TH2F("EcalIsoVsRho25DR04","Ecal Isolation of barrel photons VS Rho after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;Rho (Et/ unit area);ecalRecHitSumEtConeDR04",50,0,50,2700,-40,500);h_EcalIsoVsRho25DR04->Sumw2();
  TH2F* h_HcalIsoVsRho25DR04 = new TH2F("HcalIsoVsRho25DR04"," Hcal Isolation of barrel photons VS Rho after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;Rho (Et/ unit area);hcalTowerSumEtConeDR04",50,0,50,1525,-5,300);h_HcalIsoVsRho25DR04->Sumw2();
  TH2F* h_TrackIsoVsRho25DR04 = new TH2F("TrackIsoVsRho25DR04"," Tracker Isolation of barrel photons VS Rho after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;Rho (Et/ unit area);trkSumPtHollowConeDR04",50,0,50,2025,-5,400); h_TrackIsoVsRho25DR04->Sumw2();
  TH2F* h_CombIsoVsRho25DR04 = new TH2F("CombIsoVsRho25DR04"," Combined Isolation of barrel photons VS Rho after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;Rho (Et/ unit area);ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04",50,0,50,5200,-40,1000);h_CombIsoVsRho25DR04->Sumw2();

  TH2F* h_EcalIsoVsRhoBDR03 = new TH2F("EcalIsoVsRhoBDR03","Ecal Isolation of barrel photons VS RhoB after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;RhoB (Et/ unit area);ecalRecHitSumEtConeDR03",50,0,50,2700,-40,500);h_EcalIsoVsRhoBDR03->Sumw2();
  TH2F* h_HcalIsoVsRhoBDR03 = new TH2F("HcalIsoVsRhoBDR03"," Hcal Isolation of barrel photons VS RhoB after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;RhoB (Et/ unit area);hcalTowerSumEtConeDR03",50,0,50,1525,-5,300);h_HcalIsoVsRhoBDR03->Sumw2();
  TH2F* h_TrackIsoVsRhoBDR03 = new TH2F("TrackIsoVsRhoBDR03"," Tracker Isolation of barrel photons VS RhoB after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;RhoB (Et/ unit area);trkSumPtHollowConeDR03",50,0,50,2025,-5,400); h_TrackIsoVsRhoBDR03->Sumw2();
  TH2F* h_CombIsoVsRhoBDR03 = new TH2F("CombIsoVsRhoBDR03"," Combined Isolation of barrel photons VS RhoB after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;RhoB (Et/ unit area);ecalRecHitSumEtConeDR03+hcalTowerSumEtConeDR03+trkSumPtHollowConeDR03",50,0,50,5200,-40,1000);h_CombIsoVsRhoBDR03->Sumw2();

  TH2F* h_EcalIsoVsRhoBDR04 = new TH2F("EcalIsoVsRhoBDR04","Ecal Isolation of barrel photons VS RhoB after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;RhoB (Et/ unit area);ecalRecHitSumEtConeDR04",50,0,50,2700,-40,500);h_EcalIsoVsRhoBDR04->Sumw2();
  TH2F* h_HcalIsoVsRhoBDR04 = new TH2F("HcalIsoVsRhoBDR04"," Hcal Isolation of barrel photons VS RhoB after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;RhoB (Et/ unit area);hcalTowerSumEtConeDR04",50,0,50,1525,-5,300);h_HcalIsoVsRhoBDR04->Sumw2();
  TH2F* h_TrackIsoVsRhoBDR04 = new TH2F("TrackIsoVsRhoBDR04"," Tracker Isolation of barrel photons VS RhoB after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;RhoB (Et/ unit area);trkSumPtHollowConeDR04",50,0,50,2025,-5,400); h_TrackIsoVsRhoBDR04->Sumw2();
  TH2F* h_CombIsoVsRhoBDR04 = new TH2F("CombIsoVsRhoBDR04"," Combined Isolation of barrel photons VS RhoB after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;RhoB (Et/ unit area);ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04",50,0,50,5200,-40,1000);h_CombIsoVsRhoBDR04->Sumw2();

  TH2F* h_EcalIsoVsNVertexDR03 = new TH2F("EcalIsoVsNVertexDR03","Ecal Isolation of barrel photons VS NVertex after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;NVertex ;ecalRecHitSumEtConeDR03",50,0,50,2700,-40,500);h_EcalIsoVsNVertexDR03->Sumw2();
  TH2F* h_HcalIsoVsNVertexDR03 = new TH2F("HcalIsoVsNVertexDR03"," Hcal Isolation of barrel photons VS NVertex after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;NVertex ;hcalTowerSumEtConeDR03",50,0,50,1525,-5,300);h_HcalIsoVsNVertexDR03->Sumw2();
  TH2F* h_TrackIsoVsNVertexDR03 = new TH2F("TrackIsoVsNVertexDR03"," Tracker Isolation of barrel photons VS NVertex after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;NVertex ;trkSumPtHollowConeDR03",50,0,50,2025,-5,400); h_TrackIsoVsNVertexDR03->Sumw2();
  TH2F* h_CombIsoVsNVertexDR03 = new TH2F("CombIsoVsNVertexDR03"," Combined Isolation of barrel photons VS NVertex after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;NVertex ;ecalRecHitSumEtConeDR03+hcalTowerSumEtConeDR03+trkSumPtHollowConeDR03",50,0,50,5200,-40,1000);h_CombIsoVsNVertexDR03->Sumw2();

  TH2F* h_EcalIsoVsNVertexDR04 = new TH2F("EcalIsoVsNVertexDR04","Ecal Isolation of barrel photons VS NVertex after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;NVertex ;ecalRecHitSumEtConeDR04",50,0,50,2700,-40,500);h_EcalIsoVsNVertexDR04->Sumw2();
  TH2F* h_HcalIsoVsNVertexDR04 = new TH2F("HcalIsoVsNVertexDR04"," Hcal Isolation of barrel photons VS NVertex after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;NVertex ;hcalTowerSumEtConeDR04",50,0,50,1525,-5,300);h_HcalIsoVsNVertexDR04->Sumw2();
  TH2F* h_TrackIsoVsNVertexDR04 = new TH2F("TrackIsoVsNVertexDR04"," Tracker Isolation of barrel photons VS NVertex after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;NVertex ;trkSumPtHollowConeDR04",50,0,50,2025,-5,400); h_TrackIsoVsNVertexDR04->Sumw2();
  TH2F* h_CombIsoVsNVertexDR04 = new TH2F("CombIsoVsNVertexDR04"," Combined Isolation of barrel photons VS NVertex after Et, h/e, eta, r9<.98, CaloIdL, CaloIsoVL||R9>.8 cuts;NVertex ;ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04",50,0,50,5200,-40,1000);h_CombIsoVsNVertexDR04->Sumw2();

  TH2F* h_HcalIsoVsEcalIsoDR03 = new TH2F("HcalIsoVsHcalIsoDR03"," Hcal Isolation VS Ecal Isolation of barrel photons after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;ecalRecHitSumEtConeDR03;hcalTowerSumEtConeDR03",3250,-50,600,1525,-5,300);h_HcalIsoVsEcalIsoDR03->Sumw2();
  TH2F* h_TrackIsoVsEcalIsoDR03 = new TH2F("TrackIsoVsEcalIsoDR03"," Track Isolation VS Ecal Isolation of barrel photons after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;ecalRecHitSumEtConeDR03;trkSumPtHollowConeDR03",3250,-50,600,1525,-5,300);h_TrackIsoVsEcalIsoDR03->Sumw2();
  TH2F* h_TrackIsoVsHcalIsoDR03 = new TH2F("TrackIsoVsHcalIsoDR03"," Track Isolation VS Hcal Isolation of barrel photons after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;hcalTowerSumEtConeDR03;trkSumPtHollowConeDR03",3250,-50,600,1525,-5,300);h_TrackIsoVsHcalIsoDR03->Sumw2();
  
  TH2F* h_EcalIsoVsRho25DR03_lt6NV = new TH2F("EcalIsoVsRho25DR03_lt6NV","Ecal Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts in events with 1<NVertex<=5;Rho (Et/ unit area);ecalRecHitSumEtConeDR03",50,0,50,3250,-50,600);h_EcalIsoVsRho25DR03_lt6NV->Sumw2();
  TH2F* h_HcalIsoVsRho25DR03_lt6NV = new TH2F("HcalIsoVsRho25DR03_lt6NV"," Hcal Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts in events with 1<NVertex<=5;Rho (Et/ unit area);hcalTowerSumEtConeDR03",50,0,50,1525,-5,300);h_HcalIsoVsRho25DR03_lt6NV->Sumw2();
  TH2F* h_EcalIsoVsRho25DR03_6to10NV = new TH2F("EcalIsoVsRho25DR03_6to10NV","Ecal Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts in events with 6<NVertex<=10;Rho (Et/ unit area);ecalRecHitSumEtConeDR03",50,0,50,3250,-50,600);h_EcalIsoVsRho25DR03_6to10NV->Sumw2();
  TH2F* h_HcalIsoVsRho25DR03_6to10NV = new TH2F("HcalIsoVsRho25DR03_6to10NV"," Hcal Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts in events with 6<NVertex<=10;Rho (Et/ unit area);hcalTowerSumEtConeDR03",50,0,50,1525,-5,300);h_HcalIsoVsRho25DR03_6to10NV->Sumw2();
  TH2F* h_EcalIsoVsRho25DR03_11to15NV = new TH2F("EcalIsoVsRho25DR03_11to15NV","Ecal Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts in events with 11<NVertex<=15;Rho (Et/ unit area);ecalRecHitSumEtConeDR03",50,0,50,3250,-50,600);h_EcalIsoVsRho25DR03_11to15NV->Sumw2();
  TH2F* h_HcalIsoVsRho25DR03_11to15NV = new TH2F("HcalIsoVsRho25DR03_11to15NV"," Hcal Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts in events with 11<NVertex<=15;Rho (Et/ unit area);hcalTowerSumEtConeDR03",50,0,50,1525,-5,300);h_HcalIsoVsRho25DR03_11to15NV->Sumw2();
  TH2F* h_EcalIsoVsRho25DR03_gt15NV = new TH2F("EcalIsoVsRho25DR03_gt15NV","Ecal Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts in events with NVertex>15;Rho (Et/ unit area);ecalRecHitSumEtConeDR03",50,0,50,3250,-50,600);h_EcalIsoVsRho25DR03_gt15NV->Sumw2();
  TH2F* h_HcalIsoVsRho25DR03_gt15NV = new TH2F("HcalIsoVsRho25DR03_gt15NV"," Hcal Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts in events with NVertex>15;Rho (Et/ unit area);hcalTowerSumEtConeDR03",50,0,50,1525,-5,300);h_HcalIsoVsRho25DR03_gt15NV->Sumw2();
  


  TH2F* h_EcalIsoVsRho25DR03_ee = new TH2F("EcalIsoVsRho25DR03_ee","Ecal Isolation of barrel probe electron VS Rho25 after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);ecalRecHitSumEtConeDR03",50,0,50,3250,-50,600);h_EcalIsoVsRho25DR03_ee->Sumw2();
  TH2F* h_HcalIsoVsRho25DR03_ee = new TH2F("HcalIsoVsRho25DR03_ee"," Hcal Isolation of barrel probe electron VS Rho25 after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);hcalTowerSumEtConeDR03",50,0,50,1525,-5,300);h_HcalIsoVsRho25DR03_ee->Sumw2();
  TH2F* h_TrackIsoVsRho25DR03_ee = new TH2F("TrackIsoVsRho25DR03_ee"," Tracker Isolation of barrel probe electron VS Rho25 after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);trkSumPtHollowConeDR03",50,0,50,3025,-5,600); h_TrackIsoVsRho25DR03_ee->Sumw2();
  TH2F* h_CombIsoVsRho25DR03_ee = new TH2F("CombIsoVsRho25DR03_ee"," Combined Isolation of barrel probe electron VS Rho25 after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);ecalRecHitSumEtConeDR03+hcalTowerSumEtConeDR03+trkSumPtHollowConeDR03",50,0,50,6250,-50,1200);h_CombIsoVsRho25DR03_ee->Sumw2();
  TH2F* h_EcalIsoVsRho25DR04_ee = new TH2F("EcalIsoVsRho25DR04_ee","Ecal Isolation of barrel probe electron VS Rho25 after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);ecalRecHitSumEtConeDR04",50,0,50,3250,-50,600);h_EcalIsoVsRho25DR04_ee->Sumw2();
  TH2F* h_HcalIsoVsRho25DR04_ee = new TH2F("HcalIsoVsRho25DR04_ee"," Hcal Isolation of barrel probe electron VS Rho25 after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);hcalTowerSumEtConeDR04",50,0,50,1525,-5,300);h_HcalIsoVsRho25DR04_ee->Sumw2();
  TH2F* h_TrackIsoVsRho25DR04_ee = new TH2F("TrackIsoVsRho25DR04_ee"," Tracker Isolation of barrel probe electron VS Rho25 after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);trkSumPtHollowConeDR04",50,0,50,3025,-5,600); h_TrackIsoVsRho25DR04_ee->Sumw2();
  TH2F* h_CombIsoVsRho25DR04_ee = new TH2F("CombIsoVsRho25DR04_ee"," Combined Isolation of barrel probe electron VS Rho25 after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04",50,0,50,6250,-50,1200);h_CombIsoVsRho25DR04_ee->Sumw2();
  TH2F* h_HOverEVsRho25_ee = new TH2F("HOverEVsRho25_ee","H/E of barrel probe electron VS Rho25 after Et, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);H/E",50,0,50,41,-.05,.2);h_HOverEVsRho25_ee->Sumw2(); 
  TH2F* h_HOverEAfterIsoVsRho25_ee = new TH2F("HOverEAfterIsoVsRho25_ee","H/E of barrel probe electron VS Rho25 after Et, r9<1.0, CaloIdL, IsoVL||R9>.8 CombIso<6 cuts;Rho (Et/ unit area);H/E",50,0,50,41,-.05,.2);h_HOverEAfterIsoVsRho25_ee->Sumw2(); 

  TH2F* h_EcalIsoVsRhoBDR03_ee = new TH2F("EcalIsoVsRhoBDR03_ee","Ecal Isolation of barrel probe electron VS RhoB after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;RhoB (Et/ unit area);ecalRecHitSumEtConeDR03",50,0,50,3250,-50,600);h_EcalIsoVsRhoBDR03_ee->Sumw2();
  TH2F* h_HcalIsoVsRhoBDR03_ee = new TH2F("HcalIsoVsRhoBDR03_ee"," Hcal Isolation of barrel probe electron VS RhoB after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;RhoB (Et/ unit area);hcalTowerSumEtConeDR03",50,0,50,1525,-5,300);h_HcalIsoVsRhoBDR03_ee->Sumw2();
  TH2F* h_TrackIsoVsRhoBDR03_ee = new TH2F("TrackIsoVsRhoBDR03_ee"," Tracker Isolation of barrel probe electron VS RhoB after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;RhoB (Et/ unit area);trkSumPtHollowConeDR03",50,0,50,3025,-5,600); h_TrackIsoVsRhoBDR03_ee->Sumw2();
  TH2F* h_CombIsoVsRhoBDR03_ee = new TH2F("CombIsoVsRhoBDR03_ee"," Combined Isolation of barrel probe electron VS RhoB after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;RhoB (Et/ unit area);ecalRecHitSumEtConeDR03+hcalTowerSumEtConeDR03+trkSumPtHollowConeDR03",50,0,50,6250,-50,1200);h_CombIsoVsRhoBDR03_ee->Sumw2();
  

  TH2F* h_EcalIsoVsRho25DR03_Pt0to40 = new TH2F("EcalIsoVsRho25DR03_Pt0to40","Ecal Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);ecalRecHitSumEtConeDR03",50,0,50,3250,-50,600);h_EcalIsoVsRho25DR03_Pt0to40->Sumw2();
  TH2F* h_HcalIsoVsRho25DR03_Pt0to40 = new TH2F("HcalIsoVsRho25DR03_Pt0to40"," Hcal Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);hcalTowerSumEtConeDR03",50,0,50,1525,-5,300);h_HcalIsoVsRho25DR03_Pt0to40->Sumw2();
  TH2F* h_TrackIsoVsRho25DR03_Pt0to40 = new TH2F("TrackIsoVsRho25DR03_Pt0to40"," Tracker Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);trkSumPtHollowConeDR03",50,0,50,3025,-5,600); h_TrackIsoVsRho25DR03_Pt0to40->Sumw2();
  TH2F* h_CombIsoVsRho25DR03_Pt0to40 = new TH2F("CombIsoVsRho25DR03_Pt0to40"," Combined Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);ecalRecHitSumEtConeDR03+hcalTowerSumEtConeDR03+trkSumPtHollowConeDR03",50,0,50,6250,-50,1200);h_CombIsoVsRho25DR03_Pt0to40->Sumw2();
  TH2F* h_EcalIsoVsRho25DR03_Pt40to50= new TH2F("EcalIsoVsRho25DR03_Pt40to50","Ecal Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);ecalRecHitSumEtConeDR03",50,0,50,3250,-50,600);h_EcalIsoVsRho25DR03_Pt40to50->Sumw2();
  TH2F* h_HcalIsoVsRho25DR03_Pt40to50= new TH2F("HcalIsoVsRho25DR03_Pt40to50"," Hcal Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);hcalTowerSumEtConeDR03",50,0,50,1525,-5,300);h_HcalIsoVsRho25DR03_Pt40to50->Sumw2();
  TH2F* h_TrackIsoVsRho25DR03_Pt40to50= new TH2F("TrackIsoVsRho25DR03_Pt40to50"," Tracker Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);trkSumPtHollowConeDR03",50,0,50,3025,-5,600); h_TrackIsoVsRho25DR03_Pt40to50->Sumw2();
  TH2F* h_CombIsoVsRho25DR03_Pt40to50= new TH2F("CombIsoVsRho25DR03_Pt40to50"," Combined Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);ecalRecHitSumEtConeDR03+hcalTowerSumEtConeDR03+trkSumPtHollowConeDR03",50,0,50,6250,-50,1200);h_CombIsoVsRho25DR03_Pt40to50->Sumw2();
  TH2F* h_EcalIsoVsRho25DR03_Pt50toINF= new TH2F("EcalIsoVsRho25DR03_Pt50toINF","Ecal Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);ecalRecHitSumEtConeDR03",50,0,50,3250,-50,600);h_EcalIsoVsRho25DR03_Pt50toINF->Sumw2();
  TH2F* h_HcalIsoVsRho25DR03_Pt50toINF= new TH2F("HcalIsoVsRho25DR03_Pt50toINF"," Hcal Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);hcalTowerSumEtConeDR03",50,0,50,1525,-5,300);h_HcalIsoVsRho25DR03_Pt50toINF->Sumw2();
  TH2F* h_TrackIsoVsRho25DR03_Pt50toINF= new TH2F("TrackIsoVsRho25DR03_Pt50toINF"," Tracker Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);trkSumPtHollowConeDR03",50,0,50,3025,-5,600); h_TrackIsoVsRho25DR03_Pt50toINF->Sumw2();
  TH2F* h_CombIsoVsRho25DR03_Pt50toINF= new TH2F("CombIsoVsRho25DR03_Pt50toINF"," Combined Isolation of barrel photons VS Rho after Et, h/e, r9<1.0, CaloIdL, IsoVL||R9>.8 cuts;Rho (Et/ unit area);ecalRecHitSumEtConeDR03+hcalTowerSumEtConeDR03+trkSumPtHollowConeDR03",50,0,50,6250,-50,1200);h_CombIsoVsRho25DR03_Pt50toINF->Sumw2();
  

  bool doSinglePhoton=false;


 
  if(doRhoCorrection)cout<<"Applying Rho Pileup corrections!"<<endl;
  else cout<<"Applying No Pileup corrections!"<<endl;

  // to check duplicated events
  std::map<int, std::set<int> > allEvents;

  // start event looping
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < processNEvents; jentry++) {

    if(printLevel > 1) cout << "Get the tree contents." << endl;

    // if( (int)(((float)jentry/processNEvents)*100)%10==0 )cout << (int)(((float)jentry/processNEvents)*100) << " percent done..." << endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    if(printLevel > 1 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0)) ) {
      cout <</* endl <<*/ int(jentry) << " events processed with Run:"
	   << event->runNumber << ", Event:" << event->eventNumber << endl;
    }


    if(printLevel > 1) cout << "Initialize any global variables to be reset per event." << endl;

    InitializePerEvent();

    nCnt[0]++; // total number of events    

    if(printLevel > 1) cout << "Apply good run list." << endl;
    // uncomment this to use the Json file to flag good data (or bad depending on your outlook)  
    if(printLevel > 1)cout<<"runNumber="<<event->runNumber<<"  lumiNumber="<<event->luminosityBlockNumber<<endl;
    if(useJSON){
      if(!isInJson(event->runNumber,event->luminosityBlockNumber)) continue;
    }
    // uncomment this to print all ntuple variables
    //Print(*event);
    
    if(printLevel > 1) cout << "Check duplicated events for data only." << endl;
    if(event->isRealData){
      bool duplicateEvent = ! (allEvents[event->runNumber].insert(event->eventNumber)).second;
      if(duplicateEvent){
	cout<<"Duplicate Event! Run "<<event->runNumber<<" Event "<<event->eventNumber<<endl;
	continue;
      }
    }
    
    nCnt[1]++; // total number of events passing Json 

    if(printLevel > 1) cout << "Apply trigger selection in the event." << endl;
    bool passHLT = (useTrigger ? PassTriggers() : true);
    if(!passHLT )continue;//only accept events which pass our trigger(s)
    nCnt[2]++;// number of events that pass trigger
   
    //  Get NVertex and Rho for event
    int NVertex=0;
    for(std::vector<susy::Vertex>::iterator Vtx_it = event->vertices.begin(); Vtx_it<event->vertices.end(); Vtx_it++){
      if(    !Vtx_it->isFake() 
	     && Vtx_it->ndof>4 
	     && fabs(Vtx_it->position.z()<24.0) 
	     && sqrt(Vtx_it->position.x()*Vtx_it->position.x()+Vtx_it->position.y()*Vtx_it->position.y())<2.0 ) NVertex++;
    }
    float Rho = event->rho;float RhoB = event->rhoBarrel;float Rho25 = event->rho25;
    //cout<<"RhoB="<<RhoB<<endl;
   
    //Require at least 1 good vertex
    if(NVertex<1){
      cout<<"No Good Vertex!!!!  Run: "<<event->runNumber<<"  Event: "<<event->eventNumber<<endl;
      continue;
    }
    nCnt[3]++;// number of events that pass NVertex

    h_rho->Fill(Rho);h_rhoB->Fill(RhoB); h_rho25->Fill(Rho25);
    h_NvertexVsRho25->Fill(Rho25,NVertex);
    h_NvertexVsRhoB->Fill(RhoB,NVertex);

    if(printLevel > 0) cout << "Find loose and tight photons in the event." << endl;

    std::map<TString, std::vector<susy::Photon> >::iterator phoMap = event->photons.find("photons");
    if(phoMap != event->photons.end()) {
    
      //loop over photon collection 
      for(std::vector<susy::Photon>::iterator it_Pho = phoMap->second.begin(); it_Pho != phoMap->second.end(); it_Pho++) {
	if(printLevel > 0) cout<<"looping over photon collection"<<endl;

	//grab photon quantities
	float ecalIsoDR03 = it_Pho->ecalRecHitSumEtConeDR03;
	float hcalIsoDR03 = it_Pho->hcalTowerSumEtConeDR03();
	float trackIsoDR03 = it_Pho->trkSumPtHollowConeDR03;
	float combIsoDR03 = ecalIsoDR03 + hcalIsoDR03 + trackIsoDR03;

	float ecalIsoDR04 = it_Pho->ecalRecHitSumEtConeDR04;
	float hcalIsoDR04 = it_Pho->hcalTowerSumEtConeDR04();
	float trackIsoDR04 = it_Pho->trkSumPtHollowConeDR04;
	float combIsoDR04 = ecalIsoDR04 + hcalIsoDR04 + trackIsoDR04;

	float chargedHadronIso=it_Pho->chargedHadronIso;
	float neutralHadronIso=it_Pho->neutralHadronIso;
	float photonIso=it_Pho->photonIso;
	float PfCombIso = chargedHadronIso + neutralHadronIso + photonIso;

	float eta = std::abs(it_Pho->caloPosition.Eta());
	float hadOverEM = it_Pho->hadronicOverEm;
	float Et = it_Pho->momentum.Et();
	float R9 = it_Pho->r9;
	float sIetaIeta = it_Pho->sigmaIetaIeta;
	float sIphiIphi = it_Pho->sigmaIphiIphi;
	float Pt = it_Pho->momentum.Pt();
	float nPixel = it_Pho->nPixelSeeds;



	//----------------set up cuts-------------------
	// fiducial cuts. Look for only barrel now
	if(eta > susy::etaGapBegin && !doSinglePhoton) continue;
	nCnt[4]++;// number of photons that pass barrel eta cut

	// Spike cleaning
	//h_r9->Fill(it_Pho->r9);
	bool isSpike = (R9 > 1.0 || sIetaIeta<0.001 || sIphiIphi<0.001);
	if(isSpike && !doSinglePhoton) continue;
	nCnt[5]++;// number of photons that pass spike cut
	//h_SeedTime_afterR9->Fill(it_Pho->seedTime);
	//h_SeedTimeVSE_afterR9->Fill(it_Pho->seedTime,it_Pho->momentum.E());

	// Et cuts, 25 GeV for trailing photons. Will apply tighter for the leading one.
	if(Et < 25.0 && !doSinglePhoton) continue;
	nCnt[6]++;// number of photons that pass momentum cut
	// cuts containing EE cases for later use, but EE is not used for the time being.

	// H/E (in trigger, 0.15 for EB, 0.10 for EE)
	bool heCut = hadOverEM < 0.05 ;
	//if(!heCut)continue;
	//nCnt[7]++;// number of photons that pass  heCut

	//hlt cuts
	bool ecalIsoVLcut = ( ecalIsoDR03 < 6.0 + 0.012*Et );
	bool hcalIsoVLcut = ( hcalIsoDR03 < 4.0 + 0.005*Et );
	bool trackIsoVLcut = ( trackIsoDR03  < 4.0 + 0.002*Et   );
	bool IsoVLCut = ( ecalIsoVLcut && hcalIsoVLcut && trackIsoVLcut);
	bool CaloIdLCut = ( (hadOverEM < 0.15) && (sIetaIeta < 0.014) );
	bool R9TriggerCut = R9 > 0.8;

	bool combIsoCut = (ecalIsoDR03 + hcalIsoDR03 + trackIsoDR03)<6.0;
	bool EleTagCut = (sIetaIeta<0.011) && (nPixel>0) && heCut && combIsoCut;

	bool SinglePhoCut = ( Et>80 &&
			      eta <1.4 &&
			      nPixel==0 &&
			      R9<1.0 &&
			      //combIsoDR03<6.0 &&
			      sIetaIeta<0.011 &&
			      hadOverEM<0.05 &&
			      ( PassTrigger("HLT_Photon70_CaloIdL_HT300_v") ||
				PassTrigger("HLT_Photon70_CaloIdL_HT400_v") ||
				PassTrigger("HLT_Photon70_CaloIdXL_HT400_v")
				)
			      );
	  
	//if(!CaloIdLCut && !doSinglePhoton)continue;
	nCnt[7]++;// number of photons that pass CaloIdL cut

	//if(!IsoVLCut && !R9TriggerCut && !doSinglePhoton)continue;
	nCnt[8]++;// number of photons that pass IsoVL or R9 cut

	//fill vertex and rho h/e plots
	//h_AvgHOverEVsNVertex->Fill(NVertex,hadOverEM);
	//h_AvgHOverEVsRho25->Fill(Rho,hadOverEM);
	h_HOverEVsRho25->Fill(Rho25,hadOverEM);
	h_HOverEVsRhoB->Fill(RhoB,hadOverEM);
	h_HOverEVsNVertex->Fill(NVertex,hadOverEM);
	/*
	  if(AtLeastOnePFJet){
	  h_AvgHOverEVsNVertex_JetReq->Fill(NVertex,hadOverEM);
	  h_AvgHOverEVsRho25_JetReq->Fill(Rho,hadOverEM);
	  }
	*/
	// Fill the rest of the plots
	if(!heCut && !doSinglePhoton)continue;
	nCnt[9]++;// number of photons that pass  heCut
	/*  h_Rho_VS_NVertex_afterHE->Fill(NVertex,Rho);
	    h_Rho_VS_NVertex_Prof_afterHE->Fill(NVertex,Rho);*/
	//DR03 cone study
	/*	  h_EcalIsoDR03->Fill(ecalIsoDR03);
	  h_HcalIsoDR03->Fill(hcalIsoDR03);
	  h_TrackIsoDR03->Fill(trackIsoDR03);*/
	//pileup dependence
	/*
	  h_AvgEcalIsoVsNVertexDR04->Fill(NVertex,ecalIsoDR04);
	  h_AvgHcalIsoVsNVertexDR04->Fill(NVertex,hcalIsoDR04);
	  h_AvgTrackIsoVsNVertexDR04->Fill(NVertex,trackIsoDR04);
	  h_AvgEcalIsoVsRho25DR04->Fill(Rho,ecalIsoDR04);
	  h_AvgHcalIsoVsRho25DR04->Fill(Rho,hcalIsoDR04);
	  h_AvgTrackIsoVsRho25DR04->Fill(Rho,trackIsoDR04);
	  h_AvgEcalIsoVsRho25DR03->Fill(Rho,ecalIsoDR03);
	  h_AvgHcalIsoVsRho25DR03->Fill(Rho,hcalIsoDR03);
	  h_AvgTrackIsoVsRho25DR03->Fill(Rho,trackIsoDR03);
	  h_AvgCombIsoVsRho25DR03->Fill(Rho,combIsoDR03);
	  h_AvgCombIsoVsRho25DR04->Fill(Rho,combIsoDR04);*/

	if(doSinglePhoton && SinglePhoCut){
	  h_EcalIsoVsRho25DR03_SinglePho->Fill(Rho25,ecalIsoDR03);
	  h_HcalIsoVsRho25DR03_SinglePho->Fill(Rho25,hcalIsoDR03);
	  h_TrackIsoVsRho25DR03_SinglePho->Fill(Rho25,trackIsoDR03);
	  h_CombIsoVsRho25DR03_SinglePho->Fill(Rho25,combIsoDR03);
	  h_EcalIsoVsRhoBDR03_SinglePho->Fill(RhoB,ecalIsoDR03);
	  h_HcalIsoVsRhoBDR03_SinglePho->Fill(RhoB,hcalIsoDR03);
	  h_TrackIsoVsRhoBDR03_SinglePho->Fill(RhoB,trackIsoDR03);
	  h_CombIsoVsRhoBDR03_SinglePho->Fill(RhoB,combIsoDR03);
	}
	
	h_EcalIsoVsRho25DR03->Fill(Rho25,ecalIsoDR03);
	h_HcalIsoVsRho25DR03->Fill(Rho25,hcalIsoDR03);
	h_TrackIsoVsRho25DR03->Fill(Rho25,trackIsoDR03);
	h_CombIsoVsRho25DR03->Fill(Rho25,combIsoDR03);
	
	h_EcalIsoVsRho25DR04->Fill(Rho25,ecalIsoDR04);
	h_HcalIsoVsRho25DR04->Fill(Rho25,hcalIsoDR04);
	h_TrackIsoVsRho25DR04->Fill(Rho25,trackIsoDR04);
	h_CombIsoVsRho25DR04->Fill(Rho25,combIsoDR04);

	h_EcalIsoVsRhoBDR03->Fill(RhoB,ecalIsoDR03);
	h_HcalIsoVsRhoBDR03->Fill(RhoB,hcalIsoDR03);
	h_TrackIsoVsRhoBDR03->Fill(RhoB,trackIsoDR03);
	h_CombIsoVsRhoBDR03->Fill(RhoB,combIsoDR03);
	
	h_EcalIsoVsRhoBDR04->Fill(RhoB,ecalIsoDR04);
	h_HcalIsoVsRhoBDR04->Fill(RhoB,hcalIsoDR04);
	h_TrackIsoVsRhoBDR04->Fill(RhoB,trackIsoDR04);
	h_CombIsoVsRhoBDR04->Fill(RhoB,combIsoDR04);

	h_EcalIsoVsNVertexDR03->Fill(NVertex,ecalIsoDR03);
	h_HcalIsoVsNVertexDR03->Fill(NVertex,hcalIsoDR03);
	h_TrackIsoVsNVertexDR03->Fill(NVertex,trackIsoDR03);
	h_CombIsoVsNVertexDR03->Fill(NVertex,combIsoDR03);
	
	h_EcalIsoVsNVertexDR04->Fill(NVertex,ecalIsoDR04);
	h_HcalIsoVsNVertexDR04->Fill(NVertex,hcalIsoDR04);
	h_TrackIsoVsNVertexDR04->Fill(NVertex,trackIsoDR04);
	h_CombIsoVsNVertexDR04->Fill(NVertex,combIsoDR04);

	h_HcalIsoVsEcalIsoDR03->Fill(ecalIsoDR03,hcalIsoDR03);
	h_TrackIsoVsEcalIsoDR03->Fill(ecalIsoDR03,trackIsoDR03);
	h_TrackIsoVsHcalIsoDR03->Fill(hcalIsoDR03,trackIsoDR03);
	
	if(NVertex<=5){
	  h_EcalIsoVsRho25DR03_lt6NV->Fill(Rho25,ecalIsoDR03);
	  h_HcalIsoVsRho25DR03_lt6NV->Fill(Rho25,hcalIsoDR03);
	}
	if(NVertex>5  && NVertex<=10){
	  h_EcalIsoVsRho25DR03_6to10NV->Fill(Rho25,ecalIsoDR03);
	  h_HcalIsoVsRho25DR03_6to10NV->Fill(Rho25,hcalIsoDR03);
	}
	if(NVertex>10  && NVertex <=15){
	  h_EcalIsoVsRho25DR03_11to15NV->Fill(Rho25,ecalIsoDR03);
	  h_HcalIsoVsRho25DR03_11to15NV->Fill(Rho25,hcalIsoDR03);
	} 
	if(NVertex>15){
	  h_EcalIsoVsRho25DR03_gt15NV->Fill(Rho25,ecalIsoDR03);
	  h_HcalIsoVsRho25DR03_gt15NV->Fill(Rho25,hcalIsoDR03);
	}

	if(Pt<40){
	  h_EcalIsoVsRho25DR03_Pt0to40->Fill(Rho25,ecalIsoDR03);
	  h_HcalIsoVsRho25DR03_Pt0to40->Fill(Rho25,hcalIsoDR03);
	  h_TrackIsoVsRho25DR03_Pt0to40->Fill(Rho25,trackIsoDR03);
	  h_CombIsoVsRho25DR03_Pt0to40->Fill(Rho25,combIsoDR03);
	}
	if(Pt>=40 && Pt<50){
	  h_EcalIsoVsRho25DR03_Pt40to50->Fill(Rho25,ecalIsoDR03);
	  h_HcalIsoVsRho25DR03_Pt40to50->Fill(Rho25,hcalIsoDR03);
	  h_TrackIsoVsRho25DR03_Pt40to50->Fill(Rho25,trackIsoDR03);
	  h_CombIsoVsRho25DR03_Pt40to50->Fill(Rho25,combIsoDR03);
	}
	if(Pt>=50){
	  h_EcalIsoVsRho25DR03_Pt50toINF->Fill(Rho25,ecalIsoDR03);
	  h_HcalIsoVsRho25DR03_Pt50toINF->Fill(Rho25,hcalIsoDR03);
	  h_TrackIsoVsRho25DR03_Pt50toINF->Fill(Rho25,trackIsoDR03);
	  h_CombIsoVsRho25DR03_Pt50toINF->Fill(Rho25,combIsoDR03);
	}
	  
	ne=0;
	if(EleTagCut){
	  for(std::vector<susy::Photon>::iterator it_Ele = phoMap->second.begin(); it_Ele != phoMap->second.end(); it_Ele++) {
	    if(!isSameObject(it_Pho->caloPosition,it_Ele->caloPosition,0.1)){
	      if(InvariantMass((&*it_Pho)->momentum,(&*it_Ele)->momentum)>81 && InvariantMass((&*it_Pho)->momentum,(&*it_Ele)->momentum)<101){
		if(//Make sure other Electron passes base cuts:
		   std::abs(it_Ele->caloPosition.Eta())<susy::etaGapBegin &&
		   it_Ele->r9<1.0 &&
		   it_Ele->sigmaIetaIeta<0.014 &&
		   it_Ele->nPixelSeeds>0 &&
		   it_Ele->momentum.Et()>25/* &&
					      (it_Ele->r9>0.8 || 
					      (it_Ele->ecalRecHitSumEtConeDR03  < 6.0 + 0.012*it_Ele->momentum.Et() &&
					      it_Ele->hcalTowerSumEtConeDR03() < 4.0 + 0.005*it_Ele->momentum.Et() &&
					      it_Ele->trkSumPtHollowConeDR03   < 4.0 + 0.002*it_Ele->momentum.Et() 
					      )
					      )*/
		   ){
		  h_HOverEVsRho25_ee->Fill(Rho25,it_Ele->hadronicOverEm);
		  if(it_Ele->ecalRecHitSumEtConeDR03+it_Ele->hcalTowerSumEtConeDR03()+it_Ele->trkSumPtHollowConeDR03-(PUCorr_ECAL+PUCorr_HCAL+PUCorr_TRACK)*Rho25<6.0)h_HOverEAfterIsoVsRho25_ee->Fill(Rho25,it_Ele->hadronicOverEm);
		  if(it_Ele->hadronicOverEm<0.05){
		    ne++;
		    h_EcalIsoVsRho25DR03_ee->Fill(Rho25,it_Ele->ecalRecHitSumEtConeDR03);
		    h_HcalIsoVsRho25DR03_ee->Fill(Rho25,it_Ele->hcalTowerSumEtConeDR03());
		    h_TrackIsoVsRho25DR03_ee->Fill(Rho25,it_Ele->trkSumPtHollowConeDR03);
		    h_CombIsoVsRho25DR03_ee->Fill(Rho25,it_Ele->ecalRecHitSumEtConeDR03+it_Ele->hcalTowerSumEtConeDR03()+it_Ele->trkSumPtHollowConeDR03);
		    h_EcalIsoVsRhoBDR03_ee->Fill(RhoB,it_Ele->ecalRecHitSumEtConeDR03);
		    h_HcalIsoVsRhoBDR03_ee->Fill(RhoB,it_Ele->hcalTowerSumEtConeDR03());
		    h_TrackIsoVsRhoBDR03_ee->Fill(RhoB,it_Ele->trkSumPtHollowConeDR03);
		    h_CombIsoVsRhoBDR03_ee->Fill(RhoB,it_Ele->ecalRecHitSumEtConeDR03+it_Ele->hcalTowerSumEtConeDR03()+it_Ele->trkSumPtHollowConeDR03);
		    h_EcalIsoVsRho25DR04_ee->Fill(Rho25,it_Ele->ecalRecHitSumEtConeDR04);
		    h_HcalIsoVsRho25DR04_ee->Fill(Rho25,it_Ele->hcalTowerSumEtConeDR04());
		    h_TrackIsoVsRho25DR04_ee->Fill(Rho25,it_Ele->trkSumPtHollowConeDR04);
		    h_CombIsoVsRho25DR04_ee->Fill(Rho25,it_Ele->ecalRecHitSumEtConeDR04+it_Ele->hcalTowerSumEtConeDR04()+it_Ele->trkSumPtHollowConeDR04);
		    h_chargedHadronIsoVsRho25_ee->Fill(Rho25,it_Ele->chargedHadronIso);
		    h_neutralHadronIsoVsRho25_ee->Fill(Rho25,it_Ele->neutralHadronIso);
		    h_photonIsoVsRho25_ee->Fill(Rho25,it_Ele->photonIso);
		    h_PfCombIsoVsRho25_ee->Fill(Rho25,(it_Ele->chargedHadronIso+it_Ele->neutralHadronIso+it_Ele->photonIso));
		  }
		}
	      }
	    }
	  }   
	}
	if(ne>0)n2e++;
	if(ne>1)n3e++;

	if(printLevel > 0) cout<<"End of Photon Loop"<<endl;
	  
      }//for(it_Pho)
      
    }//if phoMap!=end (after phoMap def)
  
  } // for jentry
  
  
  cout << "---ALL DONE!------"<<endl;
  cout << " ----------------- Job Summary ----------------- " << endl;
  cout << " Total events               : " << nCnt[0] << " (" << 100*nCnt[0]/float(nCnt[0]) << "%)" << endl;
  cout << " JSON events passed         : " << nCnt[1] << " (" << 100*nCnt[1]/float(nCnt[0]) << "%)" << endl;
  cout << " HLT events passed          : " << nCnt[2] << " (" << 100*nCnt[2]/float(nCnt[0]) << "%)" << endl;
  cout << " nVertex events passed      : " << nCnt[3] << " (" << 100*nCnt[3]/float(nCnt[0]) << "%)" << endl;
  cout << " Barrel eta photons passed  : " << nCnt[4] /*<< " (" << 100*nCnt[4]/float(nCnt[0]) << "%)"*/ << endl;
  cout << " !isSpike photons passed    : " << nCnt[5] << " (" << 100*nCnt[5]/float(nCnt[4]) << "%)" << endl;
  cout << " Et>25 photons passed       : " << nCnt[6] << " (" << 100*nCnt[6]/float(nCnt[4]) << "%)" << endl;
  //cout << " CaloIdLCut photons passed  : " << nCnt[7] << " (" << 100*nCnt[7]/float(nCnt[4]) << "%)" << endl;
  //cout << " IsoVL or R9 photons passed : " << nCnt[8] << " (" << 100*nCnt[8]/float(nCnt[4]) << "%)" << endl;
  cout << " heCut photons passed       : " << nCnt[9] << " (" << 100*nCnt[9]/float(nCnt[4]) << "%)" << endl;
 
  cout << " ----------------------------------------------- " << endl;

  cout<<"Writing root output to: /data/ndpc3/b/dmorse/RA3/AnalysisOutput/Pileup/cms525v3_jec2012/hist_"<<ds<<".root"<<endl;
  //cout<<"Writing root output to: /Users/dmorse/RA3/AnalysisOutput/hist_"<<ds<<".root"<<endl;

  // close the output file

  fout->cd();
  fout->Write();
  fout->Close();
  delete fout;

}
    

void SusyEventAnalyzer::PhotonId() {
  
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  cout << "total events in files  : " << nentries << endl;

  if(processNEvents <= 0 || processNEvents > nentries) processNEvents = nentries;

  cout << "events to be processed : " << processNEvents << endl; 
  


  if(printLevel > 1) cout << "Initialize event counters." << endl;
  const int NCNT = 20;
  int nCnt[NCNT];
  for(int i=0; i<NCNT; i++) nCnt[i] = 0;

  // open hist file and define histograms

  TFile* fout = new TFile("/data/ndpc3/b/dmorse/RA3/AnalysisOutput/PhotonId/hist_"+ds+".root","RECREATE");

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->Divide(2);
  c1->cd(1)->SetPad(0.,.36,1.,1.);
  c1->cd(2)->SetPad(0.,0.,1.,.36);
  fout->cd();
  int NumBins=38;
  double EtBins[39]={0.};
  for(int i=0;i<39;i++){
    if(i<=20)EtBins[i]=i*10;//0-200
    if(i>20 && i<=30)EtBins[i]=i*20-200;//200-400
    if(i>30 && i<=33)EtBins[i]=i*50-1100;//400-600
    if(i>33)EtBins[i]=i*100-2800;//600-1000
  }

  TH1F* h_Pho_Et = new TH1F("Pho_Et","E_{T} of all reco::photon's;E_{T} (GeV);Events",100,0.,1000.);  

  TH1F* h_Pho_Et_Truth = new TH1F("Pho_Et_Truth","E_{T} of reco::photon dR matched to genPhoton;E_{T} (GeV);Events",100,0.,1000.);
  TH1F* h_Pho_Eta_Truth = new TH1F("Pho_Eta_Truth","#eta of reco::photon dR matched to genPhoton;#eta;Events",20,-2.,2.);
  TH1F* h_Pho_Phi_Truth = new TH1F("Pho_Phi_Truth","#phi of reco::photon dR matched to genPhoton;#phi;Events",20,-4.,4.);
  TH1I* h_Pho_Nvertex_Truth = new TH1I("Pho_Nvertex_Truth","NVertex of reco::photon dR matched to genPhoton;NVertex;Events",50,0,50);

  TH1F* h_Pho_Et_Truth_PF = new TH1F("Pho_Et_Truth_PF","E_{T} of reco::photon dR matched to genPhoton;E_{T} (GeV);Events",100,0.,1000.);
  TH1F* h_Pho_Eta_Truth_PF = new TH1F("Pho_Eta_Truth_PF","#eta of reco::photon dR matched to genPhoton;#eta;Events",20,-2.,2.);
  TH1F* h_Pho_Phi_Truth_PF = new TH1F("Pho_Phi_Truth_PF","#phi of reco::photon dR matched to genPhoton;#phi;Events",20,-4.,4.);
  TH1I* h_Pho_Nvertex_Truth_PF = new TH1I("Pho_Nvertex_Truth_PF","NVertex of reco::photon dR matched to genPhoton;NVertex;Events",50,0,50);
  TH1I* h_Pho_Nvertex_Truth_PF_BG = new TH1I("Pho_Nvertex_Truth_PF_BG","NVertex of reco::photon dR matched to genPhoton;NVertex;Events",50,0,50);

  TH1F* h_Pho_Et_Truth_BG = new TH1F("Pho_Et_Truth_BG","E_{T} of reco::photon dR matched to genPhoton;E_{T} (GeV);Events",100,0.,1000.);
  TH1F* h_Pho_Eta_Truth_BG = new TH1F("Pho_Eta_Truth_BG","#eta of reco::photon dR matched to genPhoton;#eta;Events",20,-2.,2.);
  TH1F* h_Pho_Phi_Truth_BG = new TH1F("Pho_Phi_Truth_BG","#phi of reco::photon dR matched to genPhoton;#phi;Events",20,-4.,4.);
  TH1I* h_Pho_Nvertex_Truth_BG = new TH1I("Pho_Nvertex_Truth_BG","NVertex of reco::photon dR matched to genPhoton;NVertex;Events",50,0,50);

  TH1F* h_Pho_Et_RA3 = new TH1F("Pho_Et_RA3","E_{T} of reco::photon passing RA3 selection;E_{T} (GeV);Events",100,0.,1000.);
  TH1F* h_Pho_Eta_RA3 = new TH1F("Pho_Eta_RA3","#eta of reco::photon passing RA3 selection;#eta;Events",20,-2.,2.);
  TH1F* h_Pho_Phi_RA3 = new TH1F("Pho_Phi_RA3","#phi of reco::photon passing RA3 selection;#phi;Events",20,-4.,4.);
  TH1I* h_Pho_Nvertex_RA3 = new TH1I("Pho_Nvertex_RA3","NVertex of reco::photon passing RA3 selection;NVertex;Events",50,0,50);

  TH1I* h_Pho_Nvertex_RA3_NoPUcorr = new TH1I("Pho_Nvertex_RA3_NoPUcorr","NVertex of reco::photon passing RA3 selection with No Rho correction;NVertex;Events",50,0,50);
  TH1I* h_Pho_Nvertex_Egamma = new TH1I("Pho_Nvertex_Egamma","NVertex of reco::photon passing Egamma selection;NVertex;Events",50,0,50);
  TH1I* h_Pho_Nvertex_Egamma_NoPUcorr = new TH1I("Pho_Nvertex_Egamma_NoPUcorr","NVertex of reco::photon passing Egamma selection with No Rho correction;NVertex;Events",50,0,50);

  TH1I* h_Pho_Nvertex_Egamma_NoTrackIso = new TH1I("Pho_Nvertex_Egamma_NoTrackIso","NVertex of reco::photon passing Egamma selection except track isolation;NVertex;Events",50,0,50);
  TH1I* h_Pho_Nvertex_Egamma_NoTrackIso_NoPUcorr = new TH1I("Pho_Nvertex_Egamma_NoTrackIso_NoPUcorr","NVertex of reco::photon passing Egamma selection except track isolation with No Rho correction;NVertex;Events",50,0,50);
  TH1I* h_Pho_Nvertex_Egamma_NoTrackIsoNoEcalIso = new TH1I("Pho_Nvertex_Egamma_NoTrackIsoNoEcalIso","NVertex of reco::photon passing Egamma selection except track or hcal isolation;NVertex;Events",50,0,50);
  TH1I* h_Pho_Nvertex_Egamma_NoTrackIsoNoEcalIso_NoPUcorr = new TH1I("Pho_Nvertex_Egamma_NoTrackIsoNoEcalIso_NoPUcorr","NVertex of reco::photon passing Egamma selection except track or hcal isolation with No Rho correction;NVertex;Events",50,0,50);
  TH1I* h_Pho_Nvertex_Egamma_NoIso = new TH1I("Pho_Nvertex_Egamma_NoIso","NVertex of reco::photon passing Egamma selection except isolation;NVertex;Events",50,0,50);
 



  TH1F* h_Pho_Et_RA3_BG = new TH1F("Pho_Et_RA3_BG","E_{T} of reco::photon passing RA3 selection;E_{T} (GeV);Events",100,0.,1000.);
  TH1F* h_Pho_Eta_RA3_BG = new TH1F("Pho_Eta_RA3_BG","#eta of reco::photon passing RA3 selection;#eta;Events",20,-2.,2.);
  TH1F* h_Pho_Phi_RA3_BG = new TH1F("Pho_Phi_RA3_BG","#phi of reco::photon passing RA3 selection;#phi;Events",20,-4.,4.);
  TH1I* h_Pho_Nvertex_RA3_BG = new TH1I("Pho_Nvertex_RA3_BG","NVertex of reco::photon passing RA3 selection;NVertex;Events",50,0,50);

  TH1I* h_Pho_Nvertex_RA3_BG_NoPUcorr = new TH1I("Pho_Nvertex_RA3_BG_NoPUcorr","NVertex of reco::photon passing RA3 selection with No Rho correction;NVertex;Events",50,0,50);
  
  TH1F* h_Pho_Et_PF = new TH1F("Pho_Et_PF","E_{T} of reco::photon passing PF selection;E_{T} (GeV);Events",100,0.,1000.);
  TH1F* h_Pho_Eta_PF = new TH1F("Pho_Eta_PF","#eta of reco::photon passing PF selection;#eta;Events",20,-2.,2.);
  TH1F* h_Pho_Phi_PF = new TH1F("Pho_Phi_PF","#phi of reco::photon passing PF selection;#phi;Events",20,-4.,4.);
  TH1I* h_Pho_Nvertex_PF = new TH1I("Pho_Nvertex_PF","NVertex of reco::photon passing PF selection;NVertex;Events",50,0,50);
  TH1I* h_Pho_Nvertex_PF_NoPUcorr = new TH1I("Pho_Nvertex_PF_NoPUcorr","NVertex of reco::photon passing PF selection;NVertex;Events",50,0,50);

  TH1F* h_Pho_Et_PF_BG = new TH1F("Pho_Et_PF_BG","E_{T} of reco::photon passing PF selection;E_{T} (GeV);Events",100,0.,1000.);
  TH1F* h_Pho_Eta_PF_BG = new TH1F("Pho_Eta_PF_BG","#eta of reco::photon passing PF selection;#eta;Events",20,-2.,2.);
  TH1F* h_Pho_Phi_PF_BG = new TH1F("Pho_Phi_PF_BG","#phi of reco::photon passing PF selection;#phi;Events",20,-4.,4.);
  TH1I* h_Pho_Nvertex_PF_BG = new TH1I("Pho_Nvertex_PF_BG","NVertex of reco::photon passing PF selection;NVertex;Events",50,0,50);
  TH1I* h_numChi0 = new TH1I("numChi0","Number of neutralinos in event",10,0,10);
  TH1I* h_numPhoFromChi0 = new TH1I("numPhoFromChi0","Number of photons from neutralinos in event",10,0,10);

  int numPhos=0,numPfPhos=0,numPhosRA3=0,numPFPhosRA3=0,numPhosTruth=0,numPhosVgamma=0;
  int numPho=0,numPhoMatch=0;
  
  if(doRhoCorrection)cout<<"Applying Rho Pileup corrections!"<<endl;
  else cout<<"Applying NO Pileup corrections!"<<endl;
 
  // start event looping
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < processNEvents; jentry++) {

    if(doNVertexCorrection && doRhoCorrection){ cout<<"Trying to apply both Rho AND Nvertex Pileup corrections!"<<endl; break;}
    
    if(printLevel > 1) cout << "Get the tree contents." << endl;
    
    // if( (int)(((float)jentry/processNEvents)*100)%10==0 )cout << (int)(((float)jentry/processNEvents)*100) << " percent done..." << endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    if(printLevel > 1 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0)) ) {
      cout <</* endl <<*/ int(jentry) << " events processed with Run:"
	   << event->runNumber << ", Event:" << event->eventNumber << endl;
    }


    if(printLevel > 1) cout << "Initialize any global variables to be reset per event." << endl;

    InitializePerEvent();

    nCnt[0]++; // total number of events    

    if(printLevel > 1) cout << "Apply good run list." << endl;
    // uncomment this to use the Json file to flag good data (or bad depending on your outlook)  
    if(printLevel > 1)cout<<"runNumber="<<event->runNumber<<"  lumiNumber="<<event->luminosityBlockNumber<<endl;
    if(event->isRealData && useJSON){
      if(!isInJson(event->runNumber,event->luminosityBlockNumber)) continue;
    }
    // uncomment this to print all ntuple variables
    //Print(*event);
    
    nCnt[1]++; // total number of events passing Json 

    if(printLevel > 1) cout << "Apply trigger selection in the event." << endl;
    bool passHLT = (useTrigger ? PassTriggers() : true);
    if(!passHLT )continue;//only accept events which pass our trigger(s)
    nCnt[2]++;// number of events that pass trigger
   
    //  Get NVertex and Rho for event
    int NVertex=0;
    for(std::vector<susy::Vertex>::iterator Vtx_it = event->vertices.begin(); Vtx_it<event->vertices.end(); Vtx_it++){
      if(    !Vtx_it->isFake() 
	     && Vtx_it->ndof>4 
	     && fabs(Vtx_it->position.z()<24.0) 
	     && sqrt(Vtx_it->position.x()*Vtx_it->position.x()+Vtx_it->position.y()*Vtx_it->position.y())<2.0 ) NVertex++;
    }
    float Rho = event->rho;float RhoB = event->rhoBarrel;float Rho25 = event->rho25;
   
    //Require at least 1 good vertex
    if(NVertex<1){
      cout<<"No Good Vertex!!!!  Run: "<<event->runNumber<<"  Event: "<<event->eventNumber<<endl;
      continue;
    }
    nCnt[3]++;// number of events that pass NVertex
  
    if(!event->isRealData &&jentry==0){
      cout<<"Neutralino Mass: " <<event->gridParams["mChi0"]<<endl;
      cout<<"Gluino Mass: " <<event->gridParams["mGluino"]<<endl;
      cout<<"Squark Mass: " <<event->gridParams["mSquark"]<<endl;
      cout<<"Cross Section: " <<event->gridParams["xsec"]<<endl;
      //cout<<"ptHat: " <<event->gridParams["ptHat"]<<endl;
    }
    
    //LumiWeight = (x-sec)(IntegratedLumi)/(# events)
    //float LumiWeight = 1;
    //float ptHatWeight=1./pow((event->gridParams["ptHat"]/15.),4.5);

    
    //----------
    std::vector<susy::Photon*>   pho_Cands_RA3_NoPUcorr;
    std::vector<susy::Photon*>   pho_Cands_RA3;
    std::vector<susy::Photon*>   pho_Cands_Egamma_NoPUcorr;
    std::vector<susy::Photon*>   pho_Cands_Egamma;
    std::vector<susy::Photon*>   pho_Cands_PF;
    std::vector<susy::Photon*>   pho_Cands_PF_NoPUcorr;
    std::vector<susy::Photon*>   pho_Cands_Truth;
    std::vector<susy::Photon*>   pho_Cands_Truth_PF;    
    std::vector<susy::Photon*>   pho_Cands_Truth_PF_NoPUcorr;
    std::vector<susy::Photon*>   pho_Cands_RA3_BG_NoPUcorr;
    std::vector<susy::Photon*>   pho_Cands_RA3_BG;
    std::vector<susy::Photon*>   pho_Cands_PF_BG;
    std::vector<susy::Photon*>   pho_Cands_Truth_BG;
    std::vector<susy::Photon*>   pho_Cands_Truth_PF_BG;
    susy::Particle* genPart = new susy::Particle; genPart->Init();
    susy::Particle* genPartBG = new susy::Particle; genPartBG->Init();

    int Chi0=0,PhoFromChi0=0;
    for(std::vector<susy::Particle>::iterator Part_it = event->genParticles.begin(); 
	Part_it != event->genParticles.end(); Part_it++){
      if(Part_it->pdgId==1000022)Chi0++;
      if(Part_it->pdgId==22/* && Part_it->motherId==1000022*/)PhoFromChi0++;
    }
    h_numChi0->Fill(Chi0);
    h_numPhoFromChi0->Fill(PhoFromChi0);
    
    if(printLevel > 1) cout << "Find loose and tight photons in the event." << endl;
    //find photons, sort by energy
    //----------
    std::map<TString, std::vector<susy::Photon> >::iterator phoMap = event->photons.find("photons");
    if(phoMap != event->photons.end()) {

      //loop over photon collection 
      for(std::vector<susy::Photon>::iterator it_Pho = phoMap->second.begin(); it_Pho != phoMap->second.end(); it_Pho++) {
	if(printLevel > 1) cout<<"looping over photon collection"<<endl;
	numPhos++;

	h_Pho_Et->Fill(it_Pho->momentum.Et());
	//if(it_Pho->chargedHadronIso!=0 || it_Pho->neutralHadronIso!=0 || it_Pho->photonIso!=0)cout<<"chargedHadronIso: "<<it_Pho->chargedHadronIso<<"  neutralHadronIso:"<<it_Pho->neutralHadronIso<<"  photonIso: "<<it_Pho->photonIso<<endl<<endl;
	//if find genMatch, put truth pho into truth histo
	bool hasMatch=false,hasMatchBG=false;
	for(std::vector<susy::Particle>::iterator Part_it = event->genParticles.begin(); 
	    Part_it != event->genParticles.end(); Part_it++){
	  if(Part_it->pdgId==22/* && Part_it->motherId==1000022*/ && isSameObject(it_Pho->caloPosition,Part_it->momentum,0.3)){
	    numPho++;
	    if(isSameObject(it_Pho->caloPosition,Part_it->momentum,0.1)){
	      numPhoMatch++;
	    }
	  }
	  if(Part_it->pdgId==22/* && Part_it->motherId==1000022*/){
	    if( isSameObject(it_Pho->caloPosition,Part_it->momentum,0.1) ){
	      genPart=&*Part_it;hasMatch=true;
	    }
	  }
	  if(Part_it->pdgId==22/* && Part_it->motherId==111*/){
	    if( isSameObject(it_Pho->caloPosition,Part_it->momentum,0.1) ){
	      genPartBG=&*Part_it;hasMatchBG=true;
	    }
	  }
	  if(hasMatch || hasMatchBG)break;
	}
	//only look at reco::photons that have truth match
	//if(hasMatch || hasMatchBG){
	 
	//do things that want all photons with no cuts here:
	    
	//----------------set up cuts-------------------
	float ecalIsoDR03=it_Pho->ecalRecHitSumEtConeDR03;
	float hcalIsoDR03=it_Pho->hcalTowerSumEtConeDR03();
	float trackIsoDR03=it_Pho->trkSumPtHollowConeDR03;
	float ecalIsoDR04=it_Pho->ecalRecHitSumEtConeDR04;
	float hcalIsoDR04=it_Pho->hcalTowerSumEtConeDR04();
	float trackIsoDR04=it_Pho->trkSumPtHollowConeDR04;
	float chargedHadronIso=it_Pho->chargedHadronIso;
	float neutralHadronIso=it_Pho->neutralHadronIso;
	float photonIso=it_Pho->photonIso;
	float PhoEt=it_Pho->momentum.Et();
	// fiducial cuts. Look for only barrel now
	bool etaCut = (std::abs(it_Pho->caloPosition.Eta()) < susy::etaGapBegin);
	    
	// Spike cleaning
	bool isSpike = (it_Pho->r9 > 1.0 || it_Pho->sigmaIetaIeta<0.001 || it_Pho->sigmaIphiIphi<0.001);
	//bool isSpikeGen = (genPart->r9 > 1.0 || genPart->sigmaIetaIeta<0.001 || genPart->sigmaIphiIphi<0.001);
	    
	// Et cuts, 25 GeV for trailing photons. Will apply tighter for the leading one.
	bool EtCut = (PhoEt > 25.0);
	    
	// H/E (in trigger, 0.15 for EB, 0.10 for EE)
	bool heCut = (it_Pho->hadronicOverEm < 0.05);
	// sigma_ietaieta (in trigger 0.014 for EB, 0.034 for EE)
	bool sIetaCut = (it_Pho->sigmaIetaIeta < 0.011);

	//extra Vgamma cuts
	bool sIetaLowCut = (it_Pho->sigmaIetaIeta > 0.001);
	bool sIphiLowCut = (it_Pho->sigmaIphiIphi > 0.001);
	bool EcalIsoEgammaCutNoPUcorr = ecalIsoDR04 < 4.2 + 0.006*PhoEt;
	bool HcalIsoEgammaCutNoPUcorr = hcalIsoDR04 < 2.2 + 0.0025*PhoEt;
	bool TrackIsoEgammaCutNoPUcorr = trackIsoDR04 < 2.0 + 0.001*PhoEt;
	bool EcalIsoEgammaCut = ecalIsoDR04 < 4.2 + 0.006*PhoEt + 0.183*Rho25;
	bool HcalIsoEgammaCut = hcalIsoDR04 < 2.2 + 0.0025*PhoEt + 0.062*Rho25;
	bool TrackIsoEgammaCut = trackIsoDR04 < 2.0 + 0.001*PhoEt + 0.0167*Rho25;
	bool EgammaIsoCutNoPUcorr=EcalIsoEgammaCutNoPUcorr && HcalIsoEgammaCutNoPUcorr && TrackIsoEgammaCutNoPUcorr;
	bool EgammaIsoCut=EcalIsoEgammaCut && HcalIsoEgammaCut && TrackIsoEgammaCut;
	//combined cuts
	//  bool combIsoCut =( ( ecalIsoDR03 + hcalIsoDR03 + trackIsoDR03 ) < 6.0 );
	bool combIsoCut_NoPUcorr =( ( ecalIsoDR03 + hcalIsoDR03 + trackIsoDR03 ) < 6.0 );
	bool combIsoCut =( ( ecalIsoDR03 - PUCorr_ECAL*Rho25 + hcalIsoDR03 - PUCorr_HCAL*Rho25 + trackIsoDR03 - PUCorr_TRACK*Rho25 ) < 6.0 );

	//bool PFIsoCut = (chargedHadronIso + neutralHadronIso + photonIso) / PhoEt < 0.1;
	bool PFIsoCut = (chargedHadronIso+neutralHadronIso+photonIso - (PUCorr_chargedHadron+PUCorr_neutralHadron+PUCorr_photon)*Rho25 )<6.0;	
	bool PFIsoCut_NoPUcorr = (chargedHadronIso+neutralHadronIso+photonIso)<6.0;
	//hlt cuts
	bool ecalIsoVLcut = ( ecalIsoDR03  < 6.0 + 0.012*PhoEt );
	bool hcalIsoVLcut = ( hcalIsoDR03 < 4.0 + 0.005*PhoEt );
	bool trackIsoVLcut = ( trackIsoDR03  < 4.0 + 0.002*PhoEt );
	bool IsoVLCut = ( ecalIsoVLcut && hcalIsoVLcut && trackIsoVLcut);
	bool CaloIdLCut = ( (it_Pho->hadronicOverEm < 0.15) && (it_Pho->sigmaIetaIeta < 0.014) );
	//Pixel cut
	bool pixelCut = (it_Pho->nPixelSeeds == 0);
	    
	bool R9Id85=it_Pho->r9>0.85;
	bool CaloId10=( it_Pho->hadronicOverEm<0.1 && it_Pho->sigmaIetaIeta < 0.014 );
	bool Iso50=(ecalIsoDR03<5.0 + 0.012*PhoEt
		    && hcalIsoDR03<5.0 + 0.005*PhoEt
		    && trackIsoDR03<5.0+ 0.002*PhoEt
		    );
	bool R9Id85orCaloId10andIso50Cut = R9Id85 || ( CaloId10 && Iso50 );
	
	bool PhoCut_NoPUcorr = (pixelCut && etaCut && EtCut && !isSpike && heCut && combIsoCut_NoPUcorr && sIetaCut /*&& CaloId10 && Iso50*/);
	bool PhoCut = (pixelCut && etaCut && EtCut && !isSpike && heCut && combIsoCut && sIetaCut /*&& CaloId10 && Iso50*/);
	bool PhoCutEgamma = (pixelCut && etaCut && EtCut && sIetaLowCut && sIphiLowCut && heCut && EgammaIsoCut && sIetaCut /*&& CaloId10 && Iso50*/);
	bool PhoCutEgamma_NoPUcorr = (pixelCut && etaCut && EtCut && sIetaLowCut && sIphiLowCut && heCut && EgammaIsoCutNoPUcorr && sIetaCut /*&& CaloId10 && Iso50*/);
	bool PhoCut_PF = (pixelCut && etaCut && EtCut && !isSpike && heCut && PFIsoCut && sIetaCut/* && CaloId10 && Iso50*/);
	bool PhoCut_PF_NoPUcorr = (pixelCut && etaCut && EtCut && !isSpike && heCut && PFIsoCut_NoPUcorr && sIetaCut/* && CaloId10 && Iso50*/);

	bool passHLT_Pho  = (useTrigger ? ( PassTrigger("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v") ) : true);

	//Now fill pho_Cands
	if(hasMatch){
	  if(1){
	    numPhosTruth++;	 
	    pho_Cands_Truth.push_back(&*it_Pho);
	    h_Pho_Nvertex_Truth->Fill(NVertex);
	  }

	  if( PhoCut_PF && passHLT_Pho){
	    pho_Cands_PF.push_back(&*it_Pho);
	    h_Pho_Nvertex_PF->Fill(NVertex);
	  }
	  if( PhoCut_PF_NoPUcorr && passHLT_Pho){
	    pho_Cands_PF_NoPUcorr.push_back(&*it_Pho);
	    h_Pho_Nvertex_PF_NoPUcorr->Fill(NVertex);
	  }

	  if(PhoCut && passHLT_Pho ){
	    numPhosRA3++;
	    pho_Cands_RA3.push_back(&*it_Pho);
	    h_Pho_Nvertex_RA3->Fill(NVertex);
	  }
	    
	  if( PhoCut_NoPUcorr && passHLT_Pho){
	    pho_Cands_RA3_NoPUcorr.push_back(&*it_Pho);
	    h_Pho_Nvertex_RA3_NoPUcorr->Fill(NVertex);
	  }

	  if(pixelCut && etaCut && EtCut && sIetaLowCut && sIphiLowCut && heCut && sIetaCut && CaloId10 && Iso50){
	    h_Pho_Nvertex_Egamma_NoIso->Fill(NVertex);
	    if(HcalIsoEgammaCut){
	      h_Pho_Nvertex_Egamma_NoTrackIsoNoEcalIso->Fill(NVertex);
	      if(HcalIsoEgammaCutNoPUcorr){
		h_Pho_Nvertex_Egamma_NoTrackIsoNoEcalIso_NoPUcorr->Fill(NVertex);
	      }
	    }//end hcal
	    if(EcalIsoEgammaCut){
	      h_Pho_Nvertex_Egamma_NoTrackIso->Fill(NVertex);
	      if(EcalIsoEgammaCutNoPUcorr){
		h_Pho_Nvertex_Egamma_NoTrackIso_NoPUcorr->Fill(NVertex);
	      }
	    }
	  }
	  
	  if( PhoCutEgamma && passHLT_Pho ){
	    numPhosVgamma++;
	    pho_Cands_Egamma.push_back(&*it_Pho);
	    h_Pho_Nvertex_Egamma->Fill(NVertex);
	  }
	    
	  if( PhoCutEgamma_NoPUcorr && passHLT_Pho){
	    pho_Cands_Egamma_NoPUcorr.push_back(&*it_Pho);
	    h_Pho_Nvertex_Egamma_NoPUcorr->Fill(NVertex);
	  }
	  
	}
	if(hasMatchBG){
	  if(1){
	    pho_Cands_Truth_BG.push_back(&*it_Pho);
	    h_Pho_Nvertex_Truth_BG->Fill(NVertex);
	  }
	    
	  if( PhoCut && passHLT_Pho ){
	    pho_Cands_RA3_BG.push_back(&*it_Pho);
	    h_Pho_Nvertex_RA3_BG->Fill(NVertex);
	  }
	    
	  if( PhoCut_NoPUcorr && passHLT_Pho){
	    pho_Cands_RA3_BG_NoPUcorr.push_back(&*it_Pho);
	    h_Pho_Nvertex_RA3_BG_NoPUcorr->Fill(NVertex);
	  }
	  
	  if( PhoCut_PF && passHLT_Pho){
	    pho_Cands_PF_BG.push_back(&*it_Pho);
	    h_Pho_Nvertex_PF_BG->Fill(NVertex);
	  }
	}
	if(printLevel > 1) cout<<"End of Photon Loop"<<endl;
	//}//if(hasMatch)
      }//for(it_Pho)
    }//phoMap

    //Now grab PF photons
    std::map<TString, std::vector<susy::Photon> >::iterator phoMapPF = event->photons.find("pfPhotonTranslator:pfphot");
    if(phoMapPF != event->photons.end()) {
      numPfPhos++;
      //loop over photon collection 
      for(std::vector<susy::Photon>::iterator it_PhoPF = phoMapPF->second.begin(); it_PhoPF != phoMapPF->second.end(); it_PhoPF++) {
	if(printLevel > 1) cout<<"looping over photon collection"<<endl;
	
	h_Pho_Et->Fill(it_PhoPF->momentum.Et());
	//if(it_PhoPF->chargedHadronIso!=0 || it_PhoPF->neutralHadronIso!=0 || it_PhoPF->photonIso!=0)cout<<"chargedHadronIso: "<<it_PhoPF->chargedHadronIso<<"  neutralHadronIso:"<<it_PhoPF->neutralHadronIso<<"  photonIso: "<<it_PhoPF->photonIso<<endl<<endl;
	//if find genMatch, put truth pho into truth histo
	bool hasMatch=false,hasMatchBG=false;
	for(std::vector<susy::Particle>::iterator Part_it = event->genParticles.begin(); 
	    Part_it != event->genParticles.end(); Part_it++){
	  if(Part_it->pdgId==22/* && Part_it->motherId==1000022*/){
	    if( isSameObject(it_PhoPF->caloPosition,Part_it->momentum,0.1) ){
	      genPart=&*Part_it;hasMatch=true;
	    }
	  }
	  if(Part_it->pdgId==22/* && Part_it->motherId==111*/){
	    if( isSameObject(it_PhoPF->caloPosition,Part_it->momentum,0.1) ){
	      genPartBG=&*Part_it;hasMatchBG=true;
	    }
	  }
	}
	//only look at reco::photons that have truth match
	//if(hasMatch || hasMatchBG){
	 
	//do things that want all photons with no cuts here:
	    
	//----------------set up cuts-------------------
	float ecalIsoDR03=it_PhoPF->ecalRecHitSumEtConeDR03;
	float hcalIsoDR03=it_PhoPF->hcalTowerSumEtConeDR03();
	float trackIsoDR03=it_PhoPF->trkSumPtHollowConeDR03;
	float chargedHadronIso=it_PhoPF->chargedHadronIso;
	float neutralHadronIso=it_PhoPF->neutralHadronIso;
	float photonIso=it_PhoPF->photonIso;
	float PhoEt=it_PhoPF->momentum.Et();
	// fiducial cuts. Look for only barrel now
	bool etaCut = (std::abs(it_PhoPF->caloPosition.Eta()) < susy::etaGapBegin);
	    
	// Spike cleaning
	bool isSpike = (it_PhoPF->r9 > 1.0);
	    
	// Et cuts, 25 GeV for trailing photons. Will apply tighter for the leading one.
	bool EtCut = (PhoEt > 25.0);
	    
	// H/E (in trigger, 0.15 for EB, 0.10 for EE)
	bool heCut = (it_PhoPF->hadronicOverEm < 0.05);
	// sigma_ietaieta (in trigger 0.014 for EB, 0.034 for EE)
	bool sIetaCut = (it_PhoPF->sigmaIetaIeta < 0.011);
	    
	    
	//combined cuts
	//  bool combIsoCut =( ( ecalIsoDR03 + hcalIsoDR03 + trackIsoDR03 ) < 6.0 );
	bool combIsoCut_NoPUcorr =( ( ecalIsoDR03 + hcalIsoDR03 + trackIsoDR03 ) < 6.0 );
	bool combIsoCut =( ( ecalIsoDR03 - PUCorr_ECAL*Rho25 + hcalIsoDR03 - PUCorr_HCAL*Rho25 + trackIsoDR03 - PUCorr_TRACK*Rho25 ) < 6.0 );
	bool PFIsoCut = (chargedHadronIso + neutralHadronIso + photonIso) / PhoEt < 0.1;
	//hlt cuts
	bool ecalIsoVLcut = ( ecalIsoDR03  < 6.0 + 0.012*PhoEt );
	bool hcalIsoVLcut = ( hcalIsoDR03 < 4.0 + 0.005*PhoEt );
	bool trackIsoVLcut = ( trackIsoDR03  < 4.0 + 0.002*PhoEt );
	bool IsoVLCut = ( ecalIsoVLcut && hcalIsoVLcut && trackIsoVLcut);
	bool CaloIdLCut = ( (it_PhoPF->hadronicOverEm < 0.15) && (it_PhoPF->sigmaIetaIeta < 0.014) );
	//Pixel cut
	bool pixelCut = (it_PhoPF->nPixelSeeds == 0);

	bool R9Id85=it_PhoPF->r9>0.85;
	bool CaloId10=( it_PhoPF->hadronicOverEm<0.1 && it_PhoPF->sigmaIetaIeta < 0.014 );
	bool Iso50=(ecalIsoDR03<5.0 + 0.012*PhoEt
		    && hcalIsoDR03<5.0 + 0.005*PhoEt
		    && trackIsoDR03<5.0+ 0.002*PhoEt
		    );
	bool R9Id85orCaloId10andIso50Cut = R9Id85 || ( CaloId10 && Iso50 );
	bool PhoCut = (pixelCut && etaCut && EtCut && !isSpike && heCut && combIsoCut && sIetaCut && CaloId10 && Iso50);
	bool PhoCut_PF = (pixelCut && etaCut && EtCut && !isSpike && heCut && PFIsoCut && sIetaCut && CaloId10 && Iso50);

	bool passHLT_Pho  = (useTrigger ? ( PassTrigger("HLT_Photon26_IsoVL_Photon18_v") 
					    || PassTrigger("HLT_Photon36_CaloIdL_Photon22_CaloIdL_v") 
					    || PassTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v") ) : true);
	    
       
	if(PhoCut)numPFPhosRA3++;
	//Now fill pho_Cands
	if(hasMatch){
	  if(1){
	    pho_Cands_Truth_PF.push_back(&*it_PhoPF);
	    h_Pho_Nvertex_Truth_PF->Fill(NVertex);
	  }
	  /*
	    if(PhoCut_PF && passHLT_Pho){
	    pho_Cands_PF.push_back(&*it_PhoPF);
	    h_Pho_Nvertex_PF->Fill(NVertex);
	    }*/
	}
	if(hasMatchBG){
	  if(1){
	    pho_Cands_Truth_PF_BG.push_back(&*it_PhoPF);
	    h_Pho_Nvertex_Truth_PF_BG->Fill(NVertex);
	  }
	    
	  if(PhoCut_PF && passHLT_Pho){
	    pho_Cands_PF_BG.push_back(&*it_PhoPF);
	    h_Pho_Nvertex_PF_BG->Fill(NVertex);
	  }
	}
	if(printLevel > 1) cout<<"End of Photon Loop"<<endl;
	//}//if(hasMatch)
      }//for(it_PhoPF)
    }//phoMapPF

    for(std::vector<susy::Photon*>::iterator pho_it_Truth = pho_Cands_Truth.begin();pho_it_Truth!=pho_Cands_Truth.end();pho_it_Truth++){
      h_Pho_Et_Truth->Fill((*pho_it_Truth)->momentum.Et());
      h_Pho_Eta_Truth->Fill((*pho_it_Truth)->momentum.Eta());
      h_Pho_Phi_Truth->Fill((*pho_it_Truth)->momentum.Phi());
    }
    for(std::vector<susy::Photon*>::iterator pho_it_Truth_PF = pho_Cands_Truth_PF.begin();pho_it_Truth_PF!=pho_Cands_Truth_PF.end();pho_it_Truth_PF++){
      h_Pho_Et_Truth_PF->Fill((*pho_it_Truth_PF)->momentum.Et());
      h_Pho_Eta_Truth_PF->Fill((*pho_it_Truth_PF)->momentum.Eta());
      h_Pho_Phi_Truth_PF->Fill((*pho_it_Truth_PF)->momentum.Phi());
    }
    for(std::vector<susy::Photon*>::iterator pho_it = pho_Cands_RA3.begin();pho_it!=pho_Cands_RA3.end();pho_it++){
      h_Pho_Et_RA3->Fill((*pho_it)->momentum.Et());
      h_Pho_Eta_RA3->Fill((*pho_it)->momentum.Eta());
      h_Pho_Phi_RA3->Fill((*pho_it)->momentum.Phi());
    }
    for(std::vector<susy::Photon*>::iterator PFpho_it = pho_Cands_PF.begin();PFpho_it!=pho_Cands_PF.end();PFpho_it++){
      h_Pho_Et_PF->Fill((*PFpho_it)->momentum.Et());
      h_Pho_Eta_PF->Fill((*PFpho_it)->momentum.Eta());
      h_Pho_Phi_PF->Fill((*PFpho_it)->momentum.Phi());
    }
    for(std::vector<susy::Photon*>::iterator pho_it_Truth_BG = pho_Cands_Truth_BG.begin();pho_it_Truth_BG!=pho_Cands_Truth_BG.end();pho_it_Truth_BG++){
      h_Pho_Et_Truth_BG->Fill((*pho_it_Truth_BG)->momentum.Et());
      h_Pho_Eta_Truth_BG->Fill((*pho_it_Truth_BG)->momentum.Eta());
      h_Pho_Phi_Truth_BG->Fill((*pho_it_Truth_BG)->momentum.Phi());
    }
    for(std::vector<susy::Photon*>::iterator pho_it_BG = pho_Cands_RA3_BG.begin();pho_it_BG!=pho_Cands_RA3_BG.end();pho_it_BG++){
      h_Pho_Et_RA3_BG->Fill((*pho_it_BG)->momentum.Et());
      h_Pho_Eta_RA3_BG->Fill((*pho_it_BG)->momentum.Eta());
      h_Pho_Phi_RA3_BG->Fill((*pho_it_BG)->momentum.Phi());
    }
    for(std::vector<susy::Photon*>::iterator PFpho_it_BG = pho_Cands_PF_BG.begin();PFpho_it_BG!=pho_Cands_PF_BG.end();PFpho_it_BG++){
      h_Pho_Et_PF_BG->Fill((*PFpho_it_BG)->momentum.Et());
      h_Pho_Eta_PF_BG->Fill((*PFpho_it_BG)->momentum.Eta());
      h_Pho_Phi_PF_BG->Fill((*PFpho_it_BG)->momentum.Phi());
    }

  }//for jentry
  TH1F* num=(TH1F*)h_Pho_Et_RA3->Rebin(NumBins,"num",EtBins);TH1F* den=(TH1F*)h_Pho_Et_Truth->Rebin(NumBins,"den",EtBins);
  TGraphAsymmErrors *h_et_RA3 = new TGraphAsymmErrors(num,den,"");
  h_et_RA3->SetMarkerSize(0.6);
  h_et_RA3->SetTitle("");
  h_et_RA3->GetXaxis()->SetTitle("E_{T}");
  h_et_RA3->GetYaxis()->SetTitle("Efficiency");
  h_et_RA3->Write("Eff_et_RA3");

  TGraphAsymmErrors *h_eta_RA3 = new TGraphAsymmErrors(h_Pho_Eta_RA3,h_Pho_Eta_Truth,"");
  h_eta_RA3->SetMarkerSize(0.6);
  h_eta_RA3->GetXaxis()->SetTitle("#eta");
  h_eta_RA3->GetYaxis()->SetTitle("Efficiency");
  h_eta_RA3->Write("Eff_eta_RA3");

  TGraphAsymmErrors *h_phi_RA3 = new TGraphAsymmErrors(h_Pho_Phi_RA3,h_Pho_Phi_Truth,"");
  h_phi_RA3->SetMarkerSize(0.6);
  h_phi_RA3->GetXaxis()->SetTitle("#phi");
  h_phi_RA3->GetYaxis()->SetTitle("Efficiency");
  h_phi_RA3->Write("Eff_phi_RA3");

  TGraphAsymmErrors *h_nvertex_RA3 = new TGraphAsymmErrors(h_Pho_Nvertex_RA3,h_Pho_Nvertex_Truth,"");
  h_nvertex_RA3->SetMarkerSize(0.6);
  h_nvertex_RA3->GetXaxis()->SetTitle("nVertex");
  h_nvertex_RA3->GetYaxis()->SetTitle("Efficiency");
  h_nvertex_RA3->Write("Eff_nvertex_RA3");
 
  TGraphAsymmErrors *h_nvertex_RA3_NoPUcorr = new TGraphAsymmErrors(h_Pho_Nvertex_RA3_NoPUcorr,h_Pho_Nvertex_Truth,"");
  h_nvertex_RA3_NoPUcorr->SetMarkerSize(0.6);
  h_nvertex_RA3_NoPUcorr->GetXaxis()->SetTitle("nVertex");
  h_nvertex_RA3_NoPUcorr->GetYaxis()->SetTitle("Efficiency");
  h_nvertex_RA3_NoPUcorr->SetLineColor(kRed);
  h_nvertex_RA3_NoPUcorr->SetMarkerColor(kRed);
  h_nvertex_RA3_NoPUcorr->Write("Eff_nvertex_RA3_NoPUcorr");

  TGraphAsymmErrors *h_nvertex_Egamma = new TGraphAsymmErrors(h_Pho_Nvertex_Egamma,h_Pho_Nvertex_Truth,"");
  h_nvertex_Egamma->SetMarkerSize(0.6);
  h_nvertex_Egamma->GetXaxis()->SetTitle("nVertex");
  h_nvertex_Egamma->GetYaxis()->SetTitle("Efficiency");
  h_nvertex_Egamma->Write("Eff_nvertex_Egamma");
 
  TGraphAsymmErrors *h_nvertex_Egamma_NoPUcorr = new TGraphAsymmErrors(h_Pho_Nvertex_Egamma_NoPUcorr,h_Pho_Nvertex_Truth,"");
  h_nvertex_Egamma_NoPUcorr->SetMarkerSize(0.6);
  h_nvertex_Egamma_NoPUcorr->GetXaxis()->SetTitle("nVertex");
  h_nvertex_Egamma_NoPUcorr->GetYaxis()->SetTitle("Efficiency");
  h_nvertex_Egamma_NoPUcorr->SetLineColor(kRed);
  h_nvertex_Egamma_NoPUcorr->SetMarkerColor(kRed);
  h_nvertex_Egamma_NoPUcorr->Write("Eff_nvertex_Egamma_NoPUcorr");

  TGraphAsymmErrors *h_nvertex_Egamma_NoIso = new TGraphAsymmErrors(h_Pho_Nvertex_Egamma_NoIso,h_Pho_Nvertex_Truth,"");
  h_nvertex_Egamma_NoIso->SetMarkerSize(0.6);
  h_nvertex_Egamma_NoIso->GetXaxis()->SetTitle("nVertex");
  h_nvertex_Egamma_NoIso->GetYaxis()->SetTitle("Efficiency");
  h_nvertex_Egamma_NoIso->Write("Eff_nvertex_Egamma_NoIso");

  TGraphAsymmErrors *h_nvertex_Egamma_NoTrackIso = new TGraphAsymmErrors(h_Pho_Nvertex_Egamma_NoTrackIso,h_Pho_Nvertex_Truth,"");
  h_nvertex_Egamma_NoTrackIso->SetMarkerSize(0.6);
  h_nvertex_Egamma_NoTrackIso->GetXaxis()->SetTitle("nVertex");
  h_nvertex_Egamma_NoTrackIso->GetYaxis()->SetTitle("Efficiency");
  h_nvertex_Egamma_NoTrackIso->Write("Eff_nvertex_Egamma_NoTrackIso");
 
  TGraphAsymmErrors *h_nvertex_Egamma_NoTrackIso_NoPUcorr = new TGraphAsymmErrors(h_Pho_Nvertex_Egamma_NoTrackIso_NoPUcorr,h_Pho_Nvertex_Truth,"");
  h_nvertex_Egamma_NoTrackIso_NoPUcorr->SetMarkerSize(0.6);
  h_nvertex_Egamma_NoTrackIso_NoPUcorr->GetXaxis()->SetTitle("nVertex");
  h_nvertex_Egamma_NoTrackIso_NoPUcorr->GetYaxis()->SetTitle("Efficiency");
  h_nvertex_Egamma_NoTrackIso_NoPUcorr->SetLineColor(kRed);
  h_nvertex_Egamma_NoTrackIso_NoPUcorr->SetMarkerColor(kRed);
  h_nvertex_Egamma_NoTrackIso_NoPUcorr->Write("Eff_nvertex_Egamma_NoTrackIso_NoPUcorr");

  TGraphAsymmErrors *h_nvertex_Egamma_NoTrackIsoNoEcalIso = new TGraphAsymmErrors(h_Pho_Nvertex_Egamma_NoTrackIsoNoEcalIso,h_Pho_Nvertex_Truth,"");
  h_nvertex_Egamma_NoTrackIsoNoEcalIso->SetMarkerSize(0.6);
  h_nvertex_Egamma_NoTrackIsoNoEcalIso->GetXaxis()->SetTitle("nVertex");
  h_nvertex_Egamma_NoTrackIsoNoEcalIso->GetYaxis()->SetTitle("Efficiency");
  h_nvertex_Egamma_NoTrackIsoNoEcalIso->Write("Eff_nvertex_Egamma_NoTrackIsoNoEcalIso");
 
  TGraphAsymmErrors *h_nvertex_Egamma_NoTrackIsoNoEcalIso_NoPUcorr = new TGraphAsymmErrors(h_Pho_Nvertex_Egamma_NoTrackIsoNoEcalIso_NoPUcorr,h_Pho_Nvertex_Truth,"");
  h_nvertex_Egamma_NoTrackIsoNoEcalIso_NoPUcorr->SetMarkerSize(0.6);
  h_nvertex_Egamma_NoTrackIsoNoEcalIso_NoPUcorr->GetXaxis()->SetTitle("nVertex");
  h_nvertex_Egamma_NoTrackIsoNoEcalIso_NoPUcorr->GetYaxis()->SetTitle("Efficiency");
  h_nvertex_Egamma_NoTrackIsoNoEcalIso_NoPUcorr->SetLineColor(kRed);
  h_nvertex_Egamma_NoTrackIsoNoEcalIso_NoPUcorr->SetMarkerColor(kRed);
  h_nvertex_Egamma_NoTrackIsoNoEcalIso_NoPUcorr->Write("Eff_nvertex_Egamma_NoTrackIsoNoEcalIso_NoPUcorr");

  //TH1F* numEt_PF=(TH1F*)h_Pho_Et_PF->Rebin(NumBins,"numEt_PF",EtBins);TH1F* denEt_PF=(TH1F*)h_Pho_Et_Truth_PF->Rebin(NumBins,"denEt_PF",EtBins);
  TGraphAsymmErrors *h_et_PF = new TGraphAsymmErrors(h_Pho_Et_PF,h_Pho_Et_Truth,"");
  h_et_PF->SetMarkerSize(0.6);
  h_et_PF->GetXaxis()->SetTitle("E_{T}");
  h_et_PF->GetYaxis()->SetTitle("Efficiency");
  h_et_PF->SetLineColor(kBlue);
  h_et_PF->SetMarkerColor(kBlue);
  h_et_PF->Write("Eff_et_PF");

  TGraphAsymmErrors *h_eta_PF = new TGraphAsymmErrors(h_Pho_Eta_PF,h_Pho_Eta_Truth,"");
  h_eta_PF->SetMarkerSize(0.6);
  h_eta_PF->GetXaxis()->SetTitle("#eta");
  h_eta_PF->GetYaxis()->SetTitle("Efficiency");
  h_eta_PF->SetLineColor(kBlue);
  h_eta_PF->SetMarkerColor(kBlue);
  h_eta_PF->Write("Eff_eta_PF");

  TGraphAsymmErrors *h_phi_PF = new TGraphAsymmErrors(h_Pho_Phi_PF,h_Pho_Phi_Truth,"");
  h_phi_PF->SetMarkerSize(0.6);
  h_phi_PF->GetXaxis()->SetTitle("#phi");
  h_phi_PF->GetYaxis()->SetTitle("Efficiency");
  h_phi_PF->SetLineColor(kBlue);
  h_phi_PF->SetMarkerColor(kBlue);
  h_phi_PF->Write("Eff_phi_PF");

  TGraphAsymmErrors *h_nvertex_PF = new TGraphAsymmErrors(h_Pho_Nvertex_PF,h_Pho_Nvertex_Truth,"");
  h_nvertex_PF->SetMarkerSize(0.6);
  h_nvertex_PF->GetXaxis()->SetTitle("nVertex");
  h_nvertex_PF->GetYaxis()->SetTitle("Efficiency");
  h_nvertex_PF->SetLineColor(kBlue);
  h_nvertex_PF->SetMarkerColor(kBlue);
  h_nvertex_PF->Write("Eff_nvertex_PF");
    
  TGraphAsymmErrors *h_nvertex_PF_NoPUcorr = new TGraphAsymmErrors(h_Pho_Nvertex_PF_NoPUcorr,h_Pho_Nvertex_Truth,"");
  h_nvertex_PF_NoPUcorr->SetMarkerSize(0.6);
  h_nvertex_PF_NoPUcorr->GetXaxis()->SetTitle("nVertex");
  h_nvertex_PF_NoPUcorr->GetYaxis()->SetTitle("Efficiency");
  h_nvertex_PF_NoPUcorr->SetLineColor(kBlue);
  h_nvertex_PF_NoPUcorr->SetMarkerColor(kBlue);
  h_nvertex_PF_NoPUcorr->Write("Eff_nvertex_PF_NoPUcorr");
    
  //Now BG
  TGraphAsymmErrors *h_et_RA3_BG = new TGraphAsymmErrors(h_Pho_Et_RA3_BG,h_Pho_Et_Truth_BG,"");
  h_et_RA3_BG->SetMarkerSize(0.6);
  h_et_RA3_BG->GetXaxis()->SetTitle("E_{T}");
  h_et_RA3_BG->GetYaxis()->SetTitle("Efficiency");
  h_et_RA3_BG->Write("Eff_et_RA3_BG");

  TGraphAsymmErrors *h_eta_RA3_BG = new TGraphAsymmErrors(h_Pho_Eta_RA3_BG,h_Pho_Eta_Truth_BG,"");
  h_eta_RA3_BG->SetMarkerSize(0.6);
  h_eta_RA3_BG->GetXaxis()->SetTitle("#eta");
  h_eta_RA3_BG->GetYaxis()->SetTitle("Efficiency");
  h_eta_RA3_BG->Write("Eff_eta_RA3_BG");

  TGraphAsymmErrors *h_phi_RA3_BG = new TGraphAsymmErrors(h_Pho_Phi_RA3_BG,h_Pho_Phi_Truth_BG,"");
  h_phi_RA3_BG->SetMarkerSize(0.6);
  h_phi_RA3_BG->GetXaxis()->SetTitle("#phi");
  h_phi_RA3_BG->GetYaxis()->SetTitle("Efficiency");
  h_phi_RA3_BG->Write("Eff_phi_RA3_BG");

  TGraphAsymmErrors *h_nvertex_RA3_BG = new TGraphAsymmErrors(h_Pho_Nvertex_RA3_BG,h_Pho_Nvertex_Truth_BG,"");
  h_nvertex_RA3_BG->SetMarkerSize(0.6);
  h_nvertex_RA3_BG->GetXaxis()->SetTitle("nVertex");
  h_nvertex_RA3_BG->GetYaxis()->SetTitle("Efficiency");
  h_nvertex_RA3_BG->Write("Eff_nvertex_RA3_BG");
 
  TGraphAsymmErrors *h_nvertex_RA3_BG_NoPUcorr = new TGraphAsymmErrors(h_Pho_Nvertex_RA3_BG_NoPUcorr,h_Pho_Nvertex_Truth_BG,"");
  h_nvertex_RA3_BG_NoPUcorr->SetMarkerSize(0.6);
  h_nvertex_RA3_BG_NoPUcorr->GetXaxis()->SetTitle("nVertex");
  h_nvertex_RA3_BG_NoPUcorr->GetYaxis()->SetTitle("Efficiency");
  h_nvertex_RA3_BG_NoPUcorr->SetLineColor(kRed);
  h_nvertex_RA3_BG_NoPUcorr->SetMarkerColor(kRed);
  h_nvertex_RA3_BG_NoPUcorr->Write("Eff_nvertex_RA3_BG_NoPUcorr");

  TGraphAsymmErrors *h_et_PF_BG = new TGraphAsymmErrors(h_Pho_Et_PF_BG,h_Pho_Et_Truth_BG,"");
  h_et_PF_BG->SetMarkerSize(0.6);
  h_et_PF_BG->GetXaxis()->SetTitle("E_{T}");
  h_et_PF_BG->GetYaxis()->SetTitle("Efficiency");
  h_et_PF_BG->SetLineColor(kBlue);
  h_et_PF_BG->SetMarkerColor(kBlue);
  h_et_PF_BG->Write("Eff_et_PF_BG");

  TGraphAsymmErrors *h_eta_PF_BG = new TGraphAsymmErrors(h_Pho_Eta_PF_BG,h_Pho_Eta_Truth_BG,"");
  h_eta_PF_BG->SetMarkerSize(0.6);
  h_eta_PF_BG->GetXaxis()->SetTitle("#eta");
  h_eta_PF_BG->GetYaxis()->SetTitle("Efficiency");
  h_eta_PF_BG->SetLineColor(kBlue);
  h_eta_PF_BG->SetMarkerColor(kBlue);
  h_eta_PF_BG->Write("Eff_eta_PF_BG");

  TGraphAsymmErrors *h_phi_PF_BG = new TGraphAsymmErrors(h_Pho_Phi_PF_BG,h_Pho_Phi_Truth_BG,"");
  h_phi_PF_BG->SetMarkerSize(0.6);
  h_phi_PF_BG->GetXaxis()->SetTitle("#phi");
  h_phi_PF_BG->GetYaxis()->SetTitle("Efficiency");
  h_phi_PF_BG->SetLineColor(kBlue);
  h_phi_PF_BG->SetMarkerColor(kBlue);
  h_phi_PF_BG->Write("Eff_phi_PF_BG");

  TGraphAsymmErrors *h_nvertex_PF_BG = new TGraphAsymmErrors(h_Pho_Nvertex_PF_BG,h_Pho_Nvertex_Truth_BG,"");
  h_nvertex_PF_BG->SetMarkerSize(0.6);
  h_nvertex_PF_BG->GetXaxis()->SetTitle("nVertex");
  h_nvertex_PF_BG->GetYaxis()->SetTitle("Efficiency");
  h_nvertex_PF_BG->SetLineColor(kBlue);
  h_nvertex_PF_BG->SetMarkerColor(kBlue);
  h_nvertex_PF_BG->Write("Eff_nvertex_PF_BG");
    
  c1->cd(1);
  h_Pho_Et_Truth->SetLineColor(kBlack);
  h_Pho_Et_RA3->SetLineColor(kRed);
  h_Pho_Et_PF->SetLineColor(kBlue);
  h_Pho_Et_Truth->Draw();
  h_Pho_Et_RA3->Draw("SAMES");
  h_Pho_Et_PF->Draw("SAMES");
  TLegend *leg = new TLegend(.45,.4,.85,.7,"","brNDC");
  leg->AddEntry(h_Pho_Et_Truth,"reco::Photon matched to genPhoton","l");
  leg->AddEntry(h_Pho_Et_RA3,"reco::Photon passing RA3 cuts (STD Iso)","l");
  leg->AddEntry(h_Pho_Et_PF,"reco::Photon passing RA3 cuts (PF Iso)","l");
  leg->SetFillColor(kWhite);
  leg->Draw("SAME");
  c1->cd(2);
  h_et_RA3->SetLineColor(kRed);
  h_et_RA3->SetMarkerColor(kRed);
  h_et_PF->SetLineColor(kBlue);
  h_et_PF->SetMarkerColor(kBlue);
  h_et_PF->Draw("AP");
  h_et_RA3->Draw("sameP");
  c1->Write("Eff_et");

  cout << "---ALL DONE!------"<<endl;
    
  cout << " ----------------- Job Summary ----------------- " << endl;
  cout << " Total events                : " << nCnt[0] << " (" << 100*nCnt[0]/float(nCnt[0]) << "%)" << endl;
  cout << " JSON events passed          : " << nCnt[1] << " (" << 100*nCnt[1]/float(nCnt[0]) << "%)" << endl;
  cout << " HLT events passed           : " << nCnt[2] << " (" << 100*nCnt[2]/float(nCnt[0]) << "%)" << endl;
  cout << " nVertex events passed       : " << nCnt[3] << " (" << 100*nCnt[3]/float(nCnt[0]) << "%)" << endl;
  cout << " Number of photons           : " << numPhos << endl;
  cout << " Number of photons from Chi0 : " << numPhosTruth << endl;
  cout << " Number of RA3 photons       : " << numPhosRA3 << endl;
  cout << " Number of Vgamma photons    : " << numPhosVgamma << endl;
  cout << " Number of PF photons        : " << numPfPhos << endl;
  cout<<" numPho: "<<numPho<<endl;
  cout<<" numPhoMatch: "<<numPhoMatch<<endl;
  cout << " ----------------------------------------------- " << endl;
    
  cout<<"Writing root output to: /data/ndpc3/b/dmorse/RA3/AnalysisOutput/PhotonId/hist_"<<ds<<".root"<<endl;
    
    
  // close the output file
  fout->cd();
  fout->Write();
  fout->Close();
  delete fout;
    
}

