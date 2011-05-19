#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include "SusyEventPrinter.h"

void Print(const TVector2& p) {
  float mod = p.Mod();
  if(mod < 1e-6) std::cout << "(0)" << std::endl;
  else {
    std::cout << "(" << p.X() << "," << p.Y() << "; mod=" << mod << "; phi=" << p.Phi() << ")" << std::endl;
  }

}

void Print(const TVector3& p) {
  float mag = p.Mag();
  if(mag < 1e-6) std::cout << "(0)" << std::endl;
  else {
    std::cout << "(" << p.Px() << "," << p.Py() << "," << p.Pz() << "; mag=" << mag << "; pt=" << p.Pt() << "; eta=" << p.Eta() << "; phi=" << p.Phi() << ")" << std::endl;
  }
}


void Print(const TLorentzVector& p) {
  float mag = p.Vect().Mag();
  if(mag < 1e-6) std::cout << "(0)" << std::endl;
  else {
    std::cout << "(" << p.Px() << "," << p.Py() << "," << p.Pz() << ";" << p.E() << "; m=" << p.M() << "; mag=" << mag << "; pt=" << p.Pt() << "; eta=" << p.Eta() << "; phi=" << p.Phi() << ")" << std::endl;
  }
}


void Print(const susy::TriggerMap& t) {
  for(susy::TriggerMap::const_iterator it = t.begin(); it != t.end(); it++) {
    std::cout << "\t" << it->first << "(prescale=" << it->second.first << ", fire=" << int(it->second.second) << ")" << std::endl;
  }
}


void Print(const susy::CorrMETData& cm) {
  std::cout << "\t\tdmEx:" << cm.dmEx << ", dmEy:" << cm.dmEy << ", dsumEt:" << cm.dsumEt << ", dSignificance:" << cm.dSignificance << std::endl;
}


void Print(const susy::MET& met) {
  std::cout << "\tsumEt : " << met.sumEt << std::endl;
  std::cout << "\tsignificance : " << met.significance << std::endl;
  std::cout << "\tmEt : "; Print(met.mEt);
  std::cout << "\tvertex : "; Print(met.vertex);
  std::cout << "\tmEtCorr ===>" << std::endl;
  for(std::vector<susy::CorrMETData>::const_iterator it = met.mEtCorr.begin(); it != met.mEtCorr.end(); it++) Print(*it);
  std::cout << std::endl;
}


void Print(const susy::Event& event) {

  std::cout << "---------- run(" << event.runNumber << "), event(" << event.eventNumber << ") ----------" << std::endl;

  std::cout << "isRealData : " << int(event.isRealData) << std::endl;
  std::cout << "luminosityBlockNumber : " << event.luminosityBlockNumber << std::endl;
  std::cout << "bunchCrossing : " << event.bunchCrossing << std::endl;
  std::cout << "avgInsRecLumi : " << event.avgInsRecLumi << std::endl;
  std::cout << "intgRecLumi : " << event.intgRecLumi << std::endl;
  std::cout << "cosmicFlag : " << int(event.cosmicFlag) << std::endl;

  std::cout << "beamSpot"; Print(event.beamSpot);

  std::cout << "l1Map size(" << event.l1Map.size() << ") =========>" << std::endl;
  Print(event.l1Map);
  std::cout << std::endl;

  std::cout << "hltMap size(" << event.hltMap.size() << ") =========>" << std::endl;
  Print(event.hltMap);
  std::cout << std::endl;

  std::cout << "metMap size(" << event.metMap.size() << ") =========>" << std::endl;
  for(std::map<TString,susy::MET>::const_iterator it = event.metMap.begin(); it != event.metMap.end(); it++) {
    std::cout << it->first << " ======>" << std::endl;
    Print(it->second);
  }
  std::cout << std::endl;

  std::cout << "vertices size(" << event.vertices.size() << ") =========>" << std::endl;
  std::cout << "tracks size(" << event.tracks.size() << ") =========>" << std::endl;
  std::cout << "superClusters size(" << event.superClusters.size() << ") =========>" << std::endl;
  std::cout << "clusters size(" << event.clusters.size() << ") =========>" << std::endl;
  std::cout << "muons size(" << event.muons.size() << ") =========>" << std::endl;
  std::cout << "electrons size(" << event.electrons.size() << ") =========>" << std::endl;
  std::cout << "photons size(" << event.photons.size() << ") =========>" << std::endl;
  std::cout << "caloJets size(" << event.caloJets.size() << ") =========>" << std::endl;
  std::cout << "pfJets size(" << event.pfJets.size() << ") =========>" << std::endl;
  std::cout << "jptJets size(" << event.jptJets.size() << ") =========>" << std::endl;
  std::cout << "generalTracks size(" << event.generalTracks.size() << ") =========>" << std::endl;
  std::cout << "simVertices size(" << event.simVertices.size() << ") =========>" << std::endl;
  std::cout << "genParticles size(" << event.genParticles.size() << ") =========>" << std::endl;
  std::cout << "gridParams size(" << event.gridParams.size() << ") =========>" << std::endl;

  std::cout << std::endl;
}
