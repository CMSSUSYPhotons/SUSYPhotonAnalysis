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


void Print(const susy::Photon& p) {

  std::cout << "\tfidBit : " << p.fidBit << std::endl;
  std::cout << "\tnPixelSeeds : " << p.nPixelSeeds << std::endl;
  std::cout << "\thadronicOverEm : " << p.hadronicOverEm << std::endl;
  std::cout << "\thadronicDepth1OverEm : " << p.hadronicDepth1OverEm << std::endl;
  std::cout << "\thadronicDepth2OverEm : " << p.hadronicDepth2OverEm << std::endl;
  std::cout << "\te1x2 : " << p.e1x2 << std::endl;
  std::cout << "\te1x5 : " << p.e1x5 << std::endl;
  std::cout << "\te2x5 : " << p.e2x5 << std::endl;
  std::cout << "\te3x3 : " << p.e3x3 << std::endl;
  std::cout << "\te5x5 : " << p.e5x5 << std::endl;
  std::cout << "\tmaxEnergyXtal : " << p.maxEnergyXtal << std::endl;
  std::cout << "\tsigmaEtaEta : " << p.sigmaEtaEta << std::endl;
  std::cout << "\tsigmaIetaIeta : " << p.sigmaIetaIeta << std::endl;
  std::cout << "\tr9 : " << p.r9 << std::endl;
  std::cout << "\tecalRecHitSumEtConeDR04 : " << p.ecalRecHitSumEtConeDR04 << std::endl;
  std::cout << "\thcalDepth1TowerSumEtConeDR04 : " << p.hcalDepth1TowerSumEtConeDR04 << std::endl;
  std::cout << "\thcalDepth2TowerSumEtConeDR04 : " << p.hcalDepth2TowerSumEtConeDR04 << std::endl;
  std::cout << "\ttrkSumPtSolidConeDR04 : " << p.trkSumPtSolidConeDR04 << std::endl;
  std::cout << "\ttrkSumPtHollowConeDR04 : " << p.trkSumPtHollowConeDR04 << std::endl;
  std::cout << "\tnTrkSolidConeDR04 : " << (int)p.nTrkSolidConeDR04 << std::endl;
  std::cout << "\tnTrkHollowConeDR04 : " << (int)p.nTrkHollowConeDR04 << std::endl;
  std::cout << "\tecalRecHitSumEtConeDR03 : " << p.ecalRecHitSumEtConeDR03 << std::endl;
  std::cout << "\thcalDepth1TowerSumEtConeDR03 : " << p.hcalDepth1TowerSumEtConeDR03 << std::endl;
  std::cout << "\thcalDepth2TowerSumEtConeDR03 : " << p.hcalDepth2TowerSumEtConeDR03 << std::endl;
  std::cout << "\ttrkSumPtSolidConeDR03 : " << p.trkSumPtSolidConeDR03 << std::endl;
  std::cout << "\ttrkSumPtHollowConeDR03 : " << p.trkSumPtHollowConeDR03 << std::endl;
  std::cout << "\tnTrkSolidConeDR03 : " << (int)p.nTrkSolidConeDR03 << std::endl;
  std::cout << "\tnTrkHollowConeDR03 : " << (int)p.nTrkHollowConeDR03 << std::endl;
  std::cout << "\tchargedHadronIso : " << p.chargedHadronIso << std::endl;
  std::cout << "\tneutralHadronIso : " << p.neutralHadronIso << std::endl;
  std::cout << "\tphotonIso : " << p.photonIso << std::endl;
  std::cout << "\tseedTime : " << p.seedTime << std::endl;
  std::cout << "\tsMaj : " << p.sMaj << std::endl;
  std::cout << "\tsMin : " << p.sMin << std::endl;
  std::cout << "\talpha : " << p.alpha << std::endl;
  std::cout << "\troundness : " << p.roundness << std::endl;
  std::cout << "\tangle : " << p.angle << std::endl;
  std::cout << "\tconvDist : " << p.convDist << std::endl;
  std::cout << "\tconvDcot : " << p.convDcot << std::endl;
  std::cout << "\tconvVtxChi2 : " << p.convVtxChi2 << std::endl;
  std::cout << "\tconvVtxNdof : " << (int)p.convVtxNdof << std::endl;
  std::cout << "\tconvVertex : "; Print(p.convVertex);
  std::cout << "\tsuperClusterIndex : " << p.superClusterIndex << std::endl;
  std::cout << "\tsuperClusterPreshowerEnergy : " << p.superClusterPreshowerEnergy << std::endl;
  std::cout << "\tsuperClusterPhiWidth : " << p.superClusterPhiWidth << std::endl;
  std::cout << "\tsuperClusterEtaWidth : " << p.superClusterEtaWidth << std::endl;
  std::cout << "\tcaloPosition : "; Print(p.caloPosition);
  std::cout << "\tvertex : "; Print(p.vertex);
  std::cout << "\tmomentum : "; Print(p.momentum);
  std::cout << "\tidPairs : ";
  for(std::map<TString,UChar_t>::const_iterator it = p.idPairs.begin(); it != p.idPairs.end(); it++ ){
    std::cout << "(" << it->first << ", " << int(it->second) << ") ";
  }
  std::cout << std::endl;

}


void Print(const susy::CaloJet& j) {
  std::cout << "\tpartonFlavour : " << j.partonFlavour << std::endl;
  std::cout << "\tjetCharge : " << j.jetCharge << std::endl;
  std::cout << "\tetaMean : " << j.etaMean << std::endl;
  std::cout << "\tphiMean : " << j.phiMean << std::endl;
  std::cout << "\tetaEtaMoment : " << j.etaEtaMoment << std::endl;
  std::cout << "\tetaPhiMoment : " << j.etaPhiMoment << std::endl;
  std::cout << "\tphiPhiMoment : " << j.phiPhiMoment << std::endl;
  std::cout << "\tmaxDistance : " << j.maxDistance << std::endl;
  std::cout << "\tjetArea : " << j.jetArea << std::endl;
  std::cout << "\tpileup : " << j.pileup << std::endl;
  std::cout << "\tnPasses : " << (int)j.nPasses << std::endl;
  std::cout << "\tnConstituents : " << (int)j.nConstituents << std::endl;
  std::cout << "\tmaxEInEmTowers : " << j.maxEInEmTowers << std::endl;
  std::cout << "\tmaxEInHadTowers : " << j.maxEInHadTowers << std::endl;
  std::cout << "\tenergyFractionHadronic : " << j.energyFractionHadronic << std::endl;
  std::cout << "\temEnergyFraction : " << j.emEnergyFraction << std::endl;
  std::cout << "\thadEnergyInHB : " << j.hadEnergyInHB << std::endl;
  std::cout << "\thadEnergyInHO : " << j.hadEnergyInHO << std::endl;
  std::cout << "\thadEnergyInHE : " << j.hadEnergyInHE << std::endl;
  std::cout << "\thadEnergyInHF : " << j.hadEnergyInHF << std::endl;
  std::cout << "\temEnergyInEB : " << j.emEnergyInEB << std::endl;
  std::cout << "\temEnergyInEE : " << j.emEnergyInEE << std::endl;
  std::cout << "\temEnergyInHF : " << j.emEnergyInHF << std::endl;
  std::cout << "\ttowersArea : " << j.towersArea << std::endl;
  std::cout << "\tn90 : " << (int)j.n90 << std::endl;
  std::cout << "\tn60 : " << (int)j.n60 << std::endl;
  std::cout << "\tfHPD : " << j.fHPD << std::endl;
  std::cout << "\tfRBX : " << j.fRBX << std::endl;
  std::cout << "\tn90Hits : " << j.n90Hits << std::endl;
  std::cout << "\tfSubDetector1 : " << j.fSubDetector1 << std::endl;
  std::cout << "\tfSubDetector2 : " << j.fSubDetector2 << std::endl;
  std::cout << "\tfSubDetector3 : " << j.fSubDetector3 << std::endl;
  std::cout << "\tfSubDetector4 : " << j.fSubDetector4 << std::endl;
  std::cout << "\trestrictedEMF : " << j.restrictedEMF << std::endl;
  std::cout << "\tnHCALTowers : " << (int)j.nHCALTowers << std::endl;
  std::cout << "\tnECALTowers : " << (int)j.nECALTowers << std::endl;
  std::cout << "\tapproximatefHPD : " << j.approximatefHPD << std::endl;
  std::cout << "\tapproximatefRBX : " << j.approximatefRBX << std::endl;
  std::cout << "\thitsInN90 : " << (int)j.hitsInN90 << std::endl;
  std::cout << "\tnumberOfHits2RPC : " << (int)j.numberOfHits2RPC << std::endl;
  std::cout << "\tnumberOfHits3RPC : " << (int)j.numberOfHits3RPC << std::endl;
  std::cout << "\tnumberOfHitsRPC : " << (int)j.numberOfHitsRPC << std::endl;
  std::cout << "\tvertex : "; Print(j.vertex);
  std::cout << "\tmomentum : "; Print(j.momentum);
  std::cout << "\tdetectorP4 : "; Print(j.detectorP4);
  std::cout << "\tjecScaleFactors : ";
  for(std::map<TString, Float_t>::const_iterator it = j.jecScaleFactors.begin();
      it != j.jecScaleFactors.end(); it++) {
    std::cout << "(" << it->first << "," << it->second << ") ";
  }
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
  for(std::vector<TVector3>::const_iterator it = event.vertices.begin(); it != event.vertices.end(); it++) {
    std::cout << "\t"; Print(*it);
  }
  std::cout << std::endl;

  std::cout << "tracks size(" << event.tracks.size() << ") =========>" << std::endl;
  std::cout << "superClusters size(" << event.superClusters.size() << ") =========>" << std::endl;
  std::cout << "clusters size(" << event.clusters.size() << ") =========>" << std::endl;
  std::cout << "muons size(" << event.muons.size() << ") =========>" << std::endl;
  std::cout << "electrons size(" << event.electrons.size() << ") =========>" << std::endl;

  std::cout << "photons size(" << event.photons.size() << ") =========>" << std::endl;
  for(std::map<TString,susy::PhotonCollection>::const_iterator it = event.photons.begin(); it != event.photons.end(); it++) {
    std::cout << it->first << " ======>" << std::endl;
    for(std::vector<susy::Photon>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
      Print(*it2);
    }
  }

  std::cout << "caloJets size(" << event.caloJets.size() << ") =========>" << std::endl;
  for(std::map<TString, susy::CaloJetCollection>::const_iterator it = event.caloJets.begin(); it != event.caloJets.end(); it++) {
    std::cout << it->first << " ======>" << std::endl;
    for(std::vector<susy::CaloJet>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
      Print(*it2);
    }
  }

  std::cout << "pfJets size(" << event.pfJets.size() << ") =========>" << std::endl;
  std::cout << "jptJets size(" << event.jptJets.size() << ") =========>" << std::endl;
  std::cout << "generalTracks size(" << event.generalTracks.size() << ") =========>" << std::endl;
  std::cout << "simVertices size(" << event.simVertices.size() << ") =========>" << std::endl;
  std::cout << "genParticles size(" << event.genParticles.size() << ") =========>" << std::endl;
  std::cout << "gridParams size(" << event.gridParams.size() << ") =========>" << std::endl;

  std::cout << std::endl;
}
