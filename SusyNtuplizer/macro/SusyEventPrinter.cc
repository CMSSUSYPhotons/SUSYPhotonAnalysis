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


void Print(const susy::Vertex& vtx) {
  std::cout << "\tchi2 : " << vtx.chi2 << std::endl;
  std::cout << "\tndof : " << vtx.ndof << std::endl;
  std::cout << "\ttracksSize : " << (int)vtx.tracksSize << std::endl;
  std::cout << "\tposition : "; Print(vtx.position);
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
  std::cout << "\tsigmaIphiIphi : " << p.sigmaIphiIphi << std::endl;
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
  std::cout << "\tchargedHadronIsoDeposit : " << p.chargedHadronIsoDeposit << std::endl;
  std::cout << "\tneutralHadronIsoDeposit : " << p.neutralHadronIsoDeposit << std::endl;
  std::cout << "\tphotonIsoDeposit : " << p.photonIsoDeposit << std::endl;
  std::cout << "\tseedTime : " << p.seedTime << std::endl;
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


void Print(const susy::PFJet& j) {
  std::cout << "\tphyDefFlavour : " << j.phyDefFlavour << std::endl;
  std::cout << "\talgDefFlavour : " << j.algDefFlavour << std::endl;
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
  std::cout << "\tchargedHadronEnergy : " << j.chargedHadronEnergy << std::endl;
  std::cout << "\tneutralHadronEnergy : " << j.neutralHadronEnergy << std::endl;
  std::cout << "\tphotonEnergy : " << j.photonEnergy << std::endl;
  std::cout << "\telectronEnergy : " << j.electronEnergy << std::endl;
  std::cout << "\tmuonEnergy : " << j.muonEnergy << std::endl;
  std::cout << "\tHFHadronEnergy : " << j.HFHadronEnergy << std::endl;
  std::cout << "\tHFEMEnergy : " << j.HFEMEnergy << std::endl;
  std::cout << "\tchargedEmEnergy : " << j.chargedEmEnergy << std::endl;
  std::cout << "\tchargedMuEnergy : " << j.chargedMuEnergy << std::endl;
  std::cout << "\tneutralEmEnergy : " << j.neutralEmEnergy << std::endl;
  std::cout << "\tchargedHadronMultiplicity : " << (int)j.chargedHadronMultiplicity << std::endl;
  std::cout << "\tneutralHadronMultiplicity : " << (int)j.neutralHadronMultiplicity << std::endl;
  std::cout << "\tphotonMultiplicity : " << (int)j.photonMultiplicity << std::endl;
  std::cout << "\telectronMultiplicity : " << (int)j.electronMultiplicity << std::endl;
  std::cout << "\tmuonMultiplicity : " << (int)j.muonMultiplicity << std::endl;
  std::cout << "\tHFHadronMultiplicity : " << (int)j.HFHadronMultiplicity << std::endl;
  std::cout << "\tHFEMMultiplicity : " << (int)j.HFEMMultiplicity << std::endl;
  std::cout << "\tchargedMultiplicity : " << (int)j.chargedMultiplicity << std::endl;
  std::cout << "\tneutralMultiplicity : " << (int)j.neutralMultiplicity << std::endl;
  std::cout << "\tvertex : "; Print(j.vertex);
  std::cout << "\tmomentum : "; Print(j.momentum);
  std::cout << "\tjecScaleFactors : ";
  for(std::map<TString, Float_t>::const_iterator it = j.jecScaleFactors.begin();
      it != j.jecScaleFactors.end(); it++) {
    std::cout << "(" << it->first << "," << it->second << ") ";
  }
  std::cout << std::endl;
  std::cout << "\tbTagDiscriminators :";
  for(std::vector<Float_t>::const_iterator it = j.bTagDiscriminators.begin();
      it != j.bTagDiscriminators.end(); it++) {
    std::cout << " " << *it;
  }
  std::cout << std::endl;
}



void Print(const susy::Particle& p) {
  std::cout << "\tstatus : " << (int)p.status << std::endl;
  std::cout << "\tmotherIndex : " << p.motherIndex << std::endl;
  std::cout << "\tpdgId : " << p.pdgId << std::endl;
  std::cout << "\tcharge : " << (int)p.charge << std::endl;
  std::cout << "\tvertex : "; Print(p.vertex);
  std::cout << "\tmomentum : "; Print(p.momentum);
  std::cout << std::endl;
}


void Print(const susy::Track& t) {
  std::cout << "\talgorithm : " << t.algorithm << std::endl;
  std::cout << "\tquality : " << t.quality << std::endl;
  std::cout << "\tnumberOfValidHits : " << (int)t.numberOfValidHits << std::endl;
  std::cout << "\tnumberOfValidTrackerHits : " << (int)t.numberOfValidTrackerHits << std::endl;
  std::cout << "\tnumberOfValidMuonHits : " << (int)t.numberOfValidMuonHits << std::endl;
  std::cout << "\tnumberOfValidPixelHits : " << (int)t.numberOfValidPixelHits << std::endl;
  std::cout << "\tnumberOfValidStripHits : " << (int)t.numberOfValidStripHits << std::endl;
  std::cout << "\tchi2 : " << t.chi2 << std::endl;
  std::cout << "\tndof : " << t.ndof << std::endl;
  std::cout << "\tcharge : " << t.charge << std::endl;
  std::cout << "\terror[5] : (";
  for(int i=0; i<5; i++) std::cout << t.error[i] << " ";
  std::cout << ")" << std::endl;
  std::cout << "\tvertex : "; Print(t.vertex);
  std::cout << "\tmomentum : "; Print(t.momentum);
  std::cout << "\textrapolatedPosition : ===> ";
  if(t.extrapolatedPositions.size() == 0) {
    std::cout << "(0)" << std::endl;
  }
  else {
    std::cout << std::endl;
    for(std::map<TString,TVector3>::const_iterator it = t.extrapolatedPositions.begin();
	it != t.extrapolatedPositions.end(); it++) {
      std::cout << "\t\t[ " << it->first << " : "; Print(it->second); std::cout << " ]" << std::endl;
    }
  }
}


void Print(const susy::PFParticle& p) {
  std::cout << "\tpdgId : " << p.pdgId << std::endl;
  std::cout << "\tcharge : " << (int)p.charge << std::endl;
  std::cout << "\tecalEnergy : " << p.ecalEnergy << std::endl;
  std::cout << "\trawEcalEnergy : " << p.rawEcalEnergy << std::endl;
  std::cout << "\thcalEnergy : " << p.hcalEnergy << std::endl;
  std::cout << "\trawHcalEnergy : " << p.rawHcalEnergy << std::endl;
  std::cout << "\tpS1Energy : " << p.pS1Energy << std::endl;
  std::cout << "\tpS2Energy : " << p.pS2Energy << std::endl;
  std::cout << "\tvertex : "; Print(p.vertex);
  std::cout << "\tpositionAtECALEntrance : "; Print(p.positionAtECALEntrance);
  std::cout << "\tmomentum : "; Print(p.momentum);
  std::cout << std::endl;
}


void Print(const susy::PUSummaryInfo& p) {
  std::cout << "\tnumInteractions : " << p.numInteractions << std::endl;
  std::cout << "\tBX : " << p.BX << std::endl;
  std::cout << "\ttrueNumInteractions : " << p.trueNumInteractions << std::endl;
  std::cout << "\tzPositions :";
  for(std::vector<float>::const_iterator it = p.zPositions.begin(); it != p.zPositions.end(); it++) {
    std::cout << " " << *it; 
  }
  std::cout << std::endl;
  std::cout << "\tsumPTLowPT :";
  for(std::vector<float>::const_iterator it = p.sumPTLowPT.begin(); it != p.sumPTLowPT.end(); it++) {
    std::cout << " " << *it; 
  }
  std::cout << std::endl;
  std::cout << "\tsumPTHighPT :";
  for(std::vector<float>::const_iterator it = p.sumPTHighPT.begin(); it != p.sumPTHighPT.end(); it++) {
    std::cout << " " << *it; 
  }
  std::cout << std::endl;
  std::cout << "\tnumTracksLowPT :";
  for(std::vector<int>::const_iterator it = p.numTracksLowPT.begin(); it != p.numTracksLowPT.end(); it++) {
    std::cout << " " << *it; 
  }
  std::cout << std::endl;
  std::cout << "\tnumTracksHighPT :";
  for(std::vector<int>::const_iterator it = p.numTracksHighPT.begin(); it != p.numTracksHighPT.end(); it++) {
    std::cout << " " << *it; 
  }
  std::cout << std::endl;
  std::cout << "\tinstLumi :";
  for(std::vector<float>::const_iterator it = p.instLumi.begin(); it != p.instLumi.end(); it++) {
    std::cout << " " << *it; 
  }
  std::cout << std::endl;
  std::cout << "\tdataMixerRun :";
  for(std::vector<unsigned int>::const_iterator it = p.dataMixerRun.begin(); it != p.dataMixerRun.end(); it++) {
    std::cout << " " << *it; 
  }
  std::cout << std::endl;
  std::cout << "\tdataMixerEvt :";
  for(std::vector<unsigned int>::const_iterator it = p.dataMixerEvt.begin(); it != p.dataMixerEvt.end(); it++) {
    std::cout << " " << *it; 
  }
  std::cout << std::endl;
  std::cout << "\tdataMixerLumiSection :";
  for(std::vector<unsigned int>::const_iterator it = p.dataMixerLumiSection.begin(); it != p.dataMixerLumiSection.end(); it++) {
    std::cout << " " << *it; 
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
  std::cout << "rho : " << event.rho << std::endl;
  std::cout << "rhoBarrel : " << event.rhoBarrel << std::endl;
  std::cout << "metFilterBit : " << event.metFilterBit << std::endl;
  std::cout << "metFilterBit break down ===> ";
  std::cout << "passCSCBeamHalo(" << event.passCSCBeamHalo() << ") ";
  std::cout << "passHcalNoise(" << event.passHcalNoise() << ") ";
  std::cout << "passEcalDeadCellTP(" << event.passEcalDeadCellTP() << ") ";
  std::cout << "passEcalDeadCellBE(" << event.passEcalDeadCellBE() << ") ";
  std::cout << "passHcalLaser(" << event.passHcalLaser() << ") ";
  std::cout << "passTrackingFailure(" << event.passTrackingFailure() << ") ";
  std::cout << "passMetFilters(" << event.passMetFilters() << ") ";
  std::cout << std::endl;

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
  for(std::vector<susy::Vertex>::const_iterator it = event.vertices.begin(); it != event.vertices.end(); it++) {
    Print(*it);
  }
  std::cout << std::endl;

  std::cout << "tracks size(" << event.tracks.size() << ") =========>" << std::endl;
  for(std::vector<susy::Track>::const_iterator it = event.tracks.begin(); it != event.tracks.end(); it++) {
    Print(*it);
  }
  std::cout << std::endl;

  std::cout << "superClusters size(" << event.superClusters.size() << ") =========>" << std::endl;
  std::cout << "clusters size(" << event.clusters.size() << ") =========>" << std::endl;
  std::cout << "muons size(" << event.muons.size() << ") =========>" << std::endl;
  std::cout << "electrons size(" << event.electrons.size() << ") =========>" << std::endl;

  std::cout << "photons size(" << event.photons.size() << ") =========>" << std::endl;
  for(std::map<TString,susy::PhotonCollection>::const_iterator it = event.photons.begin(); it != event.photons.end(); it++) {
    std::cout << it->first << " size(" << it->second.size() << ") ======>" << std::endl;
    for(std::vector<susy::Photon>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
      Print(*it2);
    }
  }
  std::cout << std::endl;

  std::cout << "caloJets size(" << event.caloJets.size() << ") =========>" << std::endl;
  for(std::map<TString, susy::CaloJetCollection>::const_iterator it = event.caloJets.begin(); it != event.caloJets.end(); it++) {
    std::cout << it->first << " size(" << it->second.size() << ") ======>" << std::endl;
    for(std::vector<susy::CaloJet>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
      Print(*it2);
    }
  }
  std::cout << std::endl;

  std::cout << "pfJets size(" << event.pfJets.size() << ") =========>" << std::endl;
  for(std::map<TString, susy::PFJetCollection>::const_iterator it = event.pfJets.begin(); it != event.pfJets.end(); it++) {
    std::cout << it->first << " size(" << it->second.size() << ") ======>" << std::endl;
    for(std::vector<susy::PFJet>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
      Print(*it2);
    }
  }
  std::cout << std::endl;


  std::cout << "jptJets size(" << event.jptJets.size() << ") =========>" << std::endl;

  std::cout << "pfParticles size(" << event.pfParticles.size() << ") =========>" << std::endl;
  for(std::map<TString, susy::PFParticleCollection>::const_iterator it = event.pfParticles.begin(); it != event.pfParticles.end(); it++) {
    std::cout << it->first << " size(" << it->second.size() << ") ======>" << std::endl;
    for(std::vector<susy::PFParticle>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
      Print(*it2);
    }
  }
  std::cout << std::endl;

  std::cout << "generalTracks size(" << event.generalTracks.size() << ") =========>" << std::endl;

  if(!event.isRealData) {
    std::cout << "pu size(" << event.pu.size() << ") =========>" << std::endl;
    for(susy::PUSummaryInfoCollection::const_iterator it = event.pu.begin(); it != event.pu.end(); it++) {
      std::cout << "\tpu : "; Print(*it);
    }
    std::cout << std::endl;

    std::cout << "simVertices size(" << event.simVertices.size() << ") =========>" << std::endl;
    for(std::vector<TVector3>::const_iterator it = event.simVertices.begin(); it != event.simVertices.end(); it++) {
      std::cout << "\tsimVertex : "; Print(*it);
    }
    std::cout << std::endl;

    std::cout << "genParticles size(" << event.genParticles.size() << ") =========>" << std::endl;
    for(std::vector<susy::Particle>::const_iterator it = event.genParticles.begin(); it != event.genParticles.end(); it++) {
      Print(*it);
    }
    std::cout << std::endl;

    std::cout << "gridParams size(" << event.gridParams.size() << ") =========>" << std::endl;
    std::cout << "\t";
    for(std::map<TString,Float_t>::const_iterator it = event.gridParams.begin(); it != event.gridParams.end(); it++) {
      std::cout << "(" << it->first << "," << it->second << ") ";
    }
    std::cout << std::endl;
  } // if(!event.isRealData)

  std::cout << std::endl;
}
