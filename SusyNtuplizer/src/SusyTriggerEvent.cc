#include "SusyTriggerEvent.h"

#include "TFile.h"
#include "TMath.h"
#include "TVector2.h"
#include "TBranch.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TCollection.h"
#include "TLeaf.h"
#include "TSystem.h"

#include <iostream>

namespace susy {

  TriggerObject::TriggerObject(Float_t _pt, Float_t _eta, Float_t _phi, Float_t _mass, Short_t _vid/* = 0*/) :
    pt(_pt),
    eta(_eta),
    phi(_phi),
    mass(_mass),
    vid(_vid)
  {
  }

  TLorentzVector
  TriggerObject::lorentzVector() const
  {
    Double_t theta(2. * TMath::ATan(TMath::Exp(-eta)));
    Double_t px(pt * TMath::Cos(phi));
    Double_t py(pt * TMath::Sin(phi));
    Double_t pz(pt / TMath::Tan(theta));
    Double_t e(TMath::Sqrt(px * px + py * py + pz * pz + mass * mass));

    return TLorentzVector(px, py, pz, e);
  }

  Double_t
  TriggerObject::deltaR(TLorentzVector const& _lv) const
  {
    Double_t dEta(_lv.Eta() - eta);
    Double_t dPhi(TVector2::Phi_mpi_pi(_lv.Phi() - phi));

    return TMath::Sqrt(dEta * dEta + dPhi * dPhi);
  }

  TriggerEvent::TriggerEvent() :
    event_(),
    filter_(),
    filterObject_(),
    object_(),
    readMode_(kFALSE),
    writeMode_(kFALSE),
    ownEventTree_(kFALSE),
    ordered_(kTRUE),
    filterIdMap_(),
    eventTree_(0),
    filterTree_(0),
    filterObjectTree_(0),
    objectTree_(0),
    susyTree_(0),
    currentTreeNumber_(-1),
    verbosity_(0)
  {
  }

  TriggerEvent::~TriggerEvent()
  {
    reset();
  }

  Bool_t
  TriggerEvent::bookTrees(TString const& _fileName)
  {
    if(readMode_ || writeMode_){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::: TriggerEvent is already initialized" << std::endl;
      return kFALSE;
    }

    TFile* file(TFile::Open(_fileName, "recreate"));

    if(!file || file->IsZombie()){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::: Cannot create " << _fileName << std::endl;
      delete file;
      return kFALSE;
    }

    ownEventTree_ = kTRUE;

    return bookTrees(file);
  }

  Bool_t
  TriggerEvent::bookTrees(TFile* _file)
  {
    if(readMode_ || writeMode_){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::bookTrees: TriggerEvent is already initialized" << std::endl;
      return kFALSE;
    }
    if(!_file || _file->IsZombie() || !(TString(_file->GetName()).Contains(".root"))){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::bookTrees: Cannot open file" << std::endl;
      return kFALSE;
    }

    writeMode_ = kTRUE;

    filterIdMap_.clear();

    _file->cd();

    eventTree_ = new TTree("triggerEventTree", "SUSY TriggerEvent");
    eventTree_->SetAutoSave(10000000);
    eventTree_->Branch("runNumber", &event_.runNumber, "runNumber/i");
    eventTree_->Branch("eventNumber", &event_.eventNumber, "eventNumber/i");
    eventTree_->Branch("filterIndex", &event_.filterBegin, "begin/L:end");
    eventTree_->Branch("objectIndex", &event_.objectBegin, "begin/L:end");

    TFile::Open(TString(_file->GetName()).ReplaceAll(".root", "_tmp_filter.root"), "recreate");
    filterTree_ = new TTree("filterTree", "SUSY Trigger Filters");
    filterTree_->Branch("id", &filter_.id, "id/s");
    filterTree_->Branch("filterObjectIndex", &filter_.filterObjectBegin, "begin/L:end");

    TFile::Open(TString(_file->GetName()).ReplaceAll(".root", "_tmp_filterObject.root"), "recreate");
    filterObjectTree_ = new TTree("filterObjectTree", "SUSY Filter Objects");
    filterObjectTree_->Branch("filterObject", &filterObject_.vid, "vid/S:key/s");

    TFile::Open(TString(_file->GetName()).ReplaceAll(".root", "_tmp_object.root"), "recreate");
    objectTree_ = new TTree("objectTree", "SUSY Trigger Objects");
    objectTree_->Branch("triggerObject", &object_.pt, "pt/F:eta:phi:mass");

    currentTreeNumber_ = 0;

    return kTRUE;
  }

  void
  TriggerEvent::initializeEvent(UInt_t _runNumber, UInt_t _eventNumber)
  {
    if(!writeMode_){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::initializeEvent: Not in write mode" << std::endl;
      return;
    }

    event_.runNumber = _runNumber;
    event_.eventNumber = _eventNumber;

    event_.filterBegin = filterTree_->GetEntries();
    event_.objectBegin = objectTree_->GetEntries();

    if(verbosity_ > 1) std::cerr << "susy::TriggerEvent::initializeEvent: Initialized event info for Run " << event_.runNumber << " Event " << event_.eventNumber << std::endl;
  }

  void
  TriggerEvent::fillObject(TriggerObject const& _obj)
  {
    // must be called right after fillFilter
    // better implementationally if this was a part of fillFilter, but we don't want CMSSW dependence in this code

    if(!writeMode_){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::fillFilter: Not in write mode" << std::endl;
      return;
    }

    object_ = _obj;
    objectTree_->Fill();
  }

  void
  TriggerEvent::fillFilter(TString const& _filterName, std::vector<Int_t> const& _vids, std::vector<UShort_t> const& _keys, std::map<UShort_t, UShort_t> const& _keyMap)
  {
    if(!writeMode_){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::fillFilter: Not in write mode" << std::endl;
      return;
    }

    UInt_t nVid(_vids.size());
    if(nVid == 0){
      if(verbosity_ > 1) std::cerr << "susy::TriggerEvent::fillFilter: Vid vector size is 0 for filter " << _filterName << std::endl;
      return;
    }

    std::map<TString, UShort_t>::iterator fItr(filterIdMap_.find(_filterName));
    if(fItr == filterIdMap_.end()){
      filter_.id = filterIdMap_.size();
      filterIdMap_[_filterName] = filter_.id;
    }
    else
      filter_.id = fItr->second;

    filter_.filterObjectBegin = filterObjectTree_->GetEntries();
    filter_.filterObjectEnd = filterObjectTree_->GetEntries() + nVid;

    filterTree_->Fill();

    // key must exist in keymap

    for(UInt_t iVid(0); iVid != nVid; ++iVid){
      filterObject_.vid = _vids[iVid];
      filterObject_.key = _keyMap.find(_keys[iVid])->second;
      filterObjectTree_->Fill();
    }

    if(verbosity_ > 2){
      std::cerr << "susy::TriggerEvent::fillFilter: Filled filter info for " << _filterName << std::endl
                << " ID " << filter_.id << " FilterObjectIndices " << filter_.filterObjectBegin << "-" << filter_.filterObjectEnd << std::endl;
      if(verbosity_ > 3){
        std::cerr << " FilterObjects:" << std::endl;
        for(UInt_t iVid(0); iVid < nVid; ++iVid)
          std::cerr << "  vid=" << _vids[iVid] << ", key=" << _keyMap.find(_keys[iVid])->second << std::endl;
      }
    }
  }

  void
  TriggerEvent::fillEvent()
  {
    if(!writeMode_){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::fillEvent: Not in write mode" << std::endl;
      return;
    }

    event_.objectEnd = objectTree_->GetEntries();
    event_.filterEnd = filterTree_->GetEntries();

    eventTree_->Fill();

    if(verbosity_ > 1){
      std::cerr << "susy::TriggerEvent::fillEvent: Finalized event info for Run " << event_.runNumber << " Event " << event_.eventNumber << std::endl
                << " FilterIndices " << event_.filterBegin << "-" << event_.filterEnd << std::endl;
    }
  }

  void
  TriggerEvent::write()
  {
    if(!writeMode_){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::write: Not in write mode" << std::endl;
      return;
    }

    if(verbosity_ > 1) std::cerr << "susy::TriggerEvent::write: Writing triggerEvents" << std::endl;

    TFile* eventFile(eventTree_->GetCurrentFile());
    eventFile->cd();
    eventTree_->Write();
    // eventTree and the file will be deleted in reset() at the end of this function

    if(verbosity_ > 1) std::cerr << "susy::TriggerEvent::write: Writing filters" << std::endl;

    TFile* filterFile(filterTree_->GetCurrentFile());
    TString filterFileName(filterFile->GetName());
    filterFile->cd();
    filterTree_->Write();
    delete filterFile;

    filterFile = TFile::Open(filterFileName);
    filterFile->GetObject("filterTree", filterTree_);
    if(!filterTree_){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::write: [ERROR] Failed writing filters" << std::endl;
      reset();
      return;
    }

    eventFile->cd();
    TTree* filterTree(filterTree_->CloneTree(-1, "fast"));
    if(!filterTree){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::write: [ERROR] Failed writing filters" << std::endl;
      reset();
      return;
    }
    filterTree->Write();
    delete filterTree;

    delete filterFile;
    filterTree_ = 0;
    if(gSystem->Unlink(filterFileName) != 0){
      TFile* tmp(TFile::Open(filterFileName, "recreate"));
      delete tmp;
    }

    if(verbosity_ > 1) std::cerr << "susy::TriggerEvent::write: Writing filterObjects" << std::endl;

    TFile* filterObjectFile(filterObjectTree_->GetCurrentFile());
    TString filterObjectFileName(filterObjectFile->GetName());
    filterObjectFile->cd();
    filterObjectTree_->Write();
    delete filterObjectFile;

    filterObjectFile = TFile::Open(filterObjectFileName);
    filterObjectFile->GetObject("filterObjectTree", filterObjectTree_);
    if(!filterObjectTree_){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::write: [ERROR] Failed writing filterObjects" << std::endl;
      reset();
      return;
    }

    eventFile->cd();
    TTree* filterObjectTree(filterObjectTree_->CloneTree(-1, "fast"));
    if(!filterObjectTree){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::write: [ERROR] Failed writing filterObjects" << std::endl;
      reset();
      return;
    }
    filterObjectTree->Write();
    delete filterObjectTree;

    delete filterObjectFile;
    filterObjectTree_ = 0;
    if(gSystem->Unlink(filterObjectFileName) != 0){
      TFile* tmp(TFile::Open(filterObjectFileName, "recreate"));
      delete tmp;
    }

    if(verbosity_ > 1) std::cerr << "susy::TriggerEvent::write: Writing objects" << std::endl;

    TFile* objectFile(objectTree_->GetCurrentFile());
    TString objectFileName(objectFile->GetName());
    objectFile->cd();
    objectTree_->Write();
    delete objectFile;

    objectFile = TFile::Open(objectFileName);
    objectFile->GetObject("objectTree", objectTree_);
    if(!objectTree_){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::write: [ERROR] Failed writing objects" << std::endl;
      reset();
      return;
    }

    eventFile->cd();
    TTree* objectTree(objectTree_->CloneTree(-1, "fast"));
    if(!objectTree){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::write: [ERROR] Failed writing objects" << std::endl;
      reset();
      return;
    }
    objectTree->Write();
    delete objectTree;

    delete objectFile;
    objectTree_ = 0;
    if(gSystem->Unlink(objectFileName) != 0){
      TFile* tmp(TFile::Open(objectFileName, "recreate"));
      delete tmp;
    }

    if(verbosity_ > 1) std::cerr << "susy::TriggerEvent::write: Writing filter table" << std::endl;

    eventFile->cd();
    TString* filterName(new TString());
    TTree* filterTable(new TTree("filterTable", "SUSY Trigger Filter Table"));
    filterTable->Branch("filterName", "TString", &filterName);
    UShort_t nIds(filterIdMap_.size());
    for(UShort_t iId(0); iId < nIds; ++iId){
      std::map<TString, UShort_t>::iterator fEnd(filterIdMap_.end());
      for(std::map<TString, UShort_t>::iterator fItr(filterIdMap_.begin()); fItr != fEnd; ++fItr){
        if(fItr->second == iId){
          if(verbosity_ > 0) std::cerr << " " << iId << ". " << fItr->first << std::endl;

          *filterName = fItr->first;
          filterTable->Fill();
          break;
        }
      }
    }
    filterTable->Write();

    delete filterTable;
    delete filterName;

    if(verbosity_ > 1) std::cerr << "susy::TriggerEvent::write: Cleanup" << std::endl;

    reset();
  }

  void
  TriggerEvent::copyEvent(TriggerEvent& _orig)
  {
    if(!writeMode_){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::copyEvent: Not in write mode" << std::endl;
      return;
    }

    if(_orig.currentTreeNumber_ != _orig.eventTree_->GetTreeNumber() && !_orig.loadTrees()){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::copyEvent: File transition failed" << std::endl;
      return;
    }

    event_.runNumber = _orig.event_.runNumber;
    event_.eventNumber = _orig.event_.eventNumber;

    event_.objectBegin = objectTree_->GetEntries();

    Long64_t objectIndex(_orig.event_.objectBegin);
    while(objectIndex != _orig.event_.objectEnd && _orig.objectTree_->GetEntry(objectIndex++) != 0){
      object_ = _orig.object_;
      objectTree_->Fill();
    }

    event_.objectEnd = objectTree_->GetEntries();

    event_.filterBegin = filterTree_->GetEntries();

    Long64_t filterIndex(_orig.event_.filterBegin);
    std::map<TString, UShort_t>::const_iterator origMapEnd(_orig.filterIdMap_.end());
    while(filterIndex != _orig.event_.filterEnd && _orig.filterTree_->GetEntry(filterIndex++) != 0){
      for(std::map<TString, UShort_t>::const_iterator origMapItr(_orig.filterIdMap_.begin()); origMapItr != origMapEnd; ++origMapItr){
	if(origMapItr->second == _orig.filter_.id){
	  std::map<TString, UShort_t>::iterator mapItr(filterIdMap_.find(origMapItr->first));
	  if(mapItr == filterIdMap_.end()){
	    filter_.id = filterIdMap_.size();
	    filterIdMap_[origMapItr->first] = filter_.id;
	  }
	  else
	    filter_.id = mapItr->second;

	  break;
	}
      }

      filter_.filterObjectBegin = filterObjectTree_->GetEntries();

      Long64_t filterObjectIndex(_orig.filter_.filterObjectBegin);
      while(filterObjectIndex != _orig.filter_.filterObjectEnd && _orig.filterObjectTree_->GetEntry(filterObjectIndex++) != 0){
	filterObject_ = _orig.filterObject_;
	filterObjectTree_->Fill();
      }

      filter_.filterObjectEnd = filterObjectTree_->GetEntries();

      filterTree_->Fill();
    }

    event_.filterEnd = filterTree_->GetEntries();

    eventTree_->Fill();
  }

  void
  TriggerEvent::reset()
  {
    if(susyTree_ && susyTree_->GetFriend("triggerEvent") == eventTree_){
      susyTree_->RemoveFriend(eventTree_);
      if(susyTree_->InheritsFrom(TChain::Class())) susyTree_->LoadTree(susyTree_->GetReadEntry());
    }

    susyTree_ = 0;
    ordered_ = kTRUE;
    currentTreeNumber_ = -1;

    if(ownEventTree_ && eventTree_){
      if(eventTree_->InheritsFrom(TChain::Class()))
        delete eventTree_;
      else
        delete eventTree_->GetCurrentFile();
    }
    eventTree_ = filterTree_ = filterObjectTree_ = objectTree_ = 0;

    filterIdMap_.clear();

    readMode_ = writeMode_ = ownEventTree_ = kFALSE;

    event_ = EventData();
    filter_ = FilterData();
    filterObject_ = FilterObject();
    object_ = TriggerObject();
  }

  Bool_t
  TriggerEvent::bindTree(TTree* _susyTree, TString const& _pattern, TString const& _replacement)
  {
    if(readMode_ || writeMode_){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::bindTree: TriggerEvent is already initialized" << std::endl;
      return kFALSE;
    }

    if(!_susyTree){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::bindTree: NULL pointer to susyTree" << std::endl;
      return kFALSE;
    }

    if(_susyTree->InheritsFrom(TChain::Class())){
      TChain* eventTree(new TChain("triggerEventTree"));

      TIter next(static_cast<TChain*>(_susyTree)->GetListOfFiles());
      TChainElement* ce(0);
      while((ce = static_cast<TChainElement*>(next()))){
        TString fileName(ce->GetTitle());
        fileName.ReplaceAll(_pattern, _replacement);
        if(eventTree->Add(fileName, 0) == 0){
          if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::bindTree: Could not open " << fileName << std::endl;
          delete eventTree;
          return kFALSE;
        }
      }
      ownEventTree_ = kTRUE;
      eventTree_ = eventTree;
    }
    else{
      TString fileName(_susyTree->GetCurrentFile()->GetName());
      fileName.ReplaceAll(_pattern, _replacement);
      TFile* source(TFile::Open(fileName));
      if(!source || source->IsZombie()){
        if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::bindTree: Could not open " << fileName << std::endl;
        delete source;
        return kFALSE;
      }
      source->GetObject("triggerEventTree", eventTree_);
      if(!eventTree_){
        if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::bindTree: triggerEventTree not found in " << fileName << std::endl;
        delete source;
        return kFALSE;
      }
    }

    return bindTree(_susyTree, eventTree_, kTRUE);
  }

  Bool_t
  TriggerEvent::bindTree(TTree* _susyTree, TTree* _triggerTree, Bool_t _ordered/* = kTRUE*/)
  {
    if(readMode_ || writeMode_){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::bindTree: TriggerEvent is already initialized" << std::endl;
      return kFALSE;
    }

    if(!_susyTree || !_triggerTree){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::bindTree: Cannot bind NULL pointers" << std::endl;
      return kFALSE;
    }

    susyTree_ = _susyTree;
    ordered_ = _ordered;

    eventTree_ = _triggerTree;

    if(ordered_){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::bindTree: Binding in synchronized mode" << std::endl;
      susyTree_->AddFriend(eventTree_, "triggerEvent");

      TObjArray* leaves(eventTree_->GetListOfLeaves());
      for(Int_t iL(0); iL != leaves->GetEntriesFast(); ++iL){
        TLeaf* leaf(static_cast<TLeaf*>(leaves->At(iL)));
        TString branchName(leaf->GetBranch()->GetName());
        susyTree_->SetBranchStatus("triggerEvent." + branchName, 1);
      }
    }
    else{
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::bindTree: Binding in unsynchronized mode" << std::endl;
      eventTree_->BuildIndex("runNumber", "eventNumber");
      eventTree_->SetBranchStatus("*", 1);
    }

    readMode_ = kTRUE;

    eventTree_->SetBranchAddress("runNumber", &event_.runNumber);
    eventTree_->SetBranchAddress("eventNumber", &event_.eventNumber);
    eventTree_->SetBranchAddress("filterIndex", &event_.filterBegin);
    eventTree_->SetBranchAddress("objectIndex", &event_.objectBegin);

    currentTreeNumber_ = -1;

    return kTRUE;
  }

  TriggerObjectCollection
  TriggerEvent::getFilterObjects(TString const& _filter, Short_t _vid/* = 0*/)
  {
    TriggerObjectCollection result;

    if(!readMode_){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::getFilterObjects: Not in read mode" << std::endl;
      return result;
    }

    if(!ordered_){
      Long64_t treeEntry(susyTree_->LoadTree(susyTree_->GetReadEntry()));
      UInt_t runNumber(susyTree_->GetBranch("runNumber")->GetLeaf("runNumber")->GetValue(treeEntry));
      UInt_t eventNumber(susyTree_->GetBranch("eventNumber")->GetLeaf("eventNumber")->GetValue(treeEntry));
      Int_t nBytes(eventTree_->GetEntryWithIndex(runNumber, eventNumber));

      if(nBytes == 0){
        if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::getFilterObjects: No event data found in nonsynchronized mode" << std::endl;
        return result;
      }
    }

    if(currentTreeNumber_ != eventTree_->GetTreeNumber() && !loadTrees()){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::getFilterObjects: File transition failed" << std::endl;
      return result;
    }

    if(_filter == "") return result;

    std::map<TString, UShort_t>::const_iterator idItr(filterIdMap_.find(_filter));
    if(idItr == filterIdMap_.end()){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::getFilterObjects: Filter " << _filter << " not found in table" << std::endl;
      return result;
    }

    UShort_t id(idItr->second);

    Long64_t objectOffset(event_.objectBegin);

    Long64_t filterIndex(event_.filterBegin);

    while(filterIndex != event_.filterEnd && filterTree_->GetEntry(filterIndex++) != 0){
      if(id != filter_.id) continue;

      if(verbosity_ > 1) std::cerr << "susy::TriggerEvent::getFilterObjects: Filter " << _filter << " found in the event " << event_.runNumber << "/" << event_.eventNumber << std::endl;

      Long64_t foIndex(filter_.filterObjectBegin);

      while(foIndex != filter_.filterObjectEnd && filterObjectTree_->GetEntry(foIndex++) != 0){
        if(_vid != 0 && filterObject_.vid != _vid) continue;

        if(verbosity_ > 2){
          std::cerr << " FilterObject vid=" << filterObject_.vid << " key=" << filterObject_.key;
          std::cerr.flush();
        }

        Long64_t objectIndex(objectOffset + filterObject_.key);

        if(objectIndex >= event_.objectEnd) continue;

        objectTree_->GetEntry(objectIndex);

        if(verbosity_ > 2) std::cerr << " pt=" << object_.pt << " eta=" << object_.eta << " phi=" << object_.phi << " mass=" << object_.mass << std::endl;

        object_.vid = filterObject_.vid;

        result.push_back(object_);
      }

      break;
    }

    return result;
  }     

  Bool_t
  TriggerEvent::hasMatchingObject(TString const& _filter, TLorentzVector const& _momentum, double _dR/* = 0.1*/, Short_t _vid/* = 0*/)
  {
    if(!readMode_){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::hasMatchingObject: Not in read mode" << std::endl;
      return kFALSE;
    }

    TriggerObjectCollection filterObjects(getFilterObjects(_filter, _vid));

    UInt_t nObj(filterObjects.size());
    for(UInt_t iObj(0); iObj < nObj; ++iObj)
      if(filterObjects[iObj].deltaR(_momentum) < _dR) return kTRUE;

    return kFALSE;
  }

  void
  TriggerEvent::Print(std::ostream& os/* = std::cout*/)
  {
    getFilterObjects(""); // initialize read

    os << "Trigger Event ======>" << std::endl;
    for(std::map<TString, UShort_t>::const_iterator fItr(filterIdMap_.begin()); fItr != filterIdMap_.end(); ++fItr){
      TriggerObjectCollection filterObjects(getFilterObjects(fItr->first));
      if(filterObjects.size() == 0) continue;
      os << "\t" << fItr->first << ": ";
      for(unsigned iO(0); iO != filterObjects.size(); ++iO){
        TriggerObject& obj(filterObjects[iO]);
        os << "(" << obj.pt << "," << obj.eta << "," << obj.phi << "," << obj.mass << ")[" << obj.vid << "] ";
      }
      os << std::endl;
    }
    os << std::endl;
  }

  Bool_t
  TriggerEvent::loadTrees()
  {
    TFile* file(eventTree_->GetCurrentFile());
    if(!file){
      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::loadTrees: Current file is NULL" << std::endl;
      return kFALSE;
    }

    filterIdMap_.clear();

    // local trees in the previous files are deleted at the file transition by ROOT
    filterTree_ = 0;
    filterObjectTree_ = 0;
    objectTree_ = 0;

    TTree* filterTable(0);
    file->GetObject("filterTable", filterTable);
    file->GetObject("filterTree", filterTree_);
    file->GetObject("filterObjectTree", filterObjectTree_);
    file->GetObject("objectTree", objectTree_);

    if(!filterTable || !filterTree_ || !filterObjectTree_ || !objectTree_){
      delete filterTable;
      delete filterTree_;
      delete filterObjectTree_;
      delete objectTree_;
      filterTree_ = 0;
      filterObjectTree_ = 0;
      objectTree_ = 0;

      if(verbosity_ > 0) std::cerr << "susy::TriggerEvent::loadTrees: Missing information on file" << std::endl;
      return kFALSE;
    }

    if(verbosity_ > 1) std::cerr << "susy::TriggerEvent::loadTrees: Loading filter table" << std::endl;

    TString* filterName(new TString());
    filterTable->SetBranchAddress("filterName", &filterName);

    UShort_t filterId(0);
    while(filterTable->GetEntry(filterId) != 0){
      if(verbosity_ > 1) std::cerr << " " << filterId << ". " << *filterName << std::endl;

      filterIdMap_[*filterName] = filterId++;
    }

    delete filterTable;
    delete filterName;

    filterTree_->SetBranchAddress("id", &filter_.id);
    filterTree_->SetBranchAddress("filterObjectIndex", &filter_.filterObjectBegin);

    filterObjectTree_->SetBranchAddress("filterObject", &filterObject_.vid);

    objectTree_->SetBranchAddress("triggerObject", &object_.pt);

    currentTreeNumber_ = eventTree_->GetTreeNumber();

    return kTRUE;
  }
}
