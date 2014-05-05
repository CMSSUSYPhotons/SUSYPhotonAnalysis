#ifndef SusyTriggerEvent_h
#define SusyTriggerEvent_h

#include <vector>
#include <map>
#include <string>

#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"

namespace susy {

  struct TriggerObject {
    Float_t pt;
    Float_t eta;
    Float_t phi;
    Float_t mass;
    Short_t vid;
    TriggerObject() : pt(0.), eta(0.), phi(0.), mass(0.), vid(0) {}
    TriggerObject(Float_t, Float_t, Float_t, Float_t, Short_t = 0);
    TriggerObject& operator=(TriggerObject const& _rhs) { pt = _rhs.pt; eta = _rhs.eta; phi = _rhs.phi; mass = _rhs.mass; vid = _rhs.vid; return *this; }
    TLorentzVector lorentzVector() const;
    Double_t deltaR(TLorentzVector const&) const;
  };

  typedef std::vector<TriggerObject> TriggerObjectCollection;

  class TriggerEvent {
  public:
    TriggerEvent();
    ~TriggerEvent();

    Bool_t bookTrees(TString const&);
    Bool_t bookTrees(TFile*);
    void initializeEvent(UInt_t, UInt_t);
    void fillObject(TriggerObject const&);
    void fillFilter(TString const&, std::vector<Int_t> const&, std::vector<UShort_t> const&, std::map<UShort_t, UShort_t> const&);
    void fillEvent();
    void write();
    void copyEvent(TriggerEvent&);

    void reset();

    Bool_t bindTree(TTree*, TString const&, TString const&);
    Bool_t bindTree(TTree*, TTree*, Bool_t = kTRUE);
    TriggerObjectCollection getFilterObjects(TString const&, Short_t = 0);
    Bool_t hasMatchingObject(TString const&, TLorentzVector const&, Double_t = 0.1, Short_t = 0);

    UInt_t getRunNumber() const { return event_.runNumber; }
    UInt_t getEventNumber() const { return event_.eventNumber; }

    void setVerbosity(Int_t _v) { verbosity_ = _v; }

  private:
    Bool_t loadTrees();

    struct EventData {
      Long64_t filterBegin;
      Long64_t filterEnd;
      Long64_t objectBegin;
      Long64_t objectEnd;
      UInt_t runNumber;
      UInt_t eventNumber;
    } event_;
    struct FilterData {
      Long64_t filterObjectBegin;
      Long64_t filterObjectEnd;
      UShort_t id; // within the current file
    } filter_;
    struct FilterObject {
      Short_t vid; // particle ID assigned within HLT
      UShort_t key; // index of the object in event
      //e.g. index of object 1 of filter A: event.objectBegin + filterObject[filter[event.filterBegin + offsetA].filterObjectBegin + 1].key)
    } filterObject_;
    TriggerObject object_;

    Bool_t readMode_;
    Bool_t writeMode_;
    Bool_t ownEventTree_;
    Bool_t ordered_;
    std::map<TString, UShort_t> filterIdMap_;
    TTree* eventTree_;
    TTree* filterTree_;
    TTree* filterObjectTree_;
    TTree* objectTree_;
    TTree* susyTree_;
    Int_t currentTreeNumber_;
    Int_t verbosity_;
  };

}

#endif
