#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "SusyTriggerEvent.h"

#include "TString.h"

#include <set>
#include <map>

class SusyTriggerNtuplizer : public edm::EDAnalyzer {
public:
  SusyTriggerNtuplizer(edm::ParameterSet const&);
  ~SusyTriggerNtuplizer() {}

private:
  void endJob();
  void analyze(edm::Event const&, edm::EventSetup const&);

  edm::InputTag triggerEventTag_;
  std::set<std::string> filterTags_;
  susy::TriggerEvent triggerEvent_;
};

SusyTriggerNtuplizer::SusyTriggerNtuplizer(edm::ParameterSet const& _ps) :
  triggerEventTag_(_ps.getUntrackedParameter<edm::InputTag>("triggerEventTag")),
  filterTags_(),
  triggerEvent_()
{
  if(_ps.existsAs<std::vector<std::string> >("filterTags", false)){
    std::vector<std::string> filterTags(_ps.getUntrackedParameter<std::vector<std::string> >("filterTags"));
    filterTags_.insert(filterTags.begin(), filterTags.end());
  }

  TString fileName(_ps.getUntrackedParameter<std::string>("fileName"));
  if(!triggerEvent_.bookTrees(fileName))
    throw cms::Exception("IOError") << "Could not book trees";
}

void
SusyTriggerNtuplizer::analyze(edm::Event const& _evt, edm::EventSetup const&)
{
  triggerEvent_.initializeEvent(_evt.id().run(), _evt.id().event());

  try{
    edm::Handle<trigger::TriggerEvent> teH;

    trigger::TriggerObjectCollection const& objects(teH->getObjects());

    std::map<unsigned short, unsigned short> keyMap;

    unsigned nF(teH->sizeFilters());
    for(unsigned iF(0); iF < nF; ++iF){
      if(filterTags_.size() != 0 && filterTags_.find(teH->filterTag(iF).label()) == filterTags_.end()) continue;

      trigger::Keys const& keys(teH->filterKeys(iF));

      for(unsigned iK(0); iK != keys.size(); ++iK){
        unsigned short key(keys[iK]);

        if(keyMap.find(key) != keyMap.end()) continue;

        trigger::TriggerObject const& obj(objects.at(key));
        triggerEvent_.fillObject(susy::TriggerObject(obj.pt(), obj.eta(), obj.phi(), obj.mass()));

        keyMap[key] = keyMap.size() - 1;
      }

      triggerEvent_.fillFilter(teH->filterTag(iF).label(), teH->filterIds(iF), keys, keyMap);
    }
  }
  catch(...){
    triggerEvent_.fillEvent();
    triggerEvent_.write();
    throw;
  }

  triggerEvent_.fillEvent();
}

void
SusyTriggerNtuplizer::endJob()
{
  triggerEvent_.write();
}

DEFINE_FWK_MODULE(SusyTriggerNtuplizer);
