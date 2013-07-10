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

class SusyTriggerNtuplizer : public edm::EDAnalyzer {
public:
  SusyTriggerNtuplizer(edm::ParameterSet const&);
  ~SusyTriggerNtuplizer() {}

private:
  void endJob();
  void analyze(edm::Event const&, edm::EventSetup const&);

  edm::InputTag triggerEventTag_;
  susy::TriggerEvent triggerEvent_;
};

SusyTriggerNtuplizer::SusyTriggerNtuplizer(edm::ParameterSet const& _ps) :
  triggerEventTag_(_ps.getUntrackedParameter<edm::InputTag>("triggerEventTag")),
  triggerEvent_()
{
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
    unsigned nObj(objects.size());
    for(unsigned iObj(0); iObj < nObj; ++iObj){
      trigger::TriggerObject const& obj(objects.at(iObj));
      triggerEvent_.fillObject(susy::TriggerObject(obj.pt(), obj.eta(), obj.phi(), obj.mass()));
    }

    unsigned nF(teH->sizeFilters());
    for(unsigned iF(0); iF < nF; ++iF){
      trigger::Vids const& vids(teH->filterIds(iF));
      trigger::Keys const& keys(teH->filterKeys(iF));
      triggerEvent_.fillFilter(teH->filterTag(iF).label(), vids, keys);
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
