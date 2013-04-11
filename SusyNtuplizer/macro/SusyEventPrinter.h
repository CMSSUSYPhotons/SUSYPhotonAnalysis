#ifndef SusyEventPrinter_H
#define SusyEventPrinter_H

#include "../src/SusyEvent.h"

void Print(const TVector2& p);
void Print(const TVector3& p);
void Print(const TLorentzVector& p);
void Print(const susy::TriggerMap& t);
void Print(const susy::MET& met);
void Print(const susy::Vertex& vtx);
void Print(const susy::Photon& p);
void Print(const susy::CaloJet& j);
void Print(const susy::PFJet& j);
void Print(const susy::Particle& p);
void Print(const susy::Track& t);
void Print(const susy::PFParticle& p);

void Print(const susy::Event& event);

#endif
