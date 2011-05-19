#ifndef SusyEventPrinter_H
#define SusyEventPrinter_H

#include "../src/SusyEvent.h"

void Print(const TVector2& p);
void Print(const TVector3& p);
void Print(const TLorentzVector& p);
void Print(const susy::TriggerMap& t);
void Print(const susy::CorrMETData& cm);
void Print(const susy::MET& met);


void Print(const susy::Event& event);

#endif
