// Original Author:  Dongwook Jang
// $Id: ana.C,v 1.8 2011/11/01 22:14:51 dwjang Exp $

void ana(TString outputName="analysis"){

  gSystem->Load("libSusyEvent.so");

  // Printing utility for ntuple variables
  gROOT->LoadMacro("SusyEventPrinter.cc+");

  // Main analysis code
  gROOT->LoadMacro("SusyEventAnalyzer.cc+");

  // chain of inputs
  TChain chain("susyTree");
  chain.Add("susyEvents.root");
  //chain->Add("dcap:///pnfs/cms/WAX/resilient/lpcpjm/SusyNtuples/cms423v2_v1/Run2011A-May10ReReco-v1/Photon/susyEvent_1_1_dLs.root");

  SusyEventAnalyzer sea(chain);

  // configuration parameters
  // any values given here will replace the default values
  sea.SetOutput(outputName);
  sea.SetPrintInterval(10000);             // print frequency
  sea.SetPrintLevel(1);                  // print level for event contents
  sea.AddHltName("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50");  // add HLT trigger path name
  sea.CopyEvents(false);                  // filter events passing final cuts
  sea.SetProcessNEvents(-1);             // number of events to be processed

  // as an example -- add the full path to your favorite Json here.  More than one can be "Include"ed
  //  sea.IncludeAJson("Cert_161079-161352_7TeV_PromptReco_Collisions11_JSON_noESpbl_v2.txt");
  //sea.IncludeAJson("anotherJSON.txt");

  TStopwatch ts;

  ts.Start();

  sea.Run();

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
