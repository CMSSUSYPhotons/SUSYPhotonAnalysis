
void ana(TString ds="relval", TString physics="ttbar") {

  gSystem->Load("libSusyEvent.so");

  gROOT->LoadMacro("SusyEventAnalyzer.cc+");

  TChain* chain = new TChain("susyTree");
  chain->Add("../susyEvents.root");

  SusyEventAnalyzer* sea = new SusyEventAnalyzer(chain);

  // configuration parameters
  // any values given here will replace the default values
  sea->SetDataset(physics+"_"+ds);                    // dataset name
  sea->SetPrintInterval(1e4);             // print frequency
  sea->SetPrintLevel(10);                  // print level for event contents
  sea->SetUseTrigger(false);
  //sea->SetHltName("HLT_Ele15_LW_L1R");    // HLT trigger path name
  sea->SetFilter(false);                  // filter events passing final cuts
  sea->SetProcessNEvents(-1);             // number of events to be processed

  TStopwatch ts;

  ts.Start();

  sea->Loop();

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
