/**
 * Arborist makes trees with weights from art ROOT.
 *
 * atm 2018/03/26
 */

#include <cassert>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TInterpreter.h>

#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "uboone/EventWeight/MCEventWeight.h"

void extract_tree(TString outfile, std::vector<std::string> filenames) {
  gInterpreter->GenerateDictionary("std::map<std::string, std::vector<double> >",
                                   "string;map;vector");

  // Output tree
  TFile fout(outfile, "recreate");
  TTree* t = new TTree("events", "");
  int eventID;
  double enu;
  std::map<std::string, std::vector<double> > weights;
  t->Branch("eventID", &eventID);
  t->Branch("Enu", &enu);
  t->Branch("weights", &weights);

  art::InputTag mct_tag { "TreeReader" };
  art::InputTag mcwgh_tag { "mcweight" };
  eventID = 0;

  std::cout << "Copying..." << std::endl;
  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {
    auto const& mct = *ev.getValidHandle<std::vector<simb::MCTruth> >(mct_tag);
    auto const& mcwgh = *ev.getValidHandle<std::vector<evwgh::MCEventWeight> >(mcwgh_tag);
    assert(mcwgh.size() == 1 && mct.size() == 1);

    enu = mct[0].GetNeutrino().Nu().E();
    weights = mcwgh[0].fWeight;

    t->Fill();

    eventID++;
  }

  std::cout << "Writing..." << std::endl;
  fout.cd();
  t->Write();
  fout.Close();

  std::cout << "Done." << std::endl;
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0]
              << " output.root input.root [input2.root ...]" << std::endl;
    return 0;
  }

  std::vector<std::string> filenames;
  for (int i=2; i<argc; i++) { 
    std::cout << "FILE: " << argv[i] << std::endl; 
    filenames.push_back(argv[i]);
  }

  extract_tree(argv[1], filenames);

  return 0;
}

