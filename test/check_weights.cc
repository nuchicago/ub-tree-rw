#include <iostream>
#include <string>
#include <vector>
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "uboone/EventWeight/MCEventWeight.h"

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " file1.root file2.root" << std::endl;
    return 0;
  }

  std::vector<std::string> f1 = { argv[1] };
  std::vector<std::string> f2 = { argv[2] };

  gallery::Event ev1(f1);
  gallery::Event ev2(f2);
  assert(ev1.numberOfEventsInFile() == ev2.numberOfEventsInFile());

  art::InputTag mlabel = { "mcweight" };

  while (!ev1.atEnd() && !ev2.atEnd()) {
    auto const& wlist1 = *ev1.getValidHandle<std::vector<evwgh::MCEventWeight> >(mlabel);
    auto const& wlist2 = *ev2.getValidHandle<std::vector<evwgh::MCEventWeight> >(mlabel);
    auto const& w1 = wlist1.at(0).fWeight;
    auto const& w2 = wlist2.at(0).fWeight;

    assert(w1.size() == w2.size());

    for (auto it : w1) {
      std::vector<double> f1 = it.second;
      std::vector<double> f2 = w2.find(it.first)->second;
      assert(f1.size() == f2.size());
      for (size_t i=0; i<f1.size(); i++) {
        if (f1[i] != f2[i]) {
          std::cerr << "Mismatch in weight " << it.first << ": "
                    << f1[i] << " != " << f2[i] << std::endl;
        }
        assert(f1[i] == f2[i]);
      }
    }

    ev1.next();
    ev2.next();
  }

  return 0;
}

