#include <algorithm>
#include <iostream>

#include "core/Event.h"
#include "io/EventReaderAMPT.h"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <AMPT_root_file or a list>\n";
    return 1;
  }
  EventReaderAMPT reader(argv[1]);

  for (int ievt = 0; ievt < 3; ++ievt) {
    Event evt;
    if (!reader.NextEvent(evt)) {
      std::cerr << "End of file at event " << ievt << "\n";
      break;
    }
    std::cout << "Event #" << ievt << " parton count = " << evt.GetPartons().size() << "\n";
    for (size_t ip = 0; ip < std::min<size_t>(evt.GetPartons().size(), 5); ++ip) {
      auto p = evt.GetPartons()[ip];
      auto pos = p->GetPosition();
      std::cout << "  [" << ip << "] pos=(" << pos[0] << "," << pos[1] << "," << pos[2] << ")"
                << " B=" << p->GetBaryonNumber() << "\n";
    }
  }
  return 0;
}
