#ifndef EVENT_WRITER_H
#define EVENT_WRITER_H

#include "Event.h"
#include <TFile.h>
#include <TTree.h>
#include <string>

class EventWriter {
public:
  EventWriter(const std::string& filename);
  ~EventWriter();

  void WriteEvent(Event* evt);
  void Close();

private:
  TFile* m_file = nullptr;
  TTree* m_tree = nullptr;
  Event* m_evtBuffer = nullptr;
  int m_eventIndex = 0;
};

#endif // EVENT_WRITER_H