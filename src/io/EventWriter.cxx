#include "io/EventWriter.h"

EventWriter::EventWriter(const std::string& filename) {
  m_file = new TFile(filename.c_str(), "RECREATE");
  m_tree = new TTree("events", "Generated events");
  m_tree->Branch("event", "Event", &m_evtBuffer);
}

EventWriter::~EventWriter() {
  Close();  // 确保析构时关闭文件
}

void EventWriter::WriteEvent(Event* evt) {
  m_evtBuffer = evt;
  m_tree->Fill();
  ++m_eventIndex;
}

void EventWriter::Close() {
  if (m_file) {
    m_file->cd();
    m_tree->Write();
    m_file->Close();
    delete m_file;
    m_file = nullptr;
  }
}