#include "core/Event.h"

unsigned int Event::sNextID = 1;

void Event::Reset(Option_t* opt) {
  // delete dynamically allocated Parton objects
  for (auto p : fPartons) {
    delete p;
  }
  // delete dynamically allocated Hadron objects
  for (auto h : fHadrons) {
    delete h;
  }
  // clear the containers
  fPartons.clear();
  fHadrons.clear();
  // call base class Clear
  TObject::Clear(opt);
}

ClassImp(Event)