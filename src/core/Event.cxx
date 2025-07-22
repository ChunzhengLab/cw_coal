#include "core/Event.h"

#include <algorithm>
#include <numeric>
#include <random>

unsigned int Event::sNextID = 1;

void Event::Reset(Option_t* opt) {
  // delete dynamically allocated Parton objects
  for (auto p : fPartons) { delete p; }
  // delete dynamically allocated Hadron objects
  for (auto h : fHadrons) { delete h; }
  // clear the containers
  fPartons.clear();
  fHadrons.clear();
  TObject::Clear(opt);
}

void Event::ShufflePartons(ShuffleLevel level) {
  double frac = 0.0;
  switch (level) {
    case ShuffleLevel::kLevel1:
      frac = 0.25;
      break;
    case ShuffleLevel::kLevel2:
      frac = 0.50;
      break;
    case ShuffleLevel::kLevel3:
      frac = 0.75;
      break;
    case ShuffleLevel::kLevel4:
      frac = 1.00;
      break;
  }
  ShufflePartons(frac);
}

void Event::ShufflePartons(double fraction) {
  size_t n = fPartons.size();
  if (n < 2) return;
  if (fraction <= 0.0) return;
  if (fraction > 1.0) fraction = 1.0;

  size_t num = static_cast<size_t>(fraction * n);
  if (num < 2) return;

  std::vector<size_t> idx(n);
  std::iota(idx.begin(), idx.end(), 0);

  static thread_local std::mt19937 rng{std::random_device{}()};
  std::shuffle(idx.begin(), idx.begin() + num, rng);

  std::vector<std::array<double, 3>> pos(num);
  for (size_t i = 0; i < num; ++i) pos[i] = fPartons[idx[i]]->GetPosition();

  std::shuffle(pos.begin(), pos.end(), rng);

  for (size_t i = 0; i < num; ++i) {
    auto p = fPartons[idx[i]];
    auto [x, y, z] = pos[i];
    p->SetPosition(x, y, z);
  }
}

ClassImp(Event)