#ifdef __CLING__
#pragma link C++ class Particle+;
#pragma link C++ class Parton+;
#pragma link C++ class Hadron+;
#pragma link C++ class Event+;

#pragma link C++ class std::vector<Particle*>+;
#pragma link C++ class std::vector<Parton*>+;
#pragma link C++ class std::vector<Hadron*>+;
#pragma link C++ class std::vector<Event*>+;

#pragma link C++ function GetMass;
#pragma link C++ function GetMassTable;
#endif