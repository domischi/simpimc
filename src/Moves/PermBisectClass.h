#ifndef PermBisectClass_H
#define PermBisectClass_H

#include "MoveClass.h"

class PermBisect : public Move
{
private:
  string species;
  int iSpecies;
  int offset;
  int nImages;
  unsigned int nLevel, nBisectBeads, nPart, nPermPart, nPermType;
  unsigned int bead0, bead1;
  RealType lambda, i4LambdaTauNBisectBeads, epsilon, logEpsilon;

  struct Cycle
  {
    RealType weight, contribution;
    int index, type;
    Ivector perm, iPerm, part;
  };
  vector<Cycle*> cycles;
  field<Cycle> all_cycles;
  Tmatrix t;

  int permType;
  Ivector permAttempt, permAccept;

  RealType constructPermTable();
  void updatePermTable();
  void BuildCycles();
  int selectCycle(RealType permTot);
  void permuteBeads(field<Bead*>& b0, field<Bead*>& b1, Cycle* c);
  void assignParticleLabels();
  void Write();

  std::vector<Bead*> affBeads;
protected:

public:
  // Constructor
  PermBisect(Path &tmpPath, RNG &tmpRNG, vector<Action*> &actionList, Input &in, IOClass &out)
    : Move(tmpPath, tmpRNG, actionList, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual int Attempt();
  virtual void Accept();
  virtual void Reject();
};


#endif
