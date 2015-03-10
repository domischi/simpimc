#ifndef BisectClass_H
#define BisectClass_H

#include "MoveClass.h"

class Bisect : public Move
{
private:
  string species;
  bool rollOver, adaptive;
  double targetRatio;
  int nImages;
  uint iSpecies;
  uint nLevel, nBisectBeads;
  uint bead0, bead1;
  double i4LambdaTauNBisectBeads, lambda;
  int refAccept, refAttempt;

  std::vector< std::shared_ptr<Bead> > affBeads;
protected:

public:
  // Constructor
  Bisect(Path &tmpPath, RNG &tmpRNG, std::vector< std::shared_ptr<Action> > &actionList, Input &in, IOClass &out)
    : Move(tmpPath, tmpRNG, actionList, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual bool Attempt();
  virtual void Accept();
  virtual void Reject();
  virtual void Reset();
  virtual void Write();

};

#endif
