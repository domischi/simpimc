#ifndef EventClass_H
#define EventClass_H

#include "Utils/config.h"       // Standard libraries
#include <sys/time.h>

class Event
{
private:

protected:

public:
  // Constructor
  Event() { timeSpent = 0.; }

  string name;
  RealType timeSpent;
  virtual void DoEvent() {};
  virtual void Write() {};
};

#endif
