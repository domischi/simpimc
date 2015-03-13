#ifndef SIMPIMC_OBSERVABLES_OBSERVABLE_CLASS_H_
#define SIMPIMC_OBSERVABLES_OBSERVABLE_CLASS_H_

#include "../event_class.h"
#include "../actions/action_class.h"

class Observable : public Event
{
private:

protected:
  Path &path;
  IO &out;
  uint32_t skip;

  struct LinearGrid
  {
    vec<double> rs;
    double r_min, r_max, dr, d_ir;
    uint32_t n_r;

    void CreateGrid(double t_r_min, double t_r_max, uint32_t t_n_r)
    {
      r_min = t_r_min;
      r_max = t_r_max;
      n_r = t_n_r;
      dr = (r_max-r_min)/(n_r-1.);
      d_ir = 1./dr;
      rs.set_size(n_r);
      for (uint32_t i=0; i<n_r; ++i)
        rs(i) = r_min + i*dr;
    };
    inline double operator() (uint32_t i) { return rs(i); };
    uint32_t ReverseMap(const double r) {
      uint32_t i = (uint32_t) nearbyint((r-r_min)*d_ir-0.5);
      if (i < 0)
        return 0;
      else
        return i;
    };
  };

  struct Histogram
  {
    vec<double> y;
    LinearGrid x;
  };

  std::string prefix;

public:
  // Constructor
  Observable(Path &tmp_path, Input &in, IO &tmp_out)
    : Event(), path(tmp_path), out(tmp_out)
  {
    name = in.GetAttribute<std::string>("name");
    type = in.GetAttribute<std::string>("type");
    skip = in.GetAttribute<uint32_t>("skip",1);
    prefix = "Observables/"+name+"/";
    out.CreateGroup(prefix);
    out.Write(prefix+"/type",type);
    first_time = 1;
    n_measure = 0;
  }

  std::string type;
  bool first_time;
  uint32_t n_measure;

  // Functions
  inline void DoEvent() {
    struct timeval time;
    gettimeofday(&time, NULL); // Start Time
    double start = time.tv_sec + (time.tv_usec / 1000000.);
    Accumulate();
    gettimeofday(&time, NULL); //END-TIME
    double end = time.tv_sec + (time.tv_usec / 1000000.);
    time_spent += end - start;
  }
  virtual void Init(Input &in) {};
  virtual void Accumulate() {};
  virtual void Reset() {};
  virtual void Write() {};
};

#endif // SIMPIMC_OBSERVABLES_OBSERVABLE_CLASS_H_