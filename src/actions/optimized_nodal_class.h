#ifndef SIMPIMC_ACTIONS_OPTIMIZED_NODAL_CLASS_H_
#define SIMPIMC_ACTIONS_OPTIMIZED_NODAL_CLASS_H_

#include "nodal_class.h"
#include <einspline/multi_bspline.h>
#include <einspline/bspline.h>

class OptimizedNodal : public Nodal
{
private:

protected:
  // Variational parameter sets
  std::vector<std::vector<double>> param_sets;
  uint32_t param_set_i, model_i;

  // Rho matrix
  virtual double GetGij(const vec<double> &r, const uint32_t slice_diff);

  // 1/(4\lambda\tau)
  virtual double Geti4LambdaTau(const uint32_t slice_diff);

  // Splines
  field<UBspline_1d_d*> rho_node_r_splines;
  virtual void SetupSpline();

public:
  // Constructor
  OptimizedNodal(Path &path, RNG &rng, Input &in, IO &out)
    : Nodal(path,rng,in,out)
  {
    Init(in);
  }

  // Functions
  virtual void Init(Input &in);
  virtual void Write();
  uint32_t GetParamSet() { return param_set_i; };
  uint32_t GetNumParamSets() { return param_sets.size(); };
  void SetParamSet(uint32_t t_param_set_i) { param_set_i = t_param_set_i; };
  void SetRandomParamSet() { SetParamSet(rng.UnifRand(param_sets.size())-1); };

};

#endif // SIMPIMC_ACTIONS_OPTIMIZED_NODAL_CLASS_H_