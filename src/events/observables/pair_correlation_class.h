#ifndef SIMPIMC_OBSERVABLES_PAIR_CORRELATION_CLASS_H_
#define SIMPIMC_OBSERVABLES_PAIR_CORRELATION_CLASS_H_

#include "observable_class.h"

/// Record the pair correlation function between two species of particles
class PairCorrelation : public Observable
{
private:
  Histogram gr; ///< Histogram representing g(r)
  std::shared_ptr<Species> species_a; ///< First species
  std::shared_ptr<Species> species_b; ///< Second species
  bool image_summation; ///< indicated whether we need to sum also over the images
  int image_max; ///< if we need to sum over images this gives the maximal elongation along one direction
  /// Accumulate the observable
  virtual void Accumulate()
  {
    path.SetMode(NEW_MODE);
    // Homogeneous
    double cofactor = path.GetSign()*path.GetImportanceWeight();
    if (species_a == species_b) {
      for (uint32_t b_i=0; b_i<path.GetNBead(); ++b_i) {
        for (uint32_t p_i=0; p_i<species_a->GetNPart(); ++p_i) {
          //for (uint32_t p_j=p_i+1; p_j<species_b->GetNPart(); ++p_j) {
          for (uint32_t p_j=p_i; p_j<species_b->GetNPart(); ++p_j) {
            vec<double> dr(path.Dr(species_a->GetBead(p_i,b_i),species_b->GetBead(p_j,b_i)));
            double mag_dr=mag(dr);
            uint32_t i = gr.x.ReverseMap(mag(dr));
            // Direct Contribution
            if (i < gr.x.n_r&&p_i!=p_j)
              gr.y(i) = gr.y(i) + 1.*cofactor;
            // Image Contribution
            if (image_summation){
              vec<double> shift(path.GetND());
              for(int i0=-image_max;i0<image_max;++i0) 
              for(int i1=-image_max;i1<image_max;++i1) 
              for(int i2=-image_max;i2<image_max;++i2){
                shift[0]=i0*path.GetL();
                shift[1]=i1*path.GetL();
                shift[2]=i2*path.GetL();
                vec<double> dr_p=dr+shift;
                uint32_t i_p = gr.x.ReverseMap(mag(dr_p));
                if (i_p < gr.x.n_r)
                  gr.y(i_p) = gr.y(i_p) + 1.*cofactor;
              }
            }
          }
        }
      }
    }
    // Homologous
    else {
      for (uint32_t b_i=0; b_i<path.GetNBead(); ++b_i) {
        for (uint32_t p_i=0; p_i<species_a->GetNPart(); ++p_i) {
          for (uint32_t p_j=0; p_j<species_b->GetNPart(); ++p_j) {
            vec<double> dr(path.Dr(species_a->GetBead(p_i,b_i),species_b->GetBead(p_j,b_i)));
            uint32_t i = gr.x.ReverseMap(mag(dr));
            // Direct Contibution
            if (i < gr.x.n_r)
              gr.y(i) = gr.y(i) + 1.*cofactor;
            // Image Contribution
            if (image_summation){
              vec<double> shift(path.GetND());
              for(int i0=-image_max;i0<image_max;++i0) 
              for(int i1=-image_max;i1<image_max;++i1) 
              for(int i2=-image_max;i2<image_max;++i2){
                shift[0]=i0*path.GetL();
                shift[1]=i1*path.GetL();
                shift[2]=i2*path.GetL();
                vec<double> dr_p=dr+shift;
                uint32_t i_p = gr.x.ReverseMap(mag(dr_p));
                if (i_p < gr.x.n_r)
                  gr.y(i_p) = gr.y(i_p) + 1.*cofactor;
              }
            }
          }
        }
      }
    }
    n_measure += 1;
  }

  /// Reset the observable's counters
  virtual void Reset()
  {
    n_measure = 0;
    gr.y.zeros();
  }

public:
  /// Constructor calls Init and sets output data_type
  PairCorrelation(Path &path, Input &in, IO &out)
    : Observable(path, in, out, "histogram")
  {
    // Read in species info
    std::string species_a_name = in.GetAttribute<std::string>("species_a");
    std::string species_b_name = in.GetAttribute<std::string>("species_b");
    species_a = path.GetSpecies(species_a_name);
    species_b = path.GetSpecies(species_b_name);

    // Read in grid info
    double r_min = in.GetAttribute<double>("r_min",0.);
    double r_max = in.GetAttribute<double>("r_max",path.GetL()/2.);
    uint32_t n_r = in.GetAttribute<double>("n_r",1000);
    gr.x.CreateGrid(r_min,r_max,n_r);
    gr.y.zeros(n_r);

    // Do we need to sum over images
    image_summation=(r_max>path.GetL()/2.0);
    if(image_summation)
      image_max=ceil((r_max-path.GetL()/2.0)/path.GetL());
    else
      image_max=0;
    // Compute rs
    vec<double> rs(gr.x.n_r-1);
    for (uint32_t i=0; i<gr.x.n_r-1; i++) {
      double r1 = gr.x(i);
      double r2 = gr.x(i+1);
      if (path.GetND() == 3)
        rs(i) = 0.75 * (r2*r2*r2*r2-r1*r1*r1*r1)/(r2*r2*r2-r1*r1*r1);
      else if (path.GetND() == 2)
        rs(i) = (r2*r2*r2-r1*r1*r1)/(r2*r2-r1*r1); // FIXME: Not sure if 2D and 1D are correct here
      else if (path.GetND() == 1)
        rs(i) = 0.5*(r2-r1);
    }

    // Write things to file
    out.Write(prefix+"/species_a", species_a_name);
    out.Write(prefix+"/species_b", species_b_name);
    out.Write(prefix+"/r_min", r_min);
    out.Write(prefix+"/r_max", r_max);
    out.Write(prefix+"/n_r", n_r);
    out.Write(prefix+"/x", rs);

    Reset();
  }

  /// Write relevant information about an observable to the output
  virtual void Write()
  {
    if (n_measure > 0) {
      // Normalize histograms
      double vol = path.GetVol();
      double norm;
      if (species_a == species_b)
        //norm = 0.5*n_measure*species_a->GetNPart()*(species_b->GetNPart()-1)*path.GetNBead()/vol;
        norm = 0.5*n_measure*species_a->GetNPart()*path.GetNBead();
      else
        //norm = n_measure*species_a->GetNPart()*species_b->GetNPart()*path.GetNBead()/vol;
        norm = n_measure*species_a->GetNPart()*path.GetNBead();
      for (uint32_t i=0; i<gr.x.n_r; i++) {
        double r1 = gr.x(i);
        double r2 = (i<(gr.x.n_r-1)) ? gr.x(i+1):(2.*gr.x(i)-gr.x(i-1));
        double r = 0.5*(r1+r2);
        double bin_vol;
        if (path.GetND() == 3)
          bin_vol = 4.*M_PI/3. * (r2*r2*r2-r1*r1*r1);
        else if (path.GetND() == 2)
          bin_vol = M_PI * (r2*r2-r1*r1);
        else if (path.GetND() == 1)
          bin_vol = r2-r1;
        gr.y(i) = gr.y(i)/(bin_vol*norm);
        //gr.y(i) = gr.y(i)/(norm);
      }

      // Write to file
      if (first_time) {
        first_time = 0;
        out.CreateExtendableDataSet("/"+prefix, "y", gr.y);
      } else {
        out.AppendDataSet("/"+prefix, "y", gr.y);
      }

      Reset();
    }
  }

};

#endif // SIMPIMC_OBSERVABLES_PAIR_CORRELATION_CLASS_H_
