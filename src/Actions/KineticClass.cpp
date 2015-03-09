#include "KineticClass.h"

void Kinetic::Init(Input &in)
{
  // Read in things
  nImages = in.getAttribute<int>("nImages");
  species = in.getAttribute<string>("species");
  speciesList.push_back(species);
  cout << "Setting up kinetic action for " << species << "..." << endl;
  path.GetSpeciesInfo(species,iSpecies);
  nPart = path.speciesList[iSpecies]->nPart;
  i4LambdaTau = 1./(4.*path.speciesList[iSpecies]->lambda*path.tau);

  // Write things to file
  out.Write("Actions/"+name+"/nImages", nImages);
  out.Write("Actions/"+name+"/species", species);

  // Setup spline
  SetupSpline();
}

// Create a spline for each possible sliceDiff
// TODO: Combine with free nodal action splines
void Kinetic::SetupSpline()
{
  // Setup grid
  Ugrid r_grid;
  if (path.PBC) {
    r_grid.start = -path.L/2.;
    r_grid.end = path.L/2.;
  } else {
    r_grid.start = -100.;
    r_grid.end = 100.;
    nImages = 0;
  }
  r_grid.num = 10000;
  double dr = (r_grid.end - r_grid.start)/(r_grid.num - 1);

  // Resize spline field
  int nSpline = path.nBead/2 + (path.nBead%2) + 1;
  rho_free_r_splines.set_size(nSpline);

  // Create splines
  for (int iSpline=0; iSpline<nSpline; ++iSpline) {
    vec<double> rho_free_r(r_grid.num), num_sum_r(r_grid.num);
    double t_i4LambdaTau = i4LambdaTau/(iSpline+1);

    // Make rho_free
    for (int i=0; i<r_grid.num; ++i) {
      double r = r_grid.start + i*dr;
      double r2 = r*r;
      double r2i4LambdaTau = r2*t_i4LambdaTau;
      rho_free_r(i) = 0.;
      if (iSpline == 0)
        num_sum_r(i) = 0.;
      for (int image=-nImages; image<=nImages; image++) {
        if (image != 0) {
          double t_r = r + image*path.L;
          double expPart = path.fexp(r2i4LambdaTau - t_r*t_r*t_i4LambdaTau);
          rho_free_r(i) += expPart;
          if (iSpline == 0 && r2 != 0.)
            num_sum_r(i) += (t_r*t_r/r2)*expPart;
        }
      }
      rho_free_r(i) = log1p(min(10.,rho_free_r(i)));
      num_sum_r(i) = log1p(min(10.,num_sum_r(i)));
    }
    BCtype_d xBC = {NATURAL, FLAT}; // fixme: Is this correct?
    UBspline_1d_d* rho_free_r_spline = create_UBspline_1d_d(r_grid, xBC, rho_free_r.memptr());
    rho_free_r_splines(iSpline) = rho_free_r_spline;
    if (iSpline == 0)
      num_sum_r_spline = create_UBspline_1d_d(r_grid, xBC, num_sum_r.memptr());
  }
}

double Kinetic::GetGaussSum(const double &r, const int sliceDiff)
{
  double gaussSum;
  eval_UBspline_1d_d(rho_free_r_splines(sliceDiff-1),r,&gaussSum);
  gaussSum = exp(0.9999*gaussSum);
  gaussSum *= exp(-(r*r*i4LambdaTau/sliceDiff));
  return gaussSum;
}

double Kinetic::GetNumSum(const double &r)
{
  double numSum;
  eval_UBspline_1d_d(num_sum_r_spline,r,&numSum);
  numSum = exp(0.9999*numSum);
  numSum *= -(r*r*i4LambdaTau/path.tau)*exp(-(r*r*i4LambdaTau));
  return numSum;
}

double Kinetic::DActionDBeta()
{
  double tot = nPart*path.nBead*path.nD/(2.*path.tau); // Constant term
  #pragma omp parallel for collapse(2) reduction(+:tot)
  for (int iP=0; iP<nPart; iP++) {
    for (int iB=0; iB<path.nBead; iB++) {
      vec<double> numSum(path.nD), gaussSum(path.nD);
      vec<double> dr(path.Dr(path(iSpecies,iP,iB),path.GetNextBead(path(iSpecies,iP,iB),1)));
      double gaussProd = 1.;
      for (int iD=0; iD<path.nD; iD++) {
        numSum(iD) = GetNumSum(dr(iD));
        gaussSum(iD) = GetGaussSum(dr(iD),1);
        gaussProd *= gaussSum(iD);
      }
      double scalarNumSum = 0.;
      for (int iD=0; iD<path.nD; iD++) {
        double numProd = 1.;
        for (int jD=0; jD<path.nD; jD++) {
          if (iD != jD)
            numProd *= gaussSum(jD);
          else
            numProd *= numSum(jD);
        }
        scalarNumSum += numProd;
      }
      tot += scalarNumSum/gaussProd;
    }
  }

  return tot;
}

double Kinetic::GetAction(const int b0, const int b1, const vector< pair<int,int> > &particles, const int level)
{
  int skip = 1<<level;
  double i4LambdaLevelTau = i4LambdaTau/skip;
  double tot = 0.;
  for (auto& p: particles) {
    int iS = p.first;
    int iP = p.second;
    if (iS == iSpecies) {
      std::shared_ptr<Bead> beadA(path(iSpecies,iP,b0));
      std::shared_ptr<Bead> beadF(path.GetNextBead(beadA,b1-b0));
      while(beadA != beadF) {
        std::shared_ptr<Bead> beadB(path.GetNextBead(beadA,skip));
        vec<double> dr(path.Dr(beadA,beadB));
        double gaussProd = 1;
        for (int iD=0; iD<path.nD; iD++)
          gaussProd *= GetGaussSum(dr(iD),skip);
        tot -= log(gaussProd);
        beadA = beadB;
      }
    }
  }

  return tot;
}

vec<double> Kinetic::GetActionGradient(const int b0, const int b1, const vector< pair<int,int> > &particles, const int level)
{
  int skip = 1<<level;
  double i4LambdaLevelTau = i4LambdaTau/skip;
  vec<double> tot;
  tot.zeros(path.nD);
  std::shared_ptr<Bead> beadA, beadB, beadC, beadF;
  for (auto& p: particles) {
    int iS = p.first;
    int iP = p.second;
    if (iS == iSpecies) {
      double gaussProd, gaussSum, dist;
      beadA = path(iSpecies,iP,b0);
      beadF = path.GetNextBead(beadA,b1-b0);
      while(beadA != beadF) {
        beadB = path.GetPrevBead(beadA,skip);
        vec<double> dr(path.Dr(beadB,beadA));
        tot -= dr;
        beadC = path.GetNextBead(beadA,skip);
        dr = path.Dr(beadA,beadC);
        tot += dr;
        beadA = beadC;
      }
    }
  }

  return 2.*i4LambdaLevelTau*tot;
}

double Kinetic::GetActionLaplacian(const int b0, const int b1, const vector< pair<int,int> > &particles, const int level)
{
  int skip = 1<<level;
  double i4LambdaLevelTau = i4LambdaTau/skip;
  double tot = 0.;
  vec<double> dr(path.nD);
  std::shared_ptr<Bead> beadA, beadF;
  for (auto& p: particles) {
    int iS = p.first;
    int iP = p.second;
    if (iS == iSpecies) {
      double gaussProd, gaussSum, dist;
      beadA = path(iSpecies,iP,b0);
      beadF = path.GetNextBead(beadA,b1-b0);
      while(beadA != beadF) {
        tot += path.nD*4.*i4LambdaLevelTau;
        beadA = path.GetNextBead(beadA,skip);
      }
    }
  }

  return tot;
}

void Kinetic::Write()
{

}
