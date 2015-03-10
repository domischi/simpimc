#include "SimulationClass.h"

void Simulation::SetupSimulation(string inFile)
{
  // Input
  in.load(inFile);

  // Build MPI Model
  procsPerGroup = in.getChild("Parallel").getAttribute<uint>("procsPerGroup",1);
  uint N = WorldComm.NumProcs();
  assert ((N % procsPerGroup) == 0);
  uint nGroups = N/procsPerGroup;
  uint myGroup = WorldComm.MyProc()/procsPerGroup;
  WorldComm.Split(myGroup, IntraComm);
  vec<int> ranks(nGroups);
  for (uint group=0; group<nGroups; group++)
    ranks(group) = group*procsPerGroup;
  WorldComm.Subset(ranks, InterComm);
  int nThreads = 1;
#if USE_OPENMP
  #pragma omp parallel
  {
    nThreads = omp_get_num_threads();
  }
#endif
  if (WorldComm.MyProc() == 0) {
    cout << "Running simpimc on " << inFile << endl;
    cout <<"# Processes: "<<N<< ", # Groups: "<<nGroups
         <<", # Processes/Group: "<<procsPerGroup
         <<", # Threads/Process: "<<nThreads<<endl;

  }

  // Output
  stringstream tmpSS;
  tmpSS << in.getChild("IO").getAttribute<string>("outputPrefix") << "." << myGroup << ".h5";
  string output = tmpSS.str();
  out.load(output);
  out.create();

  // Write input data
  out.CreateGroup("Input");
  out.Write("Input/fileName",inFile);
  string inString = in.getString();
  out.Write("Input/contents",inString);
}

void Simulation::Run()
{
#if USE_MPI
  int seed = in.getChild("RNG").getAttribute<int>("seed",(int)time(0)*(WorldComm.MyProc()+1));
#else
  int seed = in.getChild("RNG").getAttribute<int>("seed",(int)time(0));
#endif
  RNG rng(seed);
  out.CreateGroup("RNG");
  out.Write("RNG/seed",seed);

  // Algorithm
  Algorithm algorithm(WorldComm, InterComm, IntraComm);
  algorithm.Init(in, out, rng);
  algorithm.Run();

}

