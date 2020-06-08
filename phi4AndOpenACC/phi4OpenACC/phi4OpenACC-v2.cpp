//=================================================================================== 
// This is 2D phi^4 code on a torus designed to be easy to converted for openACC.
//===================================================================================

/*
v0 additions
This version is a near straight up copy of the serial version of this routine.
It is used purely to demonstrate that PGI++ can compile serial code, and to give 
some timings for:

1. Lattice Percolation
2. Swendsen Wang cluster decomposition
3. HMC trajectory
4. Phi mean and moments

v1 additions
This version is an evolution of v0, whose additions remain. In this file we shall
parallelise the HMC using some important pragmas:

1. #pragma acc data present array[...]
   This pragma tells the compiler that the data in `array` is already present on 
   the device. As such, PGI++ does not have to generate implicit copies
   of `array` each time a GPU accelerator loops is compiled

2. #pragma acc parallel loop collapse(N)
   Here we tell the compiler that there are two loops in a nested structure,
   and that each iteration of the loop is independent. We say that the loop
   may be 'flattened'.

3. #pragma acc parallel loop tile(N,M)
   Here we tell the compiler that there is spatial locality in the data that
   it can exploit. for example, one thread at (x,y) must fetch data from 
   site x+1 to perform a differential operator calculation in the +ve x-dim. 
   A different thread at (x+1,y-1) needs the exact same data to perform the
   differential operator in the +ve y-dim. The tile(N,M) clause instructs
   the compiler to store NxM lots of data in fast memory cache. You must 
   experiment with N and M.

4. #pragma acc update self (array[])
   #pragma acc update device (array[])
   These pragmas assume that `array` is defined on both the device and the host
   'update self' will copy the data on the device array to the host array.
   'update device' will copy the data on the host array to the device array.

5. Data causes. The 'data' clause in a pragma indicates that we are declaring memory 
   on the GPU. Three types of clause are listed:
    
   a) copy(A[0:N])      The 'copy' clause tells the compiler that it should copy 
                        the HOST data in the A array to the GPU at the start of 
                        the caluclation, and it is also to copy the data in the 
                        GPU back to the host at the end of the calculation.
    
   b) copyin(rhs[0:N])  The 'copyin' clause tells the compiler to copy the data 
                        from the HOST to the GPU at the start, but it should not 
                        copy the data in the GPU back to the host at the end. 
                        Can you guess what the opposite clause would be? To not 
                        copy the HOST data to the GPU at the start, but to copy 
                        from the GPU to the HOST at the end?
    
   c) create(Anew[0:N]) The 'create' clause simply creates memory on the device. 
                        No data transfer is performed for this array

v2 additions
This version is an evolution of v1, whose additions remain. In this file we shall
parallelise the cluster routines.

1. The lattice percolation is easily parallelised. We pass random numbers to the
   routine due to the unreliability of GPU RNGs

2. The Swendsen-Wang cluster identification routine requires some chicanery. 
   We take advantage of the fact that not every thread launches simultaneously
   and so some parts of the array are updated before others. This latency in
   updating the labels is parlayed into a faster method for propagating the 
   labels.

*/

#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <chrono>
using namespace std;
#include <cmath>
#include <complex>

//v1: If using the GPU, we must include the GPU specific math library
#ifdef GPU
#include <accelmath.h> //PGI math lib 
#endif
 
#define TWO_PI 6.283185307179586
#define L 128

// Parameter structure
//--------------------------------------------------------------------------
typedef struct{
  int Latsize;    //size of lattice
  double tau;     //Leapfrog integration trajectory length
  double dt;      //Leapfrog integration step size
  int nStep;      //number of HMC integration steps
  int nMeas;      //Measure every nMean iteration
  int nIter;      //Number of iterations
  int nHMC;       //Perform an HMC step every nHMC iter
  int nTherm;     //Thermalise to this trajectory
  double lambda;  //phi**4 coupling
  double musqr;   //phi**2 coupling
  bool verbose;   //Print stats to stdout
  bool debug;     //Print more stats to stdout
  int seed;       //RNG seed
  int relaxIter;  //Number of relaxation sweeps
  //v2: New parameters
  int relaxInit;  //Number of initial relaxations
  int relaxLoops; //Number of relaxation loops in GPU

} param;

//Utilities
void printLattice(const double phi[L][L]);
void hotStart(double phi[L][L], param p);
void coldStart(double phi[L][L], param p);
//v1: upgraded function to compute momenta on the device
void gaussReal_F(double field[L][L], double rands[L][L][2]);

//Measurements
double calcH(double mom[L][L], double phi[L][L],param p);
double measMag(const double phi[L][L]);
double measAbsMag(const double phi[L][L]);

//Hybrid Monte-Carlo
int hmc(double phi[L][L], param p, int iter);
void force_phi(double fphi[L][L], double phi[L][L], param p);
void update_mom(double mom[L][L], double fphi[L][L], param p, double dt);
void update_phi(double phi[L][L], double mom[L][L], param p, double dt);
void trajectory(double mom[L][L], double phi[L][L], param p);

//Cluster decomposition
void latticePercolate(bool bond[L][L][4], int label[L][L], const double phi[L][L]);
//v2: swendsenWang with no stopping condition
void swendsenWang(int label[L][L], const bool bond[L][L][4], int loops);
//v2: swendsenWang with a stopping condition
bool swendsenWangStopTest(int label[L][L], const bool bond[L][L][4]);
void flipSpins(double phi[L][L], const int label[L][L]);

// Global variables for HMC diagnostics
//--------------------------------------------------------------------------
double ave_expmdH = 0.0;
double ave_dH = 0.0;
//v2: New Statistic
int ave_relax = 0;

// v0: Global variables for timing and new timing function
//--------------------------------------------------------------------------
double start = 0.0;
double LPtime = 0.0;
double SWtime = 0.0;
double HMCtime = 0.0;
double PhiMtime = 0.0;
#include <sys/time.h>
inline double get_time() {
  
  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv, &tz);
  return 1.0*tv.tv_sec+1.0E-6*tv.tv_usec;

}

// Inlined Device Utilities
//--------------------------------------------------------------------------
// Copy lattice field (Device)
template<typename T> inline void copyField(T phi2[L][L], T phi1[L][L]) {
#pragma acc data present(phi2[0:L][0:L]) present(phi1[0:L][0:L])
  {
#pragma acc parallel loop collapse(2)    
    for(int x=0; x<L;x++)
      for(int y=0; y<L;y++)
	phi2[x][y] = phi1[x][y];
  }
}

//--------------------------------------------------------------------------

int main(int argc, char **argv) {
  
  double prog_start = get_time();

  param p;
  p.Latsize = L;
  p.lambda = 0.125; 
  p.musqr = -0.359;
  p.nStep = 15;
  p.tau = 1.0;
  p.dt = p.tau/p.nStep;
  p.nTherm = 100;
  p.nIter = 100;
  p.nMeas = 10;
  p.nHMC = 5;
  p.verbose = true;
  p.debug = false;
  p.seed = 1234;
  p.relaxIter = 10000;
  //v2: New parameters
  p.relaxInit = 1000;
  p.relaxLoops = 50;
  srand(p.seed);

  cout << "lattice size = " << L << endl;
  cout << "lambda = " << p.lambda << endl;
  cout << "musqr = " << p.musqr << endl;
  cout << "trajectory steps = "<< p.nStep << endl;
  cout << "trajectory length = " << p.tau << endl;
  cout << "trajectory step length = " << p.dt << endl;
  cout << "thermalisation iters = " << p.nTherm << endl;
  cout << "total iterations = " << p.nIter << endl;
  cout << "HMC performed every " << p.nHMC << " iters" << endl;
  cout << "relaxation sweep max = " << p.relaxIter << endl;
  //v2: New Parameters
  cout << "relaxation init = " << p.relaxInit << endl;
  cout << "relaxation loops = " << p.relaxLoops << endl;
  cout << "seed = " << p.seed << endl;
  
  double phi[L][L];
  // coldStart(phi,p);
  hotStart(phi, p);
  
  int accepted = 0;
  int measurement = 0;
  int nHMC = 0;
  double getMag = 0.0;
  double AbsPhi = 0.0;
  double Phi = 0.0;  
  double Phi2 = 0.0;
  double Phi4 = 0.0;
  bool stop = false;
  bool bond[L][L][4];
  int label[L][L];

  double vol = L*L;
  double vol2= vol*vol;
  double vol4= vol*vol*vol*vol;
  
  // v1: OpenACC Init 
#pragma acc init
  
  //Enter OpenACC parallel region
  //======================================================================
  // 'copyout' will create arrays on the device and copy host data to
  // the device. 'create' will create space on the device only. See notes above.

  //v2: Now that we are going to compute the clusters on the device, we should
  //    create the device arrays here, rather than repeatedly creating them.
#pragma acc data copyout(phi[0:L][0:L]) create(bond[0:L][0:L][0:4], label[0:L][0:L])
  {
        
    //--------------------------------------------------------------------------
    // Thermalise the field. This evolves the phi field
    // from a completely random (hot start)
    // or a completely ordered (cold start)
    // to a state of thermal equlibrium.
    for(int iter = 0; iter < p.nTherm; iter++) {
      
      if(iter%p.nHMC == 0) hmc(phi, p, iter);
      //phi has been updated on the device and the host

      latticePercolate(bond, label, phi);

      stop = false;
      swendsenWang(label, bond, p.relaxInit);
      for(int relax = 0; relax < p.relaxIter && !stop; relax += p.relaxLoops) {

	swendsenWang(label, bond, p.relaxLoops-1);
	stop = swendsenWangStopTest(label, bond);

	if(stop) {
	  //Lattice clusters identified. Flip 'em!
	  flipSpins(phi, label);
	  if(p.debug) {
	    cout << "Iter " << iter << ": relaxation sweeps: " << relax << endl;
	    for(int x=0; x<L; x++)
	      for(int y=0; y<L; y++) {
		cout << label[x][y] << " ";
		if(y==L-1) cout<<endl;
	      }
	  }	  
	}
	
	if(relax > p.relaxIter) {
	  cout << "Error in SW: clusters not identified in " << p.relaxIter << " relaxation sweeps. Exiting..." << endl;
	  exit(0);
	}
      }
    }
      
    // Thermalisation complete
    //--------------------------------------------------------------------------
    
    
    // Begin sampling the phase space of the thermalised lattice
    //--------------------------------------------------------------------------
    //reset acceptance
    accepted = 0;
    for(int iter = p.nTherm; iter < p.nIter + p.nTherm; iter++) {
      
      if(iter%p.nHMC ==  0) {
	accepted += hmc(phi, p, iter); 
	nHMC++;
      }
      
      latticePercolate(bond, label, phi);
      
      stop = false;
      swendsenWang(label, bond, p.relaxInit);
      for(int relax = 0; relax < p.relaxIter && !stop; relax += p.relaxLoops) {

	swendsenWang(label, bond, p.relaxLoops-1);
	stop = swendsenWangStopTest(label, bond);

	if(stop) {
	  //Lattice clusters identified. Flip 'em!
	  ave_relax += relax + p.relaxInit;
	  flipSpins(phi, label);
	  if(p.debug) {
	    cout << "Iter " << iter << ": relaxation sweeps: " << relax << endl;
	    for(int x=0; x<L; x++)
	      for(int y=0; y<L; y++) {
		cout << label[x][y] << " ";
		if(y==L-1) cout<<endl;
	      }
	  }	  
	}
	
	if(relax > p.relaxIter) {
	  cout << "Error in SW: clusters not identified in " << p.relaxIter << " relaxation sweeps. Exiting..." << endl;
	  exit(0);
	}
      }
      
      // Take measurements
      //--------------------------------------------------------------------------
      if((iter+1)%p.nMeas == 0) {
	
	measurement++;

	//Ensure host phi array is updated
#pragma acc update self (phi[0:L][0:L])
	
	//Get phi mean and moments
	//------------------------------------------------------------------------
	start = get_time();
	
	AbsPhi += measAbsMag(phi);
	
	getMag = measMag(phi);
	Phi += getMag;
	
	getMag *= getMag;
	Phi2 += getMag;
	
	getMag *= getMag;
	Phi4 += getMag;      
	
	double avPhi = Phi/measurement;
	double avPhi2 = Phi2/measurement;
	double avPhi4 = Phi4/measurement;
	PhiMtime += get_time() - start;
	//--------------------------------------------------------------------------
	
	// Diagnostics and observables
	//---------------------------------------------------------------------
	// 1. Try to keep the acceptance rate between 0.75 - 0.85. Do this by 
	//    varying the number of steps of the HMC integrator, keeping tau=1.0.
	// 2. <exp(-dH)> should be ~1.0.
	// 3. <dH> should be ~0.0, and ever so slightly positive. Why?
	// 4. <phi>, <phi**2>, <phi**4> and the Binder cumulant depend on
	//    musqr, lambda, and the problem size.
	if(p.verbose) {
	  cout << "measurement " << measurement << endl;
	  cout << "HMC acceptance rate = " << (double)accepted/nHMC << endl;
	  cout << "HMC <exp(-dH)>      = " << ave_expmdH/nHMC << endl;
	  cout << "HMC <dH>            = " << ave_dH/nHMC << endl;
	  cout << "CLUSTER <relax>     = " << (double)ave_relax/(iter+1-p.nTherm) << endl;
	  cout << "MEAS |phi|          = " << setprecision(12) << AbsPhi/(vol*measurement) << endl;
	  cout << "MEAS <phi>          = " << setprecision(12) << avPhi/vol << endl;
	  cout << "MEAS <phi**2>       = " << setprecision(12) << avPhi2/vol2 << endl;
	  cout << "MEAS <phi**4>       = " << setprecision(12) << avPhi4/vol4 << endl;
	  cout << "Binder Cumulant     = " << setprecision(12) << 1.0 - Phi4/(3.0*Phi2*Phi2/measurement) << endl;
	  
	}
	cout << "Timings: LP=" << LPtime << " SW=" << SWtime << " HMC=" << HMCtime << " PhiM=" << PhiMtime << " Total=" << get_time() - prog_start << endl; 
	//--------------------------------------------------------------------------
	
	// Dump data to file
	//--------------------------------------------------------------------------
	string filename;
	char fname[256];
	FILE *fp;
	filename = "phi4Data_L" + to_string(p.Latsize) + "_lambda" + to_string(p.lambda) +  "_musqr" + to_string(p.musqr) + "_tau"+ to_string(p.tau) +"_nHMCStep" + to_string(p.nStep) + "_nHMC" + to_string(p.nHMC) + "_seed" + to_string(p.seed) + ".dat";
	
	sprintf(fname, "%s", filename.c_str());	
	fp = fopen(fname, "a");	
	fprintf(fp, "%d %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
		iter+1,
		(double)accepted/nHMC,
		ave_expmdH/nHMC,
		ave_dH/nHMC,
		(double)ave_relax/(iter+1-p.nTherm),
		AbsPhi/(vol*measurement),
		avPhi/vol,
		avPhi2/vol,
		avPhi4/vol4,
		1.0 - Phi4/(3.0*Phi2*Phi2/measurement));
	fclose(fp);
	//--------------------------------------------------------------------------
      }
    }
    // End of simulation
    //--------------------------------------------------------------------------
  }
  // End of OpenACC region
  //--------------------------------------------------------------------------
  return 0;
}


//Cluster decomposition
//----------------------------------------------------------------------------
// Identify the bond structure of the Swendsen-Wang decomposition.
// v2: RNGs on the GPU do not have a solid reputation. This is an active line of 
// research in computational science. Instead, we take the performance hit of 
// generating the random numbers on the host, and trasferring them to the device.

void latticePercolate(bool bond[L][L][4], int label[L][L], const double phi[L][L]) {
  
  start = get_time();
  
  //Populate host rand arrays
  double rands1[L][L];
  double rands2[L][L];
  double rands3[L][L];
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++){
      rands1[x][y] = drand48();
      rands2[x][y] = drand48();
      rands3[x][y] = drand48();
    }
  
#pragma acc update device (bond[0:L][0:L][0:4], label[0:L][0:L])

  //v2: Notice that rands arrays are 'copied in' only. We don't care about the
  //    data they contain once we're done with the function.
#pragma acc data present(bond[0:L][0:L][0:4], label[0:L][0:L], phi[0:L][0:L]) copyin(rands1[0:L][0:L], rands2[0:L][0:L], rands3[0:L][0:L])
  {
    
    //v2: Notice that the bond array is not initialised, rather we allow the 
    //    compute threads to decide for is using an if/else. 
    //    CLASS QUESTION: Is this desireable?
#pragma acc parallel loop collapse(2)
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++){
	
	double probability = 1.0 - exp(-2.0 * phi[x][y]*phi[(x+1)%L][y]);
	if (rands1[x][y] < probability) {
	  bond[x][y][0] = true;
	  bond[(x+1)%L][y][2] = true;
	}
	else {
	  bond[x][y][0] = false;
	  bond[(x+1)%L][y][2] = false;
	} 
	
	probability = 1.0 - exp(-2.0 * phi[x][y]*phi[x][(y+1)%L]);
	if(rands2[x][y] < probability) {
	  bond[x][y][1] = true;
	  bond[x][(y+1)%L][3] = true;
	}
	else {
	  bond[x][y][1] = false;
	  bond[x][(y+1)%L][3] = false;
	}
	
	// Random spin on labels = p/m (1, 2, ... L*L):
	if(rands3[x][y] < 0.5) label[x][y] = -(y + x*L + 1);
	else label[x][y] = (y + x*L + 1);
      }
  }
  
  //Ensure to update the bond and label array for host side testing.
#pragma acc update self (bond[0:L][0:L][0:4], label[0:L][0:L])

  LPtime += get_time() - start;
  return;
}
  
  
//Use the bonds computed in latticePercolate() to identify the culsters.
//Find min label connection to local 4 point stencil. Each lattice site
//starts with a label identical to its lexicographic 2D coord. 
//We test to see if its neighbours are connected via a bond, and have a
//lower label. If so, this lower label is propagated to the site. After a finite number
//of applications of this function, each site in each cluster (as defined by the
//bond structure) will have a label unique to the cluster, equal to the
//smallest lexicographic coordinate in that cluster.
void swendsenWang(int label[L][L], const bool bond[L][L][4], int loops) {
    
  start = get_time();
  
#pragma acc update device(label[0:L][0:L], bond[0:L][0:L][0:4])

#pragma acc data present(bond[0:L][0:L][0:4], label[0:L][0:L])
  {
    //Use a quirk of the architecture to accelerate the relaxtion procedure.
    //A single relaxtion sweep is defined as a loop over x and y. By introducing
    //a new outer loop, it would appear to perform redundant steps. However, we
    //change the labels in global memory as soon as the thread is complete and the
    //streaming multiprocessors launch in an arbitrary sequence. As a result, the
    //labels are able to propagate.
    //
    //CLASS QUESTION: Why is this algorithm faster on the GPU?
    
    //v2: try using a tile(N,M). CLASS QUESTION: Why did that happen?
#pragma acc parallel loop collapse(3)
    for(int a=0; a<loops; a++) {
      for(int x=0; x<L; x++)
	for(int y=0; y<L; y++) {
	  // Find min of connection to local 4 point stencil.
	  int minLabel = label[x][y];
	  
	  if( bond[x][y][0] && (abs(minLabel) > abs(label[(x+1)%L][y])) ) {
	    minLabel = label[(x+1)%L][y];
	  }
	  
	  if( bond[x][y][1] && (abs(minLabel) > abs(label[x][(y+1)%L])) ) {
	    minLabel = label[x][(y+1)%L];
	  }
	  
	  if( bond[x][y][2] && (abs(minLabel) > abs(label[(x-1+L)%L][y])) ) {
	    minLabel = label[(x-1+L)%L][y];
	  }
	  
	  if( bond[x][y][3] && (abs(minLabel) > abs(label[x][(y-1+L)%L])) ) {
	    minLabel = label[x][(y-1+L)%L];
	  }
	  
	  label[x][y] = minLabel;
	}
    }//end loops
  }
  
#pragma acc update self(label[0:L][0:L])

  SWtime += get_time() - start;
}

bool swendsenWangStopTest(int label[L][L], const bool bond[L][L][4]) {

  start = get_time();
  
  bool stop = true;
  int newlabel[L][L];
  int minLabel;
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++) {

      minLabel = label[x][y];  
      
      if( bond[x][y][0] && (abs(minLabel) > abs(label[(x+1)%L][y])) ) {
	minLabel = label[(x+1)%L][y];
	stop = false;
      }
      
      if( bond[x][y][1] && (abs(minLabel) > abs(label[x][(y+1)%L])) ) {
	minLabel = label[x][(y+1)%L];
	stop = false;
      }
      
      if( bond[x][y][2] && (abs(minLabel) > abs(label[(x-1 + L)%L][y])) ) {
	minLabel = label[(x-1 + L)%L][y];
	stop = false;
      }
      
      if( bond[x][y][3] && (abs(minLabel) > abs(label[x][(y-1 + L)%L])) ) {
	minLabel = label[x][(y-1 + L)%L];
	stop = false;
      }      
      newlabel[x][y] =  minLabel;
    }
  
  //After sweeping through the lattice, update all the cluster labels 
  for(int x = 0; x< L; x++)
    for(int y = 0; y< L; y++)
      label[x][y] = newlabel[x][y];

  SWtime += get_time() - start;

  return stop;
}

//The clusters will all have the same spin as the spin of the lowest
//coordinate phi value. We randomly placed a minus sign on the label
//which will now serve as a 50/50 test to flip the entire cluster

void flipSpins(double phi[L][L], const int label[L][L]) {

  start = get_time();
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++) 
      if(label[x][y] < 0) phi[x][y] = -phi[x][y];

#pragma acc update device (phi[0:L][0:L])
  
  SWtime += get_time() - start;
  
}

//Hybrid Monte-Carlo
//---------------------------------------------------------------------------------
//Measure the action of the phi field. Perform an HMC trajectory on the field, and
//measure the new action. We accept/reject the new field according to the
//Metropolis test.
//
// NOTES on Metropolis:
// 1. If the candidate action is less than the old action, we always accept the
//    candidate lattice. If this happens, it means that
//    (due to Hamilton's theory of exremised action) the lattice has evolved to a
//    state which is closer to the stationary point classicial action.
// 2. If the candidate action is greater than the old action, we have evolved the lattice
//    away from the stationary point. This is allowed in a quantum field theory, but
//    we may only add it to the Markov chain with a gaussian probability, facilitated
//    by the Metropolis step.
// 3. The astute reader will note that the discrete integration we perform
//    in the trajectory() step will come with errors. The Metropolis step
//    'accounts' for these integration errors by accepting or rejecting
//    candidate configurations. (CLASS THOUGHT EXPERIMENT)

int hmc(double phi[L][L], param p, int iter) {

  start = get_time();  
  
  int accepted = 0;
  double rands[L][L][2];
  double mom[L][L];

  //v1: Make sure device phi is up to date
#pragma acc update device (phi[0:L][0:L])
  
#pragma acc data present(phi[0:L][0:L]) create(rands[0:L][0:L][0:2], mom[0:L][0:L])
  { 
    
    //populate rands array
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++) {
	rands[x][y][0] = drand48();
	rands[x][y][1] = drand48();
      }
    
    //Copy random numbers to the device
#pragma acc update device(rands[0:L][0:L][0:2])
    
    //Generate gaussian distributed momemta on the device
    gaussReal_F(mom, rands);
    
    //Perform trajectory
    double Hold = calcH(mom, phi, p);
    trajectory(mom, phi, p); // MD trajectory using Leapfrog
    double H = calcH(mom, phi, p);
    
    //record HMC diagnostics
    if(iter >= p.nTherm) {
      ave_expmdH += exp(-(H-Hold));
      ave_dH += H-Hold;
    }
  
    // Metropolis accept/reject step
    // Always accepts trajectories during first half of warm up.
    if (drand48() > exp(-(H-Hold)) && iter > p.nTherm/2-1) {
      //Keep old field. Host phi field does not need to be updated.
      //v2: Now that we are going to use phi on a device cluster function, 
      //    we should ensure that the device side array is kept up-to-date
#pragma acc update device (phi[0:L][0:L])
      accepted = 0;
    }
    else {
      //Accept candidate field
      //We must update the host phi array with the new field 
#pragma acc update self (phi[0:L][0:L])
      accepted = 1;
    }      
  }
  HMCtime += get_time() - start;
  return accepted;
  
}

//Performs the HMC trajectory via Leapfrog integration
void trajectory(double mom[L][L], double phi[L][L], param p) {

  const double dt = p.dt;
  const int nStep = p.nStep;
  
  double fphi[L][L];
#pragma acc data create(fphi[0:L][0:L])
  {
    //Initial half step:
    //P_{1/2} = P_0 - dtau/2 * fphi
    force_phi(fphi, phi, p);
    update_mom(mom, fphi, p, 0.5*dt);
    
    for(int k=1; k<nStep; k++) {
      
      //phi_{k} = phi_{k-1} + P_{k-1/2} * dt
      update_phi(phi, mom, p, dt);
      
      //P_{k+1/2} = P_{k-1/2} - fphi * dt 
      force_phi(fphi, phi, p);
      update_mom(mom, fphi,  p, dt);      
    }
    
    //Final half step.
    //phi_{n} = phi_{n-1} + P_{n-1/2} * dt
    update_phi(phi, mom, p, dt);
    force_phi(fphi, phi, p);
    update_mom(mom, fphi, p, 0.5*dt);
  }
  
  return;
}

//Compute the HMC force from the field
void force_phi(double fphi[L][L], double phi[L][L], param p) {
  
#pragma acc data present(fphi[0:L][0:L], phi[0:L][0:L])
  {
#pragma acc parallel loop tile(32,2) 
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++) {
	fphi[x][y] = 0.0;
	fphi[x][y] -= phi[(x+1)%L][y] - 2.0*phi[x][y] + phi[(x-1+L)%L][y];
	fphi[x][y] -= phi[x][(y+1)%L] - 2.0*phi[x][y] + phi[x][(y-1+L)%L];
	fphi[x][y] += (2.0*p.musqr + 4.0*p.lambda*phi[x][y]*phi[x][y]) * phi[x][y];
      }
  }
}
  
//Update the momenta from the back reaction of the field
void update_mom(double mom[L][L], double fphi[L][L], param p, double dt) {

#pragma acc data present(mom[0:L][0:L], fphi[0:L][0:L])
  {
#pragma acc parallel loop collapse(2)  
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++) 
	mom[x][y] -= fphi[x][y] * dt;
    
  }
}

//Update the field from the back reaction of the momenta
void update_phi(double phi[L][L], double mom[L][L], param p, double dt) {

#pragma acc data present(phi[0:L][0:L], mom[0:L][0:L])
  {
#pragma acc parallel loop collapse(2)  
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++) 
	phi[x][y] += mom[x][y] * dt;
    
  }
}
//---------------------------------------------------------------------------------



// Measurements
//------------------------------------------------------------------------------------
//Compute the action of the field and momenta
double calcH(double mom[L][L], double phi[L][L], param p) {

  const double musqr = p.musqr;
  const double lambda = p.lambda;
  
  double Hout = 0.0;
  //We indicate that the data is present on the device
#pragma acc data present(mom[0:L][0:L], phi[0:L][0:L])
  {
    //If we declare a variable inside a parallel region,
    //the compiler will ensure that each thread has
    //a private copy.
    double H = 0.0;
#pragma acc parallel loop tile(32,4) reduction(+:H)
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++) {
	
	H += (  0.5*(phi[x][y] - phi[(x+1)%L][y]) * (phi[x][y] - phi[(x+1)%L][y])
	      + 0.5*(phi[x][y] - phi[x][(y+1)%L]) * (phi[x][y] - phi[x][(y+1)%L])
	      + (musqr + lambda*phi[x][y]*phi[x][y])*phi[x][y]*phi[x][y] 
	      + 0.5 * mom[x][y] * mom[x][y]);
	
      }
    Hout = H;
  }
  return Hout;
}

//Compute the 'magnetisation' of the phi field
double measMag(const double phi[L][L]) {
  double mag = 0.0;
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      mag += phi[x][y];  
  return mag;
}

//Compute the absolute value of the phi field
double measAbsMag(const double phi[L][L]) {
  double mag = 0.0;
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      mag += fabs(phi[x][y]);  
  return mag;
}
//------------------------------------------------------------------------------------




// Utilities
//------------------------------------------------------------------
void printLattice(const double phi[L][L]) {
  for(int x=0; x<L; x++){
    for(int y=0; y<L; y++)
      cout <<"("<<x<<","<<y<<") = " << phi[x][y] << endl;
    cout << "\n\n";
  }
  return;
}

void coldStart(double phi[L][L], param p) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      phi[x][y] = 1.0;
  return;
}

void hotStart(double phi[L][L],param p) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++) {
      phi[x][y] = 2.0*drand48() - 1.0;
      if(drand48() < 0.5) phi[x][y] = - phi[x][y];   
    }
  return;
}  

//v1: upgraded function to compute momenta on the device
void gaussReal_F(double field[L][L], double rands[L][L][2]) {
  //normalized gaussian exp[ - phi*phi/2]  <eta^2> = 1
#pragma acc data present(field[0:L][0:L], rands[0:L][0:L][0:2])
  {
    double r, theta;
#pragma acc parallel loop tile(32,2)
    for(int x=0; x<L; x++){
      for(int y=0; y<L; y++){
	r = sqrt(-2.0*log(rands[x][y][0]));
	theta = TWO_PI*rands[x][y][1];
	field[x][y] = r*cos(theta);
      }
    }
  }
  return;
}
//------------------------------------------------------------------
