//================================================================================= 
// This is a 2d phi 4th code on a torus desinged to be easy to convet to openACC.
//=================================================================================

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
 
#define TWO_PI 6.283185307179586
#define L 64

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
} param;

//Utilities
void printLattice(const double phi[L][L]);
void hotStart(double phi[L][L], param p);
void coldStart(double phi[L][L], param p);
void gaussMom(double mom[L][L]);

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
bool swendsenWang( int label[L][L],const bool bond[L][L][4]);
void flipSpins(double phi[L][L], const int label[L][L]);

// Global variables for HMC diagnostics
//--------------------------------------------------------------------------
double ave_expmdH = 0.0;
double ave_dH = 0.0;

// Inlined Utilities
//--------------------------------------------------------------------------
// Zero lattice field.
template<typename T> inline void zeroField(T phi[L][L]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      phi[x][y] = 0.0;
}

// Copy lattice field
template<typename T> inline void copyField(T phi2[L][L],T phi1[L][L]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      phi2[x][y] = phi1[x][y];
}
//--------------------------------------------------------------------------



int main(int argc, char **argv) {
  
  param p;
  p.Latsize = L;
  p.lambda = 0.125; 
  p.musqr = -0.359;
  p.nStep = 15;
  p.tau = 1.0;
  p.dt = p.tau/p.nStep;
  p.nTherm = 1000;
  p.nIter = 10000;
  p.nMeas = 25;
  p.nHMC = 5;
  p.verbose = true;
  p.debug = false;
  p.seed = 1234;
  p.relaxIter = 1000;
  
  srand(p.seed);

  cout << "lattice size:" << L << endl;
  cout << "lambda = " << p.lambda << endl;
  cout << "musqr = " << p.musqr << endl;
  cout << "trajectory steps = "<< p.nStep << endl;
  cout << "trajectory length = " << p.tau << endl;
  cout << "trajectory step length = " << p.dt << endl;
  cout << "thermalisation iters = " << p.nTherm << endl;
  cout << "total iterations = " << p.nIter << endl;
  cout << "HMC performed every " << p.nHMC << " iters" << endl;
  cout << "relaxation sweep max = " << p.relaxIter << endl;
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
  
  double momZero[L][L];
  zeroField(momZero);


  
  //--------------------------------------------------------------------------
  // Thermalise the field. This evolves the phi field
  // from a completely random (hot start)
  // or a completely ordered (cold start)
  // to a state of thermal equlibrium.
  for(int iter = 0; iter < p.nTherm; iter++) {
    
    if(iter%p.nHMC == 0) hmc(phi, p, iter);
    
    latticePercolate(bond, label, phi);
    
    stop = false;
    for(int relax = 0; relax < p.relaxIter && !stop; relax++) {
      stop = swendsenWang(label, bond);
      if(stop) {
	//Lattice bonds identified. Flip 'em!
	flipSpins(phi, label);
	if(p.debug) cout << "Iter " << iter << ": relaxation sweeps: " << relax << endl; 
      }
    }
    if(!stop) {
      cout << "Error in SW: clusters not identified in " << p.relaxIter << " relaxation sweeps. Exiting..." << endl;
      exit(0);
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
    for(int relax = 0; relax < p.relaxIter && !stop; relax++) {
      stop = swendsenWang(label, bond);
      if(stop) {
	//Lattice bonds identified. Flip 'em!
	flipSpins(phi, label);
	if(p.debug) cout << "Iter " << iter << ": relaxation sweeps: " << relax << endl; 
      }
    }
    if(!stop) {
      cout << "Error in SW: clusters not identified in " << p.relaxIter << " relaxation sweeps. Exiting..." << endl;
      exit(0);
    }
    
    // Take measurements
    //--------------------------------------------------------------------------
    if((iter+1)%p.nMeas == 0) {
      
      measurement++;
      
      //Get phi mean and moments
      //------------------------------------------------------------------------
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
	cout << "MEAS |phi|          = " << setprecision(12) << AbsPhi/(vol*measurement) << endl;
	cout << "MEAS <phi>          = " << setprecision(12) << avPhi/vol << endl;
	cout << "MEAS <phi**2>       = " << setprecision(12) << avPhi2/vol2 << endl;
	cout << "MEAS <phi**4>       = " << setprecision(12) << avPhi4/vol4 << endl;
	cout << "Binder Cumulant     = " << setprecision(12) << 1.0 - Phi4/(3.0*Phi2*Phi2/measurement) << endl;
      }
      //--------------------------------------------------------------------------

      // Dump data to file
      //--------------------------------------------------------------------------
      string filename;
      char fname[256];
      FILE *fp;
      filename = "phi4Data_L" + to_string(p.Latsize) + "_lambda" + to_string(p.lambda) +  "_musqr" + to_string(p.musqr) + "_tau"+ to_string(p.tau) +"_nHMCStep" + to_string(p.nStep) + "_nHMC" + to_string(p.nHMC) + "_seed" + to_string(p.seed) + ".dat";

      sprintf(fname, "%s", filename.c_str());	
      fp = fopen(fname, "a");	
      fprintf(fp, "%d %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
	      iter+1,
	      (double)accepted/nHMC,
	      ave_expmdH/nHMC,
	      ave_dH/nHMC,
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
  return 0;
}



//Cluster decomposition
//----------------------------------------------------------------------------
// Identify the bond structure of the Swendsen-Wang decomposition.
void latticePercolate(bool bond[L][L][4], int label[L][L], const double phi[L][L]) {
  
  double probability; 

  //Initialise all bonds to false 
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      for(int mu  = 0; mu < 4;mu++)
	bond[x][y][mu] = false;

  //For every lattice site there are two bonds, one each in the
  //+ve x- and y-dimension.
  for(int x=0 ; x<L; x++)
    for(int y=0; y<L; y++){

      //Test if x the bond is active
      probability = 1.0 -  exp( -2.0 * phi[x][y]*phi[(x+1)%L][y]);
      if (drand48() < probability) {
	bond[x][y][0] = true;
	bond[(x+1)%L][y][2] = true;
      }

      //Test if y the bond is active
      probability = 1.0 -  exp( -2.0 * phi[x][y]*phi[x][(y+1)%L]);
      if(drand48() < probability) {
	bond[x][y][1] = true;
	bond[x][(y+1)%L][3] = true;
      }
      
      //Place a random spin on each lattice site, will use this in flipSpins()      
      if(drand48() < 0.5) label[x][y] = - (y + x*L + 1);
      else label[x][y] = (y + x*L + 1);
    }
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
bool swendsenWang(int label[L][L], const bool bond[L][L][4]) {
  
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
  
  return stop;
}

//The clusters will all have the same spin as the spin of the lowest
//coordinate phi value. We randomly placed a minus sign on the label
//which will now serve as a 50/50 test to flip the entire cluster
void flipSpins(double phi[L][L], const int label[L][L]) {
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++) {
      if(label[x][y] < 0) phi[x][y] = -phi[x][y];
    }  
  return;
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
  
  double phiOld[L][L];
  double H = 0.0, Hold = 0.0;

  //Copy the field in case of rejection
  copyField(phiOld, phi);

  
  //Step 1: Create gaussian distributed momenta
  double mom[L][L];
  zeroField(mom);  
  gaussMom(mom); 

  //Step 2: Perform trajectory
  Hold = calcH(mom, phi, p);
  trajectory(mom, phi, p); // MD trajectory using Leapfrog
  H = calcH(mom, phi, p);

  //record HMC diagnostics
  if(iter >= p.nTherm) {
    ave_expmdH += exp(-(H-Hold));
    ave_dH += H-Hold;
  }
  
  // Step 3: Metropolis accept/reject step
  // Always accepts trajectories during first half of warm up.
  if (drand48() > exp(-(H-Hold)) && iter > p.nTherm/2-1) {    

    //Keep old field
    copyField(phi, phiOld);
    return 0;
  }
  else {
    
    //Accept candidate field
    return 1;
  }
}

//Step 2: Performs the HMC trajectory via Leapfrog integration
void trajectory(double mom[L][L], double phi[L][L], param p) {

  const double dt = p.dt;
  double fphi[L][L];
  //Step 2.1: Initial half step:
  //P_{1/2} = P_0 - dtau/2 * fphi
  force_phi(fphi, phi, p);
  update_mom(mom, fphi, p, 0.5*dt);
  
  //Step 2.2:
  for(int k=1; k<p.nStep; k++) {
    
    //phi_{k} = phi_{k-1} + P_{k-1/2} * dtau
    update_phi(phi, mom, p, dt);
    
    //P_{k+1/2} = P_{k-1/2} - fphi * dtau
    force_phi(fphi, phi, p);
    update_mom(mom, fphi,  p, dt);
    
  }
  
  //step 2.3: Final half step.
  //phi_{n} = phi_{n-1} + P_{n-1/2} * dtau
  update_phi(phi, mom, p, dt);
  force_phi(fphi, phi, p);
  update_mom(mom, fphi, p, 0.5*dt);
  
  return;
}

//Compute the HMC force from the field
void force_phi(double fphi[L][L], double phi[L][L], param p) {
  
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++) {
	fphi[x][y] = 0.0;
	fphi[x][y] -= phi[(x+1)%L][y] - 2.0*phi[x][y] + phi[(x-1+L)%L][y];
	fphi[x][y] -= phi[x][(y+1)%L] - 2.0*phi[x][y] + phi[x][(y-1+L)%L];
	fphi[x][y] += (2.0*p.musqr + 4.0*p.lambda*phi[x][y]*phi[x][y]) * phi[x][y];
      }
    
}

//Update the momenta from the back reaction of the field
void update_mom(double mom[L][L], double fphi[L][L], param p, double dt) {
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++) 
      mom[x][y] -= fphi[x][y] * dt;
  
}

//Update the field from the back reaction of the momenta
void update_phi(double phi[L][L], double mom[L][L], param p, double dt) {
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++) 
      phi[x][y] += mom[x][y] * dt;
  
}
//---------------------------------------------------------------------------------



// Measurements
//------------------------------------------------------------------------------------
//Compute the action of the field and momenta
double calcH(double mom[L][L], double phi[L][L],  param p) {
  
  double Hphi = 0.0, Hmom = 0.0;  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++) {
      
      Hphi += 0.5*(phi[x][y] - phi[(x+1)%L][y])*(phi[x][y] - phi[(x+1)%L][y]);
      Hphi += 0.5*(phi[x][y] - phi[x][(y+1)%L])*(phi[x][y] - phi[x][(y+1)%L]);
      Hphi += (p.musqr + p.lambda*phi[x][y]*phi[x][y])*phi[x][y]*phi[x][y];
      
      Hmom += 0.5 * mom[x][y] * mom[x][y];
    }
  
  return Hphi + Hmom;
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
  for(int x =0;x< L;x++){
    for(int y =0;y< L;y++)
      cout <<"("<<x<<","<<y<<") = " << phi[x][y] << endl;
    cout << "\n\n";
  }
  return;
}

void coldStart(double phi[L][L], param p) {
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      phi[x][y] = 1.0;
  return;
}

void hotStart(double phi[L][L],param p) {
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++) {
      phi[x][y] = 2.0*drand48() - 1.0;
      if(drand48() < 0.5) phi[x][y] = - phi[x][y];   
    }
  return;
}  

//Step 1: Gaussian momenta generated here
void gaussMom(double field[L][L]) {
  //normalized gaussian exp[ - phi*phi/2]  <eta^2> = 1
  double r, theta;
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++){
      r = sqrt( -2.0*log(drand48()) );
      theta = TWO_PI*drand48();
      field[x][y] = r*cos(theta);
    }
  return;
}
//------------------------------------------------------------------
