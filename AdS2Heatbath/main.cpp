#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>

using namespace std;

//Library wide precision.
typedef long double Float;

#include "util.h"
#include "hyp_util.h"
#include "graph.h"

//If true, prints a bunch of numbers!
bool flags = false;

//Utilities
void copyLattice(vector<Vertex> &Lattice, bool new_to_old, param p);
void zeroField(Float *field, param p);
Float measH(vector<Vertex> Lattice, Float *mom, param p);
void gaussReal(Float *mom, param p);

//HMC routines
int hmc(vector<Vertex> &Lattice, param p, int iter);
void forceU(Float *fU, vector<Vertex> &Lattice, param p);
void update_mom(Float *mom, Float *fU, param p, double dt);
void update_phi(vector<Vertex> &Lattice, Float *mom, param p, double dt);
void trajectory(vector<Vertex> &Lattice, Float *mom, param p);

//Heatbath
void heatbath(vector<Vertex> &Lattice, param p, int iter);

Float dHAve = 0.0;
Float dHAve2 = 0.0;
Float dHAve3 = 0.0;
Float dHAve4 = 0.0;
Float dHAve5 = 0.0;
Float dHAve6 = 0.0;

Float sum = 0.0;

Float expdHAve = 0.0;
int accepted_metropolis = 0;
 
// Begin Main Program
//==============================================================
int main(int argc, char **argv) {
  
  param p;
  //Process Command line arguments
  for (int i=1; i<argc; i++){
    if(p.init(argc, argv, &i) == 0){
      continue;
    }
    printf("ERROR: Invalid option: %s\n", argv[i-1]);
    p.usage(argv);
    exit(0);
  }

  int k=0;
  if(argc > 1) p.init(argc, argv, &k);
  
  //Pseudo RNG seed
  srand48(p.seed);
  
  p.S1 = endNode(p.Levels,p) + 1 - (endNode(p.Levels-1,p) + 1);
  p.surfaceVol = p.S1*p.Lt;
  p.AdSVol = endNode(p.Levels,p) + 1;
  p.latVol = p.AdSVol;
  
  //Print paramters
  p.print();
  
  //Object to hold graph data
  vector<Vertex> Lattice(p.latVol);
  
  //-1 in Lattice indicates that node n has no connections.
  //This is used during construction to indicate if the node is yet
  //to be populated. During truncation, nodes to be removed are
  //assigned a position of -1, and all connections to that node are
  //removed.
  for(int n = 0; n < p.latVol; n++)
    for(int mu = 0; mu < p.q+2; mu++) 
      Lattice[n].nn[mu] = -1;
  
  //Construct neighbour table.
  buildGraph(Lattice, p);
  //Get the z-coords
  getComplexPositions(Lattice, p);

  //Begin (H)MC routines
  //---------------------------------------------------------------------------

  //Fake out
  //If we want to compute the lattice action, we need a momentum field.  
  Float mom[p.AdSVol];
  zeroField(mom, p);

  long int accept = 0;
  long int count = 0;
  
  //Initialise Lattice
  for(int i = 0; i<p.latVol; i++) Lattice[i].phi = 2.0*drand48() - 1.0;
  //copyLattice(Lattice, true, p);
  
  //Loop over thermalisation iterations
  for(int iter=0; iter<p.n_therm; iter++) {
    //accept = hmc(Lattice, p, iter);    
    heatbath(Lattice, p, iter);
  }

  //reset statistics
  accepted_metropolis = 0;
  //Loop over warmup iterations
  for(int iter=p.n_therm; iter<p.n_therm + p.n_meas; iter++) {

    count++;
    double Hold = measH(Lattice, mom, p);
    heatbath(Lattice, p, iter);
    double H = measH(Lattice, mom, p);
    //accept += hmc(Lattice, p, iter);
    
    cout << iter << " " << setprecision(16) << (1.0*accepted_metropolis)/(count*p.AdSVol) << " " << expdHAve/(count*p.AdSVol) << " ";
    cout << dHAve/(count*p.AdSVol) << endl;
    // cout << (2*dHAve2 - dHAve*dHAve)/(count*p.AdSVol) << endl;
    // cout << dHAve2/(count*p.AdSVol) << " ";
    // cout << dHAve3/(count*p.AdSVol) << " ";
    // cout << dHAve4/(count*p.AdSVol) << " ";
    // cout << dHAve5/(count*p.AdSVol) << " ";
    // cout << dHAve6/(count*p.AdSVol) << " ";
    // cout << sum/(count*p.AdSVol) << " ";
    // cout << 1 - dHAve + dHAve2/2 - dHAve3/6 + dHAve4/24 - dHAve5/120 + dHAve6/720 << endl;
  }
  
  return 0;
}

void heatbath(vector<Vertex> &Lattice, param p, int iter) {
  
  double phi_new = 0.0;
  double phi_new_sq = 0.0;
  double phi = 0.0;
  double phi_sq = 0.0;
  double lambda = 0.25*p.lambda;
  double msqr   = 0.50*p.msqr;
  
  double deltaH = 0.0;
  
  for (int i=0; i<p.latVol; i++) {
    
    deltaH = 0.0;
    
    phi = Lattice[i].phi;
    phi_new = phi + p.delta_phi * (2.0*drand48() - 1.0);
    phi_new_sq = phi_new*phi_new;
    phi_sq = phi*phi;
    
    //PE
    deltaH += lambda*(phi_new_sq*phi_new_sq - phi_sq*phi_sq);
    deltaH += msqr  *(phi_new_sq            - phi_sq);
    
    //KE
    for(int q=0; q<p.q; q++) {
      if(Lattice[i].nn[q] != -1) {
	deltaH += 0.5 * (phi_new_sq - phi_sq + 2*Lattice[Lattice[i].nn[q]].phi*(phi_new - phi));
      } else {
	deltaH += 0.5 * (phi_new_sq - phi_sq);
      }
    }
    
    
    if (iter >= p.n_therm) {
      expdHAve += exp(-deltaH);
      dHAve  += deltaH;
      // dHAve2 += pow(deltaH,2)/2;
      // dHAve3 += pow(deltaH,3)/6;
      // dHAve4 += pow(deltaH,4)/24;
      // dHAve5 += pow(deltaH,5)/120;
      // dHAve6 += pow(deltaH,6)/720;
      // sum += 1 - deltaH + pow(deltaH,2)/2 - pow(deltaH,3)/6 + pow(deltaH,4)/24 - pow(deltaH,5)/120 + pow(deltaH,6)/720;
    }
    
    if(deltaH < 0.0) {
      //cout << "Accepted" << endl;
      accepted_metropolis++;
      Lattice[i].phi = phi_new;
    }
    else if ( drand48() < exp(-deltaH)) {
      //cout<< "Acepted" << endl;
      accepted_metropolis++;
      Lattice[i].phi = phi_new;
    }
  }
}

// Utilities
//----------------------------------------------------------------------------

//Computes the action of the field
Float measH(vector<Vertex> Lattice, Float *mom, param p) {

  Float KE, PE, Hmom;
  Float phi_sq;
  Float phi;
  Float lambda_p = 0.25*p.lambda;
  Float msqr_p   = 0.50*p.msqr;

  for (int i=1; i<p.latVol; i++) {

    //Field action
    //------------
    phi = Lattice[i].phi;
    phi_sq = phi*phi;
    
    //PE terms
    PE += lambda_p * phi_sq*phi_sq;
    PE += msqr_p   * phi_sq;
    
    //Spatial: q=0 and q=fwdLinks+1 are on the same level. We take q=fwdLinks+1 and
    //the fwdLinks as the forward links in a forward difference operator
    //for(int q=0; q<p.q; q++) {
    for(int q=0; q<Lattice[i].fwdLinks+1; q++) {
      
      //Test if the link is connected
      if(Lattice[i].nn[q] != -1) {
	
	//Compute kinetic term
	KE += 0.5*((phi - Lattice[Lattice[i].nn[q]].phi)*
		   (phi - Lattice[Lattice[i].nn[q]].phi));
      } else {
	
	KE += 0.5*phi_sq;	
      }
    }
    
    //Momentum action
    //---------------
    Hmom += 0.5 * mom[i]*mom[i];
  }
  
  if (flags) cout << KE << " " << PE << " " << Hmom << " ";
  return KE + PE + Hmom;
}

//normalized gaussian exp[-phi*phi/2] |  <phi^2> = 1
void gaussReal(Float *field, param p) {  
  double r, theta, sum;
  for(int i=0; i<p.AdSVol; i++) {
    r = sqrt(-2.0*log(drand48()));
    theta = 2*M_PI * drand48();
    field[i] = r*cos(theta);
  }
  
  return;
}

void copyLattice(vector<Vertex> &Lattice, bool new_to_old, param p){
  
  if (new_to_old) {    
    for (int i=0; i<p.AdSVol; i++) Lattice[i].phi_old = Lattice[i].phi;
  } else {
    for (int i=0; i<p.AdSVol; i++) Lattice[i].phi = Lattice[i].phi_old;
  }
}

void zeroField(Float *field, param p){
  
  for (int i=0; i<p.AdSVol; i++) field[i] = 0.0;
}







// HMC routines
//----------------------------------------------------------------------------

int hmc(vector<Vertex> &Lattice, param p, int iter) {

  int accepted = 0;
  Float Hold = 0, H = 0;
  
  Float mom[p.AdSVol];
  gaussReal(mom, p);

  // MD trajectory using Verlet (Leapfrog)
  Hold = measH(Lattice, mom, p); if(flags) cout << endl;
  trajectory(Lattice, mom, p);
  H = measH(Lattice, mom, p); if(flags) cout << endl;
  
  if (iter+1 > p.n_therm) {
    
    // Metropolis accept/reject step.

    if( (H-Hold) != (H-Hold)) {
      cout << "Nope..." << endl;
      exit(0);
    } else {
      cout << "deltaH = " << (H-Hold) << endl;
    }
    
    dHAve += (H-Hold);
    if ( drand48() > exp(-(H-Hold)) ) {
      //Rejected
      //cout << "Rejecting " << exp(-(H-Hold)) << endl;
      copyLattice(Lattice, false, p);
    }
    else {
      //Accepted
      //cout << "Accepting " << exp(-(H-Hold)) << endl;
      accepted = 1;
      copyLattice(Lattice, true, p);
    }
  } else {
    //Auto accept
    copyLattice(Lattice, true, p);    
  } 
  return accepted;
  
}

void trajectory(vector<Vertex> &Lattice, Float *mom, param p) {

  const int n_step = p.n_step;
  const Float dt = p.tau/n_step;

  Float fU[p.AdSVol];
  
  if (flags) cout << 0 << " " << measH(Lattice, mom, p) << endl;
  
  // Implement the Leapfrog method
  //Initial half step:
  //P_{1/2} = P_0 - dtau/2 * fU  
  forceU(fU, Lattice, p);
  update_mom(mom, fU, p, 0.5*dt);
  if (flags) cout << 0.5 << " " << measH(Lattice, mom, p) << endl;
  //step loop
  for(int k=1; k<n_step; k++) {
    
    //U_{k} = U_{k-1} + P_{k-1/2} * dt
    update_phi(Lattice, mom, p, dt);
    
    //P_{k+1/2} = P_{k-1/2} - fU * dt
    forceU(fU, Lattice, p);
    update_mom(mom, fU, p, dt);

    if (flags) cout << k << " " << measH(Lattice, mom, p) << endl;
    
  } //end step loop
  
  //Final half step.
  //U_{n} = U_{n-1} + P_{n-1/2} * dt
  update_phi(Lattice, mom, p, dt);
  forceU(fU, Lattice, p);
  update_mom(mom, fU, p, 0.5*dt);

  if (flags) cout << "final" << " " << measH(Lattice, mom, p) << endl;
}

void forceU(Float *fU, vector<Vertex> &Lattice, param p) {

  zeroField(fU, p);
  const Float msqr  = 0.5*p.msqr;
  const Float lambda = 0.25*p.lambda;
  Float phi_lc = 0.0;
  
  for(int i=1; i<p.AdSVol; i++) {
    
    //if( i < endNode(p.Levels-1, p) +1) {
      
      //A convenience
      phi_lc = Lattice[i].phi;
      //cout << i << " ";
      for (int q=0; q<p.q; q++) {
	//cout << Lattice[i].nn[q] << " ";
	if(Lattice[i].nn[q] != -1) {
	  //cout << Lattice[i].nn[q] << " ";
	  fU[i] -= 1.0*(Lattice[Lattice[i].nn[q]].phi - phi_lc);
	}
	else {
	  fU[i] += 1.0*phi_lc;
	}
      }
      //cout << endl;
      fU[i] += (2.0*msqr + 4.0*lambda*phi_lc)*phi_lc*phi_lc;
    }
  //}  
}

void update_mom(Float *mom, Float *fU, param p, double dt) {

  for(int i=0; i<p.AdSVol; i++) mom[i] -= fU[i] * dt;
  
}

void update_phi(vector<Vertex> &Lattice, Float *mom, param p, double dt) {
  
  for(int i=0; i<p.AdSVol; i++) Lattice[i].phi += mom[i] * dt;
  
}
