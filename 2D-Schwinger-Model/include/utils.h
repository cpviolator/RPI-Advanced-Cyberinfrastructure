#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>
#include "latHelpers.h"
#include "fermionHelpers.h"

using namespace std;

typedef struct{
  
  //HMC
  int nstep = 25;
  double tau = 1.0;
  int iterHMC = 1000;
  int therm = 50;
  int skip = 25;
  int checkpoint = 100;
  int checkpointStart = 0;
  int maxIterCG = 1000;
  double eps = 1e-6;
  
  //physics
  int XLatsize = LX;
  int YLatsize = LX;
  double beta = 3.0;
  double m = -0.06;
  bool dynamic = true;
  
  //Smearing
  double alpha = 0.5;
  int smearIter = 1;
  
  //Measurements
  bool measPL = false; //Polyakov loops
  bool measWL = false; //Wilson loop and Creutz ratios
  bool measPC = false; //Pion
  bool measVT = false; //Vacuum trace
  
  //Wilson loop and Polyakov loop max size.
  int loopMax = LX/2;
  
} param_t;


void printParams(param_t p) {
  cout << endl;
  cout << "Physics:  XSize = "<< LX << endl;
  cout << "          YSize = "<< LY << endl;
  cout << "          Beta = "<< p.beta << endl;
  cout << "          Dynamic = " << (p.dynamic == true ? "True" : "False") << endl;
  if (p.dynamic == true) cout << "          Mass = " << p.m << endl;
  cout << "HMC:      Therm Sweeps: (" << p.therm << " accept) (" << p.therm << " accept/reject)" << endl; 
  cout << "          Data Points = " << p.iterHMC << endl;
  cout << "          Time Step = " << p.tau/p.nstep << endl;
  cout << "          Trajectory Steps " << p.nstep << endl;
  cout << "          Trajectory Length = " << p.tau << endl;
  cout << "Smearing: APE iter = " << p.smearIter << endl;
  cout << "          APE alpha = " << p.alpha << endl;
}

void constructName(string &name, param_t p) {
  name += "_LX" + to_string(p.XLatsize) + "_LY" + to_string(p.YLatsize) + "_B" + to_string(p.beta);
  if(p.dynamic == true) name += "_M"+ to_string(p.m);
  name += "_tau" + to_string(p.tau) + "_nHMCstep" + to_string(p.nstep);
}

void printLattice(Complex gauge[LX][LY][2]){
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++) {
      cout << "(x,y) = " << x << "," << y << " ";
      cout << gauge[x][y][0]<< " " << gauge[x][y][1] << endl;
    }
  return;
}

double measPlaq(Complex gauge[LX][LY][2]){
  
  double plaq = 0.0;  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++){
      plaq += real(gauge[x][y][0]*gauge[ (x+1)%LX ][y][1]*conj(gauge[x][ (y+1)%LY ][0])*conj(gauge[x][y][1]));
    }
  return plaq/(LX*LY);
}


void writeGaugeLattice(Complex gauge[LX][LY][2], string name){

  fstream outPutFile;
  outPutFile.open(name,ios::in|ios::out|ios::trunc);  
  outPutFile.setf(ios_base::fixed,ios_base::floatfield); 

  //Plaquette action header
  outPutFile << setprecision(20) <<  setw(20) << measPlaq(gauge) << endl;
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int mu=0; mu<2; mu++)
	outPutFile << setprecision(12) <<  setw(20) << arg(gauge[x][y][mu]) << endl;
  
  outPutFile.close();
  return;
  
}
void readGaugeLattice(Complex gauge[LX][LY][2], string name){

  fstream inPutFile;
  inPutFile.open(name);
  string val;
  if(!inPutFile.is_open()) {
    cout << "Error opening file " << name << endl;
    exit(0);
  }
  
  //Header check
  getline(inPutFile, val);
  double plaq = stod(val);

  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int mu=0; mu<2; mu++){
	getline(inPutFile, val);
	gauge[x][y][mu] = polar(1.0, stod(val));	  
      }
  cout << setprecision(16) << setw(20) << "Plaqette on file  = " << plaq << endl;
  cout << setprecision(16) << setw(20) << "Plaqette measured = " << measPlaq(gauge) << endl;
  double err = fabs(1.0 - plaq/measPlaq(gauge));
  if(err > 1e-12) {
    cout << "Gauge read fail!" <<  endl;
    exit(0);
  }    
  return;
}

/*===============================================================================
  Gaussian numbers with p(theta) = sqrt(beta/ 2 PI) exp( - beta* theta^2/2)
  <Gaussian^2> = 1/beta  
  Perimeter Law:  Wilson Loop = exp[ - 4 sigma L ]   sigma = - Log[ <cos(theta)> ]
  ================================================================================*/ 
void gaussStart(Complex gauge[LX][LY][2], param_t p){

  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++){
      gauge[x][y][0] = polar(1.0,sqrt(1.0/p.beta)*drand48());
      gauge[x][y][1] = polar(1.0,sqrt(1.0/p.beta)*drand48());
    }
  return;
}  

void coldStart(Complex gauge[LX][LY][2], param_t p){

  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int mu=0; mu<2; mu++)
	gauge[x][y][mu] = Complex(1.0,0.0);
  return;
}  

void gaussReal_F(double field[LX][LY][2]) {
  //normalized gaussian exp[ - phi*phi/2]  <phi|phi> = 1
  double r, theta, sum;
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++){
      r = sqrt(-2.0*log(drand48()));
      theta = TWO_PI*drand48();
      field[x][y][0] = r*cos(theta);
      
      r = sqrt(-2.0*log(drand48()));
      theta = TWO_PI*drand48();
      field[x][y][1] = r*cos(theta);
      //sum += field[x][y][0]*field[x][y][0] + field[x][y][1]*field[x][y][1];
    }
  
  //cout << "mom check: " << sum/(L*L) << endl;
  
  return;
}

void gaussReal_F(double field[LX][LY]) {
  //normalized gaussian exp[ - phi*phi/2]  <phi|phi> = 1
  double r, theta, sum;
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++){
      r = sqrt(-2.0*log(drand48()));
      theta = TWO_PI*drand48();
      field[x][y] = r*cos(theta);

      //sum += field[x][y]*field[x][y];
    }
  
  //cout << "mom check: " << sum/(L*L) << endl;
  
  return;
}


void gaussComplex_F(Complex eta[LX][LY][2], param_t p) {
  
  //normalized gaussian exp[ - eta*eta/2]  <eta|eta> = 1;
  double r1, theta1, r2, theta2, sum;
  double inv_sqrt2 = 1.0/sqrt(2);
  
  for(int x=0; x<LX; x++) {
    for(int y=0; y<LY; y++) {
      for(int s=0; s<2; s++) {
	r1 = sqrt(-2.0*log(drand48()));
	theta1 = TWO_PI*(drand48());
	r2 = sqrt(-2.0*log(drand48()));
	theta2 = TWO_PI*(drand48());
	
	eta[x][y][s] = Complex(r1*cos(theta1),r2*sin(theta2))*inv_sqrt2;
      }
    }
  }
  //cout << "GaussComplex_F: norm(eta) = " << norm2(eta)/(LX*LY*2) << endl;
  return;
}

void gaussComplex_F(Complex eta[LX][LY], param_t p) {
  
  //normalized gaussian exp[ - eta*eta/2]  <eta|eta> = 1;
  double r1, theta1, r2, theta2;
  double inv_sqrt2 = 1.0/sqrt(2);
  
  for(int x=0; x<LX; x++) {
    for(int y=0; y<LY; y++) {
      r1 = sqrt(-2.0*log(drand48()));
      theta1 = TWO_PI*(drand48());
      r2 = sqrt(-2.0*log(drand48()));
      theta2 = TWO_PI*(drand48());
      
      eta[x][y] = Complex(r1*cos(theta1),r2*sin(theta2))*inv_sqrt2;
    }
  }
  //cout << "GaussComplex_F: norm(eta) = " << norm2(eta)/(LX*LY) << endl;
  return;
  
}

//staple x is 0th, y is 1st.
//APE smearing: project back on U(1)       
void smearLink(Complex Smeared[LX][LY][2], const Complex gauge[LX][LY][2], param_t p){

  double alpha = p.alpha;
  int iter = p.smearIter;
  int xp1, xm1, yp1, ym1;
  
  Complex SmearedTmp[LX][LY][2];
  copyLat(Smeared, gauge);
  copyLat(SmearedTmp, Smeared);
  
  for(int i=0; i<iter; i++) {    
    for(int x=0; x<LX; x++) {
      xp1 = (x+1)%LX;
      xm1 = (x-1+LX)%LX;
      for(int y=0; y<LY; y++) {
	yp1 = (y+1)%LY;
	ym1 = (y-1+LY)%LY;
	
	SmearedTmp[x][y][0] += alpha * Smeared[x][y][1] * Smeared[x][yp1][0] * conj(Smeared[x][y][1]);
	SmearedTmp[x][y][0] += alpha * conj(Smeared[x][ym1][1]) * Smeared[x][ym1][0] * Smeared[xp1][ym1][1];
	SmearedTmp[x][y][1] += alpha * Smeared[x][y][0] * Smeared[xp1][y][1] * conj(Smeared[xp1][yp1][0]);
	SmearedTmp[x][y][1] += alpha * conj(Smeared[xm1][y][0]) * Smeared[xm1][y][1] * Smeared[xm1][yp1][0];
      }
    }
    
    //Project back to U(1)
    for(int x=0; x<LX; x++)
      for(int y=0; y<LY; y++)
	for(int mu=0; mu<2; mu++)
	  Smeared[x][y][mu] = polar(1.0,arg(SmearedTmp[x][y][mu]));
  }
}

void usage(char **argv) {
  
  printf("This is an exhaustive list of run-time options. Those not set at\n");
  printf("run time will be given a default value. Important input options will\n");
  printf("dumped to stdout at execution. Please ensure that they are sensible!\n\n");

  // The user is advised to study these parameters and set them accordingly.
  //-----------------------------------------------------------------------

  //HMC params
  //------------------------------------------------------------------------
  printf("--nStep <int>           The number of leapfrog steps in the HMC.\n");
  printf("--tau <double>          The length of the HMC trajectory.\n");
  printf("--iterHMC <int>         The number of HMC trajectories to compute.\n");
  printf("--therm <int>           The number of HMC trajectories the thermalise.\n");
  printf("--skip <int>            The number of HMC trajectories to skip before taking a measurement\n");  
  printf("--checkpoint <int>      Save the gauge field every <checkpoint>th HMC iter.\n");
  printf("--checkpointStart <int> Start the HMC from this trajectory.\n");
  printf("--maxIterCG <int>       The maximum number of iteratons CG will run.");
  printf("--tolCG <double>        Desired residual norm of CG solver.\n");
  
  //Physics
  //------------------------------------------------------------------------
  printf("--beta <double>         The beta value (1/g^2).\n");
  printf("--mass <double>         The fermion mass.\n");
  printf("--dynamic  <true/false> Dynamic or quenched simulation\n");  
  printf("--smearAlpha <double>   The alpha parameter in APE smearing.\n");
  printf("--smearIter <int>       The number of APE smearing hits to perform.\n");
  
  //Measurements
  //------------------------------------------------------------------------
  printf("--measPL <true/false>   Measure Polyakov loops.\n");
  printf("--measWL <true/false>   Measure Wilson loops.\n");
  printf("--measPC <true/false>   Measure Pion correlation.\n");
  
}

int init(int argc, char **argv, int *idx, param_t *p) {
  
  int ret = -1;
  int i = *idx;
  
  if( strcmp(argv[i], "--help")== 0){
    usage(argv);
    exit(0);
  }
  
  //HMC params
  //------------------------------------------------------------------------
  if( strcmp(argv[i], "--nSteps") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    
    p->nstep = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  if( strcmp(argv[i], "--tau") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    
    p->tau = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

 if( strcmp(argv[i], "--iterHMC") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    
    p->iterHMC = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }


  if( strcmp(argv[i], "--therm") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    
    p->therm = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  if( strcmp(argv[i], "--skip") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    
    p->skip = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }
  
  if( strcmp(argv[i], "--checkpoint") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    
    p->checkpoint = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }
  
  if( strcmp(argv[i], "--checkpointStart") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    
    p->checkpointStart = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  if( strcmp(argv[i], "--maxIterCG") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    
    p->maxIterCG = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  if( strcmp(argv[i], "--beta") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    
    p->beta = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }
  
  if( strcmp(argv[i], "--mass") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    
    p->m = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }
  
  if( strcmp(argv[i], "--dynamic") == 0){
    if (i+1 >= argc){
      usage(argv);
    }

    if (strcmp(argv[i+1], "true") == 0){
      p->dynamic = true;
    }else if (strcmp(argv[i+1], "false") == 0){
      p->dynamic = false;
    }else{
      fprintf(stderr, "ERROR: invalid dynamic type (true/false)\n");
      exit(1);
    }
    
    i++;
    ret = 0;
    goto out;
  }
  
  if( strcmp(argv[i], "--smearAlpha") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    
    p->alpha = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }


  if( strcmp(argv[i], "--smearIter") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    
    p->smearIter = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  if( strcmp(argv[i], "--measPL") == 0){
    if (i+1 >= argc){
      usage(argv);
    }

    if (strcmp(argv[i+1], "true") == 0){
      p->measPL = true;
    }else if (strcmp(argv[i+1], "false") == 0){
      p->measPL = false;
    }else{
      fprintf(stderr, "ERROR: invalid measPL type (true/false)\n");
      exit(1);
    }
    
    i++;
    ret = 0;
    goto out;
  }

  if( strcmp(argv[i], "--measWL") == 0){
    if (i+1 >= argc){
      usage(argv);
    }

    if (strcmp(argv[i+1], "true") == 0){
      p->measWL = true;
    }else if (strcmp(argv[i+1], "false") == 0){
      p->measWL = false;
    }else{
      fprintf(stderr, "ERROR: invalid measWL type (true/false)\n");
      exit(1);
    }
    
    i++;
    ret = 0;
    goto out;
  }

  if( strcmp(argv[i], "--measPC") == 0){
    if (i+1 >= argc){
      usage(argv);
    }

    if (strcmp(argv[i+1], "true") == 0){
      p->measPC = true;
    }else if (strcmp(argv[i+1], "false") == 0){
      p->measPC = false;
    }else{
      fprintf(stderr, "ERROR: invalid measPC type (true/false)\n");
      exit(1);
    }
    
    i++;
    ret = 0;
    goto out;
  }
  
  
  
 out:
  *idx = i;
  return ret ;
}

#endif
