#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <string>
#include <random>
#include <unistd.h>

#include "util.h"
#include "hyp_util.h"

using namespace std;

void param::usage(char **argv) {

  printf("\n\n");
  printf("This is an exhaustive list of run-time options. Those not set at\n");
  printf("run time will be given a default value. Important input options will\n");
  printf("dumped to stdout at execution. Please ensure that they are sensible!\n\n");

  // The user is advised to study these parameters and set them accordingly.
  //-----------------------------------------------------------------------
  printf("--lambda <Float>                 The phi^4 interaction term.\n");
  printf("--mSqr <Float>                   The AdS mass.\n");

  //AdS params
  printf("--q <n>                          The triangulation of the Poincare disk.\n");
  printf("--levels <n>                     The number of levels to which the Poincare Disk is generated.\n");
  printf("--centre <V/v, C/c>              Place a vertex at the centre of the Poincare disk, or the centre of a hyperbolic triangle.\n");  
  
  //MC params
  printf("--nTherm <n>                     The number of thermalisation steps (1 x metro + nCluster x Cluster)\n");
  printf("--nSkip  <n>                     The number of samples to skip between measuremnets.\n");
  printf("--nMeas  <n>                     The number of measurements to make.\n");
  printf("--nStep  <n>                     The number of HMC steps to make.\n");
  printf("--tau  <Float>                   The HMC step size.\n");  
  printf("--deltaPhi <float>               The parameter that governs how much an individual phi value is changed in the Metropolis test.\n");

  printf("--seed <long int>                Sets the seed for the RNG.\n");  
  printf("--verbosity <V/v, Q/q>           Sets the verbosity as either Verbose or Quiet.\n");  
  
}

int param::init(int argc, char **argv, int *idx) {

  int ret = -1;
  int i = *idx;

  if( strcmp(argv[i], "--help")== 0){
    usage(argv);
  }

  //Theory
  if( strcmp(argv[i], "--theory") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string theory_in(argv[i+1]);
    if (theory_in == "phi4" ||
	theory_in == "PHI4" ||
	theory_in == "Phi4") {
      theory_type = PHI4;
      
    } else if(theory_in == "ising" ||
	      theory_in == "Ising" ||
	      theory_in == "ISING") {
      theory_type = ISING;
    } else {
      cout<<"Invalid theory type ("<<theory_in<<") given. Options are "<<endl;
      cout<<"phi4: Phi fourth theory"<<endl;
      cout<<"ising: The Ising model."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }
  
  //Coupling Type
  if( strcmp(argv[i], "--couplingType") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    std::string couplingType_in(argv[i+1]);    
    if (couplingType_in == "sr"  ||
	couplingType_in == "SR"  ||
	couplingType_in == "Sr") {
      coupling_type = SR;      
    } else if(couplingType_in == "pow" ||
	      couplingType_in == "POW" ||
	      couplingType_in == "Pow") {
      coupling_type = POW;
      usePowLaw = true;
    } else if(couplingType_in == "rad" ||
	      couplingType_in == "RAD" ||
	      couplingType_in == "Rad") {
      coupling_type = RAD;
      usePowLaw = false;      
    } else {
      cout<<"Invalid coupling type "<<couplingType_in<<" given. Options are "<<endl;
      cout<<"SR: Short range (local) interacton"<<endl;
      cout<<"Pow: Long range (non-local) |x-y|^{d+s}"<<endl;
      cout<<"Rad: Long range (non-local) (cosh(dt) - cos(dtheta))^{d+s}"<<endl;    
      exit(0);
    }

    i++;
    ret = 0;
    goto out;
  } 

  //Boundary (lambda) interaction
  if( strcmp(argv[i], "--lambda") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    lambda = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  //Long-range coupling exponent
  if( strcmp(argv[i], "--sigma") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    sigma = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  //Ising (constant) coupling
  if( strcmp(argv[i], "--J") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    J = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  //Ising external field
  if( strcmp(argv[i], "--h") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    h = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  //Temporal scaling
  if( strcmp(argv[i], "--tScale") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    t_scale = atof(argv[i+1]);
    if(t_scale < 0) {
      cout<<"Invalid t_scale ("<<t_scale<<") given. Please ensure that tscale > 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //----------------//
  //   AdS params   //
  //----------------// 

  //Triangulation
  if( strcmp(argv[i], "--q") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    q = atoi(argv[i+1]);
    if(q < 7) {
      cout<<"Invalid triangulation ("<<q<<") given. Please ensure that q >= 7."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }
  
  //Levels
  if( strcmp(argv[i], "--levels") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    Levels = atoi(argv[i+1]);
    if(Levels < 3) {
      cout<<"Invalid number of levels ("<<Levels<<") given. Please ensure that Levels >= 3."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }
  
  //Vertex or body centred.
  if( strcmp(argv[i], "--centre") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string Centre(argv[i+1]);
    if (Centre == "V" || Centre == "v") {
      Vcentre = true;
    } else if (Centre == "C" || Centre == "c") {
      Vcentre = false;
    } else {
      cout<<"Invalid centre condition ("<<Centre<<") given. Use V/v for Vertexcentred or C/c for Circumcentred."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //Bulk (M) mass^2
  if( strcmp(argv[i], "--mSqr") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    msqr = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  //Bulk (M) mass^2 increment (Multishift inverter)
  if( strcmp(argv[i], "--deltaMsqr") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    delta_msqr = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  //Levels
  if( strcmp(argv[i], "--nShift") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    n_shift = atoi(argv[i+1]);
    if(n_shift < 0) {
      cout<<"Invalid number of shifts ("<<n_shift<<") given. Please ensure that n_shift >= 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }
  
  //Boundary Conditions
  if( strcmp(argv[i], "--BC") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    std::string BC(argv[i+1]);
    if (BC == "D" || BC == "d") {
      bc = true;
    } else if (BC == "N" || BC == "n") {
      bc = false;
    } else {
      cout<<"Invalid boundary condition ("<<BC<<") given. Use D/d for Dirichlet or N/n for Neumann."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //Maximum CG iterations
  if( strcmp(argv[i], "--maxIter") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    MaxIter = atoi(argv[i+1]);
    if(n_shift < 0) {
      cout<<"Invalid number of maxIter ("<<MaxIter<<") given. Please ensure that maxIter > 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //CG tolerance
  if( strcmp(argv[i], "--tol") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    tol = atof(argv[i+1]);
    if(tol < 0) {
      cout<<"Invalid tolerance ("<<tol<<") given. Please ensure that tol > 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //Source position
  if( strcmp(argv[i], "--srcPos") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    src_pos = atoi(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }


  //Bulk mass scaling factor
  if( strcmp(argv[i], "--cMsqr") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    C_msqr = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  //Bulk correlation scaling factor
  if( strcmp(argv[i], "--clatt") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    N_latt = atof(argv[i+1]);
    i++;
    ret = 0;
    goto out;
  }

  
  //--------------------//
  //     MC params      //
  //--------------------//
  
  //Metro cool-down steps.
  if( strcmp(argv[i], "--nMetroCool") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    n_metro_cool = atoi(argv[i+1]);
    if(n_metro_cool < 0) {
      cout<<"Invalid number of metro cool down steps ("<<n_metro_cool<<") given. Please ensure that n_metro_cool >= 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //Thermalisation sweeps
  if( strcmp(argv[i], "--nTherm") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    n_therm = atoi(argv[i+1]);
    if(n_therm < 0) {
      cout<<"Invalid number of thermalisation steps ("<<n_therm<<") given. Please ensure that n_therm >= 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //HMC steps
  if( strcmp(argv[i], "--nStep") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    n_step = atoi(argv[i+1]);
    if(n_therm < 0) {
      cout<<"Invalid number of HMC steps ("<<n_step<<") given. Please ensure that nStep > 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  
  //Measurements
  if( strcmp(argv[i], "--nMeas") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    n_meas = atoi(argv[i+1]);
    if(n_meas < 0) {
      cout<<"Invalid number of measurements ("<<n_meas<<") given. Please ensure that n_meas >= 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //Measurements
  if( strcmp(argv[i], "--nWrite") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    n_write = atoi(argv[i+1]);
    if(n_write <= 0) {
      cout<<"Invalid number of measurements ("<<n_write<<") given. Please ensure that n_write > 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

   //Skip steps
  if( strcmp(argv[i], "--nSkip") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    n_skip = atoi(argv[i+1]);
    if(n_skip < 0) {
      cout<<"Invalid number of skip steps ("<<n_skip<<") given. Please ensure that n_skip >= 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //Delta phi
  if( strcmp(argv[i], "--tau") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    tau= atof(argv[i+1]);
    if(tau <= 0) {
      cout<<"Invalid tau ("<<tau<<") given. Please ensure that tau > 0.0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  //HMC steps
  if( strcmp(argv[i], "--seed") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    seed = atoi(argv[i+1]);
    if(abs(seed) < 1000) {
      cout<<"Invalid seed ("<<seed<<") given. Please ensure that seed has 4 digits and is > 1000."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

  
  //Cluster steps
  if( strcmp(argv[i], "--nCluster") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    n_cluster = atoi(argv[i+1]);
    if(n_cluster < 0) {
      cout<<"Invalid number of cluster steps ("<<n_cluster<<") given. Please ensure that n_cluster >= 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }

    //Jackknife blocks
  if( strcmp(argv[i], "--nJkBlock") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    n_jkblock = atoi(argv[i+1]);
    if(n_jkblock < 0) {
      cout<<"Invalid number of jackknife blocks ("<<n_jkblock<<") given. Please ensure that n_jkblock > 0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }
  
  //Delta phi
  if( strcmp(argv[i], "--deltaPhi") == 0){
    if (i+1 >= argc){
      usage(argv);
    }
    delta_phi = atof(argv[i+1]);
    if(t_scale == 0) {
      cout<<"Invalid t_scale ("<<t_scale<<") given. Please ensure that tscale != 0.0."<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  }


  //Run-time Verbosity
  if( strcmp(argv[i], "--verbosity") == 0){
    if (i+1 >= argc){
      usage(argv);
    }  
    std::string verbose(argv[i+1]);
    if (verbose == "V" || verbose == "v") {
      verbosity = true;
    } else if(verbose == "Q" || verbose == "q") {
      verbosity = false;
    } else {
      cout<<"Invalid Verbosity condition ("<<verbose<<") given. Use v/q for verbose/quiet"<<endl;
      exit(0);
    }
    i++;
    ret = 0;
    goto out;
  } 

  
 out:
  *idx = i;
  return ret ;
  
}

void param::print() {

  bool useOMP = false;

#ifdef USE_OMP
  useOMP = true;
#endif
  
  cout<<endl;
  cout<<"********************************"<<endl;
  cout<<"*      parameter status        *"<<endl;
  cout<<"********************************"<<endl;
  cout<<endl;
  cout<<"Theory type = "<<(theory_type == PHI4 ? "Phi Fourth" : "Ising")<<endl;
  cout<<"Coupling type = "<<(coupling_type == SR ? "Short-Range" :
			     coupling_type == POW ? "Power Law" : "Radial")<<endl;
  if(theory_type == PHI4) {
    cout<<"(Phi4) Mass squared (M^2) = "<<msqr<<endl;  
    cout<<"(Phi4) coupling (Lambda) = "<<lambda<<endl;
  } else {
    cout<<"(Ising) External field strength (h) = "<<h<<endl;
    cout<<"(Ising) Constant coupling (J) = "<<J<<endl;
  }   
  if(coupling_type != SR) cout<<"Long-Range coupling sigma = "<<sigma<<endl;
  cout<<"Temporal Scaling = "<<t_scale<<endl;
  
  cout<<endl<<"* MC params *"<<endl;
  cout<<"Metropolis Cool-Down = "<<n_metro_cool<<endl;
  cout<<"Thermalisations = "<<n_therm<<endl;
  cout<<"Skip steps = "<<n_skip<<endl;
  cout<<"Measurements = "<<n_meas<<endl;
  cout<<"Cluster steps per Metropolis = "<<n_cluster<<endl;
  
  cout<<endl<<"* Lattice params *"<<endl;
  cout<<"Triangulation = "<<q<<endl;
  cout<<"Levels = "<<Levels<<endl;
  cout<<"2D AdS Boundary points = "<<S1<<endl;
  cout<<"2D AdS Slice points = "<<AdSVol<<endl;
  cout<<"TimeSlices = "<<Lt<<endl;
  cout<<"Centering = "<<(Vcentre == true ? "Vertex" : "Circum")<<" centred"<<endl;

  cout<<endl;
}


//Using the formula c(n) = (q-4)*c(n-1) - c(n-2) where c is the
//number of nodes on circumference at level n, we can construct
//the address of the end node on a given level for triangulation q:
//EN(lev,q) = SUM c(n) n=0..level
long unsigned int endNode(int lev, param &P) { 
  
  int q = P.q;
  
  //This is the case where a vertex is at the centre of the disk
  if(P.Vcentre == true) {
    //Explicit results for level <= 2.
    if(lev==0) return 0;
    if(lev==1) return q;
    if(lev==2) return q*(q-4) + q;
    
    //level >= 3
    int p=q-4; //A convenience  
    long unsigned int c[20]; //Array to hold the circumference info
    
    //initialise array
    for(int i=0; i<20; i++) c[i] = 0;
    //manually set 0,1,2 entries:
    c[0] = 0;
    c[1] = q;
    c[2] = p*q;
    
    long unsigned int EN = 0; //address of node at end of level;  
    //Loop up to lev
    for(int n=3; n<lev+1; n++) {
      c[n] = p*c[n-1] - c[n-2];
      EN += c[n];
    }
    EN += (c[1] + c[2]);
    return EN;
  }

  //This is the case where a circumcentre is at the centre of the disk.
  else {
    //Explicit results for level <= 1.
    if(lev==0) return 2;
    if(lev==1) return (q-3)*3 + 2;
    
    //level >= 2
    int p=q-4; //A convenience  
    long unsigned int c[20]; //Array to hold the circumference info
    
    //initialise array
    for(int i=0; i<20; i++) c[i] = 0;
    //manually set 0,1,2 entries:
    //NB!!! There are three nodes on level 0, but the END NODE is addressed
    //as 2. Therefore, when correcting the number of nodes on a
    //circumference, we use 3 (there are three nodes)
    //but when giving the endNode count, we use 2 (we count from 0)
    c[0] = 3;       
    c[1] = (q-3)*3;
    
    long unsigned int EN = 0; //address of node at end of level;  
    //Loop up to lev
    for(int n=2; n<lev+1; n++) {
      c[n] = p*c[n-1] - c[n-2];
      EN += c[n];
    }
    EN += (2 + c[1]); //0,1,2 nodes to add from circumference 0
    return EN;
  }
}


//Print the hyperbolic and poincare radii
void checkRadius(vector<Vertex> &Lattice, param P){

  double hyp_rad_ave = 0.0;
  double poi_rad_ave = 0.0;
  double hyp_rad_sig = 0.0;
  double poi_rad_sig = 0.0;

  int q = P.q;
  int count = 0;
  int TotNumber = (endNode(P.Levels,P)+1);
  
  for(long unsigned int n=0; n<TotNumber; n++) {
    //loop over neighbours. Break if an outer node is found.
    for(int m=0; m<q; m++) {
      if(Lattice[n].nn[m] == -1 && Lattice[n].pos != -1) {
	//cout<<"n = "<<n<<" |z| = "<<abs(Lattice[n].z)<<" s = ";
	//cout<<s(Lattice[n].z)<<endl;
	
	poi_rad_ave += abs(Lattice[n].z);
	hyp_rad_ave += s(Lattice[n].z);
	count++;
	m=q;
      }
    }
  }
  
  hyp_rad_ave /= count;
  poi_rad_ave /= count;
  
  for(long unsigned int n=0; n<TotNumber; n++) {
    for(int m=0; m<q; m++) {
      if(Lattice[n].nn[m] == -1 && Lattice[n].pos != -1) {
	hyp_rad_sig += pow(hyp_rad_ave - s(Lattice[n].z),2);
	poi_rad_sig += pow(poi_rad_ave - abs(Lattice[n].z),2);
	m=q;
      }
    }
  }

  hyp_rad_sig /= (count-1);
  poi_rad_sig /= (count-1);
  hyp_rad_sig = sqrt(hyp_rad_sig);
  poi_rad_sig = sqrt(poi_rad_sig);
  
  cout<<"HYP RAD AVE = "<<hyp_rad_ave<<endl;
  cout<<"HYP RAD SIG = "<<hyp_rad_sig<<endl;
  cout<<"POI RAD AVE = "<<poi_rad_ave<<endl;
  cout<<"POI RAD SIG = "<<poi_rad_sig<<endl;

  
  //Eyeball the output. Something out of place will
  //stick out like a sore thumb.
  //if(P.verbosity) PrintNodeTables(Lattice, P);  
}


void checkArea(const vector<Vertex> Lattice, param P) {

  double length_01 = 0.0;
  double length_02 = 0.0;
  double length_12 = 0.0;
  double equi_area = area3q(P.q);
  double ave       = 0.0;

  double sig1 = 0.0;
  double sig2 = 0.0;
  int count = 0;

  if(P.verbosity) cout<<endl<<"Checking boundary areas"<<endl;
  for(long unsigned int n=endNode(P.Levels-2,P) + 1; n<endNode(P.Levels-1,P)+1; n++) {
    for(int k=0; k<Lattice[n].fwdLinks; k++) {
      length_01 = d12(Lattice[n].z, Lattice[Lattice[n].nn[k+1]].z);
      length_02 = d12(Lattice[n].z, Lattice[Lattice[n].nn[k+2]].z);
      length_12 = d12(Lattice[Lattice[n].nn[k+1]].z, Lattice[Lattice[n].nn[k+2]].z);      
      ave += areaGeneral(P, length_01, length_02, length_12);
      count++;
    }
  }
  
  ave /= count;
  
  for(long unsigned int n=endNode(P.Levels-2,P) + 1; n<endNode(P.Levels-1,P)+1; n++) {
    for(int k=0; k<Lattice[n].fwdLinks; k++) {
      length_01 = d12(Lattice[n].z, Lattice[Lattice[n].nn[k+1]].z);
      length_02 = d12(Lattice[n].z, Lattice[Lattice[n].nn[k+2]].z);
      length_12 = d12(Lattice[Lattice[n].nn[k+1]].z, Lattice[Lattice[n].nn[k+2]].z);      
      if(P.verbosity) cout<<"n="<<n<<" area "<<k+1<<" = "<<areaGeneral(P, length_01, length_02, length_12)<<endl;
      sig1 += pow(equi_area - areaGeneral(P, length_01, length_02, length_12),2);
      sig2 += pow(ave       - areaGeneral(P, length_01, length_02, length_12),2);
    }
  }

  sig1 /= count - 1;
  sig2 /= count - 1;
  sig1 = sqrt(sig1);
  sig2 = sqrt(sig2);

  cout<<"Boundary areas"<<endl;
  cout<<"AREA EQUI = "<<equi_area<<endl;
  cout<<"AREA STD DEV W.R.T. EQUI = "<<sig1<<endl;
  cout<<"AREA AVE = "<<ave<<endl;  
  cout<<"AREA STD DEV W.R.T AVE = "<<sig2<<endl;  

}

void checkEdgeLength(const vector<Vertex> Lattice, param P) {
  
  int q = P.q;
  int Levels = P.Levels;
  double length = 0.0;
  double sig = 0.0;
  int  nn_node;
  bool Vcentre = P.Vcentre;
  double length_0 = d12(Lattice[0].z, Lattice[1].z);
  double tol = 1e-2;

  //Level 0 is specific to how the graph is centred.
  if(Vcentre) {
    if(P.verbosity) cout<<" lev =  " << 0 << endl;
    if(P.verbosity) cout<<endl<<" Node number = "<<0<<" : "<<endl;
    for(int i = 0; i < q; i++){
      nn_node = Lattice[0].nn[i];
      length = d12(Lattice[0].z, Lattice[nn_node].z);
      if(P.verbosity) cout<<" "<<Lattice[0].nn[i]<<" > "<<length<<"  ";
    }
  }
  else {
    if(P.verbosity) cout<<" lev = "<<0<<endl;
    if(P.verbosity) cout<<endl<<" Node number = "<<0<<" : "<<endl;
    for(int i = 0; i < 2; i++){
      nn_node = Lattice[0].nn[i];
      length = d12(Lattice[0].z, Lattice[nn_node].z);
      if(P.verbosity) cout << Lattice[0].nn[i] << " >  " << length<< "  ";
    }
  }
  
  for(int lev = 1; lev < Levels+1; lev++)  {
    if(P.verbosity) cout<<endl<<endl<<" lev = "<<lev<<endl;      
    for(long unsigned int n = endNode(lev-1,P) + 1;n < endNode(lev,P) + 1 ;n++) {
      if(P.verbosity) cout<<endl<<" Node number = "<<n<<":"<<endl;
      sig += pow( length_0 - d12(Lattice[n].z, Lattice[Lattice[n].nn[q-1]].z), 2);
      
      for(int i = 0; i <q; i++){
	nn_node = Lattice[n].nn[i];
	if(Lattice[n].nn[i] != -1 ) {
	  length = d12(Lattice[n].z, Lattice[nn_node].z);
	  if(P.verbosity) {
	    cout<<" to "<<Lattice[n].nn[i]<<" = "<<length<<" ";
	    if(abs(length - length_0)/length_0 > tol) cout<<"<-! "<<endl;
	    else cout<<"    "<<endl;
	  }
	}
      }
    }
  }
  sig /= endNode(Levels,P);
  sig = sqrt(sig);
  cout<<endl<<"LENGTH STD DEV = "<<sig<<endl;
  if(sig>tol) {
    cout<<"WARNING: Hypergeometric length STD_DEV has diverged over "<<tol<<endl;
    //exit(0);
  }
}

