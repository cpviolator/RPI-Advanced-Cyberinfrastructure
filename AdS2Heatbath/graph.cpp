#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>

#include "util.h"
#include "hyp_util.h"
//#include "data_io.h"
#include "graph.h"

using namespace std;

//- Get the z coordinates of every node on the Poincare disk 
void getComplexPositions(std::vector<Vertex> &Lattice, param& p){

  int q = p.q;
  int Levels = p.Levels;
  int T_offset = endNode(p.Levels,p)+1;
  
  if(p.Vcentre == true) {
    //Assume for now that the origin (level 0) is a vertex
    Lattice[0].z = 0.0;

    //Assert that node 1 is on the real axis
    complex<double> init(edgeLength(q),0.0);
    Lattice[1].z = init;
    //Rotate to create level level 1
    for(int k=1; k<q+1; k++) {
      Lattice[k].z = newVertex(init, 0.0, k-1, q);
    }
    //For every node on level >=1, the zeroth
    //nn is the (n-1)th link on the same level. If
    //n=1, the node address is the endnode value.
    for(int l=1; l<Levels+1; l++) {
      for(long unsigned int n=endNode(l-1,p)+1; n<endNode(l,p)+1; n++) {
	for(int k=0; k<q; k++) {
	  if(Lattice[n].nn[k] != -1) {
	    Lattice[Lattice[n].nn[k]].z = newVertex(Lattice[Lattice[n].nn[0]].z, Lattice[n].z, k, q);
	  }
	}
      }
    }    
  }
  
  else {
    
    double numer = sqrt(cos(M_PI*(q+6)/(6*q)) - sin(M_PI/q));
    double denom = sqrt(sin(M_PI/q) + sin(M_PI*(q+3)/(3*q)));    
    double init_mod = sqrt(norm(numer/denom));
    
    //Assume that node 0 lies on +ve real axis
    complex<double> init_0(init_mod,0.0);
    Lattice[0].z = init_0;
    
    //Assert that node 1 is node 0 rotated by 2*PI/3
    complex<double>init_1(init_mod*cos(2.0*M_PI/3.0),
			 init_mod*sin(2.0*M_PI/3.0));
    
    Lattice[1].z = init_1;
    //Rotate node 1 about node 0 to create level 0 (the equilateral triangle)
    Lattice[2].z = newVertex(init_1, init_0, 1, q);

    //For every node on level >=1, the zeroth
    //nn is the (n-1)th link on the same level. If
    //n=1, the node address is the endnode value.
    for(long unsigned int n=0; n<endNode(1,p)+1; n++) {
      for(int k=0; k<q; k++) {
	if(Lattice[n].nn[k] != -1) {
	  Lattice[Lattice[n].nn[k]].z = newVertex(Lattice[Lattice[n].nn[0]].z, Lattice[n].z, k, q);
	}
      }
    }
    for(int l=1; l<Levels+1; l++) {
      for(long unsigned int n=endNode(l-1,p)+1; n<endNode(l,p)+1; n++) {
	for(int k=0; k<q; k++) {
	  if(Lattice[n].nn[k] != -1) {
	    Lattice[Lattice[n].nn[k]].z = newVertex(Lattice[Lattice[n].nn[0]].z, Lattice[n].z, k, q);
	  }
	}
      }
    }
  }

  if(p.Lt > 1) {
    //Copy all 2D complex positions and weights along the cylinder
    for(long unsigned int n=0; n<endNode(p.Levels,p)+1; n++) 
      for(int t=1; t<p.Lt; t++) {
	Lattice[n + T_offset*t].z = Lattice[n].z;
      }
  }
}


//- For each node n, with a link to another node,
//  it checks that the neighbour table on the linked
//  node contains the original node n as a neighbour.
void connectivityCheck(vector<Vertex> &Lattice, param p){

  int q = p.q;
  int Levels = p.Levels;
  int T = p.Lt;
  int TotNumber = T*(endNode(Levels,p)+1);
  int t_offset  = 0;
  T == 1 ? t_offset = 0 : t_offset = 2;
  
  //Object to hold boolean values of graph connectivity.
  vector<Vertex> AuxLattice(TotNumber);
  //Initialise to 0.
  for(int n = 0; n <TotNumber;n++)
    for(int mu = 0; mu < q+t_offset; mu++) {
      AuxLattice[n].nn[mu] = 0;
    }
  
  for(long unsigned int n=0; n<TotNumber; n++) {
    //Check that the node is valid
    if(Lattice[n].pos != -1) {      
      for(int m=0; m<q+t_offset; m++) {
	//Check that the link is valid
	if(Lattice[n].nn[m] != -1) {
	  for(int p=0; p<q+t_offset; p++) {
	    //Loop over all links on the linked node,
	    //check if original node exists in neighbour
	    //table.
	    if( n == Lattice[ Lattice[n].nn[m] ].nn[p] ) {
	      AuxLattice[n].nn[m] = 1;
	    }
	  }
	}
      }
    }
  }

  //Eyeball the output. something out of place will
  //stick out like a sore thumb.
  //PrintNodeTables(AuxLattice, p);
}

//Truncate the graph according to the hyperbolic radius
//condition |z| < s.
void hypRadGraph(vector<Vertex> &Lattice, param &p){

  int q = p.q;
  int Levels = p.Levels;
  int T = p.Lt;
  int TotNumber = T*(endNode(Levels,p)+1);
  int t_offset  = 0;
  T == 1 ? t_offset = 0 : t_offset = 2;

  //Find radius to maximize connectivity.
  int r_min_pos = 0;
  int r_max_pos = 0;
  double r_min = 1.0;
  double r_max = 0.0;

  //Locate the node on the outer circumference with the smallest
  //radius
  for(int n=endNode(p.Levels-1,p)+1; n<endNode(p.Levels,p)+1; n++){
    if(abs(Lattice[n].z) < r_min) {
      r_min = abs(Lattice[n].z);
      r_min_pos = n;
    }
  }
  p.hyp_rad = s(Lattice[r_min_pos].z);
  cout<<"HYP_RAD = "<<p.hyp_rad<<endl;

  double hyp_rad = p.hyp_rad;
  
  for(long unsigned int n=0; n<TotNumber; n++) {
    if(s(Lattice[n].z) >= hyp_rad + 0.0) {     
      //This node must be removed. loop over its neighbours
      //and remove this specific connection.
      if(p.verbosity) cout<<"Deletng node: "<<n<<" connections: ";
      Lattice[n].pos = -1;
      for(int m=0; m<q+t_offset; m++)
	for(int k=0; k<q+t_offset; k++) {
	  if(Lattice[Lattice[n].nn[m]].nn[k] == n) {
	    Lattice[Lattice[n].nn[m]].nn[k] = -1;
	    Lattice[n].nn[m] = -1;
	    if(p.verbosity) cout<<Lattice[n].nn[m]<<" ";
	  }
	}
      if(p.verbosity) cout<<endl;
    }
  }
  
  for(int n=endNode(p.Levels-2,p)+1; n<endNode(p.Levels-1,p)+1; n++){
    if(abs(Lattice[n].z) < r_min && Lattice[n].pos != -1) {
      r_min = abs(Lattice[n].z);
      r_min_pos = n;
    }
    if(abs(Lattice[n].z) > r_max && Lattice[n].pos != -1) {
      r_max = abs(Lattice[n].z);
      r_max_pos = n;
    }
  }
  p.r_max_pos = r_max_pos;
  cout<<"R MAX POS = "<<p.r_max_pos<<endl;
  p.r_min_pos = r_min_pos;
  cout<<"R MIN POS = "<<p.r_min_pos<<endl;
  
  //Eyeball the output. Something out of place will
  //stick out like a sore thumb.
  //if(p.verbosity == "d") radiusCheck(Lattice, P); 
}


void buildGraph(vector<Vertex> &Lattice, param p) {
  
  int q = p.q;
  int Levels = p.Levels;
  int offset = endNode(Levels,p) + 1;
  
  if(p.Vcentre == true) {

    //Level 0 spatial: trivial
    Lattice[0].pos = 0;
    Lattice[0].fwdLinks = q;
    for(int mu=1; mu<q+1; mu++) {
      Lattice[0].nn[mu-1] = mu;
    }
    
    // CONVENTION: The link on the same level, going from node n to (n-1)
    //             is the 0th entry in the neighbour table. The next entries
    //             are the links going to the next higher level. The next
    //             entry is same level link going from n to (n+1). The last
    //             entries go to the lower level. This means all links in the 
    //             neighbour table are listed in anti-clockwise order.
    
    //Level 1
    for(long unsigned int n=endNode(0,p)+1; n<endNode(1,p)+1; n++){

      Lattice[n].pos = n;
      
      //This is the first node, treat it separately.
      n-1 == 0 ? Lattice[n].nn[0] = endNode(1,p) : Lattice[n].nn[0] = n-1;
      
      //Address of first new node on level l+1 from node a
      int x = endNode(1,p)+1 +(n-1)*(q-4); 
      
      //get new nodes on level 2
      for(int i=1; i<q-2; i++) {
	Lattice[n].nn[i] = x+i-1;
	Lattice[x+i-1].nn[q-1] = n;
	Lattice[x+i-1].fwdLinks = q-3;
	//By definition, the first link in this loop has two back links.
	if(i==1) {
	  Lattice[x+i-1].nn[q-2] = n;
	  Lattice[x+i-1].fwdLinks = q-4;
	  n == 1 ? Lattice[x+i-1].nn[q-1] = q : Lattice[x+i-1].nn[q-1] = n-1;
	}
      }
      Lattice[n].nn[q-2] = n%q+1;
      Lattice[n].nn[q-1] = 0;
    }
    //Fix (q-3) link on final node 
    Lattice[q].nn[q-3] = endNode(1,p)+1;
    
    //Level >=2
    for(int l=2; l<Levels+1; l++){
      
      //Get first new node on level l+1
      int x = endNode(l,p)+1;
      
      //Loop over all nodes on this level
      for(long unsigned int n=endNode(l-1,p)+1; n<endNode(l,p)+1; n++){      

	Lattice[n].pos = n;
	
	//Assign links on the same level 
	//Check if first node
	if(n == endNode(l-1,p)+1) {
	  Lattice[n].nn[0] = endNode(l,p);
	} else {
	  Lattice[n].nn[0] = n-1;
	}
	//Check if last node
	if(n == endNode(l,p)) {
	  Lattice[n].nn[Lattice[n].fwdLinks]   = n+1;
	  Lattice[n].nn[Lattice[n].fwdLinks+1] = endNode(l-1,p)+1;
	}
	else Lattice[n].nn[Lattice[n].fwdLinks+1] = n+1;
	
	//Loop over all links on node n to level l+1. If the node
	//has two back links, there are only (q-4) new links
	//to l+1, else there are (q-3).
	//By deiniftion, the first link in this loop has two back links.
	if(l<Levels) {
	  for(int i=1; i<Lattice[n].fwdLinks+1; i++) {
	    Lattice[n].nn[i] = x+i-1;
	    Lattice[x+i-1].nn[q-1] = n;
	    Lattice[x+i-1].fwdLinks = q-3;
	    if(i==1) {
	      Lattice[x+i-1].nn[q-2] = n;
	      Lattice[x+i-1].fwdLinks = q-4;
	      n == endNode(l-1,p)+1 ? Lattice[x+i-1].nn[q-1] = endNode(l-1,p) : Lattice[x+i-1].nn[q-1] = n-1;
	    }
	  }
	}
	x += Lattice[n].fwdLinks-1;

	//If this is a boundary node, correct the fwdLinks to be zero.
	if(l == Levels) Lattice[n].fwdLinks = 0;	
	//fix link q-1 on start node
	Lattice[endNode(l-1,p)+1].nn[q-1]=endNode(l-1,p);
	//fix link q-2 on start node
	Lattice[endNode(l-1,p)+1].nn[q-2]=endNode(l-2,p)+1;
	//fix link q-3 on end node
	if(n == endNode(Levels,p)) Lattice[endNode(l,p)].nn[q-3] = -1;
	else Lattice[endNode(l,p)].nn[q-3] = endNode(l,p) + 1;
      }
    }

    /*
    //Populate temporal links on t=0 disk
    for(long unsigned int n=0; n<endNode(Levels,p)+1; n++) {
      //Fwd link
      Lattice[n].nn[q  ] = n + offset;
      //Bkd link
      Lattice[n].nn[q+1] = (T-1)*offset + n;
    }
    
    //Construct disks and t links for 0 < t < T
    for(int t=1; t<T; t++)
      for(long unsigned int n=0; n<endNode(Levels,p)+1; n++) {
	for(int i=0; i<q; i++) {
	  Lattice[(t-1)*offset + n].nn[i] == -1 ?
	    Lattice[t*offset + n].nn[i] = -1 : 
	    Lattice[t*offset + n].nn[i] =
	    Lattice[(t-1)*offset + n].nn[i] + offset;
	}
	//Fwd link
	Lattice[t*offset + n].nn[q  ] = (t+1)*offset + n;    
	//Bkd link
	Lattice[t*offset + n].nn[q+1] = (t-1)*offset + n;
	//fwdLinks data
	Lattice[t*offset + n].fwdLinks = Lattice[n].fwdLinks;
	//pos data
	Lattice[t*offset + n].pos = t*offset + n;
      }
    
    
    //Correct forward t links for t = T-1
    int t=T-1;
    for(long unsigned int n=0; n<endNode(Levels,p)+1; n++) {
      Lattice[t*offset + n].nn[q] = n;
    }
    */
  }
  else {  
    //Level 0
    for(long unsigned int n=0; n<3; n++){

      Lattice[n].pos = n;
      
      //This is the first node, treat it separately.
      n == 0 ? Lattice[n].nn[0] = 2 : Lattice[n].nn[0] =  n-1;
      
      //Address of first new node on level 1 from node n
      int x = 3 + n*(q-3);
      
      //get new nodes on level 1
      for(int i=1; i<q-1; i++) {
	Lattice[n].nn[i] = x+i-1;
	Lattice[n].fwdLinks = q-2;
	
	Lattice[x+i-1].nn[q-1] = n;
	Lattice[x+i-1].fwdLinks = q-3;
	//Corrections
	if(i==1) {
	  //This node (3,7,11) has two back links to level 0:
	  //3:  0,2
	  //7:  1,0
	  //11: 2,1
	  Lattice[x+i-1].fwdLinks = q-4;
	  Lattice[x+i-1].nn[q-2] = (x+i-1)%3;
	  Lattice[x+i-1].nn[q-1] = ((x+i-1)%3 + 2)%3;	
	}
	n == 2 ? Lattice[n].nn[q-1] = 0 : Lattice[n].nn[q-1] =  n+1;
      }
    }
    //Fix (q-2) link on final node 
    Lattice[2].nn[q-2] = 3;
    
    //Level 1
    
    //Get first new node on level 2
    int x = endNode(1,p)+1;
    
    //Loop over all nodes on level 1.
    for(long unsigned int n=endNode(0,p)+1; n<endNode(1,p)+1; n++){      

      Lattice[n].pos = n;
      
      //Assign links on the same level 
      //Check if first node
      if(n == endNode(0,p)+1) {
	Lattice[n].nn[0] = endNode(1,p);
      } else {
	Lattice[n].nn[0] = n-1;
      } 
      //Check if last node
      if(n == endNode(1,p)) {
	Lattice[n].nn[Lattice[n].fwdLinks]   = n+1;
	Lattice[n].nn[Lattice[n].fwdLinks+1] = endNode(0,p)+1;
      }
      else Lattice[n].nn[Lattice[n].fwdLinks+1] = n+1;
      
      //Loop over new links
      for(int i=1; i<Lattice[n].fwdLinks+1; i++) {
	Lattice[n].nn[i] = x+i-1;
	Lattice[x+i-1].nn[q-1] = n;
	Lattice[x+i-1].fwdLinks = q-3;
	if(i==1) {
	  Lattice[x+i-1].nn[q-2] = n;
	  Lattice[x+i-1].fwdLinks = q-4;
	  n == endNode(0,p)+1 ? Lattice[x+i-1].nn[q-1] = endNode(1,p) : Lattice[x+i-1].nn[q-1] = n-1;
	}
      }
      x += Lattice[n].fwdLinks-1;
    }
    //Fix (q-3) link on final node 
    Lattice[endNode(1,p)].nn[q-3] = endNode(1,p)+1;
    
    //Level >=2
    for(int l=2; l<Levels+1; l++){
      
      //Get first new node on level l+1
      int x = endNode(l,p)+1;    
      //Loop over all nodes on this level
      for(long unsigned int n=endNode(l-1,p)+1; n<endNode(l,p)+1; n++){      
	
	//Assign links on the same level 
	//Check if first node
	if(n == endNode(l-1,p)+1) {
	  Lattice[n].nn[0] = endNode(l,p);
	} else {
	  Lattice[n].nn[0] = n-1;
	}
	//Check if last node
	if(n == endNode(l,p)) {
	  Lattice[n].nn[Lattice[n].fwdLinks]   = n+1;
	  Lattice[n].nn[Lattice[n].fwdLinks+1] = endNode(l-1,p)+1;
	}
	else Lattice[n].nn[Lattice[n].fwdLinks+1] = n+1;
	
	//Loop over all links on node n to level l+1. If the node
	//has two back links, there are only (q-4) new links
	//to l+1, else there are (q-3).
	//By deiniftion, the first link in this loop has two back links.
	if(l<Levels) {
	  for(int i=1; i<Lattice[n].fwdLinks+1; i++) {
	    Lattice[n].nn[i] = x+i-1;
	    Lattice[x+i-1].nn[q-1] = n;
	    Lattice[x+i-1].fwdLinks = q-3;
	    if(i==1) {
	      Lattice[x+i-1].nn[q-2] = n;
	      Lattice[x+i-1].fwdLinks = q-4;
	      n == endNode(l-1,p)+1 ? Lattice[x+i-1].nn[q-1] = endNode(l-1,p) : Lattice[x+i-1].nn[q-1] = n-1;
	    }
	  }
	}
	x += Lattice[n].fwdLinks-1;
	
	//fix link q-1 on start node
	Lattice[endNode(l-1,p)+1].nn[q-1]=endNode(l-1,p);
	//fix link q-2 on start node
	Lattice[endNode(l-1,p)+1].nn[q-2]=endNode(l-2,p)+1;
	//fix link q-3 on end node
	if(n == endNode(Levels,p)) Lattice[endNode(l,p)].nn[q-3] = -1;
	else Lattice[endNode(l,p)].nn[q-3] = endNode(l,p) + 1;
      }
    }

    /*
    //Populate temporal links on t=0 disk
    for(long unsigned int n=0; n<endNode(Levels,p)+1; n++) {
      //Fwd link
      Lattice[n].nn[q  ] = n + offset;
      //Bkd link
      Lattice[n].nn[q+1] = (T-1)*offset + n;
    }
    
    //Construct disks and t links for 0 < t < T
    for(int t=1; t<T; t++)
      for(long unsigned int n=0; n<endNode(Levels,p)+1; n++) {
	for(int i=0; i<q; i++) {
	  Lattice[(t-1)*offset + n].nn[i] == -1 ?
	    Lattice[t*offset + n].nn[i] = -1 : 
	    Lattice[t*offset + n].nn[i] =
	    Lattice[(t-1)*offset + n].nn[i] + offset;
	}
	//Fwd link
	Lattice[t*offset + n].nn[q  ] = (t+1)*offset + n;    
	//Bkd link
	Lattice[t*offset + n].nn[q+1] = (t-1)*offset + n;
	//fwdLinks data
	Lattice[t*offset + n].fwdLinks = Lattice[n].fwdLinks;
	//pos data
	Lattice[t*offset + n].pos = t*offset + n;
      }
    
    //Correct forward t links for t = T-1
    int t=T-1;
    for(long unsigned int n=0; n<endNode(Levels,p)+1; n++) {
      Lattice[t*offset + n].nn[q] = n;
    }
    */
  }
}
