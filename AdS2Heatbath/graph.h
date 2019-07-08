#ifndef GRAPH_H
#define GRAPH_H
#include <complex>
#include <cstring>
#include <vector>


using namespace std;

//Get the z coordinates of every node on the Poincare disk 
void getComplexPositions(std::vector<Vertex> &Lattice, param& P);

//For each node n, with a link to another node,
//it checks that the neighbour table on the linked
//node contains the original node n as a neighbour.
void connectivityCheck(std::vector<Vertex> &Lattice, param P);

//Truncate the graph according to the hyperbolic radius
//condition |z| < s.
void hypRadGraph(std::vector<Vertex> &Lattice, param &P);

//- Construct the nearest neighbour table
void buildGraph(vector<Vertex> &Lattice, param P);


#endif
