#ifndef LATHELPERS_H
#define LATHELPERS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>

using namespace std;

//Lattice utilities
//---------------------------------------------------------------------------

// Zero lattice 2D
template<typename T> inline void zeroLat(T v[LX][LY][2]) {
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int mu=0; mu<2; mu++)
	v[x][y][mu] = 0.0;
}

// Copy lattice 2D
template<typename T> inline void copyLat(T v2[LX][LY][2], const T v1[LX][LY][2]) {
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int mu=0; mu<2; mu++)
	v2[x][y][mu] = v1[x][y][mu];
}

// Zero Wilson Loop field 2D.
template<typename T> inline void zeroWL(T psi[LX/2][LY/2]) {
  for(int x=0; x<LX/2; x++)
    for(int y=0; y<LY/2; y++)
      psi[x][y] = 0.0;
}

#endif
