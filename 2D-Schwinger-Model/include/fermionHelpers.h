#ifndef FERMHELPERS_H
#define FERMHELPERS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>

using namespace std;

//Wilson Fermion Utilities
//---------------------------------------------------------------------------------
// Zero fermion field
template<typename T> inline void zeroField(T psi[LX][LY][2]) {
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int s=0; s<2; s++)
	psi[x][y][s] = 0.0;
}

// Copy fermion field
template<typename T> inline void copyField(T psi2[LX][LY][2],T psi1[LX][LY][2]) {
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int s=0; s<2; s++)
	psi2[x][y][s] = psi1[x][y][s];
}

// Inner product
template<typename T> inline T dotField(T psi1[LX][LY][2], T psi2[LX][LY][2]) {
  T scalar = (T) 0.0;
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int s=0; s<2; s++)
	scalar += conj(psi1[x][y][s])*psi2[x][y][s];
  return scalar;
}

// Norm squared 
template<typename T> inline double norm2(T psi[LX][LY][2]) {  
  double norm2 = 0.0;
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int s=0; s<2; s++)
	norm2 += (psi[x][y][s].real() * psi[x][y][s].real() + psi[x][y][s].imag() * psi[x][y][s].imag());
  
  return norm2;
}


template<typename T> inline void caxpby(T a, T X[LX][LY][2],
					T b, T Y[LX][LY][2],
					T result[LX][LY][2]){  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int s=0; s<2; s++)
	result[x][y][s] = a*X[x][y][s] + b*Y[x][y][s];
}

template<typename T> inline void axpby(double a, T X[LX][LY][2],
				       double b, T Y[LX][LY][2],
				            T result[LX][LY][2]){
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int s=0; s<2; s++)
	result[x][y][s] = a*X[x][y][s] + b*Y[x][y][s];
}

//caxpy in place 
template<typename T> inline void caxpy(T a, T X[LX][LY][2], T Y[LX][LY][2]){
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int s=0; s<2; s++)
	Y[x][y][s] += a*X[x][y][s];
}

//caxpy in result
template<typename T> inline void caxpy(T a, T X[LX][LY][2], T Y[LX][LY][2],
				       T result[LX][LY][2]){  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int s=0; s<2; s++)
	result[x][y][s] = a*X[x][y][s] + Y[x][y][s];
}

//axpy in place 
template<typename T> inline void axpy(double a, T X[LX][LY][2], T Y[LX][LY][2]){
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int s=0; s<2; s++)
	Y[x][y][s] = a*X[x][y][s] + Y[x][y][s];
}

//axpy in result
template<typename T> inline void axpy(double a, T X[LX][LY][2], T Y[LX][LY][2],
				      T result[LX][LY][2]){  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int s=0; s<2; s++)
	result[x][y][s] = a*X[x][y][s] + Y[x][y][s];
}

template<typename T> inline void xpaypbz(T X[LX][LY][2],
					 double a, T Y[LX][LY][2],
					 double b, T Z[LX][LY][2]) {  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int s=0; s<2; s++) {
	Z[x][y][s] *= b;
	Z[x][y][s] += X[x][y][s] + a*Y[x][y][s];
      }
}

template<typename T> inline void ax(double a, T X[LX][LY][2]){ 
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int s=0; s<2; s++)
	X[x][y][s] *= a;
}

template<typename T> inline void printVector(T X[LX][LY][2]){
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      cout << "(" << x << "," << y << ") = " << X[x][y][0] << " " << X[x][y][1]<<endl; 
}
    
//Staggered Fermion Utilities
//---------------------------------------------------------------------------------
// Zero fermion field
template<typename T> inline void zeroField(T psi[LX][LY]) {
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      psi[x][y] = 0.0;
}

// Copy fermion field
template<typename T> inline void copyField(T psi2[LX][LY],T psi1[LX][LY]) {
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      psi2[x][y] = psi1[x][y];
}

// Inner product 
template<typename T> inline T dotField(T psi1[LX][LY], T psi2[LX][LY]) {
  T scalar = (T) 0.0;
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      scalar += conj(psi1[x][y])*psi2[x][y];
  return scalar;
}

// Norm squared 
template<typename T> inline double norm2(T psi[LX][LY]) {
  
  double norm2 = 0.0;
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      norm2 += (psi[x][y].real() * psi[x][y].real() + psi[x][y].imag() * psi[x][y].imag());
  
  return norm2;
}


template<typename T> inline void caxpby(T a, T X[LX][LY],
					T b, T Y[LX][LY],
					T result[LX][LY]){  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      result[x][y] = a*X[x][y] + b*Y[x][y];
}

template<typename T> inline void axpby(double a, T X[LX][LY], double b,
				       T Y[LX][LY], T result[LX][LY]){
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      result[x][y] = a*X[x][y] + b*Y[x][y];
}

//Staggered fermion axpy in place 
template<typename T> inline void axpy(double a, T X[LX][LY], T Y[LX][LY]){ 
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      Y[x][y] = a*X[x][y] + Y[x][y];
}

//Staggered fermion axpy in result
template<typename T> inline void axpy(double a, T X[LX][LY], T Y[LX][LY],
				      T result[LX][LY]){ 
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      result[x][y] = a*X[x][y] + Y[x][y];
}

template<typename T> inline void xpaypbz(T X[LX][LY],
					 double a, T Y[LX][LY],
					 double b, T Z[LX][LY]) {
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++) {
      Z[x][y] *= b;
      Z[x][y] += X[x][y] + a*Y[x][y];
    }
}

template<typename T> inline void ax(double a, T X[LX][LY]){
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      X[x][y] *= a;
}

#endif
