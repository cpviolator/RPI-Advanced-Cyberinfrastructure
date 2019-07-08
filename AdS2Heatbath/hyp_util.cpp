#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <unistd.h>

#include "hyp_util.h"
#include "util.h"

using namespace std;

// Translate w to 0 
complex<double> T(complex<double> z,  complex<double> w){ //translate w to 0
  return (z - w)/(z*conj(w) + (double)1.0);
}

// Rotate at z = 0
complex<double> R(complex<double> z, complex<double> omega){
  // rotate by omega = exp [i theta] about z= 0
  return omega*z;
 }

//Reflection z accross the z1 -- z2 line
complex<double> flip(complex<double> z, complex<double> z1, complex<double> z2){
  // reflection (or flip)  z across (z1,z2)
  return T( R( conj( T(z,z1) ), z2/conj(z2) ), -z1 );
}

//Hyperbolic distance s, from origin to z
double s(complex<double> z){
  return log(((double)1.0+abs(z))/((double)1.0-abs(z)));
}

//Poincare distance |z| from origin to s
double r(double s){
  return tanh(s/2);
}

//Geodesic distance from z1 to z2
double d12(complex<double> z, complex<double> w) {
  return log ( (abs((double)1.0-conj(z)*w) + abs(z-w))/(abs((double)1.0-conj(z)*w) - abs(z-w)));  
}

//Geodesic distance from z1,t1 to z2,t2
double sigma(complex<double> z, complex<double> w, double delta_t) {
  
  double theta = atan2( (w/z).imag() , (w/z).real() );
  double r  = abs(z);
  double rp = abs(w);  
  double xi = (cosh(delta_t)*(1+r*r)*(1+rp*rp) - 4*r*rp*cos(theta)) / ((1-r*r)*(1-rp*rp));
  return acosh(xi);
  
}

//Geodesic distance from z1,t1 to z2,t2
double sigmad(double dt, double dth, double r) {
  
  double xi = (cosh(dt)*(1+r*r)*(1+r*r) - 4*r*r*cos(dth))/((1-r*r)*(1-r*r)); 
  return acosh(xi);  
}

//Geodesic distance from z1,t1 to z2,t2
double sigmaL(complex<long double> z,
	      complex<long double> w,
	      double delta_t) {
  
  long double theta = atan2( (w/z).imag() , (w/z).real() );
  long double r = abs(z);
  long double rp = abs(w);  
  long double xi = (cosh(delta_t)*(1.0L+r*r)*(1.0L+rp*rp) - 4*r*rp*cos(theta)) / ((1.0L-r*r)*(1.0L-rp*rp)); 
  return acosh(xi);   
}

// length of arc q fold triangle to origin.
double s3p(int q){  //vertex centered Arc lengeth
  return (double)2.0*acosh((double)1.0/sin(M_PI/(double)q));
}

// Area equilateral triangle with angles 2 pi/q
double area3q(int q){
  //pi - (3 * hyp_angle) = defect
  return M_PI - (double)3.0*(2.0*M_PI/(double)q);
}

// Area non-equilateral triangle with side length a,b,c
double areaGeneral(param P, double a, double b, double c) {
  //pi - (A+B+C) = defect
  
  // use general cosine law:
  // cos(C) = (cosh(c) - cosh(a)cosh(b)) / sinh(a)sinh(b)
  double C = acos( -(cosh(c) - cosh(a)*cosh(b)) / (sinh(a)*sinh(b)) );
  
  // use general sine law:
  // sin(A)/sinh(a) = sin(B)/sinh(b) = ...
  double B = asin( sinh(b)*sin(C)/sinh(c) );
  double A = asin( sinh(a)*sin(C)/sinh(c) );

  return M_PI - (A+B+C);
}

//
double centralRad(double s){
  return (sqrt( cosh(s/2.0) - (double)1.0/4.0) - 0.5*sqrt(3.0))/sqrt(cosh(s/2.0) -(double)1.0);
}

//
complex<double> DisktoUHP(complex<double> z) {
  // z = -1 -i, 1, i maps to u =  -1, 0 1, infty
  return (z + I)/((double)1.0 + I * z);
}

complex<double> UHPtoDisk(complex<double> u) {
  // u = 0, 1, infty  maps to -1 -i , 1, i  
  return (u - I)/((double)1.0 - I*u); 
}

//- Rotate z about z0 by 2*k*pi/q 
complex<double> newVertex(complex<double> z,complex<double> z0, int k, int q) {

  complex<double> w( 0.0, 2.0 * sin(k * M_PI/q) );
  complex<double> a( cos(k*M_PI/q)*((double)1.0 - norm(z0)), sin(k*M_PI/q)*((double)1.0 + norm(z0)) ); 
  w = w*z0;
  
  //cout<<"New z = "<<-(a*z - w)/(conj(w)*z - conj(a))<<endl;
  return - (a*z - w)/(conj(w)*z - conj(a)); 
}

complex<double> inversion(complex<double> z0, double r) {
  // z_image conj(z0) = r^2
  return r*2/conj(z0);
}

complex<double> squareInversion(complex<double>z0,double r1,double r2 ) {
  return inversion(inversion(z0, r1),r2);
}

double greens2D(complex<double> z, complex<double> w) {
  return -log( tanh ( log ( (abs((double)1.0-conj(z)*w) + abs(z-w))/(abs((double)1.0-conj(z)*w) - abs(z-w)) )/2 ) );    
}

double greensM2D(complex<double> z, complex<double> w, param p) {
  
  //First compute 2F1  

  double delta = p.Lt > 1 ? 1.0+sqrt(1.0+p.msqr) : 0.5+sqrt(0.25+p.msqr);
  double h = 1;
  double result = 0.0;
  double result_0 = 0.0;
  double geo = exp(-2*d12(z,w));
  double a,b,c;
  double tol = 1e-10;
  int n=0;
  bool conv = false;

  while( !conv && n < 10000 ) {    
    result_0 = result;
    a = tgamma(delta + n)/tgamma(delta);
    b = tgamma(h + n)/tgamma(h);
    c = tgamma(delta+1-h + n)/tgamma(delta+1-h);
    result += ( (a*b) / (c*tgamma(n+1)) ) * pow(geo,n);
    if( abs(result_0 - result)/result_0 < tol ) conv = true;
    n++;
    if(n%10000 == 0) cout<<n<<" 2F1 iters "<<geo<<" "<<abs(result_0 - result)/result_0<<endl; 
  }
  
  //Then compute coefficient. 
  result *= pow(geo,delta/2) * tgamma(delta) / (2*pow(M_PI,h)*tgamma(delta+1-h));

  return result;
}

double AdS2p1Prop(complex<long double> z1, complex<long double> z2,
		  double delta_t, param p){

  long double Delta = 1.0L + sqrt(1.0L + p.msqr);
  long double sigma_val = sigmaL(z1, z2, delta_t);
  
  return exp(-Delta*sigma_val) / (1.0L - exp(-2.0L*sigma_val));
}
