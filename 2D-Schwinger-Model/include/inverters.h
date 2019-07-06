#ifndef INVERTERS_H
#define INVERTERS_H

#include "dOpHelpers.h"

//===============================================================
// CG solutions to Apsi = b 
// see http://en.wikipedia.org/wiki/Conjugate_gradient_method
//===============================================================
//       x: The solution
//       b: The RHS vector
//      x0: An initial guess
//   gauge: The gauge field defining the operator
//   param: The parameter container

// Wilson g5Dg5D matrix inverter
//---------------------------------------------------------------
int Ainvpsi(Complex x[LX][LY][2], Complex b[LX][LY][2], Complex x0[LX][LY][2],
	    const Complex gauge[LX][LY][2], param_t param) {
  
  int success = 0;
  
  Complex res[LX][LY][2], p[LX][LY][2], Ap[LX][LY][2], tmp[LX][LY][2];
  double alpha, beta, denom;
  double rsq = 0, rsqNew = 0, bsqrt = 0.0, bnorm = 0.0;
  bool deflating = false;
  
  //Intialize
  zeroField(res);
  zeroField(Ap);
  zeroField(p);
  zeroField(x);
  
  // Find norm of rhs.
  bnorm = norm2(b);
  bsqrt = sqrt(bnorm);
  if(bsqrt == 0 || bsqrt != bsqrt) {
    printVector(b);
    cout << "Error in Wilson Ainvpsi: inverting on zero source... or nan!" << endl;
    exit(0);
  }
  copyField(res, b);
  
  // res = b - A*x0
  if (norm2(x0) != 0.0) {
    
    //Solve the deflated system.
    deflating = true;
    DdagDpsi(tmp, x0, gauge, param);    
    axpy(-1.0, tmp, res);
    
    cout << "using initial guess, |x0| = " << sqrt(norm2(x0))
	 << ", |b| = " << bsqrt
	 << ", |res| = " << sqrt(norm2(res)) << endl;
  }
  
  copyField(p, res);
  rsq = norm2(res);
  
  // Iterate until convergence
  int k;
  for (k=0; k<param.maxIterCG; k++) {

    // Compute Ap.
    DdagDpsi(Ap, p, gauge, param);
    
    denom = real(dotField(p, Ap));
    alpha = rsq/denom;
    
    axpy( alpha, p,  x);
    axpy(-alpha, Ap, res);
    
    // Exit if new residual is small enough
    rsqNew = norm2(res);
    //printf("CG iter %d, rsq = %g\n", k+1, rsqNew);
    if (rsqNew < param.eps*bnorm) {
      rsq = rsqNew;
      break;
    }
    
    // Update vec using new residual
    beta = rsqNew/rsq;
    rsq = rsqNew;
    
    axpy(beta, p, res, p);
    
  } // End loop over k
  
  if(k == param.maxIterCG) {
    // Failed convergence 
    printf("CG: Failed to converge iter = %d, rsq = %.16e\n", k+1, rsq); 
    success = 0; 
  } else {
    // Convergence 
    success = 1; 
  }
  
  if(deflating) {
    // x contains the solution to the deflated system b - A*x0.
    // We must add back the exact part
    axpy(1.0, x0, x);
    // x now contains the solution to the RHS b.
  }
  DdagDpsi(tmp, x, gauge, param);
  axpy(-1.0, tmp, b, res);
  
  //double truersq = real(dotField(res, res));
  //printf("CG: Converged iter = %d, rsq = %.16e, truersq = %.16e\n", k+1, rsq, truersq/(bsqrt*bsqrt));
  return success;
  
}

// let dD \equiv (d/dtheta D)
//
// d/dtheta (phi^* (DD^dag)^-1 phi) = -((DD^dag)^1 phi)^dag ([dD]*D^dag + D*[dD^dag]) ((DD^dag)^-1 phi)
//
// *****  Should optimize this to operate only on EVEN sites. ****

void forceD(double fD[LX][LY][2], Complex gauge[LX][LY][2], Complex phi[LX][LY][2],
	    Complex guess[LX][LY][2], param_t p){
  
  if(p.dynamic == true) {

    zeroLat(fD);
    
    //phip = (D^dagD)^-1 * phi
    Complex phip[LX][LY][2];
    zeroField(phip);
    
    //Ainvpsi inverts using the DdagD (g3Dg3D) operator, returns
    // phip = (D^-1 * Ddag^-1) phi = (D^-1 * g3 * D^-1 g3) phi.
    Complex guess[LX][LY][2]; //Initial guess to CG
    zeroField(guess);
    Ainvpsi(phip, phi, guess, gauge, p);
    
    //g3Dphi = g3D * phip
    Complex g3Dphi[LX][LY][2];
    zeroField(g3Dphi);
    g3Dpsi(g3Dphi, phip, gauge, p);
    
    int xp1, xm1, yp1, ym1;
    double r = 1.0;
    for(int x=0; x<LX; x++)
      for(int y=0; y<LY; y++) {

	xp1 = (x+1)%LX;
	yp1 = (y+1)%LY;
	xm1 = (x-1+LX)%LX;
	ym1 = (y-1+LY)%LY;	
	
	//mu = 0
	//upper
	// | r  1 | 
	// | 1  r |
	//lower
	// | r -1 |
	// | 1 -r |					
	fD[x][y][0] += real(I*((conj(gauge[x][y][0]) *
			       (conj(phip[xp1][y][0]) * (r*g3Dphi[x][y][0] +   g3Dphi[x][y][1]) -
				conj(phip[xp1][y][1]) * (  g3Dphi[x][y][0] + r*g3Dphi[x][y][1])))
			       -
			       (gauge[x][y][0] *
			       (conj(phip[x][y][0]) * (r*g3Dphi[xp1][y][0] -   g3Dphi[xp1][y][1]) +
				conj(phip[x][y][1]) * (  g3Dphi[xp1][y][0] - r*g3Dphi[xp1][y][1])))
			       )
			    );	
	
	//mu = 1
	//upper
	// | r -i | 
	// | i  r |
	//lower
	// | r  i |
	// | i -r |
	fD[x][y][1] += real(I*((conj(gauge[x][y][1]) *
				(conj(phip[x][yp1][0]) * (r*g3Dphi[x][y][0] - I*g3Dphi[x][y][1]) -
				 conj(phip[x][yp1][1]) * (I*g3Dphi[x][y][0] + r*g3Dphi[x][y][1])))
			       -			       
			       (gauge[x][y][1] *
				(conj(phip[x][y][0]) * (r*g3Dphi[x][yp1][0] + I*g3Dphi[x][yp1][1]) +
				 conj(phip[x][y][1]) * (I*g3Dphi[x][yp1][0] - r*g3Dphi[x][yp1][1])))
			       )
			    );
      }
  }
}


//Staggered
int Ainvpsi(Complex psi[LX][LY], Complex b[LX][LY], Complex psi0[LX][LY], const Complex gauge[LX][LY][2], param_t p) {

  int success = 0;
  
  Complex res[LX][LY] , pvec[LX][LY], Apvec[LX][LY];
  double alpha, beta, denom ;
  double rsq = 0, rsqNew = 0, bsqrt = 0.0;
  
  //Intialize  
  zeroField(res);
  zeroField(Apvec);  
  zeroField(pvec);
  
  // Find norm of rhs.
  bsqrt = real(dotField(b,b));
  bsqrt = sqrt(bsqrt);
  
  copyField(res, b); // res = b  - A psi0, for now start with phi0 = 0
  copyField(pvec, res);

  rsq = real(dotField(res,res));
  
  // Compute Ap.
  DdagDpsi(Apvec, pvec, gauge, p);

  // iterate till convergence
  int k;
  for (k=0; k<p.maxIterCG; k++) {
    
    denom = real(dotField(pvec,Apvec));
    alpha = rsq/denom;

    axpy( alpha, pvec, psi);
    axpy(-alpha, Apvec, res);
    
    // Exit if new residual is small enough
    rsqNew = real(dotField(res,res));
    
    if (sqrt(rsqNew) < p.eps*bsqrt) {
      //printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
    // Update vec using new residual
    beta = rsqNew / rsq;
    rsq = rsqNew;
    
    axpy(beta, pvec, res, pvec);
    
    // Compute the new Ap.
    DdagDpsi(Apvec, pvec, gauge,p);  
  }
  //End loop over k

  if(k == p.maxIterCG) {
    printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq); 
    success = 0; // Failed convergence 
  }
  else {
    success = 1; // Convergence 
    k++;
  }

  DdagDpsi(Apvec, psi, gauge,p);
  axpy(-1.0, Apvec, b, res);
  
  //double truersq =  real(dotField(res,res));
  //printf("CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  return success;
}

void forceD(double fD[LX][LY][2], Complex gauge[LX][LY][2], Complex phi[LX][LY],
	    Complex guess[LY][LX], param_t p) {

  if(p.dynamic == true) {

    zeroLat(fD);
    
    Complex phip[LX][LY];
    Complex Dphip[LX][LY];
    zeroField(phip);
    zeroField(Dphip);
    
    Ainvpsi(phip, phi, phip, gauge, p); // note phip = 0 for ODD
    Dpsi(Dphip, phip, gauge, p);        // restrict to Dslash, m = 0
    
    for(int x=0; x<LX; x++)
      for(int y=0; y<LY; y++){
	if( (x+y)%2 == 1) phip[x][y]  = Complex(0.0,0.0);
	if( (x+y)%2 == 0) Dphip[x][y] = Complex(0.0,0.0);
      }
    
    double eta1;
    int yp1, xp1;
    for(int x=0; x<LX; x++)
      for(int y=0; y<LY; y++) {
	
	eta1 =(1.0 - 2.0*(x%2));
	xp1 = (x+1)%LX;
	yp1 = (y+1)%LY;
	
	if( (x+y+1)%2 == 0){ 
	  fD[x][y][0] += 2.0*imag(conj(Dphip[x][y]) * gauge[x][y][0] * phip[xp1][y]);
	}
	else {
	  fD[x][y][0] += 2.0*imag(conj(Dphip[xp1][y]) * conj(gauge[x][y][0]) * phip[x][y]);
	};
	
	if( (x+y+1)%2 == 0){    
	  fD[x][y][1] += 2.0*eta1*imag(conj(Dphip[x][y]) * gauge[x][y][1] * phip[x][yp1]);
	}
	else {
	  fD[x][y][1] += 2.0*eta1*imag(conj(Dphip[x][yp1]) * conj(gauge[x][y][1]) * phip[x][y]);
	}
      }	    
  }
}

#endif
