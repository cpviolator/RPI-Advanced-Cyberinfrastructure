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
  
  zeroLat(fD);
  
  if(p.dynamic == true) {
    
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

#endif
