#ifndef DOPHELPERS_H
#define DOPHELPERS_H

// Wilson stencil
// D_{W}(n,m) = (m_{0} + 2r)\delta(n,m)
//               - 1/2 Sum [(r-\sigma_{mu}) U_{n,\mu} \delta_{n,m-\hat{\mu}} +
//                          (r+\sigma_{mu}) U^{\dagger}_{m,\mu} \delta_{n,m+\hat{\mu}}]
//
// sigma_1 = | 0  1 |  sigma_2 = | 0 -i | sigma_3 = i*sigma_1*sigma_2 = | 1  0 |
//           | 1  0 |            | i  0 |                               | 0 -1 |

void Dpsi(Complex psi2[LX][LY][2], Complex psi1[LX][LY][2],
	  const Complex gauge[LX][LY][2], param_t p){
  
  double m0 = p.m;
  double  r = 1.0;
  double constant = (2*r + m0);
  int xp1, xm1, yp1, ym1;
  
  //Sum over 0,1 directions.
  for(int x=0; x<LX; x++) {
    xp1 = (x+1)%LX;
    xm1 = (x-1+LX)%LX;
    for(int y=0; y<LY; y++) {
      yp1 = (y+1)%LY;
      ym1 = (y-1+LY)%LY;
      
      //upper
      psi2[x][y][0] = constant * psi1[x][y][0] -
	
	0.5*(     gauge[x][y][0]    * (r*psi1[xp1][y][0] - psi1[xp1][y][1]) +
	     conj(gauge[xm1][y][0]) * (r*psi1[xm1][y][0] + psi1[xm1][y][1]) +
		  
		  gauge[x][y][1]    * (r*psi1[x][yp1][0] + I*psi1[x][yp1][1]) +
	     conj(gauge[x][ym1][1]) * (r*psi1[x][ym1][0] - I*psi1[x][ym1][1]));
      
      //lower
      psi2[x][y][1] = constant * psi1[x][y][1] -
	
	0.5*(     gauge[x][y][0]    * (-psi1[xp1][y][0] + r*psi1[xp1][y][1]) -
	     conj(gauge[xm1][y][0]) * (-psi1[xm1][y][0] - r*psi1[xm1][y][1]) +

		  gauge[x][y][1]    * (-I*psi1[x][yp1][0] + r*psi1[x][yp1][1]) -
	     conj(gauge[x][ym1][1]) * (-I*psi1[x][ym1][0] - r*psi1[x][ym1][1]));
      
      
    }
  }
}

void g3psi(Complex psi2[LX][LY][2], Complex psi1[LX][LY][2]){
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++) {
      psi2[x][y][0] =  psi1[x][y][0];
      psi2[x][y][1] = -psi1[x][y][1];
    }
}

void g3psi(Complex psi1[LX][LY][2]){
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++) {
      psi1[x][y][1] *= -1.0;
    }
}

void g2psi(Complex psi2[LX][LY][2], Complex psi1[LX][LY][2]){
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++) {
      psi2[x][y][0] = -I*psi1[x][y][1];
      psi2[x][y][1] =  I*psi1[x][y][0];
    }
}

void g2psi(Complex psi1[LX][LY][2]){

  Complex tmp = 0.0;
  for(int x=0; x<LX; x++)    
    for(int y=0; y<LY; y++) {
      tmp = psi1[x][y][1];      
      psi1[x][y][1] =  I*psi1[x][y][0];
      psi1[x][y][0] = -I*tmp;
    }
}


void g1psi(Complex psi2[LX][LY][2], Complex psi1[LX][LY][2]){
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++) {
      psi2[x][y][0] = psi1[x][y][1];
      psi2[x][y][1] = psi1[x][y][0];
    }
}

void g1psi(Complex psi1[LX][LY][2]){

  Complex tmp = 0.0;
  for(int x=0; x<LX; x++)    
    for(int y=0; y<LY; y++) {
      tmp = psi1[x][y][1];      
      psi1[x][y][1] = psi1[x][y][0];
      psi1[x][y][0] = tmp;
    }
}


void g3Dpsi(Complex psi2[LX][LY][2], Complex psi1[LX][LY][2],
	    const Complex gauge[LX][LY][2], param_t p ){

  Dpsi(psi2, psi1, gauge, p);
  g3psi(psi2);
}

void Ddagpsi(Complex psi2[LX][LY][2], Complex psi1[LX][LY][2],
	     const Complex gauge[LX][LY][2], param_t p){
  
  g3psi(psi1);
  Dpsi(psi2, psi1, gauge, p);
  g3psi(psi2);  
}


void DdagDpsi(Complex psi2[LX][LY][2], Complex psi1[LX][LY][2],
	      const Complex gauge[LX][LY][2], param_t p) {
  
  //Hack for now
  Complex temp[LX][LY][2];
  Dpsi(temp, psi1, gauge, p);
  g3psi(temp);
  Dpsi(psi2, temp, gauge, p);
  g3psi(psi2);
}

/*
  For a 2D square lattice, the stencil is:

  psi2[x] = D_xy psi_y = m psi_x delta_xy 
  - eta_mu(x) [U(x,x+mu) psi[x+mu] - U*(x-mu,x) psi[x-mu]]

  The even/odd anti-hermiticity             

  eta_mu(x) [U(x,y)\deta_x+mu,y - U(y,x) delta_y+mu,
 
  1 |  0 -eta1  0 |
  - | +1    0  -1 |  , where eta0 = 1, eta1 = (-)^x = 1 - 2*(x%2)
  2 |  0 +eta1  0 |

*/

void Dpsi(Complex psi2[LX][LY], Complex psi1[LX][LY],
	  const Complex gauge[LX][LY][2], param_t p ){

  int xp1, xm1, yp1, ym1;
  double eta1;  
  
  for(int x=0; x<LX; x++) {
    eta1 =(1-2*(x%2));
    xp1 = (x+1)%LX;
    xm1 = (x-1+LX)%LX;    
    for(int y=0; y<LY; y++) {
      yp1 = (y+1)%LY;
      ym1 = (y-1+LY)%LY;

      psi2[x][y] = p.m*psi1[x][y];
      
      psi2[x][y] += - (gauge[x][y][0] * psi1[xp1][y]
		       - conj(gauge[xm1][y][0]) * psi1[xm1][y]);
      
      psi2[x][y] += - eta1*(gauge[x][y][1]*psi1[x][yp1]
			    - conj(gauge[x][ym1][1])*psi1[x][ym1]);
      
    }
  }
}

void Ddagpsi(Complex psi2[LX][LY], Complex  psi1[LX][LY],
	     const Complex gauge[LX][LY][2], param_t p ) {
  
  // For a 2D square lattice, the stencil is:
  //   1 |  0 -eta1  0 |
  //   - | +1    0  -1 |  , where eta0 = 1, eta1 = (-)^x = 1 - 2(x%L)
  //   2 |  0 +eta1  0 |

  int xp1, xm1, yp1, ym1;
  double eta1;

  for(int x=0; x<LX; x++) {
    eta1 =(1-2*(x%2));
    xp1 = (x+1)%LX;
    xm1 = (x-1+LX)%LX;    
    for(int y=0; y<LY; y++) {
      yp1 = (y+1)%LY;
      ym1 = (y-1+LY)%LY;

      psi2[x][y] = p.m*psi1[x][y];
      
      psi2[x][y] += (gauge[x][y][0] * psi1[xp1][y]
		     - conj(gauge[xm1][y][0]) * psi1[xm1][y]);
      
      psi2[x][y] += eta1*(gauge[x][y][1]*psi1[x][yp1]
			  - conj(gauge[x][ym1][1])*psi1[x][ym1]);
    }
  }
}

//=======================//
// Note: Ddag D = D Ddag //
///======================//

void DdagDpsi(Complex psi2[LX][LY], Complex  psi1[LX][LY],
	      const Complex gauge[LX][LY][2], param_t p ) {
  
  Complex psitmp[LX][LY];
  zeroField(psitmp);
  Dpsi(psitmp, psi1, gauge, p);
  Ddagpsi(psi2, psitmp, gauge, p);
}



#endif
