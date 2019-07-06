#ifndef HMC_HELPERS
#define HMC_HELPERS

//2D Routines
//-------------------------------------------------------------------------------------
void forceU(double fU[LX][LY][2], Complex gauge[LX][LY][2], param_t p) {
  
  Complex plaq0;
  Complex plaq;
  zeroLat(fU);
  int xp1, xm1, yp1, ym1;
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++) {

      xp1 = (x+1)%LX;
      xm1 = (x-1+LX)%LX;
      yp1 = (y+1)%LY;
      ym1 = (y-1+LY)%LY;
      
      plaq0 = gauge[x][y][0]*gauge[xp1][y][1]*conj(gauge[x][yp1][0])*conj(gauge[x][y][1]);
      fU[x][y][0] += p.beta*imag(plaq0);
      
      plaq =  gauge[x][ym1][0]*gauge[xp1][ym1][1]*conj(gauge[x][y][0])*conj(gauge[x][ym1][1]);
      fU[x][y][0] -= p.beta*imag(plaq);

      plaq =  gauge[x][y][1]*conj(gauge[xm1][yp1][0])*conj(gauge[xm1][y][1])*gauge[xm1][y][0];
      fU[x][y][1] += p.beta*imag(plaq);

      //This plaquette was aleady computed. We want the conjugate.
      fU[x][y][1] -= p.beta*imag(plaq0);      
    }
}

//P_{k+1/2} = P_{k-1/2} - dtau * (fU + fD)
void update_mom(double fU[LX][LY][2], double fD[LX][LY][2], double mom[LX][LY][2], double dtau){

  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int mu=0; mu<2; mu++)
	mom[x][y][mu] -= (fU[x][y][mu] - fD[x][y][mu])*dtau;
}

//U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
void update_gauge(Complex gauge[LX][LY][2], double mom[LX][LY][2], double dtau){
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int mu=0; mu<2; mu++)
	gauge[x][y][mu] *= polar(1.0, mom[x][y][mu] * dtau);
}
//----------------------------------------------------------------------------------
#endif
