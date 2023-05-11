




#include "SpectroscopicGalaxyTotal.h"

#include "IntermediateCosmology.h"
#include "LinearPowerSpectrumGR.h"
#include "SpectroscopicGalaxyTensor.h"
#include "SpectroscopicGalaxyVector.h"
#include "SpectroscopicGalaxyScalar.h"




SpectroscopicTotal::SpectroscopicTotal()
{
    //SetPowerSpectrum();
}

SpectroscopicTotal::~SpectroscopicTotal()
{
    //SetPowerSpectrum();
}




double SpectroscopicTotal::GalBispectrumTotal(int x,double fnl,double z,double k1, double k2, double mu12,double phin, double phi, double mu1)
{
	 Spectroscopic SB;
    SpectroscopicVector VB;
    SpectroscopicTensor TB;
	if(x == 1)
{
	double Bnewtonian  = SB.GalBispectrumN(1,9, fnl, z, k1,  k2,  mu12,phi, mu1);
  return Bnewtonian;
}
else if(x == 2)
{
	double Bscalar = SB.GalBispectrumg(1, 10,fnl, z, k1,  k2, mu12, phi, mu1);
	return Bscalar;
}
else if(x == 3)
{
	double Bvector = VB.GalBispectrumV(1,fnl, z, k1,  k2, mu12, phin,phi, mu1);
 return Bvector;
}
else if(x == 4)
{
	double Btensor = TB.GalBispectrumT(1,fnl, z,k1, k2,mu12, phin,  phi,mu1);
	return Btensor;
}
  else
  {
  	double Bnewtonian  = SB.GalBispectrumN(1,9,fnl,z, k1, k2, mu12,phi,  mu1);
    double Bscalar = SB.GalBispectrumg(1,10, fnl, z, k1,  k2, mu12, phi, mu1);
    double Bvector = VB.GalBispectrumV(1,fnl, z, k1,  k2, mu12, phin,phi, mu1);
    double Btensor = TB.GalBispectrumT(1,fnl, z,k1, k2,mu12, phin,  phi,mu1);
  	double tmp = Bnewtonian + Bscalar + Bvector + Btensor;
  return tmp;
  }
  
}


double SpectroscopicTotal::GalBispectrumgTotalAverage(int x, double fnl,double z,double k1, double k2,double muk,double phin)
{
    Spectroscopic SB;
    SpectroscopicVector VB;
    SpectroscopicTensor TB;

  if(x == 1)
{
  double Bnewtonian  = SB.GalBispectrumNAverage(1, 9,  fnl, z,k1,  k2,muk);
  return Bnewtonian;
}
else if(x == 2)
{
  double Bscalar = SB.GalBispectrumgAverage(1,10, fnl, z,k1,  k2, muk);
  return Bscalar;
}
else if(x == 3)
{
  double Bvector = VB.GalBispectrumgVAverage(1,  fnl, z,k1, k2, muk, phin);
 return Bvector;
}
else if(x == 4)
{
  double Btensor = TB.ReducedGalBispectrumgTAverage(1, fnl, z, k1,  k2, muk, phin);
  return Btensor;
}
  else
  {
    double Bnewtonian  = SB.GalBispectrumNAverage(1, 9,  fnl, z,k1,  k2,muk);
    double Bscalar = SB.GalBispectrumgAverage(1,10, fnl, z, k1,  k2, muk);
    double Bvector = VB.GalBispectrumgVAverage(1, fnl, z,k1, k2, muk, phin);
    double Btensor = TB.GalBispectrumgTAverage(1, fnl, z,k1, k2, muk, phin);

    double tmp = Bnewtonian + Bscalar + Bvector + Btensor;
  return tmp;
  }
  

}


double SpectroscopicTotal::ReducedGalBispectrumTotalAverage(int x, double fnl,double z,double k1, double k2,double muk,double phin)
{
   Spectroscopic SB;
    SpectroscopicVector VB;
    SpectroscopicTensor TB;
  if(x == 1)
{
  double Bnewtonian  = SB.ReducedGalBispectrumNAverage(1, 9,  fnl, z, k1, k2,muk);
  return Bnewtonian;
}
else if(x == 2)
{
  double Bscalar = SB.ReducedGalBispectrumgAverage(1,10, fnl, z,k1, k2, muk);
  return Bscalar;
}
else if(x == 3)
{
  double Bvector = VB.ReducedGalBispectrumgVAverage(1, fnl, z, k1,  k2, muk, phin);
 return Bvector;
}
else if(x == 4)
{
  double Btensor = TB.ReducedGalBispectrumgTAverage(1,  fnl, z, k1,k2,muk, phin);
  return Btensor;
}
  else
  {
    double Bnewtonian  = SB.ReducedGalBispectrumNAverage(1, 9,  fnl, z, k1, k2,muk);
    double Bscalar = SB.ReducedGalBispectrumgAverage( 1,10, fnl, z,k1, k2, muk);
    double Bvector = VB.ReducedGalBispectrumgVAverage(1, fnl, z, k1,  k2, muk, phin);
    double Btensor = TB.ReducedGalBispectrumgTAverage(1,  fnl, z, k1,k2,muk, phin);
    double tmp = Bnewtonian + Bscalar + Bvector + Btensor;
  return tmp;
  }
  

}

