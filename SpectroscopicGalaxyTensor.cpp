

#include "SpectroscopicGalaxyTensor.h"

//#include "GrowthEquation.h"
#include "MyBessel.cpp"

#include <boost/bind.hpp>
#include <stdio.h>
#include <stdlib.h>



#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include "csimpson.cxx"

#include "IntermediateCosmology.h"
#include "LinearPowerSpectrumGR.h"





SpectroscopicTensor::SpectroscopicTensor()
{
    //SetPowerSpectrum();
}

SpectroscopicTensor::~SpectroscopicTensor()
{
    //SetPowerSpectrum();
}






///////////////////////Tensor Bispectrum begins/////////////////////////
double GreenFunction(double z, double k)
{
   Cosmology C;
   const  double Mpctokm =  3.086e19;
   const  double c_light_km_s =  299792;
   const  double c_light_Mpc_s =  299792.0/Mpctokm;
   const double onegigayeartosec = 315576*pow(10,11);
   const double factor = c_light_Mpc_s*onegigayeartosec;
   double eta = C.ConformalTime(z)*factor;
    double x = k*eta;
    double jone =  3.0* jl2(1, k*eta)/(k*eta);
   //double TransferFunction  = 1.0 + 3.0*(k*eta*cos(k*eta)-sin(k*eta))/(pow(k*eta,3));
   double TransferFunction  = 1.0 + jone;
   return TransferFunction;
}

double GreenFunctionpr(double z, double k)
{
  //Cosmology C;
	IntermediateCosmology TC;


	double delta = 1.0e-5;
	double x1 = z -delta;
	double x2 = z + delta ;
	double y1 =GreenFunction(x1,k);
	double y2 = GreenFunction(x2,k);
	double deri = (y2-y1)/(x2-x1);
	double tmp = -(1.0 + z)*TC.HH(z)*deri;
	return tmp;
}


double CurylyA4(double z)
{   
	IntermediateCosmology TC;

	double tmp1 = (3.0*TC.Omegam(z) + 2.0*pow(TC.fg(z),2));
	double tmp2 = (2.0*TC.fg(z) -TC.curlyQ(z) + 2.0*TC.Hpr(z)/pow(TC.HH(z),2));
	double tmp3 = 4.0*TC.fg(z)*TC.fgpr(z)/(TC.HH(z));
	double tmp4 = -3.0*TC.Omegam(z)*(1.0 + 2.0*TC.Hpr(z)/pow(TC.HH(z),2));
	double tmp = tmp1*tmp2 + tmp3  + tmp4;
	return tmp;
}
double CurylyA5(double z)
{
	IntermediateCosmology TC;
	double tmp = (3.0*TC.Omegam(z) + 2.0*pow(TC.fg(z),2));
	return tmp;
}

double CurlyKTk1k2k3(double z,double k1, double k2, double k3, double mu12, double mu1, double mu2, double mu3, double phi, double phin)
{       
	IntermediateCosmology TC;

    double shfitsqrt  = sqrt(1.0-0.998*0.998);
	double tmp1 = (1.0-pow(mu1,2))*pow(sin(phin-phi),2);
	double tmp2A =  ((k1/k2) + mu12)*(1.0-pow(mu1,2))*pow(cos(phin-phi),2);
	double tmp2B =  mu1*mu1*(1.0 -pow(mu12,2)); 

	double tmp2C = -2.0*((k1/k2) + mu12)*mu1*sqrt(1.0-pow(mu1,2))*(sqrt(1.0-pow(mu12,2))-shfitsqrt)*cos(phin-phi);
	double tmp2  =  -pow(k2/k3,2)*(tmp2A + tmp2B + tmp2C) ;
       double tmp = (1.0-pow(mu12,2))*(tmp1 + tmp2)/(2.0*pow(k3,2));

      double everything = 3.0*TC.Omegam(z)*pow(TC.HH(z),4)*(CurylyA4(z)*GreenFunction(z, k3) + CurylyA5(z)*GreenFunctionpr(z,k3)/TC.HH(z))/pow(k3,2);
	return everything*tmp;
}

double CurlyKTk2k3k1(double z,double k1, double k2, double k3, double mu12, double mu1, double mu2, double mu3, double phi, double phin)
{
	IntermediateCosmology TC;

    double cos2A =  pow(cos(phin-phi),2) - pow(sin(phin-phi),2);
	double tmp =  -(1.0 - pow(mu1,2))*(1.0 - pow(mu12,2))*cos2A/(2.0*k3*k3);

	double everything = 3.0*TC.Omegam(z)*pow(TC.HH(z),4)*(CurylyA4(z)*GreenFunction(z, k1) + CurylyA5(z)*GreenFunctionpr(z,k1)/(TC.HH(z)))/pow(k1,2);
	return everything*tmp;
}

double CurlyKTk3k1k2(double z,double k1, double k2, double k3, double mu12, double mu1, double mu2, double mu3, double phi, double phin)
{
	IntermediateCosmology TC;

    double shfitsqrt  = sqrt(1.0-0.998*0.998);
	double tmp1 = (1.0-pow(mu1,2))*pow(sin(phin-phi),2);
	double tmp2 =  -mu12*mu12*(1.0-pow(mu1,2))*pow(cos(phin-phi),2);
	double tmp3 = 2.0*mu1*mu12*sqrt(1.0-pow(mu1,2))*(sqrt(1.0-pow(mu12,2)))*cos(phin-phi);
	double tmp4  =  -mu1*mu1*(1.0 -pow(mu12,2)); 
    double tmp = (1.0-pow(mu12,2))*(tmp1 + tmp2 + tmp3 + tmp4)/(2.0*pow(k3,2));

	double everything = 3.0*TC.Omegam(z)*pow(TC.HH(z),4)*(CurylyA4(z)*GreenFunction(z, k2) + CurylyA5(z)*GreenFunctionpr(z,k2)/(TC.HH(z)))/pow(k2,2);
	return everything*tmp;
}



double kernelTk1k2k3(int y,double fnl, double z, double k1, double k2, double k3, double mu12, double mu1, double mu2, double mu3, double phi, double phin)
{
  
   LinearPowerSpectrumGR LP;
  double KReal2 = CurlyKTk1k2k3(z,k1, k2,k3,mu12,mu1,mu2,mu3,phi, phin);
  double KReal1A = LP.Linearkernel(10, fnl,z, k1,mu1);
  double KReal1B = LP.Linearkernel(10, fnl,z, k2,mu2);
  double KImaginary2 = 0.0;
  double KImaginary1A = LP.Linearkernel(3, fnl,z, k1,mu1);
  double KImaginary1B = LP.Linearkernel(3, fnl,z, k2,mu2);

if(y == 1)
{
 double tmp1 = KReal2*(KReal1A*KReal1B - KImaginary1A*KImaginary1B);// - KImaginary2*(KReal1A*KImaginary1B + KImaginary1A*KReal1B);
return  tmp1*LP.P1P2(z,k1,k2);
}
else if(y == 2)
{
double tmp2 = KReal2*(KReal1A*KImaginary1B +  KImaginary1A*KReal1B);// + KImaginary2*(KReal1A*KReal1B - KImaginary1A*KImaginary1B);
return tmp2*LP.P1P2(z,k1,k2);
}
else 
return 0.0;
}

double kernelTk2k3k1(int y,double fnl, double z, double k1, double k2, double k3, double mu12, double mu1, double mu2, double mu3, double phi, double phin)
{
	LinearPowerSpectrumGR LP;
  double KReal2 = CurlyKTk2k3k1(z,k1, k2,k3,mu12,mu1,mu2,mu3,phi, phin);
  double KReal1A = LP.Linearkernel(10, fnl,z, k2,mu2);
  double KReal1B = LP.Linearkernel(10, fnl,z, k3,mu3);


  double KImaginary2 = 0.0;//KV231Imaginary(k1, k2, k3,  mu12,mu1, mu2,  mu3,  phi,  phin);
  double KImaginary1A = LP.Linearkernel(3, fnl,z, k2,mu2);
  double KImaginary1B = LP.Linearkernel(3, fnl,z, k3,mu3);

if(y == 1)
{
 double tmp1 = KReal2*(KReal1A*KReal1B - KImaginary1A*KImaginary1B);// - KImaginary2*(KReal1A*KImaginary1B + KImaginary1A*KReal1B);
return  tmp1*LP.P1P2(z,k2,k3);
}
else if(y == 2)
{
double tmp2 = KReal2*(KReal1A*KImaginary1B +  KImaginary1A*KReal1B);// + KImaginary2*(KReal1A*KReal1B - KImaginary1A*KImaginary1B);
return tmp2*LP.P1P2(z,k2,k3);
}
else 
return 0.0;
}


double kernelTk3k1k2(int y,double fnl, double z, double k1, double k2, double k3, double mu12, double mu1, double mu2, double mu3, double phi, double phin)
{
	LinearPowerSpectrumGR LP;
  double KReal2 = CurlyKTk3k1k2(z,k1, k2,k3,mu12,mu1,mu2,mu3,phi, phin);
  double KReal1A = LP.Linearkernel(10, fnl,z, k3,mu3);
  double KReal1B = LP.Linearkernel(10, fnl,z, k1,mu1);


  double KImaginary2 = 0.0;//KV312Imaginary(k1, k2, k3,  mu12,mu1, mu2,  mu3,  phi,  phin);
  double KImaginary1A = LP.Linearkernel(3, fnl,z, k3,mu3);
  double KImaginary1B = LP.Linearkernel(3, fnl,z, k1,mu1);

if(y == 1)
{
 double tmp1 = KReal2*(KReal1A*KReal1B - KImaginary1A*KImaginary1B);// - KImaginary2*(KReal1A*KImaginary1B + KImaginary1A*KReal1B);
return  tmp1*LP.P1P2(z,k3,k1);
}
else if(y == 2)
{
double tmp2 = KReal2*(KReal1A*KImaginary1B +  KImaginary1A*KReal1B);// + KImaginary2*(KReal1A*KReal1B - KImaginary1A*KImaginary1B);
return tmp2*LP.P1P2(z,k3,k1);
}
else 
return 0.0;

}

double SpectroscopicTensor::GalBispectrumT(int x,double fnl,double z,double k1, double k2, double mu12,double phin, double phi, double mu1)
{

	LinearPowerSpectrumGR LP;
    double r = k2/k1;
    double y = sqrt(r*r+2.0*mu12 *r + 1.0);
    double k3 =  k1*y;

    double mu2 = LP.mu2_inVT(phi,phin, mu1, mu12);
    double mu3 = LP.mu3_inVT(k1,k2, k3,phi, phin, mu1, mu12);

  //  double mu2 = mu2_in(phi,  mu1,  mu12);
   // double mu3 = mu3_in(k1,k2, k3, phi, mu1, mu12);// mu13(r, mu12)//mu23(r, mu12)
///kernelVRandIk1k2k3(y,fnl,z,k1,k2, k3, mu12,mu1,mu2, mu3, phi, phin)
  
    double tmp1 = kernelTk1k2k3(x,fnl,z,k1,k2, k3, mu12,mu1,mu2, mu3, phi, phin);
    double tmp2 = kernelTk2k3k1(x,fnl,z,k1,k2, k3, mu12,mu1,mu2, mu3, phi, phin);
    double tmp3 = kernelTk3k1k2(x,fnl,z,k1,k2, k3, mu12,mu1,mu2, mu3, phi, phin);

    return tmp1 + tmp2 + tmp3;

}


double SpectroscopicTensor::PassbyrefGalT(int x, void *params, double phin, double phi, double mu1)
{
  double *fp = (double *) params;
  double fnl = fp[0];
  double  z = fp[1];
  double k1 = fp[2];
  double k2 = fp[3];
 // double k3 = fp[4];
  double mu12 = fp[4];
  //GalBispectrumV(int x,double fnl,double z,double k1, double k2, double mu12,double phin, double phi, double mu1)
return GalBispectrumT(x, fnl, z, k1,k2,mu12,phin,phi,mu1);


}

double SpectroscopicTensor::GalBispectrumgTAverage(int x, double fnl,double z,double k1, double k2,double muk,double phin)
{
    const double EPSREL = 1e-3;
    const double EPSABS = 1e-3;
    double a[2] = {0.00001,-0.99998};
    double  b[2] = {2*M_PI,0.99998};
    double fp[5];
    fp[0] =fnl;
    fp[1] = z;
    fp[2] = k1;
    fp[3] = k2;
   // fp[4] = k3;
    fp[4] = muk;

    double tmp1 = Integrate<2>(bind(&SpectroscopicTensor::PassbyrefGalT,this,x,fp,phin,_1,_2),a,b,EPSREL,EPSABS)/(4.0*M_PI);
    //double tmp1 = Integrate(bind(&Spectroscopic::PassbyrefGalT,this,x,fp,phin,M_PI/2,_1),-1.0,1.0,EPSREL,EPSABS)/(4.0*M_PI); 
    return tmp1;
}

double SpectroscopicTensor::ReducedGalBispectrumgTAverage(int x, double fnl,double z,double k1, double k2,double muk,double phin)
{
	LinearPowerSpectrumGR LP;

    double r = k2/k1;
   double y = sqrt(r*r+2.0*muk *r + 1.0);
   double k3 =  k1*y;
   double PowerSPec_perm12 = LP.LinearPowerAverage(3,fnl, z,k1)*LP.LinearPowerAverage(3,fnl,z,k2);
   double PowerSPec_perm32 = LP.LinearPowerAverage(3,fnl, z, k3)*LP.LinearPowerAverage(3,fnl,z,k2);
   double PowerSPec_perm13 = LP.LinearPowerAverage(3,fnl,z, k1 )*LP.LinearPowerAverage(3,fnl,z,k3);

return GalBispectrumgTAverage(x, fnl,z,k1, k2,muk,phin)/(PowerSPec_perm12+PowerSPec_perm32+PowerSPec_perm13);

}









///////////////////////////////////////Tensor bispectrum ends////////////////
