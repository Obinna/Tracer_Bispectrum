
#include "SpectroscopicGalaxyVector.h"

//#include "GrowthEquation.h"
//#include "MyBessel.cpp"

#include <boost/bind.hpp>
#include <stdio.h>
#include <stdlib.h>



#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
//#include "csimpson.cxx"

#include "IntermediateCosmology.h"
#include "LinearPowerSpectrumGR.h"


SpectroscopicVector::SpectroscopicVector()
{
    //SetPowerSpectrum();
}

SpectroscopicVector::~SpectroscopicVector()
{
    //SetPowerSpectrum();
}



//////////////////////////////Tools//////////////////////////////////////



double muk13(double k1, double k2, double k3)
{
     double tmp = (k2*k2-k1*k1-k3*k3)/(2.0* k1*k3);
     return tmp;

}

double muk23(double k1, double k2, double k3)
{
    double tmp = (k1*k1-k2*k2-k3*k3)/(2.0* k2*k3);
     return tmp;

}



//////////////////Vector bispectrum begins///////////
double CurlyA2(double z)
{
	IntermediateCosmology TC;

	double tmp = 12.0*TC.Omegam(z)*pow(TC.HH(z),2)*TC.fg(z);
	return tmp;
}

double CurlyA3(double z)
{
	IntermediateCosmology TC;
	double tmp1 =  -TC.bevo(z) + 2.0*TC.curlyQ(z)+ 2.0*(1.0 - TC.curlyQ(z))/(TC.chi(z)*TC.HH(z)) + TC.Hpr(z)/pow(TC.HH(z),2);
	double tmp = -12.0*TC.Omegam(z)*pow(TC.HH(z),3)*TC.fg(z)*tmp1;
	return tmp;
}

double curlyV123(double z,double k1, double k2, double k3, double mu12, double mu1, double mu2, double mu3, double phi, double phin)
{
  // double r = k2/k1;
   double shfitsqrt  = sqrt(1.0-0.998*0.998);
   double tmp1 = (k2*k2-k1*k1)/(pow(k3,2)*k1);
   double tmp2 = sqrt(1.0-mu12*mu12) - shfitsqrt;
   double tmp3 = mu1*sqrt(1.0-mu12*mu12)- (k1/k2 + mu12)*sqrt(1.0-mu1*mu1)*cos(phin-phi);
   double tmp =  tmp1*tmp2*tmp3;
   return tmp;
}

double curlyV231(double z, double k1, double k2, double k3, double mu12, double mu1, double mu2, double mu3, double phi, double phin)
{
   //double r = k2/k1;
   double shfitsqrt  = sqrt(1.0-0.998*0.998);
   double tmp1 = (k3*k3 - k2*k2)/(pow(k3,2)*k2);
   double tmp2 = sqrt(1.0-mu1*mu1);
   double tmp3 = (sqrt(1.0-mu12*mu12)-shfitsqrt)*cos(phin-phi);
   double tmp =  tmp1*tmp2*tmp3;
   return tmp;
}


double curlyV312(double z,double k1, double k2, double k3, double mu12, double mu1, double mu2, double mu3, double phi, double phin)
{
  // double r = k2/k1;
  double mu12arg  = 1.0 - pow(mu12,2);//sqrt(abs(mu12arg) < 1e-12 ? 0 : mu12arg)
  double mu1arg = 1.0 - pow(mu1,2);//
   double sqrtmu1 = sqrt(abs(mu1arg));
    double sqrtmu12 = sqrt(abs(mu12arg) );
   double tmp = 0.0;
   double tmp1 = (k3*k3 - k1*k1)/(pow(k3,2)*k1);
   //double tmp2 = sqrt(1 - mu12*mu12) -shfitsqrt;
   double tmp2 = sqrtmu12;
   double tmp3 = mu1*(sqrtmu12)- mu12*(sqrtmu1)*cos(phin-phi);
          tmp =  tmp1*tmp2*tmp3;
  //cout<<"\t"<<"k1"<<"\t"<<k1<<"\t"<<"mu12"<<"\t"<<mu12<<"\t"<<"CurlyV312"<<"\t"<<tmp<<endl;
   return tmp;
}

double KV123Real(double z,double k1, double k2, double k3, double mu12, double mu1, double mu2, double mu3, double phi, double phin)
{
	double tmp = CurlyA2(z)*mu3*curlyV123(z,k1,k2,k3,mu12,mu1,mu2,mu3,phi,phin)/k3;
	return tmp;
}
double KV123Imaginary(double z,double k1, double k2, double k3, double mu12, double mu1, double mu2, double mu3, double phi, double phin)
{
	double tmp = CurlyA3(z)*curlyV123(z,k1,k2,k3,mu12,mu1,mu2,mu3,phi,phin)/pow(k3,2);
	return tmp;
}

double KV231Real(double z,double k1, double k2, double k3, double mu12, double mu1, double mu2, double mu3, double phi, double phin)
{
	double tmp = CurlyA2(z)*mu1*curlyV231(z,k1,k2,k3,mu12,mu1,mu2,mu3,phi,phin)/k1;
	return tmp;
}
double KV231Imaginary(double z,double k1, double k2, double k3, double mu12, double mu1, double mu2, double mu3, double phi, double phin)
{
	double tmp = CurlyA3(z)*curlyV231(z,k1,k2,k3,mu12,mu1,mu2,mu3,phi,phin)/pow(k1,2);
	return tmp;
}

double KV312Real(double z,double k1, double k2, double k3, double mu12, double mu1, double mu2, double mu3, double phi, double phin)
{
	double tmp = CurlyA2(z)*mu2*curlyV312(z,k1,k2,k3,mu12,mu1,mu2,mu3,phi,phin)/k2;
	return tmp;
}
double KV312Imaginary(double z,double k1, double k2, double k3, double mu12, double mu1, double mu2, double mu3, double phi, double phin)
{
	double tmp = CurlyA3(z)*curlyV312(z,k1,k2,k3,mu12,mu1,mu2,mu3,phi,phin)/pow(k2,2);
	return tmp;
}



////////////////////////////////////////////Bispectrum computation proper///////////////////////////////

double kernelVRandIk1k2k3(int y,double fnl, double z, double k1, double k2, double k3, double mu12, double mu1, double mu2, double mu3, double phi, double phin)
{
	LinearPowerSpectrumGR LP;
  double KReal2 = KV123Real(z,k1, k2,k3,mu12,mu1,mu2,mu3,phi, phin);
  /// SecondOrderKernelRSD(18, fnl, z, k1, k2,k3,mu12, mu1, mu2,mu3);
  double KReal1A = LP.Linearkernel(10, fnl,z, k1,mu1);
  double KReal1B = LP.Linearkernel(10, fnl,z, k2,mu2);


  double KImaginary2 = KV123Imaginary(z,k1, k2, k3,  mu12,mu1, mu2,  mu3,  phi,  phin);
  //SecondOrderKernelRSD(11, fnl, z, k1, k2,k3,mu12, mu1, mu2,mu3);
  double KImaginary1A = LP.Linearkernel(3, fnl,z, k1,mu1);
  double KImaginary1B = LP.Linearkernel(3, fnl,z, k2,mu2);

if(y == 1)
{
 double tmp1 = KReal2*(KReal1A*KReal1B - KImaginary1A*KImaginary1B) - KImaginary2*(KReal1A*KImaginary1B + KImaginary1A*KReal1B);
return  tmp1*LP.P1P2(z,k1,k2);
}
else if(y == 2)
{
double tmp2 = KReal2*(KReal1A*KImaginary1B +  KImaginary1A*KReal1B) + KImaginary2*(KReal1A*KReal1B - KImaginary1A*KImaginary1B);
return tmp2*LP.P1P2(z,k1,k2);
}
else 
return 0.0;
}

double kernelVRandIk2k3k1(int y,double fnl, double z, double k1, double k2, double k3, double mu12, double mu1, double mu2, double mu3, double phi, double phin)
{
	LinearPowerSpectrumGR LP;
  double KReal2 = KV231Real(z,k1, k2,k3,mu12,mu1,mu2,mu3,phi, phin);
  double KReal1A = LP.Linearkernel(10, fnl,z, k2,mu2);
  double KReal1B = LP.Linearkernel(10, fnl,z, k3,mu3);


  double KImaginary2 = KV231Imaginary(z,k1, k2, k3,  mu12, mu1, mu2,  mu3,  phi,  phin);
  double KImaginary1A = LP.Linearkernel(3, fnl,z, k2,mu2);
  double KImaginary1B = LP.Linearkernel(3, fnl,z, k3,mu3);

if(y == 1)
{
 double tmp1 = KReal2*(KReal1A*KReal1B - KImaginary1A*KImaginary1B) - KImaginary2*(KReal1A*KImaginary1B + KImaginary1A*KReal1B);
return  tmp1*LP.P1P2(z,k2,k3);
}
else if(y == 2)
{
double tmp2 = KReal2*(KReal1A*KImaginary1B +  KImaginary1A*KReal1B) + KImaginary2*(KReal1A*KReal1B - KImaginary1A*KImaginary1B);
return tmp2*LP.P1P2(z,k2,k3);
}
else 
return 0.0;
}


double kernelVRandIk3k1k2(int y,double fnl, double z, double k1, double k2, double k3, double mu12, double mu1, double mu2, double mu3, double phi, double phin)
{

	LinearPowerSpectrumGR LP;
  double KReal2 = KV312Real(z,k1, k2,k3,mu12,mu1,mu2,mu3,phi, phin);
  double KReal1A = LP.Linearkernel(10, fnl,z, k3,mu3);
  double KReal1B = LP.Linearkernel(10, fnl,z, k1,mu1);


  double KImaginary2 = KV312Imaginary(z,k1, k2, k3,  mu12,mu1, mu2,  mu3,  phi,  phin);
  double KImaginary1A = LP.Linearkernel(3, fnl,z, k3,mu3);
  double KImaginary1B = LP.Linearkernel(3, fnl,z, k1,mu1);

if(y == 1)
{
 double tmp1 = KReal2*(KReal1A*KReal1B - KImaginary1A*KImaginary1B) - KImaginary2*(KReal1A*KImaginary1B + KImaginary1A*KReal1B);
return  tmp1*LP.P1P2(z,k3,k1);
}
else if(y == 2)
{
double tmp2 = KReal2*(KReal1A*KImaginary1B +  KImaginary1A*KReal1B) + KImaginary2*(KReal1A*KReal1B - KImaginary1A*KImaginary1B);
return tmp2*LP.P1P2(z,k3,k1);
}
else 
return 0.0;
}


double SpectroscopicVector::GalBispectrumV(int x,double fnl,double z,double k1, double k2, double mu12,double phin, double phi, double mu1)
{
	LinearPowerSpectrumGR LP;
    double r = k2/k1;
    double y = sqrt(r*r+2.0*mu12 *r + 1.0);
    double k3 =  k1*y;

    double mu2 = LP.mu2_inVT(phi,phin, mu1, mu12);
    double mu3 = LP.mu3_inVT(k1,k2, k3,phi, phin, mu1, mu12);// mu13(r, mu12)//mu23(r, mu12)
///kernelVRandIk1k2k3(y,fnl,z,k1,k2, k3, mu12,mu1,mu2, mu3, phi, phin)
  
    double tmp1 = kernelVRandIk1k2k3(x,fnl,z,k1,k2, k3, mu12,mu1,mu2, mu3, phi, phin); // + kernelVRandIk2k3k1(x,fnl,z,k1,k2, k3, mu12,mu1,mu2, mu3, phi, phin) +  kernelVRandIk3k1k2(x,fnl,z,k1,k2, k3, mu12,mu1,mu2, mu3, phi, phin);
    double tmp2 = kernelVRandIk2k3k1(x,fnl,z,k1,k2, k3, mu12,mu1,mu2, mu3, phi, phin);   
    double tmp3 = kernelVRandIk3k1k2(x,fnl,z,k1,k2, k3, mu12,mu1,mu2, mu3, phi, phin);
    //return tmp3;

  //cout<<"\t"<<"k1"<<"\t"<<k1<<"\t"<<"tmp1"<<"\t"<<tmp1<<"\t"<<"tmp2"<<"\t"<<"\t"<<tmp2<<"\t"<<"tmp3"<<"\t"<<"\t"<<tmp3<<endl;

    return tmp1 + tmp2 + tmp3;

}




double SpectroscopicVector::PassbyrefGalV(void *params,double k1, double k2,double phi, double mu1)
{
  double *fp = (double *) params;
  double fnl = fp[0];
  double  z = fp[1];
  double x = fp[2];
  double phin = fp[3];
 // double k3 = fp[4];
  double mu12 = fp[4];
  //GalBispectrumV(int x,double fnl,double z,double k1, double k2, double mu12,double phin, double phi, double mu1)
return GalBispectrumV(x, fnl, z, k1,k2,mu12,phin,phi,mu1);


}





double SpectroscopicVector::GalBispectrumgVAverage(int x, double fnl,double z,double k1, double k2,double muk,double phin)
{
    const double EPSREL = 1e-3;
    const double EPSABS = 1e-3;
    double a[2] = {0.0,-0.9998};
    double  b[2] = {2*M_PI,1.0};

    double a3[3] = {0.00001,0.00001,-1.0};
    double  b3[3] = {2*M_PI,2*M_PI,1.0};
    double fp[4];
    fp[0] =fnl;
    fp[1] = z;
    fp[2] = x;
    fp[3] = phin;
   // fp[4] = k3;
    fp[4] = muk;
 //PassbyrefGalV(void *params,double k1, double k2,double phi, double mu1)
 //   double tmp1 = Integrate<3>(bind(&Spectroscopic::PassbyrefGalV,this,x,fp,_1,_2,_3),a3,b3,EPSREL,EPSABS)/(8.0*M_PI*M_PI); 
    double tmp1 = Integrate<2>(bind(&SpectroscopicVector::PassbyrefGalV,this,fp,k1,k2,_1,_2),a,b,EPSREL,EPSABS)/(4.0*M_PI);
    //double tmp1 = Integrate(bind(&Spectroscopic::PassbyrefGalV,this,x,fp,phi,M_PI/2.0,_1),-0.9998,1.0,EPSREL,EPSABS)/(4.0*M_PI);
  
    return tmp1;
}



double SpectroscopicVector::ReducedGalBispectrumgVAverage(int x, double fnl,double z,double k1, double k2,double muk,double phin)
{
	LinearPowerSpectrumGR LP;

   double r = k2/k1;
   double y = sqrt(r*r+2.0*muk *r + 1.0);
   double k3 =  k1*y;
   double PowerSPec_perm12 = LP.LinearPowerAverage(3,fnl, z,k1)*LP.LinearPowerAverage(3,fnl,z,k2);
   double PowerSPec_perm32 = LP.LinearPowerAverage(3,fnl, z, k3)*LP.LinearPowerAverage(3,fnl,z,k2);
   double PowerSPec_perm13 = LP.LinearPowerAverage(3,fnl,z, k1 )*LP.LinearPowerAverage(3,fnl,z,k3);

   return GalBispectrumgVAverage(x, fnl,z,k1, k2,muk,phin)/(PowerSPec_perm12+PowerSPec_perm32+PowerSPec_perm13);

}




const int Narray = 50;

//std::vector <double> karray(Narray);
std::vector <double> karray2(Narray);
//std::vector <double> phiarray(Narray);

   //double *karray = (double*)malloc(Narray * sizeof(double));
	//double *phiarray = (double*)malloc(Narray * sizeof(double));

   double *karray = (double*)malloc(Narray* sizeof(double) );
   double *phiarray = (double*)malloc(Narray* sizeof(double) );

//std::vector <double> BVtab(Narray *Narray);
//std::vector<std::vector<double> > BVtab(Narray,std::vector<double>(Narray));

/*
std::vector <const double> karray = {5.13562e-8, 1.16879e-7, 2.65998e-7, 6.05371e-7, 
 1.37773e-6, 3.1355e-6, 
 7.13592e-6, 0.0000162403, 0.0000369603, 0.000084116, 0.000191435, 
0.000435677, 0.000991533, 0.00225658, 0.00513562, 0.0116879, 
0.0265998, 0.0605371, 0.137773, 0.31355, 0.713592, 1.62403, 3.69603, 
8.4116, 19.1435};
  double karray2[Narray];


 std::vector <const double> phiarray = {0.00012618, 0.261921, 0.523717, 0.785512, 1.04731, 1.3091, 1.5709, \
1.83269, 2.09449, 2.35628, 2.61808, 2.87987, 3.14167, 3.40346, 
3.66526, 3.92705, 4.18885, 4.45064, 4.71244, 4.97424, 5.23603, 
5.49783, 5.75962, 6.02142, 6.28321};
*/
/*
const double karray[] = {5.13562e-8, 1.16879e-7, 2.65998e-7, 6.05371e-7, 
 1.37773e-6, 3.1355e-6, 
 7.13592e-6, 0.0000162403, 0.0000369603, 0.000084116, 0.000191435, 
0.000435677, 0.000991533, 0.00225658, 0.00513562, 0.0116879, 
0.0265998, 0.0605371, 0.137773, 0.31355, 0.713592, 1.62403, 3.69603, 
8.4116, 19.1435};
  
 const double phiarray[] = {0.00012618, 0.261921, 0.523717, 0.785512, 1.04731, 1.3091, 1.5709,
1.83269, 2.09449, 2.35628, 2.61808, 2.87987, 3.14167, 3.40346, 
3.66526, 3.92705, 4.18885, 4.45064, 4.71244, 4.97424, 5.23603, 
5.49783, 5.75962, 6.02142, 6.28321};
*/

//std::vector <double> BVtab(Narray *Narray);
  //const double *BVtab[Narray][Narray];
 //double karray2[Narray];
 
    double BVtabij[Narray][Narray];
// std::array<double,25> karray = {}
   // std::vector<std::vector< double> > BVtab(Narray,std::vector<double>(Narray));
double SpectroscopicVector::PassbyrefGalV2(int x, void *params, double k1, double phi, double mu1)
{
  double *fp = (double *) params;
  double fnl = fp[0];
  double  z = fp[1];
  double k2 = fp[2];
  double phin = fp[3];
  double mu12 = fp[4];
return GalBispectrumV(x, fnl, z, k1,k2,mu12,phi,phin,mu1);
}

void SpectroscopicVector::kvalues()
{
	  const  double kmin = 0.0001;
 const  double kmax = 3.0;
 const  double phimin = 0.00001;
 const  double phimax = 2.0*M_PI;
	for(int i=0;i< Narray; i++)
	{
		karray[i] = exp(((log(kmax)-log(kmin))*((double) i)/((double) Narray-1) + log(kmin) ));
        
        phiarray[i] = phimin+ i*(phimax-phimin)/((double) Narray-1);
	}
}

void SpectroscopicVector::CalculateVB(int x, double fnl,double z,double x2,double phin,double mu12)
{
    const double EPSREL = 1e2;
    const double EPSABS = 1e2;
 
    int thread_number8 = 8;
    double fp[5];
    fp[0] = fnl;
    fp[1] = z;
    //fp[2] = x2*karray[i];
    fp[3] = phin;
    fp[4] = mu12;
 //const double *BVtabint[Narray][Narray];
//#pragma omp parallel private(thread_number)
{
 //#pragma omp for schedule(static) nowait

/*tabulate karray and phiarray*/
/*	for(int i=0;i< Narray; i++)
	{
		karray[i] = exp(((log(kmax)-log(kmin))*((double) i)/((double) Narray-1) + log(kmin) ));
        
        phiarray[i] = phimin+ i*(phimax-phimin)/Narray;
	}*/
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
for(int i = 0;i< Narray;i++)
{
      //karray[i] = exp(((log(kmax)-log(kmin))*((double) i)/((double) Narray-1) + log(kmin) ));
      //phiarray[i] = phimin+ i*(phimax-phimin)/Narray;
   for(int j = 0; j < Narray; j++)
    {
//BVtab[i][j] = Integrate(bind(&Spectroscopic::GalBispectrumV,this,x,fnl,z,karray[i],karray2[i],mu12,phiarray[j],0.0,_1),-0.998,0.998,EPSREL,EPSABS)/2.0;
//BVtabij[i][j] = integrate1(bind(&Spectroscopic::GalBispectrumV,this,x,fnl,z,karray[i],karray2[i],mu12,phiarray[j],0.0,_1),-0.998,0.998)/2.0;
//BVtabij[i][j] = GalBispectrumV(x,fnl,z,karray[i],x2*karray[i],mu12,phiarray[j],0.0,1);
 //double Spectroscopic::PassbyrefGalV2(int x, void *params, double k1, double phi, double mu1);
   fp[2] = karray[i]*x2;
   BVtabij[i][j] = Integrate(bind(&SpectroscopicVector::PassbyrefGalV2,this,x,fp,karray[i],phiarray[j],_1),-0.988,0.999,EPSREL,EPSABS)/(4.0*M_PI);

//cout<<"\t"<<"k1"<<"\t"<<karray[i]<<"\t"<<"phi"<<"\t"<<phiarray[j]<<"\t"<<"bBV"<<"\t"<<"\t"<<BVtabij[i][j]<<endl;
  


   }
}




}


cout <<"\t"<< "Finished"<<"\t"<<endl;
}


 //const double karrayx = karray;

 //const double phiarrayy = phiarrayy;
//std::vector <double> BVtab(Narray *Narray);
//std::vector<std::vector<double> > vect2(Narray,std::vector<double>(Narray));
//std::vector < double> vect2;
 //vect2 = BVtabij;

 //std::vector<double> output(karray.size());
 //std::copy(karray.begin(),karray.end(), output.begin());
/*
void func()
{
	
CalculateVB(1,0.0,zs1,1.0,cos(thetam2));
double BVtab[Narray][Narray];
memcpy(BVtab, BVtabij, sizeof(BVtabij));
//return 0;
}*/

double Passintegrand(double k, double phi)
{
   cout <<"\t"<< "Now splining"<<"\t"<<endl;

   //std::vector<std::vector< double> > vect2(Narray,std::vector< double>(Narray));
   //std::vector < double> vect2;
	// double vect2[Narray][Narray] ;
   
   //const double vect2[karray.size()*phiarray.size()] =  BVtabij;
//const double BVtab = BVtabij;
  //memcpy(BVtab, BVtabij, sizeof(BVtabij));
  

	 
	//int a[3];
  //  int b[] = {1,2,3};

	double *karray2 = (double*)malloc(Narray * sizeof(double));
	double *phiarray2 = (double*)malloc(Narray * sizeof(double));

//memcpy(&vect2, &BVtabij, sizeof vect2);

	const gsl_interp2d_type *T = gsl_interp2d_bicubic;

   gsl_spline2d *spline = gsl_spline2d_alloc(T, Narray, Narray);
   gsl_interp_accel *xacc = gsl_interp_accel_alloc();
   gsl_interp_accel *yacc = gsl_interp_accel_alloc();

   // set z grid values
   double *za = (double*)malloc(Narray*Narray * sizeof(double));
 
/*
   gsl_spline2d_set(spline, za, karray[1],phiarray[1], BVtabij[1][1]);
  gsl_spline2d_set(spline,za, karray[5],phiarray[5], BVtabij[5][5]);
  gsl_spline2d_set(spline, za,karray[10],phiarray[10], BVtabij[10][10]);
  gsl_spline2d_set(spline,za, karray[15],phiarray[15], BVtabij[15][15]);
   gsl_spline2d_set(spline,za, karray[20],phiarray[20], BVtabij[20][20]);
   gsl_spline2d_set(spline,za, karray[25],phiarray[25], BVtabij[25][25]);
   
for(int i =0; i , Narray; i++)
{
	phiarray2[i] = phiarray[i];
      karray2[i] = karray[i];
}
*/
for(int i = 0; i< Narray;i++)
{
	//phiarray2[i] = phiarray[i];
    //  karray2[i] = karray[i];
   for(int j = 0; j < Narray; j++)
    {
      gsl_spline2d_set(spline, za, karray[i],phiarray[j], BVtabij[i][j]);

      //cout<<"\t"<<"k2"<<"\t"<<karray[i]<<"\t"<<"phi2"<<"\t"<<phiarray[j]<<"\t"<<"bBV"<<"\t"<<"\t"<<BVtabij[i][j]<<endl;
  
    }
    

}



  gsl_spline2d_init(spline, karray, phiarray,za, Narray, Narray);

 double tmp =  gsl_spline2d_eval(spline,k,phi,xacc, yacc);

 return tmp;


cout <<"\t"<< "Finished splining"<<"\t"<<endl;
}

double SpectroscopicVector::GalBispectrumgVAverage2(double k)
{
    const double EPSREL = 1e-8;
    const double EPSABS = 1e-8;
    
    double tmp1 = Integrate(bind(Passintegrand,k,_1),0.00001,2.0*M_PI,EPSREL,EPSABS)/(2.0*M_PI);
  
    return tmp1;
}







//////////////////Vector bispectrum ends///////////