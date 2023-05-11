#include "IntermediateCosmology.h"
//#include "MyCosmology.h"
#include "LinearPowerSpectrumGR.h"




LinearPowerSpectrumGR::LinearPowerSpectrumGR()
{
    //SetPowerSpectrum();
}

LinearPowerSpectrumGR::~LinearPowerSpectrumGR()
{
    //SetPowerSpectrum();
}


////////////////////////Geometric tools//////////////


double LinearPowerSpectrumGR::mu2_inVT(double phi,double phin, double mu1, double muk)
{

   double mu12arg  = 1.0 - pow(muk,2);//sqrt(abs(mu12arg) < 1e-12 ? 0 : mu12arg)
  double mu1arg = 1.0 - pow(mu1,2);//
   
    double sqrtmu12 = sqrt(abs(mu12arg) );
   double sqrtmu1 = sqrt(abs(mu12arg) );

   double tmp1 = mu1* muk + sqrtmu1*sqrtmu12*cos(phin-phi);
   return tmp1;
}
double LinearPowerSpectrumGR::mu3_inVT(double k1,double k2, double k3,double phi,double phin, double mu1, double muk)
{
  double tmp1 = -(k1*mu1 + k2*mu2_inVT(phi,phin,mu1, muk))/k3;
  return tmp1;
}

/////power spectrum from camb////////////////////////////////

////Power spectrum from Class








const int N = 587; ///You have to make sure this number corresponds to the number of rows in class//cat explanatory13_pk.dat | wr -l
std::vector <double> k_array(N); // define vectors to k-values
std::vector <double> pk_array(N);  // define vectors to hold pofk


void LinearPowerSpectrumGR::ReadWrite(void)
{
  ifstream inFile; //class that reads the file
    string line; ///to take care of the header comments in .dat
  inFile.open("explanatory02_pk.dat");///load class file
 
  if(!inFile)
    {
      cerr << "could not open file" << endl; ///check whether the file is readable
    }
 while (getline(inFile, line))
    {
      if (line[0] == '#') continue;///Tell compiler to ignore comments in class's .dat file

        for(int i = 1; i < N; i++)
           {
  inFile >> k_array[i] >> pk_array[i]; /// Save the content of the file in k_array and pk_array

 //cout<<"\t"<<"k="<<"\t"<<k_array[i]<<"\t"<<""<<"pk="<<"\t"<<pk_array[i]<<endl;//check whether it is doing the right thing.
          }
    }
    inFile.close();
}

double MySplinePofk(double k)
{
  Spline<double, double> CubicSpline_num(k_array,pk_array);
  double res = CubicSpline_num.interpolate(k);
  return res;
}

double LinearPowerSpectrumGR::P1P2(double z, double  k1, double k2)
  {
    Cosmology C;
    double Dz = C.D_growth_Tabulated(z);
    //double tmp =  pow(Dz,4)* MySplinePofk(k1)*MySplinePofk(k2);
    double tmp =  pow(Dz,4)* C.P_k_EH_Tabulated(k1)*C.P_k_EH_Tabulated(k2);
   return tmp;
   }
////////////////

double LinearPowerSpectrumGR::MatterPowerSpectrum(double z, double k)
{
      Cosmology C;
     double Dz = C.D_growth_Tabulated(z);
     double Pk=  pow(Dz,2)*C.P_k_EH_Tabulated(k); 
     return Pk;
}

double Firstgamma1(double z)
{
	IntermediateCosmology TC;
	double tmp1 = TC.fg(z)*(TC.bevo(z) - 2.0*TC.curlyQ(z)- TC.Hpr(z)/pow(TC.HH(z),2) -2.0*(1-TC.curlyQ(z))/(TC.chi(z)*TC.HH(z)));
	return tmp1*TC.HH(z);
}


double Firstgamma2(double z)
{
	IntermediateCosmology TC;
	double tmp1 = -TC.fg(z)*(TC.bevo(z) -3.0);
	double tmp2A = -2.0 + TC.fg(z) - TC.bevo(z) + 4.0* TC.curlyQ(z)  + TC.Hpr(z)/pow(TC.HH(z),2)+ 2.0*(1.0-TC.curlyQ(z))/(TC.chi(z)*TC.HH(z));
	double tmp2B = -1.5*TC.Omegam(z);
	double tmp2 = tmp2A*tmp2B;
        double tmp = tmp1 + tmp2;



	return tmp*pow(TC.HH(z),2);
	
}

double FirstOrdergamma1(double z)
{
	IntermediateCosmology TC;
	double tmp1 = TC.fg(z)*(TC.bevo(z) - 2.0*TC.curlyQ(z)- TC.Hpr(z)/pow(TC.HH(z),2) -2.0*(1-TC.curlyQ(z))/(TC.chi(z)*TC.HH(z)));
	return tmp1;
}
double FirstOrdergamma2(double z)
{
	IntermediateCosmology TC;
	double tmp1 = -TC.fg(z)*(TC.bevo(z) -3.0);
	double tmp2A = -2.0 + TC.fg(z) - TC.bevo(z) + 4.0* TC.curlyQ(z)  + TC.Hpr(z)/pow(TC.HH(z),2)+2.0*(1.0-TC.curlyQ(z))/(TC.chi(z)*TC.HH(z));
	double tmp2B = -1.5*TC.Omegam(z);
	double tmp2 = tmp2A*tmp2B;
        double tmp = tmp1 + tmp2;
	return tmp;
	
}


double kerneloneN(double fnl, double z,double k, double mu)
{
	IntermediateCosmology TC;
      double b1 =  TC.bias1(1.0,z);// b_one(fnl,k, z);
	double tmp1 = b1 + TC.fg(z)*mu*mu;
	return tmp1;
}

double kerneloneGRI(double z,double k, double mu)
{
	IntermediateCosmology TC;
	double tmp1 = FirstOrdergamma1(z)*mu*TC.HH(z)/k;
	return tmp1;
}

double kerneloneGRR(double z,double k, double mu)
{
	IntermediateCosmology TC;
	double tmp1 = FirstOrdergamma2(z)*pow(TC.HH(z)/k,2);
	return tmp1;
}






double LinearPowerSpectrumGR::Linearkernel(int x, double fnl, double z,double k, double mu)
{
if(x ==1)
{
 return kerneloneN(fnl,  z, k,  mu);
}
else if (x == 2)
{
return kerneloneGRR(z, k, mu);
}
else if( x == 3)
{
return kerneloneGRI(z,k,mu);
}
else
return kerneloneN(fnl,  z, k,  mu) + kerneloneGRR(z, k, mu);

}




double LinearPowerSpectrumGR::LinearPower(int x, double fnl, double z, double k, double mu)
{
       Cosmology C;
     
     double Dz = C.D_growth_Tabulated(z);
     double power_multi =  pow(Dz,2)*C.P_k_EH_Tabulated(k); 
if(x == 1)
{
double tmp =  kerneloneN(fnl, z,k,mu)*kerneloneN(fnl, z,k,mu)*power_multi;
return tmp;
}
else if(x == 2)
{
double tmp = pow(kerneloneGRR(z, k, mu),2) + pow(kerneloneGRR(z,k,mu),2);
return tmp*power_multi;
}
else
{
double tmp  = kerneloneN(fnl, z,k,mu)*kerneloneN(fnl, z,k,mu) + 2.0*kerneloneN(fnl, z,k,mu)*kerneloneGRR(z, k, mu) + pow(kerneloneGRR(z, k, mu),2) + pow(kerneloneGRI(z,k,mu),2);   
     return  tmp*power_multi;
}

}


double LinearPowerSpectrumGR::LinearPowerAverage(int x,double fnl, double z,double k)
{
 const double EPSREL = 1e-2;

  const double EPSABS = 1e-2;
 double mu_min = -0.9998;
 double mu_max = 0.9998;//LinearPower_HIHI(int x, double fnl, double z, double k, double mu)
double tmp = Integrate(bind(&LinearPowerSpectrumGR::LinearPower,this,x,fnl,z,k,_1),mu_min,mu_max,EPSREL,EPSABS)/2.0;
return tmp;
}
