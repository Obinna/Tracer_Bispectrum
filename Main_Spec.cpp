
#include<iostream>
#include <fstream>

#include <time.h>
#include <sys/time.h>

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<cmath>
#include <list>



#include <omp.h>
#include <iomanip>
#include <pthread.h>
#include <list>



#include"/home/obinna/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/rk_int.cpp"
#include"/home/obinna/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/routines.cpp"
#include"/home/obinna/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/constants.cpp"
#include"/home/obinna/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/MyCosmology.cpp"
//#include <array>



/*
#include"/Users/obinna/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/rk_int.cpp"
#include"/Users/obinna/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/routines.cpp"
#include"/Users/obinna/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/constants.cpp"
#include"/Users/obinna/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/MyCosmology.cpp"

*/
/*
#include"/Users/o-admin/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/rk_int.cpp"
#include"/Users/o-admin/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/routines.cpp"
#include"/Users/o-admin/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/constants.cpp"
#include"/Users/o-admin/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/MyCosmology.cpp"
*/


//#include"MyHIBias.cpp"
//#include"ToB_Bias.cpp"

#include "GrowthEquation.cpp"
#include "SpectroscopicGal.cpp"





using namespace std;









double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}
//#endif









int main(int argc, char **argv)
{

//const char* output = "DataSpec.dat";

//const char* output = "DataSpec1.dat";

//const char* output = "DataSpec2.dat";

//const char* output = "DataSpec3.dat";

//const char* output = "DataSpec4.dat";

//const char* output = "DataSpec5.dat";
//const char* output = "DataSpec6.dat";


//const char* output = "DataSpec7.dat";

//const char* output = "DataSpec8.dat";

//const char* output = "DataSpec9.dat";

//const char* output = "DataSpec10.dat";

//const char* output = "DataSpec11.dat";

//const char* output = "DataSpec12.dat";

//const char* output = "DataSpec13.dat";

  const char* output = "DataSpec14.dat";



ofstream myfile;
myfile.setf(ios::scientific);
myfile.open(output);


  Cosmology C;
  
/*
//  Fisher_Deriv  FD;
 C.InitializeCosmology();
//C.InitialzieSpeed2();
 C.InitializeDz();
 C.TabulateSigma_EH();

*/
C.InitializeCosmology();
C.InitializeDz();
C.InitializeSigma_EH();
C.TabulateSigma_EH();

//GrowthEquation GE;
//GE.InitialzeFloop();
//cout <<"\t"<<"F"<<"\t"<<GE.CurlyF(0.5)<<"\t"<<"G"<<"\t"<<GE.CurlyG(0.5)<<endl;



//ToBias  ToB;
//HIBias  HIB;

//HIB.InitializeAverageloop();

Spectroscopic S;
S.InitialzieSpeed4();
S.InitialzieGamma();
//S.InitialzieCurlyterms();
///S.InitialzieSpeedMatter();




  double M_min = 10e8;
  double M_max = 10e13;
  double M;
  double z;
  int nM = 1000;
  double za[5] = {10.07,5.72,3.06,1.50,0.0};
  double pa[5];
int N = 100;
 double fnlmin = -1.0*5;
 double fnlmax = 5;
vector <double > z_array(N);
vector <double > k_array(N);
vector <double > k_array2(N);
vector <double > M_array(N);
vector <double > R_array(N);
vector <double > mu1_array(N);
vector <double > mu2_array(N);
vector <double > fnl_array(N);
vector <long double> Out1(N);
vector <long double> Out2(N);
vector <long double> Out3(N);
vector <long double> Out4(N);
vector <long double> Out5(N);
vector <long double> Out6(N);
vector <long double> Out7(N);
vector <long double> Out8(N);
vector <long double> Out9(N);
vector <long double> Out10(N);
vector <long double> Out11(N);

vector <long double> Out12(N);
vector <long double> Out13(N);

vector <double> theta_array(N);
//vector <double> Out12(N);
vector <int> ell(N);


vector<vector<double> > muk_Out(N,vector<double>(N));
vector<vector<double> >  BK_Out1(N,vector<double>(N));
vector<vector<double> >  BK_Out2(N,vector<double>(N));
vector<vector<double> >  BK_Out3(N,vector<double>(N));
vector<vector<double> >  BK_Out4(N,vector<double>(N));
vector<vector<double> >  BK_Out5(N,vector<double>(N));
vector<vector<double> >  BK_Out6(N,vector<double>(N));
vector<vector<double> >  BK_Out7(N,vector<double>(N));


/*
double muk_Out[N][N];
double BK_Out1[N][N];
double BK_Out2[N][N];
double BK_Out3[N][N];
double BK_Out4[N][N];
double BK_Out5[N][N];
double BK_Out6[N][N];
*/
const  double kmin = 0.0001;
 const  double kmax = 1.0;


const  double kmin2 = 1.0;
 const  double kmax2 = 2.0;

double zmin = 0.02;
  double zmax = 5.0;


//Squeesed limit
const double k_1 =0.1;
const double k_2 =0.001;
const double k_3 =0.001;
const double k_4 =0.01;
const double k_5 =0.1;
const double x_half = 0.5;
const double x1 = 1.0;
const double x2 = 2.0;


double mumin = -0.9991;
  double mumax = 0.999;

double thetamin = 0.001;
  double thetamax = M_PI-0.0001;

const  double zs = 0.5;
const double zs1 = 1.0;
const double zs2 = 1.5;
 const double zs3 = 2.0;
const double zs4 = 2.5;


 double muin = 0.0;

///S.InitialzieSpeedM(zs);

const double thetam = 1.0*M_PI/180.0;
const double thetam2 = 177.0*M_PI/180.0;


double fnlsingle =0.0;
double fnlsingle1 =0.2;
double fnlsingle2 =0.0;
double fnlsingle3 =2.0;
double fnlsingle4 =5.0;

double smallh = 0.68;

/*
cout <<"\t"<<"K_H= zs"<<"\t"<<L.k_H(zs)<<"\t"<<"k_BAO= zs"<<"\t"<<L.k_BAO(zs)<<endl;
cout <<"\t"<<"K_H= zs1"<<"\t"<<L.k_H(zs1)<<"\t"<<"k_BAO= zs1"<<"\t"<<L.k_BAO(zs1)<<endl;
cout <<"\t"<<"K_H= zs2"<<"\t"<<L.k_H(zs2)<<"\t"<<"k_BAO= zs2"<<"\t"<<L.k_BAO(zs2)<<endl;
cout <<"\t"<<"K_H= zs3"<<"\t"<<L.k_H(zs3)<<"\t"<<"k_BAO= zs3"<<"\t"<<L.k_BAO(zs3)<<endl;
cout <<"\t"<<"K_H= zs4"<<"\t"<<L.k_H(zs4)<<"\t"<<"k_BAO= zs4"<<"\t"<<L.k_BAO(zs4)<<endl;
*/

int thread_number8 = 4;
#pragma omp parallel private(thread_number8)
{
#pragma omp for schedule(static) nowait
for(int i =0; i < N;i++)
{ 
z_array[i] = exp( (  (log(zmax)-log(zmin))*((double) i)/((double) N-1) + log(zmin) ));//zmin + i*(zmax-zmin)/N;// 
theta_array[i] = (thetamin+i*(thetamax-thetamin)/N);//exp( (  (log(thetamax)-log(thetamin))*((double) i)/((double) N-1) + log(thetamin) ));////;//
k_array[i] = exp( ((log(kmax)-log(kmin))*((double) i)/((double) N-1) + log(kmin) ));
k_array2[i] = k_array[i]/x2;//exp( ((log(kmax2)-log( kmin2))*((double) i)/((double) N-1) + log( kmin2) ));//kmin2 + i*(kmax2-kmin2)/N;//

// C.P_k_EH_Tabulated(k_array[i]);


/*
Out1[i] = S.ReducedRealSpaceBispectrum(1,0.0,zs,k_1,k_1/x1, cos(theta_array[i]));
Out2[i] = S.ReducedRealSpaceBispectrum(2,0.0,zs,k_1,k_1/x1, cos(theta_array[i]));
Out3[i] = S.ReducedRealSpaceBispectrum(3,0.0,zs,k_1,k_1/x1, cos(theta_array[i]));
Out4[i] = S.ReducedRealSpaceBispectrum(4,0.0,zs,k_1,k_1/x1, cos(theta_array[i]));
Out5[i] = S.ReducedRealSpaceBispectrum(5,0.0,zs,k_1,k_1/x1, cos(theta_array[i]));

 myfile<<theta_array[i]/M_PI<<" "<<Out1[i]<<" "<<Out2[i]<<" "<<Out3[i]<<" "<<Out4[i]<<" "<<Out5[i]<<endl;
  cout <<"i="<<i<<"\t"<<k_array2[i]<<"\t"<<" "<<theta_array[i]/M_PI<<"\t"<<" "<<Out1[i]<<"\t"<<" "<<Out2[i]<<"\t"<<" "<<Out3[i]<<"\t"<<" "<<Out4[i]<<"\t"<<" "<<Out5[i]<<endl;
*/

/*
Out1[i] = S.ReducedRealSpaceBispectrum(4,0.0,zs,k_1,k_1/x1, cos(theta_array[i]));
Out2[i] = S.ReducedRealSpaceBispectrum(5,0.0,zs,k_1,k_1/x1, cos(theta_array[i]));
Out3[i] = S.ReducedRealSpaceBispectrumfnl(4,fnlsingle3,zs,k_1,k_1/x1, cos(theta_array[i]));
Out4[i] = S.ReducedRealSpaceBispectrumfnl(5,fnlsingle3,zs,k_1,k_1/x1, cos(theta_array[i]));
Out5[i] = S.ReducedRealSpaceBispectrumfnl(6,fnlsingle3,zs,k_1,k_1/x1, cos(theta_array[i]));

 myfile<<theta_array[i]/M_PI<<" "<<Out1[i]<<" "<<Out2[i]<<" "<<Out3[i]<<" "<<Out4[i]<<" "<<Out5[i]<<endl;
  cout <<"i="<<i<<"\t"<<k_array2[i]<<"\t"<<" "<<theta_array[i]/M_PI<<"\t"<<" "<<Out1[i]<<"\t"<<" "<<Out2[i]<<"\t"<<" "<<Out3[i]<<"\t"<<" "<<Out4[i]<<"\t"<<" "<<Out5[i]<<endl;
*/


//Out1[i] = C.P_k_EH_Tabulated(k_array[i]);
//Out2[i] = S.LinearPowerAverage(1,0.0, 0.5,k_array[i]);//ReducedNewtonianBispectrumWithRSD(int x, double fnl,double z,double k1, double k2,double muk)

/*
Out1[i] =S.ReducedRealSpaceBispectrum(5,0.0,0.5,k_1,k_1/x1, cos(theta_array[i]));
Out2[i] = S.ReducedNewtonianBispectrumWithRSD(8,0.0,0.5,k_1,k_1/x1, cos(theta_array[i]));
Out3[i] =S.ReducedRealSpaceBispectrum(5,0.0,0.05,k_1,k_1/x2, cos(theta_array[i]));
Out4[i] = S.ReducedNewtonianBispectrumWithRSD(8,0.0,0.05,k_1,k_1/x2, cos(theta_array[i]));
//Out5[i] = S.ReducedGalBispectrumgAverage(8,0.0,0.5,k_1,k_1/2.0,cos(theta_array[i]));
Out5[i] =S.ReducedRealSpaceBispectrum(5,0.0,0.001,k_1,k_1/x1, cos(theta_array[i]));
Out6[i] = S.ReducedNewtonianBispectrumWithRSD(8,0.0,0.001,k_1,k_1/x1, cos(theta_array[i]));

 
 myfile<<theta_array[i]/M_PI<<" "<<Out1[i]<<" "<<Out2[i]<<" "<<Out3[i]<<" "<<Out4[i]<<" "<<Out5[i]<<" "<<Out6[i]<<endl;
  cout <<"i="<<i<<"\t"<<k_array2[i]<<"\t"<<" "<<theta_array[i]/M_PI<<"\t"<<" "<<Out1[i]<<"\t"<<" "<<Out2[i]<<"\t"<<" "<<Out3[i]<<"\t"<<" "<<Out4[i]<<"\t"<<" "<<Out5[i]<<"\t"<<" "<<Out6[i]<<endl;*/

/*
Out1[i] = S.NewtonianBispectrumWithRSD(8,0.0,0.5,k_array[i],k_array[i],cos(thetam2),0.0, 1);
Out2[i] = S.NewtonianBispectrumWithRSD(8,0.0,0.5,k_array[i],k_array[i],cos(thetam2),0.0, 0.0);
Out3[i] = S.GalBispectrumg(8,0.0,0.5,k_array[i],k_array[i],cos(thetam2),0.0, 1);
Out4[i] = S.GalBispectrumg(8,0.0,0.5,k_array[i],k_array[i],cos(thetam2),0.0, 0.0);
Out5[i] = S.GalBispectrumgAverage(8,0.0,0.5,k_array[i],k_array[i],cos(thetam2));


myfile<<k_array[i]<<" "<<abs(Out1[i])<<" "<<abs(Out2[i])<<" "<<abs(Out3[i])<<" "<<abs(Out4[i])<<" "<<abs(Out5[i])<<"\n";
cout<<"\t" <<"i="<<i<<" "<<"k="<<k_array[i]<<", "<<Out1[i]<<", "<<Out2[i]<<", "<<Out3[i]<<", "<<Out4[i]<<", "<<Out5[i]<<endl;

*/

//AveragedNewtonianBispectrumWithRSD(int x, double fnl,double z,double k1, double k2,double muk)

/*
Out1[i] = S.AveragedNewtonianBispectrumWithRSD(8,0.0,zs1,k_array[i],k_array[i],cos(thetam2));
Out2[i] = S.AveragedNewtonianBispectrumWithRSD(9,0.0,zs1,k_array[i],k_array[i],cos(thetam2));

Out3[i] = S.AveragedNewtonianBispectrumWithRSD(8,0.0,zs,k_array[i],k_array[i],cos(thetam2));
Out4[i] = S.AveragedNewtonianBispectrumWithRSD(9,0.0,zs,k_array[i],k_array[i],cos(thetam2));

Out5[i] = S.GalBispectrumgAverage(8,0.0,zs,k_array[i],k_array[i],cos(thetam2));
Out6[i] = S.GalBispectrumgAverage(8,0.0,zs1,k_array[i],k_array[i],cos(thetam2));


myfile<<k_array[i]<<" "<<abs(Out1[i])<<" "<<abs(Out2[i])<<" "<<abs(Out3[i])<<" "<<abs(Out4[i])<<" "<<abs(Out5[i])<<" "<<abs(Out6[i])<<"\n";
cout<<"\t" <<"i="<<i<<" "<<"k="<<k_array[i]<<", "<<Out1[i]<<", "<<Out2[i]<<", "<<Out3[i]<<", "<<Out4[i]<<", "<<Out4[i]<<", "<<Out5[i]<<", "<<Out6[i]<<endl;
*/


/*
Out1[i] = S.AveragedNewtonianBispectrumWithRSD(8,0.0,zs,k_array[i],k_array[i],cos(thetam2));
Out2[i] = S.AveragedNewtonianBispectrumWithRSD(8,0.0,zs1,k_array[i],k_array[i],cos(thetam2));
Out3[i] = S.AveragedNewtonianBispectrumWithRSD(8,0.0,zs2,k_array[i],k_array[i],cos(thetam2));
Out4[i] = S.AveragedNewtonianBispectrumWithRSD(8,0.0,zs3,k_array[i],k_array[i],cos(thetam2));
Out5[i] = S.AveragedNewtonianBispectrumWithRSD(8,0.0,zs4,k_array[i],k_array[i],cos(thetam2));


myfile<<k_array[i]<<" "<<abs(Out1[i])<<" "<<abs(Out2[i])<<" "<<abs(Out3[i])<<" "<<abs(Out4[i])<<" "<<abs(Out5[i])<<"\n";
cout<<"\t" <<"i="<<i<<" "<<"k="<<k_array[i]<<", "<<Out1[i]<<", "<<Out2[i]<<", "<<Out3[i]<<", "<<Out4[i]<<", "<<Out4[i]<<", "<<Out5[i]<<", "<<Out6[i]<<endl;
*/




/*
Out1[i] = S.AveragedNewtonianBispectrumWithRSD(8,0.0,zs,k_array[i],k_array[i],cos(thetam2));
Out2[i] = S.AveragedNewtonianBispectrumWithRSD(8,0.0,zs1,k_array[i],k_array[i],cos(thetam2));
Out3[i] = S.AveragedNewtonianBispectrumWithRSD(8,0.0,zs2,k_array[i],k_array[i],cos(thetam2));
Out4[i] = S.AveragedNewtonianBispectrumWithRSD(8,0.0,zs3,k_array[i],k_array[i],cos(thetam2));
Out5[i] = S.AveragedNewtonianBispectrumWithRSD(8,0.0,zs4,k_array[i],k_array[i],cos(thetam2));


myfile<<k_array[i]<<" "<<abs(Out1[i])<<" "<<abs(Out2[i])<<" "<<abs(Out3[i])<<" "<<abs(Out4[i])<<" "<<abs(Out5[i])<<"\n";
cout<<"\t" <<"i="<<i<<" "<<"k="<<k_array[i]<<", "<<Out1[i]<<", "<<Out2[i]<<", "<<Out3[i]<<", "<<Out4[i]<<", "<<Out4[i]<<", "<<Out5[i]<<endl;
*/

/*
Out1[i] = S.AveragedNewtonianBispectrumWithRSD(8,0.0,zs,k_array[i],k_array[i],cos(thetam2));
Out2[i] = S.AveragedNewtonianBispectrumWithRSD(8,0.0,zs1,k_array[i],k_array[i],cos(thetam2));
Out3[i] =  S.SquzzedbispectrumN(1.0,k_array[i], zs);
Out4[i] =  S.SquzzedbispectrumN(1.0,k_array[i], zs1);

myfile<<k_array[i]<<" "<<abs(Out1[i])<<" "<<abs(Out2[i])<<" "<<abs(Out3[i])<<" "<<abs(Out4[i])<<"\n";
cout<<"\t" <<"i="<<i<<" "<<"k="<<k_array[i]<<", "<<Out1[i]<<", "<<Out2[i]<<", "<<Out3[i]<<", "<<Out4[i]<<", "<<Out5[i]<<", "<<Out6[i]<<endl;

*/


Out1[i] =  S.Squzzedbispectrum(1.0,k_array[i], zs);
Out2[i] =  S.Squzzedbispectrum(1.0,k_array[i], zs1);
Out3[i] = S.GalBispectrumgAverage(1,0.0,zs,k_array[i],k_array[i],cos(thetam2));
Out4[i] = S.GalBispectrumgAverage(1,0.0,zs1,k_array[i],k_array[i],cos(thetam2));

myfile<<k_array[i]<<" "<<abs(Out1[i])<<" "<<abs(Out2[i])<<" "<<abs(Out3[i])<<" "<<abs(Out4[i])<<"\n";
cout<<"\t" <<"i="<<i<<" "<<"k="<<k_array[i]<<", "<<Out1[i]<<", "<<Out2[i]<<", "<<Out3[i]<<", "<<Out4[i]<<", "<<Out5[i]<<", "<<Out6[i]<<endl;

/*
Out1[i] =  S.Squzzedbispectrum(1.0,k_array[i], zs);
Out2[i] =  S.Squzzedbispectrum(1.0,k_array[i], zs1);
Out3[i] =  S.Squzzedbispectrum(1.0,k_array[i], zs2);
Out4[i] =  S.Squzzedbispectrum(1.0,k_array[i], zs3);
Out5[i] =  S.Squzzedbispectrum(1.0,k_array[i], zs4);

Out6[i] = S.GalBispectrumgAverage(1,0.0,zs,k_array[i],k_array[i],cos(thetam2));
Out7[i] = S.GalBispectrumgAverage(1,0.0,zs1,k_array[i],k_array[i],cos(thetam2));
Out8[i] = S.GalBispectrumgAverage(1,0.0,zs2,k_array[i],k_array[i],cos(thetam2));
Out9[i] = S.GalBispectrumgAverage(1,0.0,zs3,k_array[i],k_array[i],cos(thetam2));
Out10[i] = S.GalBispectrumgAverage(1,0.0,zs4,k_array[i],k_array[i],cos(thetam2));

myfile<<k_array[i]<<" "<<abs(Out1[i])<<" "<<abs(Out2[i])<<" "<<abs(Out3[i])<<" "<<abs(Out4[i])<<" "<<abs(Out5[i])<<" "<<abs(Out6[i])<<" "<<abs(Out7[i])<<" "<<abs(Out8[i])<<" "<<abs(Out9[i])<<" "<<abs(Out10[i])<<"\n";
cout <<"i="<<i<<" "<<"k="<<k_array[i]<<", "<<Out1[i]<<", "<<Out2[i]<<", "<<Out3[i]<<", "<<Out4[i]<<", "<<Out5[i]<<", "<<Out6[i]<<", "<<Out7[i]<<", "<<Out8[i]<<", "<<Out9[i]<<", "<<Out10[i]<<endl;
*/

/*
Out1[i] = S.GalBispectrumgAverage(1,0.0,zs,k_array[i],k_array[i],cos(thetam2));
Out2[i] = S.GalBispectrumgAverage(1,0.0,zs1,k_array[i],k_array[i],cos(thetam2));
Out3[i] = S.GalBispectrumgAverage(1,0.0,zs2,k_array[i],k_array[i],cos(thetam2));
Out4[i] = S.GalBispectrumgAverage(1,0.0,zs3,k_array[i],k_array[i],cos(thetam2));
Out5[i] = S.GalBispectrumgAverage(1,0.0,zs4,k_array[i],k_array[i],cos(thetam2));


myfile<<k_array[i]<<" "<<abs(Out1[i])<<" "<<abs(Out2[i])<<" "<<abs(Out3[i])<<" "<<abs(Out4[i])<<" "<<abs(Out5[i])<<"\n";
cout<<"\t" <<"i="<<i<<" "<<"k="<<k_array[i]<<", "<<Out1[i]<<", "<<Out2[i]<<", "<<Out3[i]<<", "<<Out4[i]<<", "<<Out4[i]<<", "<<Out5[i]<<endl;
*/


/*
Out1[i] = S.GalBispectrumgAverage(1,0.0,zs,k_array[i],k_array[i],cos(thetam2));
Out2[i] = S.GalBispectrumgAverage(2,0.0,zs,k_array[i],k_array[i],cos(thetam2));
Out3[i] = S.GalBispectrumgAverage(1,0.0,zs1,k_array[i],k_array[i],cos(thetam2));
Out4[i] = S.GalBispectrumgAverage(2,0.0,zs1,k_array[i],k_array[i],cos(thetam2));
Out5[i] = S.GalBispectrumgAverage(1,0.0,zs2,k_array[i],k_array[i],cos(thetam2));
Out6[i] = S.GalBispectrumgAverage(2,0.0,zs2,k_array[i],k_array[i],cos(thetam2));


myfile<<k_array[i]<<" "<<abs(Out1[i])<<" "<<abs(Out2[i])<<" "<<abs(Out3[i])<<" "<<abs(Out4[i])<<" "<<abs(Out5[i])<<" "<<abs(Out6[i])<<"\n";
cout<<"\t" <<"i="<<i<<" "<<"k="<<k_array[i]<<", "<<Out1[i]<<", "<<Out2[i]<<", "<<Out3[i]<<", "<<Out4[i]<<", "<<Out4[i]<<", "<<Out5[i]<<", "<<Out6[i]<<endl;
*/


/*
Out1[i] = S.GalBispectrumg(1,0.0,zs1,k_array[i],k_array[i],cos(thetam2),0.0, 1);
Out2[i] = S.GalBispectrumg(2,0.0,zs1,k_array[i],k_array[i],cos(thetam2),0.0, 1);
Out3[i] = S.GalBispectrumg(1,0.0,zs1,k_array[i],k_array[i],cos(thetam2),0.0, 0.0);
Out4[i] = S.GalBispectrumg(2,0.0,zs1,k_array[i],k_array[i],cos(thetam2),0.0, 0.0);
Out5[i] = S.GalBispectrumgAverage(1,0.0,zs1,k_array[i],k_array[i],cos(thetam2));
Out6[i] = S.GalBispectrumgAverage(2,0.0,zs1,k_array[i],k_array[i],cos(thetam2));


myfile<<k_array[i]<<" "<<abs(Out1[i])<<" "<<abs(Out2[i])<<" "<<abs(Out3[i])<<" "<<abs(Out4[i])<<" "<<abs(Out5[i])<<" "<<abs(Out6[i])<<"\n";
cout<<"\t" <<"i="<<i<<" "<<"k="<<k_array[i]<<", "<<Out1[i]<<", "<<Out2[i]<<", "<<Out3[i]<<", "<<Out4[i]<<", "<<Out4[i]<<", "<<Out5[i]<<", "<<Out6[i]<<endl;

*/
}
}



myfile.close();






  return 0;
}

