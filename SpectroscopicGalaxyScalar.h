

#ifndef  _Spectroscopic
#define  _Spectroscopic



class Spectroscopic
{
    public:





double  MatterBispectrum(double z,double k1, double r, double mu12);
double ReducedMatterBispectrum_Q(double z,double k1, double r, double mu12);

//double LinearPowerSpectrum(int x, double fnl, double z, double k, double mu);

double RealSpaceBispectrum(int x, double fnl,double z,double k1, double k2, double mu12);
double ReducedRealSpaceBispectrum(int x,double fnl,double z,double k1, double k2, double mu12);
double NewtonianBispectrumWithRSD(int x,double fnl,double z,double k1, double k2, double mu12,double phi, double mu1);


//double LinearPower(int x, double fnl, double z, double k, double mu);
//double LinearPowerAverage(int x,double fnl, double z,double k);
double Passbyref(int x,int l , int m, void *params ,double phi, double mu1);
double AveragedNewtonianBispectrumWithRSD(int x,int l , int m, double fnl,double z,double k1, double k2,double muk);
double ReducedNewtonianBispectrumWithRSD(int x,int l , int m, double fnl,double z,double k1, double k2,double muk);

double GalBispectrumg(int w, int x,double fnl,double z,double k1, double k2, double mu12,double phi, double mu1);

double GalBispectrumabs(int x, double fnl,double z,double k1, double k2, double mu12,double phi, double mu1);
double GalBispectrumabsQ(double fnl,double z,double k1, double k2, double mu12,double phi, double mu1);
double PassbyrefGal(int w, int x, void *params ,double phi, double mu1);

double GalBispectrumgAverage(int w, int x, double fnl,double z,double k1, double k2,double muk);
double GalBispectrumgVAverage2(double k);
double ReducedGalBispectrumgAverage(int w, int x, double fnl,double z,double k1, double k2,double muk);
  
double  RealSpaceBispectrumfnl(int x,double fnl,double z,double k1, double r, double mu12);
double ReducedRealSpaceBispectrumfnl(int x,double fnl,double z,double k1, double r, double mu12);

//double SquzzedbispectrumN(double b0,double kS, double z);
//double Squzzedbispectrum(double b0,double kS, double z);

double GalBispectrumDopplerleading(int x,double fnl,double z,double k1, double k2, double mu12,double phi, double mu1);
double PassbyrefGalDopplerleading(void *params,int x, int m, double phi, double mu1);
double GalBispectrumgDipoleleading(int x, int m, double fnl,double z,double k1, double k2,double muk);


double PassbyrefDipole(void *params , int w, double phi, double mu1);
double GalBispectrumgDipole(int w, double fnl,double z,double k1, double k2,double muk);
double Testzterms(int x,double z);






double GalBispectrumN(int w,int x,double fnl,double z,double k1, double k2, double mu12,double phi, double mu1);
double GalBispectrumNabs(int x,double fnl,double z,double k1, double k2, double mu12,double phi, double mu1);
double PassbyrefGalN(int w, int x, void *params ,double phi, double mu1);
double GalBispectrumNAverage(int w, int x, double fnl,double z,double k1, double k2,double muk);
double ReducedGalBispectrumNAverage(int w, int x, double fnl,double z,double k1, double k2,double muk);


double GalBispectrumTotal(int x,double fnl,double z,double k1, double k2, double mu12,double phin, double phi, double mu1);
double PassbyrefGalTotal(int x, void *params, double phin, double phi, double mu1);
double GalBispectrumgTotalAverage(int x, double fnl,double z,double k1, double k2,double muk,double phin);
double ReducedGalBispectrumTotalAverage(int x, double fnl,double z,double k1, double k2,double muk,double phin);
double PassbyrefGalV2(int x, void *params, double k1, double phi, double mu1);



double PassbyrefGalMulti(void *Iparams, void *params ,double phi, double mu1);
double GalBispectrumgMulti(int w,int l, int m, double fnl,double z,double k1, double k2,double muk);
double ReducedGalBispectrumgMulti(int w,int l, int m, double fnl,double z,double k1, double k2,double muk);


//////////////////////////////////Different basis///////////////////
double DiffBasisNewtonianbisepctrumintegrand(void *params,double Phi,double mu12,double mu1);
double DiffBasisImaginaryGRbisepctrumintegrand3D(void *params,double Phi,double mu12,double mu1);
double GalBispectrumgDipoleBasis(int xx, double fnl,double z,double k1, double k2,double  l1,double l2, double L);

void kvalues();
void CalculateVB(int x, double fnl,double z,double x2,double phin,double mu12);
void InitialzieSpeed4(void);
void InitialzieGamma(void);




    Spectroscopic();
 // void SetPowerSpectrum(void);
    void ReadWrite(void);
   ~Spectroscopic();










};

double CurlyIGR2(double fnl, double z, double k1, double k2, double k3,double muk, double mu1, double mu2,double mu3);
double CurlyRGR2(double fnl, double z, double k1, double k2, double k3,double muk, double mu2, double mu1,double mu3);

//double bias1(double fnl, double z);
//double bias2(double fnl, double z);
//double bias3(double fnl, double z);


double mu2_in(double phi, double mu1, double muk);
double mu3_in(double k1,double k2, double k3,double phi, double mu1, double muk);

//double P1P2(double z, double  k1, double k2);


#endif
