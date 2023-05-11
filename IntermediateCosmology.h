
#ifndef  _IntermediateCosmology
#define  _IntermediateCosmology



class IntermediateCosmology
{
    public:

double HH(double z);
double Alpha( double k,double z);
double AlphaAlpha2(double z, double  k1, double k2);
double bevo(double z);
double bevopr(double z);
double curlyQ(double z);
double curlyQpr(double z);
double curlyQzpr(double z);
double curlyQLpr(double  z);
double fg(double z);
double chi(double z);
double Omegam(double z);
double Hpr(double z);
double Omegampr(double z);
double Hpr2(double z);

double fgpr(double z);
double bias1(double fnl, double z);
double TidalBias(double z);
double bias1pr(double fnl, double z);
double bias2(double fnl, double z);
double b_one(double fnl, double k, double z);
double b_two_integrand(double fnl, double z ,double k1,double k2,double mu1);
double bias1Lpr(double fnl, double z);
double bevoLpr(double z);
double mu23(double r, double mu12);
double mu13(double r, double mu12);



    IntermediateCosmology();
    void InitialzieSpeed4(void);
    void ReadWrite(void);
   ~IntermediateCosmology();


};

double HH_in(double z);

#endif