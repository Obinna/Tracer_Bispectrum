#ifndef  _LinearPowerSpectrumGR
#define  _LinearPowerSpectrumGR



class LinearPowerSpectrumGR
{
    public:
double MatterPowerSpectrum(double z, double k);
double Linearkernel(int x, double fnl, double z,double k, double mu);

double P1P2(double z, double  k1, double k2);

double LinearPowerAverage(int x,double fnl, double z,double k);
double LinearPower(int x, double fnl, double z, double k, double mu);
double mu2_inVT(double phi,double phin, double mu1, double muk);
double mu3_inVT(double k1,double k2, double k3,double phi,double phin, double mu1, double muk);


  LinearPowerSpectrumGR();
    void ReadWrite(void);
   ~LinearPowerSpectrumGR();

    	};

    	double kerneloneN(double fnl, double z,double k, double mu);
double kerneloneGRI(double z,double k, double mu);
double kerneloneGRR(double z,double k, double mu);



#endif