
#ifndef  _SpectroscopicTotal
#define  _SpectroscopicTotal



class SpectroscopicTotal
{
    public:


double GalBispectrumTotal(int x,double fnl,double z,double k1, double k2, double mu12,double phin, double phi, double mu1);
double GalBispectrumgTotalAverage(int x, double fnl,double z,double k1, double k2,double muk,double phin);
double ReducedGalBispectrumTotalAverage(int x, double fnl,double z,double k1, double k2,double muk,double phin);




    SpectroscopicTotal();
 // void SetPowerSpectrum(void);
    
   ~SpectroscopicTotal();










};

#endif
