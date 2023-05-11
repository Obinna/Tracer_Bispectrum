
#ifndef  _SpectroscopicTensor
#define  _SpectroscopicTensor



class SpectroscopicTensor
{
    public:


double GalBispectrumT(int x,double fnl,double z,double k1, double k2, double mu12,double phin, double phi, double mu1);
double PassbyrefGalT(int x, void *params, double phin, double phi, double mu1);
double GalBispectrumgTAverage(int x, double fnl,double z,double k1, double k2,double muk,double phin);
double ReducedGalBispectrumgTAverage(int x, double fnl,double z,double k1, double k2,double muk,double phin);


  SpectroscopicTensor();
 
   ~SpectroscopicTensor();


};



#endif