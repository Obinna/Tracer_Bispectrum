
#ifndef  _SpectroscopicVector
#define  _SpectroscopicVector



class SpectroscopicVector
{
    public:


double GalBispectrumV(int x,double fnl,double z,double k1, double k2, double mu12,double phin, double phi, double mu1);
//double PassbyrefGalV(int x, void *params, double phin, double phi, double mu1);
double PassbyrefGalV(void *params,double k1, double k2,double phi, double mu1);
double GalBispectrumgVAverage(int x, double fnl,double z,double k1, double k2,double muk,double phin);


double PassbyrefGalV2(int x, void *params, double k1, double phi, double mu1);
void kvalues();
void CalculateVB(int x, double fnl,double z,double x2,double phin,double mu12);
double GalBispectrumgVAverage2(double k);
double ReducedGalBispectrumgVAverage(int x, double fnl,double z,double k1, double k2,double muk,double phin);



  SpectroscopicVector();
 
   ~SpectroscopicVector();



};



#endif