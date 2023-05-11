#include "SpectroscopicGalaxyScalar.h"

//#include "GrowthEquation.h"
//#include "MyBessel.cpp"

#include <boost/bind.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>


#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
//#include "csimpson.cxx"

#include "/home/obinna/Downloads/wignerSymbols-master/include/wignerSymbols/wignerSymbols-cpp.h"
#include "/home/obinna/Downloads/wignerSymbols-master/src/wignerSymbols-cpp.cpp"


Spectroscopic::Spectroscopic()
{
    //SetPowerSpectrum();
}

Spectroscopic::~Spectroscopic()
{
    //SetPowerSpectrum();
}



using namespace boost::math;


///////////////////////////////////Gamma terms////////////////////////////////////////////////



double Beta1_in(double fnl,double z)
{
	IntermediateCosmology IC;
  
  double tmp1 =  7.0 + 4.0* IC.bevo(z) + 4.0 *IC.curlyQ(z) - 4.0* IC.fg(z)* IC.bevo(z)+ 8.0 *IC.fg(z)*IC.curlyQ(z) - 8.0 *IC.bevo(z)*IC.curlyQ(z) + pow(IC.bevo(z),2) + 16.0 *pow(IC.curlyQ(z),2)  - 16.0*IC.curlyQLpr(z)  -2.0*IC.fgpr(z)/IC.HH(z) + IC.bevopr(z)/IC.HH(z)- 8.0*IC.curlyQpr(z)/IC.HH(z);

  double tmp2 = (IC.Hpr(z)/pow(IC.HH(z),2))*(4.0*IC.fg(z) - 6.0 - 2.0*IC.bevo(z) + 8.0*IC.curlyQ(z) + 3.0*IC.Hpr(z)/pow(IC.HH(z),2))  - IC.Hpr2(z)/pow(IC.HH(z),3);

  double tmp3 = -(2.0/(IC.chi(z)*IC.HH(z)))*(2.0 - 4.0*IC.fg(z) + 2.0*IC.bevo(z) -2.0*IC.curlyQ(z) + 4.0*IC.fg(z)*IC.curlyQ(z) - 2.0*IC.bevo(z)*IC.curlyQ(z) + 8.0*pow(IC.curlyQ(z),2) - 8.0*IC.curlyQLpr(z) +  3.0*(IC.curlyQ(z)-1.0)*IC.Hpr(z)/pow(IC.HH(z),2)  - 2.0*IC.curlyQpr(z)/IC.HH(z) );

  double tmp4 = 2.0*(1.0-IC.curlyQ(z) + 2.0*pow(IC.curlyQ(z),2) - 2.0*IC.curlyQLpr(z) )/(pow(IC.chi(z)*IC.HH(z),2) );

 

 double firstbatch = 9.0*pow(IC.Omegam(z),2)*(tmp1 + tmp2 +tmp3 +tmp4 )/4.0;

 
 double tmp6 = 13.0 - 10.0*IC.fg(z) + 8.0*IC.bevo(z) + IC.fg(z)*IC.bevo(z)- 28.0*IC.curlyQ(z) + 2.0*pow(IC.fg(z),2)- 2.0*pow(IC.bevo(z),2) + 4.0*IC.fg(z)*IC.curlyQ(z) + 8.0* IC.bevo(z)*IC.curlyQ(z)- 2.0*IC.bevopr(z)/IC.HH(z) + 8.0*IC.curlyQpr(z)/IC.HH(z);
 
 double tmp7  =  2.0*(IC.fg(z) -7.0 + 2.0*IC.bevo(z) - IC.fg(z)*IC.curlyQ(z) - 2.0*IC.bevo(z)*IC.curlyQ(z) - 2.0*IC.curlyQpr(z))/(IC.chi(z)*IC.HH(z));


 double secondbatch = 3.0*IC.Omegam(z)*IC.fg(z)*(tmp6 + tmp7)/2.0;

 
 double tmp8  =  6.0 - 9.0*IC.bevo(z)  + pow(IC.bevo(z),2) + IC.bevopr(z)/IC.HH(z) + (IC.bevo(z)-3.0)*IC.Hpr(z)/pow(IC.HH(z),2);
 
 double thirdbatch = pow(IC.fg(z),2)*tmp8  - 3.0* IC.Omegam(z)*IC.fgpr(z)/(2.0*IC.HH(z));
  
  double tmp =  firstbatch + secondbatch + thirdbatch;
 
  return tmp* pow(IC.HH(z),4);
	
}


double Beta2_in(double fnl, double z)
{
	IntermediateCosmology IC;
	double tmp1A = 1.0 - IC.bevo(z) + 2.0*IC.curlyQ(z) + 2.0*(1.0- IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) + IC.Hpr(z)/pow(IC.HH(z),2);

    double tmp1 = -4.5*pow(IC.Omegam(z),2)*tmp1A;


    double tmp2A = 2.0 + IC.fg(z)*(-1.0 - IC.bevo(z) + 2.0*IC.curlyQ(z) + 2.0*(1.0- IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) + IC.Hpr(z)/pow(IC.HH(z),2)) + IC.bevo(z) - 4.0*IC.curlyQ(z) - 2.0*(1.0-IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) - 2.0*IC.Hpr(z)/pow(IC.HH(z),2) ;
    double tmp2 = -3.0*IC.Omegam(z)*IC.fg(z)*tmp2A + 3.0* IC.Omegam(z) *IC.fgpr(z)/IC.HH(z);


	double tmp = tmp1 + tmp2;

  return tmp*pow(IC.HH(z),4);
}




double Beta3_in(double fnl, double z)
{
	IntermediateCosmology IC;


    double termzero = 9.0*pow(IC.Omegam(z),2)*(IC.fg(z) -2.0 +  2.0*IC.curlyQ(z))/4.0;

	double tmp1 =  -2.0- IC.fg(z)*(-3.0 + IC.fg(z)- 2.0*IC.bevo(z) - 3.0*IC.curlyQ(z) - 4.0*(1.0- IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) - 2.0*IC.Hpr(z)/pow(IC.HH(z),2));

	double tmp2 =  -IC.fgpr(z)/IC.HH(z) + 3.0*IC.bevo(z) + pow(IC.bevo(z),2) - 6.0*IC.bevo(z)*IC.curlyQ(z) + 4.0*IC.curlyQ(z) + 8.0*pow(IC.curlyQ(z),2) - 8.0 *IC.curlyQLpr(z)- 6.0*IC.curlyQpr(z)/IC.HH(z) + IC.bevopr(z)/IC.HH(z);


	double tmp3 = 2.0*(1.0 -IC.curlyQ(z) + 2.0*pow(IC.curlyQ(z),2)  - 2.0*IC.curlyQLpr(z))/(pow(IC.chi(z)*IC.HH(z),2));

	double tmp4 = 2.0*(-1.0- 2.0*IC.bevo(z) + 2.0*IC.bevo(z)*IC.curlyQ(z) + IC.curlyQ(z) - 6.0*pow(IC.curlyQ(z),2) + 3.0*(1.0 - IC.curlyQ(z))*IC.Hpr(z)/pow(IC.HH(z),2) + 6.0*IC.curlyQLpr(z) + 2.0*IC.curlyQpr(z)/IC.HH(z))/(IC.chi(z)*IC.HH(z));

	double tmp5 = -(3.0 + 2.0* IC.bevo(z)  - 6.0*IC.curlyQ(z)  - 3.0*IC.Hpr(z)/pow(IC.HH(z),2))*(IC.Hpr(z)/pow(IC.HH(z),2)) -IC.Hpr2(z)/pow(IC.HH(z),3);

	
	
	double firtbatch =  1.5* IC.Omegam(z)*IC.fg(z)*(tmp1 + tmp2 + tmp3 + tmp4 + tmp5);
	
	
	double tmp6 = pow(IC.fg(z),2)*(-3.0 + 2.0*IC.bevo(z)*(2.0 + (1.0-IC.curlyQ(z))/(IC.chi(z)*IC.HH(z))) - pow(IC.bevo(z),2) +  2.0* IC.bevo(z)*IC.curlyQ(z) -6.0*IC.curlyQ(z)  -  IC.bevopr(z)/IC.HH(z)- 6.0*(1.0 -IC.curlyQ(z))/(IC.chi(z)*IC.HH(z))) + 2.0*(1.0- 1.0/(IC.chi(z)*IC.HH(z))*IC.curlyQpr(z));
       
	double tmp  =  termzero + firtbatch + tmp6 ;

	return tmp* pow(IC.HH(z),3);
}




double Beta4_in(double fnl, double z)
{
	IntermediateCosmology IC;
        
   double tmpA1 =  -IC.bevo(z) + 2.0 *IC.curlyQ(z) +  2.0*(1.0-IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) + 2.0*IC.Hpr(z)/pow(IC.HH(z),2);

    double tmp1 = 1.5*IC.Omegam(z)*(tmpA1);

    double tmpA2 =  -IC.bevo(z) + 2.0 *IC.curlyQ(z) +  2.0*(1.0-IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) + IC.Hpr(z)/pow(IC.HH(z),2);

    double tmp2 = 2.0* pow(IC.fg(z),2)*tmpA2;

     double tmp = tmp1 + tmp2;
	
	return 0.0;//tmp* pow(IC.HH(z),3);

}

double Beta5_in(double fnl, double z)
{
	IntermediateCosmology IC;
    double tmpA1 =  IC.bevo(z) - 2.0 *IC.curlyQ(z) -  2.0*(1.0-IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) - IC.Hpr(z)/pow(IC.HH(z),2);

    double tmp = 3.0*IC.Omegam(z)*(tmpA1);

    return 0.0;//tmp* pow(IC.HH(z),3);

}


double Beta6_in(double fnl, double z)
{
	IntermediateCosmology IC;
    double tmpA1 =  2.0 - 2.0*IC.fg(z) + IC.bevo(z) - 4.0 *IC.curlyQ(z) -  2.0*(1.0-IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) - IC.Hpr(z)/pow(IC.HH(z),2);

    double tmp = 1.5*IC.Omegam(z)*(tmpA1);

    return tmp* pow(IC.HH(z),2);

}

double Beta7_in(double fnl, double z)
{
	IntermediateCosmology IC;
	
    double tmp = (3.0 - IC.bevo(z))*IC.fg(z);
	return tmp* pow(IC.HH(z),2);
}

double Beta8_in(double fnl, double z)
{
	IntermediateCosmology IC;
	

    double tmp1 = 3.0*IC.Omegam(z)*IC.fg(z)*(2.0 -IC.fg(z) - 2.0*IC.curlyQ(z));

    double tmp2 =  4.0 + IC.bevo(z) - pow(IC.bevo(z),2) + 4.0*IC.bevo(z)*IC.curlyQ(z) - 6.0*IC.curlyQ(z) - 4.0*pow(IC.curlyQ(z),2) + 4.0*IC.curlyQLpr(z) + 4.0*IC.curlyQpr(z)/IC.HH(z) - IC.bevopr(z)/IC.HH(z);

    double tmp3 = -2.0*(1.0-IC.curlyQ(z) + 2.0*pow(IC.curlyQ(z),2) - 2.0*IC.curlyQLpr(z))/pow(IC.chi(z)*IC.HH(z),2);

    double tmp4 = -2.0*(3.0 - 2.0*IC.bevo(z) + 2.0*IC.bevo(z)*IC.curlyQ(z) - IC.curlyQ(z) - 4.0*pow(IC.curlyQ(z),2) + 3.0*(1.0 -IC.curlyQ(z))*IC.Hpr(z)/pow(IC.HH(z),2) + 4.0*IC.curlyQLpr(z) + 2.0*IC.curlyQpr(z)/IC.HH(z));

	double tmp5  =  -(IC.Hpr(z)/pow(IC.HH(z),2))*(3.0 - 2.0*IC.bevo(z) + 4.0*IC.curlyQ(z) +3.0*IC.Hpr(z)/pow(IC.HH(z),2)) + IC.Hpr2(z)/pow(IC.HH(z),3);

        double tmp   =  tmp1 + pow(IC.fg(z),2)*( tmp2 + tmp3 + tmp4 + tmp5);

	return tmp* pow(IC.HH(z),2);
}



 double Beta9_in(double fnl, double z)
 {
 	IntermediateCosmology IC;
 	double tmp = -4.5* IC.Omegam(z)*IC.fg(z) ;//- 2.0*pow(IC.fg(z),2);
 	return 0.0;//tmp*pow(IC.HH(z),2);
 }

  double Beta10_in(double fnl, double z)
 {
 	IntermediateCosmology IC;
 	double tmp = 3.0* IC.Omegam(z)*IC.fg(z);
 	return 0.0;//tmp*pow(IC.HH(z),2);
 }



double Beta11_in(double fnl, double z)
{
	IntermediateCosmology IC;
      

    double tmp1 = -pow(IC.fg(z),2)*(-1.0 + IC.bevo(z) -2.0*IC.curlyQ(z) -  2.0*(1.0 + IC.curlyQ(z))/(IC.chi(z)*IC.HH(z))) - IC.Hpr(z)/pow(IC.HH(z),2);
	double tmp2 = 3.0* IC.Omegam(z)*IC.fg(z);
    double tmp = tmp1  + tmp2;

	return tmp* pow(IC.HH(z),2);
}

double Beta12_in(double fnl, double z)
{
	IntermediateCosmology IC;
    double tmp1 = IC.bias1(fnl,z)*(2.0 + IC.bevo(z) - 4.0 *IC.curlyQ(z) - 2.0*(1.0-IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) - IC.Hpr(z)/pow(IC.HH(z),2));

    double tmp2 = IC.bias1pr(fnl, z)/IC.HH(z) + 2.0*(2.0-1/(IC.chi(z)*IC.HH(z)))*IC.bias1Lpr(fnl,z);

	double tmp3 =  -IC.fg(z)*(IC.bias1(fnl,z)*( IC.fg(z) -3.0 + IC.bevo(z)) + IC.bias1pr(fnl, z)/IC.HH(z));


	double tmp = 1.5*IC.Omegam(z)*(tmp1 + tmp2) + tmp3;

	return tmp* pow(IC.HH(z),2);
}

double Beta13_in(double fnl, double z)
{
	IntermediateCosmology IC;
     
     double tmp1 = 9.0*pow(IC.Omegam(z),2)/4.0;

     double tmp2 = 1.5*IC.Omegam(z)*IC.fg(z)*(1.0 - 2.0*IC.fg(z) + 2.0*IC.bevo(z)- 6.0*IC.curlyQ(z)- 4.0*(1.0-IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) - 3.0*IC.Hpr(z)/pow(IC.HH(z),2)) ; 

     double tmp3 =  pow(IC.fg(z),2)*(3.0-IC.bevo(z)) ;

     double tmp =  tmp1 + tmp2 +  tmp3;

	return tmp* pow(IC.HH(z),2);
}

double Beta14_in(double fnl, double z)
{
	IntermediateCosmology IC;

     //double fnl =0.0;//
     double tmp = -1.5*IC.Omegam(z)*IC.bias1(fnl,z);
	
	return tmp*IC.HH(z);
}

double Beta15_in(double fnl, double z)
{
	IntermediateCosmology IC;
	double tmp = 2.0*pow(IC.fg(z),2);
	return tmp*IC.HH(z);
}



double Beta16_in(double fnl, double z)
{
	IntermediateCosmology IC;

      
	double tmp1 = IC.bias1(fnl,z)*(IC.fg(z) + IC.bevo(z) - 2.0 *IC.curlyQ(z) - IC.Hpr(z)/pow(IC.HH(z),2) - 2.0*(1.0-IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)));
    double tmp2 = IC.bias1pr(fnl, z)/IC.HH(z) + 2.0*(1.0-1.0/(IC.chi(z)*IC.HH(z)))*IC.bias1Lpr(fnl,z) ;
        
	double tmp = IC.fg(z)*(tmp1 + tmp2);

	return tmp*IC.HH(z);
}


double Beta17_in(double fnl, double z)
{
	IntermediateCosmology IC;
double tmp = -1.5* IC.Omegam(z)*IC.fg(z);
	return tmp*IC.HH(z);	
}

double Beta18_in(double fnl, double z)
{
	IntermediateCosmology IC;
double tmp1 = 1.5* IC.Omegam(z)*IC.fg(z);

 double tmp2 = pow(IC.fg(z),2)*(3.0 - 2.0*IC.bevo(z)  + 4.0 *IC.curlyQ(z) + 4.0*(1.0-IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) + 3.0*IC.Hpr(z)/pow(IC.HH(z),2));

 //double tmp3 = IC.fg(z)*(IC.bevo(z) - 2.0*IC.curlyQ(z) -2.0*(1.0-IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) -IC.Hpr(z)/pow(IC.HH(z),2));

  double tmp = tmp1 + tmp2;// + tmp3;

   return tmp*IC.HH(z);		
}

double Beta19_in(double fnl, double z)
{
	IntermediateCosmology IC;
        double tmp1 = IC.fg(z)*(IC.bevo(z) - 2.0*IC.curlyQ(z) -2.0*(1-IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) - IC.Hpr(z)/pow(IC.HH(z),2));

	return tmp1*IC.HH(z);

}



/////////////////////////////////////////Alphas////////////////////////

double alpha1_in(double fnl, double z)
{
	IntermediateCosmology IC;
 
 double tmp1A = 9.0 - 4.0*IC.fg(z) - IC.bevo(z) -2.0*IC.curlyQ(z) + 2.0*(1.0- IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) + IC.Hpr(z)/pow(IC.HH(z),2) + (6.0*fnl/5.0)*(-5.0 + 2.0*IC.fg(z) + IC.bevo(z) - 2.0*(1.0 - IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) - IC.Hpr(z)/pow(IC.HH(z),2));

double tmp1 = (9.0*pow(IC.Omegam(z),2)/4.0)*tmp1A;

double tmp2A = -7.0 -2.0*IC.fg(z)  + 2.0*IC.fgpr(z)/IC.HH(z) + 4.0*IC.bevo(z) - 4.0*IC.curlyQ(z) - 2.0*(1.0-IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) - 3.0*IC.Hpr(z)/pow(IC.HH(z),2) + (6.0*fnl/5.0)*(5.0 -4.0*IC.curlyQ(z) - 2.0*(1.0-IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) -4.0*IC.Hpr(z)/pow(IC.HH(z),2));

double tmp2 = 1.5*IC.Omegam(z)*IC.fg(z)*tmp2A;

 double tmp3A = -2.0 + 2.0*IC.fg(z) - IC.bevo(z) + 4.0*IC.curlyQ(z) + 2.0*(1.0 - IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) + 3.0*IC.Hpr(z)/pow(IC.HH(z),2) - 12.0*fnl/5.0;

double tmp3 =  1.5*IC.Omegam(z)*pow(IC.fg(z),2)*tmp3A;

 double tmp4 = 6.0*fnl*pow(IC.fg(z),2)*(3.0 - IC.bevo(z))/5.0 - 1.5*IC.Omegam(z)*IC.fgpr(z)*(1.0 + 6*fnl/5.0)/IC.HH(z);

 double tmp = tmp1 + tmp2 + tmp3 + tmp4;

  return tmp*pow(IC.HH(z),4);
	
}

double alpha2_in(double fnl, double z)
{
	IntermediateCosmology IC;
	double tmp1A = -1.0 + IC.bevo(z) -2.0*IC.curlyQ(z) - 2.0*(1.0- IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) - IC.Hpr(z)/pow(IC.HH(z),2);

        double tmp1 = 4.5*pow(IC.Omegam(z),2)*tmp1A;


      double tmp2A = 1.0 + 2.0*IC.fg(z) - 2.0*IC.bevo(z) + 4.0*IC.curlyQ(z) + 2.0*(1.0- IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) + 3.0*IC.Hpr(z)/pow(IC.HH(z),2) ;
       double tmp2 = 3.0*IC.Omegam(z)*IC.fg(z)*tmp2A;

      double tmp3A = -1.0 + IC.bevo(z) - 2.0*IC.curlyQ(z) - 2.0*(1.0- IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) - IC.Hpr(z)/pow(IC.HH(z),2) ;
       double tmp3 = 3.0*IC.Omegam(z)*pow(IC.fg(z),2)*tmp3A + 3.0*IC.Omegam(z)*IC.fgpr(z)/IC.HH(z);


	double tmp = tmp1 + tmp2 + tmp3 ;

  return tmp*pow(IC.HH(z),4);
}

double alpha3_in(double fnl, double z)
{
	IntermediateCosmology IC;
       double tmp1A = -1.0 + IC.bevo(z) -2.0*IC.curlyQ(z) - 2.0*(1.0- IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) - IC.Hpr(z)/pow(IC.HH(z),2);
       double tmp1B = 1.0 -2.0*fnl*(1.0 + 2.0*IC.fg(z)/IC.Omegam(z))/5.0;

        double tmp1 = 4.5*pow(IC.Omegam(z),2)*IC.fg(z)*tmp1A*tmp1B;

       
	return tmp1*pow(IC.HH(z),3);

}



double alpha4_in(double fnl,double z)
{
	IntermediateCosmology IC;
	double tmp = 3.0*IC.Omegam(z)*IC.fg(z)*(IC.bevo(z) - 2.0 *IC.curlyQ(z)  - 2.0*(1.0-IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) -IC.Hpr(z)/pow(IC.HH(z),2));
	
	return tmp*pow(IC.HH(z),3);
}

double alpha5_in(double fnl,double z)
{
	IntermediateCosmology IC;
	double tmp = -4.5*IC.Omegam(z)*IC.fg(z)*(1.0-2.0*fnl*(1.0 + 2.0*IC.fg(z)/(3.0*IC.Omegam(z)))/5.0);
	return tmp*pow(IC.HH(z),2);
}

double alpha6_in(double fnl,double z)
{
	IntermediateCosmology IC;
       double tmp1 = 3.0*IC.Omegam(z)*IC.fg(z);

	return tmp1*pow(IC.HH(z),2);
}


////////////////////////////Spline them for speed up////////////////////////////////////////////////




std::vector <double> z1_tab(k_bins2);
std::vector <double> Beta1_tab(k_bins2);
std::vector <double> Beta2_tab(k_bins2);
std::vector <double> Beta3_tab(k_bins2);
std::vector <double> Beta4_tab(k_bins2);
std::vector <double> Beta5_tab(k_bins2);
std::vector <double> Beta6_tab(k_bins2);
std::vector <double> Beta7_tab(k_bins2);
std::vector <double> Beta8_tab(k_bins2);
std::vector <double> Beta9_tab(k_bins2);
std::vector <double> Beta10_tab(k_bins2);
std::vector <double> Beta11_tab(k_bins2);
std::vector <double> Beta12_tab(k_bins2);
std::vector <double> Beta13_tab(k_bins2);
std::vector <double> Beta14_tab(k_bins2);
std::vector <double> Beta15_tab(k_bins2);
std::vector <double> Beta16_tab(k_bins2);
std::vector <double> Beta17_tab(k_bins2);
std::vector <double> Beta18_tab(k_bins2);
std::vector <double> Beta19_tab(k_bins2);
std::vector <double> Beta20_tab(k_bins2);



std::vector <double> alpha1_tab(k_bins2);
std::vector <double> alpha2_tab(k_bins2);
std::vector <double> alpha3_tab(k_bins2);
std::vector <double> alpha4_tab(k_bins2);
std::vector <double> alpha5_tab(k_bins2);
std::vector <double> alpha6_tab(k_bins2);

std::vector <double> sigma_tab(k_bins2);



void Spectroscopic::InitialzieGamma(void)
{
	   
	double  EPSREL = 1e-8;
	double z_Dz_min2;
	double z_Dz_max2;
	z_Dz_min2 = log(1e-3);
	z_Dz_max2 = log(10.0);
	const double QMIN = 1e-4;
	const double QMAX = 1e5;
	
	
	
	
#pragma omp parallel private(thread_number)
	{
#pragma omp for schedule(static) nowait
		for(int i = 0;i< k_bins2 ;i++)
		{
			z1_tab[i] = ((z_Dz_max2-z_Dz_min2)*((double) i)/((double)  k_bins2-1) + z_Dz_min2 );
			
			            Beta1_tab[i] = Beta1_in(0.0,exp(z1_tab[i]));
                        Beta2_tab[i] = Beta2_in(0.0,exp(z1_tab[i])); 
                        Beta3_tab[i] = Beta3_in(0.0,exp(z1_tab[i]));
                        Beta4_tab[i] = Beta4_in(0.0,exp(z1_tab[i]));
                        Beta5_tab[i] = Beta5_in(0.0,exp(z1_tab[i]));
                        Beta6_tab[i] = Beta6_in(0.0,exp(z1_tab[i])); 
                        Beta7_tab[i] = Beta7_in(0.0,exp(z1_tab[i]));
                        Beta8_tab[i] = Beta8_in(0.0,exp(z1_tab[i]));
                        Beta9_tab[i] = Beta9_in(0.0,exp(z1_tab[i]));
                        Beta10_tab[i] = Beta10_in(0.0,exp(z1_tab[i])); 
                        Beta11_tab[i] = Beta11_in(0.0,exp(z1_tab[i]));
                        Beta12_tab[i] = Beta12_in(0.0,exp(z1_tab[i]));
                        Beta13_tab[i] = Beta13_in(0.0,exp(z1_tab[i]));
                        Beta14_tab[i] = Beta14_in(0.0,exp(z1_tab[i])); 
                        Beta15_tab[i] = Beta15_in(0.0,exp(z1_tab[i]));
                        Beta16_tab[i] = Beta16_in(0.0,exp(z1_tab[i]));
                        Beta17_tab[i] = Beta17_in(0.0,exp(z1_tab[i]));
                        Beta18_tab[i] = Beta18_in(0.0,exp(z1_tab[i]));
                        Beta19_tab[i] = Beta19_in(0.0,exp(z1_tab[i]));
                       
                       
                        alpha1_tab[i] = alpha1_in(0.0,exp(z1_tab[i]));
                        alpha2_tab[i] = alpha2_in(0.0,exp(z1_tab[i])); 
                        alpha3_tab[i] = alpha3_in(0.0,exp(z1_tab[i]));
                        alpha4_tab[i] = alpha4_in(0.0,exp(z1_tab[i]));
                        alpha5_tab[i] = alpha5_in(0.0,exp(z1_tab[i]));
                        alpha6_tab[i] = alpha6_in(0.0,exp(z1_tab[i])); 
                       
			
//cout<<"\t"<<"GAMMA"<<"\t"<<exp(z1_tab[i])<<"\t"<<Beta1_tab[i]<<"\t"<<Beta2_tab[i]<<"\t"<<Beta3_tab[i]<<"\t"<<Beta4_tab[i]<<"\t"<<Beta5_tab[i]<<"\t"<<Beta6_tab[i]<<"\t"<<Beta7_tab[i]<<"\t"<<Beta8_tab[i]<<"\t"<<Beta9_tab[i]<<"\t"<<Beta10_tab[i]<<"\t"<<Beta11_tab[i]<<"\t"<<Beta12_tab[i]<<"\t"<<Beta13_tab[i]<<"\t"<<Beta14_tab[i]<<"\t"<<Beta15_tab[i]<<"\t"<<Beta16_tab[i]<<"\t"<<Beta17_tab[i]<<"\t"<<Beta18_tab[i]<<"\t"<<Beta19_tab[i]<<endl;	
		}
		
	}
}

double Beta1(double z)
{
	Spline<double, double> CubicSpline_g1(z1_tab,Beta1_tab);
	return CubicSpline_g1.interpolate(log(z));
   	
}

double Beta2(double z)
{
	Spline<double, double> CubicSpline_g2(z1_tab,Beta2_tab);
	return CubicSpline_g2.interpolate(log(z));
   	
}
double Beta3(double z)
{
	Spline<double, double> CubicSpline_g3(z1_tab,Beta3_tab);
	return CubicSpline_g3.interpolate(log(z));
   	
}
double Beta4(double z)
{
	Spline<double, double> CubicSpline_g4(z1_tab,Beta4_tab);
	return CubicSpline_g4.interpolate(log(z));
   	
}
double Beta5(double z)
{
	Spline<double, double> CubicSpline_g5(z1_tab,Beta5_tab);
	return CubicSpline_g5.interpolate(log(z));
   	
}
double Beta6(double z)
{
	Spline<double, double> CubicSpline_g6(z1_tab,Beta6_tab);
	return CubicSpline_g6.interpolate(log(z));
   	
}
double Beta7(double z)
{
	Spline<double, double> CubicSpline_g7(z1_tab,Beta7_tab);
	return CubicSpline_g7.interpolate(log(z));
   	
}
double Beta8(double z)
{
	Spline<double, double> CubicSpline_g8(z1_tab,Beta8_tab);
	return CubicSpline_g8.interpolate(log(z));
   	
}
double Beta9(double z)
{
	Spline<double, double> CubicSpline_g9(z1_tab,Beta9_tab);
	return CubicSpline_g9.interpolate(log(z));
   	
}
double Beta10(double z)
{
	Spline<double, double> CubicSpline_g10(z1_tab,Beta10_tab);
	return CubicSpline_g10.interpolate(log(z));
   	
}

double Beta11(double z)
{
	Spline<double, double> CubicSpline_g11(z1_tab,Beta11_tab);
	return CubicSpline_g11.interpolate(log(z));
   	
}
double Beta12(double z)
{
	Spline<double, double> CubicSpline_g12(z1_tab,Beta12_tab);
	return CubicSpline_g12.interpolate(log(z));
   	
}
double Beta13(double z)
{
	Spline<double, double> CubicSpline_g13(z1_tab,Beta13_tab);
	return CubicSpline_g13.interpolate(log(z));
   	
}
double Beta14(double z)
{
	Spline<double, double> CubicSpline_g14(z1_tab,Beta14_tab);
	return CubicSpline_g14.interpolate(log(z));
   	
}
double Beta15(double z)
{
	Spline<double, double> CubicSpline_g15(z1_tab,Beta15_tab);
	return CubicSpline_g15.interpolate(log(z));
   	

}
double Beta16(double z)
{
	Spline<double, double> CubicSpline_g16(z1_tab,Beta16_tab);
	return CubicSpline_g16.interpolate(log(z));
   	
}
double Beta17(double z)
{
	Spline<double, double> CubicSpline_g17(z1_tab,Beta17_tab);
	return CubicSpline_g17.interpolate(log(z));
   	
}
double Beta18(double z)
{
	Spline<double, double> CubicSpline_g18(z1_tab,Beta18_tab);
	return CubicSpline_g18.interpolate(log(z));
   	
}
double Beta19(double z)
{
	Spline<double, double> CubicSpline_g19(z1_tab,Beta19_tab);
	return CubicSpline_g19.interpolate(log(z));
   	
}

double alpha1(double z)
{
	Spline<double, double> CubicSpline_a1(z1_tab,alpha1_tab);
	return CubicSpline_a1.interpolate(log(z));
   	
}

double alpha2(double z)
{
	Spline<double, double> CubicSpline_a2(z1_tab,alpha2_tab);
	return CubicSpline_a2.interpolate(log(z));
   	
}
double alpha3(double z)
{
	Spline<double, double> CubicSpline_a3(z1_tab,alpha3_tab);
	return CubicSpline_a3.interpolate(log(z));
   	
}
double alpha4(double z)
{
	Spline<double, double> CubicSpline_a4(z1_tab,alpha4_tab);
	return CubicSpline_a4.interpolate(log(z));
   	
}
double alpha5(double z)
{
	Spline<double, double> CubicSpline_a5(z1_tab,alpha5_tab);
	return CubicSpline_a5.interpolate(log(z));
   	
}
double alpha6(double z)
{
	Spline<double, double> CubicSpline_a6(z1_tab,alpha6_tab);
	return CubicSpline_a6.interpolate(log(z));
   	
}


double k1dotk2(double k1, double k2, double muk)
{
    
	double tmp = k1*k2*muk;
	return tmp;
}


double S_kernel(double z, double k1, double k2, double mu)
{
  
  double tmp = pow(mu,2) - 1.0/3.0;
   return tmp;
}

double NewtonianF2(double k1, double k2, double mu12)
{
   double first_term = 5.0/7.0;
   double second_term = (k1/k2 + k2/k1)*mu12/2.;
   double third_term = (2.0/7.0)*pow(mu12,2);

   double tmp =  first_term+second_term+third_term;
   return 2.0*tmp;
}

double NewtonianG2(double k1, double k2, double mu12)
{
   double first_term = 3.0/7.0;
   double second_term = (k1/k2 + k2/k1)*mu12/(2.0);
   double third_term = (4.0/7.0)*pow(mu12,2);
   double tmp =  first_term+second_term+third_term;
   return 2.0*tmp;
}





double KRSD_Multi(double z, double fnl,double k1,double k2,double mu1, double mu2)
{
	IntermediateCosmology IC;
/*
double tmprsd1 = 2.0*IC.fg(z) *pow(mu2,2)*(b_one(fnl, k1, z)+ pow(mu1,2)*IC.fg(z));
double tmprsd2 = 2.0*IC.fg(z) *pow(mu1,2)*(b_one(fnl, k2, z)+ pow(mu2,2)*IC.fg(z));

double tmp1  =  (tmprsd1+tmprsd2)/2.0; 

double tmpdop1 = mu1*k1*IC.fg(z) *(mu2*b_one(fnl, k1, z) + pow(mu1,2)*IC.fg(z))/k2;
double tmpdop2 = mu2*k2*IC.fg(z) *(mu1*b_one(fnl, k2, z) + pow(mu2,2)*IC.fg(z))/k1;
double tmp2 = (tmpdop1+tmpdop2)/2.0;  

return tmp1 + tmp2; 
*/



 double b1 = IC.bias1(0.0,z);

 double tmp1  =  pow(IC.fg(z),2)* mu1*mu2*pow((mu1*k1 + mu2*k2),2)/(k1*k2);

 double tmp2 = (b1*IC.fg(z)/(k1*k2))*((pow(mu1,2) + pow(mu2,2))*k1*k2 + mu1*mu2*(pow(k1,2) + pow(k2,2)));

 return tmp1 + tmp2;


}
double E2(double k1, double k2, double k3,double muk)
{
   //double k3 = sqrt(k1*k1 + k2*k2 + 2.0*k1dotk2(k1, k2,  muk));
   double first_term = 3.0;
   double second_term = 2.0*(k1/k2 + k2/k1)*muk;
   double third_term = pow(muk,2);

   double tmp =  first_term + second_term + third_term;

   return pow(k1*k2/(k3*k3),2)*tmp;
}


double CurlyRGR2(double fnl, double z, double k1, double k2,double k3,double muk, double mu1, double mu2,double mu3)
{

//double k3 = sqrt(k1*k1 + k2*k2 + 2.0*k1dotk2(k1, k2,  muk));
//double mu3 = (mu1*k1 + mu2* k2)/k3;

double tmp1 = Beta1(z)  + E2(k1,k2,k3,muk)*Beta2(z);

double tmp2 = pow(k1*k2/k3,2)*(NewtonianF2(k1,k2,muk)*Beta6(z) + NewtonianG2(k1,k2,muk)*Beta7(z));

double tmp3 = mu1*mu2*k1*k2* Beta8(z);

double tmp4 =  pow(mu3*k3,2)*(Beta9(z) + E2(k1,k2,k3,muk)*Beta10(z));

double tmp5 =  k1*k2*muk*Beta11(z);

double tmp6 = (pow(k1,2) + pow(k2,2))*Beta12(z);

double tmp7 =  (pow(mu1*k1,2) + pow(mu2*k2,2))*Beta13(z); 

double tmp = (tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7)/pow(k1*k2,2);

return tmp;

}



double CurlyIGR2(double fnl, double z, double k1, double k2,double k3,double muk, double mu1, double mu2,double mu3)
{

//double k3 = sqrt(k1*k1 + k2*k2 + 2.0*k1dotk2(k1, k2,  muk));

//double mu3 = (mu1*k1 + mu2* k2)/k3;

double tmp1 = (mu1*k1 + mu2*k2)*Beta3(z) + mu3*k3*(Beta4(z) + E2(k1,k2,k3,muk) *Beta5(z));

double tmp2 = (mu1*pow(k1,3) + mu2*pow(k2,3))*Beta14(z);

double tmp3 = (mu1*k1 + mu2*k2)*k1*k2* muk*Beta15(z);

double tmp4 = (mu1*k2 + mu2*k1)*k1*k2*Beta16(z);

double tmp5 = (pow(mu1*k1,3)  + pow(mu2*k2,3))*Beta17(z);

double tmp6 = mu1*mu2*k1*k2*(mu1*k1 + mu2*k2)*Beta18(z);

double tmp7 = NewtonianG2(k1,k2,muk)*Beta19(z)*mu3*pow(k1*k2,2)/k3;


double tmp = (tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7)/pow(k1*k2,2);

return tmp;

}


double CurlyIGRSpatialone(double fnl, double z, double k1, double k2,double k3,double muk, double mu1, double mu2,double mu3)
{

//double k3 = sqrt(k1*k1 + k2*k2 + 2.0*k1dotk2(k1, k2,  muk));

//double mu3 = (mu1*k1 + mu2* k2)/k3;

double tmp1 = (mu1*k1 + mu2*k2)*Beta3(z) + mu3*k3*(Beta4(z) + E2(k1,k2,k3,muk) *Beta5(z));

return tmp1;

}

double CurlyIGRSpatialthree(double fnl, double z, double k1, double k2,double k3,double muk, double mu1, double mu2,double mu3)
{

//double k3 = sqrt(k1*k1 + k2*k2 + 2.0*k1dotk2(k1, k2,  muk));

//double mu3 = (mu1*k1 + mu2* k2)/k3;

//double tmp1 = (mu1*k1 + mu2*k2)*Beta3(z) + mu3*k3*(Beta4(z) + E2(k1,k2,k3,muk) *Beta5(z));

double tmp2 = (mu1*pow(k1,3) + mu2*pow(k2,3))*Beta14(z);

double tmp3 = (mu1*k1 + mu2*k2)*k1*k2* muk*Beta15(z);

double tmp4 = (mu1*k2 + mu2*k1)*k1*k2*Beta16(z);

double tmp5 = (pow(mu1*k1,3)  + pow(mu2*k2,3))*Beta17(z);

double tmp6 = mu1*mu2*k1*k2*(mu1*k1 + mu2*k2)*Beta18(z);

double tmp7 = NewtonianG2(k1,k2,muk)*Beta19(z)*mu3*pow(k1*k2,2)/k3;


double tmp = (tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7)/pow(k1*k2,2);

return tmp;

}










double mu23(double r, double mu12)
{
     ///double r = k2/k1;
     double y = sqrt(r*r+2.0*mu12 *r + 1.0);
     
     double tmp = -(r + mu12)/(y);
     return tmp;

}

double mu13(double r, double mu12)
{
     ///double r = k2/k1;
     double y = sqrt(r*r+2.0*mu12 *r + 1.0);
     double tmp =- (r* mu12 +1.0)/(y);
     return tmp;

}

double  Spectroscopic::MatterBispectrum(double z,double k1, double r, double mu12)
{
  LinearPowerSpectrumGR LP;

   double y = sqrt(r*r+2.0*mu12 *r + 1.0);
   double k3 = k1*y;
   double k2  = r* k1;
    

  double tmpA1 = NewtonianF2(k1, k2,mu12);
  double tmpA2 = LP.MatterPowerSpectrum(z,k1);
  double tmpA3 = LP.MatterPowerSpectrum(z,k2);
  double tmpA =  tmpA1*tmpA2*tmpA3;

  double tmpB1 = NewtonianF2(k3, k1,mu13(r, mu12));
  double tmpB2 = LP.MatterPowerSpectrum(z,k3);
  double tmpB3 = LP.MatterPowerSpectrum(z,k1);
  double tmpB = tmpB1*tmpB2*tmpB3;

  double tmpC1 = NewtonianF2(k2, k3,mu23(r, mu12));
  double tmpC2 = LP.MatterPowerSpectrum(z,k2);
  double tmpC3 = LP.MatterPowerSpectrum(z,k3);
  double tmpC =  tmpC1*tmpC2*tmpC3;

   return tmpA + tmpB + tmpC;
}




double Spectroscopic::ReducedMatterBispectrum_Q(double z,double k1, double r, double mu12)
{
 LinearPowerSpectrumGR LP;

   double y = sqrt(r*r+2.0*mu12 *r + 1.0);
   double k3 = k1*y;
   double k2  = r* k1;

   double PowerSPec_perm12 = LP.MatterPowerSpectrum(z,k1)*LP.MatterPowerSpectrum(z,k2);
   double PowerSPec_perm32 = LP.MatterPowerSpectrum(z,k3)*LP.MatterPowerSpectrum(z,k2);
   double PowerSPec_perm13 = LP.MatterPowerSpectrum(z,k1)*LP.MatterPowerSpectrum(z,k3);

   double tmp =  MatterBispectrum(z,k1,r,mu12)/(PowerSPec_perm12+PowerSPec_perm32+PowerSPec_perm13);

//cout<<"\t "<< mu1<<"\t"<<"B1"<<"\t"<<LinearPower_HIHI(1,z,fnl,k1,mu1)<<"\t"<<"B2"<<"\t"<<LinearPower_HIHI(1,fnl,z,k3,mu3)<<endl;
   return tmp;
}


double SecondOrderKernelRealSpace(int x,double fnl, double z, double k1, double k2, double k3,double mu12)
{   
IntermediateCosmology IC; 
     double HIb10 = IC.bias1(fnl,z);
    double HIb20 = IC.bias2(fnl,z);

   double HIb10k = IC.b_one(fnl,k1, z);
    double HIb20k = IC.b_two_integrand(fnl, z , k1, k2,mu12);



    if(x == 1)
       {
         return HIb10*NewtonianF2(k1,k2,mu12);
       }
    else if(x ==2)
      {
         return HIb20;
      }
    else if(x == 3)
      {
         return  IC.TidalBias(z)*S_kernel(z,k1, k2,mu12);
      }
    else if(x ==4)
      {
         return HIb10*NewtonianF2(k1,k2,mu12) + HIb20;
      }
    else if(x ==5)
      {
         return HIb10*NewtonianF2(k1,k2,mu12) + HIb20 + IC.TidalBias(z)*S_kernel(z,k1, k2,mu12);
      }

else 
      {
         return HIb10*NewtonianF2(k1,k2,mu12) + HIb20k + IC.TidalBias(z)*S_kernel(z,k1, k2,mu12);
      }

}





double SecondOrderKernelRSD(int x,double fnl, double z, double k1, double k2, double k3,double mu12, double mu1, double mu2, double mu3)
{    
	IntermediateCosmology IC;
     double HIb10 = IC.bias1(0.0,z);
    double HIb20 = IC.bias2(0.0,z);

//double k3 = sqrt(k1*k1 + k2*k2 + 2.0*k1dotk2(k1, k2,  mu12));

//double mu3 = (mu1*k1 + mu2* k2)/k3; 

 double b2k = IC.b_two_integrand(fnl, z , k1, k2,mu12);
 
    if(x == 1)
       {
         return HIb10*NewtonianF2(k1,k2,mu12);
       }
    else if(x ==2)
      {
         return HIb20;
      }
    else if(x == 3)
      {
         return  IC.TidalBias(z)*S_kernel(z,k1, k2,mu12);
      }
    else if(x ==4)
      {
         return HIb10*NewtonianF2(k1,k2,mu12) + HIb20;
      }
    else if(x ==5)
      {
         return HIb10*NewtonianF2(k1,k2,mu12) + HIb20 + IC.TidalBias(z)*S_kernel(z,k1, k2,mu12);
      }
    else if(x == 6)
      {
         return  pow(mu3,2)*IC.fg(z)*NewtonianG2(k1,k2,mu12);
 
      }
    else if(x == 7)
      {
         return  KRSD_Multi(z, fnl, k1, k2, mu1,  mu2);
 
      }
    else if(x == 8)
     {
        double tmp  = HIb10*NewtonianF2(k1,k2,mu12)
        + HIb20 + pow(mu3,2)*IC.fg(z)*NewtonianG2(k1,k2,mu12) +  KRSD_Multi(z,fnl, k1, k2, mu1, mu2);	
        return tmp;
}
  else if(x == 9)
     {
        double tmp  = HIb10*NewtonianF2(k1,k2,mu12)
        + HIb20 + pow(mu3,2)*IC.fg(z)*NewtonianG2(k1,k2,mu12) + KRSD_Multi(z, fnl, k1, k2, mu1, mu2)+IC.TidalBias(z)*S_kernel(z, k1, k2,mu12);	
        return tmp;
}
    else if(x ==10)
    {
       double tmp  = CurlyRGR2(fnl, z,  k1, k2,k3,mu12, mu1,  mu2,mu3);
       return tmp;
    }
    else if(x ==11)
    {
       double tmp  = CurlyIGR2(fnl, z,  k1, k2,k3,mu12, mu1,  mu2,mu3);
       return tmp;
    }
     else if(x ==12)
    {
       double tmp  = CurlyIGRSpatialone(fnl, z,  k1, k2,k3,mu12, mu1,  mu2,mu3);
       return tmp;
    }
    else if(x ==13)
    {
       double tmp  = CurlyIGRSpatialthree(fnl, z,  k1, k2,k3,mu12, mu1,  mu2,mu3);
       return tmp;
    }
    else
    {
      double tmp2  = HIb10*NewtonianF2(k1,k2,mu12) + HIb20 + pow(mu3,2)*IC.fg(z)*NewtonianG2(k1,k2,mu12)  + KRSD_Multi(z, fnl, k1, k2, mu1,  mu2) + CurlyRGR2(fnl, z,  k1, k2,k3,mu12, mu1,  mu2,mu3) + IC.TidalBias(z)*S_kernel(z, k1, k2,mu12);
     
    // double tmp2  = CurlyRGR2(fnl, z,  k1, k2,k3,mu12, mu1,  mu2,mu3) ;
     return tmp2;
    }

}








double mu2_in(double phi, double mu1, double muk)
{
   double tmp1 = mu1* muk + sqrt(1.0 - mu1*mu1 )*sqrt(1.0- muk*muk)*cos(phi);
   return tmp1;
}
double mu3_in(double k1,double k2, double k3,double phi, double mu1, double muk)
{
  double tmp1 = -(k1*mu1+k2*mu2_in(phi,mu1, muk))/k3;
  return tmp1;
}





double  Spectroscopic::RealSpaceBispectrum(int x,double fnl,double z,double k1, double r, double mu12)
{
	IntermediateCosmology IC;
	LinearPowerSpectrumGR LP;

     double y = sqrt(r*r + 2.0*mu12 *r + 1.0);
   double k3 = k1*y;
   double k2  = r* k1;


  //double r = k2/k1;
  //double k3 =  k1*sqrt(1.0+2.0*mu12*r+pow(r,2));
  
  double tmpA1 = SecondOrderKernelRSD(x, fnl, z, k1, k2,k3,mu12,0,0,0);
  double tmpA2 = IC.bias1(fnl,z);///Linearkernel(1, fnl,z, k1,mu1)
  double tmpA3 = IC.bias1(fnl,z);
  double tmpA =  tmpA1*tmpA2*tmpA3*LP.P1P2(z,k1,k2);

  double tmpB1 = SecondOrderKernelRSD(x, fnl, z, k1, k3,k2,mu13(r, mu12),0,0,0);
  double tmpB2 = IC.bias1(fnl,z);
  double tmpB3 = IC.bias1(fnl,z);
  double tmpB = tmpB1*tmpB2*tmpB3*LP.P1P2(z,k1,k3);

  double tmpC1 = SecondOrderKernelRSD(x, fnl, z, k2, k3,k1,mu23(r, mu12), 0,0,0);
  double tmpC2 = IC.bias1(fnl,z);
  double tmpC3 = IC.bias1(fnl,z);
  double tmpC =  tmpC1*tmpC2*tmpC3*LP.P1P2(z,k2,k3);

   return tmpA + tmpB + tmpC;
}

double Spectroscopic::ReducedRealSpaceBispectrum(int x,double fnl,double z,double k1, double r, double mu12)
{
 // double r = k2/k1;
 // double y = sqrt(r*r+2.0*mu12 *r + 1.0);
 // double k2 = r*k1; ///This ratio is the oppositte of what is in BCGS
  // double k3 =  k1*y;
	IntermediateCosmology IC;
	LinearPowerSpectrumGR LP;

  double y = sqrt(r*r+2.0*mu12 *r + 1.0);
   double k3 = k1*y;
   double k2  = r* k1;

double PowerSPec_perm12 = pow(IC.bias1(fnl,z),4)* LP.P1P2(z,k1,k2);
double PowerSPec_perm32 =pow(IC.bias1(fnl,z),4)* LP.P1P2(z,k2,k3);
double PowerSPec_perm13 = pow(IC.bias1(fnl,z),4)*LP.P1P2(z,k3,k1);
double tmp = RealSpaceBispectrum(x, fnl, z,k1,  r, mu12)/(PowerSPec_perm12+PowerSPec_perm32+PowerSPec_perm13);//See equation 155
return tmp;
}
/*

double  Spectroscopic::RealSpaceBispectrumfnl(int x,double fnl,double z,double k1, double r, double mu12)
{

     double y = sqrt(r*r + 2.0*mu12 *r + 1.0);
   double k3 = k1*y;
   double k2  = r* k1;


  //double r = k2/k1;
  //double k3 =  k1*sqrt(1.0+2.0*mu12*r+pow(r,2));
  
  double tmpA1 = SecondOrderKernelRealSpace(x, fnl, z, k1, k2,k3,mu12);
  double tmpA2 = b_one(fnl,k1, z);
  double tmpA3 = b_one(fnl,k2, z);
  double tmpA =  tmpA1*tmpA2*tmpA3*P1P2(z,k1,k2);

  double tmpB1 = SecondOrderKernelRealSpace(x, fnl, z, k1, k3,k2,mu13(r, mu12));
  double tmpB2 = b_one(fnl,k3, z);
  double tmpB3 = b_one(fnl,k1, z);
  double tmpB = tmpB1*tmpB2*tmpB3*P1P2(z,k1,k3);

  double tmpC1 = SecondOrderKernelRealSpace(x, fnl, z, k2, k3,k1,mu23(r, mu12));
  double tmpC2 = b_one(fnl,k2, z);
  double tmpC3 = b_one(fnl,k3, z);
  double tmpC =  tmpC1*tmpC2*tmpC3*P1P2(z,k2,k3);


   return tmpA + tmpB + tmpC;
}

double Spectroscopic::ReducedRealSpaceBispectrumfnl(int x,double fnl,double z,double k1, double r, double mu12)
{
 // double r = k2/k1;
 // double y = sqrt(r*r+2.0*mu12 *r + 1.0);
 // double k2 = r*k1; ///This ratio is the oppositte of what is in BCGS
  // double k3 =  k1*y;

  double y = sqrt(r*r+2.0*mu12 *r + 1.0);
   double k3 = k1*y;
   double k2  = r* k1;

double PowerSPec_perm12 = P1P2(z,k1,k2);
double PowerSPec_perm32 = P1P2(z,k2,k3);
double PowerSPec_perm13 = P1P2(z,k3,k1);
double tmp = RealSpaceBispectrumfnl(x, fnl, z,k1,  k2, mu12)/(PowerSPec_perm12+PowerSPec_perm32+PowerSPec_perm13);//See equation 155
return tmp;
}
*/
////////////////////////////////////////////Redshift space tools////////////////////////////////





double  Spectroscopic::NewtonianBispectrumWithRSD(int x,double fnl,double z,double k1, double k2, double mu12,double phi, double mu1)
{
	LinearPowerSpectrumGR LP;

   double r = k2/k1;
  double k3 =  k1*sqrt(1.0+2.0*mu12*r+pow(r,2));
  
  double mu2 = mu2_in(phi,  mu1,  mu12);
  double mu3 = mu3_in(k1,k2, k3, phi, mu1, mu12);// mu13(r, mu12)//mu23(r, mu12)
 

  double tmpA1 = SecondOrderKernelRSD(x, fnl, z, k1, k2, k3,mu12, mu1, mu2,mu3);
  double tmpA2 = LP.Linearkernel(1, fnl,z,k1,mu1);///Linearkernel(1, fnl,z, k1,mu1)
  double tmpA3 = LP.Linearkernel(1, fnl,z,k2,mu2);
  double tmpA =  tmpA1*tmpA2*tmpA3*LP.P1P2(z,k1,k2);

  double tmpB1 = SecondOrderKernelRSD(x, fnl, z, k3, k1,k2, mu13(r, mu12), mu3, mu1,mu2);
  double tmpB2 = LP.Linearkernel(1, fnl,z,k3,mu3);
  double tmpB3 = LP.Linearkernel(1, fnl,z,k1,mu1);
  double tmpB = tmpB1*tmpB2*tmpB3*LP.P1P2(z,k3,k1);

  double tmpC1 = SecondOrderKernelRSD(x, fnl, z, k2, k3,k1,mu23(r, mu12), mu2, mu3,mu1);
  double tmpC2 = LP.Linearkernel(1, fnl,z,k2,mu2);
  double tmpC3 = LP.Linearkernel(1, fnl,z,k3,mu3);
  double tmpC =  tmpC1*tmpC2*tmpC3*LP.P1P2(z,k2,k3);

   return tmpA +  tmpB + tmpC;//
}



double Spectroscopic::Passbyref(int x,int l , int m, void *params ,double phi, double mu1)
{
  double *fp = (double *) params;
  double fnl = fp[0];
  double  z = fp[1];
  double k1 = fp[2];
  double k2 = fp[3];
  //double k3 = fp[4];
  double mu12 = fp[4];
//NewtonianBispectrumWithRSD(int x,double fnl,double z,double k1, double k2, double mu12,double phi, double mu1)
  // return NewtonianBispectrumWithRSD(x,fnl,z,k1,k2,mu12,phi, mu1);

 return NewtonianBispectrumWithRSD(x,fnl,z,k1,k2,mu12,phi, mu1)*boost::math::spherical_harmonic_r(l,m, acos(mu1), phi);
 // return boost::math::spherical_harmonic_r(l,m, acos(mu1), phi);

}


double Spectroscopic::AveragedNewtonianBispectrumWithRSD(int x,int l , int m, double fnl,double z,double k1, double k2,double muk)
{
    const double EPSREL = 1e7;
    const double EPSABS = 1e7;
    double a[2] = {0.0,-0.998};
    double  b[2] = {2*M_PI,0.998};

   double fp[5];
    fp[0] =fnl;
    fp[1] = z;
    fp[2] = k1;
    fp[3] = k2;
   //fp[4] = k3;
    fp[4] = muk;

  double tmp =Integrate<2>(bind(&Spectroscopic::Passbyref,this,x,l,m,fp,_1,_2),a,b,EPSREL,EPSABS)/(sqrt(4.0*M_PI)* (2.0*l +1.0));

  //  double tmp = Integrate(bind(&Spectroscopic::Passbyref,this,x,l,m,fp,2.0*M_PI,_1),-0.998,0.998,EPSREL,EPSABS)/(4.0*M_PI*(2.0*l +1.0));
   return tmp;
}


double Spectroscopic::ReducedNewtonianBispectrumWithRSD(int x,int l , int m, double fnl,double z,double k1, double k2,double muk)
{
	LinearPowerSpectrumGR LP;
   double r = k2/k1;
 double y = sqrt(r*r+2.0*muk *r + 1.0);
   double k3 =  k1*y;
double PowerSPec_perm12 = LP.LinearPowerAverage(1,fnl, z,k1)*LP.LinearPowerAverage(1,fnl,z,k2);
double PowerSPec_perm32 = LP.LinearPowerAverage(1,fnl, z, k3)*LP.LinearPowerAverage(1,fnl,z,k2);
double PowerSPec_perm13 = LP.LinearPowerAverage(1,fnl,z, k1 )*LP.LinearPowerAverage(1,fnl,z,k3);

return  AveragedNewtonianBispectrumWithRSD(x,l,m, fnl,z,k1, k2,muk)/(PowerSPec_perm12+PowerSPec_perm32+PowerSPec_perm13);

}


double Spectroscopic::Testzterms(int x, double z)
{
    IntermediateCosmology ICC;

     double fnl =0.0;
    double s = 2.0*ICC.curlyQ(z)/5.0;
    double HIb10 = ICC.bias1(fnl,z); 
 
    double AA = ICC.bevo(z) + 1.5 * ICC.Omegam(z) - 3.0+ (2.0 - 5.0*s)*(1.0- 1.0/(ICC.chi(z)*ICC.HH(z)));
   if(x ==1)
{
     double tmp = HIb10* (AA+ ICC.fg(z))  + ICC.bias1pr(fnl,z)/ICC.HH(z) + (1.0- 1.0/(ICC.chi(z)*ICC.HH(z)))*ICC.bias1Lpr( fnl, z);
     return tmp;
}
else if(x == 2)
{
    return AA ;
}
else if(x ==3)
{
   double EE =  4.0 - 2.0* AA -1.5* ICC.Omegam(z);
   return EE;
}

else if(x ==4)
{
   double EE =  ICC.Omegam(z);
   return EE;
}

else if(x ==5)
{
   double EE =  ICC.fg(z);
   return EE;
}

else if(x ==6)
{
   double EE =  ICC.chi(z);
   return EE;
}

else if(x == 7)
{
   double EE =  s;
   return EE;
}

}

/*
std::vector<double > collection(double z)
 {
    std::vector<double> result;
    collection(std::back_inserter(result));
    return result;
}
*/

/////////////////////////////////////////////With Relativistic corrections///////////////////////////////////

/////////Terms with odd number of partial derivatives is imaginary in Fourier space


double DipoleKernel(int x,double z,double fnl,  double k1, double k2,double k3,double mu12, double mu1, double mu2,double mu3)
{

   	IntermediateCosmology IC;
      double s = 2.0*IC.curlyQ(z)/5.0;
      double HIb10 = IC.bias1(fnl,z);  

//double coef2=  IC.fg(z)*(IC.bevo(z) - 2.0*IC.curlyQ(z) -2.0*(1-IC.curlyQ(z))/(IC.chi(z)*IC.HH(z)) - IC.Hpr(z)/pow(IC.HH(z),2));

   double AA = IC.bevo(z) + 1.5 * IC.Omegam(z) - 3.0 + (2.0 - 5.0*s)*(1.0- 1.0/(IC.chi(z)*IC.HH(z)));
   double EE =  4.0 - 2.0* AA -1.5* IC.Omegam(z);
   double CC = HIb10* (AA+ IC.fg(z))  + IC.bias1pr(fnl,z)/IC.HH(z) + (1.0- 1.0/(IC.chi(z)*IC.HH(z)))*IC.bias1Lpr(fnl, z);
  

 
   double tmp1  =  -1.5*(mu1*k1/(k2*k2) + mu2*k2/(k1*k1))*IC.Omegam(z)*HIb10;

   double tmp2 =  2.0*mu12*(mu1/k2 + mu2/k1)*IC.fg(z)*IC.fg(z);

   double tmp3 =  (mu1/k1 +  mu2/k2)*IC.fg(z)*CC;

   double tmp4 = -1.5*(pow(mu1,3)*k1/(k2*k2) + pow(mu2, 3)*k2/(k1*k1)) *IC.Omegam(z)*IC.fg(z);

   double tmp5 =  mu1*mu2*(mu1/k2 +  mu2/k1)*(1.5*IC.Omegam(z) - EE*IC.fg(z))*IC.fg(z);

   double tmp6 = NewtonianG2(k1,k2,mu12)*AA*IC.fg(z)*mu3/k3;

    if(x == 1)
       {
         return IC.HH(z)*(tmp1);
       }
    else if(x ==2)
      {
         return IC.HH(z)*(tmp2);
      }
    else if(x == 3)
   {
    return IC.HH(z)*(tmp3);
    }
   else if( x ==4)
   {
    return IC.HH(z)*(tmp4);
   }
   else if(x ==5)
   { 
   return IC.HH(z)*(tmp5);
    }
    else if(x ==6)
    { 
   return IC.HH(z)*(tmp6);
    }
   else 
    {
   double tmp  = IC.HH(z)*(tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6);  
   return tmp;
}
}

double AnalyticalDipole(int x, double fnl, double z, double k1, double k2,double k3,double mu12, double mu1, double mu2,double mu3)
{

LinearPowerSpectrumGR LP;
//DipoleKernel(z, fnl, k1,k2, k3,mu12, mu1, mu2mu3)
double tmp1 = DipoleKernel(x,z, fnl, k1,k2, k3,mu12, mu1, mu2,mu3)*LP.Linearkernel(1, fnl,z,k1,mu1)*LP.Linearkernel(1, fnl,z,k2,mu2);
double tmp2 = (LP.Linearkernel(1, fnl,z,k1,mu1)*LP.Linearkernel(3, fnl,z,k2,mu2)+ LP.Linearkernel(3, fnl,z,k1,mu1)*LP.Linearkernel(1, fnl,z,k2,mu2))*SecondOrderKernelRSD(9, fnl, z, k1, k2, k3,mu12, mu1, mu2,mu3);


  return (tmp1 +tmp2)*LP.P1P2(z,k1,k2);
}


  
double Spectroscopic::GalBispectrumDopplerleading(int x,double fnl,double z,double k1, double k2, double mu12,double phi, double mu1)
{
      double r = k2/k1;
 double y = sqrt(r*r+2.0*mu12 *r + 1.0);
   double k3 =  k1*y;

  double mu2 = mu2_in(phi,  mu1,  mu12);
  double mu3 = mu3_in(k1,k2, k3, phi, mu1, mu12);// mu13(r, mu12)//mu23(r, mu12)
 double tmp1 = AnalyticalDipole(x,fnl,z,k1,k2,k3,mu12,mu1,mu2,mu3);
 double tmp2 = AnalyticalDipole(x,fnl,z,k1,k3,k2,mu13(r, mu12), mu1, mu3,mu2);
 double tmp3 = AnalyticalDipole(x,fnl,z,k2, k3,k1,mu23(r, mu12), mu2, mu3,mu1);

 return tmp1 + tmp2 + tmp3;

}



double Spectroscopic::PassbyrefGalDopplerleading(void *params,int x, int m, double phi, double mu1)
{
  double *fp = (double *) params;
  double fnl = fp[0];
  double  z = fp[1];
  double k1 = fp[2];
  double k2 = fp[3];
 // double k3 = fp[4];
  double mu12 = fp[4];

//double tmp =  -GalBispectrumDopplerleading(x,fnl, z, k1,k2,mu12,phi,mu1)*boost::math::spherical_harmonic_r(1,-1, acos(mu1), phi) + GalBispectrumDopplerleading(x,fnl, z, k1,k2,mu12,phi,mu1)*boost::math::spherical_harmonic_r(1,0, acos(mu1), phi) + GalBispectrumDopplerleading(x,fnl, z, k1,k2,mu12,phi,mu1)*boost::math::spherical_harmonic_r(1,1, acos(mu1), phi);

 ///double tmp =  GalBispectrumDopplerleading(x,fnl, z, k1,k2,mu12,phi,mu1)*boost::math::spherical_harmonic_r(1,m, acos(mu1), phi)/(4.0*M_PI*(2.0*1 +1.0));
 //return tmp;
   double coef1 =  (2./sqrt(2.))*(3./8.*M_PI);
  double coef2 =  (1./sqrt(2.))*(3./4.*M_PI);
 /* 
    if(m == -1)
       {
         return  coef1*GalBispectrumDopplerleading(x,fnl, z, k1,k2,mu12,phi,mu1)*sin( acos(mu1))*(cos(phi));
       }
else if(m == 0)
{
 return  coef2*GalBispectrumDopplerleading(x,fnl, z, k1,k2,mu12,phi,mu1)*cos(acos(mu1));
}
else if(m ==1)
{
         return  coef1*GalBispectrumDopplerleading(x,fnl, z, k1,k2,mu12,phi,mu1)*sin( acos(mu1))*(cos(phi));
       }
 else 
*/
return GalBispectrumDopplerleading(x,fnl, z, k1,k2,mu12,phi,mu1)*boost::math::spherical_harmonic_r(1,m, acos(mu1), phi);
//return  GalBispectrumDopplerleading(x,fnl, z, k1,k2,mu12,phi,mu1)*(cos(acos(mu1))+coef1*sin( acos(mu1))*(cos(phi)));
}



double Spectroscopic::GalBispectrumgDipoleleading(int x, int m, double fnl,double z,double k1, double k2,double muk)
{
    const double EPSREL = 1e6;
    const double EPSABS = 1e6;
    double a[2] = {0.0001,-0.998};
    double  b[2] = {2*M_PI,0.998};

   double fp[5];
    fp[0] =fnl;
    fp[1] = z;
   fp[2] = k1;
   fp[3] = k2;
  // fp[4] = k3;
   fp[4] = muk;  

  double tmp2 = Integrate<2>(bind(&Spectroscopic::PassbyrefGalDopplerleading,this,fp,x,m, _1,_2),a,b,EPSREL,EPSABS)/sqrt(4.0*M_PI*(2.0*1 +1.0));

  //double tmp2 = Integrate(bind(&Spectrosco  pic::PassbyrefGalDopplerleading,this,fp,x,m, 0.0,_1),-0.998,0.998,EPSREL,EPSABS)/sqrt(4.0*M_PI*(2.0*1 +1.0));
  
return tmp2;
}





///////////////////////////////////////////////////Full treatment

double ZkernelR(int y, int x,double fnl, double z, double k1, double k2,double k3,double mu12, double mu1, double mu2,double mu3)
{
	LinearPowerSpectrumGR LP;

  double KReal2 = SecondOrderKernelRSD(18, fnl, z, k1, k2,k3,mu12, mu1, mu2,mu3);

  double KRealN2 = SecondOrderKernelRSD(9, fnl, z, k1, k2,k3,mu12, mu1, mu2,mu3);
  double KRealNA = LP.Linearkernel(1, fnl,z, k1,mu1);
  double KRealNB = LP.Linearkernel(1, fnl,z, k2,mu2);

  double KReal1A = LP.Linearkernel(4, fnl,z, k1,mu1);
  double KReal1B = LP.Linearkernel(4, fnl,z, k2,mu2);


  double KImaginary2 = SecondOrderKernelRSD(11, fnl, z, k1, k2,k3,mu12, mu1, mu2,mu3);
  double KImaginary1A = LP.Linearkernel(3, fnl,z, k1,mu1);
  double KImaginary1B = LP.Linearkernel(3, fnl,z, k2,mu2);
//double tmp1;double tmp2;


   double KImaginary3 = SecondOrderKernelRSD(12, fnl, z, k1, k2,k3,mu12, mu1, mu2,mu3);


   double KImaginary4 = SecondOrderKernelRSD(13, fnl, z, k1, k2,k3,mu12, mu1, mu2,mu3);
if(y==1)
{
 double tmp1 = KReal2*(KReal1A*KReal1B - KImaginary1A*KImaginary1B) - KImaginary2*(KReal1A*KImaginary1B + KImaginary1A*KReal1B);
return  tmp1*LP.P1P2(z,k1,k2);
}
else if(y == 2)
{
 double tmp2 = KReal2*(KReal1A*KImaginary1B +  KImaginary1A*KReal1B) + KImaginary2*(KReal1A*KReal1B - KImaginary1A*KImaginary1B);
return tmp2*LP.P1P2(z,k1,k2);
}
else if(y == 3)
{
 double tmp2 = KReal2*(KReal1A*KImaginary1B +  KImaginary1A*KReal1B) + KImaginary3*(KReal1A*KReal1B - KImaginary1A*KImaginary1B);
return tmp2*LP.P1P2(z,k1,k2);
}
else if(y == 4)
{
 double tmp2 = KRealN2*(KRealNA*KImaginary1B +  KImaginary1A*KRealNB) + KImaginary4*(KRealNA*KRealNB);// - KImaginary1A*KImaginary1B);
return tmp2*LP.P1P2(z,k1,k2);
}
else 
return 0.0 ;
}





double Spectroscopic::GalBispectrumg(int w, int x,double fnl,double z,double k1, double k2, double mu12,double phi, double mu1)
{
      double r = k2/k1;
 double y = sqrt(r*r+2.0*mu12 *r + 1.0);
   double k3 =  k1*y;

  double mu2 = mu2_in(phi,  mu1,  mu12);
  double mu3 = mu3_in(k1,k2, k3, phi, mu1, mu12);// mu13(r, mu12)//mu23(r, mu12)
 double tmp1 = ZkernelR(w,x, fnl, z, k1, k2,k3,mu12, mu1, mu2,mu3);
 double tmp2 = ZkernelR(w,x, fnl, z, k1, k3,k2,mu13(r, mu12), mu1, mu3,mu2);
 double tmp3 = ZkernelR(w,x, fnl, z, k2, k3,k1,mu23(r, mu12), mu2, mu3,mu1);

 return tmp1 + tmp2 + tmp3;

}







double Spectroscopic::GalBispectrumabs(int x, double fnl,double z,double k1, double k2, double mu12,double phi, double mu1)
{
 double tmp1 = GalBispectrumg(1,x, fnl,z, k1,k2,mu12,phi, mu1);
 double tmp2 =  GalBispectrumg(2,x, fnl,z, k1,k2,mu12,phi, mu1);
//cout<<"\t "<<"\t"<<"Real"<<"\t"<<GalBispectrumg(1, fnl,z, k1,k2,k3,mu12,phi, mu1)<<"\t"<<"k3"<<"\t"<<k3<<"\t"<<"imaginary"<<"\t"<<GalBispectrumg(2, fnl,z, k1,k2,k3,mu12,phi, mu1)<<endl;
return sqrt(pow(tmp1,2) + pow(tmp2,2));
}

double Spectroscopic::PassbyrefGal(int w, int x, void *params ,double phi, double mu1)
{
  double *fp = (double *) params;
  double fnl = fp[0];
  double  z = fp[1];
  double k1 = fp[2];
  double k2 = fp[3];
 // double k3 = fp[4];
  double mu12 = fp[4];
//return GalBispectrumabs(fnl,z,k1,k2,mu12,phi, mu1);

return GalBispectrumg(w, x, fnl, z, k1,k2,mu12,phi,mu1);


}



double Spectroscopic::GalBispectrumgAverage(int w, int x, double fnl,double z,double k1, double k2,double muk)
{
    const double EPSREL = 1e-2;
    const double EPSABS = 1e-2;
    double a[2] = {0.00001,-0.998};
    double  b[2] = {2*M_PI,0.998};

   double fp[5];
    fp[0] =fnl;
    fp[1] = z;
   fp[2] = k1;
   fp[3] = k2;
  // fp[4] = k3;
   fp[4] = muk;
 
// double tmp = Integrate<2>(bind(&Spectroscopic::PassbyrefGal,this,x,fp,_1,_2),a,b,EPSREL,EPSABS)/(4.0*M_PI);
//return tmp;

  double tmp1 = Integrate<2>(bind(&Spectroscopic::PassbyrefGal,this,w,x,fp,_1,_2),a,b,EPSREL,EPSABS)/(4.0*M_PI); //GalBispectrumg(1, fnl,z, k1,k2,mu12,phi, mu1);
// double tmp2 =  Integrate<2>(bind(&Spectroscopic::PassbyrefGal,this,2,fp,_1,_2),a,b,EPSREL,EPSABS)/(4.0*M_PI);//GalBispectrumg(2, fnl,z, k1,k2,mu12,phi, mu1);
// return tmp1;
 // return sqrt(pow(tmp1,2) + pow(tmp2,2));
  
return tmp1;
}


double Spectroscopic::ReducedGalBispectrumgAverage(int w, int x, double fnl,double z,double k1, double k2,double muk)
{
	LinearPowerSpectrumGR LP;
   double r = k2/k1;
 double y = sqrt(r*r+2.0*muk *r + 1.0);
   double k3 =  k1*y;
double PowerSPec_perm12 = LP.LinearPowerAverage(3,fnl, z,k1)*LP.LinearPowerAverage(3,fnl,z,k2);
double PowerSPec_perm32 = LP.LinearPowerAverage(3,fnl, z, k3)*LP.LinearPowerAverage(3,fnl,z,k2);
double PowerSPec_perm13 = LP.LinearPowerAverage(3,fnl,z, k1 )*LP.LinearPowerAverage(3,fnl,z,k3);

return  GalBispectrumgAverage(w,x, fnl,z,k1, k2,muk)/(PowerSPec_perm12+PowerSPec_perm32+PowerSPec_perm13);

}






double Nl1l2L(double l1,double l2,double L)
{
    return (2.*l1+1)*(2.*l2+1)*(2.*L+1);
}

double Hl1l2L(double l1,double l2,double L)
{
    return WignerSymbols::wigner3j(l1, l2, L, 0, 0, 0);
}


double Spectroscopic::DiffBasisNewtonianbisepctrumintegrand(void *params,double Phi,double mu12,double mu1)
{
    double *fp = (double *) params;
   double fnl = fp[0];
   double  z = fp[1];
   double k1 = fp[2];
   double k2 = fp[3];
   double l1 = fp[4];
   double l2 = fp[5];
   double L = fp[6];
    //Theta_12 = acos(mu12)
   // Omega = acos(mu1)
       
   double Spcoef = sqrt(2.0*M_PI/(2.0*l2+ 1.0))*sqrt(2.0*M_PI/(2.0*L+ 1.0));
   double coef = Nl1l2L(l1,l2, L)*Hl1l2L(l1,l2, L);
   double sym = WignerSymbols::wigner3j(l1, l2, L, 0, 0, 0);   
   double GRBg = GalBispectrumg(1, 18, fnl, z, k1,k2,mu12,Phi,mu1);//NewtonianBispectrumWithRSD(9,fnl,z,k1,k2,mu12,Phi, mu1);
   double dimorm = 1.0/(2.0*4.0*M_PI);
   double sph2 = 0.5*sqrt(1.0/M_PI)*0.5*sqrt(1.0/M_PI);//boost::math::spherical_harmonic_r(l2,0.0, acos(mu12), phi)*boost::math::spherical_harmonic_r(L,0.0, acos(mu1), phi)
    double l1m0 = GRBg*sym * sph2* dimorm* Spcoef;
     
    return  l1m0;
}
  
double Spectroscopic::DiffBasisImaginaryGRbisepctrumintegrand3D(void *params,double Phi,double mu12,double mu1)
{
    double *fp = (double *) params;
   double fnl = fp[0];
   double  z = fp[1];
   double k1 = fp[2];
   double k2 = fp[3];
   double l1 = fp[4];
   double l2 = fp[5];
   double L = fp[6];
    //Theta_12 = acos(mu12)
   // Omega = acos(mu1)
       
   double Spcoef = sqrt(2.0*M_PI/(2.0*l2+ 1.0))*sqrt(2.0*M_PI/(2.0*L+ 1.0));
   double coef = Nl1l2L(l1,l2, L)*Hl1l2L(l1,l2, L);
   double sym = WignerSymbols::wigner3j(l1, l2, L, 0, 0, 0);   
   double GRBg = GalBispectrumg(2, 18, fnl, z, k1,k2,mu12,Phi,mu1) ;//GalBispectrumDopplerleading(7,fnl, z, k1,k2,mu12,Phi,mu1);
   double dimorm = 1.0/(2.0*4.0*M_PI);
   double sph2 = 0.5*sqrt(1.0/M_PI)*0.5*sqrt(1.0/M_PI);//boost::math::spherical_harmonic_r(l2,0.0, acos(mu12), phi)*boost::math::spherical_harmonic_r(L,0.0, acos(mu1), phi)
    double l1m0 = GRBg*sym * sph2* dimorm* Spcoef;
     
    return  l1m0;
}


double Spectroscopic::GalBispectrumgDipoleBasis(int xx, double fnl,double z,double k1, double k2,double  l1,double l2, double L)
{
    
    double a[3] = {0.0001,-0.998,-0.998};
    double  b[3] = {2*M_PI,0.998,0.998};

   double fp[7];
    fp[0] =fnl;
    fp[1] = z;
   fp[2] = k1;
   fp[3] = k2;
   fp[4] = l1;
   fp[5] = l2;
   fp[6] = L;  
if(xx == 1)
       {
const double EPSREL = 1e2;
    const double EPSABS = 1e2;
  double tmp2 = Integrate<3>(bind(&Spectroscopic::DiffBasisNewtonianbisepctrumintegrand,this,fp, _1,_2,_3),a,b,EPSREL,EPSABS);
  return tmp2;
}
else if(xx == 2)
{
const double EPSREL = 1e1;
    const double EPSABS = 1e1;
 double tmp2 = Integrate<3>(bind(&Spectroscopic::DiffBasisImaginaryGRbisepctrumintegrand3D,this,fp, _1,_2,_3),a,b,EPSREL,EPSABS);
  return tmp2;
}
else
{
cout<<"\t "<<"\t"<<"Wrong choice"<<"\t"<<endl;
return 0.0;
}
}


/////////////////////////////////Higher multipoles////////////////////////////////////////

///////////////////Gal multipoles////////////spherical_harmonic_r(unsigned n, int m, T1 theta, T2 phi)

double Spectroscopic::PassbyrefGalMulti(void *Iparams, void *params ,double phi, double mu1)
{
  double *fp = (double *) params;
  double fnl = fp[0];
  double  z = fp[1];
  double k1 = fp[2];
  double k2 = fp[3];
 // double k3 = fp[4];
  double mu12 = fp[4];
//return GalBispectrumabs(fnl,z,k1,k2,mu12,phi, mu1); 
   int *lp = (int *) Iparams;
   int w = lp[0];
   int l = lp[1];
   int m = lp[2];
//GalBispectrumg(int w, int x,double fnl,double z,double k1, double k2, double mu12,double phi, double mu1)
return GalBispectrumg(w,18,fnl, z, k1,k2,mu12,phi,mu1)*boost::math::spherical_harmonic_r(l,m, acos(mu1), phi);

}

double Spectroscopic::PassbyrefDipole(void *params , int w, double phi, double mu1)
{
  double *fp = (double *) params;
  double fnl = fp[0];
  double  z = fp[1];
  double k1 = fp[2];
  double k2 = fp[3];
 // double k3 = fp[4];
  double mu12 = fp[4];

//GalBispectrumg(int w, int x,double fnl,double z,double k1, double k2, double mu12,double phi, double mu1)
 double tmp =  GalBispectrumg(w,18,fnl, z, k1,k2,mu12,phi,mu1)*boost::math::spherical_harmonic_r(1,0, acos(mu1), phi);//
// + (GalBispectrumg(w,18,fnl, z, k1,k2,mu12,phi,mu1)*boost::math::spherical_harmonic_r(1,-1, acos(mu1), phi) -  GalBispectrumg(w,18,fnl, z, k1,k2,mu12,phi,mu1)*boost::math::spherical_harmonic_r(1,1, acos(mu1), phi))/sqrt(2.0);
  return tmp;

}

double Spectroscopic::GalBispectrumgDipole(int w, double fnl,double z,double k1, double k2,double muk)
{
    const double EPSREL = 1e7;
    const double EPSABS = 1e8;
    double a[2] = {0.00001,-0.998};
    double  b[2] = {2*M_PI,0.998};

   double fp[5];
    fp[0] =fnl;
    fp[1] = z;
   fp[2] = k1;
   fp[3] = k2;
  // fp[4] = k3;
   fp[4] = muk;
   //int x = 12;
  double tmp2 = Integrate<2>(bind(&Spectroscopic::PassbyrefDipole,this,fp, w,_1,_2),a,b,EPSREL,EPSABS)/(12.0*M_PI);
  
return tmp2;
}


double Spectroscopic::GalBispectrumgMulti(int w,int l, int m, double fnl,double z,double k1, double k2,double muk)
{
    const double EPSREL = 1e7;
    const double EPSABS = 1e8;
    double a[2] = {0.00001,-0.998};
    double  b[2] = {2*M_PI,0.998};

   double fp[5];
    fp[0] =fnl;
    fp[1] = z;
   fp[2] = k1;
   fp[3] = k2;
  // fp[4] = k3;
   fp[4] = muk;
   //int x = 12;
 
    int lparam[3];
    lparam[0] = w;  
    lparam[1] = l;
    lparam[2] = m;

/*
  int lparam2[3];
    lparam2[0] = 2;  
    lparam2[1] = l;
    lparam2[2] = m;
*/

// double tmp = Integrate<2>(bind(&Spectroscopic::PassbyrefGal,this,x,fp,_1,_2),a,b,EPSREL,EPSABS)/(4.0*M_PI);
//return tmp;

 // double tmp1 = (2.0*l +1.0)*Integrate<2>(bind(&Spectroscopic::PassbyrefGalMulti,this,lparam, fp,_1,_2),a,b,EPSREL,EPSABS)/(4.0*M_PI); //GalBispectrumg(1, fnl,z, k1,k2,mu12,phi, mu1);
  double tmp2 = Integrate<2>(bind(&Spectroscopic::PassbyrefGalMulti,this,lparam, fp,_1,_2),a,b,EPSREL,EPSABS)/(4.0*M_PI*(2.0*l +1.0)); //GalBispectrumg(1, fnl,z, k1,k2,mu12,phi, mu1);
  //return sqrt(pow(tmp1,2) + pow(tmp2,2));
  
return tmp2;
}


double Spectroscopic::ReducedGalBispectrumgMulti(int w,int l, int m, double fnl,double z,double k1, double k2,double muk)
{
	LinearPowerSpectrumGR LP;
   double r = k2/k1;
 double y = sqrt(r*r+2.0*muk *r + 1.0);
   double k3 =  k1*y;
double PowerSPec_perm12 = LP.LinearPowerAverage(3,fnl, z,k1)*LP.LinearPowerAverage(3,fnl,z,k2);
double PowerSPec_perm32 = LP.LinearPowerAverage(3,fnl, z, k3)*LP.LinearPowerAverage(3,fnl,z,k2);
double PowerSPec_perm13 = LP.LinearPowerAverage(3,fnl,z, k1 )*LP.LinearPowerAverage(3,fnl,z,k3);

return  GalBispectrumgMulti(w,l ,m,fnl,z,k1, k2,muk)/(PowerSPec_perm12+PowerSPec_perm32+PowerSPec_perm13);

}

//////////Special galaxy bispectrum in the Newtonian limit with relativictic first order terms//////////////




double ZkernelN(int y, int x, double fnl, double z, double k1, double k2,double k3,double mu12, double mu1, double mu2,double mu3)
{
	LinearPowerSpectrumGR LP;
  double KReal2 = SecondOrderKernelRSD(x, fnl, z, k1, k2,k3,mu12, mu1, mu2,mu3);
  double KReal1A = LP.Linearkernel(10, fnl,z, k1,mu1);
  double KReal1B = LP.Linearkernel(10, fnl,z, k2,mu2);


  double KImaginary2 = 0.0;//SecondOrderKernelRSD(11, fnl, z, k1, k2,k3,mu12, mu1, mu2,mu3);
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


double Spectroscopic::GalBispectrumN(int w,int x,double fnl,double z,double k1, double k2, double mu12,double phi, double mu1)
{
      double r = k2/k1;
 double y3 = sqrt(r*r+2.0*mu12 *r + 1.0);
   double k3 =  k1*y3;

  double mu2 = mu2_in(phi,  mu1,  mu12);
  double mu3 = mu3_in(k1,k2, k3, phi, mu1, mu12);// mu13(r, mu12)//mu23(r, mu12)

  
 double tmp1 = ZkernelN(w,x, fnl, z, k1, k2,k3,mu12, mu1, mu2,mu3);
 double tmp2 = ZkernelN(w,x, fnl, z, k1, k3,k2,mu13(r, mu12), mu1, mu3,mu2);
 double tmp3 = ZkernelN(w,x, fnl, z, k2, k3,k1,mu23(r, mu12), mu2, mu3,mu1);

 return tmp1 + tmp2 + tmp3;

}

double Spectroscopic::GalBispectrumNabs(int x, double fnl,double z,double k1, double k2, double mu12,double phi, double mu1)
{
 double tmp1 = GalBispectrumN(1,x, fnl,z, k1,k2,mu12,phi, mu1);
 double tmp2 =  GalBispectrumN(2,x, fnl,z, k1,k2,mu12,phi, mu1);
//cout<<"\t "<<"\t"<<"Real"<<"\t"<<GalBispectrumg(1, fnl,z, k1,k2,k3,mu12,phi, mu1)<<"\t"<<"k3"<<"\t"<<k3<<"\t"<<"imaginary"<<"\t"<<GalBispectrumg(2, fnl,z, k1,k2,k3,mu12,phi, mu1)<<endl;
return sqrt(pow(tmp1,2) + pow(tmp2,2));
}

double Spectroscopic::PassbyrefGalN(int w, int x, void *params ,double phi, double mu1)
{
  double *fp = (double *) params;
  double fnl = fp[0];
  double  z = fp[1];
  double k1 = fp[2];
  double k2 = fp[3];
 // double k3 = fp[4];
  double mu12 = fp[4];
//return GalBispectrumabs(fnl,z,k1,k2,mu12,phi, mu1);

return GalBispectrumN(w,x, fnl, z, k1,k2,mu12,phi,mu1);


}



double Spectroscopic::GalBispectrumNAverage(int w, int x, double fnl,double z,double k1, double k2,double muk)
{
    const double EPSREL = 1e-2;
    const double EPSABS = 1e-2;
    double a[2] = {0.00001,-0.998};
    double  b[2] = {2*M_PI,0.998};

   double fp[5];
    fp[0] =fnl;
    fp[1] = z;
   fp[2] = k1;
   fp[3] = k2;
  // fp[4] = k3;
   fp[4] = muk;
 
// double tmp = Integrate<2>(bind(&Spectroscopic::PassbyrefGal,this,x,fp,_1,_2),a,b,EPSREL,EPSABS)/(4.0*M_PI);
//return tmp;

  double tmp1 = Integrate<2>(bind(&Spectroscopic::PassbyrefGalN,this,w,x,fp,_1,_2),a,b,EPSREL,EPSABS)/(4.0*M_PI); //GalBispectrumg(1, fnl,z, k1,k2,mu12,phi, mu1);
// double tmp2 =  Integrate<2>(bind(&Spectroscopic::PassbyrefGal,this,2,fp,_1,_2),a,b,EPSREL,EPSABS)/(4.0*M_PI);//GalBispectrumg(2, fnl,z, k1,k2,mu12,phi, mu1);
// return tmp1;
 // return sqrt(pow(tmp1,2) + pow(tmp2,2));
  
return tmp1;
}


double Spectroscopic::ReducedGalBispectrumNAverage(int w, int x, double fnl,double z,double k1, double k2,double muk)
{
	LinearPowerSpectrumGR LP;
   double r = k2/k1;
   double y = sqrt(r*r+2.0*muk *r + 1.0);
   double k3 =  k1*y;
   double PowerSPec_perm12 = LP.LinearPowerAverage(3,fnl, z,k1)*LP.LinearPowerAverage(3,fnl,z,k2);
   double PowerSPec_perm32 = LP.LinearPowerAverage(3,fnl, z, k3)*LP.LinearPowerAverage(3,fnl,z,k2);
   double PowerSPec_perm13 = LP.LinearPowerAverage(3,fnl,z, k1 )*LP.LinearPowerAverage(3,fnl,z,k3);

return  GalBispectrumNAverage(w,x, fnl,z,k1, k2,muk)/(PowerSPec_perm12+PowerSPec_perm32+PowerSPec_perm13);

}




/*

////////////////////////////////Some comparison tools//////////////////////////////////////////////////

double Gamma7(double b0, double z)
{
        

         double fnl = 0.0;
  
        double tmp1 = bias1(b0,z)*(2.0 + IC.bevo(z) - 4.0 *curlyQ(z)  - 2.0*(1.0-curlyQ(z))/(chi(z)*IC.HH(z)) - IC.Hpr(z)/pow(IC.HH(z),2));

         double tmp2 = bias1pr(b0, z)/IC.HH(z) + 2.0*(2.0-1/(chi(z)*IC.HH(z)))*bias1Lpr(b0,z);

	double tmp3 =  fg(z)*(bias1(b0,z)*(3.0 - fg(z) -IC.bevo(z)) -bias1pr(b0, z)/IC.HH(z));


	double tmp = 1.5*IC.Omegam(z)*(tmp1 + tmp2) + tmp3;

	return tmp* pow(HH(z),2);
}

double Gamma8(double  b0,double z)
{
     
        double tmp1 = 9.0*pow(Omegam(z),2)/4.0;
        double tmp2 = 1.5*Omegam(z)*fg(z)*(3.0 - 2.0*fg(z) + 2.0*IC.bevo(z)- 3.0*curlyQ(z)- 4.0*(1.0-curlyQ(z))/(chi(z)*HH(z)) - IC.Hpr(z)/pow(HH(z),2)) ; 
        double tmp3 =  pow(fg(z),2)*(3.0-IC.bevo(z)) ;
        double tmp4  =   2.0*(pow(fg(z),2) - 3.0*Omegam(z)*fg(z)/4.0);


	double tmp =  tmp1 + tmp2 +  tmp3 + tmp4; 

	return tmp* pow(HH(z),2);
}

double Gamma9(double b0, double z)
{

     double fnl =0.0;//
     double tmp = -1.5*Omegam(z)*bias1(b0,z);
	
	return tmp*HH(z);
}

double Gamma10(double  b0, double z)
{
	double tmp = 2.0*pow(fg(z),2);
	return tmp*HH(z);
}



double Gamma11(double b0, double z)
{
	double fnl = 0.0;
	double tmp1 = bias1(b0,z)*(fg(z) + IC.bevo(z) - 2.0 *curlyQ(z) - 2.0*(1.0-curlyQ(z))/(chi(z)*HH(z)) - IC.Hpr(z)/pow(HH(z),2));
        double tmp2 = bias1pr(b0, z)/HH(z) + 2.0*(1.0-1.0/(chi(z)*HH(z)))*bias1Lpr(b0,z) ;
        //double tmp3 = -pow(fg(z),2);
	double tmp = fg(z)*(tmp1 + tmp2);// + tmp3;
	return tmp*HH(z);
}


double Gamma12(double  b0,double z)
{
double tmp = -1.5* Omegam(z)*fg(z);
	return tmp*HH(z);	
}

double Gamma13(double  b0, double z)
{
double tmp1 = 1.5* Omegam(z)*fg(z);

 double tmp2 = pow(fg(z),2)*(-3.0 + 2.0*IC.bevo(z) -4.0 *curlyQ(z) -4.0*(1.0-curlyQ(z))/(chi(z)*HH(z))- 3.0*IC.Hpr(z)/pow(HH(z),2));

 //double tmp3 = fg(z)*(IC.bevo(z) - 2.0*curlyQ(z) -2.0*(1.0-curlyQ(z))/(chi(z)*HH(z)) -Hpr(z)/pow(HH(z),2));

  double tmp = tmp1 + tmp2;// + tmp3;

   return tmp*HH(z);		
}

double Gamma14(double  b0, double z)
{
        double tmp1 = fg(z)*(IC.bevo(z) - 2.0*curlyQ(z)- Hpr(z)/pow(HH(z),2) -2.0*(1-curlyQ(z))/(chi(z)*HH(z)));

	return tmp1*HH(z);

}

double Firstgamma1(double z)
{
	double tmp1 = fg(z)*(IC.bevo(z) - 2.0*curlyQ(z)- Hpr(z)/pow(HH(z),2) -2.0*(1-curlyQ(z))/(chi(z)*HH(z)));
	return tmp1*HH(z);
}


double Firstgamma2(double z)
{
	double tmp1 = -fg(z)*(IC.bevo(z) -3.0);
	double tmp2A = -2.0 + fg(z) - IC.bevo(z) + 4.0* curlyQ(z)  + Hpr(z)/pow(HH(z),2)+ 2.0*(1.0-curlyQ(z))/(chi(z)*HH(z));
	double tmp2B = -1.5*Omegam(z);
	double tmp2 = tmp2A*tmp2B;
        double tmp = tmp1 + tmp2;



	return tmp*pow(HH(z),2);
	
}
//////////////////////////////////////////////////////comparison tools ends///////////////////////////////////////










/////////////////////////////////////////////Analytical Newtonian Bispectrum/////////////////////////////////////////


double CurlyB0( double b0, double z)
{
 double b1 = bias1(b0, z);
 double b2 = bias2(b0, z);
 

 double line1 =  20.0*pow(b1,3)/7.0 + 52.0*pow(b1,2)*IC.fg(z)/21.0 + 2.0*pow(b1,2)*b2 + 4.0*b1*b2*IC.fg(z)/3.0 + 68.0*b1*pow(IC.fg(z),2)/105.0 + 4.0*pow(b1,3)*IC.fg(z)/3.0 + 4.0*pow(b1,2)*pow(IC.fg(z),2)/3.0  + 12.0*b1*pow(IC.fg(z),3)/35.0;// + 2.0*pow(b1,3)*IC.fg(z)/3.0 + 2.0*pow(b1*IC.fg(z),2); 

//double line2 = 20.0*pow(b1,2)*IC.fg(z)/21.0 + 2.0*b1*b2*IC.fg(z)/3.0 + 62.0*b1*pow(IC.fg(z),2)/105.0 - pow(b1,3)*IC.fg(z) + 11.0*pow(b1*IC.fg(z),2)/60.0 + 41.0*b1*pow(IC.fg(z),3)/140.0 + 12.0*b2 *pow(IC.fg(z),2)/5.0 + 54.0*pow(IC.fg(z),3)/245.0 + 4.0*pow(IC.fg(z),4)/104.0;

double line2 = 2.0*b2 *pow(IC.fg(z),2)/15.0 + 12.0*pow(IC.fg(z),3)/245.0 + 4.0*pow(IC.fg(z),4)/105.0;

double line =  line1 + line2; 

 return line;


}


double CurlyB2(double b0, double z)
{
 double b1 = bias1(b0, z);
 double b2 = bias2(b0, z);
 

  double line1 =  80.0*b1*IC.fg(z) + 21.0*b1*pow(IC.fg(z),2) + 150.0*pow(b1,2) + 35.0*b1*b1*IC.fg(z) + 105.0*b1*b2 + 35.0*b2*IC.fg(z) + 18.0*pow(IC.fg(z),2);

  double line2 = 5.0*(3.0*b1 + IC.fg(z))*Gamma7(b0,z) + (5.0*b1 + 3.0 *IC.fg(z)) *Gamma8(b0,z);

  double line3 = 35.0*pow(b1,2) * IC.fg(z) + 21.0 *b1* pow(IC.fg(z),2) + 6.0 * pow(IC.fg(z),3);

  double line4 = (35.0* b1 + 7.0* IC.fg(z))*Gamma11(b0,z)  + ( 7.0* b1 + 3.0 * IC.fg(z) ) * Gamma13(b0,z);

 double line5 = (35.0 *b1 + 7.0*IC.fg(z) ) * Gamma7(b0,z) + (7.0*b1 + 3.0 *IC.fg(z))* Gamma8(b0,z);

 double line6 = (7.0*b1 + 3.0*IC.fg(z))*IC.fg(z)*Firstgamma1(z) + 7.0*(5.0*b1 + IC.fg(z))* Gamma9(b0,z) + 3.0*(7.0*b1 +IC.fg(z))*Gamma12(b0,z);

 // double line =  Firstgamma2(z)*line1;// + 

   double line = Firstgamma2(z)*line1 +  7.0*b1*line2  + Firstgamma2(z)*line3 - Firstgamma1(z)*line4 + IC.fg(z)*line5  - Firstgamma1(z)*line6;

//double line =  7.0*b1*line2  + Firstgamma2(z)*line3 - Firstgamma1(z)*line4 + IC.fg(z)*line5 - Firstgamma1(z)*line6;

  double tmp = 2.0* line/105;

 return tmp;


}

double CurlyB4(double b0, double z)
{
 double b1 = bias1(b0, z);
 double b2 = bias2(b0, z);
 

   double line1 =  5.0*(3.0*b1 + IC.fg(z))*Gamma7(b0,z) + (5.0*b1 + 3.0*IC.fg(z))*Gamma8(b0,z) -  Firstgamma1(z)*(5.0*Gamma9(b0,z) + 3.0*Gamma12(b0,z)); 
            
   double line =  2.0*Firstgamma2(z)*line1/15.0; 

   return line;


}

double Shortshort(double b0, double z)
{
 double b1 = bias1(b0, z);
 double b2 = bias2(b0, z);
  double tmp = 2.0*b2*(15.0 *pow(b1,2) + 10.0 * b1 *IC.fg(z) + 3.0 *pow(IC.fg(z),2))/15.0;
 return tmp;

}


double Shortlong(double b0, double z)
{

double b1 = bias1(b0, z);
 double b2 = bias2(b0, z);
  double tmp = 7.0*pow(Firstgamma1(z),2)*Gamma7(b0,z) + Firstgamma1(z)*(35.0 * b1 + 7.0*Gamma7(b0,z) );
return tmp;

}


double Spectroscopic::SquzzedbispectrumN(double b0,double kS, double z)
{
 Cosmology C;
  
 double kL =  kS/16.0;
 
 double tmp = Shortshort(b0, z)*P1P2(z, kS, kS) + CurlyB0(b0,z)*P1P2(z, kL, kS);
 return tmp/2.0;
  //return C.P_k_EH_Tabulated(k);

}


double Spectroscopic::Squzzedbispectrum(double b0,double kS, double z)
{
 Cosmology C;
  
 double kL = kS/16.0;
 
 double tmp1 = Shortshort(b0, z) + 5.0 *pow(Firstgamma1(z),2)/pow(kS,2);
 double tmp2 = (CurlyB0(b0,z) + CurlyB2(b0,z)/pow(kL,2) + CurlyB4(b0,z)/pow(kL,4))*P1P2(z, kL, kS);
  double tmp3 =  Shortlong(b0, z)*pow(kL*kS,2);

 double tmp =  tmp1*P1P2(z, kS, kS) + tmp2  + tmp3*P1P2(z, kL, kS);

//cout<<"\t" <<"tmp"<<" "<<"z="<<z<<"\t"<<"B2 "<<"\t"<<CurlyB2(b0,z)<<"\t"<<"b4 "<<"\t"<<CurlyB4(b0,z)<<endl;
 return tmp/2.0;
  //return C.P_k_EH_Tabulated(k);

}



double MatterPowerSpectrum(double z, double k)
{
      Cosmology C;
     double Dz = C.D_growth_Tabulated(z);
     double Pk=  pow(Dz,2)*C.P_k_EH_Tabulated(k); 
     return Pk;
}

/////power spectrum from camb////////////////////////////////

////Power spectrum from Class


const int N = 587; ///You have to make sure this number corresponds to the number of rows in class//cat explanatory13_pk.dat | wr -l
std::vector <double> k_array(N); // define vectors to k-values
std::vector <double> pk_array(N);  // define vectors to hold pofk


void Spectroscopic::ReadWrite(void)
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

double P1P2(double z, double  k1, double k2)
  {
    Cosmology C;
    double Dz = C.D_growth_Tabulated(z);
    double tmp =  pow(Dz,4)* MySplinePofk(k1)*MySplinePofk(k2);
   return tmp;
   }
////////////////

*/
