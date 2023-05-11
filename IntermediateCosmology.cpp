


#include "IntermediateCosmology.h"





IntermediateCosmology::IntermediateCosmology()
{
    //SetPowerSpectrum();
}

IntermediateCosmology::~IntermediateCosmology()
{
    //SetPowerSpectrum();
}



double IntermediateCosmology::HH(double z)
{
    return HH_in(z);
}
double HH_in(double z)
{
    Cosmology C;
    double c_light = C.H_0*3000.0;
    const  double c_light_km_s =  299792;
    // double c_light =  C.H_0*3000.0;
    return C.H(z)/(c_light*(1.0+z));
}
double Hpr_in(double z)
{
  //Cosmology C;
	double delta = 1.0e-4;
	double x1 = z -delta;
	double x2 = z + delta ;
	double y1 = HH_in(x1);
	double y2 = HH_in(x2);
	double deri = (y2-y1)/(x2-x1);
	double tmp = -(1.0 + z)*HH_in(z)*deri;
	return tmp;
}




double chi_in(double z)
{
     Cosmology C;
     double zmin1 = 1e-5;
     //double tmpchi = 0;
     double tmpchi = C.ComovingDistance(zmin1,z);
     return tmpchi;
}

double Omegam_in(double z)//Omegam
{
    Cosmology C;
    double tmpz = C.Omega_m*pow((1.0+z),3)/(pow(C.E(z),2));
    return tmpz;
}




double IntermediateCosmology::Alpha( double k,double z)
{
     
      Cosmology C;
     double g0_overg_infty = 3.0/4.0;
     double c_light =  C.H_0*3000.0;
    // double c_light =  C.H_0*3000.0;
    const  double c_light_km_s =  299792;
     double numerator = 2.0*k*k*c_light*c_light*C.D_growth_Tabulated(z)*C.Tk_tab(k);
     double denominator = 3.0*C.Omega_m*pow(C.H_0,2) ;

     double tmp= g0_overg_infty*numerator/denominator;

//cout <<"\t"<<"Alpha="<<tmp<<endl;
     return tmp;
}


  double IntermediateCosmology::AlphaAlpha2(double z, double  k1, double k2)
{
  return Alpha(k1,z)*Alpha(k2,z);

}






double IntermediateCosmology::bevo(double z)
{
	
	double delta_c = 1.686;
	//double tmp = delta_c*(HIB.HI_b_ST_10(z)-1.0);
	//return HIB.be_Sp(z)/HH(z);
       //return -2.1*sqrt(1.0+z);
    return -4.0;
}

double IntermediateCosmology::bevopr(double z)
{
	double delta = 1.0e-4;
	double x1 = z - delta;
	double x2 = z + delta ;
	double y1 = bevo(x1);
	double y2 = bevo(x2);
	double deri = (y2-y1)/(x2-x1);
	double tmp = -(1.0 + z)*HH(z)*deri;
	return 0.0;
}





double IntermediateCosmology::curlyQ(double z)
{
   //Fix later
        double s = -0.95;
	double tmp = 5.0*s/2.0;
	return tmp;
}

double IntermediateCosmology::curlyQpr(double z)
{
        double delta = 1.0e-4;
	double x1 = z - delta;
	double x2 = z + delta ;
	double y1 = curlyQ(x1);
	double y2 = curlyQ(x2);
	double deri = (y2-y1)/(x2-x1);
	double tmp = -(1.0 + z)*HH(z)*deri;
	return tmp;
}
double IntermediateCosmology::curlyQzpr(double z)
{
        double delta = 1.0e-4;
	double x1 = z - delta;
	double x2 = z + delta ;
	double y1 = curlyQ(x1);
	double y2 = curlyQ(x2);
	double deri = (y2-y1)/(x2-x1);
	double tmp = -(1.0 + z)*deri;
	return tmp;
}
double IntermediateCosmology::curlyQLpr(double  z)
{
 double tmp = 0.0;
return tmp;
}



double f22(double z)
{
  Cosmology C;
  double tmpz = C.Omega_m*pow((1.0+z),3)/(pow(C.E(z),2));
   return pow(tmpz,0.55);
}

int k_bins2 =200;
std::vector <double> z_Dz_tab2(k_bins2);
std::vector <double> chi_tab(k_bins2);
std::vector <double> f_tab(k_bins2);
std::vector <double> ISW_tab(k_bins2);
std::vector <double> OM_tab(k_bins2);
std::vector <double> jz_tab(k_bins2);
std::vector <double> Tb_tab2(k_bins2);

std::vector <double> k_arrayin(k_bins2);


std::vector <double> Hpr_tab(k_bins2);
 std::vector <double> Sig_tab1(k_bins2);
 std::vector <double> Sigv_tab1(k_bins2);



void IntermediateCosmology::InitialzieSpeed4(void)
{
	double  EPSREL = 1e-8;
	double z_Dz_min2;
	double z_Dz_max2;
	z_Dz_min2 = log(1e-5);
	z_Dz_max2 = log(10.0);
	const double QMIN = 1e-4;
	const double QMAX = 1e5;
	
	
	
	
#pragma omp parallel private(thread_number)
	{
#pragma omp for schedule(static) nowait
		for(int i = 0;i< k_bins2 ;i++)
		{
			z_Dz_tab2[i] = ((z_Dz_max2-z_Dz_min2)*((double) i)/((double)  k_bins2-1) + z_Dz_min2 );
			//k_arrayin[i] = exp( ((log(QMAX)-log(QMIN))*((double) i)/((double) k_bins2-1) + log(QMIN) ));
			Hpr_tab[i] = Hpr_in(exp(z_Dz_tab2[i]));
			f_tab[i]  = f22(exp(z_Dz_tab2[i]));
			chi_tab[i]  =   chi_in(exp(z_Dz_tab2[i]));
                        //Sig_tab1[i] =  Mysigma1(exp(z_Dz_tab2[i]));
                        //Sigv_tab1[i] =  Mysigmav(exp(z_Dz_tab2[i]));
                       // Tb_tab2[i] = HIB.T_bST(exp(z_Dz_tab2[i]));
			OM_tab[i]  =   Omegam_in(exp(z_Dz_tab2[i]));//Mysigma2(int x,double fnl,double R,double z)
			
		}
		
	}
}


double IntermediateCosmology::fg(double z)
{
	//Spline<double, double> CubicSpline_f(z_Dz_tab2,f_tab);
	//return CubicSpline_f.interpolate(log(z));
   	return f22(z);
}


double IntermediateCosmology::chi(double z)
{
	Spline<double, double> CubicSpline_chi(z_Dz_tab2,chi_tab);
	return CubicSpline_chi.interpolate(log(z));
   	
}
double IntermediateCosmology::Omegam(double z)
{
	Spline<double, double> CubicSpline_OM(z_Dz_tab2,OM_tab);
	return CubicSpline_OM.interpolate(log(z));
   	
}

double IntermediateCosmology::Hpr(double z)
{
	Spline<double, double> CubicSpline_Hpr(z_Dz_tab2,Hpr_tab);
	return CubicSpline_Hpr.interpolate(log(z));
   	
}

double IntermediateCosmology::Omegampr(double z)
{
       double delta = 1.0e-4;
	double x1 = z -delta;
	double x2 = z + delta ;
	double y1 = Omegam_in(x1);
	double y2 = Omegam_in(x2);
	double deri = (y2-y1)/(x2-x1);
	double tmp = -(1.0 + z)*HH(z)*deri;
	return tmp;
}


double IntermediateCosmology::Hpr2(double z)
{
  //Cosmology C;
	double delta = 1.0e-4;
	double x1 = z -delta;
	double x2 = z + delta ;
	double y1 = Hpr_in(x1);
	double y2 = Hpr_in(x2);
	double deri = (y2-y1)/(x2-x1);
	double tmp = -(1.0 + z)*HH(z)*deri;
	return tmp;
}










double IntermediateCosmology::fgpr(double z)
{
   double delta = 1.0e-4;
	double x1 = z -delta;
	double x2 = z + delta ;
	double y1 = f22(x1);
	double y2 = f22(x2);
	double deri = (y2-y1)/(x2-x1);
	double tmp = -(1.0 + z)*HH(z)*deri;
	return tmp;
}



double IntermediateCosmology::bias1(double fnl, double z)
{
	
	return 1.3;//sqrt(1.0+z);
}
double IntermediateCosmology::TidalBias(double z)
 {  
   
    double fnl = 0.0;
    double HIb10 = bias1(fnl,z);//HIB.HI_b_ST_10(z);
    double tmp = -2.0*(HIb10-1.0)/7.0;
    return tmp;
 }
double IntermediateCosmology::bias1pr(double fnl, double z)
{
	double delta = 1.0e-4;
	double x1 = z -delta;
	double x2 = z + delta ;
	double y1 = bias1(fnl,x1);
	double y2 = bias1(fnl,x2);
	double deri = (y2-y1)/(x2-x1);
	double tmp = -(1.0 + z)*HH(z)*deri;
	return -1.6e-4; //tmp;
}

double IntermediateCosmology::bias2(double fnl, double z)
{
	
    double b0 = 1.0;
    double b1 = sqrt(1.0+z);
      //double tmp = -0.25*b0 - 0.13*sqrt(1.0+z)*sin(0.8*z);
     // double tmp = 0.412 -2.143*b1  + 0.929*pow(b1,2) + 0.008*pow(b1 ,3);
     //double tmp = -0.1*sqrt(1.0+z);
     return -0.74;//-0.321;//tmp;
}


   double IntermediateCosmology::b_one(double fnl, double k, double z)
{
    
    Cosmology C;
    double delta_c = 1.686;
    double HIb10 = bias1(1.0,z);
    double DeltaB = 2.0* fnl *delta_c*(HIb10-1.0)/Alpha(k, z);
  // double tmp  = HIb10+ (HIb01/alpha_z(k, z));
    double tmp = HIb10 + DeltaB;
   return  tmp;
}


double IntermediateCosmology::b_two_integrand(double fnl, double z ,double k1,double k2,double mu1)
{   
   
     Cosmology C;
     double delta_c =  1.686;
   
    
   
     double HIb10 = bias1(1.0,z);  
     double HIb20 = bias2(1.0,z); 
     double Eulerian = HIb10-1.0;

     double innerbrac1 = delta_c*HIb20 + ((13.0*delta_c/21.0) -1.0)*Eulerian;
     double innerbrac2 =delta_c*HIb20-2.0*((4.0*delta_c/21.0) +1.0)*Eulerian;

     double DeltaSqB_Ist1 = fnl*innerbrac1/Alpha(k1, z);
    double DeltaSqB_Ist2 = fnl*innerbrac1/Alpha(k2, z);
     double DeltaSqB_2nd = 4.0*fnl*fnl*delta_c*innerbrac2/AlphaAlpha2(z,  k1, k2);
     double tmp =  HIb20 + (DeltaSqB_Ist1+DeltaSqB_Ist2) + DeltaSqB_2nd;

     return tmp;

} 

double IntermediateCosmology::bias1Lpr(double fnl, double z)
{
	
	return 0.0;
}
double IntermediateCosmology::bevoLpr(double z)
{
	
	return 0.0;
}



double IntermediateCosmology::mu23(double r, double mu12)
{
     ///double r = k2/k1;
     double y = sqrt(r*r+2.0*mu12 *r + 1.0);
     
     double tmp = -(r + mu12)/(y);
     return tmp;

}

double IntermediateCosmology::mu13(double r, double mu12)
{
     ///double r = k2/k1;
     double y = sqrt(r*r+2.0*mu12 *r + 1.0);
     double tmp =- (r* mu12 +1.0)/(y);
     return tmp;

}

