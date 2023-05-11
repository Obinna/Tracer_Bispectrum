#include<iostream>
#include<cmath>
//#include <boost/bind.hpp>
//#include "Common.h"
//#include "Common.cpp"
using boost::cref;
using boost::bind;

using namespace std;

//const double PI = 3.141592654;
//int thread_number = 8;

double fs(double x)
{
  return pow(cos(x), 2);
}

double integrate_S(Spline<double, double> f, double a, double b, int n = 50)
{
  double integral = 0.0;
  double lowerLimit = a, upperLimit = b;
  double h = (b - a) / static_cast<double>(n);

//#pragma omp parallel for
  for(int i = 0; i < n; i++)
    {

      upperLimit = lowerLimit + h;

      integral += ((upperLimit - lowerLimit) / 8.0) * (f.interpolate(lowerLimit) + 3.0 * f.interpolate((2.0*lowerLimit + upperLimit) / 3.0) + 
							f.interpolate(upperLimit) + 3.0 * f.interpolate((2.0*upperLimit + lowerLimit) / 3.0)
							);
      lowerLimit = upperLimit;	
    }

  return integral;
}

const int N_itr = 50;
double integrate1(double(*f)(double), double a, double b, int n = N_itr)
{
  double integral = 0.0;
  double lowerLimit = a, upperLimit = b;
  double h = (b - a) / static_cast<double>(n);
//cout<<"Calling integration routine"<<"n="<<n<<"Method used="<< "Simpson"<<endl;

#pragma omp parallel private(thread_number)
{
 #pragma omp for schedule(static) nowait
//#pragma omp parallel for
  for(int i = 0; i < n; i++)
    {

      upperLimit = lowerLimit + h;

      integral += ((upperLimit - lowerLimit) / 8.0) * (f(lowerLimit) + 3.0 * f((2.0*lowerLimit + upperLimit) / 3.0) + 
							f(upperLimit) + 3.0 * f((2.0*upperLimit + lowerLimit) / 3.0)
							);
      lowerLimit = upperLimit;	
    }
}
  return integral;
}

/*
double integrate(double f, double a, double b, int n = 10)
{

  double integral = 0.0;
  double lowerLimit = a, upperLimit = b;
  double h = (b - a) / static_cast<double>(n);


  for(int i = 0; i < n; i++)
    {

      upperLimit = lowerLimit + h;

      integral += ((upperLimit - lowerLimit) / 8.0) * (f(lowerLimit) + 3.0 * f((2.0*lowerLimit + upperLimit) / 3.0) + 
							f(upperLimit) + 3.0 * f((2.0*upperLimit + lowerLimit) / 3.0)
							);
      lowerLimit = upperLimit;	
    }

  return integral;
}*/
template<typename Function>
double integrate1(Function f, double a, double b, int n = N_itr)
{
  double integral = 0.0;
  double lowerLimit = a, upperLimit = b;
  double h = (b - a) / static_cast<double>(n);
//cout<<"Calling integration routine"<<"n="<<n<<"Method used="<< "Simpson"<<endl;

#pragma omp parallel private(thread_number)
{
 #pragma omp for schedule(static) nowait
//#pragma omp parallel for
  for(int i = 0; i < n; i++)
    {

      upperLimit = lowerLimit + h;

      integral += ((upperLimit - lowerLimit) / 8.0) * (f(lowerLimit) + 3.0 * f((2.0*lowerLimit + upperLimit) / 3.0) + 
							f(upperLimit) + 3.0 * f((2.0*upperLimit + lowerLimit) / 3.0)
							);
      lowerLimit = upperLimit;	
    }
}
  return integral;
}


double integrate2(double(*f)(double), double a, double b, int n = 50)
{

  double h = (b - a) / static_cast<double>(n);

  int i = 0, j = 1;
  double u = f(a) + f(b);
  double v = 0.0;
  
  for(;;)
    {
      u += 4.0 * f(a + i*h);
      v += 2.0 * f(a + j*h);

      i += 2;
      j += 2;

      if(i >= n - 1 && j >= n - 2)
	{
	  break;
	}
    }


  return ((u + v) * h / 3.0);
}

template<typename Function>
double integrate2(Function f, double a, double b, int n = 50)
{

  double h = (b - a) / static_cast<double>(n);

  int i = 0, j = 1;
  double u = f(a) + f(b);
  double v = 0.0;
  
  for(;;)
    {
      u += 4.0 * f(a + i*h);
      v += 2.0 * f(a + j*h);

      i += 2;
      j += 2;

      if(i >= n - 1 && j >= n - 2)
	{
	  break;
	}
    }


  return ((u + v) * h / 3.0);
}


/*
int main()
{int l =2;
  double a = 0.05;
  double b =  PI;

  cout << integrate(bind(jl,l,_1), a, b) << endl;
  cout << integrate2(fs, a, b) << endl;

  return 0;
}*/


