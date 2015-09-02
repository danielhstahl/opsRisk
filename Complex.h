#ifndef __COMPLEX_H_INCLUDED__
#define __COMPLEX_H_INCLUDED__

#include <cmath>
#include <iostream> //for debugging

class Complex {
	private:
		double real;
		double im;
	public:
		Complex();
		Complex(double, double);
		Complex multiply(double, double);
		Complex multiply(Complex);
		Complex multiply(double);
		Complex exp();
		Complex log();
		Complex add(Complex);
		Complex add(double);
		Complex subtract(Complex);
		Complex subtract(double);
		Complex power(double);
		double getReal();
		double getIm();
		Complex divide(Complex);
		Complex divide(double);
};

#endif