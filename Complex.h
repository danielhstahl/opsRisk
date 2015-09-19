#ifndef __COMPLEX_H_INCLUDED__
#define __COMPLEX_H_INCLUDED__

#include <cmath>
#include <iostream> //for debugging

class Complex {
	private:
		double real;
		double im;
	public:

		friend Complex operator+(const Complex &c1, const Complex &c2) ;
		friend Complex operator+(const double c1, const Complex &c2) ;
		friend Complex operator+(const Complex &c1, double c2) ;

		friend Complex operator-(const Complex &c1, const Complex &c2) ;
		friend Complex operator-(double c1, const Complex &c2) ;
		friend Complex operator-(const Complex &c1, double c2) ;
		friend Complex operator/(const Complex &c1, const Complex &c2) ;
		friend Complex operator/(double c1, const Complex &c2) ;
		friend Complex operator/(const Complex &c1, double c2) ;
		friend Complex operator*(const Complex &c1, const Complex &c2) ;
		friend Complex operator*(double c1, const Complex &c2);
		friend Complex operator*(const Complex &c1, double c2);
		Complex();
		Complex(double, double);
		Complex multiply(double, double) const;
		Complex multiply(const Complex&) const;
		Complex multiply(double) const;
		Complex exp();
		Complex log();
		Complex add(const Complex&) const;
		Complex add(double) const;
		Complex subtract(const Complex&) const;
		Complex subtract(double) const;
		Complex power(double) const;
		Complex divide(const Complex&) const;
		Complex divide(double) const;
	/*	Complex Complex::multiply(double x, double y) const { //consider removing
			//Complex plac=Complex(real*x-im*y, im*x+real*y);
			return Complex(real*x-im*y, im*x+real*y);
		}*/
		double getReal() const;
		double getIm() const;

};



#endif
