#include "Complex.h"

Complex::Complex() {
	real=0;
	im=0;
}
Complex::Complex(double x, double y) {
	real=x;
	im=y;
}
Complex Complex::multiply(const Complex &c) const {
	//Complex plac=Complex(real*c.getReal()-im*c.getIm(), im*c.getReal()+real*c.getIm());
	return Complex(real*c.getReal()-im*c.getIm(), im*c.getReal()+real*c.getIm());
}
Complex Complex::multiply(double c) const {
	//Complex plac=Complex(c*real, c*im);
	return Complex(c*real, c*im);
}
Complex Complex::exp() {
	//Complex plac=Complex(std::exp(real)*std::cos(im), std::exp(real)*std::sin(im));
	return Complex(std::exp(real)*std::cos(im), std::exp(real)*std::sin(im));
}
Complex Complex::log() {
	double modulus=std::sqrt(real*real+im*im);
	//Complex plac=Complex(std::log(modulus), std::atan2(im, real));
	return Complex(std::log(modulus), std::atan2(im, real));
}
Complex Complex::add(const Complex &c) const {
	//Complex plac=Complex(real+c.getReal(), im+c.getIm());
	return Complex(real+c.getReal(), im+c.getIm());
}
//template <typename T>
Complex Complex::add(double c) const {
	//Complex plac=Complex(real+c, im);
	return Complex(real+c, im);
}
//template <typename T>
Complex Complex::subtract(double c) const {
	//Complex plac=Complex(real-c, im);
	return Complex(real-c, im);
}
Complex Complex::subtract(const Complex &c) const {
	//Complex plac=Complex(real-c.getReal(), im-c.getIm());
	return Complex(real-c.getReal(), im-c.getIm());
}
Complex Complex::divide(const Complex &c) const {
	double den=c.getReal()*c.getReal()+c.getIm()*c.getIm();
	//Complex plac=Complex((real*c.getReal()+im*c.getIm())/den, (im*c.getReal()-real*c.getIm())/den);
	return Complex((real*c.getReal()+im*c.getIm())/den, (im*c.getReal()-real*c.getIm())/den);
}
Complex Complex::divide(double c) const {
	return Complex(real/c, im/c);
}
Complex Complex::power(double exponent) const {
	double modulus=std::sqrt(real*real+im*im);
	double arg=std::atan2(im, real);
	double log_re=std::log(modulus);
	double log_im=arg;
	double x_log_re=exponent*log_re;
	double x_log_im=exponent*log_im;
	double modulus_ans=std::exp(x_log_re);
	return Complex(modulus_ans*std::cos(x_log_im), modulus_ans*std::sin(x_log_im));
}

double Complex::getReal() const {
	return real;
}
double Complex::getIm() const{
	return im;
}
Complex operator+(const Complex &c1, const Complex &c2)
{
    return c1.add(c2);
}
Complex operator+(double c1, const Complex &c2)
{
    return c2.add(c1);
}
Complex operator+(const Complex &c1, double c2)
{
    return c1.add(c2);
}

Complex operator-(const Complex &c1, const Complex &c2)
{
    return c1.subtract(c2);
}
Complex operator-(double c1, const Complex &c2)
{
    return Complex(c1-c2.getReal(), -c2.getIm());
}
Complex operator-(const Complex &c1, double c2)
{
    return c1.subtract(c2);
}
Complex operator/(const Complex &c1, const Complex &c2)
{
    return c1.divide(c2);
}
Complex operator/(double c1, const Complex &c2)
{
    return c2.power(-1).multiply(c1);
}
Complex operator/(const Complex &c1, double c2)
{
    return c1.divide(c2);
}

Complex operator*(const Complex &c1, const Complex &c2)
{
    return c1.multiply(c2);
}
Complex operator*(double c1, const Complex &c2)
{
    return c2.multiply(c1);
}
Complex operator*(const Complex &c1, double c2)
{
    return c1.multiply(c2);
}
