#define _USE_MATH_DEFINES
#include <iostream>
#include <complex>
#include "FangOost.h"
#include "RungeKutta.h"
#include <fstream>
#include <cmath>
#include <ctime>
#include <chrono> //for accurate multithreading time using std::chrono
/*#include "document.h" //rapidjson
#include "writer.h" //rapidjson
#include "stringbuffer.h" //rapidjson */
template<typename T>
std::complex<T> GaussCF(const std::complex<T>& u, double mu, double sigma){
	return exp(u*mu+u*u*sigma*sigma*.5);
}
template<typename T>
std::complex<T> stableCF(const std::complex<T>& u, double alpha, double mu, double beta, double c){
	double phi=tan(alpha*.5*M_PI);
	return exp(u*mu-pow(u*std::complex<T>(0, -1)*c, alpha)*std::complex<T>(1, -beta*phi));
}
template<typename T>
std::complex<T> gammaCF(const std::complex<T>& u, double a, double b){
	return pow(1-u*b, -a);
}
template<typename T>
std::complex<T> inverseGaussianCF(const std::complex<T>& u, double mu, double lambda){
	return exp((lambda/mu)*(1-sqrt(1-(2*mu*mu*u)/lambda)));
}
template<typename T>
std::vector<std::complex<T> > DuffieODE(const auto& u, auto&& cf, const std::vector<std::complex<T> >& initialValues, double sigma, double lambda, double a, double delta, double b){ //double alpha, double mu, double beta, double c,
	double sig=sigma*sigma*.5;
	return 
	{
		initialValues[0]*initialValues[0]*sig+cf(u+initialValues[0]*delta)*lambda-lambda-initialValues[0]*a, 
		initialValues[0]*b*a
	};
}
template<typename T>
std::complex<T> expAffine(const std::vector<std::complex<T>>& vals, double v0){
	return exp(vals[0]*v0+vals[1]);
}
template<typename T>
std::complex<T> distToInvert(
	const std::complex<T>& u, 
	double t, 
	int numODE, 
	double lambda, 
	double sigma, 
	double a, 
	double b, 
	double delta,
	double v0, 
	auto&& cf
){
	return expAffine(
		/*rungekutta::compute(t, numODE, std::vector<std::complex<T> >({std::complex<T>(0, 0), std::complex<T>(0, 0)}),
			[&](double t, const std::vector<std::complex<T>>& x){
				return DuffieODE(u, cf, x, sigma, lambda, a, delta, b);
			}
		), */
		rungekutta::computeFunctional(t, numODE, std::vector<std::complex<T> >({std::complex<T>(0, 0), std::complex<T>(0, 0)}),
			[&](double t, const std::vector<std::complex<T>>& x){
				return DuffieODE(u, cf, x, sigma, lambda, a, delta, b);
			}
		),
		v0
	);
}
int main(){
	int xNum=1024;//standard parameters
	int uNum=256;//standard parameters
	double muStable=1300;
	double cStable=100;
	double alphaStable=1.1;
	double betaStable=1;
	double lambda=100;
	double rho=.9;
	double delta=rho/(muStable*lambda); //they all have the same expected value
	double b=1-rho;
	double a=.4;
	double sigma=.4;
	double t=1;
	int numODE=128;
	double v0=1;
	double xmin=0;
	double xmax=lambda*(muStable+35*cStable);
	//std::vector<std::complex<double>> inits(2);
	auto density=fangoost::computeInv(xNum, uNum, xmin, xmax, [&](const auto& u){
		return distToInvert(u, t, numODE, lambda, sigma, a, b, delta, v0, [&](const auto& uhat){
			return stableCF(uhat, alphaStable, muStable, betaStable, cStable);
		});
	});
	auto axis=fangoost::computeXRange(xNum, xmin, xmax);

	std::ofstream out("output.csv");
	for(int i=0; i<xNum; ++i){
		out<<axis[i]<<", "<<density[i]<<std::endl;
	}
	out.close();
}
