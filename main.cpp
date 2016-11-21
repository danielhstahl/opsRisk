#define _USE_MATH_DEFINES
#include <iostream>
#include "Complex.h"
#include "FangOosterlee.h"
#include "RungeKutta.h"
#include <cmath>
#include <ctime>
#include <chrono> //for accurate multithreading time using std::chrono
#include "document.h" //rapidjson
#include "writer.h" //rapidjson
#include "stringbuffer.h" //rapidjson 
Complex GaussCF(const Complex& u, double mu, double sigma){
	return exp(u*mu+u*u*sigma*sigma*.5);
}
Complex stableCF(const Complex& u, double alpha, double mu, double beta, double c){
	double phi=tan(alpha*.5*M_PI);
	return exp(u*mu-(u*Complex(0, -1)*c).pow(alpha)*Complex(1, -beta*phi));
}
Complex gammaCF(const Complex& u, double a, double b){
	return pow(1-u*b, -a);
}
Complex inverseGaussianCF(const Complex& u, double mu, double lambda){
	return exp((lambda/mu)*(1-sqrt(1-(2*mu*mu*u)/lambda)));
}
template<typename CF>
std::vector<Complex> DuffieODE(const Complex& u, const CF& cf, const std::vector<Complex>& initialValues, double sigma, double lambda, double a, double delta, double b){ //double alpha, double mu, double beta, double c,
	double sig=sigma*sigma*.5;
	return 
	{
		initialValues[0]*initialValues[0]*sig+cf(u+initialValues[0]*delta)*lambda-lambda-initialValues[0]*a, 
		initialValues[0]*b*a
	};
}
template<typename CF>
Complex distToInvert(
	const Complex& u, 
	double t, 
	int numODE, 
	double lambda, 
	double sigma, 
	double a, 
	double b, 
	double delta,
	double v0, 
	const CF &cf){
	RungeKutta rg(t, numODE);
	std::vector<Complex> inits={Complex(0, 0), Complex(0, 0)};
	inits=rg.compute(
		[&](double t, const std::vector<Complex>& x){
			return DuffieODE(u, cf, x, sigma, lambda, a, delta, b);
		},
		std::move(inits)
	);
	return exp(inits[0]*v0+inits[1]);
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
	std::vector<Complex> inits(2);
	auto runParameters=[&](std::string& parameters){
		rapidjson::Document parms;
		parms.Parse(parameters.c_str());//yield data
		parameters.clear();
		if(parms.FindMember("xNum")!=parms.MemberEnd()){
			xNum=parms["xNum"].GetInt();
		}
		if(parms.FindMember("uNum")!=parms.MemberEnd()){
			uNum=parms["uNum"].GetInt();
		}
		if(parms.FindMember("numODE")!=parms.MemberEnd()){
			numODE=parms["numODE"].GetInt();
		}
		if(parms.FindMember("muStable")!=parms.MemberEnd()){
			muStable=parms["muStable"].GetDouble();
		}
		if(parms.FindMember("cStable")!=parms.MemberEnd()){
			cStable=parms["cStable"].GetDouble();
		}
		if(parms.FindMember("alphaStable")!=parms.MemberEnd()){
			alphaStable=parms["alphaStable"].GetDouble();
		}
		if(parms.FindMember("lambda")!=parms.MemberEnd()){
			lambda=parms["lambda"].GetDouble();
		}
		if(parms.FindMember("rho")!=parms.MemberEnd()){
			rho=parms["rho"].GetDouble();
		}
		if(parms.FindMember("v0")!=parms.MemberEnd()){
			v0=parms["v0"].GetDouble();
		}
		if(parms.FindMember("t")!=parms.MemberEnd()){
			t=parms["t"].GetDouble();
		}
		if(parms.FindMember("a")!=parms.MemberEnd()){
			a=parms["a"].GetDouble();
		}
		if(parms.FindMember("sigma")!=parms.MemberEnd()){
			sigma=parms["sigma"].GetDouble();
		}
		delta=rho/(muStable*lambda); //they all have the same expected value
		b=1-rho;
		xmax=lambda*(muStable+35*cStable);
		FangOosterlee invert(uNum, xNum);
		//auto start = std::chrono::system_clock::now();
		invert.computeDistributionJSON(xmin, xmax,
			[&](Complex &u){
				return distToInvert(u, t, numODE, lambda, sigma, a, b, delta, v0,
					[&](const Complex& uhat){
						return stableCF(uhat, alphaStable, muStable, betaStable, cStable);
					});
				});
	 };
	 runParameters(std::string("{}"));
	 /*while(true){
			 std::string parameters;
			 for (parameters; std::getline(std::cin, parameters);) {
					 break;
			 }
			 runParameters(parameters);
	 }*/
}
