#define _USE_MATH_DEFINES
#include <iostream>
#include "FangOosterlee.h"
#include "RungeKutta.h"
#include "Complex.h"
#include <cmath>
#include <ctime>
#include <chrono> //for accurate multithreading time using std::chrono
#include "document.h" //rapidjson
#include "writer.h" //rapidjson
#include "stringbuffer.h" //rapidjson 
Complex GaussCF(Complex u, std::unordered_map<std::string, double> params){
	return u.multiply(params["mu"]).add(u.multiply(u).multiply(params["sigma"]*params["sigma"]*.5)).exp();
}

Complex stableCF(Complex &u, double alpha, double mu, double beta, double c){
	double phi=tan(alpha*.5*M_PI);
	Complex uC=u.multiply(Complex(0, -1));
	return u.multiply(mu).subtract(uC.multiply(c).pow(alpha).multiply(Complex(1, -beta*phi))).exp();
}
Complex gammaCF(Complex &u, double a, double b){
	return pow(1-u*b, -a);
}
Complex inverseGaussianCF(Complex &u, double mu, double lambda){
	return exp((lambda/mu)*(1-sqrt(1-(2*mu*mu*u)/lambda)));
}
template<typename CF>
std::vector<Complex> DuffieODE(Complex &u, CF &cf, std::vector<Complex> &initialValues, double sigma, double lambda, double a, double delta, double b){ //double alpha, double mu, double beta, double c,
	std::vector<Complex> vls(2);
	double sig=sigma*sigma*.5;
	Complex uBeta=u.add(initialValues[0].multiply(delta));

	vls[0]=initialValues[0].multiply(initialValues[0]).multiply(sig).add(cf(uBeta).multiply(lambda)).subtract(lambda).subtract(initialValues[0].multiply(a));
	vls[1]=initialValues[0].multiply(b*a);
	return vls;
}
template<typename CF>
Complex distToInvert(Complex &u, std::unordered_map<std::string, double> &params, std::vector<Complex> &inits, const CF &cf){
	RungeKutta rg(params["t"], (int)params["numODE"]);
	std::vector<Complex> vls=rg.compute([&](double t, std::vector<Complex> &x){
		return DuffieODE(u, cf, x, params["sigma"], params["lambda"], params["a"], params["delta"], params["b"]);
	}, inits);
	return vls[0].multiply(params["v0"]).add(vls[1]).exp();
}

int main(){
	int xNum=1024;//standard parameters
	int uNum=256;//standard parameters
	double muStable=1300;
	double cStable=100;
	double alphaStable=1.1;
	double betaStable=1;
	std::unordered_map<std::string, double> params; //im not a huge fan of this but it works
	params["lambda"]=100;
	double rho=.9;
	params["delta"]=rho/(muStable*params["lambda"]); //they all have the same expected value
	params["b"]=1-rho;
	params["a"]=.4;
	params["sigma"]=.4;
	params["t"]=1;
	params["numODE"]=128;
	params["v0"]=1;
	double xmin=0;
	double xmax=params["lambda"]*(muStable+35*cStable);
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
			params["numODE"]=parms["numODE"].GetInt();
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
			params["lambda"]=parms["lambda"].GetDouble();
		}
		if(parms.FindMember("rho")!=parms.MemberEnd()){
			rho=parms["rho"].GetDouble();
		}
		if(parms.FindMember("v0")!=parms.MemberEnd()){
			params["v0"]=parms["v0"].GetDouble();
		}
		if(parms.FindMember("t")!=parms.MemberEnd()){
			params["t"]=parms["t"].GetDouble();
		}
		if(parms.FindMember("a")!=parms.MemberEnd()){
			params["a"]=parms["a"].GetDouble();
		}
		if(parms.FindMember("sigma")!=parms.MemberEnd()){
			params["sigma"]=parms["sigma"].GetDouble();
		}
		params["delta"]=rho/(muStable*params["lambda"]); //they all have the same expected value
		params["b"]=1-rho;
		xmax=params["lambda"]*(muStable+35*cStable);
		FangOosterlee invert(uNum, xNum);
		//auto start = std::chrono::system_clock::now();
		invert.computeDistributionJSON(xmin, xmax,
			[&](Complex &u, std::unordered_map<std::string, double> &parms){
				inits[0]=Complex(0, 0);
				inits[1]=Complex(0, 0);
				return distToInvert(u, parms, inits,
					[&](Complex &uhat){
						return stableCF(uhat, alphaStable, muStable, betaStable, cStable);
					});
				}, params);
	 };
	 while(true){
			 std::string parameters;
			 for (parameters; std::getline(std::cin, parameters);) {
					 break;
			 }
			 runParameters(parameters);
	 }
}
