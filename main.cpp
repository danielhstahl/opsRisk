#define _USE_MATH_DEFINES
#include <iostream>
#include "FangOosterlee.h"
#include "RungeKutta.h"
#include "Complex.h"
#include <cmath>
#include <fstream>
#include <ctime>
//#define pi=3.14159265
Complex GaussCF(Complex u, std::map<std::string, double> params){
	return u.multiply(params["mu"]).add(u.multiply(u).multiply(params["sigma"]*params["sigma"]*.5)).exp();
}

Complex stableCF(Complex u, std::map<std::string, double> params){
	double alpha=params["alpha"];
	double phi=tan(alpha*.5*M_PI);
	Complex uC=u.multiply(Complex(0, -1));
	return u.multiply(params["mu"]).subtract(uC.multiply(params["c"]).power(alpha).multiply(Complex(1, -params["beta"]*phi))).exp();
}
std::vector<Complex> DuffieODE(double t, std::vector<Complex> initialValues, std::map<std::string, double> params, std::map<std::string, Complex> cmplParams){
	std::vector<Complex> vls(2);
	double sig=params["sigma"];
	sig=sig*sig*.5;
	double lambda=params["lambda"];
	double a=params["a"];
	//Complex uBeta=cmplParams["u"].add(initialValues[0].multiply(params["delta"])).multiply(Complex(0, -1));
	Complex uBeta=cmplParams["u"].add(initialValues[0].multiply(params["delta"]));
	vls[0]=initialValues[0].multiply(initialValues[0]).multiply(sig).add(stableCF(uBeta, params).multiply(lambda)).subtract(lambda).subtract(initialValues[0].multiply(a));
	vls[1]=initialValues[0].multiply(params["b"]*a);
	return vls;
}
Complex distToInvert(Complex u, std::map<std::string, double> params){
	std::vector<Complex> inits(2); //fortunately these will never change...
	inits[0]=Complex(0, 0);
	inits[1]=Complex(0, 0);
	std::map<std::string, Complex> cmplParam;
	cmplParam["u"]=u;
	RungeKutta rg(DuffieODE, params["t"], inits, (int)params["numODE"]);
	std::vector<Complex> vls=rg.compute(params, cmplParam);
	return vls[0].multiply(params["v0"]).add(vls[1]).exp();
}
int main(){
	FangOosterlee invert(256, 1024);
	std::map<std::string, double> params; //im not a huge fan of this but it works
	params["mu"]=1300;
	params["c"]=100;
	params["alpha"]=1.1;
	params["beta"]=1;
	params["lambda"]=100;
	double rho=.5;
	params["delta"]=rho/(params["mu"]*params["lambda"]);
	params["a"]=.4;
	params["sigma"]=.4;
	params["b"]=1-rho;
	params["t"]=1;
	params["numODE"]=128;
	params["v0"]=1;
	double xmin=-100;
	double xmax=params["lambda"]*(params["mu"]+35*params["c"]);
	//double xmax=(params["mu"]+35*params["c"]);


	clock_t t;
	t = clock();
	//std::map<std::string, std::vector<double> > results=invert.computeDistribution(distToInvert, params, xmin, xmax);
	std::map<std::string, std::vector<double> > results=invert.computeDistribution(xmin, xmax, distToInvert, params);
	t = clock() - t;
	std::cout<<"Time it took: "<<(float)t/CLOCKS_PER_SEC<<std::endl;
	std::ofstream outputCSV;
	outputCSV.open("stable.csv");
	outputCSV <<"x, y"<<std::endl;
	for(int i=0; i<1024; i++){
		outputCSV<<""<<results["x"][i]<<", "<<results["y"][i]<<std::endl;
	}
	outputCSV.close();

}
