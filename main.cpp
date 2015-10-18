#define _USE_MATH_DEFINES
#include <iostream>
#include "FangOosterlee.h"
#include "RungeKutta.h"
#include "Complex.h"
#include <cmath>
#include <fstream>
#include <ctime>
#include <chrono> //for accurate multithreading time using std::chrono
Complex GaussCF(Complex u, std::unordered_map<std::string, double> params){
	return u.multiply(params["mu"]).add(u.multiply(u).multiply(params["sigma"]*params["sigma"]*.5)).exp();
}

Complex stableCF(Complex &u, double alpha, double mu, double beta, double c){
	double phi=tan(alpha*.5*M_PI);
	Complex uC=u.multiply(Complex(0, -1));
	return u.multiply(mu).subtract(uC.multiply(c).power(alpha).multiply(Complex(1, -beta*phi))).exp();
}
std::vector<Complex> DuffieODE(double t, std::vector<Complex> &initialValues, double sigma, double lambda, double a, double delta, double b, double alpha, double mu, double beta, double c, Complex &u){
	std::vector<Complex> vls(2);
	double sig=sigma*sigma*.5;
	Complex uBeta=u.add(initialValues[0].multiply(delta));
	vls[0]=initialValues[0].multiply(initialValues[0]).multiply(sig).add(stableCF(uBeta, alpha, mu, beta, c).multiply(lambda)).subtract(lambda).subtract(initialValues[0].multiply(a));
	vls[1]=initialValues[0].multiply(b*a);
	return vls;
}
Complex distToInvert(Complex &u, std::unordered_map<std::string, double> &params, std::vector<Complex> &inits){
	RungeKutta rg(params["t"], (int)params["numODE"]);
	std::vector<Complex> vls=rg.compute(DuffieODE, inits, params["sigma"], params["lambda"], params["a"], params["delta"], params["b"], params["alpha"], params["mu"], params["beta"], params["c"],  u);
	return vls[0].multiply(params["v0"]).add(vls[1]).exp();
}
int main(){
	int xNum=4096;
	FangOosterlee invert(512, xNum);
	std::unordered_map<std::string, double> params; //im not a huge fan of this but it works
	params["mu"]=1300;
	params["c"]=100;
	params["alpha"]=1.1;
	params["beta"]=1;
	params["lambda"]=100;
	double rho=.9;
	params["delta"]=rho/(params["mu"]*params["lambda"]);
	params["a"]=.4;
	params["sigma"]=.4;
	params["b"]=1-rho;
	params["t"]=1;
	params["numODE"]=128;
	params["v0"]=1;
	double xmin=-100;
	double xmax=params["lambda"]*(params["mu"]+35*params["c"]);
	std::vector<Complex> inits(2);
	inits[0]=Complex(0, 0);
	inits[1]=Complex(0, 0);
	auto start = std::chrono::system_clock::now();
	std::unordered_map<std::string, std::vector<double> > results=invert.computeDistribution(xmin, xmax, [&](Complex &u, std::unordered_map<std::string, double> &params){return distToInvert(u, params, inits);}, params);
	auto end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);

	std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
	std::ofstream outputCSV;
	outputCSV.open("VaR9.csv");
	outputCSV <<"x, y, VaR"<<std::endl;
	for(int i=0; i<xNum; i++){
		outputCSV<<""<<results["x"][i]<<", "<<results["y"][i]<<", "<<results["VaR"][i]<<std::endl;
	}
	outputCSV.close();

}
