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
	int xNum=1024;//10x normal
	FangOosterlee invert(256, xNum); //4x normal



	//Stable params:
	double muStable=1300;
	double cStable=100;
	double alphaStable=1.1;
	double betaStable=1;

	//gamma params
	double alphaGamma=6.5;
	double betaGamma=200;

	//inverse gaussian params
	double muIG=1300;
	double lambdaIG=5000;//better be large!

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
	auto start = std::chrono::system_clock::now();
	std::unordered_map<std::string, std::vector<double> > resultsStable=invert.computeDistribution(xmin, xmax,
		[&](Complex &u, std::unordered_map<std::string, double> &parms){
			inits[0]=Complex(0, 0);
			inits[1]=Complex(0, 0);
			return distToInvert(u, parms, inits,
				[&](Complex &uhat){
					return stableCF(uhat, alphaStable, muStable, betaStable, cStable);
				});
			}, params);
	auto end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);
	std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
	std::ofstream outputCSV;
	outputCSV.open("VaR9Stable.csv");
	outputCSV <<"x, y, VaR"<<std::endl;
	for(int i=0; i<xNum; i++){
		outputCSV<<""<<resultsStable["x"][i]<<", "<<resultsStable["y"][i]<<", "<<resultsStable["VaR"][i]<<std::endl;
	}
	outputCSV.close();

	start = std::chrono::system_clock::now();
	std::unordered_map<std::string, std::vector<double> > resultsGamma=invert.computeDistribution(xmin, xmax,
		[&](Complex &u, std::unordered_map<std::string, double> &parms){
			inits[0]=Complex(0, 0);
			inits[1]=Complex(0, 0);
			return distToInvert(u, parms, inits,
				[&](Complex &uhat){
					return gammaCF(uhat, alphaGamma, betaGamma);
				});
			}, params);
	end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);

	std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
	//std::ofstream outputCSV;
	outputCSV.open("VaR9Gamma.csv");
	outputCSV <<"x, y, VaR"<<std::endl;
	for(int i=0; i<xNum; i++){
		outputCSV<<""<<resultsGamma["x"][i]<<", "<<resultsGamma["y"][i]<<", "<<resultsGamma["VaR"][i]<<std::endl;
	}
	outputCSV.close();


	start = std::chrono::system_clock::now();
	std::unordered_map<std::string, std::vector<double> > resultsIG=invert.computeDistribution(xmin, xmax,
		[&](Complex &u, std::unordered_map<std::string, double> &parms){
			inits[0]=Complex(0, 0);
			inits[1]=Complex(0, 0);
			return distToInvert(u, parms, inits,
				[&](Complex &uhat){
					return inverseGaussianCF(uhat, muIG, lambdaIG);
				});
			}, params);
	end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);

	std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
//	std::ofstream outputCSV;
	outputCSV.open("VaR9IG.csv");
	outputCSV <<"x, y, VaR"<<std::endl;
	for(int i=0; i<xNum; i++){
		outputCSV<<""<<resultsIG["x"][i]<<", "<<resultsIG["y"][i]<<", "<<resultsIG["VaR"][i]<<std::endl;
	}
	outputCSV.close();


	/*START OF .5*/
	std::unordered_map<std::string, double> params5; //im not a huge fan of this but it works

	params5["lambda"]=100;
	rho=.5;
	params5["delta"]=rho/(muStable*params["lambda"]); //they all have the same expected value
	params5["b"]=1-rho;
	params5["a"]=.4;
	params5["sigma"]=.4;

	params5["t"]=1;
	params5["numODE"]=params["numODE"];
	params5["v0"]=1;

	params5["delta"]=rho/(muStable*params5["lambda"]); //they all have the same expected value
	params5["b"]=1-rho;
	start = std::chrono::system_clock::now();
	std::unordered_map<std::string, std::vector<double> > results5Stable=invert.computeDistribution(xmin, xmax,
		[&](Complex &u, std::unordered_map<std::string, double> &parms){
			inits[0]=Complex(0, 0);
			inits[1]=Complex(0, 0);
			return distToInvert(u, parms, inits,
				[&](Complex &uhat){
					return stableCF(uhat, alphaStable, muStable, betaStable, cStable);
				});
			}, params5);
	end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);
	std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
	outputCSV.open("VaR5Stable.csv");
	outputCSV <<"x, y, VaR"<<std::endl;
	for(int i=0; i<xNum; i++){
		outputCSV<<""<<results5Stable["x"][i]<<", "<<results5Stable["y"][i]<<", "<<results5Stable["VaR"][i]<<std::endl;
	}
	outputCSV.close();


	start = std::chrono::system_clock::now();
	std::unordered_map<std::string, std::vector<double> > results5Gamma=invert.computeDistribution(xmin, xmax,
		[&](Complex &u, std::unordered_map<std::string, double> &parms){
			inits[0]=Complex(0, 0);
			inits[1]=Complex(0, 0);
			return distToInvert(u, parms, inits,
				[&](Complex &uhat){
					return gammaCF(uhat, alphaGamma, betaGamma);
				});
			}, params5);
	end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);

	std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
	outputCSV.open("VaR5Gamma.csv");
	outputCSV <<"x, y, VaR"<<std::endl;
	for(int i=0; i<xNum; i++){
		outputCSV<<""<<results5Gamma["x"][i]<<", "<<results5Gamma["y"][i]<<", "<<results5Gamma["VaR"][i]<<std::endl;
	}
	outputCSV.close();
	start = std::chrono::system_clock::now();
	std::unordered_map<std::string, std::vector<double> > results5IG=invert.computeDistribution(xmin, xmax,
		[&](Complex &u, std::unordered_map<std::string, double> &parms){
			inits[0]=Complex(0, 0);
			inits[1]=Complex(0, 0);
			return distToInvert(u, parms, inits,
				[&](Complex &uhat){
					return inverseGaussianCF(uhat, muIG, lambdaIG);
				});
			}, params5);
	end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);

	std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
	//	std::ofstream outputCSV;
	outputCSV.open("VaR5IG.csv");
	outputCSV <<"x, y, VaR"<<std::endl;
	for(int i=0; i<xNum; i++){
		outputCSV<<""<<results5IG["x"][i]<<", "<<results5IG["y"][i]<<", "<<results5IG["VaR"][i]<<std::endl;
	}
	outputCSV.close();





	/*START OF 0*/
	std::unordered_map<std::string, double> params0; //im not a huge fan of this but it works

	params0["lambda"]=100;
	rho=0;
	params0["delta"]=rho/(muStable*params["lambda"]); //they all have the same expected value
	params0["b"]=1-rho;
	params0["a"]=.4;
	params0["sigma"]=.4;

	params0["t"]=1;
	params0["numODE"]=params["numODE"];
	params0["v0"]=1;

	params0["delta"]=rho/(muStable*params0["lambda"]); //they all have the same expected value
	params0["b"]=1-rho;
	start = std::chrono::system_clock::now();
	std::unordered_map<std::string, std::vector<double> > results0Stable=invert.computeDistribution(xmin, xmax,
		[&](Complex &u, std::unordered_map<std::string, double> &parms){
			inits[0]=Complex(0, 0);
			inits[1]=Complex(0, 0);
			return distToInvert(u, parms, inits,
				[&](Complex &uhat){
					return stableCF(uhat, alphaStable, muStable, betaStable, cStable);
				});
			}, params0);
	end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);

	std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
	outputCSV.open("VaR0Stable.csv");
	outputCSV <<"x, y, VaR"<<std::endl;
	for(int i=0; i<xNum; i++){
		outputCSV<<""<<results0Stable["x"][i]<<", "<<results0Stable["y"][i]<<", "<<results0Stable["VaR"][i]<<std::endl;
	}
	outputCSV.close();


	start = std::chrono::system_clock::now();
	std::unordered_map<std::string, std::vector<double> > results0Gamma=invert.computeDistribution(xmin, xmax,
		[&](Complex &u, std::unordered_map<std::string, double> &parms){
			inits[0]=Complex(0, 0);
			inits[1]=Complex(0, 0);
			return distToInvert(u, parms, inits,
				[&](Complex &uhat){
					return gammaCF(uhat, alphaGamma, betaGamma);
				});
			}, params0);
	end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);

	std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
	//std::ofstream outputCSV;
	outputCSV.open("VaR0Gamma.csv");
	outputCSV <<"x, y, VaR"<<std::endl;
	for(int i=0; i<xNum; i++){
		outputCSV<<""<<results0Gamma["x"][i]<<", "<<results0Gamma["y"][i]<<", "<<results0Gamma["VaR"][i]<<std::endl;
	}
	outputCSV.close();


	start = std::chrono::system_clock::now();
	std::unordered_map<std::string, std::vector<double> > results0IG=invert.computeDistribution(xmin, xmax,
		[&](Complex &u, std::unordered_map<std::string, double> &parms){
			inits[0]=Complex(0, 0);
			inits[1]=Complex(0, 0);
			return distToInvert(u, parms, inits,
				[&](Complex &uhat){
					return inverseGaussianCF(uhat, muIG, lambdaIG);
				});
			}, params0);
	end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);

	std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
	//	std::ofstream outputCSV;
	outputCSV.open("VaR0IG.csv");
	outputCSV <<"x, y, VaR"<<std::endl;
	for(int i=0; i<xNum; i++){
		outputCSV<<""<<results0IG["x"][i]<<", "<<results0IG["y"][i]<<", "<<results0IG["VaR"][i]<<std::endl;
	}
	outputCSV.close();




	std::unordered_map<std::string, double> paramsLDA; //im not a huge fan of this but it works

	paramsLDA["lambda"]=100;
	rho=0;
	paramsLDA["delta"]=rho/(muStable*params["lambda"]); //they all have the same expected value
	paramsLDA["b"]=1-rho;
	paramsLDA["a"]=0;
	paramsLDA["sigma"]=0;

	paramsLDA["t"]=1;
	paramsLDA["numODE"]=params["numODE"];
	paramsLDA["v0"]=1;

	paramsLDA["delta"]=rho/(muStable*paramsLDA["lambda"]); //they all have the same expected value
	paramsLDA["b"]=1-rho;
	start = std::chrono::system_clock::now();
	std::unordered_map<std::string, std::vector<double> > resultsLDAStable=invert.computeDistribution(xmin, xmax,
		[&](Complex &u, std::unordered_map<std::string, double> &parms){
			inits[0]=Complex(0, 0);
			inits[1]=Complex(0, 0);
			return distToInvert(u, parms, inits,
				[&](Complex &uhat){
					return stableCF(uhat, alphaStable, muStable, betaStable, cStable);
				});
			}, paramsLDA);
	end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);
	std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
	outputCSV.open("VaRLDAStable.csv");
	outputCSV <<"x, y, VaR"<<std::endl;
	for(int i=0; i<xNum; i++){
		outputCSV<<""<<resultsLDAStable["x"][i]<<", "<<resultsLDAStable["y"][i]<<", "<<resultsLDAStable["VaR"][i]<<std::endl;
	}
	outputCSV.close();


	start = std::chrono::system_clock::now();
	std::unordered_map<std::string, std::vector<double> > resultsLDAGamma=invert.computeDistribution(xmin, xmax,
		[&](Complex &u, std::unordered_map<std::string, double> &parms){
			inits[0]=Complex(0, 0);
			inits[1]=Complex(0, 0);
			return distToInvert(u, parms, inits,
				[&](Complex &uhat){
					return gammaCF(uhat, alphaGamma, betaGamma);
				});
			}, paramsLDA);
	end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);

	std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
	outputCSV.open("VaRLDAGamma.csv");
	outputCSV <<"x, y, VaR"<<std::endl;
	for(int i=0; i<xNum; i++){
		outputCSV<<""<<resultsLDAGamma["x"][i]<<", "<<resultsLDAGamma["y"][i]<<", "<<resultsLDAGamma["VaR"][i]<<std::endl;
	}
	outputCSV.close();
	start = std::chrono::system_clock::now();
	std::unordered_map<std::string, std::vector<double> > resultsLDAIG=invert.computeDistribution(xmin, xmax,
		[&](Complex &u, std::unordered_map<std::string, double> &parms){
			inits[0]=Complex(0, 0);
			inits[1]=Complex(0, 0);
			return distToInvert(u, parms, inits,
				[&](Complex &uhat){
					return inverseGaussianCF(uhat, muIG, lambdaIG);
				});
			}, paramsLDA);
	end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);

	std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
	//	std::ofstream outputCSV;
	outputCSV.open("VaRLDAIG.csv");
	outputCSV <<"x, y, VaR"<<std::endl;
	for(int i=0; i<xNum; i++){
		outputCSV<<""<<resultsLDAIG["x"][i]<<", "<<resultsLDAIG["y"][i]<<", "<<resultsLDAIG["VaR"][i]<<std::endl;
	}
	outputCSV.close();




}
