#define _USE_MATH_DEFINES
#include <iostream>
#include <complex>
#include "FangOost.h"
#include "CharacteristicFunctions.h"
#include <cmath>
#include <ctime>
#include <chrono> //for accurate multithreading time using std::chrono
#include "document.h" //rapidjson
#include "writer.h" //rapidjson
#include "stringbuffer.h" //rapidjson

template<typename Container, typename Range>
void printJson(const Container& myContainer, const Range& mn, const Range& dx){
	std::cout<<"{\"y\":[";
	auto n=myContainer.size();
	for(int i=0; i<n-1;++i){
		std::cout<<myContainer[i]<<",";
	}
	std::cout<<myContainer[n-1];
	std::cout<<"],\"xmin\":"<<mn<<",\"dx\":"<<dx<<"}";
}
auto computeXMax(double lambda, double muStable, double cStable, double t){
	const double largeScale=35;//
	return lambda*(muStable+largeScale*cStable)*t;
}
auto computeDelta(double rho, double muStable, double lambda){
	return rho/(muStable*lambda);
}
auto computeB(double rho){
	return 1.0-rho;
}
int main(int argc, char* argv[]){
	int xNum=1024;//standard parameters
	int uNum=256;//standard parameters
	double muStable=1300;
	double cStable=100;
	double alphaStable=1.1;
	const double betaStable=1; //needed to be "1" to have finite expectation
	double lambda=100;
	double rho=.9;
	double a=.4;
	double sigma=.4;
	double t=1;
	int numODE=128; 
	double v0=1;
	const double xmin=0;
	if(argc>1){
		rapidjson::Document parms;
		parms.Parse(argv[1]);
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
	}
	const double xmax=computeXMax(lambda, muStable, cStable, t);
	const double delta=computeDelta(rho, muStable, lambda); //they all have the same expected value
	const double b=computeB(rho);
	auto beta=chfunctions::AlphaOrBeta(0.0, -a, sigma*sigma, lambda);
    auto alpha=chfunctions::AlphaOrBeta(0.0, a*b, 0.0, 0.0);
	auto cf=chfunctions::augCF(alphaStable, muStable, betaStable, cStable);
	const auto density=fangoost::computeInv(xNum, uNum, xmin, xmax, [&](const auto& u){
		return chfunctions::expAffine(
            rungekutta::computeFunctional(t, numODE, std::vector<std::complex<double> >({0, 0}),
                [&](double t, const std::vector<std::complex<double> >& x){
					auto cfPart=cf(u+delta*x[0]);
					return std::vector<std::complex<double> >({
						beta(x[0], cfPart),
						alpha(x[0], cfPart)
					});
                }
            ),
            v0
        );
	});

	printJson(density, xmin, fangoost::computeDX(xNum, xmin, xmax));

}
