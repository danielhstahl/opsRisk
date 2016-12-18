#define _USE_MATH_DEFINES
#include <iostream>
#include <complex>
#include "FangOost.h"
#include "CharacteristicFunctions.h"
#include <fstream>
#include <cmath>
#include <ctime>
#include <chrono> //for accurate multithreading time using std::chrono

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
	auto density=fangoost::computeInv(xNum, uNum, xmin, xmax, [&](const auto& u){
		/*return chfunctions::jumpProcess(u, t, numODE, lambda, sigma, a, b, delta, v0, [&](const auto& uhat){
			return chfunctions::stableCF(uhat, alphaStable, muStable, betaStable, cStable);
		});*/
		return chfunctions::expAffine(
            rungekutta::computeFunctional(t, numODE, std::vector<std::complex<double> >({0, 0}),
                [&](double t, const std::vector<std::complex<double> >& x){
                    return chfunctions::duffieODE(
                        u+delta*x[0],//u
                        x, //current values
                        0.0, //rho0
                        0.0, //rho1
                        -a*b, //K0
                        -a, //K1,
                        0.0, //H0
                        sigma*sigma, //H1, 
                        0.0,//l0
                        lambda, //l1
                        [&](const auto& uhat){
							return chfunctions::stableCF(uhat, alphaStable, muStable, betaStable, cStable);
						}
                    );
                }
            ),
            v0
        );
	});
	auto axis=fangoost::computeXRange(xNum, xmin, xmax);
	std::ofstream out("output.csv");
	for(int i=0; i<xNum; ++i){
		out<<axis[i]<<", "<<density[i]<<std::endl;
	}
	out.close();
}
