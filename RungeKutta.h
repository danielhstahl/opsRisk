#ifndef __RUNGEKUTTA_H_INCLUDED__
#define __RUNGEKUTTA_H_INCLUDED__

#include <cmath>
#include "Complex.h"
#include <vector>
#include <map>
#include <iostream> //debugging

class RungeKutta { 
	typedef std::vector<Complex> (*cf)(double, std::vector<Complex>, std::map<std::string, double>, std::map<std::string, Complex>); //defines cf as a pointer to a function which takes complex and outputs complex as arguments...as of now all arguments must be doubles...
	private:
		int numSteps;
		std::vector<Complex> initialValues;
		double t;
		//std::map<std::string, double> params;
		//std::map<std::string, Complex> cmplParams;
		cf funcToSolve;
		
		//void init(cf, double, std::vector<Complex>, int, std::map<std::string, double>, std::map<std::string, Complex>);
	public:
		RungeKutta(cf, double, std::vector<Complex>, int);
		/*RungeKutta(cf, double, std::vector<Complex>, int, std::map<std::string, double>);
		RungeKutta(cf, double, std::vector<Complex>, int, std::map<std::string, Complex>);
		RungeKutta(cf, double, std::vector<Complex>, int, std::map<std::string, double>, std::map<std::string, Complex>);
		*/
		std::vector<Complex> compute(std::map<std::string, double>, std::map<std::string, Complex>);
};
#endif