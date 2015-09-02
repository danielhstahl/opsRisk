#ifndef __COMPOUNDPOISSON_H_INCLUDED__
#define __COMPOUNDPOISSON_H_INCLUDED__

#include <cmath>
#include "Complex.h"
#include <vector>
#include <map>
#include <iostream> //debugging

class CompoundPoisson { 
	typedef Complex (*cf)(Complex, std::map<std::string, double>); //defines cf as a pointer to a function which takes complex and outputs complex as arguments...as of now all arguments must be doubles...
	private:
		double lambda;
		double t;
		double mu;
		double sigma;
	public:
		CompoundPoisson(double, double, double, double);
		//std::map<std::string, std::vector<double> > computeDistribution(ICharacteristicFunction*, double, double);
		Complex computeCF(cf, std::map<std::string, double>);
};
#endif