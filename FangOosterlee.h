#ifndef __FANGOOSTERLEE_H_INCLUDED__
#define __FANGOOSTERLEE_H_INCLUDED__
#define _USE_MATH_DEFINES
#include <cmath>
#include "Complex.h"
#include <vector>
#include <map>
#include <iostream> //debugging

class FangOosterlee { 
	typedef Complex (*cf)(Complex, std::map<std::string, double>); //defines cf as a pointer to a function which takes complex and outputs complex as arguments...as of now all arguments must be doubles...
	private:
		int k;
		int h;
		//double M_PI;
	public:
		FangOosterlee(int, int);
		//std::map<std::string, std::vector<double> > computeDistribution(ICharacteristicFunction*, double, double);
		std::map<std::string, std::vector<double> > computeDistribution(cf, std::map<std::string, double>, double, double);

};
#endif