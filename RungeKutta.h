#ifndef __RUNGEKUTTA_H_INCLUDED__
#define __RUNGEKUTTA_H_INCLUDED__

#include <cmath>
//#include "Complex.h"
#include <vector>
#include <unordered_map>
#include <iostream> //debugging


class RungeKutta {
	//typedef std::vector<Complex> (*cf)(double, std::vector<Complex>, std::map<std::string, double>, std::map<std::string, Complex>); //defines cf as a pointer to a function which takes complex and outputs complex as arguments...as of now all arguments must be doubles...

	private:
		int numSteps;
		//std::vector<T> initialValues;
		double t;

		//cf funcToSolve;

	public:
		RungeKutta(double,  int);
		template<typename T, typename FN, typename... ARGS>
		std::vector<T> compute(FN&& fn, std::vector<T> initialValues, ARGS&&... args){
			double h=t/numSteps;
			double hlfh=h*.5;
			double sixthh=h/6.0;
			int n=initialValues.size();
			for(int i=0; i<numSteps; i++){
				std::vector<T> newVal(n);
				std::vector<T> k1=fn(h*i, initialValues, args...);
				for(int j=0; j<n;j++){
					newVal[j]=initialValues[j]+k1[j]*hlfh;
					//newVal[j]=initialValues[j].add(k1[j].multiply(hlfh));
				}
				std::vector<T> k2=fn(h*i+hlfh, newVal, args...);
				for(int j=0; j<n;j++){
					newVal[j]=initialValues[j]+k2[j]*hlfh;
				}
				std::vector<T> k3=fn(h*i+hlfh, newVal, args...);
				for(int j=0; j<n;j++){
					newVal[j]=initialValues[j]+k3[j]*h;
				}
				std::vector<T> k4=fn(h*i+h, newVal, args...);
				for(int j=0; j<n;j++){
					initialValues[j]=initialValues[j]+(k1[j]+k2[j]*2+k3[j]*2+k4[j])*sixthh;
				}
			}
			return initialValues;
		}
};
#endif
