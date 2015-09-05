#include "RungeKutta.h"
RungeKutta::RungeKutta(cf funcToSolve_, double t_, std::vector<Complex> initialValues_, int numSteps_){
	funcToSolve=funcToSolve_;
	t=t_;
	initialValues=initialValues_;
	numSteps=numSteps_;

}
std::vector<Complex> RungeKutta::compute(std::map<std::string, double> &params, std::map<std::string, Complex> &cmplParams){
	double h=t/numSteps;
	double hlfh=h*.5;
	double sixthh=h/6.0;
	int n=initialValues.size();
	for(int i=0; i<numSteps; i++){
		std::vector<Complex> newVal(n);
		std::vector<Complex> k1=funcToSolve(h*i, initialValues, params, cmplParams);
		for(int j=0; j<n;j++){
			newVal[j]=initialValues[j].add(k1[j].multiply(hlfh));
		}
		std::vector<Complex> k2=funcToSolve(h*i+hlfh, newVal, params, cmplParams);
		for(int j=0; j<n;j++){
			newVal[j]=initialValues[j].add(k2[j].multiply(hlfh));
		}
		std::vector<Complex> k3=funcToSolve(h*i+hlfh, newVal, params, cmplParams);
		for(int j=0; j<n;j++){
			newVal[j]=initialValues[j].add(k3[j].multiply(h));
		}
		std::vector<Complex> k4=funcToSolve(h*i+h, newVal, params, cmplParams);
		for(int j=0; j<n;j++){
			initialValues[j]=initialValues[j].add(k1[j].add(k2[j].multiply(2)).add(k3[j].multiply(2)).add(k4[j]).multiply(sixthh));
		}
	}
	return initialValues;
}
