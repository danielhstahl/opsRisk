#include "CompoundPoisson.h"

CompoundPoisson::CompoundPoisson(double mu_, double sigma_, double lambda_, double t_) {
	mu=mu_; //drift
	sigma=sigma_;//volatility
	lambda=lambda_;//jump rate
	t=t_;//time
}
Complex CompoundPoisson::computeCF(Complex u, cf characteristicFunction, std::map<std::string, double> parameters) {
	return characteristicFunction(u, parameters).subtract(1.0).multiply(lambda).add(u.multiply(u).multiply(sigma*sigma*.5).add(u.multiply(mu)).multiply(t).exp();
}

