#include "FangOosterlee.h"

FangOosterlee::FangOosterlee(int k_, int h_) {
	k=k_; //u discretions
	h=h_; //x discretions
	//M_PI=3.14159265358979323846;
}
std::map<std::string, std::vector<double> > FangOosterlee::computeDistribution(cf characteristicFunction, std::map<std::string, double> parameters, double xmin, double xmax) {
	double xRange=xmax-xmin;
	double du=M_PI/xRange;
	double dx=xRange/(double)(h-1);
	double cp=2.0/xRange;
	std::vector<double> f=std::vector<double> (k);
	std::vector<double> y=std::vector<double> (h);
	std::vector<double> x=std::vector<double> (h);
	for(int j=0; j<k; j++){
		Complex u=Complex(0, du*j);
		f[j]=characteristicFunction(u, parameters).multiply(u.multiply(-xmin).exp()).getReal()*cp;
	}
	f[0]=.5*f[0];
	double exloss=0; //make these public later
	double vloss=0;
	for(int i=0;  i<h; i++){
		y[i]=0;
		x[i]=xmin+dx*i;
		for(int j=0; j<k; j++){
			y[i]=y[i]+f[j]*std::cos(du*j*dx*i);
		}
		vloss=vloss+y[i]*i*dx*dx*i;
		exloss=exloss+y[i]*i*dx;
	}
	std::map<std::string, std::vector<double> > distribution;
	distribution["x"]=x;
	distribution["y"]=y;
	return distribution;
}

