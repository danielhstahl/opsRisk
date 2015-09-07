#include "FangOosterlee.h"

FangOosterlee::FangOosterlee(int k_, int h_) {
	k=k_; //u discretions
	h=h_; //x discretions
	//M_PI=3.14159265358979323846;
}
//void extendedDist(cf chrF, Complex u, std::map<std::string, double> &parameters, double &f, double xmin, double cp){
	//f=chrF(u, parameters).multiply(u.multiply(-xmin).exp()).getReal()*cp;
//}


/*template< typename FN, typename... ARGS>
std::map<std::string, std::vector<double> > FangOosterlee::computeDistribution(double xmin, double xmax, FN&& fn, ARGS&&... args ) {
	double xRange=xmax-xmin;
	double du=M_PI/xRange;
	double dx=xRange/(double)(h-1);
	double cp=2.0/xRange;
	std::vector<double> f=std::vector<double> (k);
	std::vector<double> y=std::vector<double> (h);

	std::vector<double> x=std::vector<double> (h);
	std::vector<std::thread> thrd=std::vector<std::thread> (k);
	//parameters["cp"]=cp;
	//parameters["xmin"]=xmin;
	for(int j=0; j<k; j++){
		Complex u=Complex(0, du*j);
		//thrd[j]=std::thread(extendedDist, characteristicFunction, u,  std::ref(parameters), std::ref(f[j]),  xmin, cp);
		//f[j]=std::forward<FN>(fn), std::forward<ARGS>(u, args...).multiply(u.multiply(-xmin).exp()).getReal()*cp;
		f[j]=fn(u, args...).multiply(u.multiply(-xmin).exp()).getReal()*cp;
	}


	//for(int j=0; j<k; j++){
			//thrd[j].join();
	//}
	//for(auto& t : thrd){
		//t.join();
	//}


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
}*/
