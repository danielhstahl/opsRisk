#define _USE_MATH_DEFINES
#define _WEBSOCKETPP_CPP11_THREAD_
#define _WEBSOCKETPP_CPP11_CHRONO_
#define _WEBSOCKETPP_CPP11_TYPE_TRAITS_
#define ASIO_STANDALONE
#define ASIO_HAS_STD_ARRAY
#define ASIO_HAS_STD_ATOMIC
#define ASIO_HAS_CSTDINT
#define ASIO_HAS_STD_ADDRESSOF
#define ASIO_HAS_STD_SHARED_PTR
#define ASIO_HAS_STD_TYPE_TRAITS
#include <iostream>
#include <set>
#include <asio.hpp>
#include <websocketpp/config/asio_no_tls.hpp>
//#include <websocketpp/config/asio.hpp>
#include <websocketpp/config/core.hpp>
#include <websocketpp/server.hpp>
#include <thread>
#include <unordered_map>
#include "FangOosterlee.h"
#include "RungeKutta.h"
#include "Complex.h"
#include <cmath>
#include <ctime>
#include <chrono> //for accurate multithreading time using std::chrono
#include "document.h" //rapidjson
#include "writer.h" //rapidjson
#include "stringbuffer.h" //rapidjson

typedef websocketpp::server<websocketpp::config::asio> server;
using websocketpp::connection_hdl;
using websocketpp::lib::placeholders::_1;
using websocketpp::lib::placeholders::_2;
using websocketpp::lib::bind;


Complex GaussCF(Complex u, std::unordered_map<std::string, double>& params){
	return u.multiply(params["mu"]).add(u.multiply(u).multiply(params["sigma"]*params["sigma"]*.5)).exp();
}

Complex stableCF(Complex &u, double alpha, double mu, double beta, double c){
	double phi=tan(alpha*.5*M_PI);
	Complex uC=u.multiply(Complex(0, -1));
	return u.multiply(mu).subtract(uC.multiply(c).pow(alpha).multiply(Complex(1, -beta*phi))).exp();
}
Complex gammaCF(Complex &u, double a, double b){
	return pow(1-u*b, -a);
}
Complex inverseGaussianCF(Complex &u, double mu, double lambda){
	return exp((lambda/mu)*(1-sqrt(1-(2*mu*mu*u)/lambda)));
}
template<typename CF>
std::vector<Complex> DuffieODE(Complex &u, CF &cf, std::vector<Complex> &initialValues, double sigma, double lambda, double a, double delta, double k){ //double alpha, double mu, double beta, double c,
	std::vector<Complex> vls(2);
	double sig=sigma*sigma*.5;
	Complex uBeta=u.add(initialValues[0].multiply(delta));

	vls[0]=initialValues[0].multiply(initialValues[0]).multiply(sig).add(cf(uBeta).multiply(lambda)).subtract(lambda).subtract(initialValues[0].multiply(a*(1+k)));
	vls[1]=initialValues[0].multiply(a);
	return vls;
}
template<typename CF>
Complex distToInvert(Complex &u, std::unordered_map<std::string, double> &params, std::vector<Complex> &inits, const CF &cf){
	RungeKutta rg(params["t"], (int)params["numODE"]);
	std::vector<Complex> vls=rg.compute([&](double t, std::vector<Complex> &x){
		return DuffieODE(u, cf, x, params["sigma"], params["lambda"], params["a"], params["delta"], params["b"]);
	}, inits);
	return vls[0].multiply(params["v0"]).add(vls[1]).exp();
}

class WS{
private:

	typedef std::set<connection_hdl,std::owner_less<connection_hdl>> con_list;
	server m_server;
	con_list m_connections;
public:
	WS(){

		m_server.init_asio();
		m_server.set_open_handler(bind(&WS::on_open,this,::_1));
		m_server.set_close_handler(bind(&WS::on_close,this,::_1));
		m_server.set_message_handler(bind(&WS::on_message,this,::_1,::_2));
	}
	void myFunction(connection_hdl connection, server::message_ptr msg){
		int xNum=1024;//standard parameters
		int uNum=256;//standard parameters
		double muStable=1300;
		double cStable=100;
		double alphaStable=1.1;
		double betaStable=1;
		std::unordered_map<std::string, double> params; //im not a huge fan of this but it works
		double rho;
		double xmin;
		double xmax;
		params["lambda"]=100;
		rho=5;
		params["delta"]=rho/(muStable*params["lambda"]); //they all have the same expected value
		params["a"]=.4;
		params["k"]=1+params["delta"]*params["lambda"]*muStable/params["a"];
		params["sigma"]=.4;
		params["t"]=1;
		params["numODE"]=128;
		params["v0"]=1;
		xmin=0;
		xmax=params["lambda"]*(muStable+35*cStable);

		std::vector<Complex> inits(2);
		rapidjson::Document parms;
		parms.Parse(msg->get_payload().c_str());
		if(parms.FindMember("xNum")!=parms.MemberEnd()){
			xNum=parms["xNum"].GetInt();
		}
		if(parms.FindMember("uNum")!=parms.MemberEnd()){
			uNum=parms["uNum"].GetInt();
		}
		if(parms.FindMember("numODE")!=parms.MemberEnd()){
			params["numODE"]=parms["numODE"].GetInt();
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
			params["lambda"]=parms["lambda"].GetDouble();
		}
		if(parms.FindMember("rho")!=parms.MemberEnd()){
			rho=parms["rho"].GetDouble();
		}
		if(parms.FindMember("v0")!=parms.MemberEnd()){
			params["v0"]=parms["v0"].GetDouble();
		}
		if(parms.FindMember("t")!=parms.MemberEnd()){
			params["t"]=parms["t"].GetDouble();
		}
		if(parms.FindMember("a")!=parms.MemberEnd()){
			params["a"]=parms["a"].GetDouble();
		}
		if(parms.FindMember("sigma")!=parms.MemberEnd()){
			params["sigma"]=parms["sigma"].GetDouble();
		}
		params["delta"]=rho/(muStable*params["lambda"]); //they all have the same expected value
		//params["a"]=.4;
		if(params["a"]==0){
			params["k"]=1;
		}
		else{
			params["k"]=1+params["delta"]*params["lambda"]*muStable/params["a"];
		}

		xmax=params["lambda"]*(muStable+35*cStable);
		FangOosterlee invert(uNum, xNum);
		//auto WSsend=;
		//auto start = std::chrono::system_clock::now();
		invert.computeDistributionWS(
			[&](const std::string& message){
				m_server.send(connection,message, websocketpp::frame::opcode::text);
			},
			false,
			xmin, xmax,
			[&](Complex &u, std::unordered_map<std::string, double> &parms){
				inits[0]=Complex(0, 0);
				inits[1]=Complex(0, 0);
				return distToInvert(u, parms, inits,
					[&](Complex &uhat){
						return stableCF(uhat, alphaStable, muStable, betaStable, cStable);
					});
				}, params);
	}
	void on_open(connection_hdl hdl) {
			m_connections.insert(hdl);
	}
	void on_close(connection_hdl hdl) {
			m_connections.erase(hdl);
	}
	void on_message(connection_hdl hdl, server::message_ptr msg) {
		std::thread myThread(&WS::myFunction, this, hdl, msg);
		myThread.detach();
			/*for (auto it : m_connections) {
					m_server.send(it,msg);
			}*/
	}
	void run(uint16_t port) {
			m_server.listen(port);
			m_server.start_accept();
			m_server.run();
	}
};
int main() {
    WS server;//([&](auto& server, auto& connection));
    server.run(9002);
}
