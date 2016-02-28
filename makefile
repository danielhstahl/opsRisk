LDFLAGS=-L ../Complex -lComplex -L ../FangOosterlee -lFangOosterlee -L ../RungeKutta -lRungeKutta
INCLUDES=-I ../Complex -I ../RungeKutta -I ../FangOosterlee -I../rapidjson -I ~/Documents/thirdParty/websocketpp/ -I ~/Documents/thirdParty/asio/
opsRisk: wsmain.o
	g++ -std=c++14 -O3  -pthread wsmain.o  $(LDFLAGS) $(INCLUDES) -o opsRisk -fopenmp
wsmain.o: wsmain.cpp
	g++ -std=c++14 -O3  -pthread -c wsmain.cpp $(LDFLAGS) $(INCLUDES) -fopenmp
clean:
	     -rm *.o opsRisk
