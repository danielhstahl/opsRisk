#LDFLAGS=-L ../Complex -lComplex
INCLUDES=-I ../RungeKutta -I ../FangOost -I../rapidjson/include/rapidjson/ -I../websocketpp/ -I../asio/asio/include/ -I ../FunctionalUtilities
opsRisk: main.o
	g++ -std=c++14 -O3  -pthread main.o  $(LDFLAGS) $(INCLUDES) -o opsRisk -fopenmp
main.o: main.cpp
	g++ -std=c++14 -O3  -pthread -c main.cpp $(LDFLAGS) $(INCLUDES) -fopenmp
clean:
	     -rm *.o opsRisk
