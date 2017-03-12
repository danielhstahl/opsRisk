#LDFLAGS=-L ../Complex -lComplex
INCLUDES=-I ../RungeKutta -I ../FangOost -I ../FunctionalUtilities -I ../CharacteristicFunctions -I ../rapidjson/include/rapidjson
opsRisk: main.o
	g++ -std=c++14 -O3 $(STATIC)  -pthread main.o  $(LDFLAGS) $(INCLUDES) -o opsRisk -fopenmp
main.o: main.cpp
	g++ -std=c++14 -O3 $(STATIC) -pthread -c main.cpp $(LDFLAGS) $(INCLUDES) -fopenmp
clean:
	     -rm *.o opsRisk
