#LDFLAGS=-L ../Complex -lComplex
INCLUDES=-I ../RungeKutta -I ../FangOost -I ../FunctionalUtilities -I ../CharacteristicFunctions
opsRisk: main.o
	g++ -std=c++14 -O3  -pthread main.o  $(LDFLAGS) $(INCLUDES) -o opsRisk -fopenmp
main.o: main.cpp
	g++ -std=c++14 -O3  -pthread -c main.cpp $(LDFLAGS) $(INCLUDES) -fopenmp
clean:
	     -rm *.o opsRisk
