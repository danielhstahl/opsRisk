LDFLAGS=-L ../Complex -lComplex -L ../FangOosterlee -lFangOosterlee -L ../RungeKutta -lRungeKutta
INCLUDES=-I ../Complex -I ../RungeKutta -I ../FangOosterlee -I../rapidjson
opsRisk: main.o
	g++ -std=c++11 -O3  main.o  $(LDFLAGS) $(INCLUDES) -o opsRisk -fopenmp
main.o: main.cpp
	g++ -std=c++11 -O3  -c main.cpp $(LDFLAGS) $(INCLUDES) -fopenmp
clean:
	     -rm *.o marketRisk
