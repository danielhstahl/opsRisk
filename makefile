LDFLAGS=-L ../Complex -lComplex -L ../FangOosterlee -lFangOosterlee -L ../RungeKutta -lRungeKutta
INCLUDES=-I ../Complex -I ../RungeKutta -I ../FangOosterlee
opsRisk:
	g++ -std=c++11 -O3  main.cpp  $(LDFLAGS) $(INCLUDES) -o opsRisk -fopenmp
clean:
	     \rm *.o *~ p1
