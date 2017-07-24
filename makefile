SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

CC  = gcc-4.7 -O0 -std=gnu99
COPT  = -m64 -fPIC 

CLNFLAGS = -lconcert -lilocplex -lcplex -lpthread -lm

INCLUDE = -I/opt/ibm/ILOG/CPLEX_Studio1263/cplex/include -I/opt/ibm/ILOG/CPLEX_Studio1263/concert/include/ilconcert

LIBS = -L/opt/ibm/ILOG/CPLEX_Studio1263/cplex/lib/x86-64_linux/static_pic -L/opt/ibm/ILOG/CPLEX_Studio1263/concert/lib/x86-64_linux/static_pic

all: dwdecomp  

check.o: check.c
	$(CC) $(COPT) $(INCLUDE) $(LIBS) -std=gnu99 check.c -g -c $(CLNFLAGS)

build.o: build.c
	$(CC) $(COPT) $(INCLUDE) $(LIBS) -std=gnu99 build.c -g -c $(CLNFLAGS)

util.o: util.c
	$(CC) $(COPT) $(INCLUDE) $(LIBS) -std=gnu99 util.c -g -c $(CLNFLAGS)

branch.o: branch.c
	$(CC) $(COPT) $(INCLUDE) $(LIBS) -std=gnu99 branch.c -g -c $(CLNFLAGS)

solve.o: solve.c
	$(CC) $(COPT) $(INCLUDE) $(LIBS) -std=gnu99 solve.c -g -c $(CLNFLAGS)

debug.o: debug.c
	$(CC) $(COPT) $(INCLUDE) $(LIBS) -std=gnu99 debug.c -g -c $(CLNFLAGS)

dwdecomp: dwdecomp.c check.o dwdecomp.h build.o branch.o solve.o util.o
	$(CC) $(COPT) $(INCLUDE) $(LIBS) -std=gnu99 dwdecomp.c check.o build.o branch.o util.o solve.o -g -o dwdecomp  $(CLNFLAGS)


clean:
	rm -f *.lp *.sol *.o dwdecomp

