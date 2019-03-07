CC = gcc
CFLAGS = -fopenmp -O3
#-vec-report -fno-alias
#-static-intel
LDFLAGS = -lm

headers = readConfig.h matrix.h lattice.h operator.h stochastic.h bicg-algo.h getCorr.h readNprintPhi.h
objects = main_dis.o readConfig.o matrix.o lattice.o operator.o stochastic.o bicg-algo.o getCorr.o readNprintPhi.o

run_dis.x : $(objects)
	$(CC) $(CFLAGS) -o $@ $(objects) $(LDFLAGS)

main_dis.o : main_dis.c $(headers)
	$(CC) $(CFLAGS) -c $<

readConfig.o : readConfig.c
	$(CC) $(CFLAGS) -c $<

matrix.o : matrix.c
	$(CC) $(CFLAGS) -c $<

lattice.o : lattice.c
	$(CC) $(CFLAGS) -c $<

operator.o : operator.c matrix.h
	$(CC) $(CFLAGS) -c $<

stochastic.o : stochastic.c matrix.h
	$(CC) $(CFLAGS) -c $<

bicg-algo.o : bicg-algo.c matrix.h operator.h
	$(CC) $(CFLAGS) -c $<

getCorr.o : getCorr.c readNprintPhi.h matrix.h
	$(CC) $(CFLAGS) -c $<

readNprintPhi.o : readNprintPhi.c matrix.h
	$(CC) $(CFLAGS) -c $<

clean :
	rm -f *.o *.x
