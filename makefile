CC=g++
OS := $(shell uname)



# check for Linux and add use of OpenMP
CFLAGS=  -Wall  -Ofast -std=c++11  -flto -pipe -funit-at-a-time  -Wfatal-errors -fopenmp -lz
LDFLAGS=$(zobu) -flto -lpthread -fopenmp -lz



EXEC=bgreat

all: $(EXEC)

aligner.o: aligner.cpp aligner.h utils.h alignerGreedy.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

utils.o: utils.cpp utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

bgreat: bgreat.o   aligner.o utils.o
	$(CC) -o $@ $^ $(LDFLAGS)

bgreat.o: bgreat.cpp  aligner.h
	$(CC) -o $@ -c $< $(CFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
