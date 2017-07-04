# define the C++ compiler to use
CC = g++

# define any compile-time flags
CFLAGS = -fopenmp -std=c++0x -O3 -W -Wall # for OpenMP

# define any libraries to link into executable
LIBS = -lz

# the build target executable:
all: SAGE2

SAGE2: main.o utils.o fastAQReader.o inputReader.o readLoader.o hashTable.o economyGraph.o overlapGraph.o simplification.o copycountEstimator.o cs2.o matePair.o mergeContigs.o scaffolding.o
	$(CC) $(CFLAGS) -o SAGE2 main.o utils.o fastAQReader.o inputReader.o readLoader.o hashTable.o economyGraph.o overlapGraph.o simplification.o copycountEstimator.o cs2.o matePair.o  mergeContigs.o scaffolding.o $(LIBS)

main.o: ./main.cpp
	$(CC) $(CFLAGS) -o main.o -c ./main.cpp

utils.o: ./utils.cpp
	$(CC) $(CFLAGS) -o utils.o -c ./utils.cpp

fastAQReader.o: ./inputReader/fastAQReader.cpp
	$(CC) $(CFLAGS) -o fastAQReader.o -c ./inputReader/fastAQReader.cpp

inputReader.o: ./inputReader/inputReader.cpp
	$(CC) $(CFLAGS) -o inputReader.o -c ./inputReader/inputReader.cpp

readLoader.o: ./inputReader/readLoader.cpp
	$(CC) $(CFLAGS) -o readLoader.o -c ./inputReader/readLoader.cpp

hashTable.o: ./economyGraph/hashTable.cpp
	$(CC) $(CFLAGS) -o hashTable.o -c ./economyGraph/hashTable.cpp

economyGraph.o: ./economyGraph/economyGraph.cpp
	$(CC) $(CFLAGS) -o economyGraph.o -c ./economyGraph/economyGraph.cpp

overlapGraph.o: ./overlapGraph/overlapGraph.cpp
	$(CC) $(CFLAGS) -o overlapGraph.o -c ./overlapGraph/overlapGraph.cpp
	
simplification.o: ./overlapGraph/simplification.cpp
	$(CC) $(CFLAGS) -o simplification.o -c ./overlapGraph/simplification.cpp
	
copycountEstimator.o: ./overlapGraph/copycountEstimator.cpp
	$(CC) $(CFLAGS) -o copycountEstimator.o -c ./overlapGraph/copycountEstimator.cpp
	
cs2.o: ./overlapGraph/cs2/cs2.cpp
	$(CC) $(CFLAGS) -o cs2.o -c ./overlapGraph/cs2/cs2.cpp

matePair.o: ./matePair/matePair.cpp
	$(CC) $(CFLAGS) -o matePair.o -c ./matePair/matePair.cpp

mergeContigs.o: ./matePair/mergeContigs.cpp
	$(CC) $(CFLAGS) -o mergeContigs.o -c ./matePair/mergeContigs.cpp

scaffolding.o: ./matePair/scaffolding.cpp
	$(CC) $(CFLAGS) -o scaffolding.o -c ./matePair/scaffolding.cpp
	
clean:
	$(RM) SAGE2 *.o
