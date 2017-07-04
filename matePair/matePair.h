/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 ****************************************************************************/

#ifndef MATEPAIR_H_
#define MATEPAIR_H_

#include <deque>
#include "../utils.h"
#include "../inputReader/readLoader.h"
#include "../overlapGraph/overlapGraph.h"

class MatePairInfo // Stores mate pair information
{
public:
	uint64_t ID;
	uint8_t freq;
	bool flag;//if both mates on the same edge 0, otherwise 1
	bool type1;
	bool type2;
	int8_t library;
	MatePairInfo *next;
};

class ReadToEdgeMap //this structure stores information for mapped reads
{
public:
	uint32_t *locationForward;
	uint32_t *locationReverse;
	Edge *edge;
	ReadToEdgeMap *next;
};

class MatePair
{
public:
	MatePairInfo **matePairList;
	ReadToEdgeMap **readToEdgeList;
	uint *Mean, *standardDeviation;
	int *upperBoundOfInsert, *lowerBoundOfInsert;
	uint minimumUpperBoundOfInsert, maximumUpperBoundOfInsert;
private:
	ReadLoader *loaderObj;
	OverlapGraph *graphObj;
	int numberOfLibrary;
public:
	MatePair(OverlapGraph *graph1, ReadLoader *loader1);
	~MatePair();
	void mapMatePairsFromList(string listPath);
	void mapMatePairs(string mateFile1, string mateFile2, int library);
	void processMatePairs(deque<string> &readsArray, int library);
	void meanSdEstimation();
	void computeMeanSD(int mu, int sd, int library, int *returnMu, int *returnSD);
	void mapReadsToEdges();
	void mapReadLocations();
	int32_t* findDistanceOnEdge(Edge *edge, uint64_t read);
};

#endif /* MATEPAIR_H_ */
