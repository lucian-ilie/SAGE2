/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 ****************************************************************************/

#ifndef MERGECONTIGS_H_
#define MERGECONTIGS_H_

#include <sstream>
#include <deque>
#include <unordered_map>
#include <utility>
#include "../utils.h"
#include "matePair.h"
#include "../overlapGraph/simplification.h"

class Indexes
{
public:
	uint64_t firstNode;
	uint64_t lastNode;
	uint64_t position;
	Indexes *next;
};

class PairSupported
{
public:
	Edge *edge1;
	Edge *edge2;
};

class SupportedEdgePair
{
public:
	Edge *edge1;
	Edge *edge2;
	int support;
};

struct pairHash
{
public:
    size_t operator()(const pair<uint64_t, uint64_t> &k) const
	{
        return hash<uint64_t>()(k.first) ^ hash<uint64_t>()(k.second);
    }
};

class ContigExtender
{
	ReadLoader *loaderObj;
	MatePair *mateObj;
	OverlapGraph *graphObj;
	uint32_t contigLengthThreshold;
public:
	ContigExtender(OverlapGraph *graph1, MatePair *mate1, ReadLoader *loader1);
	~ContigExtender();
	void extendContigs(uint16_t minOverlap);
	uint64_t findPathSupports(int highTh);
	uint64_t exploreGraph(uint64_t node, int level, uint64_t &totalPaths, Edge **pathEdges, uint64_t *pathLengths, Edge ***allPairPaths, uint64_t *allPairPathLengths);
	int matchEdgeType(Edge *edge1, Edge *edge2);
	int checkFlow(Edge **edges, uint64_t top, Edge *next_edge);
	void savePath(Edge **edges, uint64_t distance, uint64_t top, uint64_t &totalPaths, Edge ***allPairPaths, uint64_t *allPairPathLengths);
	void quicksortPaths(int64_t left, int64_t right, Edge ***allPairPaths, uint64_t *allPairPathLengths);
	uint64_t partitionPaths(int64_t left, int64_t right, Edge ***allPairPaths, uint64_t *allPairPathLengths);
	void indexPaths(uint64_t num, Indexes **indx, uint64_t &totalPaths, Edge ***allPairPaths);
	int findPathsBetweenMatepairs(uint64_t mate_pair_1, uint64_t mate_pair_2, bool type1, bool type2, int library, Indexes **indx, uint64_t &totalPaths, Edge ***allPairPaths, uint64_t *allPairPathLengths, PairSupported *supportedPairs_tmp);
	uint64_t startingPosition(uint64_t firstReadOnEdge,uint64_t second_read, Indexes **indx);
	uint64_t insertPath(Edge **edges,int *support_flag, uint64_t top,uint64_t path_found, PairSupported *supportedPairs_tmp);
	int hashSupportInsert(PairSupported *supportedPairs, PairSupported sPair);
	void freeIndexedPaths(Indexes **indx);
	void quicksortListOfLongEdges(PairSupported *pairs, int *pair_support, int64_t left, int64_t right);
	uint64_t mergeByPathSupports(unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> &listSearch, unordered_multimap<uint64_t,SupportedEdgePair> &deleteList, int highTh);
	uint64_t findMatepairSupports(uint16_t minOverlap, int highTh);
	int getSupportOfEdges(Edge *edge1, Edge *edge2, uint16_t minOverlap);
	int isUniqueInEdge(Edge *edge,uint64_t read);
	uint64_t mergeByMatepairSupports(deque<SupportedEdgePair> &list, unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> &listSearch, unordered_multimap<uint64_t,SupportedEdgePair> &deleteList, int highTh);
	Edge* mergeEdgesAndUpdate(Edge *edge1, Edge *edge2, int type_of_merge);
	uint64_t mergeShortOverlaps(uint16_t minOverlap, int highTh);
	Edge** getListOfCompositeEdges1(uint64_t *number);
	Edge** getListOfFeasibleEdges1(Edge *edge, uint64_t *number);
	uint64_t mergeShortOverlap(Edge **edge_1,Edge **edge_2,int *pair_support,int32_t *gapDistance, uint64_t total_edge_pair_to_merge, unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> &listSearch, unordered_multimap<uint64_t,SupportedEdgePair> &deleteList, int highTh);
	Edge* mergeDisconnectedEdges1(Edge *edge1, Edge *edge2, int gap);
	void quicksortListOfLongEdges1(Edge **edge_1, Edge **edge_2, int * pair_support, int32_t *gapDistance, int64_t left, int64_t right);
	uint64_t mergeByPathSupportsNew(unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> &listSearch, unordered_multimap<uint64_t,SupportedEdgePair> &deleteList, int highTh);
	uint64_t findPathSupportsNew(int highTh);
};

#endif /* MERGECONTIGS_H_ */
