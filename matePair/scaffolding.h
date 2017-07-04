/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 ****************************************************************************/

#ifndef SCAFFOLDING_H_
#define SCAFFOLDING_H_

#include "../utils.h"
#include "../overlapGraph/simplification.h"
#include "mergeContigs.h"

class ScaffoldMaker
{
	ReadLoader *loaderObj;
	MatePair *mateObj;
	OverlapGraph *graphObj;
	uint32_t contigLengthThreshold;
public:
	ScaffoldMaker(OverlapGraph *graph1, MatePair *mate1, ReadLoader *loader1);
	~ScaffoldMaker();
	void makeScaffolds(uint16_t minOverlap);
	Edge** getListOfCompositeEdges(uint64_t *number);
	Edge** getListOfFeasibleEdges(Edge *edge, uint64_t *number);
	int isUniqueInEdge(Edge *edge,uint64_t read);
	void quicksortListOfLongEdges(Edge **edge_1, Edge **edge_2, int *pair_support, int32_t *gapDistance, int64_t left, int64_t right);
	Edge* mergeDisconnectedEdges(Edge *edge1,Edge *edge2,int gap);
	uint64_t mergeAccordingToSupport(Edge **edge_1,Edge **edge_2,int *pair_support,int32_t *gapDistance, uint64_t total_edge_pair_to_merge, unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> &listSearch, unordered_multimap<uint64_t,SupportedEdgePair> &deleteList, int highTh);
	uint64_t mergeAccordingToSupport2(Edge **edge_1,Edge **edge_2,int *pair_support,int32_t *gapDistance, uint64_t total_edge_pair_to_merge, unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> &listSearch, unordered_multimap<uint64_t,SupportedEdgePair> &deleteList, int highTh);
	uint64_t mergeFinal(uint16_t minOverlap, int highTh, int flag);
	uint64_t mergeFinalSmall(uint16_t minOverlap, int highTh, int flag);
};

#endif /* SCAFFOLDING_H_ */
