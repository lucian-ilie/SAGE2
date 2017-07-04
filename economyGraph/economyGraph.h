/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 ****************************************************************************/

#ifndef ECONOMYGRAPH_H_
#define ECONOMYGRAPH_H_

#include <parallel/algorithm>
#include "hashTable.h"
#include "../utils.h"

class EconomyEdge
{
public:
	uint64_t readId:40; //destination node
	uint64_t type:2;
	uint64_t mark:2;
	uint64_t length:20;
	EconomyEdge(uint64_t dest1, uint8_t type1, uint8_t mark1, uint32_t len1);
};

class ExtensionTable
{
public:
	uint64_t readId:40; // Stores the extension node ID.
	uint64_t type:2;	//	Stores orientation of the extension as forward (0) or reverse(1).
	uint64_t length:22; // Stores the length of the overhang. 
};

class EconomyGraph
{
	friend class OverlapGraph;
	uint8_t *exploredReads;
	uint8_t *markedNodes;
	EconomyEdge **economyGraphList;
	uint64_t numberOfHashMiss;
	uint16_t minOverlap;
	HashTable *hashObj;
	ReadLoader *loaderObj;
public:
	EconomyGraph(uint16_t minOvlp, HashTable *hash1);
	~EconomyGraph();
	void buildInitialOverlapGraph();
	void buildOverlapGraphEconomy();
	int insertAllEdgesOfRead(uint64_t read);
	void markTransitiveEdge(uint64_t from_ID);
	uint64_t removeTransitiveEdges(uint64_t read);
	int compareStringInBytes(uint8_t *read2bits1, uint8_t *read2bits2, uint16_t start, uint16_t readLength1, uint16_t readLength2, uint64_t read2ID);
	int compareStringInBytesPrevious(uint8_t *read2bits1, uint8_t *read2bits2, uint16_t start, uint16_t readLength1, uint16_t readLength2);
	int insertEdgeEconomy(uint64_t vertex_u, uint64_t vertex_v, uint32_t delta, uint8_t type);
	void sortEconomyGraph();
};

/* ============================================================================================
   This is the compare function which will be used for sorting using STL sort algorithm.
   This function will help to sort edges in ascending order of the length
   ============================================================================================ */
bool compareLengthBased(const EconomyEdge &first, const EconomyEdge &second);
bool compareIdBased(const EconomyEdge &first, const EconomyEdge &second);

#endif /* ECONOMYGRAPH_H_ */
