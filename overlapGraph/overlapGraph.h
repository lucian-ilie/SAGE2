/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 ****************************************************************************/

#ifndef OVERLAPGRAPH_H_
#define OVERLAPGRAPH_H_

#include "../utils.h"
#include "../economyGraph/economyGraph.h"
#include "../inputReader/readLoader.h"

class ReadOnEdge
{
public:
	uint64_t readId:40;
	uint64_t orientation:1;
	uint64_t flag:1;
	uint64_t distPrevious:11;// distance form previous
	uint64_t distNext:11;// distance to next
};

class Edge
{
public:
	uint64_t ID:40;
	uint64_t fromID:40;
	uint64_t typeOfEdge:2;
	uint64_t isReducible:8;
	uint64_t lengthOfEdge:32;
	float flow;
	ReadOnEdge *listOfReads;
	Edge *next;
	Edge *previous;
	Edge *twinEdge;
	friend ostream& operator<<(ostream &os, const Edge *hList);
	friend istream& operator>>(istream &is, Edge * & hList);
};

class OverlapGraph
{
	EconomyGraph *economyObj;
	ReadLoader *loaderObj;
	uint64_t N_90;
	uint64_t N_50;
	uint64_t N_90_N;
	uint64_t N_50_N;
	uint64_t longestContig;
	uint64_t base_pair_covered;
	uint32_t contigLengthThreshold;
public:
	Edge **graph;
	OverlapGraph(EconomyGraph *economy1, ReadLoader *loader1);
	OverlapGraph(ReadLoader *loader1);
	~OverlapGraph();
	void convertGraph();
	void insertEdgeInGraph(uint64_t from, uint64_t to, uint32_t overlap_size, uint8_t type);
	Edge* insertEdgeSimple(uint64_t vertex_u, uint64_t vertex_v, uint8_t type, ReadOnEdge *listForward, ReadOnEdge *listReverse, float flow, uint32_t delta);
	Edge* insertEdge(uint64_t vertex_u, uint64_t vertex_v, uint8_t type, ReadOnEdge *listForward, ReadOnEdge *listReverse, float flow, uint32_t delta1, uint32_t delta2);
	void insertIntoList(Edge *v);
	Edge* mergeEdges(Edge *edge1, Edge *edge2, int type_of_merge);
	int combinedEdgeType(Edge *edge1, Edge* edge2);
	int deleteEdge(Edge *edge);
	ReadOnEdge* getListOfReads(Edge *edge1, Edge *edge2);
	void saveOverlapGraphInFile(string path);
	void loadOverlapGraphFromFile(string path);
	void checkIfAllReadsPresent();
	void printGraph(string fileName, bool isScaff=false);
	void saveContigsToFile(string fileName, Edge **contig_edges, uint64_t from, uint64_t to, bool isScaff=false);
	int getDegree(uint64_t v);
	ReadOnEdge* getListOfReadsInDisconnectedEdges(Edge *edge1, Edge *edge2, int gap, uint64_t *length);
	int stringOverlapSize(string string1, string string2, int readLength1, int readLength2);
};

#endif /* OVERLAPGRAPH_H_ */
