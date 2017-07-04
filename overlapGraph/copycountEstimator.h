/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 ****************************************************************************/

#ifndef COPYCOUNTESTIMATOR_H_
#define COPYCOUNTESTIMATOR_H_

#include "../utils.h"
#include "overlapGraph.h"
#include "cs2/cs2.h"

class CopyCountEstimator
{
	OverlapGraph *graphObj;
	ReadLoader *loaderObj;
	uint32_t aStatisticsThreshold;
	uint32_t minDelta;
public:
	uint64_t genomeSize;
	CopyCountEstimator(OverlapGraph *graph1, ReadLoader *loader1);
	~CopyCountEstimator();
	uint64_t genomeSizeEstimation();
	uint64_t findGenomeSize(uint64_t previousEstimation);
	void computeMinCostFlow(string outDir);
	void computeMinCostFlow_new(string outDir);
};

#endif /* COPYCOUNTESTIMATOR_H_ */
