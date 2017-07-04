/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 ****************************************************************************/

#ifndef SIMPLIFICATION_H_
#define SIMPLIFICATION_H_

#include "../utils.h"
#include "overlapGraph.h"
#include "../matePair/matePair.h"

uint64_t contractCompositePaths(OverlapGraph *graphObj, ReadLoader *loaderObj);
uint64_t removeDeadEnds(OverlapGraph *graphObj, ReadLoader *loaderObj, int threshold);
uint64_t removeBubbles(OverlapGraph *graphObj, ReadLoader *loaderObj, int64_t closeLength);
uint64_t removeSimpleEdges(OverlapGraph *graphObj, ReadLoader *loaderObj);
uint64_t removeTransitiveEdges(OverlapGraph *graphObj, ReadLoader *loaderObj);

#endif /* SIMPLIFICATION_H_ */
