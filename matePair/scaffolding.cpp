/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 * Changes		: 
 ****************************************************************************/

#include "scaffolding.h"

extern ofstream logStream;
extern int lowThreshold, highThreshold;
extern uint64_t genomeSize, averageReadLength;
extern string outputDir, prefixName;

ScaffoldMaker::ScaffoldMaker(OverlapGraph *graph1, MatePair *mate1, ReadLoader *loader1)
{
	graphObj = graph1;
	mateObj = mate1;
	loaderObj = loader1;
	contigLengthThreshold = (int)ceilf(100*(averageReadLength/100.0));
}

ScaffoldMaker::~ScaffoldMaker()
{
	
}

/* ============================================================================================
   Merge edges not connected in the graph that are supported.
   ============================================================================================ */
void ScaffoldMaker::makeScaffolds(uint16_t minOverlap)
{
	time_t seconds_s=time(NULL);
	logStream<< "\nIn function makeScaffolds()." << endl;
	uint64_t merged=0, sum_merged=0;
	bool flg=true;
	int flag=1;
	float coverage = (loaderObj->numberOfReads*averageReadLength)/genomeSize;
	int highTh = (int)ceilf((0.10*coverage)*(100.0/averageReadLength));

	for(int iteration=1; flg; iteration++)
	{
		if(iteration>5)
			highTh = (int)ceilf((0.05*coverage)*(100.0/averageReadLength));
		merged=mergeFinal(minOverlap, highTh, flag);
		sum_merged += merged;
		if(merged==0 && highTh==(int)ceilf((0.05*coverage)*(100.0/averageReadLength)))
			flg=false;;
	}
	
	flag=2;
	flg=true;
	highTh = (int)ceilf((0.10*coverage)*(100.0/averageReadLength));
	for(int iteration=1; flg; iteration++)
	{
		if(iteration>5)
			highTh = (int)ceilf((0.05*coverage)*(100.0/averageReadLength));
		merged=mergeFinal(minOverlap, highTh, flag);
		sum_merged += merged;
		if(merged==0 && highTh==(int)ceilf((0.05*coverage)*(100.0/averageReadLength)))
			flg=false;;
	}
	
	flag=1;
	flg=true;
	highTh = 2;
	for(int iteration=1; flg; iteration++)
	{
		merged=mergeFinal(minOverlap, highTh, flag);
		sum_merged += merged;
		if(merged==0)
			flg=false;;
	}
	
	flg=true;
	highTh = 1;
	for(int iteration=1; flg; iteration++)
	{
		merged=mergeFinal(minOverlap, highTh, flag);
		sum_merged += merged;
		if(merged==0)
			flg=false;;
	}
	
	flag=2;
	flg=true;
	highTh = 1;
	for(int iteration=1; flg; iteration++)
	{
		merged=mergeFinalSmall(minOverlap, highTh, flag);
		sum_merged += merged;
		if(merged==0)
			flg=false;;
	}
	logStream<< "\nTotal number of merged scaffolds: " << sum_merged << "\n";
	logStream<< "Function makeScaffolds() in " << time(NULL)-seconds_s << " sec.\n";

}


/* ============================================================================================
   This function returns the list of all composite edges in the graph.
   ============================================================================================ */
Edge** ScaffoldMaker::getListOfCompositeEdges(uint64_t *number)
{
	Edge **listOfLongEdges, *v, *u;
	uint64_t arraySize=100, tos=0, i;
	if((listOfLongEdges=(Edge**)malloc(arraySize*sizeof(Edge*)))==NULL)
		printError(MEM_ALLOC, "listOfLongEdges");

	for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
	{
		if(graphObj->graph[i]!=NULL)
		{
			for(v=graphObj->graph[i]; v!=NULL; v=v->next)
			{
					if(i<=v->ID && v->listOfReads[0].readId!=0) //Take all composite edges.
				{
					if(i<v->ID)
						listOfLongEdges[tos++]=v;
					else
					{
						for(u=graphObj->graph[i]; u!=v; u=u->next)
						{
							if(u==v->twinEdge)// Twin edge already considered
								break;
						}
						if(u==v) // Not considered already
							listOfLongEdges[tos++]=v;
					}
					if(tos>=arraySize-5)
					{
						arraySize=arraySize<<1;
						if((listOfLongEdges=(Edge**)realloc(listOfLongEdges, arraySize*sizeof(Edge*)))==NULL)
							printError(MEM_ALLOC, "reallocating listOfLongEdges");
					}
				}
			}
		}
	}
	*number=tos;
	return listOfLongEdges;
}

/* ============================================================================================
   This function returns the list of edges that can possibly be connected to 'edge' by matepairs.
   ============================================================================================ */
Edge** ScaffoldMaker::getListOfFeasibleEdges(Edge *edge, uint64_t *number)
{
	Edge** feasibleListOfEdges;
	MatePairInfo *w;
	uint64_t k, mate_pair_1, mate_pair_2, arraySize=10, numberOfFeasibleEdges=0, dist=0, index;
	if((feasibleListOfEdges=(Edge**)malloc(arraySize*sizeof(Edge*)))==NULL)
		printError(MEM_ALLOC, "feasibleListOfEdges");
	for(index=1; index<=edge->listOfReads[0].readId; index++)
	{
		dist += edge->listOfReads[index].distPrevious;
		if(dist>mateObj->maximumUpperBoundOfInsert)
			 break;
		mate_pair_1=edge->listOfReads[index].readId;
		if(isUniqueInEdge(edge, mate_pair_1))
		{
			for(w=mateObj->matePairList[mate_pair_1]; w!=NULL; w=w->next)
			{
				mate_pair_2=w->ID;
				if(mateObj->readToEdgeList[mate_pair_2]==NULL)
					continue;
				if(mateObj->readToEdgeList[mate_pair_2]->edge==edge)
					continue;
				if(mateObj->readToEdgeList[mate_pair_2]->edge==edge->twinEdge)
					continue;
				if(isUniqueInEdge(mateObj->readToEdgeList[mate_pair_2]->edge, mate_pair_2))
				{
					for(k=0; k<numberOfFeasibleEdges; k++)
						if(feasibleListOfEdges[k]==mateObj->readToEdgeList[mate_pair_2]->edge || feasibleListOfEdges[k]==mateObj->readToEdgeList[mate_pair_2]->edge->twinEdge)
							break;
					if(k==numberOfFeasibleEdges)
					{
						feasibleListOfEdges[numberOfFeasibleEdges++]=mateObj->readToEdgeList[mate_pair_2]->edge;
						if(numberOfFeasibleEdges>=arraySize-5)
						{
							arraySize=arraySize<<1;
							if((feasibleListOfEdges=(Edge**)realloc(feasibleListOfEdges, arraySize*sizeof(Edge*)))==NULL)
								printError(MEM_ALLOC, "realloc feasibleListOfEdges");
						}
					}
				}
			}
		}
	}
	dist=0;
	for(index=1; index<=edge->twinEdge->listOfReads[0].readId; index++)
	{
		dist += edge->twinEdge->listOfReads[index].distPrevious;
		if(dist>mateObj->maximumUpperBoundOfInsert)
			 break;
		mate_pair_1=edge->twinEdge->listOfReads[index].readId;
		if(isUniqueInEdge(edge->twinEdge, mate_pair_1))
		{
			for(w=mateObj->matePairList[mate_pair_1]; w!=NULL; w=w->next)
			{
				mate_pair_2=w->ID;
				if(mateObj->readToEdgeList[mate_pair_2]==NULL)
					continue;
				if(mateObj->readToEdgeList[mate_pair_2]->edge==edge)
					continue;
				if(mateObj->readToEdgeList[mate_pair_2]->edge==edge->twinEdge)
					continue;
				if(isUniqueInEdge(mateObj->readToEdgeList[mate_pair_2]->edge, mate_pair_2))
				{
					for(k=0; k<numberOfFeasibleEdges; k++)
					{
						if(feasibleListOfEdges[k]==mateObj->readToEdgeList[mate_pair_2]->edge || feasibleListOfEdges[k]==mateObj->readToEdgeList[mate_pair_2]->edge->twinEdge)
						{
							break;
						}
					}
					if(k==numberOfFeasibleEdges)
					{
						feasibleListOfEdges[numberOfFeasibleEdges++]=mateObj->readToEdgeList[mate_pair_2]->edge;
						if(numberOfFeasibleEdges>=arraySize-5)
						{
							arraySize=arraySize<<1;
							if((feasibleListOfEdges=(Edge**)realloc(feasibleListOfEdges, arraySize*sizeof(Edge*)))==NULL)
								printError(MEM_ALLOC, "realloc feasibleListOfEdges");
						}
					}
				}
			}
		}
	}
	*number=numberOfFeasibleEdges;
	return feasibleListOfEdges;
}

/* ============================================================================================
   This function checks if read is present in edge.
   ============================================================================================ */
int ScaffoldMaker::isUniqueInEdge(Edge *edge,uint64_t read)
{
	if(mateObj->readToEdgeList[read]==NULL) //not present;
		return 0;

	if((mateObj->readToEdgeList[read]->edge==edge || mateObj->readToEdgeList[read]->edge->twinEdge==edge) && mateObj->readToEdgeList[read]->next==NULL && graphObj->graph[read]==NULL)
		return 1;
	return 0;
}

/* ============================================================================================
   Sort list of long edges according to support. If the supports are same then sort by length.
   ============================================================================================ */
void ScaffoldMaker::quicksortListOfLongEdges(Edge **edge_1, Edge **edge_2, int * pair_support, int32_t *gapDistance, int64_t left, int64_t right)
{
	int64_t i=left,j=right,pivot=*(pair_support+((left+right)>>1)),temporary,lengthSum=edge_1[(left+right)>>1]->lengthOfEdge+edge_2[(left+right)>>1]->lengthOfEdge;
	Edge *temp;
	while(i<j)
	{
		while((int64_t)*(pair_support+i)>pivot || ((int64_t)*(pair_support+i)==pivot && ((int64_t)(edge_1[i]->lengthOfEdge+edge_2[i]->lengthOfEdge))>lengthSum))
			i++;
		while((int64_t)*(pair_support+j)<pivot || ((int64_t)*(pair_support+j)==pivot && ((int64_t)(edge_1[j]->lengthOfEdge+edge_2[j]->lengthOfEdge))<lengthSum))
			j--;
		if (i<=j)
		{
			temp=*(edge_1+i); *(edge_1+i)=*(edge_1+j); *(edge_1+j)=temp;
			temp=*(edge_2+i); *(edge_2+i)=*(edge_2+j); *(edge_2+j)=temp;
			temporary=*(pair_support+i); *(pair_support+i)=*(pair_support+j); *(pair_support+j)=temporary;
			temporary=*(gapDistance+i); *(gapDistance+i)=*(gapDistance+j); *(gapDistance+j)=temporary;
			i++; j--;
		}
	}
	if (left < j )
		quicksortListOfLongEdges(edge_1,edge_2,pair_support,gapDistance,left,j);
    if (i < right)
		quicksortListOfLongEdges(edge_1,edge_2,pair_support,gapDistance,i,right);
}


/* ============================================================================================
   Two edges (e1, e2) and (e2,e3) is merged to (e1,e3) here.
   ============================================================================================ */
Edge* ScaffoldMaker::mergeDisconnectedEdges(Edge *edge1, Edge *edge2, int gap)
{
	int type=0;
	Edge *mergedEdge, *v;
	uint64_t index;
	uint32_t dist;
	int flag, i;
	ReadToEdgeMap *edg, *edg_prev, *edg_next;
	uint64_t length1, length2;
	if((edge1->typeOfEdge==0 || edge1->typeOfEdge==1) && (edge2->typeOfEdge==0 || edge2->typeOfEdge==2))
		type=0;
	else if((edge1->typeOfEdge==0 || edge1->typeOfEdge==1) && (edge2->typeOfEdge==1 || edge2->typeOfEdge==3))
		type=1;
	else if((edge1->typeOfEdge==2 || edge1->typeOfEdge==3) && (edge2->typeOfEdge==0 || edge2->typeOfEdge==2))
		type=2;
	else if((edge1->typeOfEdge==2 || edge1->typeOfEdge==3) && (edge2->typeOfEdge==1 || edge2->typeOfEdge==3))
		type=3;
	ReadOnEdge *returnList1=graphObj->getListOfReadsInDisconnectedEdges(edge1, edge2, gap, &length1);
	ReadOnEdge *returnList2=graphObj->getListOfReadsInDisconnectedEdges(edge2->twinEdge, edge1->twinEdge, gap, &length2);
	mergedEdge=graphObj->insertEdge(edge1->fromID, edge2->ID, type, returnList1, returnList2, 1.0, length1, length2);
	if(edge1->flow<=1.0)
	{
		//remove previous mappings
		for(index=1; index<=edge1->listOfReads[0].readId; index++)
		{
			if((edg=mateObj->readToEdgeList[edge1->listOfReads[index].readId])!=NULL)
			{
				if(edg->edge==edge1 || edg->edge==edge1->twinEdge) //First node
				{
					edg_next=edg->next;
					free(edg->locationForward);
					free(edg->locationReverse);
					free(edg);
					mateObj->readToEdgeList[edge1->listOfReads[index].readId] = edg_next;
				}
				else
				{
					edg_prev=mateObj->readToEdgeList[edge1->listOfReads[index].readId];
					edg = edg_prev->next;
					for(; edg!=NULL; edg_prev=edg, edg=edg_next)
					{
						edg_next=edg->next;
						if(edg->edge==edge1 || edg->edge==edge1->twinEdge)
						{
							free(edg->locationForward);
							free(edg->locationReverse);
							free(edg);
							edg_prev->next = edg_next;
							break;
						}
					}
				}
			}
		}

		graphObj->deleteEdge(edge1->twinEdge);
		graphObj->deleteEdge(edge1);
	}
	else
	{
		edge1->twinEdge->flow-=1.0;
		edge1->flow-=1.0;
	}
	if(edge2->flow<=1.0)
	{
		//remove previous mappings
		for(index=1; index<=edge2->listOfReads[0].readId; index++)
		{
			if((edg=mateObj->readToEdgeList[edge2->listOfReads[index].readId])!=NULL)
			{
				if(edg->edge==edge2 || edg->edge==edge2->twinEdge) //First node
				{
					edg_next=edg->next;
					free(edg->locationForward);
					free(edg->locationReverse);
					free(edg);
					mateObj->readToEdgeList[edge2->listOfReads[index].readId] = edg_next;
				}
				else
				{
					edg_prev=mateObj->readToEdgeList[edge2->listOfReads[index].readId];
					edg = edg_prev->next;
					for(; edg!=NULL; edg_prev=edg, edg=edg_next)
					{
						edg_next=edg->next;
						if(edg->edge==edge2 || edg->edge==edge2->twinEdge)
						{
							free(edg->locationForward);
							free(edg->locationReverse);
							free(edg);
							edg_prev->next = edg_next;
							break;
						}
					}
				}
			}
		}

		graphObj->deleteEdge(edge2->twinEdge);
		graphObj->deleteEdge(edge2);
	}
	else
	{
		edge2->twinEdge->flow-=1.0;
		edge2->flow-=1.0;
	}

	//update readToEdge and locations
	if(mergedEdge->fromID <= mergedEdge->ID)
		v = mergedEdge;
	else
		v = mergedEdge->twinEdge;

	for(index=1; index<=v->listOfReads[0].readId; index++)
	{
		flag=1;
		for(edg=mateObj->readToEdgeList[v->listOfReads[index].readId]; edg!=NULL; edg=edg->next)
		{
			if(edg->edge==v || edg->edge==v->twinEdge)
			{
				flag=0; //already present
				break;
			}
		}
		if(flag==1) // Not present in the list already.
		{
			if((edg=(ReadToEdgeMap*)malloc(sizeof(ReadToEdgeMap)))==NULL)
				printError(MEM_ALLOC, "read to edge");
			edg->edge=v;
			edg->locationForward=NULL;
			edg->locationReverse=NULL;
			edg->next=mateObj->readToEdgeList[v->listOfReads[index].readId];
			mateObj->readToEdgeList[v->listOfReads[index].readId]=edg;
		}
	}
	//
	for(i=0; i<2; i++)
	{
		if(i==0)
			v = mergedEdge;
		if(i==1)
			v = mergedEdge->twinEdge;
		dist=0;
		for(index=1; index<=v->listOfReads[0].readId; index++)
		{
			dist+=v->listOfReads[index].distPrevious;
			for(edg=mateObj->readToEdgeList[v->listOfReads[index].readId]; edg!=NULL; edg=edg->next)
			{
				if(edg->edge==v)
				{
					if(edg->locationForward==NULL) // No locationForward mapped yet. Should allocate the memory for the first time.
					{
						if((edg->locationForward=(uint32_t*)malloc((2)*sizeof(uint32_t)))==NULL)
							printError(MEM_ALLOC, "locationForward");
						edg->locationForward[0]=1;
						if(v->listOfReads[index].orientation)
							edg->locationForward[edg->locationForward[0]]=dist;
						else
							edg->locationForward[edg->locationForward[0]]=-dist;
					}
					else // Increase the allocate memory to make space for the new locationForward.
					{
						if((edg->locationForward=(uint32_t*)realloc(edg->locationForward, (edg->locationForward[0]+2)*sizeof(uint32_t)))==NULL)
							printError(MEM_ALLOC, "Reallocating locationForward");
						edg->locationForward[0]++;
						if(v->listOfReads[index].orientation)
							edg->locationForward[edg->locationForward[0]]=dist;
						else
							edg->locationForward[edg->locationForward[0]]=-dist;
					}
					break;
				}
				else if(edg->edge==v->twinEdge)
				{
					if(edg->locationReverse==NULL)// No locationReverse mapped yet. Should allocate the memory for the first time.
					{
						if((edg->locationReverse=(uint32_t*)malloc((2)*sizeof(uint32_t)))==NULL)
							printError(MEM_ALLOC, "locationReverse");
						edg->locationReverse[0]=1;
						if(v->listOfReads[index].orientation) // Check the orientation
							edg->locationReverse[edg->locationReverse[0]]=dist;
						else
							edg->locationReverse[edg->locationReverse[0]]=-dist;
					}
					else // Increase the allocate memory to make space for the new locationReverse.
					{
						if((edg->locationReverse=(uint32_t*)realloc(edg->locationReverse, (edg->locationReverse[0]+2)*sizeof(uint32_t)))==NULL)
							printError(MEM_ALLOC, "Reallocating locationReverse");
						edg->locationReverse[0]++;
						if(v->listOfReads[index].orientation) // Check the orientation
							edg->locationReverse[edg->locationReverse[0]]=dist;
						else
							edg->locationReverse[edg->locationReverse[0]]=-dist;
					}
					break;
				}
			}
		}
	}

	return mergedEdge;
}


/* ============================================================================================
   Merge edge pairs according to support.
   ============================================================================================ */
uint64_t ScaffoldMaker::mergeAccordingToSupport(Edge **edge_1,Edge **edge_2,int *pair_support,int32_t *gapDistance, uint64_t total_edge_pair_to_merge, unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> &listSearch, unordered_multimap<uint64_t,SupportedEdgePair> &deleteList, int highTh)
{
	uint64_t number_merged=0, edgeID=0, cycleFound=0;
	uint64_t i, j, k, firstEdge=0, secondEdge=0;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash>::const_iterator it;
	Edge **listOfEdges1, **listOfEdges2;
	uint64_t numOfList1, numOfList2, unique1Count=0, unique2Count=0;
	int *suppList1, *suppList2;
	bool cycleFlag;
	
	for(i=0;i<total_edge_pair_to_merge;i++) // merge according to support
	{
		if(edge_1[i]==NULL || edge_2[i]==NULL)
			continue;
		firstEdge=edge_1[i]->fromID ^ edge_1[i]->ID;
		secondEdge=edge_2[i]->fromID ^ edge_2[i]->ID;
		it = listSearch.find(make_pair(firstEdge, secondEdge));
		if(it == listSearch.end())
			continue;

		if(pair_support[i] >= highTh)
		{
			listOfEdges1 = getListOfFeasibleEdges(edge_2[i], &numOfList1);
			listOfEdges2 = getListOfFeasibleEdges(edge_1[i], &numOfList2);
			suppList1=NULL;
			suppList2=NULL;
			cycleFlag=0;
			if(numOfList1>0)
			{
				if((suppList1=(int*)malloc(numOfList1*sizeof(int)))==NULL)
					printError(MEM_ALLOC, "suppList1");
				for(j=0; j<numOfList1; j++)
				{
					suppList1[j]=0;
					firstEdge=listOfEdges1[j]->fromID ^ listOfEdges1[j]->ID;
					secondEdge=edge_2[i]->fromID ^ edge_2[i]->ID;					
					it = listSearch.find(make_pair(firstEdge, secondEdge));
					if(it != listSearch.end())
						suppList1[j]=it->second.support;
					else
					{						
						firstEdge=edge_2[i]->twinEdge->fromID ^ edge_2[i]->twinEdge->ID;
						secondEdge=listOfEdges1[j]->twinEdge->fromID ^ listOfEdges1[j]->twinEdge->ID;	
						it = listSearch.find(make_pair(firstEdge, secondEdge));
						if(it != listSearch.end())
							suppList1[j]=it->second.support;
					}
				}
			}
			if(numOfList2>0)
			{
				if((suppList2=(int*)malloc(numOfList2*sizeof(int)))==NULL)
					printError(MEM_ALLOC, "suppList2");
				for(j=0; j<numOfList2; j++)
				{
					suppList2[j]=0;
					firstEdge=edge_1[i]->fromID ^ edge_1[i]->ID;
					secondEdge=listOfEdges2[j]->fromID ^ listOfEdges2[j]->ID;
					it = listSearch.find(make_pair(firstEdge, secondEdge));
					if(it != listSearch.end())
						suppList2[j]=it->second.support;
					else
					{	
						firstEdge=listOfEdges2[j]->twinEdge->fromID ^ listOfEdges2[j]->twinEdge->ID;
						secondEdge=edge_1[i]->twinEdge->fromID ^ edge_1[i]->twinEdge->ID;
						it = listSearch.find(make_pair(firstEdge, secondEdge));
						if(it != listSearch.end())
							suppList2[j]=it->second.support;
					}
				}
			}
			unique1Count=0;
			unique2Count=0;
			
			//check second list
			for(j=0; j<numOfList2; j++)
			{
				if(listOfEdges2[j]==edge_2[i])
					continue;
				if(suppList2[j]>=highTh)
					unique2Count++;
			}
			//check in edges
			for(j=0; j<numOfList1; j++)
			{
				if(listOfEdges1[j]==edge_1[i])
					continue;
				if(suppList1[j]>=highTh)
					unique1Count++;
			}
			
			// check for cycles
			for(j=0; j<numOfList2; j++)
			{
				if(listOfEdges2[j]==edge_2[i])
					continue;
				if(suppList2[j]>=highTh)
				{
					for(k=0; k<numOfList1; k++)
					{
						if(listOfEdges1[k]==edge_1[i])
							continue;
						if(suppList1[k]>=highTh)
						{	
							if(listOfEdges2[j] == listOfEdges1[k])
							{	
								cycleFound++;
								cycleFlag=1;
							}
						}
					}				
				}
			}
			
			if(unique1Count<=1 && unique2Count<=1 && cycleFlag==0)
			{
				number_merged++;
				if(mergeDisconnectedEdges(edge_1[i], edge_2[i], gapDistance[i]))
				{
					edgeID=edge_1[i]->ID;
					auto its = deleteList.equal_range(edgeID);
					for (auto it = its.first; it != its.second; ++it) 
					{	
						firstEdge=it->second.edge1->fromID ^ it->second.edge1->ID;
						secondEdge=it->second.edge2->fromID ^ it->second.edge2->ID;
						listSearch.erase(make_pair(firstEdge, secondEdge));
					}
					edge_1[i]=NULL;
					edge_2[i]=NULL;
				}
			}
			free(listOfEdges1);
			free(listOfEdges2);
			free(suppList1);
			free(suppList2);
		}
	}
	return number_merged;
}

/* ============================================================================================
   Merge edge pairs according to support.
   ============================================================================================ */
uint64_t ScaffoldMaker::mergeAccordingToSupport2(Edge **edge_1,Edge **edge_2,int *pair_support,int32_t *gapDistance, uint64_t total_edge_pair_to_merge, unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> &listSearch, unordered_multimap<uint64_t,SupportedEdgePair> &deleteList, int highTh)
{
	uint64_t number_merged=0, edgeID=0, cycleFound=0;
	uint64_t i, j, k, firstEdge=0, secondEdge=0;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash>::const_iterator it;
	Edge **listOfEdges1, **listOfEdges2;
	uint64_t numOfList1, numOfList2, unique1Count=0, unique2Count=0;
	int *suppList1, *suppList2;
	bool cycleFlag;
	for(i=0;i<total_edge_pair_to_merge;i++) // merge according to support
	{
		if(edge_1[i]==NULL || edge_2[i]==NULL)
			continue;
		firstEdge=edge_1[i]->fromID ^ edge_1[i]->ID;
		secondEdge=edge_2[i]->fromID ^ edge_2[i]->ID;
		it = listSearch.find(make_pair(firstEdge, secondEdge));
		if(it == listSearch.end())
			continue;

		if(pair_support[i] >= highTh)
		{
			listOfEdges1 = getListOfFeasibleEdges(edge_2[i], &numOfList1);
			listOfEdges2 = getListOfFeasibleEdges(edge_1[i], &numOfList2);
			suppList1=NULL;
			suppList2=NULL;
			cycleFlag=0;
			if(numOfList1>0)
			{
				if((suppList1=(int*)malloc(numOfList1*sizeof(int)))==NULL)
					printError(MEM_ALLOC, "suppList1");
				for(j=0; j<numOfList1; j++)
				{
					suppList1[j]=0;
					firstEdge=listOfEdges1[j]->fromID ^ listOfEdges1[j]->ID;
					secondEdge=edge_2[i]->fromID ^ edge_2[i]->ID;					
					it = listSearch.find(make_pair(firstEdge, secondEdge));
					if(it != listSearch.end())
						suppList1[j]=it->second.support;
					else
					{						
						firstEdge=edge_2[i]->twinEdge->fromID ^ edge_2[i]->twinEdge->ID;
						secondEdge=listOfEdges1[j]->twinEdge->fromID ^ listOfEdges1[j]->twinEdge->ID;	
						it = listSearch.find(make_pair(firstEdge, secondEdge));
						if(it != listSearch.end())
							suppList1[j]=it->second.support;
					}
				}
			}
			if(numOfList2>0)
			{
				if((suppList2=(int*)malloc(numOfList2*sizeof(int)))==NULL)
					printError(MEM_ALLOC, "suppList2");
				for(j=0; j<numOfList2; j++)
				{
					suppList2[j]=0;
					firstEdge=edge_1[i]->fromID ^ edge_1[i]->ID;
					secondEdge=listOfEdges2[j]->fromID ^ listOfEdges2[j]->ID;
					it = listSearch.find(make_pair(firstEdge, secondEdge));
					if(it != listSearch.end())
						suppList2[j]=it->second.support;
					else
					{	
						firstEdge=listOfEdges2[j]->twinEdge->fromID ^ listOfEdges2[j]->twinEdge->ID;
						secondEdge=edge_1[i]->twinEdge->fromID ^ edge_1[i]->twinEdge->ID;
						it = listSearch.find(make_pair(firstEdge, secondEdge));
						if(it != listSearch.end())
							suppList2[j]=it->second.support;
					}
				}
			}
			unique1Count=0;
			unique2Count=0;
			
			//check second list
			for(j=0; j<numOfList2; j++)
			{
				if(listOfEdges2[j]==edge_2[i])
					continue;
				if(suppList2[j]>=highTh)
					unique2Count++;
			}
			//check in edges
			for(j=0; j<numOfList1; j++)
			{
				if(listOfEdges1[j]==edge_1[i])
					continue;
				if(suppList1[j]>=highTh)
					unique1Count++;
			}
			
			// check for cycles
			for(j=0; j<numOfList2; j++)
			{
				if(listOfEdges2[j]==edge_2[i])
					continue;
				if(suppList2[j]>=highTh)
				{
					for(k=0; k<numOfList1; k++)
					{
						if(listOfEdges1[k]==edge_1[i])
							continue;
						if(suppList1[k]>=highTh)
						{	
							if(listOfEdges2[j] == listOfEdges1[k])
							{	
								cycleFound++;
								cycleFlag=1;
							}
						}
					}				
				}
			}

			if(unique1Count<=2 && unique2Count<=2 && cycleFlag==0)
			{
				number_merged++;

				if(mergeDisconnectedEdges(edge_1[i], edge_2[i], gapDistance[i]))
				{
					edgeID=edge_1[i]->ID;
					auto its = deleteList.equal_range(edgeID);
					for (auto it = its.first; it != its.second; ++it) 
					{	
						firstEdge=it->second.edge1->fromID ^ it->second.edge1->ID;
						secondEdge=it->second.edge2->fromID ^ it->second.edge2->ID;
						listSearch.erase(make_pair(firstEdge, secondEdge));
					}
					edge_1[i]=NULL;
					edge_2[i]=NULL;
				}
			}
			free(listOfEdges1);
			free(listOfEdges2);
			free(suppList1);
			free(suppList2);
		}
	}
	return number_merged;
}

/* ============================================================================================
   This function merges long contigs that are not connected by overlaps.
   ============================================================================================ */
uint64_t ScaffoldMaker::mergeFinal(uint16_t minOverlap, int highTh, int flag)
{
	int type1a=0, type2a=0;
	uint64_t i, j, k, m, n, tos=0, number_merged=0, mate_pair_1, mate_pair_2, dist, total_edge_pair_to_merge=0, numberOfFeasibleEdges, arraySize=10000, flag_distance, index, edgeOne=0, edgeTwo=0, edgeThree=0, edgeFour=0;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> supportedListSearch;
	unordered_multimap<uint64_t,SupportedEdgePair> deleteList;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash>::const_iterator it;
	uint64_t firstEdge=0, secondEdge=0; // for calculating amount to subtract from memory usage
	int32_t *gapDistance, gapSum, support, distance_all, *distance_1, *distance_2;
	int *pair_support;
	Edge **listOfLongEdges, *edge1=NULL, *edge2=NULL, **edge_1, **edge_2, **feasibleListOfEdges;
	MatePairInfo *w;
	listOfLongEdges=getListOfCompositeEdges(&tos);
	if((edge_1=(Edge**)malloc(arraySize*sizeof(Edge*)))==NULL)
		printError(MEM_ALLOC, "edge_1");
	if((edge_2=(Edge**)malloc(arraySize*sizeof(Edge*)))==NULL)
		printError(MEM_ALLOC, "edge_2");
	if((pair_support=(int*)malloc(arraySize*sizeof(int)))==NULL)
		printError(MEM_ALLOC, "pair_support");
	if((gapDistance=(int32_t*)malloc(arraySize*sizeof(int32_t)))==NULL)
		printError(MEM_ALLOC, "gapDistance");
	for(i=0; i<tos; i++)
	{
		feasibleListOfEdges=getListOfFeasibleEdges(listOfLongEdges[i], &numberOfFeasibleEdges);
		for(j=0; j<numberOfFeasibleEdges; j++)
		{
			int64_t len1=listOfLongEdges[i]->lengthOfEdge, len2=feasibleListOfEdges[j]->lengthOfEdge;
			if((len1>200 && len2>200) || ((listOfLongEdges[i]->flow>0 && len1>contigLengthThreshold) && (feasibleListOfEdges[j]->flow>0 && len2>contigLengthThreshold)))
			{
				for(k=0; k<4; k++)
				{
					if(k==0)
					{
						edge1=listOfLongEdges[i];
						edge2=feasibleListOfEdges[j];
					}
					else if(k==1)
					{
						edge1=listOfLongEdges[i];
						edge2=feasibleListOfEdges[j]->twinEdge;
					}
					else if(k==2)
					{
						edge1=listOfLongEdges[i]->twinEdge;
						edge2=feasibleListOfEdges[j];
					}
					else if(k==3)
					{
						edge1=listOfLongEdges[i]->twinEdge;
						edge2=feasibleListOfEdges[j]->twinEdge;
					}
					dist=0;
					support=0;
					gapSum=0;
					if(edge1>edge2) 
						continue;
						
					for(index=1; index<=edge1->listOfReads[0].readId; index++)
					{
						dist += edge1->listOfReads[index].distPrevious;
						if(dist>mateObj->maximumUpperBoundOfInsert)
							break;
						mate_pair_1=edge1->listOfReads[index].readId;
						if(isUniqueInEdge(edge1, mate_pair_1)) // mate pair1 is only present in this edge
						{
							for(w=mateObj->matePairList[mate_pair_1]; w!=NULL; w=w->next)
							{
								mate_pair_2=w->ID;
								int lib=w->library;
								if(isUniqueInEdge(edge2, mate_pair_2)) // matepair 2 is only present in this edge
								{
									distance_1=mateObj->findDistanceOnEdge(edge1, mate_pair_1);
									distance_2=mateObj->findDistanceOnEdge(edge2, mate_pair_2);
									if(distance_1[0]>0 && distance_2[0]>0)
									{
										flag_distance=0;
										for(m=1; m<=(uint64_t)distance_1[0]; m++)
										{
											for(n=1; n<=(uint64_t)distance_2[0]; n++)
											{
												distance_all=10000000;
												if(w->type1==1)
													type1a=1;
												else
													type1a=-1;
												if(w->type2==1)
													type2a=1;
												else
													type2a=-1;

												if(type1a*distance_1[m]<0 && type2a*distance_2[n]<0)
													distance_all=llabs(distance_1[m])+llabs(distance_2[n]);
												if(distance_all+2*(uint64_t)averageReadLength<(uint64_t)(mateObj->upperBoundOfInsert[lib]+minOverlap))
												{
													flag_distance=1;
													break;
												}
											}
											if(flag_distance==1)
												break;
										}
										if(flag_distance==1)
										{
											support++;
											gapSum+=mateObj->Mean[lib]-(distance_all+2*averageReadLength);
										}
									}
									free(distance_1);
									free(distance_2);
								}
							}
						}
					}
					if(support>0)
					{
						edge_1[total_edge_pair_to_merge]=edge1->twinEdge;
						edge_2[total_edge_pair_to_merge]=edge2;
						pair_support[total_edge_pair_to_merge]=support;
						gapDistance[total_edge_pair_to_merge++]=gapSum/support;						
						SupportedEdgePair pairSupp;
						pairSupp.edge1 = edge1->twinEdge;
						pairSupp.edge2 = edge2;
						firstEdge = pairSupp.edge1->fromID ^ pairSupp.edge1->ID;
						secondEdge = pairSupp.edge2->fromID ^ pairSupp.edge2->ID;
						edgeOne = pairSupp.edge1->ID;
						edgeTwo = pairSupp.edge1->fromID;
						edgeThree = pairSupp.edge2->ID;
						edgeFour = pairSupp.edge2->fromID;

						pairSupp.support = support;
						it = supportedListSearch.find(make_pair(firstEdge, secondEdge));
						if(it != supportedListSearch.end())
						{	
							if(support > it->second.support)
								supportedListSearch[make_pair(firstEdge, secondEdge)] = pairSupp;
						}
						else
							supportedListSearch[make_pair(firstEdge, secondEdge)] = pairSupp;

						deleteList.insert({edgeOne,pairSupp});
						deleteList.insert({edgeTwo,pairSupp});
						deleteList.insert({edgeThree,pairSupp});
						deleteList.insert({edgeFour,pairSupp});	
							
						if(total_edge_pair_to_merge>=arraySize-5)
						{
							arraySize=arraySize<<1;
							if((edge_1=(Edge**)realloc(edge_1, arraySize*sizeof(Edge*)))==NULL)
								printError(MEM_ALLOC, "realloc edge_1");
							if((edge_2=(Edge**)realloc(edge_2, arraySize*sizeof(Edge*)))==NULL)
								printError(MEM_ALLOC, "realloc edge_2");
							if((pair_support=(int*)realloc(pair_support, arraySize*sizeof(int)))==NULL)
								printError(MEM_ALLOC, "realloc pair_support");
							if((gapDistance=(int32_t*)realloc(gapDistance, arraySize*sizeof(int32_t)))==NULL)
								printError(MEM_ALLOC, "realloc gapDistance");
						}
					}
				}
			}
		}
		free(feasibleListOfEdges);
	}
	
	uint64_t samllGap=0, totalGaps=0, gapDistanceTotal=0;
	int highGap=0;
	for(k=0; k<arraySize; k++)
	{
		if(gapDistance[k]<1000 && gapDistance[k]>0)
		{
			gapDistanceTotal += gapDistance[k];
			totalGaps++;
			if(highGap < gapDistance[k])
				highGap = gapDistance[k];
			if(gapDistance[k]<=10)
				samllGap++;
		}
	}
	
	if(total_edge_pair_to_merge>0)
	{
		quicksortListOfLongEdges(edge_1, edge_2, pair_support, gapDistance, 0, total_edge_pair_to_merge-1); // Sort edges according to the support.
		if(flag==1)
			number_merged=mergeAccordingToSupport(edge_1, edge_2, pair_support, gapDistance, total_edge_pair_to_merge, supportedListSearch, deleteList, highTh);
		if(flag==2)
			number_merged=mergeAccordingToSupport2(edge_1, edge_2, pair_support, gapDistance, total_edge_pair_to_merge, supportedListSearch, deleteList, highTh);
	}

	free(edge_1);
	free(edge_2);
	free(pair_support);
	free(gapDistance);
	free(listOfLongEdges);
	return number_merged;
}

uint64_t ScaffoldMaker::mergeFinalSmall(uint16_t minOverlap, int highTh, int flag)
{
	int type1a=0, type2a=0;
	uint64_t i, j, k, m, n, tos=0, number_merged=0, mate_pair_1, mate_pair_2, dist, total_edge_pair_to_merge=0, numberOfFeasibleEdges, arraySize=10000, flag_distance, index, edgeOne=0, edgeTwo=0, edgeThree=0, edgeFour=0;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> supportedListSearch;
	unordered_multimap<uint64_t,SupportedEdgePair> deleteList;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash>::const_iterator it;
	uint64_t firstEdge=0, secondEdge=0; // for calculating amount to subtract from memory usage
	int32_t *gapDistance, gapSum, support, distance_all, *distance_1, *distance_2;
	int *pair_support;
	Edge **listOfLongEdges, *edge1=NULL, *edge2=NULL, **edge_1, **edge_2, **feasibleListOfEdges;
	MatePairInfo *w;
	listOfLongEdges=getListOfCompositeEdges(&tos);
	if((edge_1=(Edge**)malloc(arraySize*sizeof(Edge*)))==NULL)
		printError(MEM_ALLOC, "edge_1");
	if((edge_2=(Edge**)malloc(arraySize*sizeof(Edge*)))==NULL)
		printError(MEM_ALLOC, "edge_2");
	if((pair_support=(int*)malloc(arraySize*sizeof(int)))==NULL)
		printError(MEM_ALLOC, "pair_support");
	if((gapDistance=(int32_t*)malloc(arraySize*sizeof(int32_t)))==NULL)
		printError(MEM_ALLOC, "gapDistance");
	for(i=0; i<tos; i++)
	{
		feasibleListOfEdges=getListOfFeasibleEdges(listOfLongEdges[i], &numberOfFeasibleEdges);
		for(j=0; j<numberOfFeasibleEdges; j++)
		{
			int64_t len1=listOfLongEdges[i]->lengthOfEdge, len2=feasibleListOfEdges[j]->lengthOfEdge;
			if(len1<500 && len2<500)
			{
				for(k=0; k<4; k++)
				{
					if(k==0)
					{
						edge1=listOfLongEdges[i];
						edge2=feasibleListOfEdges[j];
					}
					else if(k==1)
					{
						edge1=listOfLongEdges[i];
						edge2=feasibleListOfEdges[j]->twinEdge;
					}
					else if(k==2)
					{
						edge1=listOfLongEdges[i]->twinEdge;
						edge2=feasibleListOfEdges[j];
					}
					else if(k==3)
					{
						edge1=listOfLongEdges[i]->twinEdge;
						edge2=feasibleListOfEdges[j]->twinEdge;
					}
					dist=0;
					support=0;
					gapSum=0;
					if(edge1>edge2) 
						continue;
						
					for(index=1; index<=edge1->listOfReads[0].readId; index++)
					{
						dist += edge1->listOfReads[index].distPrevious;
						if(dist>mateObj->maximumUpperBoundOfInsert)
							break;
						mate_pair_1=edge1->listOfReads[index].readId;
						if(isUniqueInEdge(edge1, mate_pair_1)) // mate pair1 is only present in this edge
						{
							for(w=mateObj->matePairList[mate_pair_1]; w!=NULL; w=w->next)
							{
								mate_pair_2=w->ID;
								int lib=w->library;
								if(isUniqueInEdge(edge2, mate_pair_2)) // matepair 2 is only present in this edge
								{
									distance_1=mateObj->findDistanceOnEdge(edge1, mate_pair_1);
									distance_2=mateObj->findDistanceOnEdge(edge2, mate_pair_2);
									if(distance_1[0]>0 && distance_2[0]>0)
									{
										flag_distance=0;
										for(m=1; m<=(uint64_t)distance_1[0]; m++)
										{
											for(n=1; n<=(uint64_t)distance_2[0]; n++)
											{
												distance_all=10000000;
												if(w->type1==1)
													type1a=1;
												else
													type1a=-1;
												if(w->type2==1)
													type2a=1;
												else
													type2a=-1;

												if(type1a*distance_1[m]<0 && type2a*distance_2[n]<0)
													distance_all=llabs(distance_1[m])+llabs(distance_2[n]);
												if(distance_all+2*(uint64_t)averageReadLength<(uint64_t)(mateObj->upperBoundOfInsert[lib]+minOverlap))
												{
													flag_distance=1;
													break;
												}
											}
											if(flag_distance==1)
												break;
										}
										if(flag_distance==1)
										{
											support++;
											gapSum+=mateObj->Mean[lib]-(distance_all+2*averageReadLength);
										}
									}
									free(distance_1);
									free(distance_2);
								}
							}
						}
					}
					if(support>0)
					{
						edge_1[total_edge_pair_to_merge]=edge1->twinEdge;
						edge_2[total_edge_pair_to_merge]=edge2;
						pair_support[total_edge_pair_to_merge]=support;
						gapDistance[total_edge_pair_to_merge++]=gapSum/support;						
						SupportedEdgePair pairSupp;
						pairSupp.edge1 = edge1->twinEdge;
						pairSupp.edge2 = edge2;
						firstEdge = pairSupp.edge1->fromID ^ pairSupp.edge1->ID;
						secondEdge = pairSupp.edge2->fromID ^ pairSupp.edge2->ID;
						edgeOne = pairSupp.edge1->ID;
						edgeTwo = pairSupp.edge1->fromID;
						edgeThree = pairSupp.edge2->ID;
						edgeFour = pairSupp.edge2->fromID;

						pairSupp.support = support;
						it = supportedListSearch.find(make_pair(firstEdge, secondEdge));
						if(it != supportedListSearch.end())
						{	
							if(support > it->second.support)
								supportedListSearch[make_pair(firstEdge, secondEdge)] = pairSupp;
						}
						else
							supportedListSearch[make_pair(firstEdge, secondEdge)] = pairSupp;

						deleteList.insert({edgeOne,pairSupp});
						deleteList.insert({edgeTwo,pairSupp});
						deleteList.insert({edgeThree,pairSupp});
						deleteList.insert({edgeFour,pairSupp});	
							
						if(total_edge_pair_to_merge>=arraySize-5)
						{
							arraySize=arraySize<<1;
							if((edge_1=(Edge**)realloc(edge_1, arraySize*sizeof(Edge*)))==NULL)
								printError(MEM_ALLOC, "realloc edge_1");
							if((edge_2=(Edge**)realloc(edge_2, arraySize*sizeof(Edge*)))==NULL)
								printError(MEM_ALLOC, "realloc edge_2");
							if((pair_support=(int*)realloc(pair_support, arraySize*sizeof(int)))==NULL)
								printError(MEM_ALLOC, "realloc pair_support");
							if((gapDistance=(int32_t*)realloc(gapDistance, arraySize*sizeof(int32_t)))==NULL)
								printError(MEM_ALLOC, "realloc gapDistance");
						}
					}
				}
			}
		}
		free(feasibleListOfEdges);
	}
	
	uint64_t samllGap=0, totalGaps=0, gapDistanceTotal=0;
	int highGap=0;
	for(k=0; k<arraySize; k++)
	{
		if(gapDistance[k]<1000 && gapDistance[k]>0)
		{
			gapDistanceTotal += gapDistance[k];
			totalGaps++;
			if(highGap < gapDistance[k])
				highGap = gapDistance[k];
			if(gapDistance[k]<=10)
				samllGap++;
		}
	}
	
	if(total_edge_pair_to_merge>0)
	{
		quicksortListOfLongEdges(edge_1, edge_2, pair_support, gapDistance, 0, total_edge_pair_to_merge-1); // Sort edges according to the support.
		if(flag==1)
			number_merged=mergeAccordingToSupport(edge_1, edge_2, pair_support, gapDistance, total_edge_pair_to_merge, supportedListSearch, deleteList, highTh);
		if(flag==2)
			number_merged=mergeAccordingToSupport2(edge_1, edge_2, pair_support, gapDistance, total_edge_pair_to_merge, supportedListSearch, deleteList, highTh);
	}

	free(edge_1);
	free(edge_2);
	free(pair_support);
	free(gapDistance);
	free(listOfLongEdges);
	
	return number_merged;
}

