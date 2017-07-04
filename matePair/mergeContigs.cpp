/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 * Changes		: 
 *
 ****************************************************************************/

#include "mergeContigs.h"

extern ofstream logStream;
extern bool debugging;
extern int lowThreshold, highThreshold;
extern uint64_t genomeSize, averageReadLength;
extern string outputDir, prefixName;
ContigExtender::ContigExtender(OverlapGraph *graph1, MatePair *mate1, ReadLoader *loader1)
{
	loaderObj = loader1;
	mateObj = mate1;
	graphObj = graph1;
	contigLengthThreshold = (int)ceilf(100*(averageReadLength/100.0));
}

ContigExtender::~ContigExtender()
{
	
}


bool compareSupportPairs(const PairSupported &first, const PairSupported &second)
{
	if(first.edge1->ID < second.edge1->ID)
		return true;
	return false;
}

void ContigExtender::extendContigs(uint16_t minOverlap)
{
	time_t seconds_s=time(NULL);
	logStream<< "\nIn function mergeContigs()." << endl;
	uint64_t node_resolved=0, sum_node_resolved=0, stopValue=0;
	bool flg=true;
	float coverage = (loaderObj->numberOfReads*averageReadLength)/genomeSize;
	int highTh;
	highTh = (int)ceilf((0.20*coverage)*(100.0/averageReadLength));
	node_resolved=0;
	for(int iteration=1; flg; iteration++)
	{
		if(iteration==1 || iteration==2 || iteration==3)
			highTh = (int)ceilf((0.20*coverage)*(100.0/averageReadLength));
		else if((iteration==4 || iteration==5 || iteration==6) && genomeSize>1000000000)
			highTh = (int)ceilf((0.10*coverage)*(100.0/averageReadLength));
		else if((iteration>6) && genomeSize>1000000000)
			highTh = (int)ceilf((0.05*coverage)*(100.0/averageReadLength));

		node_resolved=findMatepairSupports(minOverlap, highTh);
		sum_node_resolved += node_resolved;
		if(node_resolved==0)
			flg=false;
	}
	
	highTh = (int)ceilf((0.20*coverage)*(100.0/averageReadLength));
	flg=true;
	if(genomeSize>1000000000)
		stopValue=10;
	
	for(int iteration=1; flg; iteration++)
	{
		if(iteration==1 || iteration==2 || iteration==3)
			highTh = (int)ceilf((0.20*coverage)*(100.0/averageReadLength));
		else if((iteration==4 || iteration==5 || iteration==6) && genomeSize>1000000000)
			highTh = (int)ceilf((0.10*coverage)*(100.0/averageReadLength));
		else if((iteration>6) && genomeSize>1000000000)
			highTh = (int)ceilf((0.05*coverage)*(100.0/averageReadLength));
		
		node_resolved=findPathSupports(highTh); // Merge pair of edges according to mate pair support
		sum_node_resolved += node_resolved;
		if(node_resolved<=stopValue)
			flg=false;
	}
	
	flg=true;
	if(genomeSize>1000000000)
		stopValue=10;
	
	for(int iteration=1; flg; iteration++)
	{
		if(iteration==1 || iteration==2 || iteration==3)
			highTh = (int)ceilf((0.20*coverage)*(100.0/averageReadLength));
		else if((iteration==4 || iteration==5 || iteration==6) && genomeSize>1000000000)
			highTh = (int)ceilf((0.10*coverage)*(100.0/averageReadLength));
		else if((iteration>6) && genomeSize>1000000000)
			highTh = (int)ceilf((0.05*coverage)*(100.0/averageReadLength));
		
		node_resolved=findPathSupportsNew(highTh); // Merge pair of edges according to mate pair support
		sum_node_resolved += node_resolved;
		if(node_resolved<=stopValue)
			flg=false;
	}

	node_resolved=0;
	highTh = (int)ceilf((0.05*coverage)*(100.0/averageReadLength));
	flg=true;
	for(int iteration=1; flg; iteration++)
	{
		node_resolved=mergeShortOverlaps(minOverlap, highTh);	
		sum_node_resolved += node_resolved;
		if(node_resolved==0)
			flg=false;;
	}
	logStream<< "\nTotal number of extended contigs: " << sum_node_resolved << "\n"<<flush;
	logStream<< "Function mergeContigs() in " << time(NULL)-seconds_s << " sec.\n"<<flush;
}

/* ============================================================================================
   This function will find paths from each node in the graph, one by one
   ============================================================================================ */
uint64_t ContigExtender::findPathSupports(int highTh)
{
	uint64_t i, returnValue=0;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> supportedListSearch;
	unordered_multimap<uint64_t,SupportedEdgePair> deleteList;
	bool *flag;//this is the flag whether we already found paths of the mate pairs or not.
	if((flag=(bool*)malloc((loaderObj->numberOfUniqueReads+1)*sizeof(bool)))==NULL)
		printError(MEM_ALLOC,"flag");

	for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
		flag[i]=1;

	#pragma omp parallel
	{		
		uint64_t mate_pair_1, mate_pair_2, total_supported_pairs, j, l, index, firstEdge=0, secondEdge=0, edgeOne=0, edgeTwo=0, edgeThree=0, edgeFour=0, distance=0, no_of_reads=0, last=0, arraySize = 20000, totalPaths = 0, readsHigh=0, pairSuppHigh=0;

		unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> supportedListSearchThread;
		unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash>::const_iterator it;
		unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash>::const_iterator it1;
		unordered_map<uint64_t,int> list_of_explored;
		unordered_map<uint64_t,int>::const_iterator it2;

		MatePairInfo *w;
		Edge *v;
		ReadToEdgeMap *a;
		uint64_t *list_of_reads;		
		PairSupported *supportedPairs_tmp;
		Edge **pathEdges, ***allPairPaths;
		uint64_t *pathLengths, *allPairPathLengths;
		Indexes **indx;
		if((list_of_reads=(uint64_t*)malloc(arraySize*sizeof(uint64_t)))==NULL)
			printError(MEM_ALLOC,"list_of_reads");
		if((pathLengths=(uint64_t*)malloc((20)*sizeof(uint64_t)))==NULL)
			printError(MEM_ALLOC, "pathLengths");
		if((pathEdges=(Edge**)malloc((20)*sizeof(Edge*)))==NULL)
			printError(MEM_ALLOC, "pathEdges");
		if((allPairPaths=(Edge***)malloc(2000*sizeof(Edge**)))==NULL)
			printError(MEM_ALLOC,"allPairPaths");
		if((allPairPathLengths=(uint64_t*)malloc(2000*sizeof(uint64_t)))==NULL)
			printError(MEM_ALLOC,"allPairPathLengths");
		for(j=0;j<2000;j++)
		{	
			if((allPairPaths[j]=(Edge**)malloc((10)*sizeof(Edge*)))==NULL)
				printError(MEM_ALLOC, "allPairPaths[totalPaths]");
		}		

		#pragma omp for
		for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
		{
			
			if(graphObj->graph[i]!=NULL) //for each node in the graph
			{
				for(j=0;j<arraySize;j++)
					list_of_reads[j]=0;
				for(j=0;j<20;j++)
				{	
					pathLengths[j]=0;
					pathEdges[j]=NULL;
				}
				for(j=0;j<2000;j++)
				{
					allPairPathLengths[j]=0;
					for(l=0;l<10;l++)
					{
						allPairPaths[j][l]=NULL;
					}
				}

				totalPaths=0;
				last=0;
				no_of_reads=0;

				for(v=graphObj->graph[i]; v!=NULL; v=v->next)
				{
					distance=0;
					for(index=1; index<=v->listOfReads[0].readId; index++)
					{
						distance += v->listOfReads[index].distPrevious;
						if(distance>mateObj->maximumUpperBoundOfInsert)
							break;
						if(flag[v->listOfReads[index].readId]==1)
						{
							list_of_reads[no_of_reads++]=v->listOfReads[index].readId;
							flag[v->listOfReads[index].readId]=0;
						}
					}
					distance=0;
					for(index=1; index<=v->twinEdge->listOfReads[0].readId; index++)
					{
						distance += v->twinEdge->listOfReads[index].distPrevious;
						if(distance>mateObj->maximumUpperBoundOfInsert)
							break;
						if(flag[v->twinEdge->listOfReads[index].readId]==1)
						{
							list_of_reads[no_of_reads++]=v->twinEdge->listOfReads[index].readId;
							flag[v->twinEdge->listOfReads[index].readId]=0;
						}
					}
				}
				if(readsHigh<no_of_reads)
						readsHigh=no_of_reads;
					
				for(l=0; l<no_of_reads; l++)
				{
					for(a=mateObj->readToEdgeList[list_of_reads[l]]; a!=NULL; a=a->next)
					{
						it2 = list_of_explored.find((uint64_t)a->edge->fromID);
						if(it2 == list_of_explored.end())
						{
							list_of_explored[(uint64_t)a->edge->fromID] = 1;
							last++;
							exploreGraph(a->edge->fromID, 1, totalPaths, pathEdges, pathLengths, allPairPaths, allPairPathLengths);
						}
						it2 = list_of_explored.find((uint64_t)a->edge->ID);
						if(it2 == list_of_explored.end())
						{
							list_of_explored[(uint64_t)a->edge->ID] = 1;
							last++;
							exploreGraph(a->edge->ID, 1, totalPaths, pathEdges, pathLengths, allPairPaths, allPairPathLengths);
						}
					}
				}
				list_of_explored.clear();
				
				if(totalPaths==0 || last==0 || no_of_reads==0)
					continue;

				quicksortPaths(0, totalPaths-1, allPairPaths, allPairPathLengths);
				if((indx = (Indexes**)malloc((last+2)*sizeof(Indexes*)))==NULL)
					printError(MEM_ALLOC, "indx");

				indexPaths(last, indx, totalPaths, allPairPaths);
				
				if((supportedPairs_tmp=(PairSupported*)malloc(arraySize*sizeof(PairSupported)))==NULL)
					printError(MEM_ALLOC, "supportedPairs_tmp");

				for(l=0; l<no_of_reads; l++)
				{
					mate_pair_1=list_of_reads[l];
					for(w=mateObj->matePairList[mate_pair_1]; w!=NULL; w=w->next)
					{
						mate_pair_2=w->ID;
						if(w->flag==1) 
						{
							if(graphObj->graph[mate_pair_1]==NULL && graphObj->graph[mate_pair_2]==NULL && flag[mate_pair_2]==1) 
							{
								if(totalPaths==0)
									total_supported_pairs=0;
								else
									total_supported_pairs=findPathsBetweenMatepairs(mate_pair_1, mate_pair_2, w->type1, w->type2, w->library, indx, totalPaths, allPairPaths, allPairPathLengths, supportedPairs_tmp);

								if(pairSuppHigh<total_supported_pairs)
									pairSuppHigh=total_supported_pairs;
								
								if(total_supported_pairs>0)
								{
									for(j=0; j<total_supported_pairs; j++)
									{
										SupportedEdgePair pairSupp;
										pairSupp.edge1 = supportedPairs_tmp[j].edge1;
										pairSupp.edge2 = supportedPairs_tmp[j].edge2;
										edgeOne = pairSupp.edge1->ID;
										edgeTwo = pairSupp.edge1->fromID;
										edgeThree = pairSupp.edge2->ID;
										edgeFour = pairSupp.edge2->fromID;
										firstEdge = supportedPairs_tmp[j].edge1->fromID ^ supportedPairs_tmp[j].edge1->ID;
										secondEdge = supportedPairs_tmp[j].edge2->fromID ^ supportedPairs_tmp[j].edge2->ID;
										pairSupp.support = 1;
										it = supportedListSearchThread.find(make_pair(firstEdge, secondEdge));
										if(it != supportedListSearchThread.end())
										{
											pairSupp.support = it->second.support+1;
											supportedListSearchThread.erase(make_pair(firstEdge, secondEdge));
											supportedListSearchThread[make_pair(firstEdge, secondEdge)] = pairSupp;
										}
										else
											supportedListSearchThread[make_pair(firstEdge, secondEdge)] = pairSupp;
									}
								}
							}
						}
					}
				}
				freeIndexedPaths(indx);
				free(indx);
				free(supportedPairs_tmp);
			}
		}
		
		#pragma omp critical (addToList)
		{
			for(it = supportedListSearchThread.begin(); it != supportedListSearchThread.end(); ++it)
			{
				SupportedEdgePair pairSupp;
				pairSupp.edge1 = it->second.edge1;
				pairSupp.edge2 = it->second.edge2;
				pairSupp.support = it->second.support;
				firstEdge = it->second.edge1->fromID ^ it->second.edge1->ID;
				secondEdge = it->second.edge2->fromID ^ it->second.edge2->ID;
				it1 = supportedListSearch.find(make_pair(firstEdge, secondEdge));
				// insert the pair into the unordered_map
				if(it1 != supportedListSearch.end())
				{
					pairSupp.support = it1->second.support + it->second.support;
					supportedListSearch.erase(make_pair(firstEdge, secondEdge));
					supportedListSearch[make_pair(firstEdge, secondEdge)] = pairSupp;
				}
				else
					supportedListSearch[make_pair(firstEdge, secondEdge)] = pairSupp;

				edgeOne = pairSupp.edge1->ID;
				edgeTwo = pairSupp.edge1->fromID;
				edgeThree = pairSupp.edge2->ID;
				edgeFour = pairSupp.edge2->fromID;
				// insert the edges into the unordered_multimap
				deleteList.insert({edgeOne,pairSupp});
				deleteList.insert({edgeTwo,pairSupp});
				deleteList.insert({edgeThree,pairSupp});
				deleteList.insert({edgeFour,pairSupp});
			}
			supportedListSearchThread.clear();
			for(j=0;j<2000;j++)
				free(allPairPaths[j]);
			free(allPairPaths);
			free(allPairPathLengths);
			free(pathLengths);
			free(pathEdges);
			free(list_of_reads);
		}
		
	}
	free(flag);		
	returnValue=mergeByPathSupports(supportedListSearch, deleteList, highTh);
	
	return returnValue;
}

/* ============================================================================================
   This function will find all pair of feasible path in the graph starting from node.
   Backtracking used by recursion.
   ============================================================================================ */
uint64_t ContigExtender::exploreGraph(uint64_t node, int level, uint64_t &totalPaths, Edge **pathEdges, uint64_t *pathLengths, Edge ***allPairPaths, uint64_t *allPairPathLengths)
{
	Edge *v;
	if(level>7 || totalPaths>1000)
		return 0;
	for(v=graphObj->graph[node]; v!=NULL; v=v->next)
	{
		if(v->isReducible>=3) 
			continue;
		if(level==1)
		{
			if(v->listOfReads[0].readId!=0 && v->flow>0) 
			{
				pathEdges[level]=v;
				pathLengths[level]=0;
				savePath(pathEdges, pathLengths[level], level, totalPaths, allPairPaths, allPairPathLengths);
				exploreGraph(v->ID, level+1, totalPaths, pathEdges, pathLengths, allPairPaths, allPairPathLengths);
			}
		}
		else if(matchEdgeType(pathEdges[level-1], v) && checkFlow(pathEdges,level-1,v) && (pathLengths[level-1]+pathEdges[level-1]->lengthOfEdge)<mateObj->maximumUpperBoundOfInsert) 
		{
			pathEdges[level]=v;
			pathLengths[level]=pathLengths[level-1]+pathEdges[level-1]->lengthOfEdge;
			savePath(pathEdges, pathLengths[level], level, totalPaths, allPairPaths, allPairPathLengths);
			exploreGraph(v->ID,level+1, totalPaths, pathEdges, pathLengths, allPairPaths, allPairPathLengths);
		}
	}
	return 1;
}

/* ============================================================================================
   This function check if edge1 and edge2 have proper orientation.
   ============================================================================================ */
int ContigExtender::matchEdgeType(Edge *edge1, Edge *edge2)
{
	if(((edge1->typeOfEdge==1 || edge1->typeOfEdge==3) && (edge2->typeOfEdge==2 || edge2->typeOfEdge==3))||((edge1->typeOfEdge==0 || edge1->typeOfEdge==2) && (edge2->typeOfEdge==0 || edge2->typeOfEdge==1)))
		return 1;
	return 0;
}

/* ============================================================================================
   This function checks if inserting next edge in the path violates the flow of next edge or not
   Here we added allow the edge to be used flow+2 times (as the flow in not always accurate.) 
   ============================================================================================ */
int ContigExtender::checkFlow(Edge **edges, uint64_t top, Edge *next_edge)
{
	uint64_t i;
	float flow=0;
	for(i=1;i<=top;i++)
		if(edges[i]==next_edge || edges[i]==next_edge->twinEdge)
			flow+=1;
	if(next_edge->listOfReads[0].readId==0 && flow<=next_edge->flow) // 1 extra units of flow for simple edges
		return 1;
	if(flow<next_edge->flow)
		return 1;
	return 0;
}

/* ============================================================================================
   This function will save the path found in the main memory.
   ============================================================================================ */
void ContigExtender::savePath(Edge **edges, uint64_t distance, uint64_t top, uint64_t &totalPaths, Edge ***allPairPaths, uint64_t *allPairPathLengths)
{
	uint64_t i;
	for(i=1;i<=top;i++)
		allPairPaths[totalPaths][i]=edges[i];
	allPairPaths[totalPaths][0]=edges[top];
	allPairPaths[totalPaths][top+1]=NULL;
	allPairPathLengths[totalPaths]=distance;
	totalPaths++;
}

/* ============================================================================================
   quicksortPaths() for sorting paths algorithm
   ============================================================================================ */
void ContigExtender::quicksortPaths(int64_t left, int64_t right, Edge ***allPairPaths, uint64_t *allPairPathLengths)
{
	uint64_t j;
	if( left < right )
	{
		j = partitionPaths(left,right, allPairPaths, allPairPathLengths);
		if(j>0)
			quicksortPaths(left,j-1, allPairPaths, allPairPathLengths);
		quicksortPaths(j+1,right, allPairPaths, allPairPathLengths);
	}
}


/* ============================================================================================
   This function is used to partition the array in quicksortPaths() algorithm
   ============================================================================================ */
uint64_t ContigExtender::partitionPaths(int64_t left, int64_t right, Edge ***allPairPaths, uint64_t *allPairPathLengths)
{
	int64_t i=left, j=right+1, temp_dist;
	Edge **temp_location;
	while(1)
	{
		do
			++i;
		while(i<=right && (allPairPaths[i][1]->fromID<allPairPaths[left][1]->fromID || (allPairPaths[i][1]->fromID==allPairPaths[left][1]->fromID && allPairPaths[i][0]->fromID<allPairPaths[left][0]->fromID)));
		do
			--j;
		while(allPairPaths[j][1]->fromID>allPairPaths[left][1]->fromID || (allPairPaths[j][1]->fromID==allPairPaths[left][1]->fromID && allPairPaths[j][0]->fromID>allPairPaths[left][0]->fromID));
		if( i >= j )
			break;
		temp_location=allPairPaths[i];
		allPairPaths[i]=allPairPaths[j];
		allPairPaths[j]=temp_location;
		temp_dist=allPairPathLengths[i];
		allPairPathLengths[i]=allPairPathLengths[j];
		allPairPathLengths[j]=temp_dist;
	}
	temp_location=allPairPaths[left];
	allPairPaths[left]=allPairPaths[j];
	allPairPaths[j]=temp_location;
	temp_dist=allPairPathLengths[left];
	allPairPathLengths[left]=allPairPathLengths[j];
	allPairPathLengths[j]=temp_dist;
  	return j;
}

/* ============================================================================================
   Index the string node of the paths for quick access.
   ============================================================================================ */
void ContigExtender::indexPaths(uint64_t num, Indexes **indx, uint64_t &totalPaths, Edge ***allPairPaths)
{
	uint64_t i, firstReadOnEdge=-1, second_read=-1, count=-1;
	for(i=0;i<num+2;i++)
		indx[i]=NULL;
	for(i=0; i<totalPaths; i++)
	{
		if(allPairPaths[i][1]->fromID!=firstReadOnEdge)
		{
			Indexes *v;
			if((v = (Indexes*)malloc(sizeof(Indexes)))==NULL)
				printError(MEM_ALLOC, "v");
			v->firstNode=allPairPaths[i][1]->fromID;
			v->lastNode=allPairPaths[i][0]->fromID;
			v->position=i;
			v->next=NULL;
			firstReadOnEdge=allPairPaths[i][1]->fromID;
			second_read=allPairPaths[i][0]->fromID;
			count++;
			if(indx[count]==NULL)
			{
				indx[count]=v;
			}
			else
			{
				v->next=indx[count];
				indx[count]=v;
			}
		}
		else if(allPairPaths[i][0]->fromID!=second_read)
		{
			Indexes *vv;
			if((vv = (Indexes*)malloc(sizeof(Indexes)))==NULL)
				printError(MEM_ALLOC, "vv");
			vv->firstNode=allPairPaths[i][1]->fromID;
			vv->lastNode=allPairPaths[i][0]->fromID;
			vv->next=NULL;
			vv->position=i;
			second_read=allPairPaths[i][0]->fromID;
			if(indx[count]==NULL)
			{
				indx[count]=vv;
			}
			else
			{
				vv->next=indx[count];
				indx[count]=vv;
			}
		}
	}
}

/* ============================================================================================
   This function will find all paths between mate_pair_1 and mate_pair_2 and will find
   if the all goes through (a,node), (node,b) for some a and b.
   ============================================================================================ */
int ContigExtender::findPathsBetweenMatepairs(uint64_t mate_pair_1, uint64_t mate_pair_2, bool type1, bool type2, int library, Indexes **indx, uint64_t &totalPaths, Edge ***allPairPaths, uint64_t *allPairPathLengths, PairSupported *supportedPairs_tmp)
{
	Edge **edges, *edge1=NULL, *edge2=NULL;
	ReadToEdgeMap *edge_ptr1, *edge_ptr2;
	uint64_t i, k, l, m, path_found=0, total_supported=0, x;
	int32_t *dist_on_edge1, *dist_on_edge2;
	int type1a=type1, type2a=type2;
	if(type1a==0)
		type1a=-1;
	if(type2a==0)
		type2a=-1;
	
	int *support_flag;
	if(mateObj->readToEdgeList[mate_pair_1]==NULL || mateObj->readToEdgeList[mate_pair_2]==NULL)  // not in any edge
		return 0;
	
	if((edges=(Edge**)malloc(1000*sizeof(Edge*)))==NULL)
		printError(MEM_ALLOC, "edges");
	if((support_flag=(int*)malloc(1000*sizeof(int)))==NULL)
		printError(MEM_ALLOC, "support_flag");
	for(edge_ptr1=mateObj->readToEdgeList[mate_pair_1]; edge_ptr1!=NULL; edge_ptr1=edge_ptr1->next)// find an edge containing mate_pair_1
	{
		for(edge_ptr2=mateObj->readToEdgeList[mate_pair_2]; edge_ptr2!=NULL; edge_ptr2=edge_ptr2->next)// find an edge containing mate_pair_2
		{
			if(isUniqueInEdge(edge_ptr1->edge, mate_pair_1) && isUniqueInEdge(edge_ptr2->edge, mate_pair_2))
			{
				for(x=0; x<4; x++)
				{
					if(x==0)
					{
						edge1=edge_ptr1->edge;
						edge2=edge_ptr2->edge;
					}
					else if(x==1)
					{
						edge1=edge_ptr1->edge;
						edge2=edge_ptr2->edge->twinEdge;
					}
					else if(x==2)
					{
						edge1=edge_ptr1->edge->twinEdge;
						edge2=edge_ptr2->edge;
					}
					else if(x==3)
					{
						edge1=edge_ptr1->edge->twinEdge;
						edge2=edge_ptr2->edge->twinEdge;
					}
					dist_on_edge1=mateObj->findDistanceOnEdge(edge1->twinEdge, mate_pair_1);
					dist_on_edge2=mateObj->findDistanceOnEdge(edge2, mate_pair_2);
					if(dist_on_edge1[0]==0 || dist_on_edge2[0]==0)
					{
						free(dist_on_edge1);
						free(dist_on_edge2);
						continue;
					}
					int strt;
					strt=startingPosition(edge1->ID, edge2->fromID, indx);
					if(strt==-1) //not path possible
					{
						free(dist_on_edge1);
						free(dist_on_edge2);
						continue;
					}
					for(i=strt;i<totalPaths;i++)
					{
						if(allPairPaths[i][1]->fromID>edge1->ID || (allPairPaths[i][1]->fromID==edge1->ID && allPairPaths[i][0]->fromID>edge2->fromID))
							break;
						if(edge1->ID==allPairPaths[i][1]->fromID && matchEdgeType(edge1,allPairPaths[i][1]) && allPairPaths[i][0]==edge2)
						{
							int flg=0;
							for(k=1; k<=(uint32_t)dist_on_edge1[0]; k++)
							{
								for(l=1;l<=(uint64_t)dist_on_edge2[0];l++)
								{
									int64_t distance=100000;
									if(type1a*(int)dist_on_edge1[k]<0 && type2a*(int)dist_on_edge2[l]<0)
										distance=llabs(dist_on_edge1[k])+llabs(dist_on_edge2[l])+allPairPathLengths[i];
									if(distance<mateObj->upperBoundOfInsert[library] && distance>mateObj->lowerBoundOfInsert[library])
									{
										flg=1;
										break;
									}
								}
								if(flg==1)
								break;
							}
							if(flg==1)
							{
								edges[0]=edge1;
								for(m=1; allPairPaths[i][m]!=NULL; m++)
									edges[m]=allPairPaths[i][m];
								if(insertPath(edges, support_flag, m-1, path_found, supportedPairs_tmp))
								{
									path_found++;
								}
								else
								{
									free(dist_on_edge1);
									free(dist_on_edge2);
									free(edges);
									free(support_flag);
									return 0;
								}
							}
						}
					}
					free(dist_on_edge1);
					free(dist_on_edge2);
				}
			}
		}
	}
	if(path_found>0)
	{
		for(i=0; support_flag[i]>=0; i++) // move the supported pairs (flag=1) at the beginning of the stack
		{
			if(support_flag[i]==1)
			{
				supportedPairs_tmp[total_supported]=supportedPairs_tmp[i];
				total_supported++;
			}
		}
	}
	else
	{
		total_supported=0;
	}
	free(edges);
	free(support_flag);
	return total_supported;
}

/* ============================================================================================
   This function will return that appropriate location in the array.
   ============================================================================================ */
uint64_t ContigExtender::startingPosition(uint64_t firstReadOnEdge,uint64_t second_read, Indexes **indx)
{
	uint64_t i;
	Indexes *v;
	for(i=0; indx[i]!=NULL; i++)
	{
		if(indx[i]->firstNode==firstReadOnEdge)
		{
			for(v=indx[i];v!=NULL;v=v->next)
			{
				if(v->lastNode==second_read)
				{
					return v->position;
				}
			}
		}
	}
	return -1;
}

/* ============================================================================================
   This function insert the path in the support pair list. if no pair of edge is supported it
   returns 0. Otherwise returns 1.
   ============================================================================================ */
uint64_t ContigExtender::insertPath(Edge **edges, int *support_flag, uint64_t top, uint64_t path_found, PairSupported *supportedPairs_tmp)
{
	uint64_t k,l;
	if(path_found==0) // if it is the first path found then store pair of edges (a,b) and (b,c) and mark them as supported
	{
		for(k=0; k<top+10; k++)
			support_flag[k]=-1;
		for(k=0; k<top; k++)
		{
			supportedPairs_tmp[k].edge1=edges[k];
			supportedPairs_tmp[k].edge2=edges[k+1];
			support_flag[k]=1;
		}
	}
	else // if not the first path
	{
		for(k=0; support_flag[k]>=0; k++) // unmark the pair of edges (a,b) and (b,c) not present in the new path
		{
			for(l=0;l<top;l++)
			{
				if((supportedPairs_tmp[k].edge1==edges[l] && supportedPairs_tmp[k].edge2==edges[l+1])||(supportedPairs_tmp[k].edge1==edges[l+1]->twinEdge && supportedPairs_tmp[k].edge2==edges[l]->twinEdge))
				{
					break;
				}
			}
			if(l==top) // not supported in new path
			{
				support_flag[k]=0;
			}
		}

	}
	for(k=0; support_flag[k]>=0; k++)
	{
		if(support_flag[k]==1)
		{
			break;
		}
	}
	if(support_flag[k]==-1) //no pair of edges supported. no need to continue more
	{
		return 0;
	}
	return 1;
}

/* ============================================================================================
   This function will free memory used by path index
   ============================================================================================ */
void ContigExtender::freeIndexedPaths(Indexes **indx)
{
	uint64_t i;
	Indexes *v, *v_next;
	for(i=0; indx[i]!=NULL; i++)
	{
		for(v=indx[i]; v!=NULL; v=v_next)
		{
			v_next=v->next;
			free(v);
		}
	}
}

/* ============================================================================================
   Sort list of long edges according to support. If the supports are same then sort by length.
   ============================================================================================ */
void ContigExtender::quicksortListOfLongEdges(PairSupported *pairs, int *pair_support, int64_t left, int64_t right)
{
	int64_t i=left, j=right, temporary;
	int64_t pivot=*(pair_support+((left+right)>>1));
	int64_t lengthSum=pairs[(left+right)>>1].edge1->lengthOfEdge+pairs[(left+right)>>1].edge2->lengthOfEdge;
	Edge *temp;
	while(i<j)
	{
		while((int64_t)*(pair_support+i)>pivot || ((int64_t)*(pair_support+i)==pivot && ((int64_t)(pairs[i].edge1->lengthOfEdge+pairs[i].edge2->lengthOfEdge))>lengthSum))
			i++;
		while((int64_t)*(pair_support+j)<pivot || ((int64_t)*(pair_support+j)==pivot && ((int64_t)(pairs[j].edge1->lengthOfEdge+pairs[j].edge2->lengthOfEdge))<lengthSum))
			j--;
		if (i<=j)
		{
			temp=(pairs+i)->edge1; (pairs+i)->edge1=(pairs+j)->edge1; (pairs+j)->edge1=temp;
			temp=(pairs+i)->edge2; (pairs+i)->edge2=(pairs+j)->edge2; (pairs+j)->edge2=temp;
			temporary=*(pair_support+i); *(pair_support+i)=*(pair_support+j); *(pair_support+j)=temporary;
			i++; j--;
		}
	}
	if (left < j )
		quicksortListOfLongEdges(pairs, pair_support, left, j);
    if (i < right)
		quicksortListOfLongEdges(pairs, pair_support, i, right);
}

/* ============================================================================================
   Merge pair of edges according to supports.
   ============================================================================================ */
uint64_t ContigExtender::mergeByPathSupports(unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> &listSearch, unordered_multimap<uint64_t,SupportedEdgePair> &deleteList, int highTh)
{
	uint64_t j, node_resolved=0,firstEdge=0,secondEdge=0;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash>::const_iterator it;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash>::const_iterator it1, next_it;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash>::const_iterator it3;
	unordered_map<uint64_t,int> nodeSearch;
	unordered_map<uint64_t,int>::const_iterator it2;
	Edge *v;
	deque<Edge*> inEdges;
	deque<Edge*> outEdges;
	int *inSupp, *outSupp;
	bool uniqueIn, uniqueOut, flag=false;
	uint64_t nodeToResolve;
	for(it1 = listSearch.begin(); it1 != listSearch.end(); ++it1)
	{
		if (flag)
		{
			if(next_it != listSearch.end())
			{   
				it1 = next_it;
				flag=false;
			}else
				break;
		}

		firstEdge = it1->second.edge1->fromID ^ it1->second.edge1->ID;
		secondEdge = it1->second.edge2->fromID ^ it1->second.edge2->ID;
		it3 = listSearch.find(make_pair(firstEdge, secondEdge));
		if(it3 == listSearch.end())
				continue;

		nodeToResolve=it1->second.edge1->ID;
		it2 = nodeSearch.find(nodeToResolve);
		if(it2 != nodeSearch.end())
				continue;

		if(it1->second.support >= highTh)
		{
				inEdges.clear();
				outEdges.clear();
				for(v=graphObj->graph[nodeToResolve]; v!=NULL; v=v->next)
				{
						if(v->flow>=1 && v->ID!=v->fromID)
						{
								if(v->twinEdge!=it1->second.edge1 && matchEdgeType(v->twinEdge, it1->second.edge2))
								{
										inEdges.push_back(v->twinEdge);
								}
								if(v!=it1->second.edge2 && matchEdgeType(it1->second.edge1, v))
								{
										outEdges.push_back(v);
								}
						}
				}
				inSupp=NULL;
				outSupp=NULL;
			   // calculate supports for every pair of in-out edges
			   if(inEdges.size()!=0)
				{
						if((inSupp=(int*)malloc(inEdges.size()*sizeof(int)))==NULL)
								printError(MEM_ALLOC, "inSupp");
						for(j=0; j<inEdges.size(); j++)
						{
								inSupp[j]=0;
								firstEdge=inEdges[j]->fromID ^ inEdges[j]->ID;
								secondEdge=it1->second.edge2->fromID ^ it1->second.edge2->ID;
								it = listSearch.find(make_pair(firstEdge, secondEdge));
								if(it != listSearch.end())
										inSupp[j]=it->second.support;
								else
								{
										firstEdge=it1->second.edge2->twinEdge->fromID ^ it1->second.edge2->twinEdge->ID;
										secondEdge=inEdges[j]->twinEdge->fromID ^ inEdges[j]->twinEdge->ID;
										it = listSearch.find(make_pair(firstEdge, secondEdge));
										if(it != listSearch.end())
												inSupp[j]=it->second.support;
								}
						}
				}
				if(outEdges.size()!=0)
				{
						if((outSupp=(int*)malloc(outEdges.size()*sizeof(int)))==NULL)
								printError(MEM_ALLOC, "outSupp");
						for(j=0; j<outEdges.size(); j++)
						{
								outSupp[j]=0;
								firstEdge=it1->second.edge1->fromID ^ it1->second.edge1->ID;
								secondEdge=outEdges[j]->fromID ^ outEdges[j]->ID;
								it = listSearch.find(make_pair(firstEdge, secondEdge));
								if(it != listSearch.end())
										outSupp[j]=it->second.support;
								else
								{
										firstEdge=outEdges[j]->twinEdge->fromID ^ outEdges[j]->twinEdge->ID;
										secondEdge=it1->second.edge1->twinEdge->fromID ^ it1->second.edge1->twinEdge->ID;
										it = listSearch.find(make_pair(firstEdge, secondEdge));
										if(it != listSearch.end())
												outSupp[j]=it->second.support;
								}
						}
				}
				uniqueIn=true;
				uniqueOut=true;
				for(j=0; j<outEdges.size(); j++)
				{
						if(outSupp[j]>lowThreshold)
						{       uniqueOut=false;
						}
				}
				for(j=0; j<inEdges.size(); j++)
				{
						if(inSupp[j]>lowThreshold)
						{       uniqueIn=false;
						}
				}
				if(uniqueIn && uniqueOut)
				{
						for(v=graphObj->graph[nodeToResolve]; v!=NULL; v=v->next)
						{
								it2 = nodeSearch.find(v->ID);
								if(it2 == nodeSearch.end() && v->ID!=v->fromID)
										nodeSearch[v->ID]=1;
						}
						//set next_it
						next_it=it1;
						++next_it;
						while(next_it != listSearch.end())
						{
							if (next_it->second.edge1->fromID == it1->second.edge1->ID ||
									next_it->second.edge1->ID == it1->second.edge1->ID ||
									next_it->second.edge2->fromID == it1->second.edge1->ID ||
									next_it->second.edge2->ID == it1->second.edge1->ID)
							{
								++next_it;
							}else
								break;
						}

						if(mergeEdgesAndUpdate(it1->second.edge1, it1->second.edge2, 1))
						{
								flag=true;
								node_resolved++;
								auto its = deleteList.equal_range(nodeToResolve);
								for (auto it = its.first; it != its.second; ++it)
								{
										firstEdge=it->second.edge1->fromID ^ it->second.edge1->ID;
										secondEdge=it->second.edge2->fromID ^ it->second.edge2->ID;
										listSearch.erase(make_pair(firstEdge, secondEdge));

								}
						}
				}
				free(inSupp);
				free(outSupp);
		}
	}
	return node_resolved;
}

bool compEdgePairWithSupp(SupportedEdgePair first, SupportedEdgePair second)
{
	return (first.support>second.support);
}

uint64_t ContigExtender::findMatepairSupports(uint16_t minOverlap, int highTh)
{
	uint64_t i, j, k, merged_number=0, firstEdge=0, secondEdge=0, edgeOne=0, edgeTwo=0, edgeThree=0, edgeFour=0;
	int support;
	deque<SupportedEdgePair> supportedList;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> supportedListSearch;
	unordered_multimap<uint64_t,SupportedEdgePair> deleteList;
	Edge *v;
	deque<Edge*> inEdges;
	deque<Edge*> outEdges;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash>::const_iterator it;
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
	{
		if(graphObj->graph[i]!=NULL) //for each node in the graph
		{
			inEdges.clear();
			outEdges.clear();
			// find in edges and out edges
			for(v=graphObj->graph[i]; v!=NULL; v=v->next)
			{
				if(v->flow>=1 && v->listOfReads[0].readId>0 && v->ID!=v->fromID) //composite edges with flow
				{
					if(v->typeOfEdge==0 || v->typeOfEdge==1)
						inEdges.push_back(v->twinEdge);
					else
						outEdges.push_back(v);
				}
			}
			for(j=0; j<inEdges.size(); j++)
			{
				for(k=0; k<outEdges.size(); k++)
				{
					support = getSupportOfEdges(inEdges[j]->twinEdge, outEdges[k], minOverlap);
					if(support>1)
					{
						SupportedEdgePair pairSupp;
						pairSupp.edge1 = inEdges[j];
						pairSupp.edge2 = outEdges[k];
						edgeOne = pairSupp.edge1->ID;
						edgeTwo = pairSupp.edge1->fromID;
						edgeThree = pairSupp.edge2->ID;
						edgeFour = pairSupp.edge2->fromID;
						firstEdge = inEdges[j]->fromID ^ inEdges[j]->ID;
						secondEdge = outEdges[k]->fromID ^ outEdges[k]->ID;
						pairSupp.support = support;
						supportedList.push_back(pairSupp);
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
					}
				}
			}
		}
	}
	if(supportedList.size()>0)
	{
		sort(supportedList.begin(), supportedList.end(), compEdgePairWithSupp);
		merged_number = mergeByMatepairSupports(supportedList,supportedListSearch,deleteList, highTh);
	}
		
	return merged_number;
}

int ContigExtender::getSupportOfEdges(Edge *edge1, Edge *edge2, uint16_t minOverlap)
{
	uint64_t index, dist, mate_pair_1, mate_pair_2, m, n;
	int support, flag_distance, type1a=0, type2a=0;
	MatePairInfo *w;
	int32_t *distance_1, *distance_2, distance_all;
	dist = 0;
	support = 0;
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
						}
					}
					free(distance_1);
					free(distance_2);
				}
			}
		}
	}
	return support;
}

/* ============================================================================================
   This function checks if read is present in edge.
   ============================================================================================ */
int ContigExtender::isUniqueInEdge(Edge *edge,uint64_t read)
{
	if(mateObj->readToEdgeList[read]==NULL) //not present;
		return 0;

	if((mateObj->readToEdgeList[read]->edge==edge || mateObj->readToEdgeList[read]->edge->twinEdge==edge) && mateObj->readToEdgeList[read]->next==NULL && graphObj->graph[read]==NULL)
		return 1;
	return 0;
}

uint64_t ContigExtender::mergeByMatepairSupports(deque<SupportedEdgePair> &list, unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> &listSearch, unordered_multimap<uint64_t,SupportedEdgePair> &deleteList, int highTh)
{
	uint64_t merged=0, i, j, firstEdge=0, secondEdge=0, edgeID=0;
	Edge *v;
	deque<Edge*> inEdges;
	deque<Edge*> outEdges;
	int *inSupp, *outSupp;
	bool uniqueIn, uniqueOut;
	uint64_t nodeToResolve;
	unordered_map<uint64_t,int> nodeSearch;
	unordered_map<uint64_t,int>::const_iterator it2;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash>::const_iterator it;

	for(i=0; i<list.size(); i++)
	{
		if(list[i].edge1==NULL || list[i].edge2==NULL)
			continue;
		firstEdge=list[i].edge1->fromID ^ list[i].edge1->ID;
		secondEdge=list[i].edge2->fromID ^ list[i].edge2->ID;
		it = listSearch.find(make_pair(firstEdge, secondEdge));
		if(it == listSearch.end())
			continue;
		nodeToResolve = list[i].edge1->ID;
		it2 = nodeSearch.find(nodeToResolve);
		if(it2 != nodeSearch.end())
			continue;
		if(list[i].support >= highTh)
		{
			inEdges.clear();
			outEdges.clear();
			// find in edges and out edges
			for(v=graphObj->graph[nodeToResolve]; v!=NULL; v=v->next)
			{
				if(v->flow>=1 && v->listOfReads[0].readId>0 && v->ID!=v->fromID)
				{
					if(v->twinEdge!=list[i].edge1 && matchEdgeType(v->twinEdge, list[i].edge2))
					{
						inEdges.push_back(v->twinEdge);
					}
					if(v!=list[i].edge2 && matchEdgeType(list[i].edge1, v))
					{
						outEdges.push_back(v);
					}
				}
			}
			inSupp=NULL;
			outSupp=NULL;
			// calculate supports for every pair of in-out edges
			if(inEdges.size()!=0)
			{
				if((inSupp=(int*)malloc(inEdges.size()*sizeof(int)))==NULL)
					printError(MEM_ALLOC, "inSupp");
				for(j=0; j<inEdges.size(); j++)
				{
					inSupp[j]=0;
					firstEdge=inEdges[j]->fromID ^ inEdges[j]->ID;
					secondEdge=list[i].edge2->fromID ^ list[i].edge2->ID;
					it = listSearch.find(make_pair(firstEdge, secondEdge));
					if(it != listSearch.end())
						inSupp[j]=it->second.support;
					else
					{
						firstEdge=list[i].edge2->twinEdge->fromID ^ list[i].edge2->twinEdge->ID;
						secondEdge=inEdges[j]->twinEdge->fromID ^ inEdges[j]->twinEdge->ID;
						it = listSearch.find(make_pair(firstEdge, secondEdge));
						if(it != listSearch.end())
							inSupp[j]=it->second.support;
					}
				}
			}
			if(outEdges.size()!=0)
			{
				if((outSupp=(int*)malloc(outEdges.size()*sizeof(int)))==NULL)
					printError(MEM_ALLOC, "outSupp");
				for(j=0; j<outEdges.size(); j++)
				{
					outSupp[j]=0;
					firstEdge=list[i].edge1->fromID ^ list[i].edge1->ID;
					secondEdge=outEdges[j]->fromID ^ outEdges[j]->ID;
					it = listSearch.find(make_pair(firstEdge, secondEdge));
					if(it != listSearch.end())
						outSupp[j]=it->second.support;
					else
					{
						firstEdge=outEdges[j]->twinEdge->fromID ^ outEdges[j]->twinEdge->ID;
						secondEdge=list[i].edge1->twinEdge->fromID ^ list[i].edge1->twinEdge->ID;
						it = listSearch.find(make_pair(firstEdge, secondEdge));
						if(it != listSearch.end())
							outSupp[j]=it->second.support;
					}
				}
			}

			uniqueIn=true;
			uniqueOut=true;
			for(j=0; j<outEdges.size(); j++)
			{
				if(outSupp[j]>lowThreshold)
				{	uniqueOut=false;
				}
			}
			//check in edges
			for(j=0; j<inEdges.size(); j++)
			{
				if(inSupp[j]>lowThreshold)
				{	uniqueIn=false;
				}
			}

			if(uniqueIn && uniqueOut)
			{
				for(v=graphObj->graph[nodeToResolve]; v!=NULL; v=v->next) // Insert all neighbours in the list
				{
					it2 = nodeSearch.find(v->ID);
					if(it2 == nodeSearch.end() && v->ID!=v->fromID)
						nodeSearch[v->ID]=1;
				}

				if(mergeEdgesAndUpdate(list[i].edge1, list[i].edge2, 1))
				{
					merged++;
					edgeID=list[i].edge1->ID;
					auto its = deleteList.equal_range(edgeID);
					for (auto it = its.first; it != its.second; ++it) 
					{	
						firstEdge=it->second.edge1->fromID ^ it->second.edge1->ID;
						secondEdge=it->second.edge2->fromID ^ it->second.edge2->ID;
						listSearch.erase(make_pair(firstEdge, secondEdge));
					}
				}			
			}
			free(inSupp);
			free(outSupp);
		}
		else
			break;
	}
	
	inEdges.clear();
	outEdges.clear();
	nodeSearch.clear();
	
	return merged;
}

/* ============================================================================================
   Two edges (e1, e2) and (e2,e3) is merged to (e1,e3) here.
   ============================================================================================ */
Edge* ContigExtender::mergeEdgesAndUpdate(Edge *edge1, Edge *edge2, int type_of_merge)
{
	int type=graphObj->combinedEdgeType(edge1, edge2);
	Edge *mergedEdge, *v;
	uint64_t index, dist;
	int flag, i;
	ReadToEdgeMap *edg, *edg_prev, *edg_next;
	float flow=0.0;
	if(edge1==edge2 || edge1==edge2->twinEdge || edge1->isReducible==0 || edge2->isReducible==0 || type==-1)
		return NULL;
	if(edge1->flow==0.0 && edge2->flow==0.0) // when we call this function before flow
	{
		flow=0.0;
	}
	else
	{
		if(type_of_merge==1) // only one unit of flow is used from both edges. this is used when we resolve pairs based on support
		{
			flow=1.0;
		}
		else if(type_of_merge==2) // we take minimum flow of the two edges. type_of_merge=2 when we resolve in and out tree and also when we resolve composite path
		{
			if(edge1->flow<=edge2->flow) // minimum flow in edge1
				flow=edge1->flow;
			else // minimum flow in edge2
				flow=edge2->flow;
		}
	}
	mergedEdge=graphObj->insertEdge(edge1->fromID, edge2->ID, type, graphObj->getListOfReads(edge1,edge2), graphObj->getListOfReads(edge2->twinEdge,edge1->twinEdge), flow, edge1->lengthOfEdge+edge2->lengthOfEdge, edge1->twinEdge->lengthOfEdge+edge2->twinEdge->lengthOfEdge);
	if(edge1->flow<=flow)
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
		edge1->twinEdge->flow-=flow;
		edge1->flow-=flow;
	}

	if(edge2->flow<=flow)
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
		edge2->twinEdge->flow-=flow;
		edge2->flow-=flow;
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
   This function merges long contigs that are not connected by overlaps.
   ============================================================================================ */
uint64_t ContigExtender::mergeShortOverlaps(uint16_t minOverlap, int highTh)
{
	uint64_t i, j, k, m, n, tos=0, number_merged=0, mate_pair_1, mate_pair_2, dist, total_edge_pair_to_merge=0, numberOfFeasibleEdges, arraySize=10000, flag_distance, index, edgeOne=0, edgeTwo=0, edgeThree=0, edgeFour=0;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> supportedListSearch;
	unordered_multimap<uint64_t,SupportedEdgePair> deleteList;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash>::const_iterator it;
	uint64_t firstEdge=0, secondEdge=0; // for calculating amount to subtract from memory usage
	int32_t *gapDistance, gapSum, support, distance_all, *distance_1, *distance_2;
	int *pair_support, type1a=0, type2a=0;
	Edge **listOfLongEdges, *edge1=NULL, *edge2=NULL, **edge_1, **edge_2, **feasibleListOfEdges;
	MatePairInfo *w;
	listOfLongEdges=getListOfCompositeEdges1(&tos);
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
		feasibleListOfEdges=getListOfFeasibleEdges1(listOfLongEdges[i], &numberOfFeasibleEdges);
		for(j=0; j<numberOfFeasibleEdges; j++)
		{
			int64_t len1=listOfLongEdges[i]->lengthOfEdge, len2=feasibleListOfEdges[j]->lengthOfEdge;
			if((listOfLongEdges[i]->flow>0 && len1>contigLengthThreshold) && (feasibleListOfEdges[j]->flow>0 && len2>contigLengthThreshold))
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

					if(support>0 && gapSum<=0)
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
	if(total_edge_pair_to_merge>0)
	{
		quicksortListOfLongEdges1(edge_1, edge_2, pair_support, gapDistance, 0, total_edge_pair_to_merge-1); // Sort edges according to the support.
		number_merged=mergeShortOverlap(edge_1, edge_2, pair_support, gapDistance, total_edge_pair_to_merge, supportedListSearch, deleteList, highTh);
	}

	free(edge_1);
	free(edge_2);
	free(pair_support);
	free(gapDistance);
	free(listOfLongEdges);
	return number_merged;
}

/* ============================================================================================
   This function returns the list of all composite edges in the graph.
   ============================================================================================ */
Edge** ContigExtender::getListOfCompositeEdges1(uint64_t *number)
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
Edge** ContigExtender::getListOfFeasibleEdges1(Edge *edge, uint64_t *number)
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
   Merge edge pairs according to support.
   ============================================================================================ */
uint64_t ContigExtender::mergeShortOverlap(Edge **edge_1,Edge **edge_2,int *pair_support,int32_t *gapDistance, uint64_t total_edge_pair_to_merge, unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> &listSearch, unordered_multimap<uint64_t,SupportedEdgePair> &deleteList, int highTh)
{
	uint64_t number_merged=0, edgeID=0;
	uint64_t i, j, firstEdge=0, secondEdge=0;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash>::const_iterator it;
	Edge *edge1_p, *edge2_p, *edge1_t_p, *edge2_t_p;
	Edge **listOfEdges1, **listOfEdges2;
	uint64_t numOfList1, numOfList2;
	int *suppList1, *suppList2;
	bool unique1, unique2;
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
			listOfEdges1 = getListOfFeasibleEdges1(edge_2[i], &numOfList1);
			listOfEdges2 = getListOfFeasibleEdges1(edge_1[i], &numOfList2);
			suppList1=NULL;
			suppList2=NULL;
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
			unique1=true;
			unique2=true;
			//check second list
			for(j=0; j<numOfList2; j++)
			{
				if(listOfEdges2[j]==edge_2[i])
					continue;
				if(suppList2[j] >= highTh)
				{	unique2=false;
				}
			}
			//check in edges
			for(j=0; j<numOfList1; j++)
			{
				if(listOfEdges1[j]==edge_1[i])
					continue;
				if(suppList1[j] >= highTh)
				{	unique1=false;
				}
			}
			if(unique1 && unique2)
			{
				number_merged++;
				edge1_p=edge_1[i];
				edge2_p=edge_2[i];
				edge1_t_p=edge_1[i]->twinEdge;
				edge2_t_p=edge_2[i]->twinEdge;
				if(mergeDisconnectedEdges1(edge_1[i], edge_2[i], gapDistance[i]))
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
   Two edges (e1, e2) and (e2,e3) is merged to (e1,e3) here.
   ============================================================================================ */
Edge* ContigExtender::mergeDisconnectedEdges1(Edge *edge1, Edge *edge2, int gap)
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
   Sort list of long edges according to support. If the supports are same then sort by length.
   ============================================================================================ */
void ContigExtender::quicksortListOfLongEdges1(Edge **edge_1, Edge **edge_2, int * pair_support, int32_t *gapDistance, int64_t left, int64_t right)
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
		quicksortListOfLongEdges1(edge_1,edge_2,pair_support,gapDistance,left,j);
    if (i < right)
		quicksortListOfLongEdges1(edge_1,edge_2,pair_support,gapDistance,i,right);
}


/* ============================================================================================
   This function will find paths from each node in the graph, one by one
   ============================================================================================ */
uint64_t ContigExtender::findPathSupportsNew(int highTh)
{
	uint64_t i, returnValue=0;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> supportedListSearch;
	unordered_multimap<uint64_t,SupportedEdgePair> deleteList;
	bool *flag;//this is the flag whether we already found paths of the mate pairs or not.
	if((flag=(bool*)malloc((loaderObj->numberOfUniqueReads+1)*sizeof(bool)))==NULL)
		printError(MEM_ALLOC,"flag");

	for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
		flag[i]=1;

	#pragma omp parallel
	{		
		uint64_t mate_pair_1, mate_pair_2, total_supported_pairs, j, l, index, firstEdge=0, secondEdge=0, edgeOne=0, edgeTwo=0, edgeThree=0, edgeFour=0, distance=0, no_of_reads=0, last=0, arraySize = 20000, totalPaths = 0, readsHigh=0, pairSuppHigh=0;

		unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> supportedListSearchThread;
		unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash>::const_iterator it;
		unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash>::const_iterator it1;
		unordered_map<uint64_t,int> list_of_explored;
		unordered_map<uint64_t,int>::const_iterator it2;

		MatePairInfo *w;
		Edge *v;
		ReadToEdgeMap *a;
		uint64_t *list_of_reads;		
		PairSupported *supportedPairs_tmp;
		Edge **pathEdges, ***allPairPaths;
		uint64_t *pathLengths, *allPairPathLengths;
		Indexes **indx;
		if((list_of_reads=(uint64_t*)malloc(arraySize*sizeof(uint64_t)))==NULL)
			printError(MEM_ALLOC,"list_of_reads");
		if((pathLengths=(uint64_t*)malloc((20)*sizeof(uint64_t)))==NULL)
			printError(MEM_ALLOC, "pathLengths");
		if((pathEdges=(Edge**)malloc((20)*sizeof(Edge*)))==NULL)
			printError(MEM_ALLOC, "pathEdges");
		if((allPairPaths=(Edge***)malloc(2000*sizeof(Edge**)))==NULL)
			printError(MEM_ALLOC,"allPairPaths");
		if((allPairPathLengths=(uint64_t*)malloc(2000*sizeof(uint64_t)))==NULL)
			printError(MEM_ALLOC,"allPairPathLengths");
		for(j=0;j<2000;j++)
		{	
			if((allPairPaths[j]=(Edge**)malloc((10)*sizeof(Edge*)))==NULL)
				printError(MEM_ALLOC, "allPairPaths[totalPaths]");
		}		

		#pragma omp for
		for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
		{
			
			if(graphObj->graph[i]!=NULL) //for each node in the graph
			{
				for(j=0;j<arraySize;j++)
					list_of_reads[j]=0;
				for(j=0;j<20;j++)
				{	
					pathLengths[j]=0;
					pathEdges[j]=NULL;
				}
				for(j=0;j<2000;j++)
				{
					allPairPathLengths[j]=0;
					for(l=0;l<10;l++)
					{
						allPairPaths[j][l]=NULL;
					}
				}

				totalPaths=0;
				last=0;
				no_of_reads=0;

				for(v=graphObj->graph[i]; v!=NULL; v=v->next)
				{
					distance=0;
					for(index=1; index<=v->listOfReads[0].readId; index++)
					{
						distance += v->listOfReads[index].distPrevious;
						if(distance>mateObj->maximumUpperBoundOfInsert)
							break;
						if(flag[v->listOfReads[index].readId]==1)
						{
							list_of_reads[no_of_reads++]=v->listOfReads[index].readId;
							flag[v->listOfReads[index].readId]=0;
						}
					}
					distance=0;
					for(index=1; index<=v->twinEdge->listOfReads[0].readId; index++)
					{
						distance += v->twinEdge->listOfReads[index].distPrevious;
						if(distance>mateObj->maximumUpperBoundOfInsert)
							break;
						if(flag[v->twinEdge->listOfReads[index].readId]==1)
						{
							list_of_reads[no_of_reads++]=v->twinEdge->listOfReads[index].readId;
							flag[v->twinEdge->listOfReads[index].readId]=0;
						}
					}
				}
				if(readsHigh<no_of_reads)
						readsHigh=no_of_reads;
					
				for(l=0; l<no_of_reads; l++)
				{
					for(a=mateObj->readToEdgeList[list_of_reads[l]]; a!=NULL; a=a->next)
					{
						it2 = list_of_explored.find((uint64_t)a->edge->fromID);
						if(it2 == list_of_explored.end())
						{
							list_of_explored[(uint64_t)a->edge->fromID] = 1;
							last++;
							exploreGraph(a->edge->fromID, 1, totalPaths, pathEdges, pathLengths, allPairPaths, allPairPathLengths);
						}
						it2 = list_of_explored.find((uint64_t)a->edge->ID);
						if(it2 == list_of_explored.end())
						{
							list_of_explored[(uint64_t)a->edge->ID] = 1;
							last++;
							exploreGraph(a->edge->ID, 1, totalPaths, pathEdges, pathLengths, allPairPaths, allPairPathLengths);
						}
					}
				}
				list_of_explored.clear();
				
				if(totalPaths==0 || last==0 || no_of_reads==0)
					continue;

				quicksortPaths(0, totalPaths-1, allPairPaths, allPairPathLengths);
				if((indx = (Indexes**)malloc((last+2)*sizeof(Indexes*)))==NULL)
					printError(MEM_ALLOC, "indx");

				indexPaths(last, indx, totalPaths, allPairPaths);
				
				if((supportedPairs_tmp=(PairSupported*)malloc(arraySize*sizeof(PairSupported)))==NULL)
					printError(MEM_ALLOC, "supportedPairs_tmp");

				for(l=0; l<no_of_reads; l++)
				{
					mate_pair_1=list_of_reads[l];
					for(w=mateObj->matePairList[mate_pair_1]; w!=NULL; w=w->next)
					{
						mate_pair_2=w->ID;
						if(w->flag==1) 
						{
							if(graphObj->graph[mate_pair_1]==NULL && graphObj->graph[mate_pair_2]==NULL && flag[mate_pair_2]==1) 
							{
								if(totalPaths==0)
									total_supported_pairs=0;
								else
									total_supported_pairs=findPathsBetweenMatepairs(mate_pair_1, mate_pair_2, w->type1, w->type2, w->library, indx, totalPaths, allPairPaths, allPairPathLengths, supportedPairs_tmp);

								if(pairSuppHigh<total_supported_pairs)
									pairSuppHigh=total_supported_pairs;
								
								if(total_supported_pairs>0)
								{
									for(j=0; j<total_supported_pairs; j++)
									{
										SupportedEdgePair pairSupp;
										pairSupp.edge1 = supportedPairs_tmp[j].edge1;
										pairSupp.edge2 = supportedPairs_tmp[j].edge2;
										edgeOne = pairSupp.edge1->ID;
										edgeTwo = pairSupp.edge1->fromID;
										edgeThree = pairSupp.edge2->ID;
										edgeFour = pairSupp.edge2->fromID;
										firstEdge = supportedPairs_tmp[j].edge1->fromID ^ supportedPairs_tmp[j].edge1->ID;
										secondEdge = supportedPairs_tmp[j].edge2->fromID ^ supportedPairs_tmp[j].edge2->ID;
										pairSupp.support = 1;
										it = supportedListSearchThread.find(make_pair(firstEdge, secondEdge));
										if(it != supportedListSearchThread.end())
										{
											pairSupp.support = it->second.support+1;
											supportedListSearchThread.erase(make_pair(firstEdge, secondEdge));
											supportedListSearchThread[make_pair(firstEdge, secondEdge)] = pairSupp;
										}
										else
											supportedListSearchThread[make_pair(firstEdge, secondEdge)] = pairSupp;
									}
								}
							}
						}
					}
				}
				freeIndexedPaths(indx);
				free(indx);
				free(supportedPairs_tmp);
			}
		}
		
		#pragma omp critical (addToList)
		{
			for(it = supportedListSearchThread.begin(); it != supportedListSearchThread.end(); ++it)
			{
				SupportedEdgePair pairSupp;
				pairSupp.edge1 = it->second.edge1;
				pairSupp.edge2 = it->second.edge2;
				pairSupp.support = it->second.support;
				firstEdge = it->second.edge1->fromID ^ it->second.edge1->ID;
				secondEdge = it->second.edge2->fromID ^ it->second.edge2->ID;
				it1 = supportedListSearch.find(make_pair(firstEdge, secondEdge));
				// insert the pair into the unordered_map
				if(it1 != supportedListSearch.end())
				{
					pairSupp.support = it1->second.support + it->second.support;
					supportedListSearch.erase(make_pair(firstEdge, secondEdge));
					supportedListSearch[make_pair(firstEdge, secondEdge)] = pairSupp;
				}
				else
					supportedListSearch[make_pair(firstEdge, secondEdge)] = pairSupp;

				edgeOne = pairSupp.edge1->ID;
				edgeTwo = pairSupp.edge1->fromID;
				edgeThree = pairSupp.edge2->ID;
				edgeFour = pairSupp.edge2->fromID;
				// insert the edges into the unordered_multimap
				deleteList.insert({edgeOne,pairSupp});
				deleteList.insert({edgeTwo,pairSupp});
				deleteList.insert({edgeThree,pairSupp});
				deleteList.insert({edgeFour,pairSupp});
			}
			supportedListSearchThread.clear();
			for(j=0;j<2000;j++)
				free(allPairPaths[j]);
			free(allPairPaths);
			free(allPairPathLengths);
			free(pathLengths);
			free(pathEdges);
			free(list_of_reads);
		}
		
	}
	free(flag);		
	returnValue=mergeByPathSupportsNew(supportedListSearch, deleteList, highTh);
	
	return returnValue;
}

/* ============================================================================================
   Merge pair of edges according to supports.
   ============================================================================================ */
uint64_t ContigExtender::mergeByPathSupportsNew(unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash> &listSearch, unordered_multimap<uint64_t,SupportedEdgePair> &deleteList, int highTh)
{
	uint64_t j, node_resolved=0,firstEdge=0,secondEdge=0;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash>::const_iterator it;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash>::const_iterator it1, next_it;
	unordered_map<pair<uint64_t,uint64_t>,SupportedEdgePair,pairHash>::const_iterator it3;
	unordered_map<uint64_t,int> nodeSearch;
	unordered_map<uint64_t,int>::const_iterator it2;
	Edge *v;
	deque<Edge*> inEdges;
	deque<Edge*> outEdges;
	int *inSupp, *outSupp;
	bool uniqueIn, uniqueOut, flag=false;
	uint64_t nodeToResolve;
	for(it1 = listSearch.begin(); it1 != listSearch.end(); ++it1)
	{
		if (flag)
		{
			if(next_it != listSearch.end())
			{   
				it1 = next_it;
				flag=false;
			}else
				break;
		}

		firstEdge = it1->second.edge1->fromID ^ it1->second.edge1->ID;
		secondEdge = it1->second.edge2->fromID ^ it1->second.edge2->ID;
		it3 = listSearch.find(make_pair(firstEdge, secondEdge));
		if(it3 == listSearch.end())
				continue;

		nodeToResolve=it1->second.edge1->ID;
		it2 = nodeSearch.find(nodeToResolve);
		if(it2 != nodeSearch.end())
				continue;

		if(it1->second.support >= highTh)
		{
				inEdges.clear();
				outEdges.clear();
				for(v=graphObj->graph[nodeToResolve]; v!=NULL; v=v->next)
				{
						if(v->flow>=1 && v->ID!=v->fromID)
						{
								if(v->twinEdge!=it1->second.edge1 && matchEdgeType(v->twinEdge, it1->second.edge2))
								{
										inEdges.push_back(v->twinEdge);
								}
								if(v!=it1->second.edge2 && matchEdgeType(it1->second.edge1, v))
								{
										outEdges.push_back(v);
								}
						}
				}
				inSupp=NULL;
				outSupp=NULL;
			   // calculate supports for every pair of in-out edges
			   if(inEdges.size()!=0)
				{
						if((inSupp=(int*)malloc(inEdges.size()*sizeof(int)))==NULL)
								printError(MEM_ALLOC, "inSupp");
						for(j=0; j<inEdges.size(); j++)
						{
								inSupp[j]=0;
								firstEdge=inEdges[j]->fromID ^ inEdges[j]->ID;
								secondEdge=it1->second.edge2->fromID ^ it1->second.edge2->ID;
								it = listSearch.find(make_pair(firstEdge, secondEdge));
								if(it != listSearch.end())
										inSupp[j]=it->second.support;
								else
								{
										firstEdge=it1->second.edge2->twinEdge->fromID ^ it1->second.edge2->twinEdge->ID;
										secondEdge=inEdges[j]->twinEdge->fromID ^ inEdges[j]->twinEdge->ID;
										it = listSearch.find(make_pair(firstEdge, secondEdge));
										if(it != listSearch.end())
												inSupp[j]=it->second.support;
								}
						}
				}
				if(outEdges.size()!=0)
				{
						if((outSupp=(int*)malloc(outEdges.size()*sizeof(int)))==NULL)
								printError(MEM_ALLOC, "outSupp");
						for(j=0; j<outEdges.size(); j++)
						{
								outSupp[j]=0;
								firstEdge=it1->second.edge1->fromID ^ it1->second.edge1->ID;
								secondEdge=outEdges[j]->fromID ^ outEdges[j]->ID;
								it = listSearch.find(make_pair(firstEdge, secondEdge));
								if(it != listSearch.end())
										outSupp[j]=it->second.support;
								else
								{
										firstEdge=outEdges[j]->twinEdge->fromID ^ outEdges[j]->twinEdge->ID;
										secondEdge=it1->second.edge1->twinEdge->fromID ^ it1->second.edge1->twinEdge->ID;
										it = listSearch.find(make_pair(firstEdge, secondEdge));
										if(it != listSearch.end())
												outSupp[j]=it->second.support;
								}
						}
				}
				uniqueIn=true;
				uniqueOut=true;
				for(j=0; j<outEdges.size(); j++)
				{
						if(outSupp[j]>lowThreshold)
						{       uniqueOut=false;
					   }
				}
				for(j=0; j<inEdges.size(); j++)
				{
						if(inSupp[j]>lowThreshold)
						{       uniqueIn=false;
						}
				}
				if((uniqueIn && uniqueOut) || (uniqueIn==true && uniqueOut==false) || (uniqueIn==false && uniqueOut==true))
				{
						for(v=graphObj->graph[nodeToResolve]; v!=NULL; v=v->next)
						{
								it2 = nodeSearch.find(v->ID);
								if(it2 == nodeSearch.end() && v->ID!=v->fromID)
										nodeSearch[v->ID]=1;
						}
						//set next_it
						next_it=it1;
						++next_it;
						while(next_it != listSearch.end())
						{
							if (next_it->second.edge1->fromID == it1->second.edge1->ID ||
									next_it->second.edge1->ID == it1->second.edge1->ID ||
									next_it->second.edge2->fromID == it1->second.edge1->ID ||
									next_it->second.edge2->ID == it1->second.edge1->ID)
							{
								++next_it;
							}else
								break;
						}

						if(mergeEdgesAndUpdate(it1->second.edge1, it1->second.edge2, 1))
						{
								flag=true;
								node_resolved++;
								auto its = deleteList.equal_range(nodeToResolve);
								for (auto it = its.first; it != its.second; ++it)
								{
										firstEdge=it->second.edge1->fromID ^ it->second.edge1->ID;
										secondEdge=it->second.edge2->fromID ^ it->second.edge2->ID;
										listSearch.erase(make_pair(firstEdge, secondEdge));

								}
						}
				}
				free(inSupp);
				free(outSupp);
		}
	}
	return node_resolved;
}
