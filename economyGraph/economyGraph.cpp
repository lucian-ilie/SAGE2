/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 ****************************************************************************/

#include "economyGraph.h"
extern ofstream logStream;

EconomyEdge::EconomyEdge(uint64_t dest1, uint8_t type1, uint8_t mark1, uint32_t len1)
{
	readId = dest1;
	type = type1;
	mark = mark1;
	length = len1;
}

EconomyGraph::EconomyGraph(uint16_t minOvlp, HashTable *hash1)
{
	exploredReads=NULL;
	markedNodes=NULL;
	economyGraphList=NULL;
	numberOfHashMiss=0;
	minOverlap = minOvlp;
	hashObj = hash1;
	loaderObj = hash1->loaderObj;
}

EconomyGraph::~EconomyGraph()
{
	free(economyGraphList);
}

/* ============================================================================================
   This function builds unique extensions in the overlap graph from hash table.
   ============================================================================================ */
void EconomyGraph::buildInitialOverlapGraph()
{
	logStream<< "In function buildInitialOverlapGraph().\n";
	logStream.flush();
	
	time_t seconds_s=time(NULL);
	uint64_t i, connectionsLimit=300;
	uint16_t hashStringLength = hashObj->hashStringLength;

	ExtensionTable *rightExtension = new ExtensionTable[loaderObj->numberOfUniqueReads+1];
	ExtensionTable *leftExtension = new ExtensionTable[loaderObj->numberOfUniqueReads+1];

	if((exploredReads=(uint8_t*)malloc((loaderObj->numberOfUniqueReads+1)*sizeof(uint8_t)))==NULL)
		printError(MEM_ALLOC, "exploredReads");
	
	if((economyGraphList=(EconomyEdge**)malloc((loaderObj->numberOfUniqueReads+1)*sizeof(EconomyEdge*)))==NULL)
		printError(MEM_ALLOC, "graphEconomy");

	#pragma omp parallel for
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++) // Initialize all flags
	{
		exploredReads[i]=0;
		rightExtension[i].readId=0;
		leftExtension[i].readId=0;
		economyGraphList[i]=NULL;
	}
	
	#pragma omp parallel
	{
		Read *r1, *r2, *previousRead;
		uint64_t k,read2, prevIDRight=0, prevIDLeft=0, connections=0;
		uint64_t *uint64Window;
		int64_t index,j;
		int markAmbigRight=0, markAmbigLeft=0, type=0, itsAmbigRight=0, itsAmbigLeft=0, prevTypeRight=0, prevLengthRight=0, prevTypeLeft=0, prevLengthLeft=0, markFirstRight=0;
		#pragma omp for
		for(i=1; i<=loaderObj->numberOfUniqueReads; i++) //for each R_i find unique extension on each side.
		{		
			r1 = loaderObj->getRead(i);	
			itsAmbigRight=0;
			itsAmbigLeft=0;
			connections=0;
			// find the unique extensions of the read
			for(j=0; j<=(int64_t)(r1->length-hashStringLength); j++)
			{			
				uint64Window=get64Bit2Int(r1->readInt, j, hashStringLength);
				index=hashObj->hashTableSearch(uint64Window, numberOfHashMiss);
				free((uint64_t*) uint64Window);
				if(index!=-1) // found in hash table
				{
					markAmbigRight=0;
					markAmbigLeft=0;
					markFirstRight=0;
					for(k=1; k<=hashObj->hashTableList[index][0].readId; k++)
					{	
						read2=hashObj->hashTableList[index][k].readId;
						r2 = loaderObj->getRead(read2);
						type=hashObj->hashTableList[index][k].type;
						if(type==0 && read2!=i && j<=(uint16_t)(r1->length-minOverlap) && compareStringInBytes(r1->readInt, r2->readInt, j, r1->length, r2->length, read2)) // check if the overlap extends to the end of the read
						{
							connections++;
							if(rightExtension[i].readId==0) //first read on this side of the extension
							{
								rightExtension[i].readId=read2;
								rightExtension[i].type=0;							
								rightExtension[i].length=r2->length-(r1->length-j);
								prevIDRight=read2;
								prevTypeRight=0;
								prevLengthRight=j;
								markAmbigRight=1;
								markFirstRight=1;
							}
							else // not the first read to extend this side
							{
								previousRead=loaderObj->getRead(prevIDRight);
								if(prevTypeRight==0) // the previous read is in the forward direction
								{
									//check that the previous alignment is consistent
									if(compareStringInBytesPrevious(previousRead->readInt, r2->readInt, j-prevLengthRight, previousRead->length, r2->length))
									{	
										if(markAmbigRight==1)
										{	
											if(r2->length > previousRead->length)
											{												
												if(markFirstRight==1)
												{	rightExtension[i].readId=read2;
													rightExtension[i].type=0;							
													rightExtension[i].length=r2->length-(r1->length-j);
													prevIDRight=read2;
													prevTypeRight=0;
													prevLengthRight=j;
												}
												else
												{	
													prevIDRight=read2;
													prevTypeRight=0;
													prevLengthRight=j;
												}
													
											}
										}
										else
										{											
											prevIDRight=read2;
											prevTypeRight=0;
											prevLengthRight=j;
											markAmbigRight=1;
										}
									}
									else // this in an ambiguous extension
										itsAmbigRight=1;
								}
								else // the previous read is in the reverse direction
								{
									//check that the previous alignment is consistent
									if(compareStringInBytesPrevious(previousRead->readReverseInt, r2->readInt, j-prevLengthRight, previousRead->length, r2->length))
									{	
										if(markAmbigRight==1)
										{	
											if(r2->length > previousRead->length)
											{												
												if(markFirstRight==1)
												{	rightExtension[i].readId=read2;
													rightExtension[i].type=0;							
													rightExtension[i].length=r2->length-(r1->length-j);
													prevIDRight=read2;
													prevTypeRight=0;
													prevLengthRight=j;
												}
												else
												{	
													prevIDRight=read2;
													prevTypeRight=0;
													prevLengthRight=j;
												}
													
											}
										}
										else
										{											
											prevIDRight=read2;
											prevTypeRight=0;
											prevLengthRight=j;
											markAmbigRight=1;
										}
									}
									else // this in an ambiguous extension
										itsAmbigRight=1;
								}
							}
						}					
						else if(type==2 && read2!=i && j<=(uint16_t)(r1->length-minOverlap) && compareStringInBytes(r1->readInt, r2->readReverseInt, j, r1->length, r2->length, read2)) // check if the overlap extends to the end of the read in reverse
						{
							connections++;
							if(rightExtension[i].readId==0)
							{
								rightExtension[i].readId=read2;
								rightExtension[i].type=1;
								rightExtension[i].length=r2->length-(r1->length-j);
								prevIDRight=read2;
								prevTypeRight=1;
								prevLengthRight=j;
								markAmbigRight=1;
								markFirstRight=1;
							}
							else  // not the first read to extend this side
							{
								previousRead=loaderObj->getRead(prevIDRight);
								if(prevTypeRight==0) // the previous read is in the forward direction
								{
									//check that the previous alignment is consistent
									if(compareStringInBytesPrevious(previousRead->readInt, r2->readReverseInt, j-prevLengthRight, previousRead->length, r2->length))
									{	
										if(markAmbigRight==1)
										{	
											if(r2->length > previousRead->length)
											{												
												if(markFirstRight==1)
												{	rightExtension[i].readId=read2;
													rightExtension[i].type=1;							
													rightExtension[i].length=r2->length-(r1->length-j);
													prevIDRight=read2;
													prevTypeRight=1;
													prevLengthRight=j;
												}
												else
												{	
													prevIDRight=read2;
													prevTypeRight=1;
													prevLengthRight=j;
												}													
											}
										}
										else
										{											
											prevIDRight=read2;
											prevTypeRight=1;
											prevLengthRight=j;
											markAmbigRight=1;
										}
									}
									else // this in an ambiguous extension
										itsAmbigRight=1;
								}
								else// the previous read is in the reverse direction
								{
									//check that the previous alignment is consistent
									if(compareStringInBytesPrevious(previousRead->readReverseInt, r2->readReverseInt, j-prevLengthRight, previousRead->length, r2->length))
									{	
										if(markAmbigRight==1)
										{	
											if(r2->length > previousRead->length)
											{												
												if(markFirstRight==1)
												{	rightExtension[i].readId=read2;
													rightExtension[i].type=1;							
													rightExtension[i].length=r2->length-(r1->length-j);
													prevIDRight=read2;
													prevTypeRight=1;
													prevLengthRight=j;
												}
												else
												{	
													prevIDRight=read2;
													prevTypeRight=1;
													prevLengthRight=j;
												}													
											}
										}
										else
										{											
											prevIDRight=read2;
											prevTypeRight=1;
											prevLengthRight=j;
											markAmbigRight=1;
										}
									}
									else // this in an ambiguous extension
										itsAmbigRight=1;
								}
							}

						}				
						else if(type==1 && read2!=i && j>=(uint16_t)(minOverlap-hashStringLength) && compareStringInBytes(r1->readReverseInt, r2->readReverseInt, r1->length-j-hashStringLength, r1->length, r2->length, read2)) // check if the overlap extends to the end of the read in reverse
						{
							connections++;
							if(leftExtension[i].readId==0)
							{
								leftExtension[i].readId=read2;
								leftExtension[i].type=0;
								leftExtension[i].length=r2->length-j-hashStringLength;
								prevIDLeft=read2;
								prevTypeLeft=0;
								prevLengthLeft=r1->length-j-hashStringLength;
								markAmbigLeft=1;
							}
							else  // not the first read to extend this side
							{
								previousRead=loaderObj->getRead(prevIDLeft);
								if(prevTypeLeft==0) // the previous read is in the forward direction
								{
									//check that the previous alignment is consistent
									if(compareStringInBytesPrevious(r2->readReverseInt, previousRead->readReverseInt, prevLengthLeft-(r1->length-j-hashStringLength), r2->length, previousRead->length))
									{	
										if(markAmbigLeft==1)
										{
											if(r2->length > previousRead->length)
											{
												leftExtension[i].readId=read2;
												leftExtension[i].type=0;
												leftExtension[i].length=r2->length-j-hashStringLength;							
												prevIDLeft=read2;
												prevTypeLeft=0;
												prevLengthLeft=r1->length-j-hashStringLength;
											}
										}
										else
										{	
											leftExtension[i].readId=read2;
											leftExtension[i].type=0;
											leftExtension[i].length=r2->length-j-hashStringLength;																
											prevIDLeft=read2;
											prevTypeLeft=0;
											prevLengthLeft=r1->length-j-hashStringLength;
											markAmbigLeft=1;
										}
									}
									else // this in an ambiguous extension
										itsAmbigLeft=1;
								}
								else // the previous read is in the reverse direction
								{
									//check that the previous alignment is consistent
									if(compareStringInBytesPrevious(r2->readReverseInt, previousRead->readInt, prevLengthLeft-(r1->length-j-hashStringLength), r2->length, previousRead->length))
									{	
										if(markAmbigLeft==1)
										{
											if(r2->length > previousRead->length)
											{
												leftExtension[i].readId=read2;
												leftExtension[i].type=0;
												leftExtension[i].length=r2->length-j-hashStringLength;							
												prevIDLeft=read2;
												prevTypeLeft=0;
												prevLengthLeft=r1->length-j-hashStringLength;
											}
										}
										else
										{	
											leftExtension[i].readId=read2;
											leftExtension[i].type=0;
											leftExtension[i].length=r2->length-j-hashStringLength;																
											prevIDLeft=read2;
											prevTypeLeft=0;
											prevLengthLeft=r1->length-j-hashStringLength;
											markAmbigLeft=1;
										}
									}
									else // this in an ambiguous extension
										itsAmbigLeft=1;
								}
							}
						}
						else if(type==3 && read2!=i && j>=(uint16_t)(minOverlap-hashStringLength) && compareStringInBytes(r1->readReverseInt, r2->readInt, r1->length-j-hashStringLength, r1->length, r2->length, read2)) // check if the overlap extends to the end of the read
						{
							connections++;
							if(leftExtension[i].readId==0) //first read on this side of the extension
							{
								leftExtension[i].readId=read2;
								leftExtension[i].type=1;
								leftExtension[i].length=r2->length-j-hashStringLength;
								prevIDLeft=read2;
								prevTypeLeft=1;
								prevLengthLeft=r1->length-j-hashStringLength;
								markAmbigLeft=1;
							}
							else // not the first read to extend this side
							{
								previousRead=loaderObj->getRead(prevIDLeft);
								if(prevTypeLeft==0) // the previous read is in the forward direction
								{
									//check that the previous alignment is consistent
									if(compareStringInBytesPrevious(r2->readInt, previousRead->readReverseInt, prevLengthLeft-(r1->length-j-hashStringLength), r2->length, previousRead->length))
									{	
										if(markAmbigLeft==1)
										{
											if(r2->length > previousRead->length)
											{
												leftExtension[i].readId=read2;
												leftExtension[i].type=1;
												leftExtension[i].length=r2->length-j-hashStringLength;							
												prevIDLeft=read2;
												prevTypeLeft=1;
												prevLengthLeft=r1->length-j-hashStringLength;
											}
										}
										else
										{	
											leftExtension[i].readId=read2;
											leftExtension[i].type=1;
											leftExtension[i].length=r2->length-j-hashStringLength;																
											prevIDLeft=read2;
											prevTypeLeft=1;
											prevLengthLeft=r1->length-j-hashStringLength;
											markAmbigLeft=1;
										}
									}
									else // this in an ambiguous extension
										itsAmbigLeft=1;
								}
								else // the previous read is in the reverse direction
								{
									//check that the previous alignment is consistent
									if(compareStringInBytesPrevious(r2->readInt, previousRead->readInt, prevLengthLeft-(r1->length-j-hashStringLength), r2->length, previousRead->length))
									{	
										if(markAmbigLeft==1)
										{
											if(r2->length > previousRead->length)
											{
												leftExtension[i].readId=read2;
												leftExtension[i].type=1;
												leftExtension[i].length=r2->length-j-hashStringLength;							
												prevIDLeft=read2;
												prevTypeLeft=1;
												prevLengthLeft=r1->length-j-hashStringLength;
											}
										}
										else
										{	
											leftExtension[i].readId=read2;
											leftExtension[i].type=1;
											leftExtension[i].length=r2->length-j-hashStringLength;																
											prevIDLeft=read2;
											prevTypeLeft=1;
											prevLengthLeft=r1->length-j-hashStringLength;
											markAmbigLeft=1;
										}
									}
									else // this in an ambiguous extension
										itsAmbigLeft=1;
								}
							}
						}
					}
				}
			}				

			if(connections>connectionsLimit)
				exploredReads[i]=5;

			if(itsAmbigRight==1 || itsAmbigLeft==1)
			{
				leftExtension[i].length=0;
				rightExtension[i].length=0;
			}
		}
	}
	
	uint64_t contained=0, containedSize=0;
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
	{
		if(exploredReads[i]!=6)
		{
			// find if there is a reciprocal extension to the left and right
			if((leftExtension[i].length!=0 && (rightExtension[leftExtension[i].readId].readId == i || leftExtension[leftExtension[i].readId].readId == i)) && (rightExtension[i].length!=0 && (rightExtension[rightExtension[i].readId].readId == i || leftExtension[rightExtension[i].readId].readId == i)))
			{	
				if(exploredReads[leftExtension[i].readId]!=4)
				{	if(leftExtension[i].type == 0)
						insertEdgeEconomy(i, leftExtension[i].readId, leftExtension[i].length, 0);
					else
						insertEdgeEconomy(i, leftExtension[i].readId, leftExtension[i].length, 1);
				}
				if(exploredReads[rightExtension[i].readId]!=4)
				{	if(rightExtension[i].type == 0)
						insertEdgeEconomy(i, rightExtension[i].readId, rightExtension[i].length, 3);
					else
						insertEdgeEconomy(i, rightExtension[i].readId, rightExtension[i].length, 2);
				}
				contained++;
				exploredReads[i]=4;
			}	
		}
		else
			containedSize++;
	}
	
	delete [] rightExtension;
	delete [] leftExtension;
	
	logStream<< "     Total contained by extension: " << contained << "\n";
	logStream<< "          Total contained by size: " << containedSize << "\n";
	logStream<< "            Total left to explore: " << loaderObj->numberOfUniqueReads - contained - containedSize  << "\n";
	logStream<< "Function buildInitialOverlapGraph() in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
}

/* ============================================================================================
   This function creates the overlap graph from hash table.
   ============================================================================================ */
void EconomyGraph::buildOverlapGraphEconomy()
{
	
	logStream<< "In function buildOverlapGraphEconomy().\n";
	logStream.flush();
	uint64_t i=0, totalEdgeInserted=0, transitiveEdgesRemoved=0;
	uint64_t *queue, start=0, end=0, read1=0, read2=0, read3=0, index1=0, index2=0;

	time_t seconds_s=time(NULL);
	if((markedNodes=(uint8_t*)malloc((loaderObj->numberOfUniqueReads+1)*sizeof(uint8_t)))==NULL)
		printError(MEM_ALLOC, "markedNodes");
	if((queue=(uint64_t *)malloc((loaderObj->numberOfUniqueReads+1)*sizeof(uint64_t)))==NULL)
		printError(MEM_ALLOC, "queue");
	
	#pragma omp parallel for
	for(i=0; i<=loaderObj->numberOfUniqueReads; i++) // Initialize all flags
		markedNodes[i]=0;

	for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
	{
		if(exploredReads[i]==0)
		{
			start=0; // Initialize queue start and end
			end=0;
			queue[end++]=i; // Put the read in the queue
			while(start<end) // This loop will explore all connected component starting from read i
			{
				read1=queue[start++];
				if(exploredReads[read1]==0)
					totalEdgeInserted+=insertAllEdgesOfRead(read1); // Explore current node
				if(economyGraphList[read1]!=NULL) // Read has some edges (required only for the first read when a new queue starts
				{
					if(exploredReads[read1]==1) // Explore neighbours first
					{
						for(index1=1; index1<=economyGraphList[read1][0].readId; index1++) // Explore all neighbours
						{
							read2=economyGraphList[read1][index1].readId;
							if(exploredReads[read2]==0) // Not explored
							{
								queue[end++]=read2; // Put in the queue
								totalEdgeInserted+=insertAllEdgesOfRead(read2);
							}
						}
						markTransitiveEdge(read1); // Mark transitive edges
					}
					if(exploredReads[read1]==2)
					{
						for(index1=1; index1<=economyGraphList[read1][0].readId; index1++) // Then explore all neighbour's neighbours
						{
							read2=economyGraphList[read1][index1].readId;
							if(exploredReads[read2]==1)
							{
								for(index2=1; index2<=economyGraphList[read2][0].readId; index2++) // Explore all neighbours neighbours
								{
									read3=economyGraphList[read2][index2].readId;
									if(exploredReads[read3]==0) // Not explored
									{
										queue[end++]=read3; // Put in the queue
										totalEdgeInserted+=insertAllEdgesOfRead(read3);
									}
								}
								markTransitiveEdge(read2); // Mark transitive edge
							}
						}
						transitiveEdgesRemoved+=removeTransitiveEdges(read1); // Remove the transitive edges
					}
				}
			}
		}
	}
	free(markedNodes);
	free(exploredReads);
	free(queue);

	logStream<< "     Total edges inserted: " << totalEdgeInserted << "\n";
	logStream<< "Total number of hash miss: " << numberOfHashMiss << "\n";
	logStream<< "  Transitive edge removed: " << transitiveEdgesRemoved << "\n";
	logStream<< "Function buildOverlapGraphEconomy() in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
}

/* ============================================================================================
   This function inserts all the edges of a read and marks it as explored. At the end it sorts
   the edges according to their length.
   ============================================================================================ */
int EconomyGraph::insertAllEdgesOfRead(uint64_t read)
{
	uint64_t j, k, read1=read, totalEdgeInserted=0;
	uint64_t *uint64Window;
	int64_t index;
	if(exploredReads[read1]!=0) // already explored
		return 0;
	else
		exploredReads[read1]=1; // mark as the node is explored
	Read *r1 = loaderObj->getRead(read1);
	uint16_t hashStringLength = hashObj->hashStringLength;
	for(j=0; j<=(uint16_t)(r1->length-hashStringLength); j++)
	{
		uint64Window=get64Bit2Int(r1->readInt, j, hashStringLength);
		index=hashObj->hashTableSearch(uint64Window, numberOfHashMiss);
		free((uint64_t*) uint64Window);
		if(index!=-1) // found in hash table
		{
			for(k=1; k<=hashObj->hashTableList[index][0].readId; k++)
			{
				uint64_t read2=hashObj->hashTableList[index][k].readId;
				int type=hashObj->hashTableList[index][k].type;
				Read *r2 = loaderObj->getRead(read2);
				int32_t ovlpLength=-1;
				int8_t typeOfEdgeToInsert=-1;
				if(exploredReads[read2])
					continue;
				if(type==0 && j<=(uint16_t)(r1->length-minOverlap) && compareStringInBytesPrevious(r1->readInt, r2->readInt, j, r1->length, r2->length))
				{
					ovlpLength=r2->length-(r1->length-j);
					typeOfEdgeToInsert=3;
				}
				else if(type==1 && j>=(uint16_t)(minOverlap-hashStringLength) && compareStringInBytesPrevious(r1->readReverseInt, r2->readReverseInt, r1->length-j-hashStringLength, r1->length, r2->length))
				{
					ovlpLength=r2->length-j-hashStringLength;
					typeOfEdgeToInsert=0;
				}
				else if(type==2 && j<=(uint16_t)(r1->length-minOverlap) && compareStringInBytesPrevious(r1->readInt, r2->readReverseInt, j, r1->length, r2->length))
				{
					ovlpLength=r2->length-(r1->length-j);
					typeOfEdgeToInsert=2;
				}
				else if(type==3 && j>=(uint16_t)(minOverlap-hashStringLength) && compareStringInBytesPrevious(r1->readReverseInt, r2->readInt, r1->length-j-hashStringLength, r1->length, r2->length))
				{
					ovlpLength=r2->length-j-hashStringLength;
					typeOfEdgeToInsert=1;
				}
				if(ovlpLength!=-1)
				{
					totalEdgeInserted+=insertEdgeEconomy(read1, read2, ovlpLength, typeOfEdgeToInsert);
				}
			}
		}
	}
	// sort in decreasing order of overlap length
	if(economyGraphList[read1]!=NULL && economyGraphList[read1][0].readId>1) // Read has some edges
		sort(economyGraphList[read1]+1, economyGraphList[read1]+economyGraphList[read1][0].readId+1, compareLengthBased); //sort according to length
	return totalEdgeInserted*2; // counting the twin edges
}

/* ============================================================================================
   Mark transitive edges using Gene Myers algorithm.
   ============================================================================================ */
void EconomyGraph::markTransitiveEdge(uint64_t from_ID)
{
	uint64_t i, j, a, b;
	uint8_t type1, type2;
	for(i=1; i<=economyGraphList[from_ID][0].readId; i++)
		markedNodes[economyGraphList[from_ID][i].readId]=1; // In play
	for(i=1; i<=economyGraphList[from_ID][0].readId; i++)
	{
		a=economyGraphList[from_ID][i].readId;
		if(markedNodes[a]==1) // If In play
		{
			for(j=1; j<=economyGraphList[a][0].readId; j++)
			{
				b=economyGraphList[a][j].readId;
				if(markedNodes[b]==1)
				{
					type1=economyGraphList[from_ID][i].type;
					type2=economyGraphList[a][j].type;
					if((type1==0 || type1==2) && (type2==0 || type2==1)) // check edge orientation
						markedNodes[b]=2; //Eliminated
					else if((type1==1||type1==3) && (type2==2 || type2==3)) //check edge orientation
						markedNodes[b]=2; //Eliminated
				}
			}
		}
	}
	for(i=1; i<=economyGraphList[from_ID][0].readId; i++)
	{
		if(markedNodes[economyGraphList[from_ID][i].readId]==2) //Eliminated
			economyGraphList[from_ID][i].mark=1; //marked for deletion
	}

	for(i=1; i<=economyGraphList[from_ID][0].readId; i++)
		markedNodes[economyGraphList[from_ID][i].readId]=0; // Mark as vacant.
	markedNodes[from_ID]=0; // Mark as vacant.
	exploredReads[from_ID]=2;
}

uint64_t EconomyGraph::removeTransitiveEdges(uint64_t read)
{
	uint64_t index, numberRemoved, newListSize=0;
	EconomyEdge *newList=NULL;
	for(index=1; index<=economyGraphList[read][0].readId; index++)
	{
		if(!(economyGraphList[read][index].mark))
		{
			if(newListSize==0) // First Edge of read
			{
				if((newList=(EconomyEdge*)malloc(2*sizeof(EconomyEdge)))==NULL)
					printError(MEM_ALLOC, "malloc newList");
			}
			else
			{
				if((newList=(EconomyEdge*)realloc(newList, (newListSize+2)*sizeof(EconomyEdge)))==NULL)
					printError(MEM_ALLOC, "realloc newList");
			}
			newList[++newListSize]=economyGraphList[read][index];
		}
	}
	numberRemoved=economyGraphList[read][0].readId-newListSize;
	free(economyGraphList[read]);
	economyGraphList[read]=newList;
	economyGraphList[read][0].readId=newListSize;
	return numberRemoved;
}

/* ============================================================================================
   Performs string comparison using 64 bit ints.
   ============================================================================================ */
int EconomyGraph::compareStringInBytes(uint8_t *read2bits1, uint8_t *read2bits2, uint16_t start, uint16_t readLength1, uint16_t readLength2, uint64_t read2ID)
{
	uint64_t *inta, *intb;
	uint16_t hashStringLength = hashObj->hashStringLength;
	uint16_t start1=start+hashStringLength, start2=hashStringLength, length;
	if(readLength2-start2<=readLength1-start1)
	{
		while(start2<readLength2)
		{
			length=min(64, readLength2-start2);
			inta=get64Bit2Int(read2bits1, start1, length);
			intb=get64Bit2Int(read2bits2, start2, length);
			if(inta[0]!=intb[0] || inta[1]!=intb[1] )
			{
				free((uint64_t *) inta);
				free((uint64_t *) intb);
				return 0;
			}
			free((uint64_t *) inta);
			free((uint64_t *) intb);
			start1+=length;
			start2+=length;
		}
		exploredReads[read2ID]=6; //this is a contained read.
		return 0; 
	}	
	else
	{
		while(start1<readLength1)
		{
			length=min(64, readLength1-start1);
			inta=get64Bit2Int(read2bits1, start1, length);
			intb=get64Bit2Int(read2bits2, start2, length);
			if(inta[0]!=intb[0] || inta[1]!=intb[1] )
			{
				free((uint64_t *) inta);
				free((uint64_t *) intb);
				return 0;
			}
			free((uint64_t *) inta);
			free((uint64_t *) intb);
			start1+=length;
			start2+=length;
		}
		return 1;
	}
}

/* ============================================================================================
   Performs string comparison using 64 bit ints.
   ============================================================================================ */
int EconomyGraph::compareStringInBytesPrevious(uint8_t *read2bits1, uint8_t *read2bits2, uint16_t start, uint16_t readLength1, uint16_t readLength2)
{
	uint64_t *inta, *intb;
	uint16_t hashStringLength = hashObj->hashStringLength;
	uint16_t start1=start+hashStringLength, start2=hashStringLength, length;
	if(readLength2-start2<=readLength1-start1)
	{
		while(start2<readLength2)
		{
			length=min(64, readLength2-start2);
			inta=get64Bit2Int(read2bits1, start1, length);
			intb=get64Bit2Int(read2bits2, start2, length);
			if(inta[0]!=intb[0] || inta[1]!=intb[1] )
			{
				free((uint64_t *) inta);
				free((uint64_t *) intb);
				return 0;
			}
			free((uint64_t *) inta);
			free((uint64_t *) intb);
			start1+=length;
			start2+=length;
		}
		return 1; 
	}	
	else
	{
		while(start1<readLength1)
		{
			length=min(64, readLength1-start1);
			inta=get64Bit2Int(read2bits1, start1, length);
			intb=get64Bit2Int(read2bits2, start2, length);
			if(inta[0]!=intb[0] || inta[1]!=intb[1] )
			{
				free((uint64_t *) inta);
				free((uint64_t *) intb);
				return 0;
			}
			free((uint64_t *) inta);
			free((uint64_t *) intb);
			start1+=length;
			start2+=length;
		}
		return 1;
	}
}

/* ============================================================================================
   This function inserts edges in the economy graph.
   ============================================================================================ */
int EconomyGraph::insertEdgeEconomy(uint64_t vertex_u, uint64_t vertex_v, uint32_t delta, uint8_t type)
{
	if(vertex_u==vertex_v) //do not want to add edge (a, a)
		return 0;

	uint8_t typeForward64=type, typeReverse64=reverseEdgeType(type);
	
	EconomyEdge edgeuv(vertex_v , typeForward64, 0, delta);
	int delta2 = loaderObj->getRead(vertex_u)->length-(loaderObj->getRead(vertex_v)->length-delta);
	EconomyEdge edgevu(vertex_u, typeReverse64, 0, delta2);
	if(economyGraphList[vertex_u]==NULL) // First Edge of vertex_u
	{
		if((economyGraphList[vertex_u]=(EconomyEdge*)malloc(2*sizeof(EconomyEdge)))==NULL)
			printError(MEM_ALLOC, "graphEconomy[vertex_u]");
		economyGraphList[vertex_u][0].readId = 0;
	}
	else
	{
		if((economyGraphList[vertex_u]=(EconomyEdge*)realloc(economyGraphList[vertex_u], (economyGraphList[vertex_u][0].readId+2)*sizeof(EconomyEdge)))==NULL)
			printError(MEM_ALLOC, "graphEconomy[vertex_u]");
	}

	if(economyGraphList[vertex_v]==NULL) // First Edge of vertex_v
	{
		if((economyGraphList[vertex_v]=(EconomyEdge*)malloc(2*sizeof(EconomyEdge)))==NULL)
			printError(MEM_ALLOC, "graphEconomy[vertex_v]");
		economyGraphList[vertex_v][0].readId = 0;
	}
	else
	{
		if((economyGraphList[vertex_v]=(EconomyEdge*)realloc(economyGraphList[vertex_v], (economyGraphList[vertex_v][0].readId+2)*sizeof(EconomyEdge)))==NULL)
			printError(MEM_ALLOC, "graphEconomy[vertex_v]");
	}
	economyGraphList[vertex_u][++economyGraphList[vertex_u][0].readId] = edgeuv;
	economyGraphList[vertex_v][++economyGraphList[vertex_v][0].readId] = edgevu;
	return 1;
}

//for comparison
//first length, then id, then type
bool compareLengthBased(const EconomyEdge &first, const EconomyEdge &second)
{
	if(first.length > second.length)
		return true;
	if(first.length < second.length)
		return false;
	//first.length == second.length
	if(first.readId > second.readId)
		return true;
	if(first.readId < second.readId)
		return false;
	//first.length == second.length && first.destination == second.destination
	if(first.type > second.type)
		return true;
	if(first.type < second.type)
		return false;
	//first.length == second.length && first.destination == second.destination && first.type == second.type
	return false;
}

//for comparison
//first id, then type, then length
bool compareIdBased(const EconomyEdge &first, const EconomyEdge &second)
{
	if(first.readId < second.readId)
		return true;
	if(first.readId > second.readId)
		return false;
	//first.destination == second.destination
	if(first.type < second.type)
		return true;
	if(first.type > second.type)
		return false;
	//first.destination == second.destination && first.type == second.type
	if(first.length < second.length)
		return true;
	if(first.length > second.length)
		return false;
	//first.destination == second.destination && first.type == second.type && first.length == second.length
	return false;
}


void EconomyGraph::sortEconomyGraph()
{
	uint64_t i;
	logStream<< "\nIn function sortEdgesEconomy().\n";
	logStream.flush();
	time_t second_s=time(NULL);
	#pragma omp parallel for
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
	{
		if(economyGraphList[i]!=NULL)
		{
			if(economyGraphList[i][0].readId>1)
				sort(economyGraphList[i]+1, economyGraphList[i]+economyGraphList[i][0].readId+1, compareIdBased); //sort according to id
		}
	}
	logStream<< "Function sortEdgesEconomy() in " << time(NULL)-second_s << " sec.\n";
	logStream.flush();
}

