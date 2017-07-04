/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 ****************************************************************************/

#include "matePair.h"

extern ofstream logStream;
extern uint64_t averageReadLength;

MatePair::MatePair(OverlapGraph *graph1, ReadLoader *loader1)
{
	graphObj = graph1;
	loaderObj = loader1;
	numberOfLibrary = 0;
	readToEdgeList = NULL;
	minimumUpperBoundOfInsert = 1000000;
	maximumUpperBoundOfInsert = 0;
	Mean = NULL;
	standardDeviation = NULL;
	lowerBoundOfInsert = NULL;
	upperBoundOfInsert = NULL;
	uint64_t i;
	if((matePairList=(MatePairInfo**)malloc((loaderObj->numberOfUniqueReads+1)*sizeof(MatePairInfo*)))==NULL)
		printError(MEM_ALLOC, "matePairList");
	#pragma omp parallel for
	for(i=0; i<=loaderObj->numberOfUniqueReads; i++)
		matePairList[i]=NULL;
}

MatePair::~MatePair()
{
	uint64_t i;
	MatePairInfo *w, *w_next;
	#pragma omp parallel for private(i, w, w_next)
	for(i=0; i<=loaderObj->numberOfUniqueReads; i++)
	{
		for(w=matePairList[i]; w!=NULL; w=w_next)
		{
			w_next = w->next;
			free(w);
		}
	}
	free(matePairList);

	free(Mean);
	free(standardDeviation);
	free(upperBoundOfInsert);
	free(lowerBoundOfInsert);

	ReadToEdgeMap *edg, *edg_next;
	#pragma omp parallel for private(i, edg, edg_next)
	for(i=0; i<=loaderObj->numberOfUniqueReads; i++)
	{
		if(readToEdgeList[i]!=NULL) // Nodes
		{
			for(edg=readToEdgeList[i]; edg!=NULL; edg=edg_next)
			{
				edg_next=edg->next;
				free(edg->locationForward);
				free(edg->locationReverse);
				free(edg);
			}
		}
	}
	free(readToEdgeList);
}

void MatePair::mapMatePairsFromList(string listPath)
{
	logStream<< "\nIn function mapMatePairsFromList().\n";
	logStream.flush();
	time_t seconds_s=time(NULL);

	string line;
	ifstream fin(listPath.c_str());
	if(fin.is_open()==false)
		printError(OPEN_FILE, listPath);
	string var, val, var1, val1;
	uint32_t mateFile=0;
	size_t pos;
	while(getline(fin, line)>0)
	{
		if(line.empty() || line[0]=='#')
			continue;
		pos = line.find('=');
		if(pos==line.npos)
		{
			cout<< "[ERROR] List of input files in a wrong format!\n\n";
			exit(0);
		}
		var = trim(line.substr(0, pos));
		val = trim(line.substr(pos+1));
		if(mateFile%2==0 && var=="f1")
		{
			var1 = var;
			val1 = val;
		}
		else if(mateFile%2==0 && var=="f")
		{
			mapMatePairs(val, "", (mateFile/2)+1);
			mateFile++;
		}
		else if(mateFile%2==1 && var=="f2")
		{
			mapMatePairs(val1, val, (mateFile/2)+1);
		}
		else
		{
			cout<< "[ERROR] List of input files in a wrong format!\n\n";
			exit(0);
		}
		mateFile++;
	}

	logStream<< "\nNumber of datasets: " << mateFile/2 << "\n";
	logStream<< "Function mapMatePairsFromList() in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
}

/* ============================================================================================
   This function maps pair of reads. It will load the reads from file again.
   ============================================================================================ */
void MatePair::mapMatePairs(string mateFile1, string mateFile2, int library)
{
	logStream<< "\nIn function mapMatePairs().\n";
	time_t seconds_s=time(NULL);
	uint64_t readInFile=0;
	deque<string> readsArray;
	int64_t arrayMem=0;
	InputReader *readerObj = new InputReader(mateFile1, mateFile2);
	logStream << "\tReading from file: " << mateFile1 << " " << (mateFile2!=""? "& "+mateFile2 : "" ) << endl;
	while(readerObj->getNextRead(readInFile))
	{
		readsArray.push_back(readerObj->read.sequence);
		readInFile++;
		arrayMem += sizeof(readerObj->read.sequence)+readerObj->read.sequence.capacity()*sizeof(char);
		if(readInFile%1000000==0)
		{
			arrayMem += sizeof(readsArray);
			processMatePairs(readsArray, library);
			readsArray.clear();
			arrayMem=0;
		}
	}
	if(readsArray.size()!=0)
	{
		arrayMem += sizeof(readsArray);
		processMatePairs(readsArray, library);
		readsArray.clear();
	}

	delete readerObj;
	numberOfLibrary=library;
	logStream << "\tReads in file: " << readInFile << "\n" 
			  << "Function mapMatePairs() in " << time(NULL)-seconds_s << " sec." << endl;

}

void MatePair::processMatePairs(deque<string> &readsArray, int library)
{
	uint64_t i;
	string read1, read2;
	int64_t id1, id2;
	bool type1, type2, flag_insert;
	MatePairInfo *u, *v, *w;
	#pragma omp parallel default(none) shared(readsArray,library)\
	private(i,read1,read2,id1,id2,type1,type2,flag_insert,w,u,v)
	{
		#pragma omp for
		for(i=0; i<readsArray.size()-1; i+=2)
		{
			read1 = readsArray[i];
			read2 = readsArray[i+1];
			if(isGoodRead(read1, loaderObj->minOverlap) && isGoodRead(read2, loaderObj->minOverlap))
			{
				id1=loaderObj->getIdOfRead(read1);
				id2=loaderObj->getIdOfRead(read2);
				if(id1>0)
					type1=1;
				else
					type1=0;
				if(id2>0)
					type2=1;
				else
					type2=0;
				id1=llabs(id1);
				id2=llabs(id2);
				flag_insert=1;
				for(w=matePairList[id1]; w!=NULL; w=w->next)  
				{
					if(w->ID==(uint64_t)id2 && w->type1==type1 && w->type2==type2 && w->library==library)
					{
						flag_insert=0;
						w->freq++;
						break;
					}
				}
				if(flag_insert==1) // Insert if not already present
				{
					if((u = (MatePairInfo*)malloc(sizeof(MatePairInfo)))==NULL)
						printError(MEM_ALLOC, "u");
					u->ID=id2;
					u->freq=1;
					u->flag=1;
					u->type1=type1;
					u->type2=type2;
					u->library=library;
					u->next=matePairList[id1];
					matePairList[id1]=u;
				}
				flag_insert=1;
				for(w=matePairList[id2]; w!=NULL; w=w->next) 
				{
					if(w->ID==(uint64_t)id1 && w->type1==type2 && w->type2==type1 && w->library==library)
					{
						flag_insert=0;
						w->freq++;
						break;
					}
				}
				if(flag_insert==1) // Insert if not already present.
				{
					if((v = (MatePairInfo*)malloc(sizeof(MatePairInfo)))==NULL)
						printError(MEM_ALLOC, "v");
					v->ID=id1;
					v->freq=1;
					v->flag=1;
					v->type1=type2;
					v->type2=type1;
					v->library=library;
					v->next=matePairList[id2];
					matePairList[id2]=v;
				}
			}
		}
	}
}

/* ============================================================================================
   This function will estimate mean and standard deviation from matepairs on the same edge
   ============================================================================================ */
void MatePair::meanSdEstimation()
{
	logStream<< "\nIn function meanSdEstimation().\n";
	logStream.flush();
	time_t second_s=time(NULL);
	int counter, i;
	if((Mean=(uint*)malloc((numberOfLibrary+1)*sizeof(uint)))==NULL)
		printError(MEM_ALLOC, "Mean");
	if((standardDeviation=(uint*)malloc((numberOfLibrary+1)*sizeof(uint)))==NULL)
		printError(MEM_ALLOC, "standardDeviation");
	if((upperBoundOfInsert=(int*)malloc((numberOfLibrary+1)*sizeof(int)))==NULL)
		printError(MEM_ALLOC, "upperBoundOfInsert");
	if((lowerBoundOfInsert=(int*)malloc((numberOfLibrary+1)*sizeof(int)))==NULL)
		printError(MEM_ALLOC, "lowerBoundOfInsert");
	for(i=0; i<=numberOfLibrary; i++)
	{
		Mean[i]=0;
		standardDeviation[i]=0;
	}
	mapReadsToEdges();
	mapReadLocations();
	for(counter=1; counter<=numberOfLibrary; counter++)
	{
		int mu=5000, sd=5000, returnMu, returnSD; // Initial mean set to 5000
		logStream<< "\n\tComputing mean and SD of library " << counter << "\n";
		for(i=1;i<=10;i++) 
		{
			logStream<< "\tEstimation " << i << ": ";
			computeMeanSD(mu,sd,counter,&returnMu,&returnSD);
			if(abs(mu-returnMu)<=returnMu/100  && abs(sd-returnSD)<=returnSD/100)
			{
				mu=returnMu;
				sd=returnSD;
				logStream<< "Mean: " << mu + averageReadLength << " \tSD: " << sd << " (Final estimation)\n";

				break;
			}
			else
			{
				mu=returnMu;
				sd=returnSD;
				logStream<< "Mean: " << mu + averageReadLength << " \tSD: " << sd << "\n";

				}
		}
		Mean[counter]=returnMu + averageReadLength;
		standardDeviation[counter]=returnSD;
	}
	minimumUpperBoundOfInsert=1000000;
	maximumUpperBoundOfInsert=0;
	for(i=1; i<=numberOfLibrary; i++)
	{
		lowerBoundOfInsert[i]=Mean[i]-3*standardDeviation[i];
		if(lowerBoundOfInsert[i]<0) 
			lowerBoundOfInsert[i]=0;
		upperBoundOfInsert[i]=Mean[i]+3*standardDeviation[i];

		if(minimumUpperBoundOfInsert>=(uint)upperBoundOfInsert[i])
			minimumUpperBoundOfInsert=upperBoundOfInsert[i];
		if(maximumUpperBoundOfInsert<=(uint)upperBoundOfInsert[i])
			maximumUpperBoundOfInsert=upperBoundOfInsert[i];

		logStream<< "\n\tLibaray " << i << ": Mean: " << Mean[i] << " \tSD: " << standardDeviation[i] << " \tLB: " << lowerBoundOfInsert[i] << " \tUB: " << upperBoundOfInsert[i];
	}
	maximumUpperBoundOfInsert=(3*(int)ceilf(averageReadLength/100.0))*maximumUpperBoundOfInsert;
	logStream<< "\n\n\tMinimum upper bound: " << minimumUpperBoundOfInsert << "\n"
			<< "\tMaximum upper bound: " << maximumUpperBoundOfInsert << "\n"
			<< "Function meanSdEstimation() in " << time(NULL)-second_s <<" sec.\n";
	logStream.flush();
}

/* ============================================================================================
   Next estimation of mean and standard deviation from the previous estimation.
   ============================================================================================ */
void MatePair::computeMeanSD(int mu, int sd, int library, int *returnMu, int *returnSD)
{
	uint64_t i, mate_pair_1, mate_pair_2, distance;
	int32_t *distance1, *distance2;
	int64_t sum=0, counter=0;
	long double sqrd_error=0;
	MatePairInfo *a;
	ReadToEdgeMap *u, *v;
	int64_t sum_local, counter_local;
	long double sqrd_error_local;
	#pragma omp parallel default(none) shared(mu,library,sum,counter,sqrd_error)\
	private(i,mate_pair_1,a,mate_pair_2,u,v,distance1,distance2,distance,sum_local,counter_local,sqrd_error_local)
	{
		sum_local=0;
		counter_local=0;
		sqrd_error_local=0;
		#pragma omp for
		for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
		{
			mate_pair_1=i;
			for(a=matePairList[mate_pair_1]; a!=NULL; a=a->next)
			{
				mate_pair_2=a->ID;
				if(a->flag==0 && mate_pair_1<mate_pair_2 && a->library==library)
				{
					for(u=readToEdgeList[mate_pair_1]; u!=NULL; u=u->next)
					{
						for(v=readToEdgeList[mate_pair_2]; v!=NULL; v=v->next)
						{
							if(u->edge==v->edge || u->edge==v->edge->twinEdge) // mate_pairs on same edge
							{
								distance1 = findDistanceOnEdge(v->edge, mate_pair_1);
								distance2 = findDistanceOnEdge(v->edge, mate_pair_2);
								distance=llabs(llabs(distance1[1])-llabs(distance2[1]));
								if(distance1[0]==1 && distance2[0]==1 && distance<(uint64_t)(4*mu))
								{
									sqrd_error_local+=(mu-distance)*(mu-distance);
									counter_local++;
									sum_local+=distance;
									free(distance1);
									free(distance2);
									break;
								}
								free(distance1);
								free(distance2);
							}
						}
					}
				}
			}
		}
		#pragma omp critical (update_meanSD)
		{
			sum += sum_local;
			counter += counter_local;
			sqrd_error += sqrd_error_local;
		}
	}
	logStream<< "Mate-pairs considered: " << counter << " \t";
	mu=sum/counter;
	sd=sqrt(sqrd_error/(counter-1));
	*returnMu=mu;
	*returnSD=sd;
}


/* ===========================================================================================
   This function will remove unnecessary mate pairs (Mate pairs that are on the same edge).
   ============================================================================================ */
void MatePair::mapReadsToEdges()
{
	logStream<< "\nIn function mapReadsToEdges().\n";
	logStream.flush();
	time_t seconds_s=time(NULL);
	uint64_t i, flag, index;
	if(readToEdgeList==NULL)
	{
		if((readToEdgeList = (ReadToEdgeMap**)malloc((loaderObj->numberOfUniqueReads+1)*sizeof(ReadToEdgeMap*)))==NULL)
			printError(MEM_ALLOC, "Read to Edge mapping");
		#pragma omp parallel for default(none) private(i)
		for(i=0; i<=loaderObj->numberOfUniqueReads; i++)
			readToEdgeList[i]=NULL;
	}
	Edge *v;
	MatePairInfo *x;
	ReadToEdgeMap *edg, *edg_temp;
	#pragma omp parallel default(none) private(i,edg,edg_temp)
	{
		#pragma omp for
		for(i=0; i<=loaderObj->numberOfUniqueReads; i++) // delete previous mapping of read to edge
		{
			if(readToEdgeList[i]!=NULL) // Nodes
			{
				for(edg=readToEdgeList[i]; edg!=NULL; edg=edg_temp)
				{
					edg_temp=edg->next;
					free(edg->locationForward);
					free(edg->locationReverse);
					free(edg);
				}
			}
			readToEdgeList[i]=NULL;
		}
	}
	#pragma omp parallel default(none) private(i,v,index,edg,edg_temp,flag)
	{
		#pragma omp for
		for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
		{
			if(graphObj->graph[i]!=NULL)
			{
				for(v=graphObj->graph[i]; v!=NULL; v=v->next)
				{
					if(i<=v->ID)
					{
						for(index=1; index<=v->listOfReads[0].readId; index++)
						{
							flag=1;
							for(edg_temp=readToEdgeList[v->listOfReads[index].readId]; edg_temp!=NULL; edg_temp=edg_temp->next)
							{
								if(edg_temp->edge==v || edg_temp->edge==v->twinEdge)
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
								edg->next=readToEdgeList[v->listOfReads[index].readId];
								readToEdgeList[v->listOfReads[index].readId]=edg;
							}
						}
					}
				}
			}
		}
	}
	ReadToEdgeMap *a,*b;
	#pragma omp parallel default(none) private(i,x,a,b)
	{
		#pragma omp for
		for(i=1; i<=loaderObj->numberOfUniqueReads; i++) // for each read i
		{
			for(x=matePairList[i]; x!=NULL; x=x->next) // for each matepair x of read i
			{
				if(x->flag==1)
				{
					for(a=readToEdgeList[i]; a!=NULL; a=a->next)
					{
						for(b=readToEdgeList[x->ID]; b!=NULL; b=b->next)
						{
							if(a->edge==b->edge || a->edge==b->edge->twinEdge) //both same edge
								x->flag=0; //we do not delete it, just mark it
						}
					}
				}
			}
		}
	}
	logStream<< "Function mapReadsToEdges() in " << time(NULL)-seconds_s <<  " sec.\n";
	logStream.flush();
}


/* ============================================================================================
   This function will map each read on every edge and location of the read on the edges.
   ============================================================================================ */
void MatePair::mapReadLocations()
{
	logStream<< "\nIn function mapReadLocations().\n";
	logStream.flush();
	time_t seconds_s=time(NULL);
	uint64_t i, index;
	uint32_t dist;
	ReadToEdgeMap *edge_ptr;
	Edge *v;
	#pragma omp parallel default(none) private(i,v,dist,index,edge_ptr)
	{
		#pragma omp for
		for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
		{
			for(v=graphObj->graph[i]; v!=NULL; v=v->next)
			{
				if(v->listOfReads[0].readId==0)//if simple edge
					continue;
				dist=0;
				for(index=1; index<=v->listOfReads[0].readId; index++)
				{
					dist+=v->listOfReads[index].distPrevious;
					for(edge_ptr=readToEdgeList[v->listOfReads[index].readId]; edge_ptr!=NULL; edge_ptr=edge_ptr->next)
					{
						if(edge_ptr->edge==v)
						{
							if(edge_ptr->locationForward==NULL) 
							{
								if((edge_ptr->locationForward=(uint32_t*)malloc((2)*sizeof(uint32_t)))==NULL)
									printError(MEM_ALLOC, "locationForward");
								edge_ptr->locationForward[0]=1;
								if(v->listOfReads[index].orientation)
									edge_ptr->locationForward[edge_ptr->locationForward[0]]=dist;
								else
									edge_ptr->locationForward[edge_ptr->locationForward[0]]=-dist;
							}
							else 
							{
								if((edge_ptr->locationForward=(uint32_t*)realloc(edge_ptr->locationForward, (edge_ptr->locationForward[0]+2)*sizeof(uint32_t)))==NULL)
									printError(MEM_ALLOC, "Reallocating locationForward");
								edge_ptr->locationForward[0]++;
								if(v->listOfReads[index].orientation)
									edge_ptr->locationForward[edge_ptr->locationForward[0]]=dist;
								else
									edge_ptr->locationForward[edge_ptr->locationForward[0]]=-dist;
							}
							break;
						}
						else if(edge_ptr->edge==v->twinEdge)
						{
							if(edge_ptr->locationReverse==NULL)
							{
								if((edge_ptr->locationReverse=(uint32_t*)malloc((2)*sizeof(uint32_t)))==NULL)
									printError(MEM_ALLOC, "locationReverse");
								edge_ptr->locationReverse[0]=1;
								if(v->listOfReads[index].orientation) 
									edge_ptr->locationReverse[edge_ptr->locationReverse[0]]=dist;
								else
									edge_ptr->locationReverse[edge_ptr->locationReverse[0]]=-dist;
							}
							else 
							{
								if((edge_ptr->locationReverse=(uint32_t*)realloc(edge_ptr->locationReverse, (edge_ptr->locationReverse[0]+2)*sizeof(uint32_t)))==NULL)
									printError(MEM_ALLOC, "Reallocating locationReverse");
								edge_ptr->locationReverse[0]++;
								if(v->listOfReads[index].orientation)
									edge_ptr->locationReverse[edge_ptr->locationReverse[0]]=dist;
								else
									edge_ptr->locationReverse[edge_ptr->locationReverse[0]]=-dist;
							}
							break;
						}
					}
				}
			}
		}
	}
	logStream<< "Function mapReadLocations() in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
}


/* ============================================================================================
   This function finds distance on edge.
   ============================================================================================ */
int32_t* MatePair::findDistanceOnEdge(Edge *edge, uint64_t read)
{
	ReadToEdgeMap *edge_ptr;
	uint64_t i;
	int32_t *distance;
	for(edge_ptr=readToEdgeList[llabs(read)]; edge_ptr!=NULL; edge_ptr=edge_ptr->next)
	{
		if(edge_ptr->edge==edge)
		{
			if(edge_ptr->locationForward!=NULL)
			{
				if((distance=(int32_t*)malloc((edge_ptr->locationForward[0]+1)*sizeof(int32_t)))==NULL)
					printError(MEM_ALLOC, "distance");
				for(i=0; i<=edge_ptr->locationForward[0]; i++)
					distance[i]=edge_ptr->locationForward[i];
				return distance;
			}
		}
		else if(edge_ptr->edge==edge->twinEdge)
		{
			if(edge_ptr->locationReverse!=NULL)
			{
				if((distance=(int32_t*)malloc((edge_ptr->locationReverse[0]+1)*sizeof(int32_t)))==NULL)
					printError(MEM_ALLOC, "distance");
				for(i=0; i<=edge_ptr->locationReverse[0]; i++)
					distance[i]=edge_ptr->locationReverse[i];
				return distance;
			}
		}
	}
	if((distance=(int32_t*)malloc(1*sizeof(int32_t)))==NULL)
		printError(MEM_ALLOC, "distance"); // Not found.
	distance[0]=0;
	return distance;
}
