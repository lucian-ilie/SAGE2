/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 ****************************************************************************/

#include "copycountEstimator.h"

extern ofstream logStream;
extern bool saveAll;

CopyCountEstimator::CopyCountEstimator(OverlapGraph *graph1, ReadLoader *loader1)
{
	graphObj = graph1;
	loaderObj = loader1;
	genomeSize=0;
	aStatisticsThreshold = 3;
	minDelta = 1000;
}

CopyCountEstimator::~CopyCountEstimator()
{
	
}

/* ============================================================================================
   This function estimates the genome size.
   ============================================================================================ */
uint64_t CopyCountEstimator::genomeSizeEstimation()
{
	logStream<< "\nIn function genomeSizeEstimation().\n";
	logStream.flush();
	time_t seconds_s=time(NULL);
	uint64_t currentGenomeSize, previousGenomeSize=0;
	uint16_t i;
	for(i=1;i<=10;i++) 
	{
		currentGenomeSize=findGenomeSize(previousGenomeSize);
		if(previousGenomeSize==currentGenomeSize) 
		{
			logStream<< "\tGenome Size: " << currentGenomeSize << " (Final estimation)\n";
			break;
		}
		else
		{
			previousGenomeSize=currentGenomeSize;
			logStream<< "\tGenome Size: " << currentGenomeSize << " (Estimation " << i << ")\n";
		}
	}
	logStream<< "Function genomeSizeEstimation in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
	genomeSize = currentGenomeSize;
	return currentGenomeSize;
}

/* ============================================================================================
   Find next estimation from previous estimation.
   ============================================================================================ */
uint64_t CopyCountEstimator::findGenomeSize(uint64_t previousEstimation)
{
	uint64_t i, delta, k, deltaSum=0, kSum=0, currentGenomeSize, index;
	float a_statistic;
	Edge *v;
	Read *r1;
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
	{
		for(v=graphObj->graph[i]; v!=NULL; v=v->next)
		{
			if(i<v->ID)
			{
				delta=v->lengthOfEdge-1;
				k=0;
				if(v->listOfReads[0].readId!=0)
				{
					for(index=1; index<=v->listOfReads[0].readId; index++)
					{
						r1 = loaderObj->getRead(v->listOfReads[index].readId);
						k += r1->frequency;
					}
				}
				if(previousEstimation!=0)
				{
					a_statistic=(float)((float)delta*(float)((float)loaderObj->numberOfReads/(float)previousEstimation)-(float)((float)k*log((float)2)));
					if(a_statistic>=aStatisticsThreshold && delta>=minDelta)
					{
						deltaSum += delta;
						kSum += k;
					}
				}
				else if(v->lengthOfEdge>500)  // Only take edges longer than 500 for the first estimation
				{
					deltaSum += delta;
					kSum += k;
				}
			}
		}
	}
	currentGenomeSize=(uint64_t)((float)loaderObj->numberOfReads/(float)kSum*(float)deltaSum);
	return currentGenomeSize;
}

/* ============================================================================================
   This function generates flow in the graph.
   ============================================================================================ */
void CopyCountEstimator::computeMinCostFlow(string outDir)
{
	logStream<< "\nIn function computeMinCostFlow().\n";
	logStream.flush();
	time_t seconds_s=time(NULL), seconds_start_CS2;
	uint64_t edgesUsedMoreThanOnce=0, simple_edge=0, edgesUsedExactlyOnce=0, nodes_in_graph=0, edges_in_graph=0;
	uint64_t numberOfEdges=0, numberOfNodes=0, cnt, sinkToSourceCap;
	uint64_t *node_index, *inverse_node_index, *from, *to, *flow, i, j, k, delta, index;
	int xi, cost_scale_factor=10, flow_lb, lb_1=0, lb_2=0, lb_3=0, ub_1=0, ub_2=0, ub_3=0;
	int lb=0, ub=100, first_UB=3, second_UB=5, third_UB=1000;
	uint64_t n_1in, n_1out, n_2in, n_2out,node_1=0, node_2=0, node_3=0, node_4=0;
	uint64_t u_1in=0, u_1out=0, u_2in=0, u_2out=0, v_1in=0, v_1out=0, v_2in=0, v_2out=0;
	float cost, cost_1, cost_2, cost_3, a_statistics[25];
	double cost_of_cs2;
	Edge *v;
	Read *r1;
	string cs2InputFileName = outDir+"input.cs2";
	string cs2OutputFileName = outDir+"output.cs2";

	for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
	{
		if(graphObj->graph[i]!=NULL)
			nodes_in_graph++;
		for(v=graphObj->graph[i]; v!=NULL; v=v->next)
		{
			if(i>=v->ID)
			{
				edges_in_graph++;
				if(v->listOfReads[0].readId==0)
					simple_edge++;
				else
				{
					delta=v->lengthOfEdge-1;
					k=0;
					for(index=1; index<=v->listOfReads[0].readId; index++)
					{
						r1 = loaderObj->getRead(v->listOfReads[index].readId);
						k += r1->frequency;
					}
					for(j=1;j<=20;j++)
						a_statistics[j]=(float)((float)delta*(float)((float)loaderObj->numberOfReads/(float)genomeSize)-(float)((float)k*log((float)((float)(j+1)/(float)j))));
					if(a_statistics[1]>=aStatisticsThreshold && delta>=minDelta)
						edgesUsedExactlyOnce++;
					else
						edgesUsedMoreThanOnce++;
				}
			}
		}
	}
	logStream<<  "\t          Number of nodes: " << nodes_in_graph << "\n"
			<< "\t          Number of edges: " << edges_in_graph << "\n"
			<<  "\t             Simple edges: " << simple_edge << "\n"
			<<  "\t  Edges used exactly once: " << edgesUsedExactlyOnce << "\n"
			<<  "\tEdges used more than once: " << edgesUsedMoreThanOnce << "\n\n";

	numberOfNodes=nodes_in_graph*4+2;
	numberOfEdges=1+6*nodes_in_graph+2*(edges_in_graph-(edgesUsedMoreThanOnce))+6*(edgesUsedMoreThanOnce);

	if((from = (uint64_t *)malloc(numberOfEdges*sizeof(uint64_t)))==NULL)
		printError(MEM_ALLOC, "from");
	if((to = (uint64_t *)malloc(numberOfEdges*sizeof(uint64_t)))==NULL)
		printError(MEM_ALLOC, "to");
	if((flow = (uint64_t *)malloc(numberOfEdges*sizeof(uint64_t)))==NULL)
		printError(MEM_ALLOC, "flow");

	if((node_index=(uint64_t *)malloc((loaderObj->numberOfUniqueReads+1)*sizeof(uint64_t)))==NULL)
		printError(MEM_ALLOC, "node index");

	#pragma omp parallel for
	for(i=0; i<=loaderObj->numberOfUniqueReads; i++)
		node_index[i]=0;
	if((inverse_node_index=(uint64_t *)malloc((numberOfNodes+4)*sizeof(uint64_t)))==NULL)
		printError(MEM_ALLOC, "inverse_node index");

	#pragma omp parallel for
	for(i=0;i<numberOfNodes+4;i++)
		inverse_node_index[i]=0;

	sinkToSourceCap = numberOfNodes * 100;

	ofstream fp_overlap(cs2InputFileName.c_str());
	if(fp_overlap.is_open()==false)
		printError(OPEN_FILE, cs2InputFileName);

	fp_overlap<< std::fixed;
	fp_overlap<< "p min " << numberOfNodes << " " << numberOfEdges << "\n";
	fp_overlap<< "n 1 0\n";
	fp_overlap<< "n 2 0\n";
	fp_overlap<< "a 2 1 1 " << sinkToSourceCap << " 0\n";
	cnt=3;
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++) // edges for nodes
	{
		if(graphObj->graph[i]!=NULL)
		{
			n_1in=cnt;
			n_1out=cnt+1;
			n_2in=cnt+2;
			n_2out=cnt+3;

			node_index[i]=cnt;
			inverse_node_index[cnt]=i;
			inverse_node_index[cnt+1]=i;
			inverse_node_index[cnt+2]=i;
			inverse_node_index[cnt+3]=i;

			lb=0;
			ub=100;
			cost=1000000;
			fp_overlap<< "a 1 " << n_1in << " " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";//source to node with lower bound 0 upper bound 100 and cost 1000000
			fp_overlap<< "a " << n_1out << " 2 " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";//node to sink with lower bound 0 upper bound 100 and cost 1000000
			fp_overlap<< "a 1 " << n_2in << " " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";//source to node with lower bound 0 upper bound 100 and cost 1000000
			fp_overlap<< "a " << n_2out << " 2 " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";//node to sink with lower bound 0 upper bound 100 and cost 1000000
			lb=0;
			ub=100;
			cost=1;// every read can be used any number of times for real reads. should be used at least once for corrected reads.
			fp_overlap<< "a " << n_1in << " " << n_1out << " " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";//node_in to node_out with lower bound 0 upper bound 100 and cost 1
			fp_overlap<< "a " << n_2in << " " << n_2out << " " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";//node_in to node_out with lower bound 0 upper bound 100 and cost 1
			cnt=cnt+4;
		}
	}
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++) // edges for edges
	{
		for(v=graphObj->graph[i]; v!=NULL; v=v->next)
		{
			if(i>=v->ID)
			{
				j=node_index[v->fromID];
				u_1in=j;
				u_1out=j+1;
				u_2in=j+2;
				u_2out=j+3;
				j=node_index[v->ID];
				v_1in=j;
				v_1out=j+1;
				v_2in=j+2;
				v_2out=j+3;
				delta=v->lengthOfEdge-1;
				k=0;
				if(v->listOfReads[0].readId==0) // these are simple edges
				{
					k=0;
					a_statistics[1]=0;
				}
				else	//composite edges
				{
					for(index=1; index<=v->listOfReads[0].readId; index++)
					{
						r1 = loaderObj->getRead(v->listOfReads[index].readId);
						k += r1->frequency;
					}
					for(j=1; j<=20; j++)
						a_statistics[j]=(float)((float)delta*(float)((float)loaderObj->numberOfReads/(float)genomeSize)-(float)((float)k*log((float)((float)(j+1)/(float)j))));
				}
				if(v->typeOfEdge==0)
				{
					node_1=v_1out;
					node_2=u_1in;
					node_3=u_2out;
					node_4=v_2in;
				}
				else if(v->typeOfEdge==1)
				{
					node_1=u_2out;
					node_2=v_1in;
					node_3=v_2out;
					node_4=u_1in;
				}
				else if(v->typeOfEdge==2)
				{
					node_1=u_1out;
					node_2=v_2in;
					node_3=v_1out;
					node_4=u_2in;
				}
				else if(v->typeOfEdge==3)
				{
					node_1=u_1out;
					node_2=v_1in;
					node_3=v_2out;
					node_4=u_2in;
				}
				if(v->listOfReads[0].readId==0) // simple edge: could be used any number of times
				{
					lb=0;
					ub=100;
					cost=10000;
					fp_overlap<< "a " << node_1 << " " << node_2 << " " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";
					fp_overlap<< "a " << node_3 << " " << node_4 << " " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";
				}
				else if(a_statistics[1]>=aStatisticsThreshold && delta>=minDelta) //these edges should traversed exactly once;
				{
					lb=1;
					ub=1;
					cost=1;
					fp_overlap<< "a " << node_1 << " " << node_2 << " " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";
					fp_overlap<< "a " << node_3 << " " << node_4 << " " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";
				}
				else // these are composite edges and should be used 1 to any-number of times times
				{
					flow_lb=0; // default flow lower bound. To be in the safe side we can set it to 0 for real reads
					if(delta>minDelta)
					{
						for(j=1; j<=19; j++)
						{
							if (a_statistics[j]<-aStatisticsThreshold && a_statistics[j+1]>aStatisticsThreshold)
							{
								flow_lb=j+1;
								break;
							}
						}
					}
					else
					{
						if(v->listOfReads[0].readId>35) // composite edges with many reads should have flow at least 1
							flow_lb=1;
					}
					cost_1=0;
					cost_2=0;
					cost_3=0;
					uint64_t n=loaderObj->numberOfReads, N=genomeSize;
					for(index=1; index<=v->listOfReads[0].readId; index++)
					{
						r1 = loaderObj->getRead(v->listOfReads[index].readId);
						xi = r1->frequency;
						cost_1=cost_1+(-(xi*log((float)(first_UB)))-(n-xi)*log(N-first_UB)-(-(xi*log((float)(1)))-(n-xi)*log((float)(N-1))));
						cost_2=cost_2+(-(xi*log((float)(second_UB)))-(n-xi)*log((float)(N-second_UB))-(-(xi*log((float)(first_UB)))-(n-xi)*log((float)(N-first_UB))));
						cost_3=cost_3+(-(xi*log((float)(third_UB)))-(n-xi)*log(N-third_UB)-(-(xi*log((float)(second_UB)))-(n-xi)*log((float)(N-second_UB))));
					}
					lb_1=min(flow_lb, first_UB);
					lb_2=max(0, min(flow_lb, second_UB)-first_UB);
					lb_3=max(0, min(flow_lb,third_UB)-second_UB);
					ub_1=first_UB;
					ub_2=second_UB-first_UB;
					ub_3=third_UB-second_UB;

					cost_1=cost_1/(ub_1-1);
					cost_2=cost_2/ub_2;
					cost_3=cost_3/ub_3;

					fp_overlap<< "a " << node_1 << " " << node_2 << " " << lb_1 << " " << ub_1 << " " << cost_1*cost_scale_factor << "\n";
					fp_overlap<< "a " << node_3 << " " << node_4 << " " << lb_1 << " " << ub_1 << " " << cost_1*cost_scale_factor << "\n";

					fp_overlap<< "a " << node_1 << " " << node_2 << " " << lb_2 << " " << ub_2 << " " << cost_2*cost_scale_factor << "\n";
					fp_overlap<< "a " << node_3 << " " << node_4 << " " << lb_2 << " " << ub_2 << " " << cost_2*cost_scale_factor << "\n";

					fp_overlap<< "a " << node_1 << " " << node_2 << " " << lb_3 << " " << ub_3 << " " << cost_3*cost_scale_factor << "\n";
					fp_overlap<< "a " << node_3 << " " << node_4 << " " << lb_3 << " " << ub_3 << " " << cost_3*cost_scale_factor << "\n";
				}
			}
		}
	}
	fp_overlap.close();
	seconds_start_CS2=time(NULL);
	logStream<< "\n\tCalling CS2()\n";
	cost_of_cs2=main_cs2(cs2InputFileName.c_str(), cs2OutputFileName.c_str()); // cs2 is called and the output is written to file
	logStream<< "\tCost of CS2 : " << cost_of_cs2 << "\n";
	logStream<< "\tFunction CS2() in " << time(NULL)-seconds_start_CS2 << " sec.\n\n";
	ifstream fp_flow(cs2OutputFileName.c_str());
	if(fp_flow.is_open()==false)
		printError(OPEN_FILE, cs2OutputFileName);
	for(i=0; i<numberOfEdges; i++)
		fp_flow >> from[i] >> to[i] >> flow[i];
	fp_flow.close();

	for(i=0; i<numberOfEdges; i++) // for flow in the edges.
	{
		uint64_t from_node, to_node;
		if((from[i]==1 && to[i]==2) || (from[i]==2 && to[i]==1)) // between source and sink
			logStream<< "\t         Flow from T to S: " << flow[i] << "\n";
		else if(from[i]>2 && to[i]>2) // other nodes
		{
			from_node=inverse_node_index[from[i]];
			to_node=inverse_node_index[to[i]];
			if(from_node!=to_node && to_node>=1 && from_node>=1)
			{
				for(v=graphObj->graph[from_node]; v!=NULL; v=v->next)
				{
					if(v->ID==to_node)
					{
						v->flow+=flow[i];
						break;
					}
				}
			}
		}
	}
	uint64_t flow_different=0;
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++) //take average of flow in the twin edges
	{
		for(v=graphObj->graph[i]; v!=NULL; v=v->next)
		{
			float averageFlow = (v->flow+v->twinEdge->flow)/2;
			if(v->flow!=v->twinEdge->flow)
			{
				flow_different++;
			}
			v->flow = averageFlow;
			v->twinEdge->flow = averageFlow;
		}
	}
	free(from);
	free(to);
	free(flow);
	free(node_index);
	free(inverse_node_index);
	logStream<< "\t  Flow different in edges: " << flow_different << "\n";
	logStream<< "Function computeMinCostFlow() in " << time(NULL)-seconds_s << " sec.\n";
	if(saveAll==false)
	{
		deleteFile(cs2InputFileName);
		deleteFile(cs2OutputFileName);
	}
}

/* ============================================================================================
   This function generates flow in the graph. New implementation.
   ============================================================================================ */
void CopyCountEstimator::computeMinCostFlow_new(string outDir)
{
	logStream<< "\nIn function computeMinCostFlow().\n";
	logStream.flush();
	time_t seconds_s=time(NULL), seconds_start_CS2;
	uint64_t edgesUsedMoreThanOnce=0, simple_edge=0, edgesUsedExactlyOnce=0, nodes_in_graph=0, edges_in_graph=0;
	uint64_t numberOfEdges=0, numberOfNodes=0, cnt, sinkToSourceCap;
	uint64_t *node_index, *inverse_node_index, *from, *to, *flow, i, j, k, delta, index;
	int xi, cost_scale_factor=10, flow_lb, lb_1=0, lb_2=0, lb_3=0, ub_1=0, ub_2=0, ub_3=0;
	int lb=0, ub=1000, first_UB=3, second_UB=5, third_UB=1000;
	uint64_t n_1in, n_1out, n_2in, n_2out,node_1=0, node_2=0, node_3=0, node_4=0;
	uint64_t u_1in=0, u_1out=0, u_2in=0, u_2out=0, v_1in=0, v_1out=0, v_2in=0, v_2out=0;
	float cost, cost_1, cost_2, cost_3, a_statistics[25];
	double cost_of_cs2;
	Edge *v;
	Read *r1;
	string cs2InputFileName = outDir+"input.cs2";
	string cs2OutputFileName = outDir+"output.cs2";

	for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
	{
		if(graphObj->graph[i]!=NULL)
			nodes_in_graph++;
		for(v=graphObj->graph[i]; v!=NULL; v=v->next)
		{
			if(i>=v->ID)
			{
				edges_in_graph++;
				if(v->listOfReads[0].readId==0)
					simple_edge++;
				else
				{
					delta=v->lengthOfEdge-1;
					k=0;
					for(index=1; index<=v->listOfReads[0].readId; index++)
					{
						r1 = loaderObj->getRead(v->listOfReads[index].readId);
						k += r1->frequency;
					}
					for(j=1;j<=20;j++)
						a_statistics[j]=(float)((float)delta*(float)((float)loaderObj->numberOfReads/(float)genomeSize)-(float)((float)k*log((float)((float)(j+1)/(float)j))));
					if(a_statistics[1]>=aStatisticsThreshold && delta>=minDelta)
						edgesUsedExactlyOnce++;
					else
						edgesUsedMoreThanOnce++;
				}
			}
		}
	}
	logStream<<  "\t          Number of nodes: " << nodes_in_graph << "\n"
			<< "\t          Number of edges: " << edges_in_graph << "\n"
			<<  "\t             Simple edges: " << simple_edge << "\n"
			<<  "\t  Edges used exactly once: " << edgesUsedExactlyOnce << "\n"
			<<  "\tEdges used more than once: " << edgesUsedMoreThanOnce << "\n\n";

	numberOfNodes=nodes_in_graph*4+2;
	numberOfEdges=1+6*nodes_in_graph+2*(edges_in_graph-(edgesUsedMoreThanOnce))+6*(edgesUsedMoreThanOnce);

	if((from = (uint64_t *)malloc(numberOfEdges*sizeof(uint64_t)))==NULL)
		printError(MEM_ALLOC, "from");
	if((to = (uint64_t *)malloc(numberOfEdges*sizeof(uint64_t)))==NULL)
		printError(MEM_ALLOC, "to");
	if((flow = (uint64_t *)malloc(numberOfEdges*sizeof(uint64_t)))==NULL)
		printError(MEM_ALLOC, "flow");

	if((node_index=(uint64_t *)malloc((loaderObj->numberOfUniqueReads+1)*sizeof(uint64_t)))==NULL)
		printError(MEM_ALLOC, "node index");
	for(i=0; i<=loaderObj->numberOfUniqueReads; i++)
		node_index[i]=0;
	if((inverse_node_index=(uint64_t *)malloc((numberOfNodes+4)*sizeof(uint64_t)))==NULL)
		printError(MEM_ALLOC, "inverse_node index");
	for(i=0;i<numberOfNodes+4;i++)
		inverse_node_index[i]=0;

	sinkToSourceCap = numberOfNodes * 1000;

	ofstream fp_overlap(cs2InputFileName.c_str());
	if(fp_overlap.is_open()==false)
		printError(OPEN_FILE, cs2InputFileName);

	fp_overlap<< std::fixed;
	fp_overlap<< "p min " << numberOfNodes << " " << numberOfEdges << "\n";
	fp_overlap<< "n 1 0\n";
	fp_overlap<< "n 2 0\n";
	fp_overlap<< "a 2 1 1 " << sinkToSourceCap << " 0\n";
	cnt=3;
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++) // edges for nodes
	{
		if(graphObj->graph[i]!=NULL)
		{
			n_1in=cnt;
			n_1out=cnt+1;
			n_2in=cnt+2;
			n_2out=cnt+3;

			node_index[i]=cnt;
			inverse_node_index[cnt]=i;
			inverse_node_index[cnt+1]=i;
			inverse_node_index[cnt+2]=i;
			inverse_node_index[cnt+3]=i;

			lb=0;
			ub=1000;
			cost=1000000;
			fp_overlap<< "a 1 " << n_1in << " " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";//source to node with lower bound 0 upper bound 100 and cost 1000000
			fp_overlap<< "a " << n_1out << " 2 " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";//node to sink with lower bound 0 upper bound 100 and cost 1000000
			fp_overlap<< "a 1 " << n_2in << " " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";//source to node with lower bound 0 upper bound 100 and cost 1000000
			fp_overlap<< "a " << n_2out << " 2 " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";//node to sink with lower bound 0 upper bound 100 and cost 1000000
			lb=0;
			ub=1000;
			cost=1;// every read can be used any number of times for real reads. should be used at least once for corrected reads.
			fp_overlap<< "a " << n_1in << " " << n_1out << " " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";//node_in to node_out with lower bound 0 upper bound 100 and cost 1
			fp_overlap<< "a " << n_2in << " " << n_2out << " " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";//node_in to node_out with lower bound 0 upper bound 100 and cost 1
			cnt=cnt+4;
		}
	}
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++) // edges for edges
	{
		for(v=graphObj->graph[i]; v!=NULL; v=v->next)
		{
			if(i>=v->ID)
			{
				j=node_index[v->fromID];
				u_1in=j;
				u_1out=j+1;
				u_2in=j+2;
				u_2out=j+3;
				j=node_index[v->ID];
				v_1in=j;
				v_1out=j+1;
				v_2in=j+2;
				v_2out=j+3;
				delta=v->lengthOfEdge-1;
				k=0;
				if(v->listOfReads[0].readId==0) // these are simple edges
				{
					k=0;
					a_statistics[1]=0;
				}
				else	//composite edges
				{
					for(index=1; index<=v->listOfReads[0].readId; index++)
					{
						r1 = loaderObj->getRead(v->listOfReads[index].readId);
						k += r1->frequency;
					}
					for(j=1; j<=20; j++)
						a_statistics[j]=(float)((float)delta*(float)((float)loaderObj->numberOfReads/(float)genomeSize)-(float)((float)k*log((float)((float)(j+1)/(float)j))));
				}
				if(v->typeOfEdge==0)
				{
					node_1=v_1out;
					node_2=u_1in;
					node_3=u_2out;
					node_4=v_2in;
				}
				else if(v->typeOfEdge==1)
				{
					node_1=u_2out;
					node_2=v_1in;
					node_3=v_2out;
					node_4=u_1in;
				}
				else if(v->typeOfEdge==2)
				{
					node_1=u_1out;
					node_2=v_2in;
					node_3=v_1out;
					node_4=u_2in;
				}
				else if(v->typeOfEdge==3)
				{
					node_1=u_1out;
					node_2=v_1in;
					node_3=v_2out;
					node_4=u_2in;
				}
				if(v->listOfReads[0].readId==0) // simple edge: could be used any number of times
				{
					lb=0;
					ub=1000;
					cost=10000;
					fp_overlap<< "a " << node_1 << " " << node_2 << " " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";
					fp_overlap<< "a " << node_3 << " " << node_4 << " " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";
				}
				else if(a_statistics[1]>=aStatisticsThreshold && delta>=minDelta) //these edges should traversed exactly once;
				{
					lb=1;
					ub=1;
					cost=1;
					fp_overlap<< "a " << node_1 << " " << node_2 << " " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";
					fp_overlap<< "a " << node_3 << " " << node_4 << " " << lb << " " << ub << " " << cost*cost_scale_factor << "\n";
				}
				else // these are composite edges and sould be used 1 to any-number of times times
				{
					//flow_lb=0; // default flow lower bound. To be in the safe side we can set it to 0 for real reads
					flow_lb=1; //Ehsan
					if(delta>minDelta)
					{
						for(j=1; j<=19; j++)
						{
							if (a_statistics[j]<-aStatisticsThreshold && a_statistics[j+1]>aStatisticsThreshold)
							{
								flow_lb=j+1;
								break;
							}
						}
					}
					cost_1=0;
					cost_2=0;
					cost_3=0;
					uint64_t n=loaderObj->numberOfReads, N=genomeSize;
					for(index=1; index<=v->listOfReads[0].readId; index++)
					{
						r1 = loaderObj->getRead(v->listOfReads[index].readId);
						xi = r1->frequency;
						cost_1=cost_1+(-(xi*log((float)(first_UB)))-(n-xi)*log(N-first_UB)-(-(xi*log((float)(1)))-(n-xi)*log((float)(N-1))));
						cost_2=cost_2+(-(xi*log((float)(second_UB)))-(n-xi)*log((float)(N-second_UB))-(-(xi*log((float)(first_UB)))-(n-xi)*log((float)(N-first_UB))));
						cost_3=cost_3+(-(xi*log((float)(third_UB)))-(n-xi)*log(N-third_UB)-(-(xi*log((float)(second_UB)))-(n-xi)*log((float)(N-second_UB))));
					}
					lb_1=min(flow_lb, first_UB);
					lb_2=max(0, min(flow_lb, second_UB)-first_UB);
					lb_3=max(0, min(flow_lb,third_UB)-second_UB);
					ub_1=first_UB;
					ub_2=second_UB-first_UB;
					ub_3=third_UB-second_UB;

					cost_1=cost_1/(ub_1-1);
					cost_2=cost_2/ub_2;
					cost_3=cost_3/ub_3;

					fp_overlap<< "a " << node_1 << " " << node_2 << " " << lb_1 << " " << ub_1 << " " << cost_1*cost_scale_factor << "\n";
					fp_overlap<< "a " << node_3 << " " << node_4 << " " << lb_1 << " " << ub_1 << " " << cost_1*cost_scale_factor << "\n";

					fp_overlap<< "a " << node_1 << " " << node_2 << " " << lb_2 << " " << ub_2 << " " << cost_2*cost_scale_factor << "\n";
					fp_overlap<< "a " << node_3 << " " << node_4 << " " << lb_2 << " " << ub_2 << " " << cost_2*cost_scale_factor << "\n";

					fp_overlap<< "a " << node_1 << " " << node_2 << " " << lb_3 << " " << ub_3 << " " << cost_3*cost_scale_factor << "\n";
					fp_overlap<< "a " << node_3 << " " << node_4 << " " << lb_3 << " " << ub_3 << " " << cost_3*cost_scale_factor << "\n";
				}
			}
		}
	}
	fp_overlap.close();
	seconds_start_CS2=time(NULL);
	logStream<< "\n\tCalling CS2()\n";
	cost_of_cs2=main_cs2(cs2InputFileName.c_str(), cs2OutputFileName.c_str()); // cs2 is called and the output is written to file
	logStream<< "\tCost of CS2 : " << cost_of_cs2 << "\n";
	logStream<< "\tFunction CS2() in " << time(NULL)-seconds_start_CS2 << " sec.\n\n";
	ifstream fp_flow(cs2OutputFileName.c_str());
	if(fp_flow.is_open()==false)
		printError(OPEN_FILE, cs2OutputFileName);
	for(i=0; i<numberOfEdges; i++)
		fp_flow >> from[i] >> to[i] >> flow[i];
	fp_flow.close();

	for(i=0; i<numberOfEdges; i++) // for flow in the edges.
	{
		uint64_t from_node, to_node;
		if((from[i]==1 && to[i]==2) || (from[i]==2 && to[i]==1)) // between sourse and sink
			logStream<< "\t         Flow from T to S: " << flow[i] << "\n";
		else if(from[i]>2 && to[i]>2) // other nodes
		{
			from_node=inverse_node_index[from[i]];
			to_node=inverse_node_index[to[i]];
			if(from_node!=to_node && to_node>=1 && from_node>=1)
			{
				for(v=graphObj->graph[from_node]; v!=NULL; v=v->next)
				{
					if(v->ID==to_node)
					{
						v->flow+=flow[i];
						break;
					}
				}
			}
		}
	}
	uint64_t flow_different=0;
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++) //take average of flow in the twin edges
	{
		for(v=graphObj->graph[i]; v!=NULL; v=v->next)
		{
			float averageFlow = (v->flow+v->twinEdge->flow)/2;
			if(v->flow!=v->twinEdge->flow)
			{
				flow_different++;
			}
			v->flow = averageFlow;
			v->twinEdge->flow = averageFlow;
		}
	}
	free(from);
	free(to);
	free(flow);
	free(node_index);
	free(inverse_node_index);
	logStream<< "\t  Flow different in edges: " << flow_different << "\n";
	logStream<< "Function computeMinCostFlow() in " << time(NULL)-seconds_s << " sec.\n";
	if(saveAll==false)
	{
		deleteFile(cs2InputFileName);
		deleteFile(cs2OutputFileName);
	}
}
