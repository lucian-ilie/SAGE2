/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 ****************************************************************************/

#include "overlapGraph.h"

extern ofstream logStream;
extern uint64_t genomeSize, averageReadLength;

ostream& operator<<(ostream &os, const Edge *p)
{
	os << p->fromID << "\t" << p->ID << "\t" << p->typeOfEdge << "\t" << p->isReducible << "\t" << p->lengthOfEdge << "\t" << p->flow << "\t" << p->listOfReads[0].readId << "\n";
	for(uint64_t i=1; i<=p->listOfReads[0].readId; i++)
	{
		os << p->listOfReads[i].readId << "\t" << p->listOfReads[i].orientation << "\t" << p->listOfReads[i].flag << "\t" << p->listOfReads[i].distPrevious << "\t" << p->listOfReads[i].distNext << "\n";
	}
	return os;
}

OverlapGraph::OverlapGraph(EconomyGraph *economy1, ReadLoader *loader1)
{
	economyObj = economy1;
	loaderObj = loader1;
	N_90 = 0;
	N_50 = 0;
	N_90_N = 0;
	N_50_N = 0;
	longestContig = 0;
	base_pair_covered = 0;
	contigLengthThreshold = (int)ceilf(100*(averageReadLength/100.0));
	uint64_t i;
	if((graph=(Edge**)malloc((loaderObj->numberOfUniqueReads+1)*sizeof(Edge*)))==NULL)
		printError(MEM_ALLOC, "graph");

	#pragma omp parallel for
	for(i=0; i<=loaderObj->numberOfUniqueReads; i++)
		graph[i]=NULL;
}

OverlapGraph::OverlapGraph(ReadLoader *loader1)
{
	economyObj = NULL;
	loaderObj = loader1;
	N_90 = 0;
	N_50 = 0;
	N_90_N = 0;
	N_50_N = 0;
	longestContig = 0;
	base_pair_covered = 0;
	contigLengthThreshold = (int)ceilf(100*(averageReadLength/100.0));
	uint64_t i;
	if((graph=(Edge**)malloc((loaderObj->numberOfUniqueReads+1)*sizeof(Edge*)))==NULL)
		printError(MEM_ALLOC, "graph");

	#pragma omp parallel for
	for(i=0; i<=loaderObj->numberOfUniqueReads; i++)
		graph[i]=NULL;
}

OverlapGraph::~OverlapGraph()
{
	Edge *edge_v, *edge_v_next;
	uint64_t i;
	#pragma omp parallel for private(i, edge_v, edge_v_next)
	for(i=0; i<=loaderObj->numberOfUniqueReads; i++)
	{
		if(graph[i]!=NULL)
		{
			for(edge_v=graph[i]; edge_v!=NULL; edge_v=edge_v_next)
			{
				edge_v_next=edge_v->next;
				deleteEdge(edge_v); 
			}
		}
	}
	free(graph); 
}

/* ============================================================================================
   Converts the economy graph to overlap graph.
   ============================================================================================ */
void OverlapGraph::convertGraph()
{
	logStream<< "\nIn function convertGraph().\n";
	logStream.flush();
	time_t second_s=time(NULL);
	uint64_t i;
	uint64_t j, nodea, nodeb;
	uint32_t length;
	uint8_t type;
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++) // Total execution time O(nd)
	{
		nodea=i;
		if(economyObj->economyGraphList[i]!=NULL) // Delete the edges of a node in O(d) time.
		{
			for(j=1; j<=economyObj->economyGraphList[i][0].readId; j++)
			{
				nodeb=economyObj->economyGraphList[i][j].readId;
				if(j>1 && economyObj->economyGraphList[i][j-1].readId==economyObj->economyGraphList[i][j].readId && economyObj->economyGraphList[i][j-1].type==economyObj->economyGraphList[i][j].type)
					continue;
				if(nodea<nodeb)
				{
					type=economyObj->economyGraphList[i][j].type;
					length=economyObj->economyGraphList[i][j].length;
					insertEdgeInGraph(nodea, nodeb, length, type);
				}
			}
			free(economyObj->economyGraphList[i]);
		}
	}
	logStream<< "Function convertGraph() in " << time(NULL)-second_s << " sec.\n";
	logStream.flush();
}

/* ============================================================================================
   This function inserts edges in the overlap graph.
   ============================================================================================ */
void OverlapGraph::insertEdgeInGraph(uint64_t from, uint64_t to, uint32_t overlap_size, uint8_t type)
{
	ReadOnEdge *listForward=NULL, *listReverse=NULL;
	if((listForward=(ReadOnEdge*)malloc((1)*sizeof(ReadOnEdge)))==NULL)
		printError(MEM_ALLOC, "listForward");
	if((listReverse=(ReadOnEdge*)malloc((1)*sizeof(ReadOnEdge)))==NULL)
		printError(MEM_ALLOC, "listReverse");
	listForward[0].readId = 0;
	listReverse[0].readId = 0;
	insertEdgeSimple(from, to, type, listForward, listReverse, 0.0, overlap_size);
}

/* ============================================================================================
   Insert edge in the overlap graph. It does not check for duplicates and also does not consider
   transitive edges.
   ============================================================================================ */
Edge* OverlapGraph::insertEdgeSimple(uint64_t u, uint64_t v, uint8_t type, ReadOnEdge *listForward, ReadOnEdge *listReverse, float flow, uint32_t delta)
{
	Edge *ef, *er; // forward and reverse edges
	if((ef = (Edge*)malloc(sizeof(Edge)))==NULL)
		printError(MEM_ALLOC, "v");
	if((er = (Edge*)malloc(sizeof(Edge)))==NULL)
		printError(MEM_ALLOC, "u");
	uint32_t u_l = loaderObj->getRead(u)->length;
	uint32_t v_l = loaderObj->getRead(v)->length;
	uint32_t u_d = delta;
	uint32_t v_d = u_l - (v_l - u_d);
	ef->fromID=u;					er->fromID=v;
	ef->ID=v;						er->ID=u;
	ef->typeOfEdge=type;			er->typeOfEdge=reverseEdgeType(type);
	ef->isReducible=1;				er->isReducible=1;
	ef->lengthOfEdge=u_d;			er->lengthOfEdge=v_d;
	ef->flow=flow;					er->flow=flow;
	ef->listOfReads=listForward;	er->listOfReads=listReverse;
	ef->next=NULL;					er->next=NULL;
	ef->previous=NULL;				er->previous=NULL;
	ef->twinEdge=er;				er->twinEdge=ef;
	insertIntoList(ef);				insertIntoList(er);
	return ef;
}

/* ============================================================================================
   Insert edge in the overlap graph. It does not check for duplicates and also does not consider
   transitive edges.
   ============================================================================================ */
Edge* OverlapGraph::insertEdge(uint64_t u, uint64_t v, uint8_t type, ReadOnEdge *listForward, ReadOnEdge *listReverse, float flow, uint32_t delta1, uint32_t delta2)
{
	Edge *ef, *er; // forward and reverse edges
	if((ef = (Edge*)malloc(sizeof(Edge)))==NULL)
		printError(MEM_ALLOC, "v");
	if((er = (Edge*)malloc(sizeof(Edge)))==NULL)
		printError(MEM_ALLOC, "u");
	uint32_t u_d = delta1;
	uint32_t v_d = delta2;
	ef->fromID=u;					er->fromID=v;
	ef->ID=v;						er->ID=u;
	ef->typeOfEdge=type;			er->typeOfEdge=reverseEdgeType(type);
	ef->isReducible=1;				er->isReducible=1;
	ef->lengthOfEdge=u_d;			er->lengthOfEdge=v_d;
	ef->flow=flow;					er->flow=flow;
	ef->listOfReads=listForward;	er->listOfReads=listReverse;
	ef->next=NULL;					er->next=NULL;
	ef->previous=NULL;				er->previous=NULL;
	ef->twinEdge=er;				er->twinEdge=ef;
	insertIntoList(ef);				insertIntoList(er);
	return ef;
}

/* ============================================================================================
   This function inserts v in the list of edges of v->from_ID. Double linked list.
   ============================================================================================ */
void OverlapGraph::insertIntoList(Edge *v) // insert in linear time;
{
	if(graph[v->fromID]==NULL) //first time
	{
		graph[v->fromID]=v;
	}
	else //insert in the top
	{
		v->next=graph[v->fromID];
		graph[v->fromID]->previous=v;
		graph[v->fromID]=v;
	}
}

/* ============================================================================================
   Two edges (e1, e2) and (e2,e3) is merged to (e1,e3) here.
   ============================================================================================ */
Edge* OverlapGraph::mergeEdges(Edge *edge1, Edge *edge2, int type_of_merge)
{
	int type=combinedEdgeType(edge1, edge2);
	Edge *return_value;
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
	return_value=insertEdge(edge1->fromID, edge2->ID, type, getListOfReads(edge1,edge2), getListOfReads(edge2->twinEdge,edge1->twinEdge), flow, edge1->lengthOfEdge+edge2->lengthOfEdge, edge1->twinEdge->lengthOfEdge+edge2->twinEdge->lengthOfEdge);
	if(edge1->flow<=flow)
	{
		deleteEdge(edge1->twinEdge);
		deleteEdge(edge1);
	}
	else
	{
		edge1->twinEdge->flow-=flow; edge1->flow-=flow;
	}

	if(edge2->flow<=flow)
	{
		deleteEdge(edge2->twinEdge);
		deleteEdge(edge2);
	}
	else
	{
		edge2->twinEdge->flow-=flow; edge2->flow-=flow;
	}
	return return_value;
}

/* ============================================================================================
   This function will return type of edge if you merge edge1 and edge2
   ============================================================================================ */
int OverlapGraph::combinedEdgeType(Edge *edge1, Edge* edge2)
{
	if((edge1->typeOfEdge==0 && edge2->typeOfEdge==0) || (edge1->typeOfEdge==1 && edge2->typeOfEdge==2)) return 0;
	else if((edge1->typeOfEdge==0 && edge2->typeOfEdge==1) || (edge1->typeOfEdge==1 && edge2->typeOfEdge==3)) return 1;
	else if((edge1->typeOfEdge==2 && edge2->typeOfEdge==0) || (edge1->typeOfEdge==3 && edge2->typeOfEdge==2)) return 2;
	else if((edge1->typeOfEdge==2 && edge2->typeOfEdge==1) || (edge1->typeOfEdge==3 && edge2->typeOfEdge==3)) return 3;
	else return -1;
}

/* ============================================================================================
   This function deletes the edge between vertex_a and vertex_b.
   ============================================================================================ */
int OverlapGraph::deleteEdge(Edge *edge) // Delete in constant time.
{
	if(edge->previous==NULL && edge->next==NULL) //Only one element
	{
		graph[edge->fromID]=NULL;
	}
	else if(edge->previous==NULL && edge->next!=NULL) //First element
	{
		graph[edge->fromID]=edge->next;
		edge->next->previous=NULL;
	}
	else if(edge->previous!=NULL && edge->next==NULL) //Last element;
	{
		edge->previous->next=NULL;
	}
	else // Middle element
	{
		edge->previous->next=edge->next;
		edge->next->previous=edge->previous;
	}
	free(edge->listOfReads);
	free(edge); // free the edge
	return 1;
}

/* ============================================================================================
   Merges the list of reads in two edges.
   ============================================================================================ */
ReadOnEdge* OverlapGraph::getListOfReads(Edge *edge1, Edge *edge2)
{
	ReadOnEdge *returnList;
	uint64_t node=edge1->ID;
	uint8_t orientation=0, flag=0;
	uint32_t distFromPrev, distToNext;
	uint32_t i;
	if(edge1->typeOfEdge==1 || edge1->typeOfEdge==3)
		orientation=1;

	if((returnList=(ReadOnEdge*)malloc((edge1->listOfReads[0].readId+edge2->listOfReads[0].readId+2)*sizeof(ReadOnEdge)))==NULL)
		printError(MEM_ALLOC, "returnList");
	returnList[0].readId=edge1->listOfReads[0].readId+edge2->listOfReads[0].readId+1;

	if(edge1->listOfReads[0].readId==0)
		 distFromPrev = edge1->lengthOfEdge; // If edge1 is simple edge
	else
		distFromPrev = edge1->listOfReads[edge1->listOfReads[0].readId].distNext; // If edge1 is composite edge

	if(edge2->listOfReads[0].readId==0)
		distToNext = edge2->lengthOfEdge; // If edge2 is simple edge
	else
		distToNext = edge2->listOfReads[1].distPrevious; // If edge2 is composite edge

	for(i=1; i<=edge1->listOfReads[0].readId; i++) // Copy list from edge1
		returnList[i]=edge1->listOfReads[i];

	// Insert current node
	returnList[edge1->listOfReads[0].readId+1].readId = node;
	returnList[edge1->listOfReads[0].readId+1].orientation = orientation;
	returnList[edge1->listOfReads[0].readId+1].flag = flag;
	returnList[edge1->listOfReads[0].readId+1].distPrevious = distFromPrev;
	returnList[edge1->listOfReads[0].readId+1].distNext = distToNext;

	for(i=1; i<=edge2->listOfReads[0].readId; i++) // Copy list from edge2
		returnList[edge1->listOfReads[0].readId+i+1]=edge2->listOfReads[i];
	return returnList; // Return the new list
}

void OverlapGraph::saveOverlapGraphInFile(string path)
{
	logStream<< "\nIn function saveOverlapGraphInFile().\n";
	time_t seconds_s=time(NULL);
	logStream<< "\tSaving in the file : " << path << "\n";
	logStream.flush();
	ofstream fout(path.c_str());
	if(fout.is_open()==false)
		printError(OPEN_FILE, path);
	//save other info first
	fout<< genomeSize << "\n";
	fout<< loaderObj->numberOfReads << "\n";
//	fout<< loaderObj->numberOfUniqueReads << "\n";
	fout<< averageReadLength << "\n";

	Edge* u;
	for(uint64_t i=1; i<=loaderObj->numberOfUniqueReads; i++)
	{
		for(u=graph[i]; u!=NULL && u->next!=NULL; u=u->next);
		for( ; u!=NULL; u=u->previous)
		{
			if(i<=u->ID)
			{
				fout<< u << "\n";
				fout<< u->twinEdge << "\n";
			}
		}
	}
	fout.close();
	logStream<< "Function saveOverlapGraphInFile in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
}

void OverlapGraph::loadOverlapGraphFromFile(string path)
{
	logStream<< "\nIn function loadOverlapGraphFromFile().\n";
	time_t seconds_s=time(NULL);
	logStream<< "\tLoading from the file : " << path << "\n";
	logStream.flush();

	ifstream fin(path.c_str());
	if(fin.is_open()==false)
		printError(OPEN_FILE, path);
	//load other info first
	fin >> genomeSize;
	fin >> loaderObj->numberOfReads;
//	fin >> loaderObj->numberOfUniqueReads;
	fin >> averageReadLength;

	uint64_t from1, from1Twin, to1, to1Twin, listSize1, listSize1Twin, readID1, j;
	uint16_t type1, type1Twin, reducible1, reducible1Twin, orientation1, flag1;
	uint16_t distPrevious1, distNext1;
	uint32_t delta1, delta1Twin;
	float flow1, flow1Twin;
	ReadOnEdge *listOfReads1, *listOfReads1Twin;
	while(fin >> from1 >> to1 >> type1 >> reducible1 >> delta1 >> flow1 >> listSize1)
	{
		//read edge
		if((listOfReads1=(ReadOnEdge*)malloc((listSize1+1)*sizeof(ReadOnEdge)))==NULL)
			printError(MEM_ALLOC, "listOfReads1");
		listOfReads1[0].readId = listSize1;
		for(j=1; j<=listSize1; j++)
		{
			fin>> readID1 >> orientation1 >> flag1 >> distPrevious1 >> distNext1;
			listOfReads1[j].readId = readID1;
			listOfReads1[j].orientation = orientation1;
			listOfReads1[j].flag = flag1;
			listOfReads1[j].distPrevious = distPrevious1;
			listOfReads1[j].distNext = distNext1;
		}
		//read twin edge
		fin >> from1Twin >> to1Twin >> type1Twin >> reducible1Twin >> delta1Twin >> flow1Twin >> listSize1Twin;
		if((listOfReads1Twin=(ReadOnEdge*)malloc((listSize1Twin+1)*sizeof(ReadOnEdge)))==NULL)
			printError(MEM_ALLOC, "listOfReads1Twin");
		listOfReads1Twin[0].readId = listSize1Twin;
		for(j=1; j<=listSize1Twin; j++)
		{
			fin>> readID1 >> orientation1 >> flag1 >> distPrevious1 >> distNext1;
			listOfReads1Twin[j].readId = readID1;
			listOfReads1Twin[j].orientation = orientation1;
			listOfReads1Twin[j].flag = flag1;
			listOfReads1Twin[j].distPrevious = distPrevious1;
			listOfReads1Twin[j].distNext = distNext1;
		}
		//insert edges into lists
		Edge *v, *u;
		if((v = (Edge*)malloc(sizeof(Edge)))==NULL)
			printError(MEM_ALLOC, "edge v");
		if((u = (Edge*)malloc(sizeof(Edge)))==NULL)
			printError(MEM_ALLOC, "edge u");
		v->fromID=from1;				u->fromID=from1Twin;
		v->ID=to1;						u->ID=to1Twin;
		v->typeOfEdge=type1;			u->typeOfEdge=type1Twin;
		v->isReducible=reducible1;		u->isReducible=reducible1Twin;
		v->lengthOfEdge=delta1;			u->lengthOfEdge=delta1Twin;
		v->flow=flow1;					u->flow=flow1Twin;
		v->listOfReads=listOfReads1;	u->listOfReads=listOfReads1Twin;
		v->next=NULL;					u->next=NULL;
		v->previous=NULL;				u->previous=NULL;
		v->twinEdge=u;					u->twinEdge=v;
		insertIntoList(v);				insertIntoList(u);
	}
	fin.close();
	logStream<< "Function loadOverlapGraphFromFile in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
}

void OverlapGraph::checkIfAllReadsPresent()
{
	logStream<< "\nIn function checkIfAllReadsPresent().\n";
	logStream.flush();
	time_t seconds_s=time(NULL);
	uint8_t *flag;
	uint64_t i, readsMissing=0, index, readsPresent=0;
	Edge *v;
	Read *r1;
	if((flag=(uint8_t*)malloc((loaderObj->numberOfUniqueReads+1)*sizeof(uint8_t)))==NULL)
		printError(MEM_ALLOC, "Flag");

	#pragma omp parallel for
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
		flag[i]=0;
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
	{
		if(graph[i]!=NULL)
		{
			flag[i]=1;
			for(v=graph[i]; v!=NULL; v=v->next)
			{
				for(index=1; index<=v->listOfReads[0].readId; index++)
				{
					flag[v->listOfReads[index].readId]=1;
				}
			}
		}
	}
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
	{
		r1 = loaderObj->getRead(i);
		if(flag[i])
		{
			readsPresent+=r1->frequency;
		}
		else
		{
			readsMissing+=r1->frequency;
		}
	}
	logStream<< "\tNumber of missing reads: " << readsMissing << "\n";
	logStream<< "\tNumber of present reads: " << readsPresent << "\n";
	//set the new number of reads
	loaderObj->numberOfReads = readsPresent;
	free(flag);
	logStream<< "Function checkIfAllReadsPresent in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
}

//for comparison
bool compareEdgeLengthBased(const Edge *first, const Edge *second)
{
	if(first->lengthOfEdge > second->lengthOfEdge)
		return true;
	return false;
}

/* ============================================================================================
   This function prints the overlap graph in overlap_graph.gdl file. The graph can be viewed by
   aisee (free software available at http://www.aisee.com/)
   ============================================================================================ */
void OverlapGraph::printGraph(string fileName, bool isScaff)
{
	logStream<< "\nIn function printGraph().\n";
	logStream.flush();
	time_t seconds_s=time(NULL);
	Edge **contig_edges, *v, *u;
	uint64_t i, lengthOfArray=1000, numberOfContigsFound=0;
	uint16_t color, thickness, flg;

	if((contig_edges=(Edge **)malloc(lengthOfArray*sizeof(Edge *)))==NULL)
		printError(MEM_ALLOC, "contig_edges");

	ofstream fpo_overlap((fileName+".gdl").c_str());
	if(fpo_overlap.is_open()==false)
		printError(OPEN_FILE, fileName+".gdl");

	fpo_overlap<< "graph: {\nlayoutalgorithm :forcedir\nfdmax:704\ntempmax:254\ntempmin:0\ntemptreshold:3\ntempscheme:3\ntempfactor:1.08\nrandomfactor:100\ngravity:0.0\nrepulsion:161\nattraction:43\nignore_singles:yes\nnode.fontname:\"helvB10\"\nnode.shape:box\nnode.width:80\nnode.height:20\nnode.borderwidth:1\nnode.bordercolor:31\n";
	for(i=1; i<=loaderObj->numberOfUniqueReads;i++) // Nodes
		if(graph[i]!=NULL)
			fpo_overlap<< "node: { title:\"" << i << "\" label: \"" << i << "\"}\n";
	fpo_overlap<< "edge.arrowstyle:none\n";

	for(i=1; i<=loaderObj->numberOfUniqueReads; i++) //Edges
	{
		for(v=graph[i]; v!=NULL; v=v->next)
		{
			if(i>=v->ID)
			{
				if(i==v->ID) // do not double count loops (a,a)
				{
					flg=0;
					for(u=graph[i]; u!=v; u=u->next)
					{
						if(u==v->twinEdge)// twin edge already considered
						{
							flg=1; break;
						}
					}
					if(flg==1)
						continue;
				}
				if(v->listOfReads[0].readId!=0)
					thickness=5;//Thick composite edges
				else
					thickness=1;//Thin simple edge
				if(v->typeOfEdge==1)
					color=3; //color code for aisee 3=green
				else if(v->typeOfEdge==2)
					color=1; //1=blue
				else
					color=2; //2=red
				if(N_90!=0 && v->lengthOfEdge>N_90 )
					color=12; //edges longer than N90 12=darkmagenta
				if(!v->isReducible)
					color=31; //edges that are not reducible are color in black
				uint64_t k=0;
				k=v->listOfReads[0].readId;
				if(v->typeOfEdge==0)
					fpo_overlap<< "edge: { source:\"" << i << "\" target:\"" << v->ID << "\" thickness: " << thickness << " backarrowstyle:solid color: " << color << " label: \"Flow: " << v->flow << " Delta: " << v->lengthOfEdge << ", k: " << k << "\" }\n";
				else if(v->typeOfEdge==1)
					fpo_overlap<< "edge: { source:\"" << i << "\" target:\"" << v->ID << "\" thickness: " << thickness << " backarrowstyle:solid arrowstyle:solid color: " << color << " label: \"Flow: " << v->flow << " Delta: " << v->lengthOfEdge << ", k: " << k << "\" }\n";
				else if(v->typeOfEdge==2)
					fpo_overlap<< "edge: { source:\"" << i << "\" target:\"" << v->ID << "\" thickness: " << thickness << " color: " << color << " label: \"Flow: " << v->flow << " Delta: " << v->lengthOfEdge << ", k: " << k << "\" }\n";
				else if(v->typeOfEdge==3)
					fpo_overlap<< "edge: { source:\"" << i << "\" target:\"" << v->ID << "\" thickness: " << thickness << " arrowstyle:solid color: " << color << " label: \"Flow: " << v->flow << " Delta: " << v->lengthOfEdge << ", k: " << k << "\" }\n";
				if(v->lengthOfEdge>=contigLengthThreshold)
				{
					++numberOfContigsFound;
					contig_edges[numberOfContigsFound]=v;
					if(numberOfContigsFound>=lengthOfArray-5)
					{
						lengthOfArray += 1000; // Increase the array size as needed.
						if((contig_edges=(Edge **)realloc(contig_edges,lengthOfArray*sizeof(Edge *)))==NULL)
							printError(MEM_ALLOC,"reallocating contig_edges failed");
					}
				}
			}
		}
	}
	fpo_overlap<< "}}";
	fpo_overlap.close();
	sort(contig_edges+1, contig_edges+numberOfContigsFound+1, compareEdgeLengthBased); // Quicksort the contigs according to their lengths.
	saveContigsToFile(fileName, contig_edges, 1, numberOfContigsFound, isScaff); //Save contigs in a fasta file.
	free(contig_edges);
	logStream<< "Function printGraph() in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
}

int OverlapGraph::getDegree(uint64_t v)
{
	int deg = 0;
	Edge *e;
	for(e=graph[v]; e!=NULL; e=e->next)
		deg++;
	return deg;
}

/* ============================================================================================
   Sort contigs according to their lengths.
   ============================================================================================ */
void OverlapGraph::saveContigsToFile(string fileName, Edge **contig_edges, uint64_t from, uint64_t to, bool isScaff)
{
	logStream<< "\tWriting to file: " << fileName << ".fasta" << "\n\n";
	logStream.flush();
	uint64_t i, j, startIndex, contig_length, contig_length_tmp, N_X=0, N50=0, N90=0, N50_count=0, N90_count=0, gapCount=0, nCount=0, length=0, nextDistance=0, index=0, counter=0, prevDistance=0;
	base_pair_covered=0;
	Edge *v;
	Read *r1, *r2, *rp, *rn;
	string contig_string, str, returnString;

	ofstream fpo_contigs((fileName+".fasta").c_str());
	if(fpo_contigs.is_open()==false)
		printError(OPEN_FILE, fileName+".fasta");

	for(i=from; i<=to; i++)
	{
		if(contig_edges[i]->lengthOfEdge>=contigLengthThreshold)
			base_pair_covered+=contig_edges[i]->lengthOfEdge;
	}
	logStream<< "\tTable: Partial list of long " << (isScaff ? "scaffolds" : "contigs") << " in printGraph().\n";
	logStream<< "\t---------------------------------------------------------------------------\n";
	logStream<< "\t| # |           EDGE          |    LENGTH   |  FLOW | GAP COUNT | N COUNT |\n";
	logStream<< "\t---------------------------------------------------------------------------\n";
	logStream.flush();
	for(i=from; i<=to; i++)
	{
		v=contig_edges[i];

		rn = loaderObj->getRead(v->ID);
		rp = loaderObj->getRead(v->fromID);
		length=v->lengthOfEdge+rp->length;
		contig_string.assign(length, 'N');
		counter=0;
		
		if(v->typeOfEdge==0 || v->typeOfEdge==1)
			str=bytesToChars(rp->readReverseInt, rp->length);
		else
			str=bytesToChars(rp->readInt, rp->length);
		
		for(j=0; j<rp->length; j++)
			contig_string[counter++]=str[j]; // Copy the first read

		if(v->listOfReads[0].readId==0)
			nextDistance=v->lengthOfEdge;
		else
		{
			for(index=1; index<=v->listOfReads[0].readId; index++) // Copy intermediate reads.
			{
				if(v->listOfReads[index].orientation) // Forward orientation
				{
					rn = loaderObj->getRead(v->listOfReads[index].readId);
					str=bytesToChars(rn->readInt, rn->length);
				}
				else // Reverse orientation.
				{
					rn = loaderObj->getRead(v->listOfReads[index].readId);
					str=bytesToChars(rn->readReverseInt, rn->length);
				}
				prevDistance=v->listOfReads[index].distPrevious;
				if(prevDistance>rn->length) // there is a gap
				{
					gapCount += 1;
					nCount += prevDistance-rn->length;
					for(j=0; j<prevDistance-rn->length; j++) // fill in the gaps with N's
						contig_string[counter++]='N';
					for(j=0; j<rn->length; j++) // add the read
						contig_string[counter++]=str[j];
				}
				else//Not a gap, append the string from read
				{
					for(j=(rn->length-prevDistance); j<rn->length; j++)
						contig_string[counter++]=str[j];
				}
				nextDistance=v->listOfReads[index].distNext;
				rp = rn;
			}
		}
		
		if(v->typeOfEdge==1 || v->typeOfEdge==3)
		{
			rn = loaderObj->getRead(v->ID);
			str=bytesToChars(rn->readInt, rn->length);
		}
		else
		{
			rn = loaderObj->getRead(v->ID);
			str=bytesToChars(rn->readReverseInt, rn->length);
		}
		
		for(j=rn->length-nextDistance; j<rn->length; j++)
			contig_string[counter++]=str[j];// Copy the last read.
	
		if(getDegree(v->fromID)==1 && getDegree(v->ID)==1)
		{
			contig_length = contig_string.size();
			startIndex = 0;
		}
		else if(getDegree(v->fromID)==1)
		{
			r1 = loaderObj->getRead(v->fromID);
			r2 = loaderObj->getRead(v->ID);
			contig_length = r1->length + (v->lengthOfEdge - ((r2->length)/2));
			startIndex = 0;
		}
		else if(getDegree(v->ID)==1)
		{
			r1 = loaderObj->getRead(v->fromID);
			r2 = loaderObj->getRead(v->ID);
			contig_length = (r1->length - ((r1->length)/2)) + v->lengthOfEdge;
			startIndex = (r1->length - ((r1->length)/2));
		}
		else
		{
			r1 = loaderObj->getRead(v->fromID);
			r2 = loaderObj->getRead(v->ID);
			contig_length = (r1->length - ((r1->length)/2)) + (v->lengthOfEdge - ((r2->length)/2));
			startIndex = (r1->length - ((r1->length)/2));
		}

		if(i==from)
			longestContig=contig_length;
		if(contig_length>=contigLengthThreshold) //Put contigs in contig file in sorted order
		{
			fpo_contigs<< ">" << (isScaff ? "scaffold" : "contig") << "_" << i << " flow=" << v->flow << " Edge length=" << contig_length << " String length: " << contig_string.size() << " e(" << v->fromID << ", " << v->ID << ") gapCount=" << gapCount << " nCount=" << nCount << "\n";
			contig_length_tmp = contig_length;
			while(contig_length_tmp>=60)
			{
				fpo_contigs<< contig_string.substr(startIndex, 60) << "\n";
				startIndex += 60;
				contig_length_tmp -= 60;
			}
			if(contig_length_tmp>0)
			{
				fpo_contigs<< contig_string.substr(startIndex, contig_length_tmp) << "\n";
			}
		}
		if(i<=10) // Print only large 100 contigs in output file. This slows down the program if there are many contigs.
			logStream<< "\t|" << setw(2) << i << " |" << setw(11) << contig_edges[i]->fromID << ", " << setw(11) << contig_edges[i]->ID << " |" << setw(9) << contig_length << " BP |" << setw(6) << contig_edges[i]->flow << " |" << setw(10) << gapCount << " |" << setw(8) << nCount << " |\n";
		N_X+=contig_length;
		if(!N50_count && N_X>=genomeSize/2)
		{
			 N50=contig_length;
			 N50_count=i;
		}
		if(!N90_count && N_X>=genomeSize*9/10)
		{
			N90=contig_length;
			N90_count=i;
		}
	}
	logStream<< "\t---------------------------------------------------------------------------\n\n";
	N_50=N50;
	N_90=N90;
	N_50_N=N50_count;
	N_90_N=N90_count;
	float genomeCovered=0;
	if(genomeSize)
		genomeCovered=(float)((float)base_pair_covered/(float)genomeSize)*100;
	else
		genomeCovered=0;
	fpo_contigs.close();
	logStream<< "\t        Genome Size: " << genomeSize << " BP\n";
	logStream<< "\t    Number of Reads: " << loaderObj->numberOfReads << "\n";
	logStream<< "\t     Genome Covered: " << base_pair_covered << " BP (" << genomeCovered << "%)\n";
	if(isScaff)
	{
		logStream<< "\tNumber of Scaffolds: " << to << "\n";
		logStream<< "\t   Longest Scaffold: " << longestContig << "\n";
	}
	else
	{
		logStream<< "\t  Number of Contigs: " << to << "\n";
		logStream<< "\t     Longest Contig: " << longestContig << "\n";
	}
	logStream<< "\t                N50: " << N50 << " BP (" << N50_count << ")\n";
	logStream<< "\t                N90: " << N90 << " BP (" << N90_count << ")\n";
	logStream.flush();
}

/* ============================================================================================
   Merges the list of reads in two edges. when the edges does not have common node.
   ============================================================================================ */
ReadOnEdge* OverlapGraph::getListOfReadsInDisconnectedEdges(Edge *edge1, Edge *edge2, int gap, uint64_t *length)
{
	uint64_t node1, orient1, distPrev1, distNext1, node2, orient2, distPrev2=0, distNext2=0;
	ReadOnEdge *returnList;
	uint64_t i;
	int index=1, flag;
	string string1, string2;
	
	node1=edge1->ID;
	Read *r1 = loaderObj->getRead(node1);
	int readLength1 = r1->length;
	if(edge1->typeOfEdge==1 || edge1->typeOfEdge==3)
	{
		orient1=1;
		string1=bytesToChars(r1->readInt, r1->length);
	}
	else
	{
		orient1=0;
		string1=bytesToChars(r1->readReverseInt, r1->length);
	}
	distPrev1=edge1->listOfReads[edge1->listOfReads[0].readId].distNext; // If edge1 is composite edge

	node2=edge2->fromID;
	r1 = loaderObj->getRead(node2);
	int readLength2 = r1->length;
	if(edge2->typeOfEdge==2 || edge2->typeOfEdge==3)
	{
		orient2=1;
		string2=bytesToChars(r1->readInt, r1->length);
	}
	else
	{
		orient2=0;
		string2=bytesToChars(r1->readReverseInt, r1->length);
	}
	distNext2=edge2->listOfReads[1].distPrevious; // If edge2 is composite edge

	if(node1==node2) // only one read will be inserted
	{
		flag=1;
		distNext1=distNext2;
		*length=edge1->lengthOfEdge+edge2->lengthOfEdge;
	}
	else // both the nodes should be inserted
	{
		flag=2;
		int len=stringOverlapSize(string1, string2, readLength1, readLength2);
		if(len>0) // node1 and node2 overlap
		{
			distNext1=len;
			distPrev2=len;
			*length=edge1->lengthOfEdge+edge2->lengthOfEdge+len;
		}
		else // node1 and node2 do not overlap
		{
			if(gap<0) // Concatenate the string
			{
				distNext1=readLength1;
				distPrev2=readLength2;
				*length=edge1->lengthOfEdge+edge2->lengthOfEdge+distNext1;

			}
			else // Add N's between the reads.
			{
				distNext1=gap+readLength1;
				distPrev2=gap+readLength2;
				*length=edge1->lengthOfEdge+edge2->lengthOfEdge+distNext1;
			}
		}
	}

	if((returnList=(ReadOnEdge*)malloc((edge1->listOfReads[0].readId+edge2->listOfReads[0].readId+flag+1)*sizeof(ReadOnEdge)))==NULL)
		printError(MEM_ALLOC, "returnList");
	returnList[0].readId=edge1->listOfReads[0].readId+edge2->listOfReads[0].readId+flag;
	for(i=1; i<=edge1->listOfReads[0].readId; i++) // Copy list from edge1
		returnList[index++]=edge1->listOfReads[i];

	returnList[index].readId = node1; // Insert current node
	returnList[index].orientation = orient1;
	returnList[index].distPrevious = distPrev1;
	returnList[index++].distNext = distNext1;
	if(flag==2)
	{
		returnList[index].readId = node2; // Insert current node
		returnList[index].orientation = orient2;
		returnList[index].distPrevious = distPrev2;
		returnList[index++].distNext = distNext2;
	}

	for(i=1; i<=edge2->listOfReads[0].readId; i++) // Copy list from edge2
		returnList[index++]=edge2->listOfReads[i];
	return returnList; // Return the new list
}

int OverlapGraph::stringOverlapSize(string string1, string string2, int readLength1, int readLength2)
{
	int i, j, misMatchCount;
	for(i=0; i<(readLength1-10); i++)
	{
		misMatchCount=0;
		for(j=0; j<(readLength2-i); j++)
		{
			if(string1[i+j]!=string2[j])
			{
				misMatchCount++;
				if(misMatchCount>((readLength2-i)/10)+1)
				{
					break;
				}
			}
		}
		if(j==(readLength2-i)) // at least 10 base pairs are same
		{
			return i;
		}
	}
	return -1; // do not overlap
}
