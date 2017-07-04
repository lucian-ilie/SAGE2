/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 ****************************************************************************/

#include "simplification.h"

extern ofstream logStream;

/* ============================================================================================
   This function removes composite path in the graph
   ============================================================================================ */
uint64_t contractCompositePaths(OverlapGraph *graphObj, ReadLoader *loaderObj)
{
	logStream<< "\nIn function contractCompositePaths().\n";
	logStream.flush();
	time_t seconds_s=time (NULL);
	uint64_t i, how_many_removed=0, unable_to_remove=0;
	int flag;
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
	{
		if(graphObj->graph[i]!=NULL && graphObj->graph[i]->next!=NULL && graphObj->graph[i]->next->next==NULL) // Nodes with exactly two neighbors
		{
			flag=0;
			Edge *u;
			for(u=graphObj->graph[graphObj->graph[i]->ID]; u!=NULL; u=u->next) // This loop will avoid multiple edges between 2 nodes
			{
				if(u->ID==graphObj->graph[i]->next->ID)
				{
					flag=1; break;
				}
			}
			if(flag==1)
				continue;
			if(graphObj->graph[i]->flow!=graphObj->graph[i]->next->flow)
			{
				logStream<< "In-flow and out-flow different at node "<< i << ". flow: " << graphObj->graph[i]->flow << " vs " << graphObj->graph[i]->next->flow << "\n";
				continue;
			}
			if(graphObj->mergeEdges(graphObj->graph[i]->twinEdge, graphObj->graph[i]->next, 2))
				how_many_removed++;
			else
				unable_to_remove++;
		}
	}
	logStream<< "\t   Nodes removed: "<< how_many_removed << "\n"
			<<  "\tUnable to remove: "<< unable_to_remove << "\n"
			<<  "Function contractCompositePaths() in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
	return how_many_removed;
}

/* ===========================================================================================
   This function will remove dead ends in the graph. These are the reads with errors.
   ============================================================================================ */
uint64_t removeDeadEnds(OverlapGraph *graphObj, ReadLoader *loaderObj, int threshold)
{
	time_t seconds_s=time (NULL);
	logStream<< "In function removeDeadEnds().\n";
	logStream.flush();
	uint64_t i, deleted_reads_with_error=0, flag=0, loops=0;
	Edge *v, *v_next;
	int in_edge, out_edge, count=0;
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
	{
		if(graphObj->graph[i]!=NULL)
		{
			in_edge=0;
			out_edge=0;
			flag=0;
			for(v=graphObj->graph[i]; v!=NULL; v=v->next)
			{
				if(v->listOfReads[0].readId!=0) //composite edge
				{
					count=0;
					count=v->listOfReads[0].readId;
					if(count>threshold) //composite edges with only few reads in it
					{
						flag=1;
						break;
					}
				}
				if(v->fromID==v->ID) // if there is a loop
				{
						flag=1;
						loops++;
						break;
				}

				if(v->typeOfEdge==0 || v->typeOfEdge==1)
					in_edge++;
				else
					out_edge++;
			}
			if(flag==0 && ((in_edge==0 && out_edge>0)||(in_edge>0 && out_edge==0))) // only one type of edge
			{
				for(v=graphObj->graph[i]; v!=NULL; v=v_next)
				{
					v_next=v->next;
					graphObj->deleteEdge(v->twinEdge);
					graphObj->deleteEdge(v);
				}
				deleted_reads_with_error++;
			}
		}
	}
	logStream<< "\tSimple edges deleted: "<< deleted_reads_with_error << "\n"
			 << "\t         Total loops: "<< loops << "\n"
			<< "Function removeDeadEnds() in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
	return deleted_reads_with_error;
}

/* ===========================================================================================
   This function will remove bubbles in the graph. These are the reads with errors.
   ============================================================================================ */
uint64_t removeBubbles(OverlapGraph *graphObj, ReadLoader *loaderObj, int64_t closeLength)
{
	time_t seconds_s=time (NULL);
	logStream<< "In function removeBubbles().\n";
	logStream.flush();
	uint64_t i, deleted_reads_with_error=0;
	Edge *v;
	int in_edge,out_edge;
	Edge *in_edge_p=NULL, *out_edge_p=NULL, *edge_other=NULL;
	uint64_t a, b;
	int64_t delta_1, delta_2, no_1, no_2;
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
	{
		if(graphObj->graph[i]!=NULL)
		{
			in_edge=0;
			out_edge=0;
			for(v=graphObj->graph[i]; v!=NULL; v=v->next)
			{
				if(v->typeOfEdge==0 || v->typeOfEdge==1)
				{
					in_edge++;
					in_edge_p=v->twinEdge;

				}
				else
				{
					out_edge++;
					out_edge_p=v;
				}
			}
			if(in_edge==1 && out_edge==1)
			{
				delta_1=in_edge_p->lengthOfEdge+out_edge_p->lengthOfEdge;
				a=in_edge_p->fromID;
				b=out_edge_p->ID;
				for(v=graphObj->graph[a]; v!=NULL; v=v->next)
				{
					if(v->ID==b)
					{
						edge_other=v;
						delta_2=v->lengthOfEdge;
						no_1=1;
						no_2=0;
						if(in_edge_p->listOfReads[0].readId!=0)
							no_1+=in_edge_p->listOfReads[0].readId;
						if(out_edge_p->listOfReads[0].readId!=0)
							no_1+=out_edge_p->listOfReads[0].readId;
						if(edge_other->listOfReads[0].readId!=0)
							no_2+=edge_other->listOfReads[0].readId;
						if(llabs(delta_1-delta_2)<closeLength) //similar length
						{
							if(no_1<(no_2/2))
							{
								graphObj->deleteEdge(in_edge_p->twinEdge);
								graphObj->deleteEdge(in_edge_p);
								graphObj->deleteEdge(out_edge_p->twinEdge);
								graphObj->deleteEdge(out_edge_p);
								deleted_reads_with_error++;

							}
							if(no_2<(no_1/2))
							{
								graphObj->deleteEdge(edge_other->twinEdge);
								graphObj->deleteEdge(edge_other);
								deleted_reads_with_error++;
							}
						}
						break;
					}
				}
			}
		}
	}
	logStream<< "\tSimple edges deleted: "<< deleted_reads_with_error << "\n"
			<< "Function removeBubbles() in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
	return deleted_reads_with_error;
}

/* ===========================================================================================
   This function will remove bubbles in the graph. These are the reads with errors.
   ============================================================================================ */
uint64_t removeTransitiveEdges(OverlapGraph *graphObj, ReadLoader *loaderObj)
{
	time_t seconds_s=time (NULL);
	logStream<< "In function removeTransitiveEdges().\n";
	logStream.flush();
	uint64_t i, deleted_transitive_edges=0;
	Edge *v;
	int in_edge,out_edge;
	Edge *in_edge_p=NULL, *out_edge_p=NULL, *edge_other=NULL;
	uint64_t a, b;
	int64_t delta_1, delta_2, no_1, no_2;
	for(i=1; i<=loaderObj->numberOfUniqueReads; i++)
	{
		if(graphObj->graph[i]!=NULL)
		{
			in_edge=0;
			out_edge=0;
			for(v=graphObj->graph[i]; v!=NULL; v=v->next)
			{
				if(v->typeOfEdge==0 || v->typeOfEdge==1)
				{
					in_edge++;
					in_edge_p=v->twinEdge;

				}
				else
				{
					out_edge++;
					out_edge_p=v;
				}
			}
			if(in_edge==1 && out_edge==1)
			{
				delta_1=in_edge_p->lengthOfEdge+out_edge_p->lengthOfEdge;
				a=in_edge_p->fromID;
				b=out_edge_p->ID;
				for(v=graphObj->graph[a]; v!=NULL; v=v->next)
				{
					if(v->ID==b)
					{
						edge_other=v;
						delta_2=v->lengthOfEdge;
						no_1=0;
						no_2=0;
						if(in_edge_p->listOfReads[0].readId!=0)
							no_1+=in_edge_p->listOfReads[0].readId;
						if(out_edge_p->listOfReads[0].readId!=0)
							no_1+=out_edge_p->listOfReads[0].readId;
						if(edge_other->listOfReads[0].readId!=0)
							no_2+=edge_other->listOfReads[0].readId;
						if(llabs(delta_1-delta_2)<50) //similar length
						{
							if(no_2==0 && no_1==0)
							{
								graphObj->deleteEdge(edge_other->twinEdge);
								graphObj->deleteEdge(edge_other);
								deleted_transitive_edges++;
							}
						}
						break;
					}
				}
			}
		}
	}
	logStream<< "\tTransitive edges deleted: "<< deleted_transitive_edges << "\n"
			<< "Function removeTransitiveEdges() in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
	return deleted_transitive_edges;
}
