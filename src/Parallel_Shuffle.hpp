/* 
 * File	  : Parallel_Shuffle.cpp 
 * Purpose: Shuffle edges parallely (MPI program) on arbitrary id Galib format graph
 * Output : Write the output graph in adjacency list (Galib) or edge list format
 * Author : Md Hasanuzzaman Bhuiyan
 * Email  : mhb@vbi.vt.edu
 * Date	  : May 2, 2015
 */


#ifndef SHUFFLE_H_
#define SHUFFLE_H_

#include <map>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include <cstring>
#include <sstream>

#include "utility.hpp"
#include "utility2.hpp"
#include "TSet.hpp"
#include "mpi.h"

#ifndef NULL
	#define  NULL 0
#endif

#define P1_Q2 0
#define P1_P2 1

#define PENDING_ARRAY_SIZE 		64
#define BZ_ARRAY_SIZE 			64
#define EDGE_UPDATE_LIST_SIZE 	 8

#define MSG_LEN 	 512
#define MAX_ATTEMPT 1000

#define EDGE_FACTOR 1.20
#define NODE_FACTOR 1.20

//#define DEBUG_DETAIL

using namespace std;

typedef  long LI;
typedef  long Vertex;
typedef  long Vsize;
typedef  long Ssize;
typedef  char Byte;

typedef struct
{
	Vertex u;
	Vertex v;
} Edge;

typedef struct
{
	Vertex u;
	Vertex v;
	int	  iv;	// -1 if it is iq1 of (p1,q1), -2 if it is iq2 of (p2,q2), or processor index for eul array to search eul[iv] array to update iy2
} BzEdge;

typedef struct
{
	Vertex dx2, dy2, diy2, dx1, dy1, type;
} EdgeUpdate;


// ==================== WeightedGraph class =========================

template <typename DType>
class WeightedGraph 
{
	public:
		int 		 myRank;		// rank of this processor
		int 		 noProcs;	// no. of total processors
		
		Vsize     	 n;			// # of total vertices in the graph
		Vsize     	 m;			// # of total edges    in the graph
		Ssize		 noShuffle;	// # of total shuffles in the  graph
		double		 visitRate;	// fraction of edges participating in edge switches, 0 < visitRate <=1
		
		Vsize     	 ln;	// # of local vertices in this partition
		Vsize     	 lm;	// # of local edges    in this partition
		Ssize		 lt;	// # of shuffles in this partition going to take place in this step
		Ssize		 tt;	// # of total shuffles already done in this partition
		
		Set<LI> 	*nlist;	// neighbor list. nlist[i] is the neighbor list of vertex i (where i is the mapped id of a vertex)
		Set<LI> 	*index;	// contains the index of the adjacent edges in the edge list array
		//Set<Byte>   *llist;	// label    list. llist[i] is the label  list of the adjacent edges of vertex i
		//Set<DType>	*wlist;	// weight   list. wlist[i] is the weight list of the adjacent edges of vertex i
		
		map<LI, Vertex> mapTable;	// original id to mapped id for the local nodes
		LI			*name;		// mapped   id to original id for the local nodes
		
		Set<Edge>	edge;	// contains the edges of this partition for efficient uniform random selection in O(1) time
		
		Set<Edge>	pending;// contains the pending edges of this partition that may be added in near future!
		Set<BzEdge>	bz;		// contains the busy edges of this partition that can not be involved in another edge switch
		Set<EdgeUpdate> *eul; // eul -> edgeUpdateList
		
		double		*q;		// array of probabilities, q[i] represents the probability of picking edges from processor P_i
							// q[i] = |E_i|/|E|
		double		*cq;	// cumulative probability array, used to choose a processor to select second edge of an edge switch
		
		// BEGIN: these variables are used for computation in member functions only
		Vertex p1, q1, p2, q2, mp1, mp2, mq1, mq2, iq1, iq2, ni, nj, nk, toss, Pi, Pj, Pk, largeSize, smallSize, mediumSize, discardSize;	// one edge is (p1,q1) and the other one is (p2,q2)
		LI noSteps, *send, *recv, maxDeg, maxBufSize;
		int tag, stepCompMsgCount;
		MPI_Status  status;
		MPI_Request request;
		// END: these variables are used for computation in member functions only
		
		WeightedGraph()
		{
			nlist   = NULL;
			index   = NULL;
			//llist 	= NULL;
			//wlist 	= NULL;
			q	    = NULL;
			cq		= NULL;
			send	= NULL;
			recv	= NULL;
			eul 	= NULL;
			name	= NULL;
			noProcs = stepCompMsgCount = 0;
			n = m = ln = lm = 0;
			noShuffle = lt = tt = 0;
			discardSize	= 1;
			smallSize	= 2;
			mediumSize	= 3;
			largeSize 	= 4;
		}
		
		~WeightedGraph()
		{
			FreeGraph();
		}
		void ReadGraph(const char *gname);
		void ServeAReadMsg(LI*, LI&, LI&);
		void ServeGraphReadMsgsTillSend(LI*, LI&, LI&, LI&);
		void ServeGraphReadMsgs(LI*, LI&, LI&, LI&);
		void WriteGraph(const char *oname, int write_flag);
		void ServeAWriteMsg(LI*);
		void ServeWriteMsgs(LI*, LI&);
		void ServeWriteMsgsTillSend(LI*, LI&);
		void VerifyReadGraph();
		void Initialize(double fraction, Ssize noofsteps, int rank, int size)
		{
			visitRate 	= fraction;
			noSteps		= noofsteps;
			myRank 		= rank;
			noProcs 	= size;
			
			q  	 = new double[noProcs];	// array containing the probability of picking edge from each processor
			cq	 = new double[noProcs];	// cumulative probability array
			send = new LI[largeSize];
			recv = new LI[largeSize];
			
			eul  = new Set<EdgeUpdate> [noProcs];
			for(int i = 0; i<noProcs; i++)
				eul[i].init(EDGE_UPDATE_LIST_SIZE);
			
			pending.init(PENDING_ARRAY_SIZE);
			bz.init(BZ_ARRAY_SIZE);
		}
		
		// return processor id to which this vertex will belong
		int Hash_Function(Vertex vertexId)
		{
			return vertexId%noProcs;
		}
		
		// randomly select one edge and assign the corresponding variables
		// return TRUE if edge selected successfully, otherwise FALSE
		void Select_Edge(Vertex &u, Vertex &v, Vertex &mu, Vertex &mv, Vertex &iv)
		{
			if(lm > 0)
			{
				Vertex f = ((Vsize)rand())%(lm<<1);	// f = x % 2lm
				Vertex e = f%lm;				// randomly select one edge, e -> index of the selected edge in "Edge" array
				toss = f/lm;
				//Vertex e = ((Vsize)rand())%lm;
				//toss = rand()%2;
				mu = edge[e].u;					// mu -> mapped id of u
				u  = name[mu];					//  u -> real id
				iv = edge[e].v;					// iv -> index of v in u's adjacency list
				v  = nlist[mu][iv];				//  v -> real id
				if(Hash_Function(v) == myRank)
					mv = mapTable.at(v);				// mv -> mapped id of v
												// the selected edge is (u,v), where u and v are real node ids
				#ifdef DEBUG_DETAIL
					printf("Proc %d selected edge (%ld, %ld) = (%ld, %ld) at index %ld with toss = %ld\n", myRank, u, v, edge[e].u, edge[e].v, e, toss);
				#endif
			}
			else
			{
				cout << "\n\nNumber of local edge lm became 0 at processor " << myRank << "\n\nTERMINATING the program\n\n";
				Abnormally_Terminate();
			}
		}
		
		// checks  whether loop will be created after edge switching btwn edges (u1,v1) and (u2,v2)
		// returns TRUE if loop is created, FALSE otherwise
		bool Create_Loop(Vertex u1, Vertex v1, Vertex u2, Vertex v2)
		{
			return (u1==v2 || u2==v1);
		}
		
		// checks  whether edge switching btwn edges (u1,v1) and (u2,v2) are useless [i.e., the new edges will remain same as before]
		// returns TRUE if edge switches are useless, FALSE otherwise
		bool Useless_Switch(Vertex u1, Vertex v1, Vertex u2, Vertex v2)
		{
			return (u1==u2 || v1==v2);
		}
		
		// checks  whether edge switching btwn local edges (p1,q1) and (p2,q2) are useless or create loop
		// this function is written for ONLY LOCAL switches, where both edges belong to the same processor
		// returns TRUE if edge switches are useless or create loop, FALSE otherwise
		bool Create_Loop_Useless_Switch_Local()
		{
			return (p1==p2 || q1==q2 || p1==q2 || p2==q1);
		}
		
		// checks whether edge switching between edges (u1,v1) and (u2,v2) are useless or create loop
		// this function is for only GLOBAL switching
		// checking u1==u2 is not required since two edges are selected from two different processors
		// returns TRUE if edge switches are useless or create loop, FALSE otherwise
		bool Create_Loop_Useless_Switch_Global(Vertex u1, Vertex v1, Vertex u2, Vertex v2)
		{
			return (v1==v2 || u2==v1 || u1==v2);
		}
		
		// checks whether parallel edge will be created: possibility of (u,v)
		// u, v -> real id of nodes; mu, mv -> mapped id of u and v respectively
		// returns TRUE if edge switching creates parallel edge, FALSE otherwise
		bool Create_Parallel_Edge(Vertex mu, Vertex v)
		{
			return nlist[mu].member(v);
		}
		
		
		void Insert_In_EUL(int index, EdgeUpdate a)
		{
			eul[index].Dinsert(a);
			
			#ifdef DEBUG_DETAIL
				printf("\tProc %d added an eul edge (dx1,dy1,dx2,dy2,diy2,type) = (%ld,%ld,%ld,%ld,%ld,%ld) in the eul[%d] array ... eul[%ld].Size = %ld\n", myRank, a.dx1, a.dy1, a.dx2, a.dy2, a.diy2, a.type, index, index, eul[index].Size());
			#endif
		}
		
		bool Get_Edge_From_EUL(int index, Vertex cx1, Vertex cy1, EdgeUpdate &temp)
		{
			bool val = false;
			LI eulSize = eul[index].Size(), counter;
			
			for (counter = 0; counter < eulSize; counter++)
			{
				if(eul[index][counter].dx1 == cx1 && eul[index][counter].dy1 == cy1)
				{
					temp = eul[index][counter];
					if(eul[index].decrementSize() == -1)
						cout << "Error inside Get_Edge_From_EUL\n";
					eul[index][counter] = eul[index][eul[index].Size()];
					val = true;
					break;
				}
			}
			
			if(val == false)
			{
				eulSize = eul[index].Size();
				cout << "\n\nError ! Error !! Error !!!\n\n(" << cx1 << "," << cy1 << ") edge not found in the EUL array at index " << index << " at proc " << myRank << " while trying to remove it\n";
				cout << "Size of eul[" << index << "] array : " << eulSize << "\ncontent of eul[" << index << "] array : ";
					
				for (counter = 0; counter < eulSize; counter++)
					cout << "(" << eul[index][counter].dx2 << "," << eul[index][counter].dy2 << "," << eul[index][counter].diy2 << "," << eul[index][counter].dx1 << "," << eul[index][counter].dy1 << eul[index][counter].type << "), ";
					
				cout << "\n\n";
			}
			
			#ifdef DEBUG_DETAIL
				if(val == true)
					printf("\tProc %d removed eul edge (%ld, %ld) from eul[%d] array ... eul[%ld].Size = %ld\n", myRank, cx1, cy1, index, index, eul[index].Size());
			#endif
			
			return val;
		}
		
		// insert the edge (x,y) in the pending array, as this is a potential new edge
		// and return the noProcs of the pending array
		LI Insert_In_Pending(Vertex x, Vertex y)
		{
			pending.Dinsert((Edge){x,y});
			
			#ifdef DEBUG_DETAIL
				printf("\tProc %d added a potential edge (%ld, %ld) in the pending array ... pendSize = %ld\n", myRank, x, y, pending.Size());
			#endif
			
			return pending.Size();
		}
		
		
		// checks whether adding the potential new edge (x,y) is safe or not 
		// ensuring whether the same edge is waiting or not in the pending array to get added later by some other edge switch operation
		// return TRUE if it is in the pending array, meaning we can NOT add this edge in the graph as it is NOT SAFE
		// return FALSE, if there is no pending edge (x,y), so it is SAFE to add the edge (x,y)
		bool Potential_Parallel_Edge(Vertex x, Vertex y)
		{
			LI pendingSize = pending.Size();
			
			for (LI counter = 0; counter < pendingSize; counter++)
				if(pending[counter].u == x && pending[counter].v == y)
				{
					#ifdef DEBUG_DETAIL
						printf("\tContinue because of possibility of POTENTIAL parallel edge (%ld, %ld) at proc %d\n", x, y, myRank);
					#endif
					return true;
				}
			
			return false;
		}
		
		
		// remove the edge (x,y) from the pending array: potential edge (x,y) was inserted into this array while there was a possibility of creating the new edge (x,y)
		// if NOT found in the pending array, then something is definitely WRONG: 
		//		- may be the edge was not inserted but the remove function has been called !
		//		- may be the edge getting deleted for different edge switch !
		// return whether the egde was removed successfully or not: TRUE - successful removal, FALSE - unsuccessful
		bool Remove_From_Pending(Vertex x, Vertex y)
		{
			bool val = false;
			LI pendSize = pending.Size(), counter;
			
			for (counter = 0; counter < pendSize; counter++)
			{
				if(pending[counter].u == x && pending[counter].v == y)
				{
					if(pending.decrementSize() == -1)
						cout << "Error inside Remove_From_Pending\n";
					pending[counter] = pending[pending.Size()];
					val = true;
					break;
				}
			}
			
			if(val == false)
			{
				pendSize = pending.Size();
				cout << "\n\nError ! Error !! Error !!!\n\n(" << x << "," << y << ") edge not found in the pending array at proc " << myRank << " while trying to remove it\n";
				cout << "Size of pending array : " << pendSize << "\ncontent of pending array : ";
					
				for (counter = 0; counter < pendSize; counter++)
					cout << "(" << pending[counter].u << "," << pending[counter].v << "), ";
					
				cout << "\n\n";
			}
			
			#ifdef DEBUG_DETAIL
				if(val == true)
					printf("\tProc %d removed the pending edge (%ld, %ld) from pending array ... pendSize = %ld\n", myRank, x, y, pending.Size());
			#endif
			
			return val;
		}
		
		// new edge will be (u,v) or (v,u)
		// return TRUE if it is safe to add the new edge, FALSE otherwise
		
		bool Can_Add_Edge(Vertex u, Vertex mu, Vertex v, Vertex mv)
		{
			if(u < v)	// new edge will be (u,v)
			{
				if(Create_Parallel_Edge(mu, v) || Potential_Parallel_Edge(u, v))
				{	
					#ifdef DEBUG_DETAIL
						printf("\tContinue because of possibility of parallel edge (%ld, %ld) at proc %d\n", u, v, myRank);
					#endif
					return false;
				}
			}
			else if(Create_Parallel_Edge(mv, u) || Potential_Parallel_Edge(v, u))	// new edge will be (v,u)
			{	
				#ifdef DEBUG_DETAIL
					printf("\tContinue because of possibility of parallel edge (%ld, %ld) at proc %d\n", v, u, myRank);
				#endif
				return false;
			}

			return true;
		}
		
		// insert the edge (x,y) in the bz array, as this edge got involved in edge switching
		// and return the noProcs of the bz array
		LI Insert_In_Busy(Vertex x, Vertex y, int z)
		{
			bz.Dinsert((BzEdge){x,y,z});
			
			#ifdef DEBUG_DETAIL
				printf("\tProc %d added the edge (%ld, %ld) in the bz array ... bzSize = %ld\n", myRank, x, y, bz.Size());
			#endif
			return bz.Size();
		}
		
		
		// checks whether edge (x,y) is busy or not: return TRUE if busy, otherwise return FALSE
		bool Edge_Busy(Vertex x, Vertex y)
		{
			LI bzSize = bz.Size();
			
			for (LI counter = 0; counter < bzSize; counter++)
			{
				if(bz[counter].u == x && bz[counter].v == y)
				{
					#ifdef DEBUG_DETAIL
						printf("\tCONTINUE because of busy edge (%ld, %ld) at proc %d\n", x, y, myRank);
					#endif
					return true;
				}
			}
			
			return false;
		}
		
		
		void Update_iq_If_Busy(Vertex x, Vertex y, LI iy)
		{
			LI bzSize = bz.Size();
			
			for (LI counter = 0; counter < bzSize; counter++)
			{
				if(bz[counter].u == x && bz[counter].v == y)
				{
					if(bz[counter].iv == -1)		// busy as (p1,q1) either in global or local shuffle
						iq1 = iy;
					else if(bz[counter].iv == -2)	// busy as (p2,q2) in local-shuffle
						iq2 = iy;
					else if(bz[counter].iv >= 0)	// busy as (x2,y2) while serving tag=101 msg in a global shuffle
					{
						int index = bz[counter].iv;
						LI eulSize = eul[index].Size(), counter;
			
						for (counter = 0; counter < eulSize; counter++)
						{
							if(eul[index][counter].dx2 == x && eul[index][counter].dy2 == y)
							{
								eul[index][counter].diy2 = iy;
								break;
							}
						}
					}
					else
					{
						printf("\n\nError at proc %d inside Update_iq_If_Busy function ...\n\n", myRank);
						Abnormally_Terminate();
					}
					break;
				}
			}
		}
		
		// remove the edge (x,y) from the busy array: edge (x,y) was inserted into this array while it was busy in edge switching
		// if NOT found in the bz array, then something is definitely WRONG: 
		//		- may be the edge was not inserted but the remove function has been called !
		//		- may be the edge getting deleted for different edge switch !
		// return whether the egde was removed successfully or not: TRUE - successful removal, FALSE - unsuccessful
		bool Remove_From_Busy(Vertex x, Vertex y)
		{
			bool val = false;
			LI bzSize = bz.Size(), counter;
			
			for (counter = 0; counter < bzSize; counter++)
			{
				if(bz[counter].u == x && bz[counter].v == y)
				{
					if(bz.decrementSize() == -1)
						cout << "Error inside Remove_From_Busy\n";
					bz[counter] = bz[bz.Size()];
					val = true;
					break;
				}
			}
			
			if(val == false)
			{
				bzSize = bz.Size();
				cout << "\n\nError ! Error !! Error !!!\n\n(" << x << "," << y << ") edge not found in the bz array at proc " << myRank << " while trying to remove it\n";
				cout << "size of bz array : " << bzSize << "\ncontent of bz array : ";
					
				for (counter = 0; counter < bzSize; counter++)
					cout << "(" << bz[counter].u << "," << bz[counter].v << "), ";
					
				cout << "\n\n";
			}
			
			#ifdef DEBUG_DETAIL
				if(val == true)
					printf("\tProc %d removed the busy edge (%ld, %ld) from the bz array ... bzSize = %ld\n", myRank, x, y, bz.Size());
			#endif
			return val;
		}
		
		
		// Replace edge (u1,v1) by (u1,v2)
		void Replace_u1_v1_by_u1_v2(Vertex mu1, Vertex iv1, Vertex v2)
		{
			#ifdef DEBUG_DETAIL
				printf("Proc %d: edge (%ld, %ld) at index (%ld, %ld) getting replaced by (%ld, %ld)\n", myRank, name[mu1], nlist[mu1][iv1], mu1, iv1, name[mu1], v2);
			#endif
			nlist[mu1][iv1] = v2;	// replace (u1,v1) by (u1,v2)
		}
		
		// Replace edge (u1,v1) by (v2,u1) : remove edge (u1,v1) and add edge (v2,u1) in the same processor
		void Replace_u1_v1_by_v2_u1(Vertex u1, Vertex iv1, Vertex mu1, Vertex mv2)
		{
			#ifdef DEBUG_DETAIL
				printf("Proc %d: edge (%ld, %ld) getting replaced by (%ld, %ld)\n", myRank, name[mu1], nlist[mu1][iv1], name[mv2], u1);
			#endif
			nlist[mv2].Dinsert(u1);		// add    edge (v2,u1)
			nlist[mu1].decrementSize();	// remove edge (u1,v1)
			index[mu1].decrementSize();
			nlist[mu1][iv1] = nlist[mu1][Deg(mu1)];
			edge[index[mu1][Deg(mu1)]] = (Edge){mv2,Deg(mv2)-1};
			index[mv2].Dinsert(index[mu1][Deg(mu1)]);
			
			Update_iq_If_Busy(name[mu1], nlist[mu1][Deg(mu1)], iv1);
			
			if(lm != edge.Size())
			{
				printf("\n\nAfter replacing (u1,v1) by (v2,u1) ... \n\nProc %d: lm != edge array size, lm = %ld, edge array size = %ld\n\n", myRank, myRank, lm, edge.Size());
				Abnormally_Terminate();
			}
		}
		
		// Remove edge (u,v): mu -> mapped id of u, iv -> index of v in u's adjacency list
		void Remove_Edge_u_v(Vertex mu, Vertex iv)
		{
			#ifdef DEBUG_DETAIL
				printf("Proc %d: edge (%ld, %ld) is getting deleted\n", myRank, name[mu], nlist[mu][iv]);
			#endif
			
			nlist[mu].decrementSize();
			index[mu].decrementSize();
			nlist[mu][iv] = nlist[mu][Deg(mu)];
			
			Update_iq_If_Busy(name[mu], nlist[mu][Deg(mu)], iv);
			
			lm--;
			if(lm <= 0)
			{
				cout << "\n\nlocal edge lm became 0 at processor " << myRank << "\n\n";
				Abnormally_Terminate();
			}
			edge.decrementSize();
			edge[index[mu][Deg(mu)]] = edge[lm];
			index[edge[lm].u][edge[lm].v] = index[mu][Deg(mu)];
		}
		
		// Add a new edge (u,v) : mu -> mapped id of u, v -> real id
		void Add_Edge_u_v(Vertex mu, Vertex v)
		{
			nlist[mu].Dinsert(v);
			index[mu].Dinsert(lm);
			edge.Dinsert((Edge){mu,Deg(mu)-1});
			lm++;
			
			#ifdef DEBUG_DETAIL
				printf("Proc %d: Adding new edge (%ld, %ld)\n", myRank, name[mu], v);
			#endif
			
			if(lm != edge.Size())
			{
				printf("\n\nProc %d: After edge add\n\nProc %d: lm != edge array size, lm = %ld, edge array size = %ld\n\n", myRank, myRank, lm, edge.Size());
				Abnormally_Terminate();
			}
			
		}
		
		void Send_Abnormal_Termination_Signal()
		{
			LI temp[1];
			for(int a = 0; a < noProcs; a++)
				if(a != myRank)
					MPI_Send(temp, 1, MPI_LONG, a, 271, MPI_COMM_WORLD);
		}
		
		void Abnormally_Terminate()
		{
			Send_Abnormal_Termination_Signal();
			MPI_Finalize();
		}
		
		void Select_2nd_Edge_For_Msg_With_Tag_101 (Vertex x1, Vertex y1, Vertex &x2, Vertex &y2, Vertex &mx2, Vertex &my2, Vertex &iy2, Vertex &k);
		
		// Processor j (this processor) receives acknowledgement signal 105 from processor i to proceed the pending works of the edge switch
		// Processor i selected edge (x1,y1) and Processor j selected edge (x2,y2)
		// Processor j will now update the edge (x2,y2) with (x2,y1) or (y1,x2)
		void Serve_Msg_With_Tag_105(Vertex x1, Vertex y1, Vertex x2, Vertex y2, Vertex mx2, Vertex iy2, Vertex type)
		{
			if(type == 2)			// replace (x2,y2) by (x2,y1)
			{	
				Replace_u1_v1_by_u1_v2(mx2, iy2, y1);
				Remove_From_Pending(x2,y1);
			}
			else if(type == 3)		// remove edge (x2,y2)
				Remove_Edge_u_v(mx2, iy2);
			else if(type == 4)		// replace (x2,y2) by (y1,x2) : remove edge (x2,y2) and add edge (y1,x2)
			{
				Replace_u1_v1_by_v2_u1(x2, iy2, mx2, mapTable.at(y1));
				Remove_From_Pending(y1,x2);
			}
			else if(type == 22)		// replace (x2,y2) by (x2,x1)
			{	
				Replace_u1_v1_by_u1_v2(mx2, iy2, x1);
				Remove_From_Pending(x2,x1);
			}
			else if(type == 32)		// replace (x2,y2) by (y1,y2) or (y2,y1)
			{	
				if(y1 < y2)
				{
					Replace_u1_v1_by_v2_u1(y2, iy2, mx2, mapTable.at(y1));	// remove (x2,y2), add (y1,y2)
					Remove_From_Pending(y1,y2);
				}
				else
				{
					Replace_u1_v1_by_v2_u1(y1, iy2, mx2, mapTable.at(y2));	// remove (x2,y2), add (y2,y1)
					Remove_From_Pending(y2,y1);
				}
			}
			else
			{
				cout << "\n\nERROR ! ERROR !! ERROR !!!\n\nWhile serving Msg with tag = 5, edge extracted from eul array with type = " << type << " at processor " << myRank << "\n\n";
			}
			
			Remove_From_Busy(x2,y2);
		}
		
		void Shuffle();
		void Local_Shuffle();
		void Global_Shuffle();
		void Send_Step_Complete_Msg();
		bool Serve_Recv_Messages();
		Vsize Deg(Vertex v) {return nlist[v].Size();}
		Vsize Size() {return n;}
		Vertex &operator()(Vertex v, Vsize i) {return nlist[v][i];}
		Set<LI> &operator[](Vertex v) {return nlist[v];}
		void FreeGraph();
};
// ======================= END: WeightedGraph class ===========================



/* =========================== ReadGraph function =============================
 *
 * Reads the Galib format input graph and loads the partitioned data
 *
============================================================================ */



template <typename DType>	// DType is the data type for weight
void WeightedGraph<DType>::ServeAReadMsg(LI *recvBuffer, LI &mapIndex, LI &edgeCounter)
{
	LI nextNodeIndex = 1;
	for(LI j=0; j<recvBuffer[0]; j++)
	{
		mapTable[recvBuffer[nextNodeIndex]] = mapIndex;
		name[mapIndex] = recvBuffer[nextNodeIndex++];
		LI tmpDeg = recvBuffer[nextNodeIndex++];
		LI tmpArrSize = Max(1,(LI)ceil(tmpDeg * NODE_FACTOR));
		
		/***************    Newly added with initial nlist size equal to degree in Galib format    ******************/
		LI galibDeg = recvBuffer[nextNodeIndex++];
		tmpArrSize  = galibDeg;
		/*																											*/
		
		nlist[mapIndex].init(tmpArrSize);
		index[mapIndex].init(tmpArrSize);
			
		for(LI x=0; x < tmpDeg; x++)
		{
			nlist[mapIndex].insert(recvBuffer[nextNodeIndex++]);
			index[mapIndex].insert(edgeCounter);
			edge.insert((Edge){mapIndex, Deg(mapIndex)-1});
			edgeCounter++;
		}
		mapIndex++;
	}
}


template <typename DType>	// DType is the data type for weight
void WeightedGraph<DType>::ServeGraphReadMsgsTillSend(LI *recvBuffer, LI &mapIndex, LI &ackCounter, LI &edgeCounter)
{
	int sendFlag = 0;
	while(!sendFlag)	// until true meaning safe to use the send buffer again
	{
		MPI_Test(&request, &sendFlag, &status);
		int recvFlag = 1;
		while(recvFlag)
		{
			MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &recvFlag, &status);
			if(recvFlag)
				ServeGraphReadMsgs(recvBuffer, mapIndex, ackCounter, edgeCounter);
		}
	}
}


template <typename DType>	// DType is the data type for weight
void WeightedGraph<DType>::ServeGraphReadMsgs(LI *recvBuffer, LI &mapIndex, LI &ackCounter, LI &edgeCounter)
{
	MPI_Recv(recvBuffer, maxBufSize, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	if(status.MPI_TAG == 85)
		ackCounter++;
	else
		ServeAReadMsg(recvBuffer, mapIndex, edgeCounter);
}

template <typename DType>	// DType is the data type for weight
void WeightedGraph<DType>::ReadGraph(const char *gname) // 1st param: original graph name
{	
	int 	temp, flag, sendFlag;
	LI 		i, j, deg, mapIndex, u, v, fileSize, byteQuota, byteStart, byteEnd, nextBytePos, wordCount, proc, mxDegree = 0, ackCounter = 0;
	DType  	wt, number;
	string  line;
	istringstream istr;
	
    ifstream ifp(gname);
	ifp >> n;
	ifp.seekg(0, ifp.end);
    fileSize  = ifp.tellg();
	
	byteQuota = (LI)ceil(double(fileSize)/noProcs);
	byteStart = (myRank == 0) ? 0 : myRank * byteQuota;
	byteEnd   = byteStart + byteQuota;	// for last processor, it may exceed fileSize, hence need to use ifp.good()
	
	ifp.clear();
	ifp.seekg(byteStart, ifp.beg);
	
	if(myRank > 0 && ifp.good())	// may fall in broken line
		getline(ifp, line);
	
	while(ifp.good())
	{
		byteStart = ifp.tellg();
		getline(ifp, line);
		istr.clear();
		istr.str(line);
		
		wordCount = 0;
		while(istr >> number)
			wordCount++;
		if(wordCount == 2)	// Found the beginning of an adjacency list of a vertex in Galib format: <vertexID> <Degree>
			break;
	}
	
	// sending msg containing byte address up to which the previous processor should read
	if(myRank > 0)	// except first processor
		MPI_Send(&byteStart, 1, MPI_LONG, myRank-1, 99, MPI_COMM_WORLD);	// sending the start address of graph read of this processor to previous processor
	
	if(myRank != noProcs-1)	// except last processor
		MPI_Recv(&byteEnd, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	else
		byteEnd = Min(byteEnd, fileSize);
	
	/*
	ifp.clear();
	ifp.seekg(byteStart, ifp.beg);
	ifp >> u >> deg;
	printf("Proc %d will start reading from %ld and ending before %ld, u = %ld, deg = %ld\n", myRank, byteStart, byteEnd, u, deg);
	*/
	
	ifp.clear();
	ifp.seekg(byteStart, ifp.beg);
	nextBytePos = byteStart;
	
	LI *nodeCount = new LI[noProcs];		// no. of node count at each processor
	for(i=0; i<noProcs; i++)
		nodeCount[i] = 0;
	
	LI *edgeCount = new LI[noProcs];		// no. of node count at each processor
	for(i=0; i<noProcs; i++)
		edgeCount[i] = 0;
	
	while(nextBytePos < byteEnd)
	{
		ifp >> u >> deg;
		
		if(!ifp.good() || ifp.tellg() >= byteEnd)	// may reach EOF
			break;
		mxDegree = Max(mxDegree, deg);
		proc = Hash_Function(u);
		nodeCount[proc]++;
		for(i=0; i < deg; i++)
		{
			ifp >> v >> wt >> temp;
			if(u < v)
				edgeCount[proc]++;
		}
		nextBytePos = 1 + ifp.tellg();	// adding 1 required for the new-line/carriage-return
	}
	
	/*
	for(i=0; i<noProcs; i++)
		printf("(%ld, %ld, %ld, %ld) ", myRank, i, nodeCount[i], edgeCount[i]);
	printf("\n\n");
	*/
	
	int *recvCount = new int[noProcs];
	for(i = 0; i < noProcs; i++)
		recvCount[i] = 1;
	
	MPI_Reduce_scatter(nodeCount, &ln, recvCount, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);	// ln = local no of nodes
	MPI_Reduce_scatter(edgeCount, &lm, recvCount, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);	// lm = local no of edges
	MPI_Allreduce(&lm, &m, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);						// m  = total no of edges
	MPI_Allreduce(&mxDegree, &maxDeg, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);			// maxDeg = maximum degree of a node
	
	FreeMem(nodeCount);
	FreeMem(edgeCount);
	FreeMem(recvCount);
	
	if (visitRate < 0.9999)
		noShuffle = Ssize( - ((double)m / 2.0) * log(1.0 - visitRate));	// natural logarithm / ln / e-based log
	else 
		noShuffle = Ssize( visitRate * ((double)m / 2.0) * log((double)m));	// natural logarithm / ln / e-based log

	//if(myRank == 0)
		//cout << "Proc " << myRank << ": number of nodes = " << n << ", number of edges = " << m << ", noShuffle = " << noShuffle << endl;
	
	LI edgeAllocSize = (LI)ceil(lm * EDGE_FACTOR);
	
	nlist = new Set<LI>	[ln];
	index = new Set<LI> [ln];
	name  = new LI		[ln];
	edge.init(edgeAllocSize);
	//llist = new Set<Byte>  [ln];
	//wlist = new Set<DType> [ln];		

	maxBufSize 		 = Max (maxDeg+10, MSG_LEN);		//	maximum buffer size, extra 5 because of node id, degree, etc.
	LI *recvBuffer	 = new LI[maxBufSize];
	LI *sendBufIndex = new LI[noProcs];
	LI **sendBuffer  = new LI*[noProcs];
	
	for(i=0; i<noProcs; i++)
	{
		sendBuffer[i] 	 = new LI[maxBufSize];	// buffering messages that will be sent to processor <i>
		sendBuffer[i][0] = 0;	// how many number of nodes' adjacency lists belong to the buffer that will be sent processor <i>
		sendBufIndex[i]  = 1;	// next available index/slot in the buffer for processor <i>
	}
	
	ifp.clear();
	ifp.seekg(byteStart, ifp.beg);
	nextBytePos 	= byteStart;
	mapIndex 		= 0;
	LI edgeCounter 	= 0;
	
	while(nextBytePos < byteEnd)
	{
		ifp >> u >> deg;
		if(!ifp.good() || ifp.tellg() >= byteEnd)	// may reach EOF
			break;
		proc = Hash_Function(u);	// <u> belongs to processor <proc>
		
		if(proc == myRank)
		{
			mapTable[u]  	= mapIndex;
			name[mapIndex] 	= u;
			nlist[mapIndex].init(deg);
			index[mapIndex].init(deg);
			
			for(i=0; i < deg; i++)
			{
				ifp >> v >> wt >> temp;
				if(u < v)
				{
					nlist[mapIndex].insert(v);
					index[mapIndex].insert(edgeCounter);
					edge.insert((Edge){mapIndex, Deg(mapIndex)-1});
					edgeCounter++;
				}
			}
			mapIndex++;
		}
		else
		{
			if(maxBufSize - sendBufIndex[proc] < deg + 5)	// may not have enough space
			{
				MPI_Isend(sendBuffer[proc], sendBufIndex[proc], MPI_LONG, proc, 81, MPI_COMM_WORLD, &request);
				ServeGraphReadMsgsTillSend(recvBuffer, mapIndex, ackCounter, edgeCounter);
				
				sendBuffer[proc][0] = 0;
				sendBufIndex[proc]  = 1;
			}
			
			sendBuffer[proc][0]++;					// one more node <u>
			sendBuffer[proc][sendBufIndex[proc]++] = u;
			LI degIndex = sendBufIndex[proc]++;	// degree (assuming u < v) will be fill up later at degIndex
			sendBuffer[proc][degIndex] = 0;
			
			/***************    Newly added with initial nlist size equal to degree in Galib format    ******************/
			sendBuffer[proc][sendBufIndex[proc]++] = deg;	// put original degree in Galib format
			/*																											*/
			
			for(i=0; i < deg; i++)
			{
				ifp >> v >> wt >> temp;
				if(u < v)
				{
					sendBuffer[proc][sendBufIndex[proc]++] = v;
					sendBuffer[proc][degIndex]++;	// degree assuming u < v
				}
			}
			
			flag = 1;
			while(flag)
			{
				MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
				if(flag)
					ServeGraphReadMsgs(recvBuffer, mapIndex, ackCounter, edgeCounter);
			}
		}
		
		nextBytePos = 1 + ifp.tellg();	// adding 1 required for the end of line
	}
	
	for(i=0; i<noProcs; i++)	// send the remaining adjacency lists
	{
		if(sendBuffer[i][0] && i != myRank)
		{
			MPI_Isend(sendBuffer[i], sendBufIndex[i], MPI_LONG, i, 81, MPI_COMM_WORLD, &request);
			ServeGraphReadMsgsTillSend(recvBuffer, mapIndex, ackCounter, edgeCounter);
		}
	}
	
	for(i=0; i<noProcs; i++)	// send acknowledgement
	{
		if(i != myRank)
		{
			MPI_Isend(sendBuffer[i], 1, MPI_LONG, i, 85, MPI_COMM_WORLD, &request);		// acknowledgement
			ServeGraphReadMsgsTillSend(recvBuffer, mapIndex, ackCounter, edgeCounter);
		}
	}
	
	while(ackCounter != noProcs-1)
		ServeGraphReadMsgs(recvBuffer, mapIndex, ackCounter, edgeCounter);
	
	ifp.close();
	
	if(mapIndex != ln)
		printf("Error at proc %d, ln = %ld, mapIndex = %ld\n", myRank, ln, mapIndex);
	
	if(edgeCounter != lm)
		printf("Error at proc %d, lm = %ld, edgeCounter = %ld\n", myRank, lm, edgeCounter);
	
	for(i=0; i<noProcs; i++)
		FreeMem(sendBuffer[i]);
	FreeMem(sendBuffer);
	FreeMem(sendBufIndex);
	FreeMem(recvBuffer);
}

// ===================== END: ReadGraph function ==============================


// ====================== VerifyReadGraph function ============================

template <typename DType>	// DType is the data type for weight
void WeightedGraph<DType>::VerifyReadGraph( )
{
	char *str = new char[100];
	str[0] = 0;
	
	stringstream ss;
	ss << myRank;
	const char *ch = ss.str().c_str();
	strcpy(str, "initread-");
	strcat(str, ch);
	strcat(str, ".txt");
	
	FILE *fp = fopen(str, "at");
	
	fprintf(fp, "\nThe details of processor id %d of total %d processors is given below:\n\n", myRank, noProcs);
	fprintf(fp, "Total nodes: %ld\nTotal edges: %ld\nLocal nodes: %ld\nLocal edges: %ld\nTotal shuffle: %ld\nvisit rate: %.2lf\nno of steps = %ld\n\n", n, m, ln, lm, noShuffle, visitRate, noSteps);
	
	LI i, u, v, sizeu;
	fprintf(fp, "Edge list of proc %d: total %ld edges, maxSize = %ld\n", myRank, lm, edge.getMaxSize());
	for (i=0; i<lm; i++)
		fprintf(fp, "\t(%ld, %ld) = (%ld, %ld)\n", edge[i].u, edge[i].v, name[edge[i].u], nlist[edge[i].u][edge[i].v]);
	fprintf(fp, "\n\n");
	
	fprintf(fp, "Adjacency list of proc %d: total %ld nodes\n", myRank, ln);
 	for(u=0; u < ln; u++)
	{
		sizeu = Deg(u);
		fprintf(fp, "%ld (%ld) %ld %ld\n", name[u], u, sizeu, nlist[u].getMaxSize());
		
		for (i=0; i<sizeu; i++)			// scanning all the adjacent nodes
		{
			v = nlist[u][i];
			fprintf(fp, "\t%ld\n", v);
		}
		for (; i < nlist[u].getMaxSize(); i++)
			fprintf(fp, "\t\t%ld", nlist[u][i]);
		fprintf(fp, "\n");
	}
	
	fclose(fp);
	delete []str;
}
// =================== END: VerifyReadGraph function ==========================




template <typename DType>	// DType is the data type for weight
void WeightedGraph<DType>::ServeAWriteMsg(LI *recvBuffer)
{
	LI nextNodeIndex = 1;
	for(LI j=0; j<recvBuffer[0]; j++)
	{
		LI uVertex 	= recvBuffer[nextNodeIndex++];
		LI tmpDeg 	= recvBuffer[nextNodeIndex++];
		for(LI x=0; x < tmpDeg; x++)
		{
			LI v = mapTable.at(recvBuffer[nextNodeIndex++]);
			nlist[v].insert(uVertex);
		}
	}
}



template <typename DType>	// DType is the data type for weight
void WeightedGraph<DType>::ServeWriteMsgs(LI *recvBuffer, LI &ackCounter)
{
	int recvFlag = 1;
	while(recvFlag)
	{
		MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &recvFlag, &status);
		if(recvFlag)
		{
			MPI_Recv(recvBuffer, maxBufSize, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			if(status.MPI_TAG == 75)
				ackCounter++;
			else
				ServeAWriteMsg(recvBuffer);
		}
	}
}


template <typename DType>	// DType is the data type for weight
void WeightedGraph<DType>::ServeWriteMsgsTillSend(LI *recvBuffer, LI &ackCounter)
{
	int sendFlag = 0;
	while(!sendFlag)	// until true meaning safe to use the send buffer again
	{
		MPI_Test(&request, &sendFlag, &status);
		ServeWriteMsgs(recvBuffer, ackCounter);
	}
}


// ==================== WriteGraph function ===================================

template <typename DType>	// DType is the data type for weight
void WeightedGraph<DType>::WriteGraph(const char *oname, int write_flag) // 1st param: output file/graph name
{
	LI u, sizeu, i, j, e, ackCounter = 0, nextNodeIndex;
	int flag;
	char *str = new char[strlen(oname) + 10];
	
	stringstream ss;
	ss << myRank;
	const char *ch = ss.str().c_str();
	strcpy(str, ch);
	strcat(str, oname);
	
	FILE *fp = fopen(str, "wt");
	
	if(!fp)
	{
		cout << "Can not open file: " << oname << endl;
		Abnormally_Terminate();
	}
	
	if(write_flag == 2)		// 2 - edge list
	{
		for(u=0; u < ln; u++)
		{
			sizeu = Deg(u);
			for (i=0; i<sizeu; i++)			// scanning all the adjacent nodes
				fprintf(fp, "%ld\t%ld\n", name[u], nlist[u][i]);
		}
	}
	else if(write_flag == 1)	// 1 - Galib format
	{
		LI *nodeDegree	 = new LI[ln];
		LI *recvBuffer	 = new LI[maxBufSize];
		LI *degIndex 	 = new LI[noProcs];
		LI *sendBufIndex = new LI[noProcs];
		LI **sendBuffer  = new LI*[noProcs];
		for(i=0; i<noProcs; i++)
		{
			sendBuffer[i] 	 = new LI[maxBufSize];	// buffering messages that will be sent to processor <i>
			sendBuffer[i][0] = 0;	// how many number of nodes' adjacency lists belong to the buffer that will be sent processor <i>
			sendBufIndex[i]  = 1;	// next available index/slot in the buffer for processor <i>
		}
	
		for(i=0; i<ln; i++)
			nodeDegree[i] = Deg(i);
			
		for(u=0; u<ln; u++)
		{
			LI uName = name[u];
			LI uDeg  = nodeDegree[u];
			
			if(uDeg > 0)
			{
				for(i=0; i<noProcs; i++)
				{
					if(i != myRank)
					{
						if(maxBufSize - sendBufIndex[i] < uDeg + 5)	// may not have enough space
						{
							MPI_Isend(sendBuffer[i], sendBufIndex[i], MPI_LONG, i, 71, MPI_COMM_WORLD, &request);
							ServeWriteMsgsTillSend(recvBuffer, ackCounter);
							sendBuffer[i][0] = 0;
							sendBufIndex[i]  = 1;
						}
						
						sendBuffer[i][0]++;
						sendBuffer[i][sendBufIndex[i]++] = uName;
						degIndex[i] = sendBufIndex[i];
						sendBuffer[i][sendBufIndex[i]++] = 0;
					}
				}
			}
			
			for(i=0; i<uDeg; i++)
			{
				LI v 	= nlist[u][i];
				LI proc = Hash_Function(v);
				
				if(proc != myRank)
				{
					sendBuffer[proc][sendBufIndex[proc]++] = v;
					sendBuffer[proc][degIndex[proc]]++;
				}
				else
					nlist[mapTable.at(v)].insert(uName);
			}
			
			ServeWriteMsgs(recvBuffer, ackCounter);
		}
		
		for(i=0; i<noProcs; i++)
		{
			if(sendBuffer[i][0] && i != myRank)
			{
				MPI_Isend(sendBuffer[i], sendBufIndex[i], MPI_LONG, i, 71, MPI_COMM_WORLD, &request);
				ServeWriteMsgsTillSend(recvBuffer, ackCounter);
			}
		}
		
		for(i=0; i<noProcs; i++)
		{
			if(i != myRank)
			{
				MPI_Isend(sendBuffer[i], 1, MPI_LONG, i, 75, MPI_COMM_WORLD, &request);		// acknowledgement
				ServeWriteMsgsTillSend(recvBuffer, ackCounter);
			}
		}
		
		while(ackCounter != noProcs-1)
		{
			MPI_Recv(recvBuffer, maxBufSize, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			if(status.MPI_TAG == 75)
				ackCounter++;
			else
				ServeAWriteMsg(recvBuffer);
		}
		
		if(myRank == 0)
			fprintf(fp, "%ld\n", n);
		
		for(u=0; u < ln; u++)
		{
			sizeu = Deg(u);
			fprintf(fp, "%ld\t%ld\n", name[u], sizeu);
			for (i=0; i<sizeu; i++)			// scanning all the adjacent nodes
				fprintf(fp, "\t%ld\t1\t0\n", nlist[u][i]);
		}
		
		for(i=0; i<noProcs; i++)
			FreeMem(sendBuffer[i]);
		FreeMem(sendBuffer);
		FreeMem(sendBufIndex);
		FreeMem(degIndex);
		FreeMem(recvBuffer);
		FreeMem(nodeDegree);
	}
	else 
		cout << "Incorrect write_flag input ...\nEither provide 0 (no write) or 1 (Galib or adjacency format) or 2 (edge list format) as input for write_flag\n";
	
	fclose(fp);
	FreeMem(str);
}


// =================== END: WriteGraph function ===============================


// ======================= FreeGraph function =================================
template <typename DType>
void WeightedGraph<DType>::FreeGraph()
{
	if (nlist)
		for (Vsize i=0; i < ln; i++) 
			nlist[i].destroy();
	FreeMem(nlist);

	if (index)
		for (Vsize i=0; i < ln; i++) 
			index[i].destroy();
	FreeMem(index);

	if (eul)
		for (int j=0; j < noProcs; j++) 
			eul[j].destroy();
	FreeMem(eul);

	FreeMem(send);
	FreeMem(recv);
	bz.destroy();
	pending.destroy();
	edge.destroy();
	FreeMem(q);
	FreeMem(cq);
}
// ====================== END: FreeGraph function =============================



// processor j (this processor) receives the edge (x1,y1) from processor i
// it picks the edge (x2,y2) from its own memory so that no parallel edge will be created
// the edge switch may be a cross switch with probability 0.5 or a straight switch with prob 0.5
template <typename DType>	// DType is the data type for weight
void WeightedGraph<DType>::Select_2nd_Edge_For_Msg_With_Tag_101(Vertex x1, Vertex y1, Vertex &x2, Vertex &y2, Vertex &mx2, Vertex &my2, Vertex &iy2, Vertex &k)
{
	int source = status.MPI_SOURCE;
	
	for(nj = 0; nj < MAX_ATTEMPT; nj++)
	{
		k   = -1;
		tag =  102;
	
		Select_Edge(x2, y2, mx2, my2, iy2);	// x2,y2 -> real ids, the edge is (x2,y2)
											// mx2 -> mapped id of x2
											// my2 -> mapped id of y2 (only if it belongs to this processor)
											// iy2 -> index position of y2 in x2's adjacency list
		
		if(Edge_Busy(x2,y2) || Create_Loop_Useless_Switch_Global(x1, y1, x2, y2))
		{
			#ifdef DEBUG_DETAIL
				printf("\tProc %d going to re-select another edge bcz of parallel edge possibility\n", myRank);
			#endif
			continue;
		}
		
		if(toss == P1_Q2)	// (x1,y1) and (x2,y2) will get replaced by [(x1,y2) or (y2,x1)] and [(x2,y1) or (y1,x2)]
		{
			if(x2 < y1)	// the new edge is going to be (x2,y1)
			{
				if(Create_Parallel_Edge(mx2,y1) || Potential_Parallel_Edge(x2,y1))
				{
					#ifdef DEBUG_DETAIL
						printf("\tProc %d going to re-select another edge bcz of parallel edge possibility\n", myRank);
					#endif
					continue;
				}
				
				if(y2 < x1)
				{
					k = Hash_Function(y2);
					
					#ifdef DEBUG_DETAIL
						printf("\tProc %d: node %ld is located in proc %ld\n", myRank, y2, k);
					#endif
					
					if(k == myRank)
					{
						if(Create_Parallel_Edge(my2, x1) || Potential_Parallel_Edge(y2,x1))
							continue;
						tag = 106;								// signal Pi to remove (x1,y1) since both updates are done at Pj
						Add_Edge_u_v(my2, x1);					// add edge (y2,x1)
						Replace_u1_v1_by_u1_v2(mx2, iy2, y1);	// (p1,q1) getting replaced by (p1,q2) - successful swap
					}
					else if(k != source)
						tag = 108;
				}
				
				if(tag != 106)					// (x2,y2) will get replaced by (x2,y1)
				{
					Insert_In_Busy(x2,y2,source);		// make (x2,y2) busy
					Insert_In_Pending(x2,y1);	// make the edge (x2,y1) a potential edge
					Insert_In_EUL(source, (EdgeUpdate){x2, y2, iy2, x1, y1, 2});	// (x2,y2) will get replaced by (x2,y1)
				}
			}
			else // the new edges are going to be (y1,x2) and (x1,y2) since x1 < y1 < x2 < y2
			{
				k = Hash_Function(y1);		// y1 is located in processor k
			
				#ifdef DEBUG_DETAIL
					printf("\tProc %d: node %ld is located in proc %ld\n", myRank, y1, k);
				#endif
				
				if(k != myRank)
				{
					if(k == source)	// y1 is located in the sending processor Pi
						tag = 143;
					else			// y1 is located in some processor, Pk, in between the sending proc. Pi and this processor Pj
						tag = 148;
					
					Insert_In_Busy(x2,y2,source);		// make (x2,y2) busy
					Insert_In_EUL(source, (EdgeUpdate){x2, y2, iy2, x1, y1, 3});	// (x2,y2) will get deleted
				}
				// else, y1 is located in this processor
				else if(Create_Parallel_Edge(mapTable.at(y1), x2) || Potential_Parallel_Edge(y1,x2)) // checking for parallel edge (y1,x2)
					continue;
				else
				{
					Insert_In_Busy(x2,y2,source);		// make (x2,y2) busy
					Insert_In_Pending(y1,x2);	// make the edge (y1,x2) a potential edge
					Insert_In_EUL(source, (EdgeUpdate){x2, y2, iy2, x1, y1, 4});	// (x2,y2) will get replaced by (y1,x2)
				}
			}
		}
		
		else	// (x1,y1) and (x2,y2) will get replaced by [(x1,x2) or (x2,x1)] and [(y1,y2) or (y2,y1)]
		{
			LI miny, maxy;
			
			if(y1 < y2)
			{
				miny = y1;
				maxy = y2;
			}
			else
			{
				miny = y2;
				maxy = y1;
			}
			
			k = Hash_Function(miny);
			
			#ifdef DEBUG_DETAIL
				printf("\tProc %d: miny = %ld, maxy = %ld, node %ld is located in proc %ld\n", myRank, miny, maxy, miny, k);
			#endif
			
			if(x2 < x1)		// (x2,x1) at this processor Pj
			{
				if(Create_Parallel_Edge(mx2, x1) || Potential_Parallel_Edge(x2,x1))
					continue;
				
				tag = 28;
				
				if(k != myRank)
				{		
					Insert_In_Busy(x2,y2,source);		// make (x2,y2) busy
					Insert_In_Pending(x2,x1);	// make the edge (x2,x1) a potential edge
					Insert_In_EUL(source, (EdgeUpdate){x2, y2, iy2, x1, y1, 22});	// (x2,y2) will get replaced by (x2,x1)
				}
				else	// else, min(y1,y2) is located in this processor
				{
					LI mminy = mapTable.at(miny);
					if(Create_Parallel_Edge(mminy, maxy) || Potential_Parallel_Edge(miny,maxy)) // checking for parallel edge (y1,y2) or (y2,y1)
						continue;
					else
					{
						Replace_u1_v1_by_u1_v2(mx2, iy2, x1);	// (x2,y2) getting replaced by (x2,x1) - successful swap
						Add_Edge_u_v(mminy, maxy);				// Add (y1,y2) or (y2,y1)
						tag = 106;	// both the edges are updated at this processor
					}
				}
			}
			else	// (x1, x2) at processor Pi
			{
				tag = 38;
				
				if(k != myRank)
				{	
					Insert_In_Busy(x2,y2,source);		// make (x2,y2) busy
					Insert_In_EUL(source, (EdgeUpdate){x2, y2, iy2, x1, y1, 3});	// (x2,y2) will be removed
				}
				else	// else, min(y1,y2) is located in this processor
				{
					LI mminy = mapTable.at(miny);
					if(Create_Parallel_Edge(mminy, maxy) || Potential_Parallel_Edge(miny,maxy)) // checking for parallel edge (y1,y2) or (y2,y1)
						continue;
					else
					{
						tag = 36;
						Insert_In_Busy(x2,y2,source);			// make (x2,y2) busy
						Insert_In_Pending(miny,maxy);
						Insert_In_EUL(source, (EdgeUpdate){x2, y2, iy2, x1, y1, 32});	// (x2,y2) will be replaced by (miny,maxy)
					}
				}
			}
		}
		
		break;
	}
	
	if(nj == MAX_ATTEMPT)
	{
		cout << "\n\nMaximum number of failed attempt made to select the 2nd edge to serve msg with tag = 1 at processor " << myRank << "\n";
		cout << "\nTotal edge switching done at processor " << myRank << " is: " << ni << "\n";
		cout << "\nNumber of local edges at processor " << myRank << " is " << lm << "\n\n";
		Abnormally_Terminate();
	}
}


template <typename DType>	// DType is the data type for weight
void WeightedGraph<DType>::Send_Step_Complete_Msg()
{
	#ifdef DEBUG_DETAIL
		cout << "\n  ********************************************\n  Proc " << myRank << " sending TERMINATING signal to processor ";
	#endif
		
	for(int a = 0; a < noProcs; a++)
		if(a != myRank)
		{
			MPI_Send(send, discardSize, MPI_LONG, a, 228, MPI_COMM_WORLD);
			
			#ifdef DEBUG_DETAIL
				cout << a << ", ";
			#endif
		}
	
	#ifdef DEBUG_DETAIL
		cout << "\n  ********************************************\n";
	#endif
}


template <typename DType>	// DType is the data type for weight
bool WeightedGraph<DType>::Serve_Recv_Messages()
{
	LI  a, b, x1, y1, x2, y2, mx2, my2, iy2, k, nx1, ny1, nx2, ny2;
	EdgeUpdate temp;
	
	#ifdef DEBUG_DETAIL
		if(status.MPI_TAG == 28 || status.MPI_TAG == 38 || status.MPI_TAG == 108 || status.MPI_TAG == 148)
			printf("Proc %d received (%ld, %ld, %ld, %ld) from proc %d with tag = %d\n", myRank, recv[0], recv[1], recv[2], recv[3], status.MPI_SOURCE, status.MPI_TAG);
		else if(status.MPI_TAG == 115 || status.MPI_TAG == 133)
			printf("Proc %d received (%ld, %ld, %ld) from proc %d with tag = %d\n", myRank, recv[0], recv[1], recv[2], status.MPI_SOURCE, status.MPI_TAG);			
		else
			printf("Proc %d received (%ld, %ld) from proc %d with tag = %d\n", myRank, recv[0], recv[1], status.MPI_SOURCE, status.MPI_TAG);			
	#endif

	// receive message with tag = 101
	// processor Pj (this processor) receives the edge (x1,y1) from processor Pi
	// it picks the edge (x2,y2) from its own memory so that no parallel edge will be created
	// edge switching between edge (x1,y1) and (x2,y2)
	// processor Pj will now send the (x2,y2) to processor Pi or to an intermediate processor Pk
	// processor Pj will update the edge (x2,y2) after getting acknowledgement from processor Pi or Pk so that no parallel edge will be created
	
	if(status.MPI_TAG == 101)
	{				
		x1 = recv[0];
		y1 = recv[1];
		
		Select_2nd_Edge_For_Msg_With_Tag_101(x1, y1, x2, y2, mx2, my2, iy2, k);
		
		if(tag == 108 || tag == 28)
		{
			send[0] = y2;
			send[1] = x1;
			send[2] = y1;
			send[3] = status.MPI_SOURCE;
			MPI_Send(send, largeSize, MPI_LONG, k, tag, MPI_COMM_WORLD);
			
			#ifdef DEBUG_DETAIL
				printf("Proc %d is sending (%ld, %ld, %ld, %ld) to proc %ld with tag = %d\n", myRank, send[0], send[1], send[2], send[3], k, tag);
			#endif
		}
		else if(tag == 148 || tag == 38)
		{
			send[0] = y1;
			send[1] = x2;
			send[2] = y2;
			send[3] = status.MPI_SOURCE;
			MPI_Send(send, largeSize, MPI_LONG, k, tag, MPI_COMM_WORLD);
			
			#ifdef DEBUG_DETAIL
				printf("Proc %d is sending (%ld, %ld, %ld, %ld) to proc %ld with tag = %d\n", myRank, send[0], send[1], send[2], send[3], k, tag);
			#endif
		}
		else	// tag = 102 or 106 or 143 or 36
		{
			send[0] = x2;
			send[1] = y2;
			MPI_Send(send, smallSize, MPI_LONG, status.MPI_SOURCE, tag, MPI_COMM_WORLD);
			#ifdef DEBUG_DETAIL
				printf("Proc %d is sending (%ld, %ld) to proc %d with tag = %d\n", myRank, send[0], send[1], status.MPI_SOURCE, tag);
			#endif
		}
		
		#ifdef DEBUG_DETAIL
			printf("Proc %d served msg with tag = %ld\n\n", myRank, (LI)status.MPI_TAG);
		#endif
	}
	
	else if(status.MPI_TAG == 102)	// one edge update at Pi and other one at Pj
	{
		p2  = recv[0];
		q2  = recv[1];
		
		send[0] = p1;
		send[1] = q1;
		
		if(p1 < q2)
		{
			if(Create_Parallel_Edge(mp1, q2) || Potential_Parallel_Edge(p1, q2))	// check for parallel edge possibility for (p1,q2)
				tag = 113;
			else
			{
				tag = 105;
				Replace_u1_v1_by_u1_v2(mp1, iq1, q2);				// (p1,q1) getting replaced by (p1,q2) - successful swap
			}
		}
		else 
		{	
			mq2 = mapTable.at(q2);
			if(Create_Parallel_Edge(mq2, p1) || Potential_Parallel_Edge(q2, p1))	// check for parallel edge possibility for (q2,p1)
				tag = 113;
			else
			{
				tag = 105;
				Replace_u1_v1_by_v2_u1(p1, iq1, mp1, mq2);	// (p1,q1) getting replaced by (q2,p1) - successful swap
			}
		}
		
		MPI_Send(send, smallSize, MPI_LONG, Pj, tag, MPI_COMM_WORLD);	// send acknowledgement/discard message
		Remove_From_Busy(p1,q1);
		
		#ifdef DEBUG_DETAIL
			printf("Proc %d is sending (%ld, %ld) to proc %ld with tag = %d\n", myRank, send[0], send[1], Pj, tag);
			printf("Proc %d served msg with tag = %ld\n", myRank, (LI)status.MPI_TAG);
		#endif
		
		#ifdef DEBUG_DETAIL
			fprintf(tempFP, "Proc %d is sending (%ld, %ld) to proc %ld with tag = %d\n", myRank, send[0], send[1], Pj, tag);
			fprintf(tempFP, "Proc %d served msg with tag = %ld\n", myRank, (LI)status.MPI_TAG);
		#endif
		
		return true;
	}
	
	else if(status.MPI_TAG == 105)	// successful swap - finish the pending works at Pj
	{	
		Get_Edge_From_EUL(status.MPI_SOURCE, recv[0], recv[1], temp);
		Serve_Msg_With_Tag_105(temp.dx1, temp.dy1, temp.dx2, temp.dy2, mapTable.at(temp.dx2), temp.diy2, temp.type);
		#ifdef DEBUG_DETAIL
			printf("Proc %d served msg with tag = %ld\n", myRank, (LI)status.MPI_TAG);
		#endif
	}
	
	else if(status.MPI_TAG == 115)	// successful swap - finish the pending works at Pj
	{	
		Get_Edge_From_EUL(recv[2], recv[0], recv[1], temp);
		Serve_Msg_With_Tag_105(temp.dx1, temp.dy1, temp.dx2, temp.dy2, mapTable.at(temp.dx2), temp.diy2, temp.type);
		#ifdef DEBUG_DETAIL
			printf("Proc %d served msg with tag = %ld\n", myRank, (LI)status.MPI_TAG);
		#endif
	}
	
	else if(status.MPI_TAG == 106)	// both the edges were updated at the other processor Pj
	{
		Remove_Edge_u_v(mp1,iq1);
		Remove_From_Busy(p1,q1);
		#ifdef DEBUG_DETAIL
			printf("Proc %d served msg with tag = %ld\n", myRank, (LI)status.MPI_TAG);
		#endif
		return true;
	}
	
	else if(status.MPI_TAG == 153)	// discard at Pi
	{
		Remove_From_Busy(p1,q1);
		tag = 113;
		#ifdef DEBUG_DETAIL
			printf("Proc %d served msg with tag = %ld\n", myRank, (LI)status.MPI_TAG);
		#endif
		return true;
	}
	
	else if(status.MPI_TAG == 108)	// checking for possibility of new edge (y2,x1) in a processor Pk
	{
		ny2 		  = recv[0];
		nx1 = send[0] = recv[1];
		ny1 = send[1] = recv[2];
		Pi  = send[2] = recv[3];
		
		LI mny2 = mapTable.at(ny2);
		
		if(Create_Parallel_Edge(mny2, nx1) || Potential_Parallel_Edge(ny2,nx1))
		{
			MPI_Send(send, discardSize, MPI_LONG, Pi, 153, MPI_COMM_WORLD);			// Pk sends discard to Pi
			MPI_Send(send, mediumSize, MPI_LONG, status.MPI_SOURCE, 133, MPI_COMM_WORLD);	// Pk sends discard to Pj
			#ifdef DEBUG_DETAIL
				printf("Proc %d is sending a discard msg with tag = 153 to proc %ld\n", myRank, Pi);
				printf("Proc %d is sending a discard msg (%ld, %ld, %ld) with tag = 133 to proc %d\n", myRank, send[0], send[1], send[2], status.MPI_SOURCE);
			#endif
		}
		else
		{
			Add_Edge_u_v(mny2,nx1);
			MPI_Send(send, discardSize, MPI_LONG, Pi, 106, MPI_COMM_WORLD);			// success, send positive acknowledgement
			MPI_Send(send, mediumSize, MPI_LONG, status.MPI_SOURCE, 115, MPI_COMM_WORLD);	// send acknowledgement message
			#ifdef DEBUG_DETAIL
				printf("Proc %d is sending an ack. msg with tag = 106 to proc %ld\n", myRank, Pi);
				printf("Proc %d is sending an ack. msg (%ld, %ld, %ld) with tag = 115 to proc %d\n", myRank, send[0], send[1], send[2], status.MPI_SOURCE);
			#endif
		}
		#ifdef DEBUG_DETAIL
			printf("Proc %d served msg with tag = %ld\n", myRank, (LI)status.MPI_TAG);
		#endif
	}
	
	else if(status.MPI_TAG == 113)	// discard in processor Pj
	{
		Get_Edge_From_EUL(status.MPI_SOURCE, recv[0], recv[1], temp);
		
		if(temp.type == 2)
			Remove_From_Pending(temp.dx2, temp.dy1);
		else if(temp.type == 4)
			Remove_From_Pending(temp.dy1, temp.dx2);
		else if(temp.type == 22)
			Remove_From_Pending(temp.dx2, temp.dx1);
		else if(temp.type == 32)
		{
			if(temp.dy1 < temp.dy2)
				Remove_From_Pending(temp.dy1, temp.dy2);
			else
				Remove_From_Pending(temp.dy2, temp.dy1);
		}
		else if(temp.type != 3)
		{
			cout << "\n\nERROR ! ERROR !! ERROR !!!\n\nWhile DISCARDING Msg with tag = 13, edge extracted from eul array with type = " << temp.type << " at processor " << myRank << "\n\n";
		}
	
		Remove_From_Busy(temp.dx2, temp.dy2);
		#ifdef DEBUG_DETAIL
			printf("Proc %d served msg with tag = %ld\n", myRank, (LI)status.MPI_TAG);
		#endif
	}
	
	else if(status.MPI_TAG == 133)	// discard in processor Pj
	{
		Get_Edge_From_EUL(recv[2], recv[0], recv[1], temp);
		
		if(temp.type == 2)
			Remove_From_Pending(temp.dx2, temp.dy1);
		else if(temp.type == 4)
			Remove_From_Pending(temp.dy1, temp.dx2);
		else if(temp.type == 22)
			Remove_From_Pending(temp.dx2, temp.dx1);
		else if(temp.type == 32)
		{
			if(temp.dy1 < temp.dy2)
				Remove_From_Pending(temp.dy1, temp.dy2);
			else
				Remove_From_Pending(temp.dy2, temp.dy1);
		}
		else if(temp.type != 3)
		{
			cout << "\n\nERROR ! ERROR !! ERROR !!!\n\nWhile DISCARDING Msg with tag = 13, edge extracted from eul array with type = " << temp.type << " at processor " << myRank << "\n\n";
		}
	
		Remove_From_Busy(temp.dx2, temp.dy2);
		#ifdef DEBUG_DETAIL
			printf("Proc %d served msg with tag = %ld\n", myRank, (LI)status.MPI_TAG);
		#endif
	}
	
	else if(status.MPI_TAG == 114)	// discard at proc Pi
	{
		send[0] = p1;
		send[1] = q1;
		
		tag = 113;
		MPI_Send(send, smallSize, MPI_LONG, Pj, tag, MPI_COMM_WORLD);	// send discard message to processor Pj
		Remove_From_Busy(p1,q1);
		#ifdef DEBUG_DETAIL
			printf("Proc %d is sending discard msg to proc %ld with tag = 113\n", myRank, Pj);
			printf("Proc %d served msg with tag = %ld\n", myRank, (LI)status.MPI_TAG);
		#endif
		return true;
	}
	
	else if(status.MPI_TAG == 143)	// both edge update at Pi
	{
		p2 = recv[0];
		q2 = recv[1];
		
		send[0] = p1;
		send[1] = q1;
		
		// check for parallel edge possibility for (p1,q2) and (q1, p2)
		if(Create_Parallel_Edge(mp1, q2) || Potential_Parallel_Edge(p1, q2) || Create_Parallel_Edge(mq1, p2) || Potential_Parallel_Edge(q1, p2))	
			tag = 113;
		else
		{
			tag = 105;								// successful swap
			Replace_u1_v1_by_u1_v2(mp1, iq1, q2);	// (p1,q1) getting replaced by (p1,q2) 
			Add_Edge_u_v(mq1, p2);					// add new edge (q1,p2)
		}
	
		MPI_Send(send, smallSize, MPI_LONG, Pj, tag, MPI_COMM_WORLD);	// send acknowledgement/discard message
		
		Remove_From_Busy(p1,q1);
		#ifdef DEBUG_DETAIL
			printf("Proc %d is sending msg (%ld, %ld) with tag = %d to proc %ld\n", myRank, send[0], send[1], tag, Pj);
			printf("Proc %d served msg with tag = %ld\n", myRank, (LI)status.MPI_TAG);
		#endif
		return true;
	}
	
	else if(status.MPI_TAG == 145)	// one edge may be updated at this processor Pi and other one at Pk
	{
		p2 = recv[0];
		q2 = recv[1];
		
		if(Create_Parallel_Edge(mp1,q2) || Potential_Parallel_Edge(p1,q2))	// discard
		{
			send[0] = q1;
			send[1] = p2;
			MPI_Send(send, smallSize, MPI_LONG, status.MPI_SOURCE, 123, MPI_COMM_WORLD);	// send discard to intermediate processor Pk
			
			#ifdef DEBUG_DETAIL
				printf("Proc %d is sending discard msg (%ld, %ld) with tag = 123 to proc %d\n", myRank, send[0], send[1], status.MPI_SOURCE);
			#endif
			
			send[0] = p1;
			send[1] = q1;
			tag = 113;
			MPI_Send(send, smallSize, MPI_LONG, Pj, tag, MPI_COMM_WORLD);					// send discard to other processor Pj
			#ifdef DEBUG_DETAIL
				printf("Proc %d is sending discard msg (%ld, %ld) with tag = 113 to proc %ld\n", myRank, send[0], send[1], Pj);
			#endif
		}
		else	// successful
		{
			Replace_u1_v1_by_u1_v2(mp1, iq1, q2);											// (p1,q1) getting replaced by (p1,q2)
			
			send[0] = q1;
			send[1] = p2;
			MPI_Send(send, smallSize, MPI_LONG, status.MPI_SOURCE, 125, MPI_COMM_WORLD);	// send acknowledgement to intermediate processor Pk
			#ifdef DEBUG_DETAIL
				printf("Proc %d is sending ack. msg (%ld, %ld) with tag = 125 to proc %d\n", myRank, send[0], send[1], status.MPI_SOURCE);
			#endif
			
			send[0] = p1;
			send[1] = q1;
			tag = 105;
			MPI_Send(send, smallSize, MPI_LONG, Pj, tag, MPI_COMM_WORLD);					// send acknowledgement to other processor Pj
			#ifdef DEBUG_DETAIL
				printf("Proc %d is sending ack. msg (%ld, %ld) with tag = 105 to proc %ld\n", myRank, send[0], send[1], Pj);
			#endif
		}
		
		Remove_From_Busy(p1,q1);
		#ifdef DEBUG_DETAIL
			printf("Proc %d served msg with tag = %ld\n", myRank, (LI)status.MPI_TAG);
		#endif
		return true;
	}
	
	else if(status.MPI_TAG == 148)	// checking for possibility of new edge (y1,x2) in a processor Pk
	{
		ny1 		  = recv[0];
		nx2 = send[0] = recv[1];
		ny2 = send[1] = recv[2];
		Pi  		  = recv[3];
		
		if(Create_Parallel_Edge(mapTable.at(ny1), nx2) || Potential_Parallel_Edge(ny1,nx2))
		{
			MPI_Send(send, discardSize, MPI_LONG, Pi, 114, MPI_COMM_WORLD);		// DISCARD bcz adding new edge (y1,x2) can create parallel edge, send discard msg to Pi
			#ifdef DEBUG_DETAIL
				printf("Proc %d is sending discard msg with tag = 114 to proc %ld\n", myRank, Pi);
			#endif
		}
		else
		{
			Insert_In_Pending(ny1,nx2);
			MPI_Send(send, smallSize, MPI_LONG, Pi, 145, MPI_COMM_WORLD);		// success, send positive acknowledgement
			#ifdef DEBUG_DETAIL
				printf("Proc %d is sending msg (%ld, %ld) with tag = 145 to proc %ld\n", myRank, send[0], send[1], Pi);
			#endif
		}
		#ifdef DEBUG_DETAIL
			printf("Proc %d served msg with tag = %ld\n", myRank, (LI)status.MPI_TAG);
		#endif
	}
	
	else if(status.MPI_TAG == 123)	// discard at intermediate processor Pk
	{
		Remove_From_Pending(recv[0], recv[1]);
		#ifdef DEBUG_DETAIL
			printf("Proc %d served msg with tag = %ld\n", myRank, (LI)status.MPI_TAG);
		#endif
	}
	
	else if(status.MPI_TAG == 125)	// discard at intermediate processor Pk
	{
		Add_Edge_u_v(mapTable.at(recv[0]), recv[1]);
		Remove_From_Pending(recv[0], recv[1]);
		#ifdef DEBUG_DETAIL
			printf("Proc %d served msg with tag = %ld\n", myRank, (LI)status.MPI_TAG);
		#endif
	}
	
	else if(status.MPI_TAG == 28)	// checking for possibility of new edge (y1,y2) or (y2,y1) in some processor Pk or in Pi
	{
		ny2 		  = recv[0];
		nx1 = send[0] = recv[1];
		ny1 = send[1] = recv[2];
		Pi  = send[2] = recv[3];
		
		LI miny, maxy, mminy;
		
		if(ny1 < ny2)
		{
			miny = ny1;
			maxy = ny2;
		}
		else
		{
			miny = ny2;
			maxy = ny1;
		}
		
		mminy = mapTable.at(miny);
		
		#ifdef DEBUG_DETAIL
			printf("\tProc %d: miny = %ld, maxy = %ld, mminy = %ld\n", myRank, miny, maxy, mminy);
		#endif
		
		if(Create_Parallel_Edge(mminy, maxy) || Potential_Parallel_Edge(miny, maxy))	// discard
		{
			if(myRank == Pi)
			{
				Remove_From_Busy(p1,q1);
				tag = 113;
				MPI_Send(send, smallSize, MPI_LONG, Pj, tag, MPI_COMM_WORLD);	// Pi sends discard to Pj
				#ifdef DEBUG_DETAIL
					printf("\tProc %d is sending discard msg (%ld, %ld) with tag = 113 to proc %ld\n", myRank, send[0], send[1], Pj);
				#endif
				return true;
			}
			else
			{
				MPI_Send(send, discardSize, MPI_LONG, Pi, 153, MPI_COMM_WORLD);			// Pk sends discard to Pi
				MPI_Send(send, mediumSize, MPI_LONG, status.MPI_SOURCE, 133, MPI_COMM_WORLD);	// Pk sends discard to Pj
				#ifdef DEBUG_DETAIL
					printf("\tProc %d is sending discard msg with tag = 153 to proc %ld\n", myRank, Pi);
					printf("\tProc %d is sending discard msg (%ld, %ld, %ld) with tag = 133 to proc %d\n", myRank, send[0], send[1], send[2], status.MPI_SOURCE);
				#endif
			}			
		}
		else	// success
		{
			if(myRank == Pi)
			{
				Replace_u1_v1_by_v2_u1(maxy, iq1, mp1, mminy);	// (p1,q1) getting replaced by (mminy,maxy) - successful swap
				Remove_From_Busy(p1,q1);
				tag = 105;
				MPI_Send(send, smallSize, MPI_LONG, Pj, tag, MPI_COMM_WORLD);	// Pi sends ack to Pj
				return true;
			}
			else
			{
				Add_Edge_u_v(mminy,maxy);
				MPI_Send(send, discardSize, MPI_LONG, Pi, 106, MPI_COMM_WORLD);			// Pk sends ack to Pi; Pi will remove (x1,y1)
				MPI_Send(send, mediumSize, MPI_LONG, status.MPI_SOURCE, 115, MPI_COMM_WORLD);	// Pk sends ack to Pj; Pj will replace (x2,y2) by (x2,x1)
				#ifdef DEBUG_DETAIL
					printf("\tProc %d is sending ack. msg with tag = 106 to proc %ld\n", myRank, Pi);
					printf("\tProc %d is sending ack. msg (%ld, %ld, %ld) with tag = 115 to proc %d\n", myRank, send[0], send[1], send[2], status.MPI_SOURCE);
				#endif
			}
		}
		
		#ifdef DEBUG_DETAIL
			printf("Proc %d served msg with tag = 28\n", myRank);
		#endif
	}
	
	else if(status.MPI_TAG == 38)	// checking for possibility of new edge (y1,y2) or (y2,y1) in some processor Pk or in Pi
	{
		ny1 		  = recv[0];
		nx2 = send[0] = recv[1];
		ny2 = send[1] = recv[2];
		Pi  = send[2] = recv[3];
		
		LI miny, maxy, mminy;
		
		if(ny1 < ny2)
		{
			miny = ny1;
			maxy = ny2;
		}
		else
		{
			miny = ny2;
			maxy = ny1;
		}
		
		mminy = mapTable.at(miny);
		
		#ifdef DEBUG_DETAIL
			printf("\tProc %d: miny = %ld, maxy = %ld, mminy = %ld\n", myRank, miny, maxy, mminy);
		#endif
		
		if(Create_Parallel_Edge(mminy, maxy) || Potential_Parallel_Edge(miny, maxy))	// discard
		{
			if(myRank == Pi)
			{
				send[0] = p1;
				send[1] = q1;
				tag = 113;
				MPI_Send(send, smallSize, MPI_LONG, Pj, tag, MPI_COMM_WORLD);	// Pi sends discard to Pj
				Remove_From_Busy(p1,q1);
				#ifdef DEBUG_DETAIL
					printf("\tProc %d is sending discard msg (%ld, %ld) with tag = 113 to proc %ld\n", myRank, send[0], send[1], Pj);
				#endif
				return true;
			}
			else
			{
				MPI_Send(send, discardSize, MPI_LONG, Pi, 114, MPI_COMM_WORLD);			// Pk sends discard to Pi		
				#ifdef DEBUG_DETAIL
					printf("\tProc %d is sending discard msg with tag = 114 to proc %ld\n", myRank, Pi);
				#endif
			}
		}
		else	// success
		{
			if(myRank == Pi)
			{
				if(Create_Parallel_Edge(mp1,nx2) || Potential_Parallel_Edge(p1,nx2))
					tag = 113;	// discard
				else
				{
					Add_Edge_u_v(mminy,maxy);				// add (q1,q2) or (q2,q1)
					Replace_u1_v1_by_u1_v2(mp1, iq1, nx2);	// (p1,q1) getting replaced by (p1,p2) 
					tag = 105;
				}
				
				send[0] = p1;
				send[1] = q1;
				MPI_Send(send, smallSize, MPI_LONG, Pj, tag, MPI_COMM_WORLD);	// Pi sends ack to Pj
				Remove_From_Busy(p1,q1);
				#ifdef DEBUG_DETAIL
					printf("\tProc %d is sending msg (%ld, %ld) with tag = %d to proc %ld\n", myRank, send[0], send[1], tag, Pj);
				#endif
				return true;
			}
			else
			{
				Insert_In_Pending(miny,maxy);
				MPI_Send(send, smallSize, MPI_LONG, Pi, 32, MPI_COMM_WORLD);
				#ifdef DEBUG_DETAIL
					printf("\tProc %d is sending msg (%ld, %ld) with tag = 32 to proc %ld\n", myRank, send[0], send[1], Pj);
				#endif
			}
		}
		
		#ifdef DEBUG_DETAIL
			printf("Proc %d served msg with tag = 38\n", myRank);
		#endif
	}
	
	else if(status.MPI_TAG == 32)	// at processor Pi
	{
		p2 = recv[0];
		q2 = recv[1];
		
		send[0] = Min(q1,q2);
		send[1] = Max(q1,q2);
		
		if(Create_Parallel_Edge(mp1,p2) || Potential_Parallel_Edge(p1,p2))
		{
			MPI_Send(send, smallSize, MPI_LONG, status.MPI_SOURCE, 123, MPI_COMM_WORLD);	// send discard msg to intermediate processor Pk
			#ifdef DEBUG_DETAIL
				printf("\tProc %d is sending msg (%ld, %ld) with tag = 123 to proc %d\n", myRank, send[0], send[1], status.MPI_SOURCE);
			#endif
			
			send[0] = p1;
			send[1] = q1;
			tag = 113;	
			MPI_Send(send, smallSize, MPI_LONG, Pj, tag, MPI_COMM_WORLD);	// send discard to other processor Pj
			#ifdef DEBUG_DETAIL
				printf("\tProc %d is sending discard msg (%ld, %ld) with tag = %d to proc %ld\n", myRank, send[0], send[1], tag, Pj);
			#endif
		}
		else
		{
			Replace_u1_v1_by_u1_v2(mp1, iq1, p2);	// replace (p1,q1) by (p1,p2)
			MPI_Send(send, smallSize, MPI_LONG, status.MPI_SOURCE, 125, MPI_COMM_WORLD);	// send acknowledgement at intermediate processor Pk
			#ifdef DEBUG_DETAIL
				printf("\tProc %d is sending ack. msg (%ld, %ld) with tag = 125 to proc %d\n", myRank, send[0], send[1], status.MPI_SOURCE);
			#endif
			
			send[0] = p1;
			send[1] = q1;
			tag = 105;	
			MPI_Send(send, smallSize, MPI_LONG, Pj, tag, MPI_COMM_WORLD);	// send ack. at other processor Pj
			#ifdef DEBUG_DETAIL
				printf("\tProc %d is sending ack. msg (%ld, %ld) with tag = %d to proc %ld\n", myRank, send[0], send[1], tag, Pj);
			#endif
		}
		
		Remove_From_Busy(p1,q1);
		
		#ifdef DEBUG_DETAIL
			printf("Proc %d served msg with tag = 32\n", myRank);
		#endif
		
		return true;
	}
	
	else if(status.MPI_TAG == 36)	// at processor Pi
	{
		p2 = recv[0];
		q2 = recv[1];
		
		send[0] = p1;
		send[1] = q1;
		
		if(Create_Parallel_Edge(mp1,p2) || Potential_Parallel_Edge(p1,p2))
			tag = 113;
		else
		{
			tag = 105;
			Replace_u1_v1_by_u1_v2(mp1, iq1, p2);	// replace (p1,q1) by (p1,p2)
		}
		
		Remove_From_Busy(p1,q1);
		MPI_Send(send, smallSize, MPI_LONG, Pj, tag, MPI_COMM_WORLD);	// send ack. or discard to other processor Pj
		
		#ifdef DEBUG_DETAIL
			printf("\tProc %d is sending msg (%ld, %ld) with tag = %d to proc %ld\n", myRank, send[0], send[1], tag, Pj);
			printf("Proc %d served msg with tag = 36\n", myRank);
		#endif
		
		return true;
	}	
	
	else if(status.MPI_TAG == 201)	// edge add request (y1,y2) or (y2,y1) for a local edge switch in some other processor <status.MPI_SOURCE>
	{
		a = recv[0];
		b = recv[1];
		
		LI ma = mapTable.at(a);
		
		#ifdef DEBUG_DETAIL
			printf("Proc %d: node %ld is mapped to nlist[%ld]\n", myRank, a, ma);
		#endif
		
		if(Can_Add_Edge(a, ma, b, 0) == false)	// add new edge (a,b)
		{
			MPI_Send(recv, discardSize, MPI_LONG, status.MPI_SOURCE, 206, MPI_COMM_WORLD);	// send discard signal
			
			#ifdef DEBUG_DETAIL
				cout << "\tProc " << myRank << " sending discard signal to processor " << status.MPI_SOURCE << " with tag = 206 for the edge (" << a << "," << b << ")\n";
			#endif
		}
		else
		{
			Add_Edge_u_v(ma, b);	// add the new edge (a,b)
			MPI_Send(recv, discardSize, MPI_LONG, status.MPI_SOURCE, 205, MPI_COMM_WORLD);	// send ack signal
			
			#ifdef DEBUG_DETAIL
				cout << "\tProc " << myRank << " sending ack signal to processor " << status.MPI_SOURCE << " with tag = 205 for the edge (" << a << "," << b << ")\n";
			#endif
		}
						
		#ifdef DEBUG_DETAIL
			cout << "\tProc " << myRank << " served msg with tag = 201\n";
		#endif
	}
				
	else if(status.MPI_TAG == 205)	// ack of local edge switch request : add request was to add (y1,y2) or (y2,y1)
	{
		if(p1 < p2)
		{
			Replace_u1_v1_by_u1_v2(mp1, iq1, p2);	// replacing edge (p1,q1) by new edge (p1,p2)
			Remove_Edge_u_v(mp2,iq2);				// delete (p2,q2)
			Remove_From_Pending(p1,p2);				// remove (p1,p2) from pending
		}
		else
		{
			Replace_u1_v1_by_u1_v2(mp2, iq2, p1);	// replacing edge (p2,q2) by new edge (p2,p1)
			Remove_Edge_u_v(mp1,iq1);				// delete (p1,q1)
			Remove_From_Pending(p2,p1);				// remove (p2,p1) from pending
		}
		
		Remove_From_Busy(p1,q1);
		Remove_From_Busy(p2,q2);
		return true;
	}
		
	else if(status.MPI_TAG == 206)	// discard of local edge switch request : add request was to add (y1,y2) or (y2,y1)
	{
		if(p1 < p2)
			Remove_From_Pending(p1,p2);
		else
			Remove_From_Pending(p2,p1);
				
		Remove_From_Busy(p1,q1);
		Remove_From_Busy(p2,q2);
		tag = 113;
		return true;
	}
	
	else if(status.MPI_TAG == 211)	// edge add request (a,b) for a local edge switch in some other processor <status.MPI_SOURCE>
	{
		a = recv[0];
		b = recv[1];
		
		LI ma = mapTable.at(a);
		
		if(Can_Add_Edge(a, ma, b, 0) == false)	// parallel edge (a,b)?
		{
			MPI_Send(recv, discardSize, MPI_LONG, status.MPI_SOURCE, 216, MPI_COMM_WORLD);	// send discard signal
			
			#ifdef DEBUG_DETAIL
				cout << "\tProc " << myRank << " sending discard signal to processor " << status.MPI_SOURCE << " with tag = 216 for the edge (" << a << "," << b << ")\n";
			#endif
		}
		else
		{
			Add_Edge_u_v(ma, b);	// add the new edge (a,b)
			MPI_Send(recv, discardSize, MPI_LONG, status.MPI_SOURCE, 215, MPI_COMM_WORLD);	// send ack signal
			
			#ifdef DEBUG_DETAIL
				cout << "\tProc " << myRank << " sending ack signal to processor " << status.MPI_SOURCE << " with tag = 215 for the edge (" << a << "," << b << ")\n";
			#endif
		}
		
		#ifdef DEBUG_DETAIL
			printf("\tProc %d served msg with tag = %d\n", myRank, status.MPI_TAG);
		#endif
	}
	
	else if(status.MPI_TAG == 215)	// ack of local edge switch request : add request was to add (q1,p2)
	{
		Replace_u1_v1_by_u1_v2(mp1,iq1,q2);		// replacing edge (p1,q1) by new edge (p1,q2)
		Remove_Edge_u_v(mp2,iq2);				// delete (p2,q2)
		Remove_From_Pending(p1,q2);				// remove (p1,q2) from pending
		Remove_From_Busy(p1,q1);
		Remove_From_Busy(p2,q2);
		#ifdef DEBUG_DETAIL
			printf("\tProc %d served msg with tag = %d\n", myRank, status.MPI_TAG);
		#endif
		return true;
	}
	
	else if(status.MPI_TAG == 216)	// discard of local edge switch request : add request was to add (y1,y2) or (y2,y1)
	{
		Remove_From_Pending(p1,q2);
		Remove_From_Busy(p1,q1);
		Remove_From_Busy(p2,q2);
		tag = 113;
		#ifdef DEBUG_DETAIL
			printf("\tProc %d served msg with tag = %d\n", myRank, status.MPI_TAG);
		#endif
		return true;
	}
	
	else if(status.MPI_TAG == 221)	// edge add request (y1,y2) or (y2,y1) for a local edge switch in some other processor <status.MPI_SOURCE>
	{
		a = recv[0];
		b = recv[1];
		
		LI ma = mapTable.at(a);
		
		if(Can_Add_Edge(a, ma, b, 0) == false)	// add new edge (a,b)
		{
			MPI_Send(recv, discardSize, MPI_LONG, status.MPI_SOURCE, 226, MPI_COMM_WORLD);	// send discard signal
			
			#ifdef DEBUG_DETAIL
				cout << "\tProc " << myRank << " sending discard signal to processor " << status.MPI_SOURCE << " with tag = 226 for the edge (" << a << "," << b << ")\n";
			#endif
		}
		else
		{
			Add_Edge_u_v(ma, b);	// add the new edge (a,b)
			MPI_Send(recv, discardSize, MPI_LONG, status.MPI_SOURCE, 225, MPI_COMM_WORLD);	// send ack signal
			
			#ifdef DEBUG_DETAIL
				cout << "\tProc " << myRank << " sending ack signal to processor " << status.MPI_SOURCE << " with tag = 225 for the edge (" << a << "," << b << ")\n";
			#endif
		}
		
		#ifdef DEBUG_DETAIL
			printf("\tProc %d served msg with tag = %d\n", myRank, status.MPI_TAG);
		#endif
	}
	
	else if(status.MPI_TAG == 225)	// ack of local edge switch request : add request was to add (q1,p2)
	{
		Replace_u1_v1_by_u1_v2(mp2, iq2, q1);	// replacing edge (p2,q2) by new edge (p2,q1)
		Remove_Edge_u_v(mp1,iq1);				// delete (p1,q1)
		Remove_From_Pending(p2,q1);				// remove (p2,q1) from pending
		Remove_From_Busy(p1,q1);
		Remove_From_Busy(p2,q2);
		#ifdef DEBUG_DETAIL
			printf("\tProc %d served msg with tag = %d\n", myRank, status.MPI_TAG);
		#endif
		return true;
	}
	
	else if(status.MPI_TAG == 226)	// discard of local edge switch request : add request was to add (y1,y2) or (y2,y1)
	{
		Remove_From_Pending(p2,q1);
		Remove_From_Busy(p1,q1);
		Remove_From_Busy(p2,q2);
		tag = 113;
		#ifdef DEBUG_DETAIL
			printf("\tProc %d served msg with tag = %d\n", myRank, status.MPI_TAG);
		#endif
		return true;
	}
	
	else if(status.MPI_TAG == 228) // some other processor finished all of its' edge switches! so, increment the counter
	{	
		stepCompMsgCount++;
		#ifdef DEBUG_DETAIL
			printf("\tProc %d served msg with tag = %d\n", myRank, status.MPI_TAG);
		#endif
	}
	
	else if(status.MPI_TAG == 271)
	{
		//MPI_Abort(MPI_COMM_WORLD, 271);
		Abnormally_Terminate();
	}
	
	return false;
}


template <typename DType>	// DType is the data type for weight
void WeightedGraph<DType>::Shuffle()
{
	stepCompMsgCount = tag = 0;
	for(ni = 0; ni < lt; ni++)
	{
		Pj = BRSearch(cq, 0, noProcs-1, (double)rand()/RAND_MAX);
		
		if(Pj == myRank)
			Local_Shuffle();
		else
			Global_Shuffle();
			
		if(tag == 113)
		{
			ni--;
			tag = 0;
		}
	}
	
	Send_Step_Complete_Msg();
	
	while(stepCompMsgCount != noProcs-1)
	{
		MPI_Recv(recv, largeSize, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		Serve_Recv_Messages();
	}
}	


template <typename DType>	// DType is the data type for weight
void WeightedGraph<DType>::Local_Shuffle()
{
	LI minq, maxq, mminq;
	
	for(nk = 0; nk < MAX_ATTEMPT; nk++)
	{
		Select_Edge(p1, q1, mp1, mq1, iq1);	// p1, q1 -> real ids, the edge is (p1,q1)
											// mp1 -> mapped id of p1
											// mq1 -> mapped id of q1 (only if it belongs to this processor)
											// iq1 -> index position of q1 in p1's adjacency list
			
		Select_Edge(p2, q2, mp2, mq2, iq2);	// p2, q2 -> real ids, the edge is (p2,q2)
											// mp2 -> mapped id of p2
											// mq2 -> mapped id of q2 (only if it belongs to this processor)
											// iq2 -> index position of q2 in p2's adjacency list
		
		if(Edge_Busy(p1,q1) || Edge_Busy(p2,q2))
			continue;
		
		#ifdef DEBUG_DETAIL
			printf("Local shuffle #%ld: toss = %ld : selected edges are (%ld, %ld) and (%ld, %ld) at Proc %d\n", ni, toss, p1, q1, p2, q2, myRank);
		#endif
		
		// one edge is (p1, q1) and the other one is (p2, q2)
		// loop and useless edge switch checking
		if(Create_Loop_Useless_Switch_Local())
		{
			#ifdef DEBUG_DETAIL
				printf("\tProc %d: Discard because of possibility of loop or useless edge switch ... select another pair\n", myRank);
			#endif
			continue;
		}
		
		// parallel edge checking : possibility of [(p1,q2) or (q2,p1)] and [(p2,q1) or (q1,p2)]
		if(toss == P1_Q2)	// (p1,q1) will get replaced by (p1,q2) or (q2,p1) and (p2,q2) will get replaced by (p2,q1) or (q1,p2)
		{
			if(p1 < q2)
			{
				if(Can_Add_Edge(p1, mp1, q2, 0) == false)		// check for parallel edge (p1,q2)
					continue;
				
				if(p2 < q1)
				{
					if(Can_Add_Edge(p2, mp2, q1, 0) == false)	// check for parallel edge (p2,q1)
						continue;
				
					#ifdef DEBUG_DETAIL
						printf("Local shuffle #%ld: toss = %ld : taking place between (%ld, %ld) and (%ld, %ld) ar Proc %d\n", ni, toss, p1, q1, p2, q2, myRank);
					#endif
					
					Replace_u1_v1_by_u1_v2(mp1, iq1, q2);	// (p1,q1) is getting replaced by (p1,q2)
					Replace_u1_v1_by_u1_v2(mp2, iq2, q1);	// (p2,q2) is getting replaced by (p2,q1)
				}
				else
				{
					
					Pk = Hash_Function(q1);
					
					if(Pk == myRank)								// (q1,p2) is in this processor
					{
						if(Can_Add_Edge(q1, mq1, p2, 0) == false)	// check for parallel edge (q1,p2)
							continue;
						
						Replace_u1_v1_by_u1_v2(mp1, iq1, q2);		// (p1,q1) is getting replaced by (p1,q2)
						Replace_u1_v1_by_v2_u1(p2, iq2, mp2, mq1);	// (p2,q2) is getting replaced by (q1,p2)
					}
					else	// (q1,p2) is in processor Pk
					{
						Insert_In_Busy(p1,q1,-1);
						Insert_In_Busy(p2,q2,-2);
						
						Insert_In_Pending(p1,q2);
						
						send[0] = q1;
						send[1] = p2;
						
						MPI_Send(send, smallSize, MPI_LONG, Pk, 211, MPI_COMM_WORLD);
						
						#ifdef DEBUG_DETAIL
							printf("Proc %d is sending edge add request (%ld, %ld) to proc %ld with tag = 211\n", myRank, send[0], send[1], Pk);
						#endif
						
						bool  temp  = false;
						while(temp == false)
						{
							MPI_Recv(recv, largeSize, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
							temp = Serve_Recv_Messages();
						}
					}
				}
			}
			else	// q2 < p1, p2 < q2 < p1 < q1
			{		// new edges are (q2,p1) and (p2,q1)
				if(Can_Add_Edge(p2, mp2, q1, 0) == false)	// check for parallel edge (p2,q1)
					continue;
				
				Pk = Hash_Function(q2);
				
				if(Pk == myRank)
				{
					if(Can_Add_Edge(q2, mq2, p1, 0) == false)	// check for parallel edge (q2,p1)
						continue;
					
					Replace_u1_v1_by_v2_u1(p1, iq1, mp1, mq2);	// replacing edge (p1,q1) by new edge (q2,p1)
					Replace_u1_v1_by_u1_v2(mp2, iq2, q1);		// (p2,q2) is getting replaced by (p2,q1)
				}
				else
				{
					Insert_In_Busy(p1,q1,-1);
					Insert_In_Busy(p2,q2,-2);
					
					Insert_In_Pending(p2,q1);
					
					send[0] = q2;
					send[1] = p1;
					
					MPI_Send(send, smallSize, MPI_LONG, Pk, 221, MPI_COMM_WORLD);
						
					#ifdef DEBUG_DETAIL
						printf("\tProc %d is sending (%ld, %ld) to proc %ld with tag = 221\n", myRank, send[0], send[1], Pk);
					#endif
					
					bool  temp  = false;
					while(temp == false)
					{
						MPI_Recv(recv, largeSize, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
						temp = Serve_Recv_Messages();
					}
				}
			}
 		}
		else	// P1_P2
		{
			
			if(Can_Add_Edge(p1, mp1, p2, mp2) == false)	// check for parallel edge (p1,p2) or (p2,p1)			
				continue;
			
			// the second edge would be (minq, maxq)
			if(q1 < q2)
			{
				minq  = q1;
				maxq  = q2;
				mminq = mq1;
			}
			else
			{
				minq  = q2;
				maxq  = q1;
				mminq = mq2;
			}
			
			Pk = Hash_Function(minq);
			
			if(Pk == myRank)	// (minq,maxq) is in this processor
			{
			
				if(Can_Add_Edge(minq, mminq, maxq, 0) == false)	// check for parallel edge (minq,maxq) [(q1,q2) or (q2,q1)]
					continue;
				
				if(p1 < p2)
				{
					Replace_u1_v1_by_u1_v2(mp1, iq1, p2);	// Replacing edge (p1,q1) by new edge (p1,p2)
					Remove_Edge_u_v(mp2,iq2);				// delete (p2,q2)
				}
				else
				{
					Replace_u1_v1_by_u1_v2(mp2, iq2, p1);	// Replacing edge (p2,q2) by new edge (p2,p1)
					Remove_Edge_u_v(mp1,iq1);				// delete (p1,q1)
				}
				
				Add_Edge_u_v(mminq, maxq);	// add edge (q1,q2) or (q2,q1)
			}
			else	// (minq, maxq) is in processor Pk
			{
				Insert_In_Busy(p1,q1,-1);
				Insert_In_Busy(p2,q2,-2);
				
				if(p1 < p2)
					Insert_In_Pending(p1,p2);
				else
					Insert_In_Pending(p2,p1);
				
				send[0] = minq;
				send[1] = maxq;
				
				MPI_Send(send, smallSize, MPI_LONG, Pk, 201, MPI_COMM_WORLD);
				
				#ifdef DEBUG_DETAIL
					printf("\tProc %d is sending (%ld, %ld) to proc %ld with tag = 201\n", myRank, send[0], send[1], Pk);
				#endif
				
				bool  temp  = false;
				while(temp == false)
				{
					MPI_Recv(recv, largeSize, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
					temp = Serve_Recv_Messages();
				}
			}
		}
		
		break;
	}
	
	if(nk == MAX_ATTEMPT)
	{
		cout << "\n\nMAXIMUM number of attempt made for a local edge switch at proc " << myRank << "\n";
		cout << "\nTotal edge switching done at processor " << myRank << " is: " << ni << "\n";
		cout << "\nNumber of local edges at processor " << myRank << " is " << lm << "\n\n";
		Abnormally_Terminate();
	}
}


template <typename DType>	// DType is the data type for weight
void WeightedGraph<DType>::Global_Shuffle()
{
	// Edge switching going on between this processor (processor <myRank> or Pi) and processor Pj	

	for(nk = 0; nk < MAX_ATTEMPT; nk++)
	{
		#ifdef DEBUG_DETAIL
			printf("\nGlobal shuffle #%ld between processor %d and %ld\n\n", ni, myRank, Pj);
		#endif
		
		Select_Edge(p1, q1, mp1, mq1, iq1);	// p1,q1 -> real ids, the edge is (p1,q1)
											//   mp1 -> mapped id of p1
											//   mq1 -> mapped id of q1 (only if it belongs to this processor)
											//   iq1 -> index position of q1 in p1's adjacency list											
		
		if(Edge_Busy(p1,q1) == true)		// (p1,q1) is busy, select another edge
			continue;
		
		Insert_In_Busy(p1,q1,-1);				// make (p1,q1) busy
			
		send[0] = p1;
		send[1] = q1;
		
		tag = 101;
		
		MPI_Send(send, smallSize, MPI_LONG, Pj, tag, MPI_COMM_WORLD);
		
		#ifdef DEBUG_DETAIL
			printf("Global shuffle #%ld: proc %d sending (%ld, %ld) to proc %ld with tag = %d\n", ni, myRank, send[0], send[1], Pj, tag);
		#endif
		
		bool  temp  = false;
		while(temp == false)
		{
			MPI_Recv(recv, largeSize, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			temp = Serve_Recv_Messages();
		}
		
		break;
	}
	
	if(nk == MAX_ATTEMPT)
	{
		cout << "\n\nMAXIMUM number of attempt made for a global edge switch at proc " << myRank << "\n";
		cout << "\nTotal edge switching done at processor " << myRank << " is: " << ni << "\n";
		cout << "\nNumber of local edges at processor " << myRank << " is " << lm << "\n\n";
		Abnormally_Terminate();
	}
}

#endif /* #ifndef SHUFFLE_H_ */