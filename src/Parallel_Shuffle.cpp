/* 
 * File	  : Parallel_Shuffle.cpp 
 * Purpose: Shuffle edges parallely (MPI program) on arbitrary id Galib format graph
 * Output : Write the output graph in adjacency list (Galib) or edge list format
 * Author : Md Hasanuzzaman Bhuiyan
 * Email  : mhb@vbi.vt.edu
 * Date	  : May 2, 2015
 */

#include "Parallel_Shuffle.hpp"
#include "Serial_Multinomial.hpp"
#include "Timer.hpp"
#include <sys/time.h>
#include "mpi.h"


using namespace std;

typedef  double WType;

struct rec{
    double 	value;
    int   	index;
} in, out;

void MPI_Shuffle(char *inputGraph, double shuffleFraction, char *outputGraph, LI noSteps, int writeFlag, int rank, int noProcs)
{

	/************************************** Initialization and Read Graph *************************************************************/

	double timeTaken;
	LI noSwitchAtThisStep, multAtThisProc;
	WeightedGraph <WType> g;
	
	g.Initialize(shuffleFraction, noSteps, rank, noProcs);
	
	if(g.myRank == 0)
		printf("Reading the input graph: %s\n", inputGraph);
	
	Timer trd(1);
	
	g.ReadGraph(inputGraph);
	
	in.value = trd.getsec();
	in.index = g.myRank;
	
	//if(g.myRank == 0)
		//printf("Processor %d finished reading ... waiting for others to finish their reading ...\n", g.myRank);
	
	/************************************** Print READ time *************************************************************/
	
	// MIN-READ
	MPI_Reduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
	if(g.myRank == 0)
		printf("Reading input graph completed ...\n\nREAD Time:\n==========\n\nMinimum read time = %.2lf seconds at processor = %d\n", out.value, out.index);
	
	// MAX-READ
	MPI_Reduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
	if(g.myRank == 0)
		printf("Maximum read time = %.2lf seconds at processor = %d\n", out.value, out.index);
	
	// AVG-READ
	MPI_Reduce(&in.value, &timeTaken, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(g.myRank == 0)
		printf("Average read time = %.2lf seconds\n", timeTaken/g.noProcs);
	
	/************************************** Allocation and Initialization *************************************************************/
	
	
	LI  	*temp 		= new LI[g.noProcs];
	double 	*prob		= new double[g.noProcs];
	int 	*recvcounts = new int[g.noProcs];
	for(int i=0; i<g.noProcs; i++)
	{
		recvcounts[i] = 1;
		prob[i] = 0.0;
	}
	
	/************************************** Shuffle *************************************************************/
	
	if(g.myRank == 0)
		printf("\nSwitching edges ...\n");
	
	Timer t(1);
	
	for(LI step = 0; step < g.noSteps; step++)
	{
		prob[g.myRank] = g.lm*1.0L/g.m;		// probability value update
		MPI_Allreduce(prob, g.q, g.noProcs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		
		g.cq[0] = g.q[0];					// cumulative probability value update
		for(int i = 1; i < g.noProcs; i++)
			g.cq[i] = g.cq[i-1] + g.q[i];
		
		noSwitchAtThisStep 	= (step < g.noShuffle%g.noSteps)? g.noShuffle/g.noSteps + 1 : g.noShuffle/g.noSteps;
		multAtThisProc 		= (g.myRank < noSwitchAtThisStep%g.noProcs) ? noSwitchAtThisStep/g.noProcs + 1 : noSwitchAtThisStep/g.noProcs;
		
		Multinomial(multAtThisProc, g.q, g.noProcs, temp);	// multinomial distribution
		MPI_Reduce_scatter(temp, &g.lt, recvcounts, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
		g.tt += g.lt;
		g.Shuffle();		// edge switch
	}
	
	in.value = t.getsec();
	
	/************************************** Print SHUFFLE time *************************************************************/
	
	// MIN-SHUFFLE
	MPI_Reduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
	if(g.myRank == 0)
		printf("Switching completed ...\n\nSHUFFLE Time:\n=============\n\nMinimum shuffle time = %.2lf seconds at processor = %d\n", out.value, out.index);
	
	// MAX-SHUFFLE
	MPI_Reduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
	if(g.myRank == 0)
	{
		printf("\n******************************************************************************\n");
		printf("Maximum shuffle time = %.2lf seconds at processor = %d", out.value, out.index);
		printf("\n******************************************************************************\n\n");
	}
	
	// AVG-SHUFFLE
	MPI_Reduce(&in.value, &timeTaken, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(g.myRank == 0)
		printf("Average shuffle time = %.2lf seconds\n", timeTaken/g.noProcs);
	
	/************************************** Write Output Graph *************************************************************/
	
	MPI_Barrier(MPI_COMM_WORLD);	// mandatory barrier ... otherwise, during edge switch some processor may receive message which are supposed to receive during graph write
	
	if(writeFlag != 0)
	{
		if(g.myRank == 0)
			printf("\nWriting output graph ...\n");

		Timer tw(1);
		
		g.WriteGraph(outputGraph, writeFlag);
		
		in.value = tw.getsec();
		
		/************************************** Print Write time *************************************************************/
		
		// MIN-WRITE
		MPI_Reduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
		if(g.myRank == 0)
			cout << "Writing output graph completed ...\n\nWrite Time:\n===========\n\nMinimum write time = " << out.value << " seconds at processor = " << out.index << endl;
		
		// MAX-WRITE
		MPI_Reduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
		if(g.myRank == 0)
			cout << "Maximum write time = " << out.value << " seconds at processor = " << out.index << endl;
		
		// AVG-WRITE
		MPI_Reduce(&in.value, &timeTaken, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if(g.myRank == 0)
			cout << "Average write time = " << timeTaken/g.noProcs << " seconds" << endl;
	}
	
	if(g.myRank == 0)
		cout << "\nTerminating the program successfully ...\nWaiting to finalize ...\n";
	
	FreeMem(temp);
	FreeMem(recvcounts);
	FreeMem(prob);
}


/*
	This program takes total 5 command line arguments as input
			1st arg : input graph file name (in Galib format)
			2nd arg : fraction of shuffle, range (0.0, 1.0], t = 0.5, 1.0 means 50% and 100% edge switch respectively
			3rd arg : output graph file name
			4th arg : number of steps
			5th arg : write the output file?
				 0 means: No
				 1 means: Yes, in Galib format
				 2 means: Yes, in edge list format
*/

int main(int argc, char **argv)
{	
	int rank, noProcs;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &noProcs);
	
	if(argc != 6)
	{
		if(rank == 0)
		{
			cout << "\nWrong command line argument input\n\nUsage: " << argv[0] << " <input_graph_file_name> <fraction-of-shuffle> <output_graph_file_name> <noSteps> <writeFlag>" << endl;
			cout << "Total parameters passed " << argc << endl;
			for(int i = 1; i<=argc; i++)
				cout << argv[i] << endl;
		}
		MPI_Finalize();
	}
	
	srand((unsigned)time(NULL));
	MPI_Shuffle(argv[1], atof(argv[2]), argv[3], atol(argv[4]), atoi(argv[5]), rank, noProcs);
	
	MPI_Finalize();
	return 0;
}
