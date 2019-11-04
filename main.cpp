
#include <iostream>
#include "PatternMatcher.h"
#include "MatchingTester.h"
#include "mpi.h"

using namespace std;
int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int rank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	vector<string> genomes =
	{
		"./files/genome_1mb",
		"genome_50mb.fa",
		"genome_100mb.fa",
		"genome_250mb.fa",
		"genome_1gb.fa",
		// "genome_3gb.fa"
	};
	vector<string> patterns =
	{
		"./files/pattern_32.fa",
		"pattern_128.fa",
		"pattern_256.fa",
		"pattern_521.fa"
		"pattern_1024.fa"
	};
	for (string genome : genomes)
	{
		for (string pattern : patterns)
		{
			PatternMatcher pm;
			if (!rank)
			{
				cout << "Files: " << endl;
				cout << genome << endl;
				cout << pattern << endl;
				cout << "Number of nodes: " << size << endl;
				cout << "Loading files..." << endl;
				pm.loadTextChunk(genome);
				pm.loadPattern(pattern);
				cout << "Files Loaded. Matching started..." << endl;
			}
			MPI_Barrier(MPI_COMM_WORLD);

			MatchingTester::testNaive(pm);
			MatchingTester::testNaiveParallel(pm);
			MatchingTester::testNaiveParallelOpenMP(pm);
			MatchingTester::testBoyerMoore(pm);
			MatchingTester::testBoyerMooreParallel(pm);

			MatchingTester::testCodedNaive(pm);
			MatchingTester::testCodedNaiveParallel(pm);
			MatchingTester::testCodedNaiveParallelOpenMP(pm);
			MatchingTester::testCodedBoyerMoore(pm);
			MatchingTester::testCodedBoyerMooreParallel(pm);
			MatchingTester::testCodedBoyerMooreParallelOpenMP(pm);

			if(!rank)
			{
				cout << "Matching for files done." << endl << "==================================================" << endl;
			}
		}
	}
	MPI_Finalize();
	return 0;
}