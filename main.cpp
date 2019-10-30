
#include <iostream>
#include "PatternMatcher.h"
#include "MatchingTester.h"
#include "mpi.h";

using namespace std;
int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	string hi;
	hi = "hi";
	string text;
	text = "./files/1mb.fa";
	string pattern;
	pattern = "./files/pattern_8.fa";
	PatternMatcher pm;
	if (!rank)
	{
		cout << "Loading files..." << endl;
		// pm.loadText(text);
		pm.loadTextChunk(text);
		pm.loadPattern(pattern);
		cout << "Files Loaded. Matching started..." << endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);


	MatchingTester::testNaive(pm);

	// MatchingTester::testCodedNaiveParallelOpenMP(pm);

	// MatchingTester::testCodedBoyerMooreOpenMP(pm);

	MatchingTester::testCodedBoyerMooreParallelOpenMP(pm);

	if (!rank)
	{
		cout << "Matching over." << endl;
	}
	MPI_Finalize();
	return 0;
}