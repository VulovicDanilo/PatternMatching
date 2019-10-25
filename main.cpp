
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
	string genome = "./files/50mb.fa";
	string pattern = "./files/pattern_128.fa";
	if (!rank) cout << "Loading genomes..." << endl;
	PatternMatcher pm(genome, pattern);

	if (!rank) cout << "Genomes loaded. Matching started..." << endl;

	MatchingTester::testNaive(pm);

	MatchingTester::testCodedNaive(pm);

	MatchingTester::testCodedNaiveOpenMP(pm);

	//MatchingTester::testBoyerMoore(pm);

	//MatchingTester::testKMP(pm);

	//MatchingTester::testNaiveOpenMP(pm);

	// MatchingTester::testNaiveParallel(pm);

	//MatchingTester::testNaiveParallelOpenMP(pm);

	//MatchingTester::testSmithWaterman(pm);


	if (!rank) cout << "Matching over." << endl;
	MPI_Finalize();
	return 0;
}