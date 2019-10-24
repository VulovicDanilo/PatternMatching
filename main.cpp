
#include <iostream>
#include "PatternMatcher.h"
#include "MatchingTester.h"
#include "mpi.h";

using namespace std;
int main(int argc, char* argv[])
{
	string genome = "../../files/1mb.fa";
	string pattern = "../../files/pattern_8.fa";
	PatternMatcher pm(genome, pattern);

	MPI_Init(&argc, &argv);

	MatchingTester::testNaive(pm);

	MatchingTester::testNaiveOpenMP(pm);

	MatchingTester::testNaiveParallel(pm);

	MatchingTester::testNaiveParallelOpenMP(pm);

	MatchingTester::testSmithWaterman(pm);

	MPI_Finalize();
	return 0;
}