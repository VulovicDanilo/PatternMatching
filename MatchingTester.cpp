#pragma once

#include "MatchingTester.h"
#include <iterator>
using namespace std;

bool MatchingTester::test(PatternMatcher& pm, vector<int> matches)
{
	set<int> s(matches.begin(), matches.end());
	matches.assign(s.begin(), s.end());
	bool correct = true;
	for (unsigned int i = 0; (i < matches.size() && correct); i++)
	{
		int index = matches[i];
		if (index <= pm.textLength() - pm.patternLength())
		{
			for (int j = 0; j < pm.patternLength(); j++)
				if (pm.text[index + j] != pm.pattern[j])
					correct = false;
		}
		else
		{
			correct = false;
		}
	}
	if (correct)
	{
		cout << "Found patterns are correct (Matches: " << matches.size() << ")" << endl;
	}
	else
	{
		cout << "Incorrect patterns (Matches: " << matches.size() << ")" << endl;
	}
	return correct;
}

bool MatchingTester::testNaive(PatternMatcher& pm)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		vector<int> matches = pm.naive();
		return MatchingTester::test(pm, matches);
	}
}

bool MatchingTester::testBoyerMoore(PatternMatcher& pm)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		vector<int> matches = pm.boyer_moore();
		return MatchingTester::test(pm, matches);
	}
}

bool MatchingTester::testBoyerMooreParallel(PatternMatcher& pm)
{
	int rank;
	vector<int> matches = pm.boyer_mooreParallel();
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		return MatchingTester::test(pm, matches);
	}
}

bool MatchingTester::testKMP(PatternMatcher& pm)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		vector<int> matches = pm.kmp();
		return MatchingTester::test(pm, matches);
	}
}

bool MatchingTester::testSmithWaterman(PatternMatcher& pm)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		vector<int> matches = pm.smith_waterman();
		return true;
	}
}

bool MatchingTester::testCodedNaive(PatternMatcher& pm)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		vector<int> matches = pm.coded_naive();
		return MatchingTester::test(pm, matches);
	}
}

bool MatchingTester::testCodedNaiveOpenMP(PatternMatcher& pm)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		vector<int> matches = pm.coded_naiveOpenMP();
		return MatchingTester::test(pm, matches);
	}
}

bool MatchingTester::testCodedNaiveParallel(PatternMatcher& pm)
{
	int rank;
	vector<int> matches = pm.coded_naiveParallel();
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		return MatchingTester::test(pm, matches);
	}
}

bool MatchingTester::testCodedNaiveParallelOpenMP(PatternMatcher & pm)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<int> matches = pm.coded_naiveParallelOpenMP();
	if (rank == 0)
	{
		return MatchingTester::test(pm, matches);
	}
}

bool MatchingTester::testCodedBoyerMoore(PatternMatcher& pm)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		vector<int> matches = pm.coded_boyer_moore();
		return MatchingTester::test(pm, matches);
	}
}

bool MatchingTester::testCodedBoyerMooreOpenMP(PatternMatcher & pm)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		vector<int> matches = pm.coded_boyer_mooreOpenMP();
		return MatchingTester::test(pm, matches);
	}
}

bool MatchingTester::testCodedBoyerMooreParallel(PatternMatcher& pm)
{
	int rank;
	vector<int> matches = pm.coded_boyer_mooreParallel();
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		return MatchingTester::test(pm, matches);
	}
}

bool MatchingTester::testCodedBoyerMooreParallelOpenMP(PatternMatcher& pm)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<int> matches = pm.coded_boyer_mooreParallelOpenMP();
	if (rank == 0)
	{
		return MatchingTester::test(pm, matches);
	}
}


bool MatchingTester::testNaiveParallel(PatternMatcher& pm)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<int> matches = pm.naiveParallel();
	if (rank == 0)
	{
		return MatchingTester::test(pm, matches);
	}
}

bool MatchingTester::testNaiveOpenMP(PatternMatcher& pm)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		vector<int> matches = pm.naiveOpenMP();
		return MatchingTester::test(pm, matches);
	}
}

bool MatchingTester::testNaiveParallelOpenMP(PatternMatcher& pm)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<int> matches = pm.naiveParallelOpenMP();
	if (rank == 0)
	{
		return MatchingTester::test(pm, matches);
	}
}
