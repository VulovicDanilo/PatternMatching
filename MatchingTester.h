#pragma once

#include <string>
#include <vector>
#include "PatternMatcher.h"

class MatchingTester
{
private:

	static bool test(PatternMatcher& pm, std::vector<int> matches);
public:
	MatchingTester() {}

	static bool testNaive(PatternMatcher& pm);
	static bool testNaiveParallel(PatternMatcher& pm);
	static bool testNaiveOpenMP(PatternMatcher& pm);
	static bool testNaiveParallelOpenMP(PatternMatcher& pm);

	static bool testBoyerMoore(PatternMatcher& pm);
	static bool testKMP(PatternMatcher& pm);
	static bool testSmithWaterman(PatternMatcher& pm);

	static bool testCodedNaive(PatternMatcher& pm);
	static bool testCodedNaiveOpenMP(PatternMatcher& pm);
	static bool testCodedNaiveParallel(PatternMatcher& pm);
};