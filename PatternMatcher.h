#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include <omp.h>
#include <set>
#include "Timer.h"
#include "CodedPattern.h"

class PatternMatcher
{
private:
	std::vector<char> text;
	std::vector<char> pattern;

	std::vector<unsigned char> encode(std::vector<char> text);
	char getCodedValue(char base);

	friend class MatchingTester;

public:
	PatternMatcher() {	}
	PatternMatcher(std::string textFileName, std::string patternFileName)
	{
		std::ifstream tfs(textFileName), pfs(patternFileName);
		if (!tfs.fail() && !pfs.fail())
		{
			text.assign((std::istreambuf_iterator<char>(tfs)),
				(std::istreambuf_iterator<char>()));
			pattern.assign((std::istreambuf_iterator<char>(pfs)),
				(std::istreambuf_iterator<char>()));

			text.erase(std::remove(text.begin(), text.end(), '\n'), text.end());
			pattern.erase(std::remove(pattern.begin(), pattern.end(), '\n'), pattern.end());

			
		}
		else
		{
			std::cout << "there is a problem with one of file streams" << std::endl;
		}
	}

	inline size_t textLength() { return text.size(); }
	inline size_t patternLength() { return pattern.size(); }

	std::vector<int> naive();
	std::vector<int> boyer_moore();
	std::vector<int> kmp();
	std::vector<int> smith_waterman();
	std::vector<int> needleman_wunsch();

	std::vector<int> naiveOpenMP();
	std::vector<int> naiveParallel();
	std::vector<int> naiveParallelOpenMP();

	static std::vector<int> naive(std::vector<char> text, std::vector<char> pattern);
	static std::vector<int> boyer_moore(std::vector<char> text, std::vector<char> pattern);
	static std::vector<int> naiveOpenMP(std::vector<char> text, std::vector<char> pattern);
	static std::vector<int> smithWaterman(std::vector<char> text, std::vector<char> pattern);

	std::vector<int> coded_naive();

	std::vector<int> coded_naiveOpenMP();
	std::vector<int> coded_naiveParallel();
	std::vector<int> coded_naiveParallelOpenMP();

	static std::vector<int> coded_naive(std::vector<char> text, std::vector<char> pattern);
};
