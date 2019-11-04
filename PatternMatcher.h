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
#include <stdio.h>

class PatternMatcher
{
private:
	std::vector<char> text;
	std::vector<char> pattern;

	std::vector<unsigned char> codedText;
	std::vector<unsigned char> codedPattern;

	std::vector<unsigned char> encode(std::vector<char> text);
	std::vector<unsigned char> encode2(std::vector<char> text);
	std::vector<unsigned char> encode3(std::vector<char> text);
	void encodeParallel3(std::vector<char> text);

	char getCodedValue(char base);

	friend class MatchingTester;

public:
	PatternMatcher() 
	{	
		text = std::vector<char>();
		pattern = std::vector<char>();
	}
	PatternMatcher(std::string textFileName, std::string patternFileName)
	{
		this->loadText(textFileName);
		this->loadPattern(patternFileName);
	}

	void loadText(std::string textFileName);
	void loadPattern(std::string patternFileName);
	void loadTextChunk(std::string textFileName);

	inline std::size_t textLength() { return text.size(); }
	inline std::size_t patternLength() { return pattern.size(); }

	std::vector<int> naive();
	std::vector<int> boyer_moore();
	std::vector<int> kmp();
	std::vector<int> smith_waterman();
	std::vector<int> needleman_wunsch();

	std::vector<int> naiveOpenMP();
	std::vector<int> naiveParallel();
	std::vector<int> naiveParallelOpenMP();

	std::vector<int> boyer_mooreParallel();

	static std::vector<int> naive(std::vector<char> text, std::vector<char> pattern);
	static std::vector<int> boyer_moore(std::vector<char> text, std::vector<char> pattern);
	static std::vector<int> naiveOpenMP(std::vector<char> text, std::vector<char> pattern);
	static std::vector<int> smithWaterman(std::vector<char> text, std::vector<char> pattern);

	std::vector<int> coded_naive();
	std::vector<int> coded_boyer_moore();

	std::vector<int> coded_naiveOpenMP();
	std::vector<int> coded_naiveParallel();
	std::vector<int> coded_naiveParallelOpenMP();

	std::vector<int> coded_boyer_mooreOpenMP();
	std::vector<int> coded_boyer_mooreParallel();
	std::vector<int> coded_boyer_mooreParallelOpenMP();

	static std::vector<int> coded_naive(std::vector<unsigned char> text, std::vector<unsigned char> pattern);
	static std::vector<int> coded_naiveOpenMP(std::vector<unsigned char> text, std::vector<unsigned char> pattern);
	static std::vector<int> coded_boyer_moore(std::vector<unsigned char> text, std::vector<unsigned char> pattern);
	static std::vector<int> coded_boyer_mooreOpenMP(std::vector<unsigned char> codedText, std::vector<unsigned char> codedPattern);

};
