#pragma once

#include <chrono>
#include <string>
#include <vector>
#include <algorithm>

class Timer
{
private:
	std::chrono::time_point<std::chrono::steady_clock> startTime;
	std::chrono::time_point<std::chrono::steady_clock> finishTime;

public:
	Timer() {}

	void start() { startTime = std::chrono::high_resolution_clock::now(); }
	void stop() { finishTime = std::chrono::high_resolution_clock::now(); }

	double elapsed() 
	{
		std::chrono::duration<double> duration = finishTime - startTime;
		return duration.count();
	}
};