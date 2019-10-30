#pragma once

#include <chrono>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>

class Timer
{
private:
	std::chrono::time_point<std::chrono::system_clock> startTime;
	std::chrono::time_point<std::chrono::system_clock> finishTime;

public:
	Timer() {}

	void start() { startTime = std::chrono::system_clock::now(); }
	void stop() { finishTime = std::chrono::system_clock::now(); }

	double elapsed() 
	{
		std::chrono::duration<double> duration = (std::chrono::duration<double>) (finishTime - startTime);
		return duration.count();
	}
};