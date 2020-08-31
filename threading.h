#pragma once

#include <iostream>
#include <chrono>
#include <thread>
#include <mutex>

//timing struct for printing its lifetime (measuring time of its existence)
struct Timer
{
	std::string name;
	std::chrono::time_point<std::chrono::system_clock> start;
	std::chrono::duration<float> duration;

	Timer(std::string name) : name(name)
	{
		start = std::chrono::high_resolution_clock::now();
	}

	~Timer()
	{
		duration = std::chrono::high_resolution_clock::now() - start;
		
		std::cout << "Timer " << name << ": " << 1000*duration.count() << " ms\n"; 
	}
};

//simple for loop on multiple threads
template<typename func, typename... args>
double threadProbabilityFor(int numOfThreads, int repeats, func f, args... a)
{
	std::vector<double> outputParts;
	outputParts.reserve(numOfThreads);

	std::mutex output_mutex;

	auto forPart = [&](int begin, int partRepeats){
		
		double outputPart = 0;

		for (int i = begin; i < partRepeats + begin; ++i) outputPart += f(i, a...);

		std::lock_guard<std::mutex> guard(output_mutex);
		outputParts.push_back(outputPart);

	};
	
	int part 	  = repeats / numOfThreads;
	int remainder = repeats % numOfThreads;

	std::vector<std::thread> threads;
  	threads.reserve(numOfThreads);

	for (int i = 0; i < numOfThreads; ++i)
	{
		int extra = i < remainder ? i : remainder;
		int begin = i*part + extra;
		threads.push_back(std::thread(forPart, begin, part + (i < remainder)));
	}

	for (int i = 0; i < numOfThreads; ++i) threads[i].join();

	double output = 0;

	for(double outputPart : outputParts) output += outputPart;

	return output;
}