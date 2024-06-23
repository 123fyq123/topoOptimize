#pragma once
#include<iostream>
#include<chrono>
#include<thread>

class Timer
{
public:
	Timer()
	{
		start = std::chrono::steady_clock::now();
	}
    void out(const std::string& s){
        end = std::chrono::steady_clock::now();
        duration = end - start;
        float ms = duration.count() * 1000.0f;
        std::cout << s + " Timer took " << ms << "ms" << std::endl;
    }
	~Timer()
	{
		end = std::chrono::steady_clock::now();
		duration = end - start;
		float ms = duration.count() * 1000.0f;
		std::cout << "Timer took " << ms << "ms" << std::endl;
	}
private:
	std::chrono::time_point<std::chrono::steady_clock> start, end;
	std::chrono::duration<float> duration;
};

//void Function() {
//	Timer();
//	for (int i = 0; i < 100; ++i) {
//		std::cout << "Do something!" << std::endl;
//	}
//}
//
//int main() {
//
//	Function();
//	std::cin.get();
//}