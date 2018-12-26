#pragma once
#include <chrono>

template<typename _ty = std::chrono::microseconds>
class Timer
{
public:
	Timer(bool startNow = true)
	{
		if (startNow) StartTimer();
	}
	void StartTimer(const std::string& _printThis = std::string())
	{
		if (!_printThis.empty()) std::cout << _printThis << std::endl;
		start = std::chrono::system_clock::now();
	}
	double EndTimer(const std::string& _printThis, bool restartFlag)
	{
		auto seconds = EndTimer(_printThis);
		if (restartFlag)
			StartTimer();
		return seconds;
	}
	double EndTimer(const std::string& _printThis = std::string())
	{
		end = std::chrono::system_clock::now();
		duration = std::chrono::duration_cast<_ty>(end - start);
		seconds = (double(duration.count()) * _ty::period::num / _ty::period::den);
		std::string s = " elapsed time:  " + std::to_string(seconds) + "s";
		if (!_printThis.empty()) std::cout << _printThis << s << std::endl;
		//else std::cout << s << std::endl;
		return seconds;
	}
private:
	std::chrono::time_point<std::chrono::system_clock> start, end;
	_ty duration;
	double seconds;
};