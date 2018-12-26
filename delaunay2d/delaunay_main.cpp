#include "delaunay_2d.h"
#include <iostream>
#include <chrono>
#include <fstream>

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

void save_as_obj(std::string path, const Eigen::MatrixX2d& Pts2d, const Eigen::MatrixX3i& Tri)
{
	std::ofstream file(path);
	if (file)
	{
		for (int i = 0; i < Pts2d.rows(); ++i)
		{
			file << "v " << Pts2d(i, 0) << " " << Pts2d(i, 1) << " 0" << std::endl;
		}
		for (int i = 0; i < Tri.rows(); ++i)
		{
			file << "f " << Tri(i, 0) + 1 << " " << Tri(i, 1) + 1 << " " << Tri(i, 2) + 1 << std::endl;
		}
	}
}

int main()
{
	Eigen::MatrixXd pts;
	pts.setRandom(300000, 2);

	Timer<> timer;
	auto res = Delaunay2D(pts);
	timer.EndTimer("TimerEnd: ");

	std::cout << res.GetTri().rows() << std::endl;

	//save_as_obj("test.obj", pts, res.GetTri());
	//std::ofstream file("out.txt");
	//file << res.GetMatlabCmd() << std::endl;
	//std::cout << res.GetMatlabCmd() << std::endl;

/*
	Eigen::MatrixXd out1;
	Eigen::VectorXi out2;
	std::cout << pts << std::endl << std::endl;
	sort(pts, col, 1, out1, out2);
	std::cout << out1 << std::endl;
	std::cout << out2 << std::endl;*/

	return 0;
}