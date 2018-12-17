#include "delaunay_2d.h"
#include <iostream>

int main()
{
	Eigen::Matrix<double, 50, 2> pts;
	pts.setRandom();

	auto res = Delaunay2D(pts);

	std::cout << res.GetTri() << std::endl;
	std::cout << res.GetMatlabCmd() << std::endl;

	return 0;
}