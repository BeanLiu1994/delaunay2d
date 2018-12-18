#include "delaunay_2d.h"
#include <iostream>

int main()
{
	Eigen::MatrixXd pts;
	pts.setRandom(15000, 2);

	auto res = Delaunay2D(pts);

	//std::cout << res.GetTri() << std::endl;
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