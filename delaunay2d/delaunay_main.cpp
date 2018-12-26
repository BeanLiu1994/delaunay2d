//#define DELAUNAY_2D_USE_EIGEN

#include "delaunay_2d.h"
#include "Timer.h"
#include <iostream>
#include <fstream>


void save_as_obj(std::string path, const std::vector<Delaunay2D::Point>& Pts2d, const std::vector<Delaunay2D::Triangle>& Tri)
{
	std::ofstream file(path);
	if (file)
	{
		for (int i = 0; i < Pts2d.size(); ++i)
		{
			file << "v " << Pts2d[i][0] << " " << Pts2d[i][1] << " 0" << std::endl;
		}
		for (int i = 0; i < Tri.size(); ++i)
		{
			file << "f " << Tri[i][0] + 1 << " " << Tri[i][1] + 1 << " " << Tri[i][2] + 1 << std::endl;
		}
	}
}

void test()
{
	std::vector<Delaunay2D::Point> pts(700 * 800);
	for (int i = 0; i < pts.size(); ++i)
	{
		pts[i][0] = rand();
		pts[i][1] = rand();
	}

	Timer<> timer;
	auto res = Delaunay2D::Delaunay2D(pts);
	timer.EndTimer("TimerEnd: ");

	std::cout << res.GetTri().size() << std::endl;

	save_as_obj("test1.obj", pts, res.GetTri());
}

#ifdef DELAUNAY_2D_USE_EIGEN
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

void test_eigen()
{
	Eigen::MatrixXd pts;
	pts.setRandom(400000, 2);

	Timer<> timer;
	auto res = Delaunay2D::Delaunay2D(pts);
	timer.EndTimer("TimerEnd: ");

	std::cout << res.GetTri_Eigen().rows() << std::endl;

	save_as_obj("test2.obj", pts, res.GetTri_Eigen());
}
#endif



int main()
{
	test();
#ifdef DELAUNAY_2D_USE_EIGEN
	test_eigen();
#endif

	//std::ofstream file("out.txt");
	//file << res.GetMatlabCmd() << std::endl;
	//std::cout << res.GetMatlabCmd() << std::endl;
	
	return 0;
}