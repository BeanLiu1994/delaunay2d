//#define DELAUNAY_2D_USE_EIGEN

#include "delaunay_2d.h"
#include "Timer.h"
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>

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

std::vector<Delaunay2D::Triangle> Run_std(const std::vector<Delaunay2D::Point>& pts)
{
	Timer<> timer;
	auto res = Delaunay2D::Delaunay2D(pts);
	timer.EndTimer("TimerEnd: ");

	return res.GetTri();
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

Eigen::MatrixX3i Run_eigen(const Eigen::MatrixX2d& pts)
{
	Timer<> timer;
	auto res = Delaunay2D::Delaunay2D(pts);
	timer.EndTimer("TimerEnd: ");

	return res.GetTri_Eigen();
}
#endif


//random
int test_0()
{
	std::vector<Delaunay2D::Point> pts(70 * 80);
	for (int i = 0; i < pts.size(); ++i)
	{
		pts[i][0] = rand();
		pts[i][1] = rand();
	}

	save_as_obj("random.obj", pts, Run_std(pts));
#ifdef DELAUNAY_2D_USE_EIGEN
	Eigen::MatrixX2d pts_eigen;
	pts_eigen.resize(pts.size(), 2);
	for (int i = 0; i < pts.size(); ++i)
	{
		pts_eigen(i, 0) = pts[i][0];
		pts_eigen(i, 1) = pts[i][1];
	}

	Run_eigen(pts_eigen);
	save_as_obj("random_eigen.obj", pts_eigen, Run_eigen(pts_eigen));
#endif

	//std::ofstream file("out.txt");
	//file << res.GetMatlabCmd() << std::endl;
	//std::cout << res.GetMatlabCmd() << std::endl;
	
	return 0;
}

int test_1()
{
	std::vector<Delaunay2D::Point> pts(20 * 20);
	for (int i = 0; i < 20; ++i)
		for (int j = 0; j < 20; ++j)
		{
			pts[i * 20 + j][0] = i;
			pts[i * 20 + j][1] = j;
		}

	save_as_obj("grid.obj", pts, Run_std(pts));
	return 0;
}

int test_2()
{
	std::vector<Delaunay2D::Point> pts(19 + 19);
	for (int i = 0; i < 19; ++i)
	{
		pts[i][0] = 0;
		pts[i][1] = i - 9;
	}
	for (int i = 19; i < 38; ++i)
	{
		pts[i][0] = i - 18;
		pts[i][1] = 0;
	}

	save_as_obj("axis.obj", pts, Run_std(pts));
	return 0;
}

std::vector<Delaunay2D::Point> generateCircle(int n_points, double theta_initial, const Delaunay2D::Point& center, double radius)
{
	std::vector<Delaunay2D::Point> ret(n_points);
	double delta = 2 * M_PI / n_points;
	double theta = theta_initial;
	for (int i = 0; i < n_points; ++i)
	{
		ret[i][0] = center[0] + radius * std::cos(theta);
		ret[i][1] = center[1] + radius * std::sin(theta);
		theta += delta;
	}
	return ret;
}
int test_3()
{
	std::vector<Delaunay2D::Point> pts;
	auto append = [&pts](const std::vector<Delaunay2D::Point>&rhs) 
	{
		pts.insert(pts.begin(), rhs.begin(), rhs.end());
	};
	append(generateCircle(13, 0, { 0.0,10.0 }, 3));
	append(generateCircle(13, 0.02, { 10.0,10.0 }, 3));
	append(generateCircle(13, 0, { 10.0,0.0 }, 3));
	append(generateCircle(13, 0.05, { 0.0,0.0 }, 3));
	append(generateCircle(13, 0, { -10.0,0.0 }, 1));
	append(generateCircle(13, 0.05, { -10.0,-10.0 }, 0.5));
	append(generateCircle(13, 0.05, { 0.0,-10.0 }, 0.8));
	append(generateCircle(13, 0.05, { -10.0,10.0 }, 1.2));
	append(generateCircle(13, 0.05, { 10.0,-10.0 }, 2.2));

	save_as_obj("circles.obj", pts, Run_std(pts));
	return 0;
}

int main()
{
	//test_0();
	//test_1();
	//test_2();
	test_3();
	return 0;
}