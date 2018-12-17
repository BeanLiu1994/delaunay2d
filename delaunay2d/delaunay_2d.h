#pragma once

#include <Eigen/Eigen>
#include <string>
#include <utility>

typedef std::pair<int, int> Edge;

class Delaunay2D
{
	std::string MatlabCmd;
	void assign(const Eigen::MatrixX2d& PtsToInsert);
	Eigen::MatrixX3i Tri;
public:
	Delaunay2D(const Eigen::MatrixX2d& PtsToInsert);
	Eigen::MatrixX3i GetTri() const;
	std::string GetMatlabCmd() const;
};

inline bool isInsideTriangleCircumCircle(const Eigen::RowVector2d& pts_now,
	const Eigen::RowVector2d& p0, const Eigen::RowVector2d& p1, const Eigen::RowVector2d& p2)
{
	// ½â·½³Ì
	double p0_v = p0.dot(p0);
	double p1_v = p1.dot(p1);
	double p2_v = p2.dot(p2);

	Eigen::Matrix2d A;
	A.row(0) = p0 - p1;
	A.row(1) = p0 - p2;
	Eigen::Vector2d b(p0_v - p1_v, p0_v - p2_v);
	Eigen::RowVector2d center = (A.transpose() * A).colPivHouseholderQr().solve(A.transpose() * b).array() / 2;

	return (pts_now - center).norm() <= (p1 - center).norm();
}