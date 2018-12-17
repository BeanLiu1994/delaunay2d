#include "delaunay_2d.h"
#include <array>
#include <iostream>
#include <algorithm>
#include <string>

void error(std::string err_info)
{
	std::cerr << err_info << std::endl;
	throw std::runtime_error(err_info);
}


Delaunay2D::Delaunay2D(const Eigen::MatrixX2d& PtsToInsert)
{
	assign(PtsToInsert);
}

void Delaunay2D::assign(const Eigen::MatrixX2d& PtsToInsert)
{
	std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> Tri;

	int n = PtsToInsert.rows();
	if (n < 3)
	{
		error("Delaunay reqiures more than 3 points.");
	}

	MatlabCmd += "myline = @(p1,p2,opt_title,opt_content)line([p1(1);p2(1)], [p1(2);p2(2)], opt_title, opt_content);\n";
	MatlabCmd += "figure; hold on; axis equal; grid on;title('this algorithm')\n";

	Eigen::RowVector2d max_bound = PtsToInsert.colwise().maxCoeff().array();
	Eigen::RowVector2d min_bound = PtsToInsert.colwise().minCoeff().array();

	Eigen::Matrix<double, -1, 2, Eigen::RowMajor> Convert;
	Convert.resize(n + 4, 2);
	Convert.topRows(n) = PtsToInsert;

	// 生成全包三角形(两个)
	Convert.row(n + 0) = Eigen::RowVector2d((max_bound[0] + min_bound[0]) / 2, (3 * max_bound[1] - min_bound[1]) / 2);
	Convert.row(n + 1) = Eigen::RowVector2d((3 * min_bound[0] - max_bound[0]) / 2, (max_bound[1] + min_bound[1]) / 2);
	Convert.row(n + 2) = Eigen::RowVector2d((3 * max_bound[0] - min_bound[0]) / 2, (max_bound[1] + min_bound[1]) / 2);
	Convert.row(n + 3) = Eigen::RowVector2d((max_bound[0] + min_bound[0]) / 2, (3 * min_bound[1] - max_bound[1]) / 2);

	MatlabCmd += "p=[";
	for (int i = 0; i < n + 4; ++i)
	{
		MatlabCmd += std::to_string(Convert.row(i)[0]) + ", " + std::to_string(Convert.row(i)[1]) + ";\n";
	}
	MatlabCmd += "];\nplot(p(1:end-4,1),p(1:end-4,2),'o');\nplot(p(end-3:end,1),p(end-3:end,2),'ro');\n";

	// 遍历所有三角形
	Tri.emplace_back(n + 0, n + 1, n + 2);
	Tri.emplace_back(n + 3, n + 1, n + 2);
	for (int ith_pts = 0; ith_pts < n; ++ith_pts)
	{
		std::vector<int> ToDelete;
		std::vector<Edge> BadEdges;
		for (int ith_tri = 0; ith_tri < Tri.size(); ++ith_tri)
		{
			if (isInsideTriangleCircumCircle(
				Convert.row(ith_pts),
				Convert.row(Tri[ith_tri][0]),
				Convert.row(Tri[ith_tri][1]),
				Convert.row(Tri[ith_tri][2])
			))
			{
				ToDelete.push_back(ith_tri);

				std::array<int, 3> tri_now_sorted;
				std::copy(Tri[ith_tri].data(), Tri[ith_tri].data() + 3, tri_now_sorted.begin());
				std::sort(tri_now_sorted.begin(), tri_now_sorted.end());

				BadEdges.emplace_back(Tri[ith_tri][0], Tri[ith_tri][1]);
				BadEdges.emplace_back(Tri[ith_tri][1], Tri[ith_tri][2]);
				BadEdges.emplace_back(Tri[ith_tri][0], Tri[ith_tri][2]);
			}
		}
		std::sort(ToDelete.begin(), ToDelete.end(), std::greater<int>());
		ToDelete.erase(std::unique(ToDelete.begin(), ToDelete.end()), ToDelete.end());
		for (const auto& m : ToDelete)
			Tri.erase(Tri.begin() + m);

		ToDelete.clear();
		// 删除重复边，得到多边形边界，此时必为一个凸的，而插入点将分割它。
		for (int ith_edge = 0; ith_edge < BadEdges.size() - 1; ++ith_edge)
		{
			for (int cmp_edge = ith_edge + 1; cmp_edge < BadEdges.size(); ++cmp_edge)
			{
				if (BadEdges[ith_edge] == BadEdges[cmp_edge])
				{
					ToDelete.push_back(ith_edge);
					ToDelete.push_back(cmp_edge);
				}
			}
		}
		std::sort(ToDelete.begin(), ToDelete.end(), std::greater<int>());
		ToDelete.erase(std::unique(ToDelete.begin(), ToDelete.end()), ToDelete.end());

		for (const auto& m : ToDelete)
		{
			MatlabCmd += "myline(p(1+" + std::to_string(BadEdges[m].first) + ",:),p(1+" + std::to_string(BadEdges[m].second) + ",:), 'Color', [0.8,0.8,0.8]);\n";
			BadEdges.erase(BadEdges.begin() + m);
		}
		MatlabCmd += "pause(0.05);\n";

		for (const auto& m : BadEdges)
		{
			Tri.emplace_back(m.first, m.second, ith_pts);
			MatlabCmd += "myline(p(1+" + std::to_string(Tri.back()[0]) + ",:),p(1+" + std::to_string(Tri.back()[1]) + ",:), 'Color', [0,0,1]);\n";
			MatlabCmd += "myline(p(1+" + std::to_string(Tri.back()[0]) + ",:),p(1+" + std::to_string(Tri.back()[2]) + ",:), 'Color', [0,0,1]);\n";
			MatlabCmd += "myline(p(1+" + std::to_string(Tri.back()[1]) + ",:),p(1+" + std::to_string(Tri.back()[2]) + ",:), 'Color', [0,0,1]);\n";
		}
		MatlabCmd += "pause(0.05);\n";
	}

	MatlabCmd += "figure;axis equal;grid on;res = delaunayTriangulation(p);triplot(res.ConnectivityList, p(:,1), p(:,2));title('official algorithm from matlab');\n";
	//删除辅助三角形
	Tri.erase(std::remove_if(Tri.begin(), Tri.end(), [n](const decltype(Tri)::value_type &t) {	return (t.array() >= n).any();	}), Tri.end());

	Eigen::MatrixX3i ret;
	ret.resize(Tri.size(), 3);
	for (int i = 0; i < Tri.size(); ++i)
	{
		ret.row(i) = Tri[i];
	}

	this->Tri = ret;
}

Eigen::MatrixX3i Delaunay2D::GetTri() const
{
	return Tri;
}

std::string Delaunay2D::GetMatlabCmd() const
{
	return MatlabCmd;
}

