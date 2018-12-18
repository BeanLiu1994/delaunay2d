#include "delaunay_2d.h"
#include <array>
#include <iostream>
#include <algorithm>
#include <string>
#include <functional>

void error(std::string err_info)
{
	std::cerr << err_info << std::endl;
	throw std::runtime_error(err_info);
}

void RemoveDuplicate(std::vector<int>& ToProcess)
{
	std::sort(ToProcess.begin(), ToProcess.end(), std::greater<int>());
	ToProcess.erase(std::unique(ToProcess.begin(), ToProcess.end()), ToProcess.end());
}

Delaunay2D::Delaunay2D(const Eigen::MatrixX2d& PtsToInsert)
{
	assign(PtsToInsert);
}

void Delaunay2D::assign(const Eigen::MatrixX2d& PtsToInsert)
{
	int n = PtsToInsert.rows();
	if (n < 3)
	{
		error("Delaunay reqiures more than 3 points.");
	}

	Eigen::RowVector2d max_bound = PtsToInsert.colwise().maxCoeff().array();
	Eigen::RowVector2d min_bound = PtsToInsert.colwise().minCoeff().array();

	Eigen::MatrixX2d out_sorted;
	Eigen::VectorXi out_index;
	sort(PtsToInsert, row, 0, out_sorted, out_index);

	Eigen::Matrix<double, -1, 2, Eigen::RowMajor> Convert;
	Convert.resize(n + 4, 2);
	Convert.topRows(n) = out_sorted;

	// 生成全包三角形(两个)
	Convert.row(n + 0) = Eigen::RowVector2d((max_bound[0] + min_bound[0]) / 2, (3 * max_bound[1] - min_bound[1]) / 2);
	Convert.row(n + 1) = Eigen::RowVector2d((3 * min_bound[0] - max_bound[0]) / 2, (max_bound[1] + min_bound[1]) / 2);
	Convert.row(n + 2) = Eigen::RowVector2d((3 * max_bound[0] - min_bound[0]) / 2, (max_bound[1] + min_bound[1]) / 2);
	Convert.row(n + 3) = Eigen::RowVector2d((max_bound[0] + min_bound[0]) / 2, (3 * min_bound[1] - max_bound[1]) / 2);

#ifdef GenMatlabCmd
	MatlabCmd += "myline = @(p1,p2,opt_title,opt_content)line([p1(1);p2(1)], [p1(2);p2(2)], opt_title, opt_content);\n";
	MatlabCmd += "figure; hold on; axis equal; grid on;title('this algorithm')\n";
	MatlabCmd += "p=[";
	for (int i = 0; i < n + 4; ++i)
	{
		MatlabCmd += std::to_string(Convert.row(i)[0]) + ", " + std::to_string(Convert.row(i)[1]) + ";\n";
	}
	MatlabCmd += "];\nplot(p(1:end-4,1),p(1:end-4,2),'o');\nplot(p(end-3:end,1),p(end-3:end,2),'ro');\n";
#endif

	std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> FinalTri, TmpTri;

	TmpTri.emplace_back(n + 0, n + 1, n + 2);
	TmpTri.emplace_back(n + 3, n + 1, n + 2);
	
	// 遍历所有三角形
	for (int ith_pts = 0; ith_pts < n; ++ith_pts)
	{
		std::vector<int> ToDelete;
		std::vector<Edge> BadEdges;
		for (int ith_tri = 0; ith_tri < TmpTri.size(); ++ith_tri)
		{
			auto cmp_result = isInsideTriangleCircumCircle(
				Convert.row(ith_pts),
				Convert.row(TmpTri[ith_tri][0]),
				Convert.row(TmpTri[ith_tri][1]),
				Convert.row(TmpTri[ith_tri][2])
			);
			if (cmp_result == WeakAccept)
			{
				ToDelete.push_back(ith_tri);

				std::array<int, 3> tri_now_sorted;
				std::copy(TmpTri[ith_tri].data(), TmpTri[ith_tri].data() + 3, tri_now_sorted.begin());
				std::sort(tri_now_sorted.begin(), tri_now_sorted.end());

				BadEdges.emplace_back(TmpTri[ith_tri][0], TmpTri[ith_tri][1]);
				BadEdges.emplace_back(TmpTri[ith_tri][1], TmpTri[ith_tri][2]);
				BadEdges.emplace_back(TmpTri[ith_tri][0], TmpTri[ith_tri][2]);
			}
			else if (cmp_result == StrongReject)
			{
				// this is now a stable triangle.
				ToDelete.push_back(ith_tri);
				FinalTri.push_back(TmpTri[ith_tri]);
			}
		}
		RemoveDuplicate(ToDelete);
		for (const auto& m : ToDelete)
			TmpTri.erase(TmpTri.begin() + m);

		ToDelete.clear();
		// 删除重复边，得到多边形边界，此时必为一个凸的，而插入点将分割它。
		for (int ith_edge = 0; ith_edge < (int)BadEdges.size() - 1; ++ith_edge)
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
		RemoveDuplicate(ToDelete);
		for (const auto& m : ToDelete)
		{
#ifdef GenMatlabCmd
			MatlabCmd += "myline(p(1+" + std::to_string(BadEdges[m].first) + ",:),p(1+" + std::to_string(BadEdges[m].second) + ",:), 'Color', [0.8,0.8,0.8]);\n";
#endif
			BadEdges.erase(BadEdges.begin() + m);
		}
#ifdef GenMatlabCmd
		MatlabCmd += "pause(0.05);\n";
#endif

		for (const auto& m : BadEdges)
		{
			TmpTri.emplace_back(m.first, m.second, ith_pts);
#ifdef GenMatlabCmd
			MatlabCmd += "myline(p(1+" + std::to_string(TmpTri.back()[0]) + ",:),p(1+" + std::to_string(TmpTri.back()[1]) + ",:), 'Color', [0,0,1]);\n";
			MatlabCmd += "myline(p(1+" + std::to_string(TmpTri.back()[0]) + ",:),p(1+" + std::to_string(TmpTri.back()[2]) + ",:), 'Color', [0,0,1]);\n";
			MatlabCmd += "myline(p(1+" + std::to_string(TmpTri.back()[1]) + ",:),p(1+" + std::to_string(TmpTri.back()[2]) + ",:), 'Color', [0,0,1]);\n";
#endif
		}
#ifdef GenMatlabCmd
		MatlabCmd += "pause(0.05);\n";
#endif
	}

	FinalTri.insert(FinalTri.end(), TmpTri.begin(), TmpTri.end());
#ifdef GenMatlabCmd
	MatlabCmd += "figure;axis equal;grid on;res = delaunayTriangulation(p);triplot(res.ConnectivityList, p(:,1), p(:,2));title('official algorithm from matlab');\n";
#endif
	//删除辅助三角形
	FinalTri.erase(std::remove_if(FinalTri.begin(), FinalTri.end(), [n](const decltype(FinalTri)::value_type &t) {	return (t.array() >= n).any();	}), FinalTri.end());

	this->Tri.resize(FinalTri.size(), 3);
	for (int i = 0; i < FinalTri.size(); ++i)
	{
		for (int j = 0; j < Tri.cols(); ++j)
			this->Tri.row(i)[j] = out_index(FinalTri[i][j]);
	}

}

Eigen::MatrixX3i Delaunay2D::GetTri() const
{
	return Tri;
}

#ifdef GenMatlabCmd
std::string Delaunay2D::GetMatlabCmd() const
{
	return MatlabCmd;
}
#endif
