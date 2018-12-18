#include "delaunay_2d.h"
#include <array>
#include <iostream>
#include <algorithm>
#include <string>
#include <map>
#include <set>
#include <deque>
#include <functional>
#include <numeric>

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

// 面片相邻
std::vector<std::vector<int>> SolveAdjTri(const Eigen::MatrixX3i& Tri)
{
	// 假设每边最多两个三角形共享
	std::vector<std::vector<int>> result(Tri.rows());
	std::map<std::pair<int, int>, int> edge_set;

	std::array<int, 6> _ind_for_sorted_tri{ 0,1,0,2,1,2 };

	for (int i = 0; i < Tri.rows(); ++i)
	{
		std::array<int, 3> tri_now_sorted{ Tri(i,0), Tri(i,1) ,Tri(i,2) };
		std::sort(tri_now_sorted.begin(), tri_now_sorted.end());
		for (int j = 0; j < 3; ++j)
		{
			auto key = std::make_pair(tri_now_sorted[_ind_for_sorted_tri[2 * j]], tri_now_sorted[_ind_for_sorted_tri[2 * j + 1]]);
			if (edge_set.count(key))
			{
				result[i].push_back(edge_set[key]);
				result[edge_set[key]].push_back(i);
			}
			else
			{
				edge_set[key] = i;
			}
		}
	}

	//std::vector<Eigen::VectorXi> ret(Tri.rows());
	//for (int i = 0; i < ret.size(); ++i)
	//{
	//	ret[i].resize(result[i].size());
	//	for (int j = 0; j < result[i].size(); ++j)
	//	{
	//		ret[i][j] = result[i][j];
	//	}
	//}
	return result;
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
	MatlabCmd += "figure;res = delaunayTriangulation(p);triplot(res.ConnectivityList, p(:,1), p(:,2));title('official algorithm from matlab');axis equal;grid on;\n";
#endif
	//删除辅助三角形
	FinalTri.erase(std::remove_if(FinalTri.begin(), FinalTri.end(), [n](const decltype(FinalTri)::value_type &t) {	return (t.array() >= n).any();	}), FinalTri.end());

	this->Tri.resize(FinalTri.size(), 3);
	for (int i = 0; i < FinalTri.size(); ++i)
	{
		for (int j = 0; j < Tri.cols(); ++j)
			this->Tri.row(i)[j] = out_index(FinalTri[i][j]);
	}

	SameOrientation();

#ifdef GenMatlabCmd
	MatlabCmd += "figure; hold on; axis equal; grid on;title('final layout')\n";
	MatlabCmd += "p_ori=[";
	for (int i = 0; i < n; ++i)
	{
		MatlabCmd += std::to_string(PtsToInsert.row(i)[0]) + ", " + std::to_string(PtsToInsert.row(i)[1]) + ";\n";
	}
	MatlabCmd += "];\nplot(p_ori(:,1),p_ori(:,2),'o');\n"; 
	MatlabCmd += "tri=[";
	for (int i = 0; i < this->Tri.rows(); ++i)
	{
		MatlabCmd += std::to_string(1 + this->Tri.row(i)[0]) + ", " + std::to_string(1 + this->Tri.row(i)[1]) + ", " + std::to_string(1 + this->Tri.row(i)[2]) + ";\n";
	}
	MatlabCmd += "];\n";
	MatlabCmd += "for i = 1:size(p_ori, 1)\ntext(p_ori(i, 1), p_ori(i, 2), num2str(i));\nend\n";
	MatlabCmd += "for i = 1:size(tri, 1)\n";
	MatlabCmd += "myline(p_ori(tri(i, 1), :), p_ori(tri(i, 2), :), 'Color', [0, 0, 1]);\n"; 
	MatlabCmd += "myline(p_ori(tri(i, 1), :), p_ori(tri(i, 3), :), 'Color', [0, 0, 1]);\n";
	MatlabCmd += "myline(p_ori(tri(i, 2), :), p_ori(tri(i, 3), :), 'Color', [0, 0, 1]);\nend\n";
#endif
}

// 其实做成图遍历并调整就行
void Delaunay2D::SameOrientation()
{
	enum check_status { not_visited, visited_need_flip, visited_normal };
	std::vector<check_status> flip_status(Tri.rows(), not_visited);
	auto Adj = SolveAdjTri(Tri);
	std::deque<int> visiting_queue;
	// intial condition
	visiting_queue.push_back(0);
	flip_status[0] = visited_normal; // 固定该面片为正常
	// visiting all triangles
	while (!visiting_queue.empty())
	{
		int now_visiting = visiting_queue.front();
		visiting_queue.pop_front();

		for (int i = 0; i < Adj[now_visiting].size(); ++i)
		{
			int now_comparing = Adj[now_visiting][i];

			// update this by checking edge order, 如果一样就是需要反向
			std::array<int, 2> now_visiting_ind, now_comparing_ind;
			int now_visiting_ind_pos = 0, now_comparing_ind_pos = 0;
			for (int j = 0; j < 3; ++j)
				for (int k = 0; k < 3; ++k)
					if(Tri(now_visiting, j) == Tri(now_comparing, k))
					{
						now_visiting_ind[now_visiting_ind_pos++] = j;
						now_comparing_ind[now_comparing_ind_pos++] = k;
					}
			check_status for_update;
			if ((now_visiting_ind[1] - now_visiting_ind[0] + 3) % 3 == (now_comparing_ind[1] - now_comparing_ind[0] + 3) % 3)
			{
				if(flip_status[now_visiting] == visited_normal)
					for_update = visited_need_flip;
				else
					for_update = visited_normal;
			}
			else
			{
				for_update = flip_status[now_visiting];
			}

			if (flip_status[now_comparing] == not_visited)
			{
				flip_status[now_comparing] = for_update;
				visiting_queue.push_back(now_comparing);
			}
			else
			{
				if (for_update != flip_status[now_comparing])
					error("may be Mobius, check your input triangle.");//对本例不会出现这种问题
			}
		}
	}
	// flip all that need to flip
	for (int i = 0; i < Tri.rows(); ++i)
	{
		if (flip_status[i] == visited_need_flip)
		{
			using std::swap;
			swap(Tri(i, 1), Tri(i, 2));
		}
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
