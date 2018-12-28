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
#define DBCHECK

inline std::string EdgeToMatlabCmd(const std::vector<Delaunay2D::Edge>& ori)
{
	std::string cmd = "\nblist = [";
	for (const auto& m : ori)
	{
		cmd += std::to_string(m.first) + ", " + std::to_string(m.second) + "; ";
	}
	cmd += "]+1;\nfor i = 1:size(blist, 1)\n	myline(p(blist(i, 1), :), p(blist(i, 2), :), 'Color', [1, 0, 0]);\nend";
	return cmd;
}

void warning(std::string err_info)
{
	std::cerr << err_info << std::endl;
}
void error(std::string err_info)
{
	warning(err_info);
	throw std::runtime_error(err_info);
}

void warning(bool judge, std::string err_info)
{
	if (judge) warning(err_info);
}

void error(bool judge, std::string err_info)
{
	if (judge) error(err_info);
}

namespace Delaunay2D
{
	// 面片相邻
	std::vector<std::vector<int>> SolveAdjTri(const std::vector<Triangle>& Tri)
	{
		// 假设每边最多两个三角形共享
		std::vector<std::vector<int>> result(Tri.size());
		std::map<std::pair<int, int>, int> edge_set;

		std::array<int, 6> _ind_for_sorted_tri{ 0,1,0,2,1,2 };

		for (int i = 0; i < Tri.size(); ++i)
		{
			std::array<int, 3> tri_now_sorted{ Tri[i][0], Tri[i][1] ,Tri[i][2] };
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


	bool check_circular(const std::vector<Edge>& BadEdges)
	{
		std::map<int, std::vector<int>> mp;
		for (int i = 0; i < BadEdges.size(); ++i)
		{
			mp[BadEdges[i].first].push_back(BadEdges[i].second);
			mp[BadEdges[i].second].push_back(BadEdges[i].first);
		}
		for (auto iter = mp.begin(); iter != mp.end(); ++iter)
		{
			if (iter->second.size() != 2)
				return false;
		}
		return true;
	}


#ifdef DELAUNAY_2D_USE_EIGEN
	Delaunay2D::Delaunay2D(const Eigen::MatrixX2d& PtsToInsert)
	{
		assign(PtsToInsert);
	}
#endif

	Delaunay2D::Delaunay2D(const std::vector<Point>& PtsToInsert)
	{
		assign(PtsToInsert);
	}

#ifdef DELAUNAY_2D_USE_EIGEN
	void Delaunay2D::assign(const Eigen::MatrixX2d& PtsToInsert)
	{
		auto Pts_Input = MatToStdVec(PtsToInsert);
		assign(Pts_Input);
	}
#endif

	void Delaunay2D::assign(const std::vector<Point>& Pts)
	{
		if (Pts.size() < 3)
		{
			error("Delaunay reqiures more than 3 points.");
		}

		Point max_bound{ Pts[0][0], Pts[0][1] }, min_bound{ Pts[0][0], Pts[0][1] };
		for (const auto& m : Pts)
		{
			max_bound[0] = std::max(m[0], max_bound[0]);
			max_bound[1] = std::max(m[1], max_bound[1]);
			min_bound[0] = std::min(m[0], min_bound[0]);
			min_bound[1] = std::min(m[1], min_bound[1]);
		}

		std::vector<int> index;
		std::vector<Point> sorted(Pts);
		bool hasDuplicateElem = RemoveDuplicateAndSort(sorted, index, std::less<Point>());
		if (hasDuplicateElem)
			warning("DuplicateElem Detected, output triangulation may not use all the input points.");

		int n = sorted.size();

		std::vector<Point> Convert(n + 4);
		std::copy(sorted.begin(), sorted.end(), Convert.begin());

		// 生成全包三角形(两个)
		Convert.at(n + 0) = Point{ (max_bound[0] + min_bound[0]) / 2, (3 * max_bound[1] - min_bound[1]) / 2 };
		Convert.at(n + 1) = Point{ (3 * min_bound[0] - max_bound[0]) / 2, (max_bound[1] + min_bound[1]) / 2 };
		Convert.at(n + 2) = Point{ (3 * max_bound[0] - min_bound[0]) / 2, (max_bound[1] + min_bound[1]) / 2 };
		Convert.at(n + 3) = Point{ (max_bound[0] + min_bound[0]) / 2, (3 * min_bound[1] - max_bound[1]) / 2 };

#ifdef GenMatlabCmd
		MatlabCmd += "myline = @(p1,p2,opt_title,opt_content)line([p1(1);p2(1)], [p1(2);p2(2)], opt_title, opt_content);\n";
		MatlabCmd += "figure; hold on; axis equal; grid on;title('this algorithm')\n";
		MatlabCmd += "p=[";
		for (int i = 0; i < n + 4; ++i)
		{
			MatlabCmd += std::to_string(Convert[i][0]) + ", " + std::to_string(Convert[i][1]) + ";\n";
		}
		MatlabCmd += "];\nplot(p(1:end-4,1),p(1:end-4,2),'o');\nplot(p(end-3:end,1),p(end-3:end,2),'ro');\n";
#endif

		std::vector<std::array<int, 3>> TmpTri;

		TmpTri.emplace_back(std::array<int, 3>{ n + 0, n + 1, n + 2 });
		TmpTri.emplace_back(std::array<int, 3>{ n + 3, n + 1, n + 2 });

		// 遍历所有三角形
		for (int ith_pts = 0; ith_pts < n; ++ith_pts)
		{
			std::vector<int> ToDelete;
			std::vector<Edge> BadEdges;
			for (int ith_tri = 0; ith_tri < TmpTri.size(); ++ith_tri)
			{
				auto cmp_result = isInsideTriangleCircumCircle(
					Convert.at(ith_pts),
					Convert.at(TmpTri[ith_tri][0]),
					Convert.at(TmpTri[ith_tri][1]),
					Convert.at(TmpTri[ith_tri][2])
				);
				if (cmp_result == WeakAccept)
				{
					ToDelete.push_back(ith_tri);

					std::sort(TmpTri[ith_tri].data(), TmpTri[ith_tri].data() + 3);

					BadEdges.emplace_back(TmpTri[ith_tri][0], TmpTri[ith_tri][1]);
					BadEdges.emplace_back(TmpTri[ith_tri][1], TmpTri[ith_tri][2]);
					BadEdges.emplace_back(TmpTri[ith_tri][0], TmpTri[ith_tri][2]);
				}
				else if (cmp_result == StrongReject)
				{
					// this is now a stable triangle.
					ToDelete.push_back(ith_tri);
					Tri.push_back(TmpTri[ith_tri]);
				}
			}
			RemoveDuplicateAndSort(ToDelete);
			for (const auto& m : ToDelete)
				TmpTri.erase(TmpTri.begin() + m);

			ToDelete.clear();
			// 删除重复边，得到多边形边界，此时必为一个凸的(??)，而插入点将分割它。
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
			auto ori = BadEdges;
			RemoveDuplicateAndSort(ToDelete);
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

#ifdef DBCHECK
			error(!check_circular(BadEdges), "BadEdges are bad." + EdgeToMatlabCmd(ori));
#endif
			for (const auto& m : BadEdges)
			{
				TmpTri.emplace_back(std::array<int, 3>{m.first, m.second, ith_pts});
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

		Tri.insert(Tri.end(), TmpTri.begin(), TmpTri.end());
#ifdef GenMatlabCmd
		MatlabCmd += "figure;res = delaunayTriangulation(p);triplot(res.ConnectivityList, p(:,1), p(:,2));title('official algorithm from matlab');axis equal;grid on;\n";
#endif
		//删除辅助三角形
		Tri.erase(std::remove_if(Tri.begin(), Tri.end(),
			[n](const decltype(Tri)::value_type &t) {
			return std::any_of(t.begin(), t.end(), [&n](const typename std::decay<decltype(t)>::type::value_type& v) {return v >= n; });
		}
		), Tri.end());

		for (int i = 0; i < Tri.size(); ++i)
		{
			for (int j = 0; j < Tri[0].size(); ++j)
				Tri[i][j] = index[Tri[i][j]];
		}

		ForceSameOrientation();

#ifdef GenMatlabCmd
		MatlabCmd += "figure; hold on; axis equal; grid on;title('final layout')\n";
		MatlabCmd += "p_ori=[";
		for (int i = 0; i < n; ++i)
		{
			MatlabCmd += std::to_string(Pts[i][0]) + ", " + std::to_string(Pts[i][0]) + ";\n";
		}
		MatlabCmd += "];\nplot(p_ori(:,1),p_ori(:,2),'o');\n";
		MatlabCmd += "tri=[";
		for (int i = 0; i < this->Tri.size(); ++i)
		{
			MatlabCmd += std::to_string(1 + this->Tri[i][0]) + ", " + std::to_string(1 + this->Tri[i][1]) + ", " + std::to_string(1 + this->Tri[i][2]) + ";\n";
		}
		MatlabCmd += "];\n";
		MatlabCmd += "for i = 1:size(p_ori, 1)\ntext(p_ori(i, 1), p_ori(i, 2), num2str(i));\nend\n";
		MatlabCmd += "for i = 1:size(tri, 1)\n";
		MatlabCmd += "myline(p_ori(tri(i, 1), :), p_ori(tri(i, 2), :), 'Color', [0, 0, 1]);\n";
		MatlabCmd += "myline(p_ori(tri(i, 1), :), p_ori(tri(i, 3), :), 'Color', [0, 0, 1]);\n";
		MatlabCmd += "myline(p_ori(tri(i, 2), :), p_ori(tri(i, 3), :), 'Color', [0, 0, 1]);\nend\n";
#endif
	}



	// 做成图遍历并调整
	void Delaunay2D::ForceSameOrientation()
	{
		enum check_status { not_visited, visited_need_flip, visited_normal };
		std::vector<check_status> flip_status(Tri.size(), not_visited);
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
						if (Tri[now_visiting][j] == Tri[now_comparing][k])
						{
							now_visiting_ind[now_visiting_ind_pos++] = j;
							now_comparing_ind[now_comparing_ind_pos++] = k;
						}
				check_status for_update;
				if ((now_visiting_ind[1] - now_visiting_ind[0] + 3) % 3 == (now_comparing_ind[1] - now_comparing_ind[0] + 3) % 3)
				{
					if (flip_status[now_visiting] == visited_normal)
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
#ifdef DBCHECK
					if (for_update != flip_status[now_comparing])
						error("may be Mobius, check your input points.");//对本例不会出现这种问题,如果出错说明上面的求解有bug
#endif
				}
			}
		}
		// flip all that need to flip
		for (int i = 0; i < Tri.size(); ++i)
		{
			if (flip_status[i] == visited_need_flip)
			{
				using std::swap;
				swap(Tri[i][1], Tri[i][2]);
			}
		}
	}

	const std::vector<Triangle>& Delaunay2D::GetTri() const
	{
		return Tri;
	}

#ifdef DELAUNAY_2D_USE_EIGEN
	Eigen::MatrixX3i Delaunay2D::GetTri_Eigen() const
	{
		Eigen::MatrixX3i ret;
		ret.resize(Tri.size(), 3);
		for (int i = 0; i < Tri.size(); ++i)
		{
			ret.row(i) = Eigen::RowVector3i(Tri[i][0], Tri[i][1], Tri[i][2]);
		}
		return ret;
	}
#endif

#ifdef GenMatlabCmd
	std::string Delaunay2D::GetMatlabCmd() const
	{
		return MatlabCmd;
	}
#endif


}

