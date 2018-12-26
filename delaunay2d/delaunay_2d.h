#pragma once

#include <Eigen/Eigen>
#include <string>
#include <utility>
#include <numeric>
//#define GenMatlabCmd

typedef std::pair<int, int> Edge;

class Delaunay2D
{
#ifdef GenMatlabCmd
	std::string MatlabCmd;
#endif
	void assign(const Eigen::MatrixX2d& PtsToInsert);
	void assign(const std::vector<std::pair<double, double>>& Pts);
	void ForceSameOrientation();
	Eigen::MatrixX3i Tri;
public:
	Delaunay2D(const Eigen::MatrixX2d& PtsToInsert);
	Eigen::MatrixX3i GetTri() const;
#ifdef GenMatlabCmd
	std::string GetMatlabCmd() const;
#endif
};


enum CircleCmpStatus
{
	StrongReject,
	WeakReject,
	WeakAccept
};


template<typename scalar>
inline scalar dot(const std::pair<scalar, scalar>& pts1, const std::pair<scalar, scalar>& pts2)
{
	return pts1.first * pts2.first + pts1.second * pts2.second;
}

template<typename scalar>
inline scalar dot(const std::pair<scalar, scalar>& pts)
{
	return dot(pts, pts);
}

template<typename scalar>
inline std::pair<scalar, scalar> operator-(const std::pair<scalar, scalar>& pts1, const std::pair<scalar, scalar>& pts2)
{
	return std::make_pair(pts1.first - pts2.first, pts1.second - pts2.second);
}


//compared by a sorted points (by x dim), here we treat x axis specially.
//由于本作按照x排好序了，比较的时候会额外比较x的情况
template<typename scalar>
CircleCmpStatus isInsideTriangleCircumCircle(const std::pair<scalar, scalar>& pts_now,
	const std::pair<scalar, scalar>& p0, const std::pair<scalar, scalar>& p1, const std::pair<scalar, scalar>& p2)
{
	// 解方程
	double p0_v = dot(p0);
	std::pair<scalar, scalar> delta_p_val(p0_v - dot(p1), p0_v - dot(p2));

	std::pair<scalar, scalar>
		delta_vec_p01 = (p0 - p1),
		delta_vec_p02 = (p0 - p2);

	double checker = 0.5 / (delta_vec_p01.first * delta_vec_p02.second - delta_vec_p01.second * delta_vec_p02.first);

	std::pair<scalar, scalar> center(
		checker * (delta_p_val.first * delta_vec_p02.second - delta_p_val.second * delta_vec_p01.second),
		checker * (delta_p_val.second * delta_vec_p01.first - delta_p_val.first * delta_vec_p02.first));

	double r = dot(p1 - center);

	// center 与 pts 的 x 比较, 假定已经按照x轴排序了
	// compare pts and center on x axis. assume that points are sorted by x axis.
	if (std::pow(pts_now.first - center.first, 2) > r)
		return StrongReject;

	if (dot(pts_now - center) > r)
		return WeakReject;
	
	return WeakAccept;
}


// returns if no duplicate element is removed.
template<typename ElemTy, typename Pred = std::greater<ElemTy>>
bool RemoveDuplicateAndSort(std::vector<ElemTy>& ToProcess, Pred pred = Pred())
{
	size_t n = ToProcess.size();
	std::sort(ToProcess.begin(), ToProcess.end(), pred);
	auto last = std::unique(ToProcess.begin(), ToProcess.end());
	ToProcess.erase(last, ToProcess.end());
	return (last - ToProcess.begin()) == n;
}

template<typename ElemTy, typename Pred = std::greater<ElemTy>>
bool RemoveDuplicateAndSort(std::vector<ElemTy>& ToProcess, std::vector<int>& index, Pred pred = Pred())
{
	size_t n = ToProcess.size();
	index.resize(ToProcess.size());
	std::iota(index.begin(), index.end(), 0);

	std::vector<std::pair<ElemTy, int>> combined;
	combined.reserve(ToProcess.size());
	for (int i = 0; i < ToProcess.size(); ++i)
		combined.emplace_back(ToProcess[i], index[i]);

	std::sort(combined.begin(), combined.end(),
		[&pred](std::pair<ElemTy, int>& v1, std::pair<ElemTy, int>& v2) {return pred(v1.first, v2.first); }
	);

	long long finish = std::unique(combined.begin(), combined.end(),
		[](std::pair<ElemTy, int>& v1, std::pair<ElemTy, int>& v2) {return v1.first == v2.first; }
	) - combined.begin();

	ToProcess.erase(ToProcess.begin() + finish, ToProcess.end());
	index.erase(index.begin() + finish, index.end());
	for (int i = 0; i < finish; ++i)
	{
		ToProcess[i] = combined[i].first;
		index[i] = combined[i].second;
	}
	return (finish) == n;
}

template<typename scalar>
std::vector<std::pair<scalar, scalar>> MatToStdVec(const Eigen::Matrix<scalar, -1, 2>& In)
{
	std::vector<std::pair<scalar, scalar>> ret;
	ret.reserve(In.rows());
	for (int i = 0; i < In.rows(); ++i)
	{
		ret.emplace_back(In(i, 0), In(i, 1));
	}
	return ret;
}